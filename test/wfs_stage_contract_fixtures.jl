# Test-only concrete implementations of the public prepared WFS stage protocol.
# These fixtures intentionally contain no production WFS-family behavior.

struct ContractPlaneDevice <: AdaptiveOpticsSim.AbstractPlaneDevice
    identifier::Int
end

struct ContractDeviceArray{T,N,A<:Array{T,N}} <: AbstractArray{T,N}
    storage::A
    device::ContractPlaneDevice
end

Base.size(storage::ContractDeviceArray) = size(storage.storage)
Base.getindex(storage::ContractDeviceArray, indices...) =
    getindex(storage.storage, indices...)
Base.setindex!(storage::ContractDeviceArray, value, indices...) =
    setindex!(storage.storage, value, indices...)
Base.IndexStyle(::Type{<:ContractDeviceArray}) = IndexLinear()
AdaptiveOpticsSim.array_backend_selector(::Type{<:ContractDeviceArray}) =
    CPUBackend()
AdaptiveOpticsSim.plane_device(storage::ContractDeviceArray) = storage.device

struct ContractRateModel{T<:AbstractFloat,M}
    scale::T
    opd_gain::T
    output_metadata::M
end

struct ContractRatePlan{M,I,O}
    model::M
    input::I
    output::O
end

function _contract_require_optical_domains(input, output)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "contract fixture input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:optical_formation, :device,
            "contract fixture input and output devices differ"))
    return nothing
end

function AdaptiveOpticsSim.prepare_wfs_optical_formation(
    model::ContractRateModel, input, output::IntensityMap)
    AdaptiveOpticsSim.validate_wfs_optical_input(input)
    AdaptiveOpticsSim.validate_wfs_optical_products(output)
    output.metadata === model.output_metadata ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "contract fixture output metadata does not match its model"))
    input.metadata.dimensions == output.metadata.dimensions ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "contract fixture input and output dimensions differ"))
    _contract_require_optical_domains(input, output)
    return ContractRatePlan(model, input, output)
end

@inline function _contract_require_rate_binding(output::IntensityMap,
    input::PupilFunction, plan::ContractRatePlan)
    output.metadata === plan.output.metadata &&
        output.values === plan.output.values ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "contract fixture rate output does not match its prepared storage"))
    input.metadata === plan.input.metadata &&
        input.amplitude === plan.input.amplitude &&
        input.opd === plan.input.opd ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "contract fixture pupil input does not match its prepared storage"))
    return nothing
end

@inline function _contract_require_rate_binding(output::IntensityMap,
    input::ElectricField, plan::ContractRatePlan)
    output.metadata === plan.output.metadata &&
        output.values === plan.output.values ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "contract fixture rate output does not match its prepared storage"))
    input.metadata === plan.input.metadata && input.values === plan.input.values ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "contract fixture electric-field input does not match its prepared storage"))
    return nothing
end

function AdaptiveOpticsSim.form_wfs_optical_products!(output::IntensityMap,
    input::PupilFunction, plan::ContractRatePlan)
    _contract_require_rate_binding(output, input, plan)
    scale = plan.model.scale
    opd_gain = plan.model.opd_gain
    @. output.values = scale * (abs2(input.amplitude) +
        opd_gain * abs(input.opd))
    return output
end

function AdaptiveOpticsSim.form_wfs_optical_products!(output::IntensityMap,
    input::ElectricField, plan::ContractRatePlan)
    _contract_require_rate_binding(output, input, plan)
    scale = plan.model.scale
    @. output.values = scale * abs2(input.values)
    return output
end

function AdaptiveOpticsSim.validate_wfs_optical_formation_binding(
    output::IntensityMap, input, plan::ContractRatePlan)
    _contract_require_rate_binding(output, input, plan)
    return nothing
end

struct ContractBundleRateModel{M<:Tuple}
    models::M
end

struct ContractBundleRatePlan{P<:Tuple,I,O}
    plans::P
    input::I
    output::O
end

@inline _contract_prepare_rate_plans(::Tuple{}, input, ::Tuple{}) = ()

function _contract_prepare_rate_plans(::Tuple{}, input, outputs::Tuple)
    throw(WFSPreparationError(:optical_formation, :plane_count,
        "contract fixture has more optical outputs than models"))
end

function _contract_prepare_rate_plans(models::Tuple, input, ::Tuple{})
    throw(WFSPreparationError(:optical_formation, :plane_count,
        "contract fixture has fewer optical outputs than models"))
end

function _contract_prepare_rate_plans(models::Tuple, input, outputs::Tuple)
    plan = prepare_wfs_optical_formation(first(models), input, first(outputs))
    return (plan,
        _contract_prepare_rate_plans(Base.tail(models), input,
            Base.tail(outputs))...)
end

function AdaptiveOpticsSim.prepare_wfs_optical_formation(
    model::ContractBundleRateModel, input, output::OpticalProductBundle)
    AdaptiveOpticsSim.validate_wfs_optical_input(input)
    AdaptiveOpticsSim.validate_wfs_optical_products(output)
    plans = _contract_prepare_rate_plans(model.models, input, output.products)
    return ContractBundleRatePlan(plans, input, output)
end

@inline _contract_form_rate_plans!(::Tuple{}, ::Tuple{}, input) = nothing

@inline function _contract_form_rate_plans!(plans::Tuple, outputs::Tuple,
    input)
    form_wfs_optical_products!(first(outputs), input, first(plans))
    return _contract_form_rate_plans!(Base.tail(plans), Base.tail(outputs),
        input)
end

function AdaptiveOpticsSim.form_wfs_optical_products!(
    output::OpticalProductBundle, input, plan::ContractBundleRatePlan)
    output.products === plan.output.products ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "contract fixture optical bundle does not match its prepared leaves"))
    _contract_form_rate_plans!(plan.plans, output.products, input)
    return output
end

function AdaptiveOpticsSim.validate_wfs_optical_formation_binding(
    output::OpticalProductBundle, input, plan::ContractBundleRatePlan)
    output.products === plan.output.products || throw(WFSPreparationError(
        :optical_formation, :prepared_binding,
        "contract fixture optical bundle does not match its prepared leaves"))
    @inbounds for index in eachindex(plan.plans)
        component_input = AdaptiveOpticsSim.four_pupil_bundle_input(input,
            index)
        AdaptiveOpticsSim.validate_wfs_optical_formation_binding(
            output[index], component_input, plan.plans[index])
    end
    return nothing
end

struct ContractDetectorBinding{D,I,P,O}
    detector::D
    input::I
    plan::P
    observation::O
end

function contract_detector_binding(detector::Detector, input::IntensityMap;
    units=:detector_signal, layout=:detector_frame,
    normalized_to_photon_rate=nothing)
    plan = prepare_detector_acquisition(detector, input;
        normalized_to_photon_rate=normalized_to_photon_rate)
    observation = WFSObservation(output_frame(detector); units, layout)
    return ContractDetectorBinding(detector, input, plan, observation)
end

struct ContractDetectorAcquisitionModel{B<:Tuple}
    bindings::B
end

struct ContractDetectorAcquisitionPlan{B<:Tuple,I,O}
    bindings::B
    inputs::I
    observations::O
end

@inline _contract_contains_product(product::IntensityMap,
    target::IntensityMap) = product.metadata === target.metadata &&
    product.values === target.values

@inline _contract_contains_product(bundle::OpticalProductBundle,
    target::IntensityMap) = _contract_contains_product(bundle.products, target)

@inline _contract_contains_product(::Tuple{}, target::IntensityMap) = false

@inline function _contract_contains_product(products::Tuple,
    target::IntensityMap)
    return _contract_contains_product(first(products), target) ||
        _contract_contains_product(Base.tail(products), target)
end

@inline function _contract_require_observation_binding(binding,
    observation::WFSObservation)
    observation.storage === binding.observation.storage &&
        observation.metadata === binding.observation.metadata &&
        observation.units === binding.observation.units ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "contract fixture observation does not match its detector binding"))
    output_frame(binding.detector) === observation.storage ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "contract fixture detector output storage was replaced"))
    return nothing
end

@inline function _contract_require_single_binding(bindings::Tuple,
    observations::WFSObservation)
    length(bindings) == 1 || throw(WFSPreparationError(:acquisition,
        :plane_count, "multiple detector bindings require a tuple of observations"))
    _contract_require_observation_binding(first(bindings), observations)
    return nothing
end

@inline _contract_require_tuple_bindings(::Tuple{}, ::Tuple{}) = nothing

function _contract_require_tuple_bindings(::Tuple{}, observations::Tuple)
    throw(WFSPreparationError(:acquisition, :plane_count,
        "contract fixture has more observations than detector bindings"))
end

function _contract_require_tuple_bindings(bindings::Tuple, ::Tuple{})
    throw(WFSPreparationError(:acquisition, :plane_count,
        "contract fixture has fewer observations than detector bindings"))
end

@inline function _contract_require_tuple_bindings(bindings::Tuple,
    observations::Tuple)
    _contract_require_observation_binding(first(bindings), first(observations))
    return _contract_require_tuple_bindings(Base.tail(bindings),
        Base.tail(observations))
end

@inline _contract_require_binding_outputs(bindings::Tuple,
    observation::WFSObservation) =
    _contract_require_single_binding(bindings, observation)

@inline _contract_require_binding_outputs(bindings::Tuple,
    observations::Tuple) =
    _contract_require_tuple_bindings(bindings, observations)

@inline _contract_require_binding_inputs(::Tuple{}, products) = nothing

@inline function _contract_require_binding_inputs(bindings::Tuple, products)
    binding = first(bindings)
    _contract_contains_product(products, binding.input) ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "contract fixture detector input is absent from the optical products"))
    return _contract_require_binding_inputs(Base.tail(bindings), products)
end

@inline _contract_require_unique_detectors(::Tuple{}) = nothing

@inline function _contract_require_detector_absent(detector, ::Tuple{})
    return nothing
end

@inline function _contract_require_detector_absent(detector, bindings::Tuple)
    detector.state === first(bindings).detector.state &&
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "one detector state cannot appear in multiple acquisition bindings"))
    return _contract_require_detector_absent(detector, Base.tail(bindings))
end

@inline function _contract_require_unique_detectors(bindings::Tuple)
    binding = first(bindings)
    _contract_require_detector_absent(binding.detector, Base.tail(bindings))
    return _contract_require_unique_detectors(Base.tail(bindings))
end

function AdaptiveOpticsSim.prepare_wfs_acquisition(
    model::ContractDetectorAcquisitionModel, optical_products, observations)
    AdaptiveOpticsSim.validate_wfs_optical_products(optical_products)
    AdaptiveOpticsSim.validate_wfs_observations(observations)
    _contract_require_unique_detectors(model.bindings)
    _contract_require_binding_inputs(model.bindings, optical_products)
    _contract_require_binding_outputs(model.bindings, observations)
    return ContractDetectorAcquisitionPlan(model.bindings, optical_products,
        observations)
end

@inline _contract_require_live_detector_outputs(::Tuple{}) = nothing

@inline function _contract_require_live_detector_outputs(bindings::Tuple)
    binding = first(bindings)
    output_frame(binding.detector) === binding.observation.storage ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "contract fixture detector output storage was replaced after preparation"))
    return _contract_require_live_detector_outputs(Base.tail(bindings))
end

@inline function _contract_require_rng_arity(bindings::Tuple,
    ::AbstractRNG)
    length(bindings) == 1 || throw(WFSPreparationError(:acquisition,
        :detector_mapping,
        "multiple detector bindings require one concrete RNG per detector"))
    return nothing
end

@inline function _contract_require_rng_arity(bindings::Tuple, rngs::Tuple)
    length(bindings) == length(rngs) || throw(WFSPreparationError(
        :acquisition, :detector_mapping,
        "detector binding and RNG counts must match"))
    return nothing
end

@inline _contract_require_rng_type(::AbstractRNG) = nothing

function _contract_require_rng_type(rng)
    throw(WFSPreparationError(:acquisition, :detector_mapping,
        "each detector binding requires a concrete AbstractRNG"))
end

@inline _contract_require_rng_types(::Tuple{}) = nothing

@inline function _contract_require_rng_types(rngs::Tuple)
    _contract_require_rng_type(first(rngs))
    return _contract_require_rng_types(Base.tail(rngs))
end

@inline _contract_require_rng_types(rng::AbstractRNG) = nothing

@inline function _contract_capture_bindings!(bindings::Tuple,
    rng::AbstractRNG)
    length(bindings) == 1 || throw(WFSPreparationError(:acquisition,
        :detector_mapping,
        "multiple detector bindings require one concrete RNG per detector"))
    binding = first(bindings)
    result = capture!(binding.detector, binding.input, binding.plan, rng)
    result === binding.observation.storage ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "contract fixture detector returned unexpected output storage"))
    return nothing
end

@inline _contract_capture_bindings!(::Tuple{}, ::Tuple{}) = nothing

function _contract_capture_bindings!(::Tuple{}, rngs::Tuple)
    throw(WFSPreparationError(:acquisition, :detector_mapping,
        "contract fixture has more RNGs than detector bindings"))
end

function _contract_capture_bindings!(bindings::Tuple, ::Tuple{})
    throw(WFSPreparationError(:acquisition, :detector_mapping,
        "contract fixture has fewer RNGs than detector bindings"))
end

@inline function _contract_capture_bindings!(bindings::Tuple, rngs::Tuple)
    binding = first(bindings)
    result = capture!(binding.detector, binding.input, binding.plan,
        first(rngs))
    result === binding.observation.storage ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "contract fixture detector returned unexpected output storage"))
    return _contract_capture_bindings!(Base.tail(bindings), Base.tail(rngs))
end

function AdaptiveOpticsSim.acquire_wfs_observation!(observations,
    optical_products, plan::ContractDetectorAcquisitionPlan, rng)
    _contract_require_rng_arity(plan.bindings, rng)
    _contract_require_rng_types(rng)
    _contract_require_binding_inputs(plan.bindings, optical_products)
    _contract_require_binding_outputs(plan.bindings, observations)
    _contract_require_live_detector_outputs(plan.bindings)
    _contract_capture_bindings!(plan.bindings, rng)
    return observations
end

function AdaptiveOpticsSim.validate_wfs_acquisition_binding(observations,
    optical_products, plan::ContractDetectorAcquisitionPlan)
    _contract_require_binding_inputs(plan.bindings, optical_products)
    _contract_require_binding_outputs(plan.bindings, observations)
    _contract_require_live_detector_outputs(plan.bindings)
    return nothing
end

struct ContractSumEstimator{U,K}
    units::U
    kind::K
end

struct ContractSumEstimatorPlan{I,M}
    input::I
    measurement::M
end

@inline _contract_observation_count(::WFSObservation) = 1
@inline _contract_observation_count(observations::Tuple) = length(observations)

function _contract_require_sum_shape(input::WFSObservation,
    storage::Base.RefValue)
    return nothing
end

function _contract_require_sum_shape(input::WFSObservation,
    storage::AbstractVector)
    length(storage) == 1 || throw(WFSPreparationError(:estimation, :shape,
        "one observation requires one measurement element"))
    return nothing
end

function _contract_require_sum_shape(input::Tuple, storage::AbstractVector)
    length(storage) == length(input) || throw(WFSPreparationError(:estimation,
        :shape, "measurement length must match observation count"))
    return nothing
end

function _contract_require_sum_shape(input, storage)
    throw(WFSPreparationError(:estimation, :shape,
        "contract sum estimator requires Ref or vector output"))
end

function _contract_require_measurement_semantics(model,
    measurement::WFSMeasurement)
    isequal(measurement.units, model.units) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "measurement units are incompatible with the prepared estimator"))
    isequal(measurement.metadata.kind, model.kind) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "measurement kind is incompatible with the prepared estimator"))
    return nothing
end

function _contract_require_measurement_domain(input,
    measurement::WFSMeasurement)
    typeof(input.metadata.backend) === typeof(measurement.metadata.backend) ||
        throw(WFSPreparationError(:estimation, :backend,
            "estimator input and measurement backends differ"))
    input.metadata.device == measurement.metadata.device ||
        throw(WFSPreparationError(:estimation, :device,
            "estimator input and measurement devices differ"))
    return nothing
end

@inline _contract_require_measurement_domains(::Tuple{},
    ::WFSMeasurement) = nothing

@inline function _contract_require_measurement_domains(inputs::Tuple,
    measurement::WFSMeasurement)
    _contract_require_measurement_domain(first(inputs), measurement)
    return _contract_require_measurement_domains(Base.tail(inputs),
        measurement)
end

@inline _contract_require_measurement_domains(input,
    measurement::WFSMeasurement) =
    _contract_require_measurement_domain(input, measurement)

function _contract_require_cpu_reduction(measurement::WFSMeasurement)
    typeof(measurement.metadata.backend) === CPUBackend ||
        throw(WFSPreparationError(:estimation, :backend,
            "contract reduction estimators are CPU-only"))
    measurement.metadata.device == AdaptiveOpticsSim.HostPlaneDevice() ||
        throw(WFSPreparationError(:estimation, :device,
            "contract reduction estimators require host storage"))
    return nothing
end

function _contract_require_reduction_numeric_type(input,
    measurement::WFSMeasurement)
    output_type = measurement.metadata.numeric_type
    output_type <: AbstractFloat || throw(WFSPreparationError(
        :estimation, :estimator,
        "contract reduction estimators require real floating-point output"))
    input.metadata.numeric_type === output_type ||
        throw(WFSPreparationError(:estimation, :estimator,
            "reduction input and output numeric types must match"))
    return nothing
end

@inline _contract_require_reduction_numeric_types(::Tuple{},
    ::WFSMeasurement) = nothing

@inline function _contract_require_reduction_numeric_types(inputs::Tuple,
    measurement::WFSMeasurement)
    _contract_require_reduction_numeric_type(first(inputs), measurement)
    return _contract_require_reduction_numeric_types(Base.tail(inputs),
        measurement)
end

@inline _contract_require_reduction_numeric_types(input,
    measurement::WFSMeasurement) =
    _contract_require_reduction_numeric_type(input, measurement)

function AdaptiveOpticsSim.prepare_wfs_estimation(model::ContractSumEstimator,
    input, measurement::WFSMeasurement)
    AdaptiveOpticsSim.validate_wfs_observations(input)
    AdaptiveOpticsSim.validate_wfs_measurement(measurement)
    _contract_require_measurement_semantics(model, measurement)
    _contract_require_cpu_reduction(measurement)
    _contract_require_measurement_domains(input, measurement)
    _contract_require_reduction_numeric_types(input, measurement)
    _contract_require_sum_shape(input, measurement.storage)
    return ContractSumEstimatorPlan(input, measurement)
end

@inline function _contract_require_observation_storage(
    actual::WFSObservation, prepared::WFSObservation)
    actual.storage === prepared.storage &&
        actual.metadata === prepared.metadata && actual.units === prepared.units ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "estimator observation does not match prepared storage"))
    return nothing
end

@inline _contract_require_observation_storage(::Tuple{}, ::Tuple{}) = nothing

@inline function _contract_require_observation_storage(actual::Tuple,
    prepared::Tuple)
    _contract_require_observation_storage(first(actual), first(prepared))
    return _contract_require_observation_storage(Base.tail(actual),
        Base.tail(prepared))
end

@inline function _contract_require_measurement_storage(actual::WFSMeasurement,
    prepared::WFSMeasurement)
    actual.storage === prepared.storage &&
        actual.metadata === prepared.metadata && actual.units === prepared.units ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "estimator measurement does not match prepared storage"))
    return nothing
end


@inline function _contract_write_sums!(storage::Base.RefValue,
    observation::WFSObservation)
    storage[] = sum(observation.storage)
    return storage
end

@inline function _contract_write_sums!(storage::AbstractVector,
    observation::WFSObservation)
    @inbounds storage[1] = sum(observation.storage)
    return storage
end

@inline _contract_write_tuple_sums!(storage, ::Tuple{}, index::Int) = storage

@inline function _contract_write_tuple_sums!(storage,
    observations::Tuple, index::Int)
    @inbounds storage[index] = sum(first(observations).storage)
    return _contract_write_tuple_sums!(storage, Base.tail(observations),
        index + 1)
end

@inline _contract_write_sums!(storage::AbstractVector,
    observations::Tuple) = _contract_write_tuple_sums!(storage,
    observations, 1)

function AdaptiveOpticsSim.estimate_wfs_measurement!(
    measurement::WFSMeasurement, input,
    plan::ContractSumEstimatorPlan)
    _contract_require_observation_storage(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    _contract_write_sums!(measurement.storage, input)
    return measurement
end

AdaptiveOpticsSim.wfs_measurement_path(::ContractSumEstimatorPlan) =
    AcquiredObservationPath()

function AdaptiveOpticsSim.validate_wfs_estimation_binding(
    measurement::WFSMeasurement, input, plan::ContractSumEstimatorPlan)
    _contract_require_observation_storage(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    return nothing
end

struct ContractCopyEstimator{U,K}
    units::U
    kind::K
end

struct ContractCopyEstimatorPlan{I,M}
    input::I
    measurement::M
end

function AdaptiveOpticsSim.prepare_wfs_estimation(model::ContractCopyEstimator,
    input::WFSObservation, measurement::WFSMeasurement)
    AdaptiveOpticsSim.validate_wfs_observation(input)
    AdaptiveOpticsSim.validate_wfs_measurement(measurement)
    _contract_require_measurement_semantics(model, measurement)
    input.metadata.dimensions == measurement.metadata.dimensions ||
        throw(WFSPreparationError(:estimation, :shape,
            "copy estimator input and output dimensions differ"))
    input.metadata.numeric_type === measurement.metadata.numeric_type ||
        throw(WFSPreparationError(:estimation, :shape,
            "copy estimator input and output numeric types differ"))
    typeof(input.metadata.backend) === typeof(measurement.metadata.backend) ||
        throw(WFSPreparationError(:estimation, :backend,
            "copy estimator input and output backends differ"))
    input.metadata.device == measurement.metadata.device ||
        throw(WFSPreparationError(:estimation, :device,
            "copy estimator input and output devices differ"))
    return ContractCopyEstimatorPlan(input, measurement)
end

function AdaptiveOpticsSim.estimate_wfs_measurement!(
    measurement::WFSMeasurement, input::WFSObservation,
    plan::ContractCopyEstimatorPlan)
    _contract_require_observation_storage(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    copyto!(measurement.storage, input.storage)
    return measurement
end

AdaptiveOpticsSim.wfs_measurement_path(::ContractCopyEstimatorPlan) =
    AcquiredObservationPath()

function AdaptiveOpticsSim.validate_wfs_estimation_binding(
    measurement::WFSMeasurement, input, plan::ContractCopyEstimatorPlan)
    _contract_require_observation_storage(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    return nothing
end

struct ContractDirectEstimator{T<:AbstractFloat,U,K}
    scale::T
    units::U
    kind::K
end

struct ContractDirectEstimatorPlan{M,I,O}
    model::M
    input::I
    measurement::O
end

@inline _contract_direct_result_type(input::PupilFunction) =
    eltype(input.opd)
@inline _contract_direct_result_type(input::ElectricField) =
    typeof(real(zero(eltype(input.values))))

function _contract_require_direct_numeric_type(input,
    measurement::WFSMeasurement)
    expected = _contract_direct_result_type(input)
    measurement.metadata.numeric_type === expected ||
        throw(WFSPreparationError(:estimation, :estimator,
            "direct-estimator output numeric type is incompatible with its input"))
    return nothing
end

function AdaptiveOpticsSim.prepare_wfs_estimation(
    model::ContractDirectEstimator, input,
    measurement::WFSMeasurement)
    AdaptiveOpticsSim.validate_wfs_optical_input(input)
    AdaptiveOpticsSim.validate_wfs_measurement(measurement)
    _contract_require_measurement_semantics(model, measurement)
    _contract_require_cpu_reduction(measurement)
    _contract_require_measurement_domain(input, measurement)
    _contract_require_direct_numeric_type(input, measurement)
    measurement.storage isa AbstractVector ||
        throw(WFSPreparationError(:estimation, :shape,
            "direct contract estimator requires vector storage"))
    length(measurement.storage) == 2 ||
        throw(WFSPreparationError(:estimation, :shape,
            "direct contract estimator requires two output values"))
    return ContractDirectEstimatorPlan(model, input, measurement)
end

@inline function _contract_require_direct_input(actual::PupilFunction,
    prepared::PupilFunction)
    actual.metadata === prepared.metadata && actual.amplitude === prepared.amplitude &&
        actual.opd === prepared.opd || throw(WFSPreparationError(:estimation,
        :prepared_binding, "direct estimator pupil input was replaced"))
    return nothing
end

@inline function _contract_require_direct_input(actual::ElectricField,
    prepared::ElectricField)
    actual.metadata === prepared.metadata && actual.values === prepared.values ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "direct estimator field input was replaced"))
    return nothing
end

function AdaptiveOpticsSim.estimate_wfs_measurement!(
    measurement::WFSMeasurement, input::PupilFunction,
    plan::ContractDirectEstimatorPlan)
    _contract_require_direct_input(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    scale = plan.model.scale
    @inbounds begin
        measurement.storage[1] = scale * sum(input.opd)
        measurement.storage[2] = scale * sum(input.amplitude)
    end
    return measurement
end

function AdaptiveOpticsSim.estimate_wfs_measurement!(
    measurement::WFSMeasurement, input::ElectricField,
    plan::ContractDirectEstimatorPlan)
    _contract_require_direct_input(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    scale = plan.model.scale
    @inbounds begin
        measurement.storage[1] = scale * sum(real, input.values)
        measurement.storage[2] = scale * sum(imag, input.values)
    end
    return measurement
end

AdaptiveOpticsSim.wfs_measurement_path(::ContractDirectEstimatorPlan) =
    DirectMeasurementPath()

function AdaptiveOpticsSim.validate_wfs_estimation_binding(
    measurement::WFSMeasurement, input, plan::ContractDirectEstimatorPlan)
    _contract_require_direct_input(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    return nothing
end

struct ContractDirectCopyEstimator{U,K}
    units::U
    kind::K
end

struct ContractDirectCopyEstimatorPlan{I,M}
    input::I
    measurement::M
end

function AdaptiveOpticsSim.prepare_wfs_estimation(
    model::ContractDirectCopyEstimator, input::ElectricField,
    measurement::WFSMeasurement)
    AdaptiveOpticsSim.validate_wfs_optical_input(input)
    AdaptiveOpticsSim.validate_wfs_measurement(measurement)
    _contract_require_measurement_semantics(model, measurement)
    input.metadata.dimensions == measurement.metadata.dimensions ||
        throw(WFSPreparationError(:estimation, :shape,
            "direct-copy input and measurement dimensions differ"))
    typeof(input.metadata.backend) === typeof(measurement.metadata.backend) ||
        throw(WFSPreparationError(:estimation, :backend,
            "direct-copy input and measurement backends differ"))
    input.metadata.device == measurement.metadata.device ||
        throw(WFSPreparationError(:estimation, :device,
            "direct-copy input and measurement devices differ"))
    _contract_require_direct_numeric_type(input, measurement)
    return ContractDirectCopyEstimatorPlan(input, measurement)
end

function AdaptiveOpticsSim.estimate_wfs_measurement!(
    measurement::WFSMeasurement, input::ElectricField,
    plan::ContractDirectCopyEstimatorPlan)
    _contract_require_direct_input(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    @. measurement.storage = abs2(input.values)
    return measurement
end

AdaptiveOpticsSim.wfs_measurement_path(::ContractDirectCopyEstimatorPlan) =
    DirectMeasurementPath()

function AdaptiveOpticsSim.validate_wfs_estimation_binding(
    measurement::WFSMeasurement, input, plan::ContractDirectCopyEstimatorPlan)
    _contract_require_direct_input(input, plan.input)
    _contract_require_measurement_storage(measurement, plan.measurement)
    return nothing
end

struct ContractPackedAcquisition{R,T<:AbstractFloat}
    regions::R
    duration::T
end

struct ContractPackedAcquisitionPlan{M,I,O,V}
    model::M
    inputs::I
    observation::O
    views::V
end

function _contract_two_products(bundle::OpticalProductBundle)
    return _contract_two_products(bundle.products)
end

function _contract_two_products(products::Tuple{A,B}) where {A,B}
    return products
end

function _contract_two_products(products)
    throw(WFSPreparationError(:acquisition, :plane_count,
        "packed contract acquisition requires exactly two optical products"))
end

function _contract_validate_regions(regions, dimensions)
    length(regions) == 2 || throw(WFSPreparationError(:acquisition,
        :plane_count, "packed contract acquisition requires two regions"))
    coverage = falses(dimensions)
    for region in regions
        rows, cols = region
        first(rows) >= 1 && last(rows) <= dimensions[1] &&
            first(cols) >= 1 && last(cols) <= dimensions[2] ||
            throw(WFSPreparationError(:acquisition, :detector_mapping,
                "packed observation region is outside its storage"))
        @inbounds for col in cols, row in rows
            coverage[row, col] && throw(WFSPreparationError(:acquisition,
                :detector_mapping, "packed observation regions overlap"))
            coverage[row, col] = true
        end
    end
    all(coverage) || throw(WFSPreparationError(:acquisition,
        :detector_mapping,
        "packed observation regions must cover the declared storage"))
    return nothing
end

@inline _contract_require_packed_cell_measure(::CellIntegratedMeasure) =
    nothing

function _contract_require_packed_cell_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:acquisition, :radiometry,
        "packed contract acquisition requires cell-integrated photon rates"))
end

function AdaptiveOpticsSim.prepare_wfs_acquisition(
    model::ContractPackedAcquisition, optical_products,
    observation::WFSObservation)
    AdaptiveOpticsSim.validate_wfs_optical_products(optical_products)
    AdaptiveOpticsSim.validate_wfs_observation(observation)
    isfinite(model.duration) && model.duration > zero(model.duration) ||
        throw(WFSPreparationError(:acquisition, :duration,
            "packed acquisition duration must be finite and positive"))
    inputs = _contract_two_products(optical_products)
    observation.storage isa AbstractMatrix ||
        throw(WFSPreparationError(:acquisition, :shape,
            "packed observation storage must be a matrix"))
    observation.metadata.layout == model.regions ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "packed observation layout does not match acquisition regions"))
    _contract_validate_regions(model.regions, size(observation.storage))
    first_input = first(inputs)
    second_input = first(Base.tail(inputs))
    _contract_require_packed_cell_measure(first_input.metadata.spatial_measure)
    _contract_require_packed_cell_measure(second_input.metadata.spatial_measure)
    output_type = observation.metadata.numeric_type
    output_type <: AbstractFloat || throw(WFSPreparationError(
        :acquisition, :detector_mapping,
        "packed observation storage must use a real floating-point type"))
    first_input.metadata.numeric_type === output_type &&
        second_input.metadata.numeric_type === output_type &&
        typeof(model.duration) === output_type ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "packed input, duration, and observation numeric types must match"))
    first_region = first(model.regions)
    second_region = first(Base.tail(model.regions))
    size(first_input.values) == map(length, first_region) ||
        throw(WFSPreparationError(:acquisition, :shape,
            "first optical product does not fit its packed region"))
    size(second_input.values) == map(length, second_region) ||
        throw(WFSPreparationError(:acquisition, :shape,
            "second optical product does not fit its packed region"))
    typeof(first_input.metadata.backend) ===
        typeof(observation.metadata.backend) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "packed input and observation backends differ"))
    first_input.metadata.device == observation.metadata.device ||
        throw(WFSPreparationError(:acquisition, :device,
            "packed input and observation devices differ"))
    typeof(second_input.metadata.backend) ===
        typeof(observation.metadata.backend) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "packed input and observation backends differ"))
    second_input.metadata.device == observation.metadata.device ||
        throw(WFSPreparationError(:acquisition, :device,
            "packed input and observation devices differ"))
    first_view = view(observation.storage, first_region...)
    second_view = view(observation.storage, second_region...)
    return ContractPackedAcquisitionPlan(model, inputs, observation,
        (first_view, second_view))
end

function AdaptiveOpticsSim.acquire_wfs_observation!(
    observation::WFSObservation, optical_products,
    plan::ContractPackedAcquisitionPlan, rng)
    AdaptiveOpticsSim.validate_wfs_acquisition_binding(observation,
        optical_products, plan)
    inputs = _contract_two_products(optical_products)
    first_view = first(plan.views)
    second_view = first(Base.tail(plan.views))
    first_values = first(inputs).values
    second_values = first(Base.tail(inputs)).values
    duration = plan.model.duration
    @. first_view = first_values * duration
    @. second_view = second_values * duration
    return observation
end


function AdaptiveOpticsSim.validate_wfs_acquisition_binding(
    observation::WFSObservation, optical_products,
    plan::ContractPackedAcquisitionPlan)
    inputs = _contract_two_products(optical_products)
    first(inputs).values === first(plan.inputs).values &&
        first(Base.tail(inputs)).values ===
            first(Base.tail(plan.inputs)).values ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "packed optical input storage was replaced"))
    observation.storage === plan.observation.storage &&
        observation.metadata === plan.observation.metadata ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "packed observation storage was replaced"))
    return nothing
end

function run_contract_stages!(optical_output, optical_input, optical_plan,
    observation, acquisition_plan, rng, measurement, estimator_plan)
    form_wfs_optical_products!(optical_output, optical_input, optical_plan)
    acquire_wfs_observation!(observation, optical_output, acquisition_plan,
        rng)
    estimate_wfs_measurement!(measurement, observation, estimator_plan)
    return measurement
end
