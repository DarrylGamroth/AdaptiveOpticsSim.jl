struct TestIlluminationDefinition{C}
    combination::C
end

struct TestIlluminationEvaluator{D,C}
    destination::D
    combination::C
end

@inline AdaptiveOpticsSim.illumination_combination(
    evaluator::TestIlluminationEvaluator) = evaluator.combination

function AdaptiveOpticsSim.prepare_illumination_evaluator(
    definition::TestIlluminationDefinition, destination, boundary)
    return TestIlluminationEvaluator(destination, definition.combination)
end

function AdaptiveOpticsSim.validate_illumination_evaluator_binding(
    evaluator::TestIlluminationEvaluator, destination, boundary)
    evaluator.destination === destination || throw(PlantPreparationError(
        :illumination, :prepared_binding,
        "test illumination evaluator lost its destination binding"))
    return nothing
end

@inline function AdaptiveOpticsSim.evaluate_illumination!(destination,
    evaluator::TestIlluminationEvaluator, model_time, rng::AbstractRNG)
    return destination
end

struct MissingCombinationIlluminationDefinition end
struct MissingCombinationIlluminationEvaluator{D}
    destination::D
end

AdaptiveOpticsSim.prepare_illumination_evaluator(
    ::MissingCombinationIlluminationDefinition, destination, boundary) =
    MissingCombinationIlluminationEvaluator(destination)

function AdaptiveOpticsSim.validate_illumination_evaluator_binding(
    evaluator::MissingCombinationIlluminationEvaluator, destination,
    boundary)
    return nothing
end

struct InvalidResultIlluminationDefinition end
struct InvalidResultIlluminationEvaluator{D}
    destination::D
end

AdaptiveOpticsSim.prepare_illumination_evaluator(
    ::InvalidResultIlluminationDefinition, destination, boundary) =
    InvalidResultIlluminationEvaluator(destination)
AdaptiveOpticsSim.illumination_combination(
    ::Type{<:InvalidResultIlluminationEvaluator}) = SingleIllumination()
AdaptiveOpticsSim.validate_illumination_evaluator_binding(
    evaluator::InvalidResultIlluminationEvaluator, destination,
    boundary) = nothing
AdaptiveOpticsSim.evaluate_illumination!(destination,
    evaluator::InvalidResultIlluminationEvaluator, model_time,
    rng::AbstractRNG) = nothing

mutable struct MutableIlluminationDefinition
    value::Float64
end

struct MutableEvaluatorIlluminationDefinition end

mutable struct MutableIlluminationEvaluator{D}
    destination::D
end

AdaptiveOpticsSim.prepare_illumination_evaluator(
    ::MutableEvaluatorIlluminationDefinition, destination, boundary) =
    MutableIlluminationEvaluator(destination)

struct TimedPupilIllumination{T<:AbstractFloat}
    offset_m::T
    time_scale_m_per_s::T
    noise_scale_m::T
end

mutable struct TimedPupilIlluminationState{T,S,A}
    evaluations::Int
    last_time_s::T
    support_template::S
    amplitude_template::A
end

struct PreparedTimedPupilIllumination{P,S,B,D}
    params::P
    state::S
    backend::B
    device::D
end

AdaptiveOpticsSim.illumination_combination(
    ::Type{<:PreparedTimedPupilIllumination}) =
    ExclusiveIlluminationSelection()

function AdaptiveOpticsSim.prepare_illumination_evaluator(
    definition::TimedPupilIllumination, destination::PupilFunction,
    ::PupilFunctionIlluminationEntry)
    T = eltype(destination.opd)
    params = TimedPupilIllumination(T(definition.offset_m),
        T(definition.time_scale_m_per_s), T(definition.noise_scale_m))
    state = TimedPupilIlluminationState(0, zero(T),
        copy(destination.support), copy(destination.amplitude))
    return PreparedTimedPupilIllumination(params, state,
        backend(destination), destination.metadata.device)
end

function AdaptiveOpticsSim.validate_illumination_evaluator_binding(
    evaluator::PreparedTimedPupilIllumination,
    destination::PupilFunction, ::PupilFunctionIlluminationEntry)
    typeof(backend(destination)) === typeof(evaluator.backend) || throw(
        PlantPreparationError(:illumination, :backend,
            "timed pupil evaluator backend changed"))
    destination.metadata.device == evaluator.device || throw(
        PlantPreparationError(:illumination, :device,
            "timed pupil evaluator device changed"))
    size(evaluator.state.support_template) == size(destination.support) ||
        throw(PlantPreparationError(:illumination, :shape,
            "timed pupil support template shape changed"))
    size(evaluator.state.amplitude_template) == size(destination.amplitude) ||
        throw(PlantPreparationError(:illumination, :shape,
            "timed pupil amplitude template shape changed"))
    return nothing
end

function AdaptiveOpticsSim.evaluate_illumination!(
    destination::PupilFunction,
    evaluator::PreparedTimedPupilIllumination, model_time,
    rng::AbstractRNG)
    state = evaluator.state
    params = evaluator.params
    T = eltype(destination.opd)
    copyto!(destination.support, state.support_template)
    copyto!(destination.amplitude, state.amplitude_template)
    state.evaluations += 1
    state.last_time_s = T(model_time)
    value = params.offset_m + params.time_scale_m_per_s * T(model_time) +
        params.noise_scale_m * randn(rng, T)
    fill!(destination.opd, value)
    return destination
end

struct NativeDetectorIlluminationPathModel{T<:AbstractFloat}
    rate::T
end

struct TimedPupilIlluminationPathModel{T<:AbstractFloat}
    params::TimedPupilIllumination{T}
end

struct IlluminationFrameAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
    quantum_efficiency::T
end

struct IlluminationIdentityExecution{P}
    product::P
end

for model in (NativeDetectorIlluminationPathModel,
    TimedPupilIlluminationPathModel, IlluminationFrameAcquisitionModel)
    @eval AdaptiveOpticsSim.plant_model_definition_style(
        ::Type{<:$model}) = ColdPlantModelDefinition()
end

function AdaptiveOpticsSim.validate_path_execution_binding(
    execution::IlluminationIdentityExecution, input, result)
    execution.product === input === result || throw(PlantPreparationError(
        :path, :prepared_binding,
        "illumination identity execution lost its product binding"))
    return nothing
end

function AdaptiveOpticsSim.execute_path!(result, input,
    execution::IlluminationIdentityExecution)
    AdaptiveOpticsSim.validate_path_execution_binding(execution, input,
        result)
    return result
end

function illumination_detector_input(telescope::Telescope,
    source::AbstractSource)
    T = eltype(pupil_reflectivity(telescope))
    values = similar(pupil_reflectivity(telescope), T,
        telescope.params.resolution, telescope.params.resolution)
    fill!(values, zero(T))
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=MetricCoordinates(),
        sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(source))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

function AdaptiveOpticsSim.prepare_path_executor(
    model::NativeDetectorIlluminationPathModel,
    definition::OpticalPathDefinition, source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere)
    destination = illumination_detector_input(telescope, source)
    entry = prepare_illumination_entry(
        UniformIntensityIllumination(model.rate;
            combination=SingleIllumination()),
        destination, DetectorInputIlluminationEntry();
        visibility=(downstream_path=path_id(definition),
            starts_at=:detector_input))
    execution = IlluminationIdentityExecution(destination)
    return PreparedPathExecutor(definition, source, telescope, atmosphere,
        destination, destination, execution;
        materialization=entry,
        optical_model=(kind=:native_uniform_detector_illumination,
            photon_rate=model.rate),
        propagation_model=:detector_input_identity,
        model_revisions=UInt(1))
end

function AdaptiveOpticsSim.prepare_path_executor(
    model::TimedPupilIlluminationPathModel,
    definition::OpticalPathDefinition, source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere)
    pupil = PupilFunction(telescope; backend=backend(telescope))
    visibility = Dict{Symbol,Any}(
        :downstream_path => path_id(definition).name,
        :starts_at => :pupil_function)
    entry = prepare_illumination_entry(model.params, pupil,
        PupilFunctionIlluminationEntry(); visibility)
    imaging = prepare_direct_imaging(pupil, source; zero_padding=2)
    return PreparedPathExecutor(definition, source, telescope, atmosphere,
        pupil, direct_imaging_output(imaging), imaging;
        materialization=entry,
        optical_model=(kind=:timed_pupil_illumination_direct_imaging,
            illumination=(offset_m=model.params.offset_m,
                time_scale_m_per_s=model.params.time_scale_m_per_s,
                noise_scale_m=model.params.noise_scale_m)),
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1))
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::IlluminationFrameAcquisitionModel,
    definition::AcquisitionDefinition, path::PreparedPathExecutor)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure_s),
        noise=NoiseNone(), qe=T(model.quantum_efficiency),
        response_model=NullFrameResponse(), T=T, backend=path.key.backend)
    execution = AdaptiveOpticsSim.FrameAcquisitionExecution(detector,
        path.result)
    metadata = (kind=:detector_frame, units=:detected_electrons,
        geometry=path.result.metadata,
        detector=detector_export_metadata(detector))
    products = AdaptiveOpticsSim.AcquisitionProducts(
        execution.observation; metadata)
    return prepare_full_optical_provider(execution, products)
end

function illumination_test_plant(path_models::NamedTuple;
    reverse_declarations::Bool=false, run_seed=0x7400)
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = KolmogorovAtmosphere(telescope; r0=T(0.2), L0=T(25), T=T)
    sources = (
        alpha=Source(band=:custom, wavelength=T(0.75e-6),
            photon_irradiance=T(3), T=T),
        beta=Source(band=:custom, wavelength=T(0.85e-6),
            photon_irradiance=T(4), T=T),
    )
    paths = (
        OpticalPathDefinition(:alpha, sources.alpha, path_models.alpha),
        OpticalPathDefinition(:beta, sources.beta, path_models.beta),
    )
    acquisitions = (
        AcquisitionDefinition(:alpha_frame, :alpha,
            IlluminationFrameAcquisitionModel(T(0.25), T(0.5))),
        AcquisitionDefinition(:beta_frame, :beta,
            IlluminationFrameAcquisitionModel(T(0.5), T(0.75))),
    )
    ordered_paths = reverse_declarations ? reverse(paths) : paths
    ordered_acquisitions = reverse_declarations ?
        reverse(acquisitions) : acquisitions
    definition = PlantDefinition(; telescope, atmosphere,
        paths=ordered_paths, acquisitions=ordered_acquisitions)
    return prepare_plant(definition; run_seed)
end

function illumination_materialization(path::PreparedPathExecutor)
    return AdaptiveOpticsSim.path_materialization(path)
end

function illumination_test_payloads(telescope, source)
    pupil = PupilFunction(telescope)
    field = ElectricField(pupil, source)
    intensity = illumination_detector_input(telescope, source)
    return (; pupil, field, intensity)
end

function illumination_evaluation_allocations(entry, model_time, rng)
    evaluate_illumination!(entry, model_time, rng)
    return @allocated evaluate_illumination!(entry, model_time, rng)
end

function captured_illumination_preparation_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_illumination_preparation_error(f, component::Symbol,
    reason::Symbol)
    error = captured_illumination_preparation_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function illumination_contract_intensity_map(values::AbstractArray;
    kind=FocalPlane(), coordinate_domain=AngularCoordinates(),
    spectral=MonochromaticChannel(convert(eltype(values), 1.0e-6)),
    normalization=PhotonRateNormalization(),
    spatial_measure=CellIntegratedMeasure(),
    coherence=IncoherentIntensityAddition())
    sampling = (one(eltype(values)), one(eltype(values)))
    metadata = OpticalPlaneMetadata(kind, values;
        coordinate_domain, sampling, spectral, normalization,
        spatial_measure, coherence)
    return IntensityMap(metadata, values)
end

function illumination_selection_execution_allocations(selection,
    epoch::AtmosphereEpoch)
    execute_acquisition_selection!(selection, epoch)
    return @allocated execute_acquisition_selection!(selection, epoch)
end

@testset "Typed calibration-illumination entry seam" begin
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T)
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(3), T=T)
    payloads = illumination_test_payloads(telescope, source)
    rng = Xoshiro(0x7410)

    boundary_payloads = (
        (PupilFunctionIlluminationEntry(), payloads.pupil),
        (ElectricFieldIlluminationEntry(), payloads.field),
        (IntensityMapIlluminationEntry(), payloads.intensity),
        (ExternalOpticsResultIlluminationEntry(), payloads.field),
        (ExternalOpticsResultIlluminationEntry(), payloads.intensity),
        (DetectorInputIlluminationEntry(), payloads.intensity),
    )
    for (boundary, payload) in boundary_payloads
        entry = prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), payload,
            boundary; visibility=(consumer=:test_path,))
        @test illumination_entry_boundary(entry) === boundary
        @test illumination_destination(entry) === payload
        @test illumination_combination(entry) isa SingleIllumination
        @test illumination_visibility(entry) == (consumer=:test_path,)
        @test @inferred(evaluate_illumination!(entry, T(0.125), rng)) ===
            payload
    end

    coherent_entry = prepare_illumination_entry(
        TestIlluminationDefinition(CoherentFieldCombination()),
        payloads.field, ElectricFieldIlluminationEntry();
        visibility=(consumer=:coherent_path,))
    incoherent_entry = prepare_illumination_entry(
        TestIlluminationDefinition(IncoherentIntensityAddition()),
        payloads.intensity, IntensityMapIlluminationEntry();
        visibility=(consumer=:incoherent_path,))
    exclusive_entry = prepare_illumination_entry(
        TestIlluminationDefinition(ExclusiveIlluminationSelection()),
        payloads.pupil, PupilFunctionIlluminationEntry();
        visibility=(consumer=:selected_path,))
    @test illumination_combination(coherent_entry) isa
        CoherentFieldCombination
    @test illumination_combination(incoherent_entry) isa
        IncoherentIntensityAddition
    @test illumination_combination(exclusive_entry) isa
        ExclusiveIlluminationSelection

    visibility = Dict(:downstream_path => :science,
        :surface => :detector_input)
    visibility_entry = prepare_illumination_entry(
        TestIlluminationDefinition(SingleIllumination()),
        payloads.intensity, DetectorInputIlluminationEntry(); visibility)
    visibility[:surface] = :mutated
    exposed_visibility = illumination_visibility(visibility_entry)
    exposed_visibility[:surface] = :also_mutated
    @test illumination_visibility(visibility_entry) == Dict(
        :downstream_path => :science, :surface => :detector_input)
    @test visibility_entry.visibility == Dict(
        :downstream_path => :science, :surface => :detector_input)

    contract_snapshot = visibility_entry.contract
    @test contract_snapshot.metadata == payloads.intensity.metadata
    @test AdaptiveOpticsSim.validate_illumination_entry_binding(
        visibility_entry, payloads.intensity) === visibility_entry
    @test !applicable(AdaptiveOpticsSim.PreparedIlluminationEntry,
        illumination_entry_boundary(visibility_entry),
        illumination_evaluator(visibility_entry),
        illumination_destination(visibility_entry),
        illumination_combination(visibility_entry),
        illumination_visibility(visibility_entry),
        contract_snapshot)
    @test !applicable(
        AdaptiveOpticsSim.validate_illumination_evaluator_binding,
        illumination_evaluator(visibility_entry),
        illumination_destination(visibility_entry),
        illumination_entry_boundary(visibility_entry),
        contract_snapshot)

    native_models = (
        alpha=NativeDetectorIlluminationPathModel(T(8)),
        beta=NativeDetectorIlluminationPathModel(T(12)),
    )
    native_plant = illumination_test_plant(native_models;
        run_seed=0x7420)
    native_selection = prepare_acquisition_selection(native_plant,
        (:beta_frame, :alpha_frame))
    @test @inferred(execute_acquisition_selection_at!(native_selection,
        T(0.01))) === native_selection
    alpha_path = prepared_path(native_plant, :alpha)
    beta_path = prepared_path(native_plant, :beta)
    @test all(==(T(8)), path_result(alpha_path).values)
    @test all(==(T(12)), path_result(beta_path).values)
    @test all(==(T(1)), acquisition_observation(
        prepared_acquisition(native_plant, :alpha_frame)))
    @test all(==(T(4.5)), acquisition_observation(
        prepared_acquisition(native_plant, :beta_frame)))
    @test illumination_entry_boundary(
        illumination_materialization(alpha_path)) isa
        DetectorInputIlluminationEntry
    @test path_input(alpha_path) === path_result(alpha_path)

    timed_models = (
        alpha=TimedPupilIlluminationPathModel(TimedPupilIllumination(
            T(1e-9), T(2e-9), T(0.25e-9))),
        beta=TimedPupilIlluminationPathModel(TimedPupilIllumination(
            T(-1e-9), T(3e-9), T(0.5e-9))),
    )
    timed_plant = illumination_test_plant(timed_models;
        run_seed=0x7430)
    reordered_plant = illumination_test_plant(timed_models;
        reverse_declarations=true, run_seed=0x7430)
    timed_selection = prepare_acquisition_selection(timed_plant,
        (:alpha_frame, :beta_frame))
    reordered_selection = prepare_acquisition_selection(reordered_plant,
        (:beta_frame, :alpha_frame))
    execute_acquisition_selection_at!(timed_selection, T(0.02))
    execute_acquisition_selection_at!(reordered_selection, T(0.02))
    for id in (:alpha, :beta)
        path = prepared_path(timed_plant, id)
        reordered_path = prepared_path(reordered_plant, id)
        @test path_input(path).opd == path_input(reordered_path).opd
        @test path_result(path).values ≈ path_result(reordered_path).values
        entry = illumination_materialization(path)
        evaluator = illumination_evaluator(entry)
        @test evaluator.state.evaluations == 1
        @test evaluator.state.last_time_s == T(0.02)
        @test illumination_combination(entry) isa
            ExclusiveIlluminationSelection
        @test illumination_visibility(entry)[:downstream_path] === id
    end
    alpha_entry = illumination_materialization(
        prepared_path(timed_plant, :alpha))
    alpha_before = copy(illumination_destination(alpha_entry).opd)
    assert_illumination_preparation_error(
        () -> materialize_path_input!(prepared_path(timed_plant, :alpha),
            current_epoch(plant_atmosphere(timed_plant.definition))),
        :illumination, :rng_owner)
    @test materialize_path_input!(prepared_path(timed_plant, :alpha),
        current_epoch(plant_atmosphere(timed_plant.definition)),
        Xoshiro(0x7431)) === illumination_destination(alpha_entry)
    @test illumination_destination(alpha_entry).opd != alpha_before

    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            MissingCombinationIlluminationDefinition(), payloads.intensity,
            DetectorInputIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :missing_combination)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(CoherentFieldCombination()),
            payloads.intensity, IntensityMapIlluminationEntry();
            visibility=(consumer=:x,)),
        :illumination, :combination)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(IncoherentIntensityAddition()),
            payloads.field, ElectricFieldIlluminationEntry();
            visibility=(consumer=:x,)),
        :illumination, :combination)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(:undeclared), payloads.intensity,
            DetectorInputIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :combination)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), payloads.field,
            PupilFunctionIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :entry_payload)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), payloads.pupil,
            DetectorInputIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :entry_payload)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()),
            payloads.intensity, DetectorInputIlluminationEntry();
            visibility=nothing),
        :illumination, :visibility)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(MutableIlluminationDefinition(1),
            payloads.intensity, DetectorInputIlluminationEntry();
            visibility=(consumer=:x,)),
        :illumination, :mutable_definition)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            MutableEvaluatorIlluminationDefinition(), payloads.intensity,
            DetectorInputIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :mutable_evaluator)

    invalid_result_entry = prepare_illumination_entry(
        InvalidResultIlluminationDefinition(), payloads.intensity,
        DetectorInputIlluminationEntry(); visibility=(consumer=:x,))
    assert_illumination_preparation_error(
        () -> evaluate_illumination!(invalid_result_entry, zero(T), rng),
        :illumination, :evaluator_result)
    assert_illumination_preparation_error(
        () -> evaluate_illumination!(visibility_entry, T(Inf), rng),
        :illumination, :model_time)

    invalid_plane = illumination_contract_intensity_map(zeros(T, 2, 2);
        kind=PupilPlane(), coordinate_domain=MetricCoordinates())
    invalid_spectral = illumination_contract_intensity_map(zeros(T, 2, 2);
        spectral=UnspecifiedSpectralCoordinate())
    invalid_measure = illumination_contract_intensity_map(zeros(T, 2, 2);
        spatial_measure=PointSampledMeasure())
    invalid_normalization = illumination_contract_intensity_map(
        zeros(T, 2, 2); normalization=UnspecifiedNormalization())
    invalid_coherence = illumination_contract_intensity_map(zeros(T, 2, 2);
        coherence=CoherentFieldCombination())
    for (payload, reason) in (
        (invalid_plane, :entry_plane),
        (invalid_spectral, :spectral_sampling),
        (invalid_measure, :radiometry),
        (invalid_normalization, :radiometry),
        (invalid_coherence, :combination),
    )
    assert_illumination_preparation_error(
            () -> prepare_illumination_entry(
                TestIlluminationDefinition(SingleIllumination()), payload,
                DetectorInputIlluminationEntry(); visibility=(consumer=:x,)),
            :illumination, reason)
    end
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), zeros(T, 2, 2),
            IntensityMapIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :entry_payload)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()),
            payloads.intensity, :unsupported_boundary;
            visibility=(consumer=:x,)),
        :illumination, :entry_payload)
    malformed_support = falses(1, 1)
    malformed_pupil = PupilFunction{
        typeof(payloads.pupil.metadata),typeof(malformed_support),
        typeof(payloads.pupil.amplitude),typeof(payloads.pupil.opd),
        typeof(backend(payloads.pupil)),
    }(payloads.pupil.metadata, malformed_support,
        payloads.pupil.amplitude, payloads.pupil.opd,
        AdaptiveOpticsSim.aperture_revision(payloads.pupil))
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), malformed_pupil,
            PupilFunctionIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :shape)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()),
            invalid_coherence, IntensityMapIlluminationEntry();
            visibility=(consumer=:x,)),
        :illumination, :combination)
    incoherent_field_values = copy(payloads.field.values)
    incoherent_field_metadata = OpticalPlaneMetadata(PupilPlane(),
        incoherent_field_values;
        coordinate_domain=MetricCoordinates(),
        sampling=payloads.field.metadata.sampling,
        origin=payloads.field.metadata.origin,
        orientation=payloads.field.metadata.orientation,
        spectral=payloads.field.metadata.spectral,
        normalization=payloads.field.metadata.normalization,
        spatial_measure=payloads.field.metadata.spatial_measure,
        coherence=IncoherentIntensityAddition())
    incoherent_field = ElectricField(incoherent_field_metadata,
        incoherent_field_values)
    assert_illumination_preparation_error(
        () -> prepare_illumination_entry(
            TestIlluminationDefinition(SingleIllumination()), incoherent_field,
            ElectricFieldIlluminationEntry(); visibility=(consumer=:x,)),
        :illumination, :combination)

    normalized_input = illumination_contract_intensity_map(zeros(T, 2, 2);
        normalization=DimensionlessNormalization())
    @test prepare_illumination_entry(
        TestIlluminationDefinition(SingleIllumination()), normalized_input,
        DetectorInputIlluminationEntry(); visibility=(consumer=:x,)) isa
        AdaptiveOpticsSim.PreparedIlluminationEntry

    uniform_entry = prepare_illumination_entry(
        UniformIntensityIllumination(T(7);
            combination=SingleIllumination()), payloads.intensity,
        DetectorInputIlluminationEntry(); visibility=(consumer=:x,))
    @test evaluate_illumination!(uniform_entry, zero(T), rng) ===
        payloads.intensity
    @test all(==(T(7)), payloads.intensity.values)
    @test_throws InvalidConfiguration UniformIntensityIllumination(-one(T);
        combination=SingleIllumination())
    @test_throws InvalidConfiguration UniformIntensityIllumination(T(Inf);
        combination=SingleIllumination())

    if coverage_instrumented()
        @test_skip "illumination allocation assertions are disabled under coverage instrumentation"
    else
        @test illumination_evaluation_allocations(uniform_entry, zero(T),
            rng) == 0
        @test illumination_evaluation_allocations(alpha_entry, T(0.02),
            rng) == 0
        @test illumination_selection_execution_allocations(native_selection,
            current_epoch(plant_atmosphere(native_plant.definition))) == 0
    end
end
