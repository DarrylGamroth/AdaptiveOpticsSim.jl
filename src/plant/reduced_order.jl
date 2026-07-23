#
# Deterministic command-responsive reduced-order AO acquisition provider.
#
# This is deliberately a bounded linear direct-measurement model. It preserves
# plant-time command causality and acquisition timing without making
# diffraction, detector, or science-performance claims.
#

const _REDUCED_ORDER_COMPONENT = :reduced_order

@noinline function _reduced_order_error(reason::Symbol,
    message::AbstractString)
    throw(PlantPreparationError(_REDUCED_ORDER_COMPONENT, reason,
        String(message)))
end

function _reduced_order_float_vector(values::AbstractVector,
    label::AbstractString)
    isempty(values) && _reduced_order_error(:shape,
        "$label must not be empty")
    source_type = eltype(values)
    source_type <: Real || _reduced_order_error(:numeric_type,
        "$label must contain real values")
    E = float(source_type)
    isconcretetype(E) && E <: AbstractFloat ||
        _reduced_order_error(:numeric_type,
            "$label must have a concrete real numeric element type")
    copied = Vector{E}(values)
    all(isfinite, copied) || _reduced_order_error(:nonfinite,
        "$label must contain only finite values")
    return copied
end

function _reduced_order_float_vector(values, label::AbstractString)
    _reduced_order_error(:shape,
        "$label must be an AbstractVector; got $(typeof(values))")
end

"""
    HarmonicDisturbanceModel(amplitudes, frequencies_hz;
        offsets=zeros(...), phases_rad=zeros(...))

Deterministic time-correlated modal disturbance

`offset + amplitude * sin(2π * frequency_hz * t + phase_rad)`

evaluated at explicit plant time. Arrays are copied into immutable cold
configuration ownership. Repeated preparation with the same values produces
an exactly replayable disturbance sequence at the same timestamps.
"""
struct HarmonicDisturbanceModel{T<:AbstractFloat,V<:Vector{T}}
    offsets::V
    amplitudes::V
    frequencies_hz::V
    phases_rad::V
end

function HarmonicDisturbanceModel(amplitudes::AbstractVector,
    frequencies_hz::AbstractVector; offsets=nothing, phases_rad=nothing)
    amplitude_values = _reduced_order_float_vector(amplitudes,
        "harmonic disturbance amplitudes")
    T = eltype(amplitude_values)
    frequencies = Vector{T}(frequencies_hz)
    n_modes = length(amplitude_values)
    length(frequencies) == n_modes || _reduced_order_error(:shape,
        "harmonic disturbance frequencies must match the amplitude count")
    all(isfinite, frequencies) && all(value -> value >= zero(T),
        frequencies) || _reduced_order_error(:frequency,
        "harmonic disturbance frequencies must be finite and nonnegative")
    offset_values = offsets === nothing ? zeros(T, n_modes) :
        Vector{T}(offsets)
    phase_values = phases_rad === nothing ? zeros(T, n_modes) :
        Vector{T}(phases_rad)
    length(offset_values) == n_modes || _reduced_order_error(:shape,
        "harmonic disturbance offsets must match the amplitude count")
    length(phase_values) == n_modes || _reduced_order_error(:shape,
        "harmonic disturbance phases must match the amplitude count")
    all(isfinite, offset_values) || _reduced_order_error(:nonfinite,
        "harmonic disturbance offsets must be finite")
    all(isfinite, phase_values) || _reduced_order_error(:nonfinite,
        "harmonic disturbance phases must be finite")
    return HarmonicDisturbanceModel(offset_values, amplitude_values,
        frequencies, phase_values)
end

function HarmonicDisturbanceModel(amplitudes, frequencies_hz; kwargs...)
    _reduced_order_error(:shape,
        "harmonic disturbance amplitudes and frequencies must be AbstractVector values")
end

"""
    ReducedOrderCommandResponse(endpoint, operator;
        units, sign_convention, basis, basis_revision)

Calibrated mapping from one effective command endpoint into reduced residual
coordinates. A vector operator binds a scalar command; a matrix binds a
fixed-shape array command after column-major vectorization. The operator
contains the physical response sign. Expected command semantics are checked
against the prepared endpoint before an event loop is constructed.
"""
struct ReducedOrderCommandResponse{
    O<:Union{Vector{<:AbstractFloat},Matrix{<:AbstractFloat}},
}
    endpoint::CommandEndpointID
    operator::O
    units::CommandUnit
    sign_convention::CommandSignConvention
    basis::CommandBasis
    basis_revision::CommandBasisRevision
end

function _reduced_order_operator(operator::Union{
    AbstractVector,AbstractMatrix})
    source_type = eltype(operator)
    source_type <: Real || _reduced_order_error(:numeric_type,
        "command-response operator must contain real values")
    E = float(source_type)
    isconcretetype(E) && E <: AbstractFloat ||
        _reduced_order_error(:numeric_type,
            "command-response operator must have a concrete real numeric element type")
    copied = ndims(operator) == 1 ? Vector{E}(operator) :
        Matrix{E}(operator)
    isempty(copied) && _reduced_order_error(:shape,
        "command-response operator must not be empty")
    all(isfinite, copied) || _reduced_order_error(:nonfinite,
        "command-response operator must contain only finite values")
    return copied
end

function ReducedOrderCommandResponse(endpoint,
    operator::Union{AbstractVector,AbstractMatrix}; units,
    sign_convention, basis::CommandBasis, basis_revision)
    return ReducedOrderCommandResponse(
        _as_command_endpoint_id(endpoint),
        _reduced_order_operator(operator),
        _as_command_unit(units),
        _as_command_sign_convention(sign_convention),
        basis,
        _as_command_basis_revision(basis_revision),
    )
end

function ReducedOrderCommandResponse(endpoint, operator; kwargs...)
    _reduced_order_error(:shape,
        "command-response operator must be an AbstractVector or AbstractMatrix; got $(typeof(operator))")
end

@inline _require_reduced_order_response(
    response::ReducedOrderCommandResponse) = response

function _require_reduced_order_response(response)
    _reduced_order_error(:response_type,
        "reduced-order responses must contain ReducedOrderCommandResponse values; got $(typeof(response))")
end

function _reduced_order_response_tuple(responses::Tuple)
    foreach(_require_reduced_order_response, responses)
    return responses
end

function _reduced_order_response_tuple(responses::NamedTuple)
    return _reduced_order_response_tuple(values(responses))
end

function _reduced_order_response_tuple(responses::AbstractVector)
    return _reduced_order_response_tuple(Tuple(responses))
end

function _reduced_order_response_tuple(responses)
    _reduced_order_error(:response_type,
        "reduced-order responses must be a Tuple, NamedTuple, or AbstractVector")
end

function _require_unique_reduced_order_responses(responses::Tuple)
    isempty(responses) && _reduced_order_error(:empty_responses,
        "a command-responsive reduced-order model requires at least one command response")
    @inbounds for right in 2:length(responses)
        endpoint = responses[right].endpoint
        for left in 1:(right - 1)
            responses[left].endpoint == endpoint &&
                _reduced_order_error(:duplicate_endpoint,
                    "reduced-order endpoint $endpoint has more than one response operator")
        end
    end
    return responses
end

@inline function _require_reduced_order_omission(value::Symbol)
    isempty(String(value)) && _reduced_order_error(:omitted_effects,
        "reduced-order omitted effects must be nonempty Symbols")
    return value
end

function _require_reduced_order_omission(value)
    _reduced_order_error(:omitted_effects,
        "reduced-order omitted effects must be Symbols; got $(typeof(value))")
end

function _reduced_order_omissions(values::Tuple)
    isempty(values) && _reduced_order_error(:omitted_effects,
        "reduced-order omitted effects must be declared")
    @inbounds for (index, value) in pairs(values)
        _require_reduced_order_omission(value)
        value in values[1:(index - 1)] &&
            _reduced_order_error(:omitted_effects,
                "reduced-order omitted effects must be unique")
    end
    return values
end

_reduced_order_omissions(values::AbstractVector) =
    _reduced_order_omissions(Tuple(values))

function _reduced_order_omissions(values)
    _reduced_order_error(:omitted_effects,
        "reduced-order omitted effects must be a Tuple or AbstractVector")
end

function _checked_reduced_order_revision(value::Integer)
    value > 0 || _reduced_order_error(:calibration_revision,
        "reduced-order calibration revision must be positive")
    value <= typemax(UInt32) || _reduced_order_error(
        :calibration_revision,
        "reduced-order calibration revision exceeds UInt32 range")
    return UInt32(value)
end

_checked_reduced_order_revision(::Bool) =
    _reduced_order_error(:calibration_revision,
        "reduced-order calibration revision must be an integer, not Bool")

function _checked_reduced_order_revision(value)
    _reduced_order_error(:calibration_revision,
        "reduced-order calibration revision must be an Integer; got $(typeof(value))")
end

"""
    LinearReducedOrderAcquisitionModel(disturbance, path_projection,
        sensor_operator, command_responses; ...)

Cold linear reduced-order AO acquisition model. `path_projection` maps shared
disturbance modes into path residual coordinates, each calibrated command
response adds in those coordinates, and `sensor_operator` maps the residual to
one direct WFS measurement vector.

The operating envelope and omitted effects are mandatory because this model
does not establish full optical or detector fidelity.
"""
struct LinearReducedOrderAcquisitionModel{
    T<:AbstractFloat,
    D<:HarmonicDisturbanceModel{T},
    P<:Matrix{T},
    S<:Matrix{T},
    R<:Tuple,
    K,U,Q,V,E,
}
    disturbance::D
    path_projection::P
    sensor_operator::S
    command_responses::R
    measurement_kind::K
    measurement_units::U
    residual_kind::Q
    residual_units::V
    calibration_revision::UInt32
    operating_envelope::E
    omitted_effects::Tuple
end

function _convert_reduced_order_response(
    response::ReducedOrderCommandResponse, ::Type{T}) where {
    T<:AbstractFloat,
}
    operator = ndims(response.operator) == 1 ?
        Vector{T}(response.operator) : Matrix{T}(response.operator)
    return ReducedOrderCommandResponse(response.endpoint, operator;
        units=response.units, sign_convention=response.sign_convention,
        basis=response.basis, basis_revision=response.basis_revision)
end

function LinearReducedOrderAcquisitionModel(
    disturbance::HarmonicDisturbanceModel{T},
    path_projection::AbstractMatrix,
    sensor_operator::AbstractMatrix,
    command_responses;
    measurement_kind,
    measurement_units,
    residual_kind,
    residual_units,
    calibration_revision,
    operating_envelope,
    omitted_effects,
) where {T<:AbstractFloat}
    projection = Matrix{T}(path_projection)
    sensor = Matrix{T}(sensor_operator)
    isempty(projection) && _reduced_order_error(:shape,
        "reduced-order path projection must not be empty")
    isempty(sensor) && _reduced_order_error(:shape,
        "reduced-order sensor operator must not be empty")
    all(isfinite, projection) || _reduced_order_error(:nonfinite,
        "reduced-order path projection must contain only finite values")
    all(isfinite, sensor) || _reduced_order_error(:nonfinite,
        "reduced-order sensor operator must contain only finite values")
    size(projection, 2) == length(disturbance.amplitudes) ||
        _reduced_order_error(:disturbance_dimensions,
            "path-projection columns must match disturbance mode count")
    size(sensor, 2) == size(projection, 1) ||
        _reduced_order_error(:sensor_dimensions,
            "sensor-operator columns must match residual coordinate count")
    responses = _require_unique_reduced_order_responses(
        _reduced_order_response_tuple(command_responses))
    converted = map(response ->
        _convert_reduced_order_response(response, T), responses)
    for response in converted
        size(response.operator, 1) == size(projection, 1) ||
            _reduced_order_error(:response_dimensions,
                "every command-response row count must match residual coordinates")
    end
    measurement_kind === nothing && _reduced_order_error(:measurement_kind,
        "reduced-order measurement kind must be declared")
    measurement_units === nothing && _reduced_order_error(
        :measurement_units,
        "reduced-order measurement units must be declared")
    residual_kind === nothing && _reduced_order_error(:residual_kind,
        "reduced-order residual kind must be declared")
    residual_units === nothing && _reduced_order_error(:residual_units,
        "reduced-order residual units must be declared")
    operating_envelope === nothing && _reduced_order_error(
        :operating_envelope,
        "reduced-order operating envelope must be declared")
    omissions = _reduced_order_omissions(omitted_effects)
    return LinearReducedOrderAcquisitionModel(
        disturbance, projection, sensor, converted,
        deepcopy(measurement_kind), deepcopy(measurement_units),
        deepcopy(residual_kind), deepcopy(residual_units),
        _checked_reduced_order_revision(calibration_revision),
        deepcopy(operating_envelope), omissions)
end

plant_model_definition_style(
    ::Type{<:LinearReducedOrderAcquisitionModel}) =
    ColdPlantModelDefinition()

struct PreparedHarmonicDisturbance{
    T<:AbstractFloat,
    V<:AbstractVector{T},
}
    offsets::V
    amplitudes::V
    frequencies_hz::V
    phases_rad::V
end

mutable struct HarmonicDisturbanceState{
    T<:AbstractFloat,
    V<:AbstractVector{T},
}
    modes::V
    timestamp::PlantTimestamp
    initialized::Bool
end

struct PreparedReducedOrderCommandResponse{
    O<:Union{AbstractVector,AbstractMatrix},
}
    endpoint::CommandEndpointID
    operator::O
    units::CommandUnit
    sign_convention::CommandSignConvention
    basis::CommandBasis
    basis_revision::CommandBasisRevision
end

"""
Prepared single-writer linear reduced-order provider.

The mutable disturbance and residual values are separate from immutable
operators and metadata. Event-loop preparation resolves command identities to
fixed endpoint slots before any sample is evaluated.
"""
struct PreparedLinearReducedOrderProvider{
    HD<:PreparedHarmonicDisturbance,
    S<:HarmonicDisturbanceState,
    P<:AbstractMatrix,
    H<:AbstractMatrix,
    R<:Tuple,
    V<:AbstractVector,
    B<:AbstractArrayBackend,
    PD<:AbstractPlaneDevice,
    E,
}
    disturbance::HD
    state::S
    path_projection::P
    sensor_operator::H
    command_responses::R
    residual::V
    command_workspace::V
    backend::B
    device::PD
    calibration_revision::UInt32
    operating_envelope::E
    omitted_effects::Tuple
end

@inline acquisition_provider_style(
    ::Type{<:PreparedLinearReducedOrderProvider}) =
    CommandResponsiveReducedOrderProviderStyle()
@inline acquisition_provider_payload_work(
    ::Type{<:PreparedLinearReducedOrderProvider}) =
    :linear_reduced_order_direct_measurement

function _reduced_order_backend_copy(values::AbstractArray{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    destination = allocate_array(selector, T, size(values)...)
    copyto!(destination, values)
    return destination
end

function _prepare_harmonic_disturbance(
    model::HarmonicDisturbanceModel{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    offsets = _reduced_order_backend_copy(model.offsets, selector)
    amplitudes = _reduced_order_backend_copy(model.amplitudes, selector)
    frequencies = _reduced_order_backend_copy(model.frequencies_hz, selector)
    phases = _reduced_order_backend_copy(model.phases_rad, selector)
    prepared = PreparedHarmonicDisturbance(offsets, amplitudes, frequencies,
        phases)
    modes = similar(offsets)
    fill!(modes, zero(T))
    return prepared, HarmonicDisturbanceState(modes, zero(PlantTimestamp),
        false)
end

function _prepare_reduced_order_response(
    response::ReducedOrderCommandResponse,
    selector::AbstractArrayBackend)
    operator = _reduced_order_backend_copy(response.operator, selector)
    return PreparedReducedOrderCommandResponse(response.endpoint, operator,
        response.units, response.sign_convention, response.basis,
        response.basis_revision)
end

function _linear_reduced_order_metadata(
    model::LinearReducedOrderAcquisitionModel, measurement_count::Int)
    return (
        kind=:direct_wfs_measurement,
        units=deepcopy(model.measurement_units),
        geometry=(dimensions=(measurement_count,), layout=:dense),
        sampling=(quantity=:exposure_averaged_measurement,
            hold=:zero_order),
        model=(
            family=:linear_reduced_order_ao,
            calibration_revision=model.calibration_revision,
            residual_kind=deepcopy(model.residual_kind),
            residual_units=deepcopy(model.residual_units),
            residual_metric=:root_mean_square,
            operating_envelope=deepcopy(model.operating_envelope),
            omitted_effects=model.omitted_effects,
        ),
        semantics=:complete_acquisition,
    )
end

function prepare_acquisition_provider(
    model::LinearReducedOrderAcquisitionModel{T},
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
) where {T<:AbstractFloat}
    selector = getfield(path.key, :backend)
    projection = _reduced_order_backend_copy(model.path_projection,
        selector)
    sensor = _reduced_order_backend_copy(model.sensor_operator, selector)
    responses = map(response ->
        _prepare_reduced_order_response(response, selector),
        model.command_responses)
    disturbance, state = _prepare_harmonic_disturbance(model.disturbance,
        selector)
    residual = similar(projection, T, size(projection, 1))
    command_workspace = similar(residual)
    fill!(residual, zero(T))
    fill!(command_workspace, zero(T))
    measurement_values = allocate_array(selector, T, size(sensor, 1))
    fill!(measurement_values, zero(T))
    measurement = WFSMeasurement(measurement_values;
        units=deepcopy(model.measurement_units),
        kind=deepcopy(model.measurement_kind))
    products = AcquisitionProducts(nothing, measurement;
        metadata=_linear_reduced_order_metadata(model, length(
            measurement_values)))
    implementation = PreparedLinearReducedOrderProvider(
        disturbance, state, projection, sensor, responses, residual,
        command_workspace, selector, plane_device(measurement_values),
        model.calibration_revision, deepcopy(model.operating_envelope),
        model.omitted_effects)
    return PreparedAcquisitionProvider(implementation, products)
end

function validate_acquisition_provider_binding(
    implementation::PreparedLinearReducedOrderProvider,
    path_result,
    products::AcquisitionProducts)
    products.observation === nothing || _reduced_order_error(
        :product_contract,
        "linear reduced-order provider must not publish a fictitious observation")
    measurement = _require_reduced_order_measurement(products.measurement)
    storage = measurement_storage(measurement)
    length(storage) == size(implementation.sensor_operator, 1) ||
        _reduced_order_error(:sensor_dimensions,
            "reduced-order measurement length changed after preparation")
    typeof(backend(storage)) === typeof(implementation.backend) ||
        _reduced_order_error(:backend,
            "reduced-order measurement backend changed after preparation")
    plane_device(storage) == implementation.device ||
        _reduced_order_error(:device,
            "reduced-order measurement device changed after preparation")
    size(implementation.path_projection, 2) ==
        length(implementation.state.modes) ||
        _reduced_order_error(:disturbance_dimensions,
            "reduced-order disturbance storage changed after preparation")
    return nothing
end

@inline _require_reduced_order_measurement(
    measurement::WFSMeasurement) = measurement

function _require_reduced_order_measurement(measurement)
    _reduced_order_error(:product_contract,
        "linear reduced-order provider requires a WFSMeasurement product; got $(typeof(measurement))")
end

function execute_acquisition_provider!(
    ::AcquisitionProducts, path_result,
    ::PreparedLinearReducedOrderProvider, rngs)
    _reduced_order_error(:requires_event_timeline,
        "linear reduced-order acquisition requires a prepared event loop with effective command state")
end

abstract type _PreparedReducedOrderEventResponse end

struct _PreparedScalarReducedOrderEventResponse{
    R<:PreparedReducedOrderCommandResponse,
} <: _PreparedReducedOrderEventResponse
    response::R
    endpoint_slot::UInt32
end

struct _PreparedArrayReducedOrderEventResponse{
    R<:PreparedReducedOrderCommandResponse,
} <: _PreparedReducedOrderEventResponse
    response::R
    endpoint_slot::UInt32
end

struct PreparedLinearReducedOrderEventProvider{
    P<:PreparedLinearReducedOrderProvider,
    R<:Tuple,
}
    provider::P
    responses::R
end

function _reduced_order_endpoint_slot(endpoints,
    id::CommandEndpointID)
    @inbounds for index in eachindex(endpoints)
        command_endpoint_id(endpoints[index]) == id && return UInt32(index)
    end
    _reduced_order_error(:unknown_command_endpoint,
        "reduced-order response references unknown endpoint $id")
end

function _require_reduced_order_schema_semantics(
    response::PreparedReducedOrderCommandResponse,
    schema::PlantCommandSchema)
    command_units(schema) == response.units ||
        _reduced_order_error(:command_units,
            "reduced-order response units do not match endpoint $(response.endpoint)")
    command_sign_convention(schema) == response.sign_convention ||
        _reduced_order_error(:command_sign_convention,
            "reduced-order response sign convention does not match endpoint $(response.endpoint)")
    command_basis(schema) == response.basis ||
        _reduced_order_error(:command_basis,
            "reduced-order response basis does not match endpoint $(response.endpoint)")
    command_basis_revision(schema) == response.basis_revision ||
        _reduced_order_error(:command_basis_revision,
            "reduced-order response basis revision does not match endpoint $(response.endpoint)")
    return nothing
end

function _prepare_reduced_order_event_response(
    response::PreparedReducedOrderCommandResponse{<:AbstractVector},
    endpoints)
    slot = _reduced_order_endpoint_slot(endpoints, response.endpoint)
    endpoint = @inbounds endpoints[Int(slot)].endpoint
    schema = command_schema(endpoint)
    isempty(command_dimensions(schema)) || _reduced_order_error(
        :command_dimensions,
        "vector command-response operator requires a scalar endpoint")
    command_numeric_type(schema) === eltype(response.operator) ||
        _reduced_order_error(:command_numeric_type,
            "reduced-order response numeric type does not match endpoint $(response.endpoint)")
    _require_reduced_order_schema_semantics(response, schema)
    length(response.operator) > 0 || _reduced_order_error(
        :response_dimensions,
        "scalar command-response operator must not be empty")
    return _PreparedScalarReducedOrderEventResponse(response, slot)
end

function _prepare_reduced_order_event_response(
    response::PreparedReducedOrderCommandResponse{<:AbstractMatrix},
    endpoints)
    slot = _reduced_order_endpoint_slot(endpoints, response.endpoint)
    binding = @inbounds endpoints[Int(slot)]
    endpoint = binding.endpoint
    schema = command_schema(endpoint)
    dimensions = command_dimensions(schema)
    isempty(dimensions) && _reduced_order_error(:command_dimensions,
        "matrix command-response operator requires an array endpoint")
    size(response.operator, 2) == prod(dimensions) ||
        _reduced_order_error(:command_dimensions,
            "command-response columns do not match endpoint $(response.endpoint) dimensions")
    command_numeric_type(schema) === eltype(response.operator) ||
        _reduced_order_error(:command_numeric_type,
            "reduced-order response numeric type does not match endpoint $(response.endpoint)")
    typeof(backend(endpoint)) ===
        typeof(backend(response.operator)) ||
        _reduced_order_error(:command_backend,
            "reduced-order response backend does not match endpoint $(response.endpoint)")
    plane_device(binding.initial_command) ==
        plane_device(response.operator) ||
        _reduced_order_error(:command_device,
            "reduced-order response device does not match endpoint $(response.endpoint)")
    _require_reduced_order_schema_semantics(response, schema)
    return _PreparedArrayReducedOrderEventResponse(response, slot)
end

function prepare_linear_reduced_order_event_provider(
    provider::PreparedLinearReducedOrderProvider, endpoints)
    length(provider.command_responses) == length(endpoints) ||
        _reduced_order_error(:path_visibility,
            "the current all-path visibility contract requires one explicit reduced-order response for every prepared command endpoint")
    responses = map(response ->
        _prepare_reduced_order_event_response(response, endpoints),
        provider.command_responses)
    @inbounds for endpoint in endpoints
        id = command_endpoint_id(endpoint)
        any(response -> response.response.endpoint == id, responses) ||
            _reduced_order_error(:path_visibility,
                "reduced-order provider omits currently visible endpoint $id")
    end
    return PreparedLinearReducedOrderEventProvider(provider, responses)
end

function _evaluate_harmonic_disturbance!(
    prepared::PreparedHarmonicDisturbance{T},
    state::HarmonicDisturbanceState{T},
    timestamp::PlantTimestamp) where {T<:AbstractFloat}
    state.initialized && timestamp < state.timestamp &&
        _reduced_order_error(:time_regression,
            "reduced-order disturbance time regressed")
    time_seconds = plant_time_seconds(timestamp, T)
    @. state.modes = prepared.offsets + prepared.amplitudes *
        sin(T(2π) * prepared.frequencies_hz * time_seconds +
            prepared.phases_rad)
    state.timestamp = timestamp
    state.initialized = true
    return state.modes
end

@inline function _apply_reduced_order_response!(
    residual, workspace,
    response::_PreparedScalarReducedOrderEventResponse,
    command_applications)
    command = effective_command(
        @inbounds command_applications[Int(response.endpoint_slot)])
    @. workspace = response.response.operator * command
    residual .+= workspace
    return residual
end

@inline function _apply_reduced_order_response!(
    residual, workspace,
    response::_PreparedArrayReducedOrderEventResponse,
    command_applications)
    command = effective_command(
        @inbounds command_applications[Int(response.endpoint_slot)])
    mul!(workspace, response.response.operator, vec(command))
    residual .+= workspace
    return residual
end

@inline _apply_reduced_order_responses!(
    residual, workspace, ::Tuple{}, command_applications) = residual

@inline function _apply_reduced_order_responses!(
    residual, workspace, responses::Tuple, command_applications)
    _apply_reduced_order_response!(residual, workspace, first(responses),
        command_applications)
    return _apply_reduced_order_responses!(residual, workspace,
        Base.tail(responses), command_applications)
end

function evaluate_linear_reduced_order_sample!(
    destination::AbstractVector,
    prepared::PreparedLinearReducedOrderEventProvider,
    timestamp::PlantTimestamp,
    command_applications)
    provider = prepared.provider
    length(destination) == size(provider.sensor_operator, 1) ||
        _reduced_order_error(:sensor_dimensions,
            "direct-measurement sample length changed after preparation")
    disturbance = _evaluate_harmonic_disturbance!(provider.disturbance,
        provider.state, timestamp)
    mul!(provider.residual, provider.path_projection, disturbance)
    _apply_reduced_order_responses!(provider.residual,
        provider.command_workspace, prepared.responses,
        command_applications)
    mul!(destination, provider.sensor_operator, provider.residual)
    return destination
end

@inline reduced_order_disturbance(
    provider::PreparedLinearReducedOrderProvider) = provider.state.modes
@inline reduced_order_residual(
    provider::PreparedLinearReducedOrderProvider) = provider.residual
@inline reduced_order_sample_timestamp(
    provider::PreparedLinearReducedOrderProvider) =
    provider.state.initialized ? provider.state.timestamp : nothing

@inline reduced_order_disturbance(provider::PreparedAcquisitionProvider) =
    reduced_order_disturbance(provider.implementation)
@inline reduced_order_residual(provider::PreparedAcquisitionProvider) =
    reduced_order_residual(provider.implementation)
@inline reduced_order_sample_timestamp(
    provider::PreparedAcquisitionProvider) =
    reduced_order_sample_timestamp(provider.implementation)

@inline reduced_order_disturbance(owner::PreparedAcquisitionOwner) =
    reduced_order_disturbance(owner.provider)
@inline reduced_order_residual(owner::PreparedAcquisitionOwner) =
    reduced_order_residual(owner.provider)
@inline reduced_order_sample_timestamp(owner::PreparedAcquisitionOwner) =
    reduced_order_sample_timestamp(owner.provider)

function reduced_order_residual_rms(value)
    residual = reduced_order_residual(value)
    return sqrt(sum(abs2, residual) / length(residual))
end
