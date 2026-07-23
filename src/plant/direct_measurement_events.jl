const _DIRECT_MEASUREMENT_EVENT_COMPONENT = :direct_measurement_acquisition

@noinline function _direct_measurement_event_error(reason::Symbol,
    message::AbstractString)
    throw(PlantScheduleError(_DIRECT_MEASUREMENT_EVENT_COMPONENT, reason,
        String(message)))
end

"""
    DirectMeasurementAcquisitionDefinition(exposure_duration;
        readout_duration=PlantDuration(0),
        readiness_delay=PlantDuration(0))

Virtual-time lifecycle for an intentional direct measurement. The measurement
is the time average of the most recently sampled instantaneous product over
the half-open exposure interval. `readout_duration` and `readiness_delay`
retain the same completion and retrigger semantics as a frame acquisition,
without introducing a fictitious detector or photon-rate product.
"""
struct DirectMeasurementAcquisitionDefinition <:
    AbstractAcquisitionLifecycleDefinition
    exposure_duration::PlantDuration
    readout_duration::PlantDuration
    readiness_delay::PlantDuration

    function DirectMeasurementAcquisitionDefinition(
        exposure_duration::PlantDuration,
        readout_duration::PlantDuration,
        readiness_delay::PlantDuration)
        iszero(exposure_duration) && _direct_measurement_event_error(
            :invalid_definition,
            "direct-measurement exposure duration must be > 0 ns")
        return new(exposure_duration, readout_duration, readiness_delay)
    end
end

DirectMeasurementAcquisitionDefinition(
    exposure_duration::PlantDuration;
    readout_duration::PlantDuration=zero(PlantDuration),
    readiness_delay::PlantDuration=zero(PlantDuration),
) = DirectMeasurementAcquisitionDefinition(exposure_duration,
    readout_duration, readiness_delay)

mutable struct _DirectMeasurementAcquisitionBinding end

"""
Prepared direct-measurement lifecycle and caller-owned integration storage.

Construct with `prepare_direct_measurement_acquisition`. The instantaneous
sample, exposure integral, and published measurement are distinct buffers.
"""
struct PreparedDirectMeasurementAcquisition{
    M<:WFSMeasurement,
    T<:AbstractFloat,
    V<:AbstractVector{T},
} <: AbstractPreparedAcquisitionLifecycle
    binding::_DirectMeasurementAcquisitionBinding
    measurement::M
    definition::DirectMeasurementAcquisitionDefinition
    exposure_seconds::T
    instantaneous_sample::V
    integrated_sample::V
end

@enum DirectMeasurementAcquisitionStatus::UInt8 begin
    DirectMeasurementReady = 0x01
    DirectMeasurementExposureActive = 0x02
    DirectMeasurementReadoutPending = 0x03
    DirectMeasurementReadoutComplete = 0x04
end

"""Separately owned, single-writer direct-measurement lifecycle state."""
mutable struct DirectMeasurementAcquisitionState <:
    AbstractAcquisitionLifecycleState
    binding::_DirectMeasurementAcquisitionBinding
    status::DirectMeasurementAcquisitionStatus
    sequence::UInt64
    exposure_start::PlantTimestamp
    exposure_close::PlantTimestamp
    integrated_through::PlantTimestamp
    readout_complete::PlantTimestamp
    readiness::PlantTimestamp
end

function DirectMeasurementAcquisitionState(
    prepared::PreparedDirectMeasurementAcquisition)
    origin = zero(PlantTimestamp)
    return DirectMeasurementAcquisitionState(prepared.binding,
        DirectMeasurementReady, UInt64(0), origin, origin, origin, origin,
        origin)
end

@inline direct_measurement_acquisition_status(
    state::DirectMeasurementAcquisitionState) = state.status
@inline direct_measurement_acquisition_sequence(
    state::DirectMeasurementAcquisitionState) = state.sequence
@inline exposure_start_timestamp(
    state::DirectMeasurementAcquisitionState) = state.exposure_start
@inline exposure_close_timestamp(
    state::DirectMeasurementAcquisitionState) = state.exposure_close
@inline integrated_through_timestamp(
    state::DirectMeasurementAcquisitionState) = state.integrated_through
@inline readout_complete_timestamp(
    state::DirectMeasurementAcquisitionState) = state.readout_complete
@inline acquisition_readiness_timestamp(
    state::DirectMeasurementAcquisitionState) = state.readiness

"""
    prepare_direct_measurement_acquisition(measurement, definition)

Bind one vector-valued floating-point `WFSMeasurement` to a fixed
direct-measurement lifecycle. Preparation allocates the instantaneous sample
and exposure-integral buffers; repeated lifecycle execution reuses them.
"""
@inline _require_direct_measurement_storage(
    storage::AbstractVector) = storage

function _require_direct_measurement_storage(storage)
    _direct_measurement_event_error(:shape,
        "direct measurement storage must be an AbstractVector")
end

function prepare_direct_measurement_acquisition(
    measurement::WFSMeasurement,
    definition::DirectMeasurementAcquisitionDefinition)
    validate_wfs_measurement(measurement)
    storage = _require_direct_measurement_storage(
        measurement_storage(measurement))
    T = eltype(storage)
    T <: AbstractFloat || _direct_measurement_event_error(:numeric_type,
        "direct measurement storage must use an AbstractFloat element type")
    instantaneous = similar(storage)
    integrated = similar(storage)
    fill!(instantaneous, zero(T))
    fill!(integrated, zero(T))
    exposure_seconds = plant_duration_seconds(
        definition.exposure_duration, T)
    return PreparedDirectMeasurementAcquisition(
        _DirectMeasurementAcquisitionBinding(), measurement, definition,
        exposure_seconds, instantaneous, integrated)
end

@inline function _require_direct_measurement_binding(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState)
    state.binding === prepared.binding ||
        _direct_measurement_event_error(:foreign_state,
            "direct-measurement state belongs to another prepared lifecycle")
    return nothing
end

function begin_exposure!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp)
    _require_direct_measurement_binding(prepared, state)
    state.status == DirectMeasurementReady ||
        _direct_measurement_event_error(:acquisition_busy,
            "direct measurement is not ready for a new exposure")
    state.sequence == typemax(UInt64) &&
        _direct_measurement_event_error(:sequence_overflow,
            "direct-measurement sequence overflow")
    close_timestamp = timestamp + prepared.definition.exposure_duration
    readout_timestamp = close_timestamp +
        prepared.definition.readout_duration
    readiness_timestamp = readout_timestamp +
        prepared.definition.readiness_delay
    fill!(prepared.integrated_sample,
        zero(eltype(prepared.integrated_sample)))
    state.sequence += UInt64(1)
    state.exposure_start = timestamp
    state.exposure_close = close_timestamp
    state.integrated_through = timestamp
    state.readout_complete = readout_timestamp
    state.readiness = readiness_timestamp
    state.status = DirectMeasurementExposureActive
    return state
end

"""
    accumulate_direct_measurement_interval!(prepared, state, start, stop)

Integrate the currently held instantaneous sample over `[start, stop)`.
Sampling changes are composed by first integrating through the sample
timestamp, then replacing the held sample.
"""
function accumulate_direct_measurement_interval!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    start::PlantTimestamp,
    stop::PlantTimestamp)
    _require_direct_measurement_binding(prepared, state)
    state.status == DirectMeasurementExposureActive ||
        _direct_measurement_event_error(:exposure_inactive,
            "direct-measurement integration requires an active exposure")
    start == state.integrated_through ||
        _direct_measurement_event_error(:noncontiguous_interval,
            "direct-measurement integration must continue from prior progress")
    start <= stop <= state.exposure_close ||
        _direct_measurement_event_error(:invalid_interval,
            "direct-measurement interval lies outside the active exposure")
    start == stop && return state
    T = eltype(prepared.integrated_sample)
    duration = plant_duration_seconds(stop - start, T)
    @. prepared.integrated_sample +=
        prepared.instantaneous_sample * duration
    state.integrated_through = stop
    return state
end

function close_exposure!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp)
    _require_direct_measurement_binding(prepared, state)
    state.status == DirectMeasurementExposureActive ||
        _direct_measurement_event_error(:exposure_inactive,
            "cannot close an inactive direct-measurement exposure")
    timestamp == state.exposure_close ||
        _direct_measurement_event_error(:close_timestamp,
            "direct-measurement close timestamp does not match its schedule")
    state.integrated_through == timestamp ||
        _direct_measurement_event_error(:incomplete_integration,
            "direct measurement was not integrated through exposure close")
    state.status = DirectMeasurementReadoutPending
    return state
end

function complete_readout!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp,
    ::AbstractRNG)
    _require_direct_measurement_binding(prepared, state)
    state.status == DirectMeasurementReadoutPending ||
        _direct_measurement_event_error(:readout_state,
            "direct measurement is not pending readout")
    timestamp == state.readout_complete ||
        _direct_measurement_event_error(:readout_timestamp,
            "direct-measurement readout timestamp does not match its schedule")
    output = measurement_storage(prepared.measurement)
    @. output = prepared.integrated_sample / prepared.exposure_seconds
    state.status = DirectMeasurementReadoutComplete
    return prepared.measurement
end

function mark_acquisition_ready!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp)
    _require_direct_measurement_binding(prepared, state)
    state.status == DirectMeasurementReadoutComplete ||
        _direct_measurement_event_error(:readiness_state,
            "direct measurement has not completed readout")
    timestamp == state.readiness ||
        _direct_measurement_event_error(:readiness_timestamp,
            "direct-measurement readiness timestamp does not match its schedule")
    state.status = DirectMeasurementReady
    return state
end
