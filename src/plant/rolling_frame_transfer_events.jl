abstract type _RollingShutterEventMode end
struct _RollingExposureEventMode <: _RollingShutterEventMode end
struct _GlobalResetEventMode <: _RollingShutterEventMode end

@inline _rolling_shutter_event_mode(::RollingExposure) =
    _RollingExposureEventMode()
@inline _rolling_shutter_event_mode(::GlobalResetExposure) =
    _GlobalResetEventMode()

"""
    RollingShutterAcquisitionDefinition(exposure_duration;
        readiness_delay=PlantDuration(0))

Immutable virtual-time definition for one rolling-shutter frame. The detector's
prepared `RollingShutter` timing model supplies the row-band size, line time,
and rolling-exposure versus global-reset semantics. The definition supplies the
nominal per-band exposure and the delay between complete-frame readout and the
next accepted frame start.
"""
struct RollingShutterAcquisitionDefinition <:
    AbstractDetectorAcquisitionEventDefinition
    exposure_duration::PlantDuration
    readiness_delay::PlantDuration

    function RollingShutterAcquisitionDefinition(
        exposure_duration::PlantDuration, readiness_delay::PlantDuration)
        iszero(exposure_duration) && _detector_acquisition_event_error(
            :invalid_definition,
            "rolling-shutter exposure duration must be > 0 ns")
        return new(exposure_duration, readiness_delay)
    end
end

RollingShutterAcquisitionDefinition(exposure_duration::PlantDuration;
    readiness_delay::PlantDuration=zero(PlantDuration)) =
    RollingShutterAcquisitionDefinition(exposure_duration, readiness_delay)

mutable struct _RollingShutterAcquisitionBinding end

struct PreparedRollingShutterAcquisition{
    D<:Detector,
    P<:DetectorAcquisitionPlan,
    RP<:FrameReadoutProducts,
    M<:_RollingShutterEventMode,
    T<:AbstractFloat,
} <: AbstractPreparedDetectorAcquisition
    binding::_RollingShutterAcquisitionBinding
    detector::D
    plan::P
    readout_products::RP
    definition::RollingShutterAcquisitionDefinition
    mode::M
    exposure_seconds::T
    detection_efficiency::T
    line_duration::PlantDuration
    row_group_size::Int
    row_count::Int
    band_count::Int
end

mutable struct RollingShutterAcquisitionState <:
    AbstractDetectorAcquisitionEventState
    binding::_RollingShutterAcquisitionBinding
    status::DetectorAcquisitionStatus
    sequence::UInt64
    frame_start::PlantTimestamp
    integrated_through::PlantTimestamp
    readout_complete::PlantTimestamp
    readiness::PlantTimestamp
    opened_bands::Int
    closed_bands::Int
end

function RollingShutterAcquisitionState(
    prepared::PreparedRollingShutterAcquisition)
    origin = zero(PlantTimestamp)
    return RollingShutterAcquisitionState(prepared.binding,
        DetectorAcquisitionReady, UInt64(0), origin, origin, origin, origin,
        0, 0)
end

@inline detector_acquisition_status(state::RollingShutterAcquisitionState) =
    state.status
@inline detector_acquisition_sequence(state::RollingShutterAcquisitionState) =
    state.sequence
@inline exposure_start_timestamp(state::RollingShutterAcquisitionState) =
    state.frame_start
@inline integrated_through_timestamp(state::RollingShutterAcquisitionState) =
    state.integrated_through
@inline readout_complete_timestamp(state::RollingShutterAcquisitionState) =
    state.readout_complete
@inline acquisition_readiness_timestamp(
    state::RollingShutterAcquisitionState) = state.readiness
@inline rolling_band_count(prepared::PreparedRollingShutterAcquisition) =
    prepared.band_count
@inline rolling_opened_band_count(state::RollingShutterAcquisitionState) =
    state.opened_bands
@inline rolling_closed_band_count(state::RollingShutterAcquisitionState) =
    state.closed_bands

function rolling_band_rows(prepared::PreparedRollingShutterAcquisition,
    index::Integer)
    1 <= index <= prepared.band_count ||
        _detector_acquisition_event_error(:invalid_band_index,
            "rolling-shutter band index is outside the prepared frame")
    first_row = (Int(index) - 1) * prepared.row_group_size + 1
    last_row = min(Int(index) * prepared.row_group_size, prepared.row_count)
    return first_row:last_row
end

@inline function _rolling_frame_readout_offset(
    prepared::PreparedRollingShutterAcquisition)
    return prepared.definition.exposure_duration +
        prepared.line_duration * prepared.band_count
end

@inline function _rolling_band_open_timestamp(
    ::_RollingExposureEventMode,
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, index::Int)
    return state.frame_start + prepared.line_duration * (index - 1)
end

@inline function _rolling_band_open_timestamp(
    ::_GlobalResetEventMode,
    ::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, ::Int)
    return state.frame_start
end

@inline function rolling_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, index::Integer)
    _require_rolling_shutter_event_binding(prepared, state)
    rolling_band_rows(prepared, index)
    return _rolling_band_open_timestamp(prepared.mode, prepared, state,
        Int(index))
end

@inline function rolling_band_close_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, index::Integer)
    _require_rolling_shutter_event_binding(prepared, state)
    rolling_band_rows(prepared, index)
    return _rolling_band_close_timestamp(prepared.mode, prepared, state,
        Int(index))
end

@inline function _rolling_band_close_timestamp(
    ::_RollingExposureEventMode,
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, index::Int)
    return _rolling_band_open_timestamp(prepared.mode, prepared, state,
        index) + prepared.definition.exposure_duration
end

@inline function _rolling_band_close_timestamp(
    ::_GlobalResetEventMode,
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, index::Int)
    return state.frame_start + prepared.definition.exposure_duration +
        prepared.line_duration * (index - 1)
end

@inline function _next_rolling_band_open_timestamp(
    ::_RollingExposureEventMode,
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    state.status == DetectorExposureActive || return nothing
    state.opened_bands < prepared.band_count || return nothing
    return _rolling_band_open_timestamp(prepared.mode, prepared, state,
        state.opened_bands + 1)
end

@inline _next_rolling_band_open_timestamp(
    ::_GlobalResetEventMode, ::PreparedRollingShutterAcquisition,
    ::RollingShutterAcquisitionState) = nothing

@inline function next_rolling_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    _require_rolling_shutter_event_binding(prepared, state)
    return _next_rolling_band_open_timestamp(prepared.mode, prepared, state)
end

@inline function next_rolling_band_close_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorExposureActive || return nothing
    state.closed_bands < prepared.band_count || return nothing
    return rolling_band_close_timestamp(prepared, state,
        state.closed_bands + 1)
end

@inline _require_rolling_timing(timing::RollingShutter) = timing

function _require_rolling_timing(timing::AbstractFrameTimingModel)
    _detector_acquisition_event_error(:unsupported_timing,
        "rolling-shutter event acquisition requires a RollingShutter detector timing model")
end

function _require_rolling_sensor(::CMOSSensor)
    return nothing
end

function _require_rolling_sensor(sensor::FrameSensorType)
    _detector_acquisition_event_error(:unsupported_sensor,
        "rolling-shutter event acquisition requires a CMOS-family frame " *
        "sensor; got $(typeof(sensor))")
end

function prepare_rolling_shutter_acquisition(det::Detector,
    map::IntensityMap, definition::RollingShutterAcquisitionDefinition;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    plan = prepare_detector_acquisition(det, map;
        normalized_to_photon_rate=normalized_to_photon_rate)
    return prepare_rolling_shutter_acquisition(det, map, plan, definition)
end

function prepare_rolling_shutter_acquisition(det::Detector,
    map::IntensityMap, plan::DetectorAcquisitionPlan,
    definition::RollingShutterAcquisitionDefinition)
    timing = _require_rolling_timing(det.params.timing_model)
    _require_rolling_sensor(det.params.sensor)
    _require_detector_event_idle(det)
    _require_prepared_acquisition(det, map, plan)
    T = eltype(det.state.frame)
    exposure_seconds = plant_duration_seconds(definition.exposure_duration, T)
    isequal(det.params.integration_time, exposure_seconds) ||
        _detector_acquisition_event_error(:exposure_duration,
            "detector integration time must exactly match the prepared rolling-shutter duration")
    line_duration = _quantized_plant_duration(timing.line_time)
    row_count = size(det.state.frame, 1)
    band_count = cld(row_count, timing.row_group_size)
    detection_efficiency = plan.rate_scale * plan.quantum_efficiency
    isfinite(detection_efficiency) && detection_efficiency >= zero(T) ||
        _detector_acquisition_event_error(:detection_efficiency,
            "prepared detector rate scaling and quantum efficiency are not representable")
    return PreparedRollingShutterAcquisition(
        _RollingShutterAcquisitionBinding(), det, plan,
        det.state.readout_products, definition,
        _rolling_shutter_event_mode(timing.exposure_mode), exposure_seconds,
        detection_efficiency, line_duration, timing.row_group_size,
        row_count, band_count)
end

@inline function _require_rolling_shutter_event_binding(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    state.binding === prepared.binding ||
        _detector_acquisition_event_error(:foreign_state,
            "rolling-shutter state belongs to another prepared acquisition")
    det = prepared.detector
    plan = prepared.plan
    det.params === plan.detector_params && det.state === plan.detector_state &&
        det.state.frame === plan.detector_frame &&
        det.state.readout_products === prepared.readout_products &&
        typeof(backend(det)) === typeof(plan.detector_backend) ||
        _detector_acquisition_event_error(:prepared_binding,
            "rolling-shutter detector storage changed after preparation")
    return nothing
end

@inline function _require_rolling_shutter_progress(
    prepared::PreparedRollingShutterAcquisition,
    ::RollingShutterAcquisitionState)
    det = prepared.detector
    !det.state.readout_ready && iszero(det.state.integrated_time) ||
        _detector_acquisition_event_error(:detector_state_changed,
            "rolling-shutter detector state changed outside its event owner")
    return nothing
end

@inline _initial_rolling_opened_bands(
    ::_RollingExposureEventMode, ::Int) = 1
@inline _initial_rolling_opened_bands(
    ::_GlobalResetEventMode, band_count::Int) = band_count

function begin_exposure!(prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorAcquisitionReady ||
        _detector_acquisition_event_error(:retrigger,
            "cannot begin a rolling-shutter frame while its acquisition is busy")
    det = prepared.detector
    iszero(det.state.integrated_time) && det.state.readout_ready ||
        _detector_acquisition_event_error(:detector_busy,
            "detector is not ready to begin a rolling-shutter frame")
    state.sequence != typemax(UInt64) ||
        _detector_acquisition_event_error(:sequence_overflow,
            "rolling-shutter acquisition sequence is exhausted")
    iszero(state.sequence) || timestamp >= state.readiness ||
        _detector_acquisition_event_error(:time_regression,
            "rolling-shutter frame start precedes prior readiness")

    readout_complete = timestamp + _rolling_frame_readout_offset(prepared)
    readiness = readout_complete + prepared.definition.readiness_delay
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = false
    state.status = DetectorExposureActive
    state.sequence += UInt64(1)
    state.frame_start = timestamp
    state.integrated_through = timestamp
    state.readout_complete = readout_complete
    state.readiness = readiness
    state.opened_bands = _initial_rolling_opened_bands(prepared.mode,
        prepared.band_count)
    state.closed_bands = 0
    return nothing
end

@inline function _rolling_active_rows(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    first_band = state.closed_bands + 1
    last_band = state.opened_bands
    first_row = (first_band - 1) * prepared.row_group_size + 1
    last_row = min(last_band * prepared.row_group_size, prepared.row_count)
    return first_row:last_row
end

@inline function _next_rolling_transition_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    next_open = next_rolling_band_open_timestamp(prepared, state)
    next_close = next_rolling_band_close_timestamp(prepared, state)
    next_open === nothing && return next_close
    next_close === nothing && return next_open
    return min(next_open, next_close)
end

function accumulate_rolling_exposure_interval!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, start::PlantTimestamp,
    stop::PlantTimestamp, rng::AbstractRNG)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "rolling integration requires an active frame")
    _require_rolling_shutter_progress(prepared, state)
    start == state.integrated_through ||
        _detector_acquisition_event_error(:noncontiguous_interval,
            "rolling integration must begin at the prior interval boundary")
    start < stop || _detector_acquisition_event_error(:invalid_interval,
        "rolling integration interval must have positive duration")
    next_transition = _next_rolling_transition_timestamp(prepared, state)
    next_transition === nothing || stop <= next_transition ||
        _detector_acquisition_event_error(:missed_band_transition,
            "rolling integration crosses a pending row-band transition")

    det = prepared.detector
    T = eltype(det.state.frame)
    interval_seconds = plant_duration_seconds(stop - start, T)
    capture_signal_pipeline!(det, prepared.plan.input_values, rng,
        interval_seconds, prepared.detection_efficiency, false,
        prepared.exposure_seconds)
    accumulate_incremental_charge_generation!(det, rng, interval_seconds)
    rows = _rolling_active_rows(prepared, state)
    @views det.state.accum_buffer[rows, :] .+= det.state.frame[rows, :]
    advance_thermal!(det, interval_seconds)
    state.integrated_through = stop
    return nothing
end

function open_next_rolling_band!(
    prepared::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_RollingExposureEventMode},
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "row-band open requires an active rolling frame")
    _require_rolling_shutter_progress(prepared, state)
    expected = next_rolling_band_open_timestamp(prepared, state)
    expected === nothing && _detector_acquisition_event_error(
        :rolling_open_complete, "all rolling row bands are already open")
    timestamp == expected || _detector_acquisition_event_error(
        :rolling_open_not_due,
        "row-band open does not match the next prepared timestamp")
    state.integrated_through == timestamp ||
        _detector_acquisition_event_error(:incomplete_integration,
            "active row bands must be integrated through the next open boundary")
    state.opened_bands += 1
    return rolling_band_rows(prepared, state.opened_bands)
end

function open_next_rolling_band!(
    ::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_GlobalResetEventMode},
    ::RollingShutterAcquisitionState, ::PlantTimestamp)
    _detector_acquisition_event_error(:global_reset,
        "global-reset rolling shutter opens every row band at frame start")
end

function close_next_rolling_band!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "row-band close requires an active rolling frame")
    _require_rolling_shutter_progress(prepared, state)
    expected = next_rolling_band_close_timestamp(prepared, state)
    expected === nothing && _detector_acquisition_event_error(
        :rolling_close_complete, "all rolling row bands are already closed")
    timestamp == expected || _detector_acquisition_event_error(
        :rolling_close_not_due,
        "row-band close does not match the next prepared timestamp")
    state.integrated_through == timestamp ||
        _detector_acquisition_event_error(:incomplete_integration,
            "active row bands must be integrated through the close boundary")
    state.closed_bands += 1
    rows = rolling_band_rows(prepared, state.closed_bands)
    state.closed_bands == prepared.band_count &&
        (state.status = DetectorReadoutPending)
    return rows
end

function complete_readout!(prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorReadoutPending ||
        _detector_acquisition_event_error(:invalid_transition,
            "rolling readout completion requires every row band to be closed")
    _require_rolling_shutter_progress(prepared, state)
    timestamp == state.readout_complete ||
        _detector_acquisition_event_error(:readout_not_due,
            "rolling readout completion does not match its prepared timestamp")
    det = prepared.detector
    copyto!(det.state.frame, det.state.accum_buffer)
    finalize_incremental_capture!(det, rng, prepared.exposure_seconds)
    output = write_output!(det)
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    state.status = DetectorReadoutComplete
    return output
end

function mark_acquisition_ready!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_rolling_shutter_event_binding(prepared, state)
    state.status == DetectorReadoutComplete ||
        _detector_acquisition_event_error(:invalid_transition,
            "rolling readiness requires complete frame readout")
    timestamp == state.readiness ||
        _detector_acquisition_event_error(:readiness_not_due,
            "rolling readiness does not match its prepared timestamp")
    det = prepared.detector
    iszero(det.state.integrated_time) && !det.state.readout_ready ||
        _detector_acquisition_event_error(:detector_state_changed,
            "rolling detector state changed before readiness")
    det.state.readout_ready = true
    state.status = DetectorAcquisitionReady
    return output_frame(det)
end

"""
    FrameTransferAcquisitionDefinition(exposure_duration;
        readout_duration=PlantDuration(0))

Exact event definition for an EMCCD frame-transfer acquisition. The sensor's
`FrameTransferAcquisition` supplies the image-to-storage transfer duration.
One prepared storage frame permits image-area integration to overlap storage
readout without aliasing either frame.
"""
struct FrameTransferAcquisitionDefinition <:
    AbstractDetectorAcquisitionEventDefinition
    exposure_duration::PlantDuration
    readout_duration::PlantDuration

    function FrameTransferAcquisitionDefinition(
        exposure_duration::PlantDuration, readout_duration::PlantDuration)
        iszero(exposure_duration) && _detector_acquisition_event_error(
            :invalid_definition,
            "frame-transfer exposure duration must be > 0 ns")
        return new(exposure_duration, readout_duration)
    end
end

FrameTransferAcquisitionDefinition(exposure_duration::PlantDuration;
    readout_duration::PlantDuration=zero(PlantDuration)) =
    FrameTransferAcquisitionDefinition(exposure_duration, readout_duration)

mutable struct _FrameTransferAcquisitionBinding end

struct PreparedFrameTransferAcquisition{
    D<:Detector,
    P<:DetectorAcquisitionPlan,
    RP<:FrameReadoutProducts,
    A<:AbstractMatrix,
    T<:AbstractFloat,
} <: AbstractPreparedDetectorAcquisition
    binding::_FrameTransferAcquisitionBinding
    detector::D
    plan::P
    readout_products::RP
    definition::FrameTransferAcquisitionDefinition
    storage_frame::A
    exposure_seconds::T
    detection_efficiency::T
    transfer_duration::PlantDuration
end

@enum _FrameTransferImageStatus::UInt8 begin
    _FrameTransferImageReady = 0x01
    _FrameTransferImageActive = 0x02
    _FrameTransferPending = 0x03
end

@enum _FrameTransferStorageStatus::UInt8 begin
    _FrameTransferStorageEmpty = 0x01
    _FrameTransferStorageReadout = 0x02
end

mutable struct FrameTransferAcquisitionState <:
    AbstractDetectorAcquisitionEventState
    binding::_FrameTransferAcquisitionBinding
    image_status::_FrameTransferImageStatus
    storage_status::_FrameTransferStorageStatus
    sequence::UInt64
    image_sequence::UInt64
    storage_sequence::UInt64
    product_sequence::UInt64
    exposure_start::PlantTimestamp
    exposure_close::PlantTimestamp
    integrated_through::PlantTimestamp
    transfer_complete::PlantTimestamp
    storage_readout_complete::PlantTimestamp
    product_ready::PlantTimestamp
end

function FrameTransferAcquisitionState(
    prepared::PreparedFrameTransferAcquisition)
    origin = zero(PlantTimestamp)
    return FrameTransferAcquisitionState(prepared.binding,
        _FrameTransferImageReady, _FrameTransferStorageEmpty, UInt64(0),
        UInt64(0), UInt64(0), UInt64(0), origin, origin, origin, origin,
        origin, origin)
end

@inline frame_transfer_storage_capacity(
    ::PreparedFrameTransferAcquisition) = 1
@inline detector_acquisition_sequence(state::FrameTransferAcquisitionState) =
    state.sequence
@inline frame_transfer_image_sequence(state::FrameTransferAcquisitionState) =
    state.image_sequence
@inline frame_transfer_storage_sequence(
    state::FrameTransferAcquisitionState) = state.storage_sequence
@inline frame_transfer_product_sequence(
    state::FrameTransferAcquisitionState) = state.product_sequence
@inline exposure_start_timestamp(state::FrameTransferAcquisitionState) =
    state.exposure_start
@inline exposure_close_timestamp(state::FrameTransferAcquisitionState) =
    state.exposure_close
@inline integrated_through_timestamp(state::FrameTransferAcquisitionState) =
    state.integrated_through
@inline frame_transfer_complete_timestamp(
    state::FrameTransferAcquisitionState) = state.transfer_complete
@inline readout_complete_timestamp(state::FrameTransferAcquisitionState) =
    state.storage_readout_complete
@inline acquisition_product_ready_timestamp(
    state::FrameTransferAcquisitionState) = state.product_ready

@inline _require_frame_transfer_mode(mode::FrameTransferAcquisition) = mode

function _require_frame_transfer_mode(mode::AbstractEMCCDAcquisitionMode)
    _detector_acquisition_event_error(:unsupported_acquisition_mode,
        "frame-transfer event acquisition requires FrameTransferAcquisition; got $(typeof(mode))")
end

@inline function _require_frame_transfer_sensor(sensor::EMCCDSensor)
    return _require_frame_transfer_mode(sensor.acquisition_mode)
end

function _require_frame_transfer_sensor(sensor::FrameSensorType)
    _detector_acquisition_event_error(:unsupported_sensor,
        "frame-transfer event acquisition requires an EMCCD sensor; got $(typeof(sensor))")
end

function prepare_frame_transfer_acquisition(det::Detector,
    map::IntensityMap, definition::FrameTransferAcquisitionDefinition;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    plan = prepare_detector_acquisition(det, map;
        normalized_to_photon_rate=normalized_to_photon_rate)
    return prepare_frame_transfer_acquisition(det, map, plan, definition)
end

function prepare_frame_transfer_acquisition(det::Detector,
    map::IntensityMap, plan::DetectorAcquisitionPlan,
    definition::FrameTransferAcquisitionDefinition)
    _require_global_shutter_acquisition(det)
    mode = _require_frame_transfer_sensor(det.params.sensor)
    _require_detector_event_idle(det)
    _require_prepared_acquisition(det, map, plan)
    T = eltype(det.state.frame)
    exposure_seconds = plant_duration_seconds(definition.exposure_duration, T)
    isequal(det.params.integration_time, exposure_seconds) ||
        _detector_acquisition_event_error(:exposure_duration,
            "detector integration time must exactly match the prepared frame-transfer duration")
    transfer_duration = _quantized_plant_duration(mode.transfer_time)
    storage_frame = similar(det.state.frame)
    fill!(storage_frame, zero(eltype(storage_frame)))
    for aliased in (det.state.frame, det.state.accum_buffer,
            output_frame(det))
        Base.mightalias(storage_frame, aliased) &&
            _detector_acquisition_event_error(:storage_alias,
                "frame-transfer storage must not alias detector image or output storage")
    end
    detection_efficiency = plan.rate_scale * plan.quantum_efficiency
    isfinite(detection_efficiency) && detection_efficiency >= zero(T) ||
        _detector_acquisition_event_error(:detection_efficiency,
            "prepared detector rate scaling and quantum efficiency are not representable")
    return PreparedFrameTransferAcquisition(
        _FrameTransferAcquisitionBinding(), det, plan,
        det.state.readout_products, definition, storage_frame,
        exposure_seconds, detection_efficiency, transfer_duration)
end

@inline function _require_frame_transfer_event_binding(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState)
    state.binding === prepared.binding ||
        _detector_acquisition_event_error(:foreign_state,
            "frame-transfer state belongs to another prepared acquisition")
    det = prepared.detector
    plan = prepared.plan
    det.params === plan.detector_params && det.state === plan.detector_state &&
        det.state.frame === plan.detector_frame &&
        det.state.readout_products === prepared.readout_products &&
        typeof(backend(det)) === typeof(plan.detector_backend) ||
        _detector_acquisition_event_error(:prepared_binding,
            "frame-transfer detector storage changed after preparation")
    return nothing
end

@inline function _require_frame_transfer_progress(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState)
    det = prepared.detector
    expected_readout_ready = frame_transfer_storage_empty(state) &&
        frame_transfer_image_ready(state)
    iszero(det.state.integrated_time) &&
        det.state.readout_ready == expected_readout_ready ||
        _detector_acquisition_event_error(:detector_state_changed,
            "frame-transfer detector state changed outside its event owner")
    return nothing
end

@inline frame_transfer_image_ready(state::FrameTransferAcquisitionState) =
    state.image_status == _FrameTransferImageReady
@inline frame_transfer_storage_empty(state::FrameTransferAcquisitionState) =
    state.storage_status == _FrameTransferStorageEmpty
@inline frame_transfer_readout_pending(state::FrameTransferAcquisitionState) =
    state.storage_status == _FrameTransferStorageReadout

function begin_exposure!(prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp)
    _require_frame_transfer_event_binding(prepared, state)
    _require_frame_transfer_progress(prepared, state)
    frame_transfer_image_ready(state) ||
        _detector_acquisition_event_error(:retrigger,
            "frame-transfer image area is not ready for another exposure")
    state.sequence != typemax(UInt64) ||
        _detector_acquisition_event_error(:sequence_overflow,
            "frame-transfer acquisition sequence is exhausted")
    iszero(state.sequence) || timestamp >= state.transfer_complete ||
        _detector_acquisition_event_error(:time_regression,
            "frame-transfer exposure begins before the image area became available")
    future_close = timestamp + prepared.definition.exposure_duration
    future_transfer = future_close + prepared.transfer_duration
    if frame_transfer_readout_pending(state)
        future_transfer > state.storage_readout_complete ||
            _detector_acquisition_event_error(:storage_capacity,
                "next frame transfer would reach the one-frame storage " *
                "area before prior readout completes")
    end
    det = prepared.detector
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = false
    state.sequence += UInt64(1)
    state.image_sequence = state.sequence
    state.image_status = _FrameTransferImageActive
    state.exposure_start = timestamp
    state.exposure_close = future_close
    state.integrated_through = timestamp
    state.transfer_complete = future_transfer
    return nothing
end

function accumulate_exposure_interval!(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, start::PlantTimestamp,
    stop::PlantTimestamp, rng::AbstractRNG)
    _require_frame_transfer_event_binding(prepared, state)
    state.image_status == _FrameTransferImageActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "frame-transfer integration requires an active image area")
    _require_frame_transfer_progress(prepared, state)
    start == state.integrated_through ||
        _detector_acquisition_event_error(:noncontiguous_interval,
            "frame-transfer integration must begin at the prior boundary")
    start < stop || _detector_acquisition_event_error(:invalid_interval,
        "frame-transfer integration interval must have positive duration")
    stop <= state.exposure_close ||
        _detector_acquisition_event_error(:interval_after_close,
            "frame-transfer integration extends beyond exposure close")
    det = prepared.detector
    T = eltype(det.state.frame)
    interval_seconds = plant_duration_seconds(stop - start, T)
    exposure_start = start == state.exposure_start
    capture_signal_pipeline!(det, prepared.plan.input_values, rng,
        interval_seconds, prepared.detection_efficiency, exposure_start,
        prepared.exposure_seconds)
    accumulate_incremental_charge_generation!(det, rng, interval_seconds)
    det.state.accum_buffer .+= det.state.frame
    advance_thermal!(det, interval_seconds)
    state.integrated_through = stop
    return nothing
end

function close_exposure!(prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp)
    _require_frame_transfer_event_binding(prepared, state)
    state.image_status == _FrameTransferImageActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "frame-transfer close requires an active image area")
    _require_frame_transfer_progress(prepared, state)
    timestamp == state.exposure_close ||
        _detector_acquisition_event_error(:close_not_due,
            "frame-transfer close does not match its prepared timestamp")
    state.integrated_through == timestamp ||
        _detector_acquisition_event_error(:incomplete_integration,
            "frame-transfer exposure must be integrated through close")
    state.image_status = _FrameTransferPending
    return nothing
end

function complete_frame_transfer!(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp)
    _require_frame_transfer_event_binding(prepared, state)
    _require_frame_transfer_progress(prepared, state)
    state.image_status == _FrameTransferPending ||
        _detector_acquisition_event_error(:invalid_transition,
            "image-to-storage transfer requires a closed image-area exposure")
    frame_transfer_storage_empty(state) ||
        _detector_acquisition_event_error(:storage_capacity,
            "frame-transfer storage is still occupied by prior readout")
    timestamp == state.transfer_complete ||
        _detector_acquisition_event_error(:transfer_not_due,
            "frame transfer does not match its prepared timestamp")
    det = prepared.detector
    copyto!(prepared.storage_frame, det.state.accum_buffer)
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    state.storage_status = _FrameTransferStorageReadout
    state.storage_sequence = state.image_sequence
    state.storage_readout_complete = timestamp +
        prepared.definition.readout_duration
    state.image_status = _FrameTransferImageReady
    return prepared.storage_frame
end

function complete_readout!(prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _require_frame_transfer_event_binding(prepared, state)
    _require_frame_transfer_progress(prepared, state)
    frame_transfer_readout_pending(state) ||
        _detector_acquisition_event_error(:invalid_transition,
            "frame-transfer readout requires one occupied storage frame")
    timestamp == state.storage_readout_complete ||
        _detector_acquisition_event_error(:readout_not_due,
            "frame-transfer readout does not match its prepared timestamp")
    det = prepared.detector
    copyto!(det.state.frame, prepared.storage_frame)
    finalize_incremental_capture!(det, rng, prepared.exposure_seconds)
    output = write_output!(det)
    fill!(prepared.storage_frame, zero(eltype(prepared.storage_frame)))
    state.product_sequence = state.storage_sequence
    state.product_ready = timestamp
    state.storage_sequence = UInt64(0)
    state.storage_status = _FrameTransferStorageEmpty
    det.state.readout_ready = frame_transfer_image_ready(state)
    return output
end
