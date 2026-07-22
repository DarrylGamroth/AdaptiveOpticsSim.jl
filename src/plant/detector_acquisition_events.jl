const _DETECTOR_ACQUISITION_EVENT_COMPONENT = :detector_acquisition_events

abstract type AbstractDetectorAcquisitionEventDefinition end
abstract type AbstractPreparedDetectorAcquisition end
abstract type AbstractDetectorAcquisitionEventState end

@noinline function _detector_acquisition_event_error(reason::Symbol,
    message::String)
    throw(DetectorAcquisitionError(_DETECTOR_ACQUISITION_EVENT_COMPONENT,
        reason, message))
end

"""
    GlobalShutterAcquisitionDefinition(exposure_duration;
        readout_duration=PlantDuration(0),
        readiness_delay=PlantDuration(0))

Immutable virtual-time definition for one global-shutter frame acquisition.
`exposure_duration` is the half-open integration interval. `readout_duration`
places readout completion after exposure close, and `readiness_delay` places the
next-acquisition readiness transition after readout completion. The initial
event model rejects a trigger received while the detector is not ready.
"""
struct GlobalShutterAcquisitionDefinition <:
    AbstractDetectorAcquisitionEventDefinition
    exposure_duration::PlantDuration
    readout_duration::PlantDuration
    readiness_delay::PlantDuration

    function GlobalShutterAcquisitionDefinition(
        exposure_duration::PlantDuration,
        readout_duration::PlantDuration,
        readiness_delay::PlantDuration)
        iszero(exposure_duration) && _detector_acquisition_event_error(
            :invalid_definition,
            "global-shutter exposure duration must be > 0 ns")
        return new(exposure_duration, readout_duration, readiness_delay)
    end
end

GlobalShutterAcquisitionDefinition(exposure_duration::PlantDuration;
    readout_duration::PlantDuration=zero(PlantDuration),
    readiness_delay::PlantDuration=zero(PlantDuration)) =
    GlobalShutterAcquisitionDefinition(exposure_duration, readout_duration,
        readiness_delay)

abstract type _DetectorEventReadoutStyle end
struct _PostExposureDetectorReadout <: _DetectorEventReadoutStyle end
struct _ScheduledUpTheRampReadout <: _DetectorEventReadoutStyle end

@inline _detector_event_readout_style(::FrameSensorType) =
    _PostExposureDetectorReadout()

@inline _detector_event_readout_style(
    ::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling}) =
    _ScheduledUpTheRampReadout()

mutable struct _GlobalShutterAcquisitionBinding end

struct _PreparedGlobalShutterAcquisitionToken end
const _PREPARED_GLOBAL_SHUTTER_ACQUISITION_TOKEN =
    _PreparedGlobalShutterAcquisitionToken()

"""
Prepared detector/map binding and fixed event schedule for one global-shutter
acquisition owner. Construct with [`prepare_global_shutter_acquisition`](@ref).
"""
struct PreparedGlobalShutterAcquisition{
    D<:Detector,
    P<:DetectorAcquisitionPlan,
    RP<:FrameReadoutProducts,
    R<:_DetectorEventReadoutStyle,
    T<:AbstractFloat,
} <: AbstractPreparedDetectorAcquisition
    binding::_GlobalShutterAcquisitionBinding
    detector::D
    plan::P
    readout_products::RP
    definition::GlobalShutterAcquisitionDefinition
    readout_style::R
    exposure_seconds::T
    detection_efficiency::T
    read_offsets::Memory{PlantDuration}

    function PreparedGlobalShutterAcquisition(
        ::_PreparedGlobalShutterAcquisitionToken,
        binding::_GlobalShutterAcquisitionBinding, detector::D, plan::P,
        readout_products::RP, definition::GlobalShutterAcquisitionDefinition,
        readout_style::R, exposure_seconds::T, detection_efficiency::T,
        read_offsets::Memory{PlantDuration}) where {
        D<:Detector,P<:DetectorAcquisitionPlan,
        RP<:FrameReadoutProducts,R<:_DetectorEventReadoutStyle,
        T<:AbstractFloat,
    }
        return new{D,P,RP,R,T}(binding, detector, plan, readout_products,
            definition, readout_style, exposure_seconds,
            detection_efficiency, read_offsets)
    end
end

@enum DetectorAcquisitionStatus::UInt8 begin
    DetectorAcquisitionReady = 0x01
    DetectorExposureActive = 0x02
    DetectorReadoutPending = 0x03
    DetectorReadoutComplete = 0x04
end

"""Separately owned, single-writer state for one prepared acquisition."""
mutable struct GlobalShutterAcquisitionState <:
    AbstractDetectorAcquisitionEventState
    binding::_GlobalShutterAcquisitionBinding
    status::DetectorAcquisitionStatus
    sequence::UInt64
    exposure_start::PlantTimestamp
    exposure_close::PlantTimestamp
    integrated_through::PlantTimestamp
    readout_complete::PlantTimestamp
    readiness::PlantTimestamp
    next_read_index::Int
end

function GlobalShutterAcquisitionState(
    prepared::PreparedGlobalShutterAcquisition)
    origin = zero(PlantTimestamp)
    return GlobalShutterAcquisitionState(prepared.binding,
        DetectorAcquisitionReady, UInt64(0), origin, origin, origin, origin,
        origin, 1)
end

@inline detector_acquisition_status(state::GlobalShutterAcquisitionState) =
    state.status
@inline detector_acquisition_sequence(state::GlobalShutterAcquisitionState) =
    state.sequence
@inline exposure_start_timestamp(state::GlobalShutterAcquisitionState) =
    state.exposure_start
@inline exposure_close_timestamp(state::GlobalShutterAcquisitionState) =
    state.exposure_close
@inline integrated_through_timestamp(state::GlobalShutterAcquisitionState) =
    state.integrated_through
@inline readout_complete_timestamp(state::GlobalShutterAcquisitionState) =
    state.readout_complete
@inline acquisition_readiness_timestamp(state::GlobalShutterAcquisitionState) =
    state.readiness

@inline nondestructive_read_count(
    prepared::PreparedGlobalShutterAcquisition) = length(prepared.read_offsets)

@inline function nondestructive_read_offset(
    prepared::PreparedGlobalShutterAcquisition, index::Integer)
    1 <= index <= length(prepared.read_offsets) ||
        _detector_acquisition_event_error(:invalid_read_index,
            "nondestructive-read index is outside the prepared schedule")
    return @inbounds prepared.read_offsets[index]
end

@inline function next_nondestructive_read_timestamp(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        state.status == DetectorReadoutPending || return nothing
    state.next_read_index <= length(prepared.read_offsets) || return nothing
    return state.exposure_start +
        @inbounds(prepared.read_offsets[state.next_read_index])
end

function _empty_detector_read_offsets()
    return Memory{PlantDuration}(undef, 0)
end

function _even_detector_read_offsets(exposure::PlantDuration, n_reads::Int)
    exposure.nanoseconds >= n_reads - 1 ||
        _detector_acquisition_event_error(:unrepresentable_read_schedule,
            "up-the-ramp exposure must contain at least one nanosecond per read interval")
    offsets = Memory{PlantDuration}(undef, n_reads)
    denominator = Int128(n_reads - 1)
    exposure_nanoseconds = Int128(exposure.nanoseconds)
    @inbounds for read_index in 1:n_reads
        numerator = exposure_nanoseconds * Int128(read_index - 1)
        offsets[read_index] = PlantDuration(Int64(numerator ÷ denominator))
    end
    return offsets
end

function _quantized_plant_duration(seconds::Real)
    isfinite(seconds) && seconds >= zero(seconds) ||
        _detector_acquisition_event_error(:invalid_read_duration,
            "detector read duration must be finite and nonnegative")
    scaled = big(seconds) * big(_PLANT_NANOSECONDS_PER_SECOND)
    scaled <= typemax(Int64) || _detector_acquisition_event_error(
        :unrepresentable_read_duration,
        "detector read duration exceeds the plant-time range")
    nanoseconds = round(Int64, scaled, RoundNearest)
    seconds > zero(seconds) && iszero(nanoseconds) &&
        _detector_acquisition_event_error(:unrepresentable_read_duration,
            "positive detector read duration is below plant-time resolution")
    return PlantDuration(nanoseconds)
end

@inline function _prepare_detector_read_offsets(
    ::_PostExposureDetectorReadout, ::FrameSensorType,
    ::Detector, ::GlobalShutterAcquisitionDefinition)
    return _empty_detector_read_offsets()
end

function _prepare_detector_read_offsets(::_ScheduledUpTheRampReadout,
    sensor::HgCdTeAvalancheArraySensor, det::Detector,
    definition::GlobalShutterAcquisitionDefinition)
    mode = sensor.sampling_mode
    offsets = _even_detector_read_offsets(definition.exposure_duration,
        mode.n_reads)
    T = eltype(det.state.frame)
    physical_read_duration = _quantized_plant_duration(sampling_read_time(sensor,
        size(det.state.frame), det.params.readout_window, T))
    @inbounds for read_index in 2:length(offsets)
        spacing = offsets[read_index] - offsets[read_index - 1]
        physical_read_duration <= spacing ||
            _detector_acquisition_event_error(:read_spacing,
                "up-the-ramp read duration exceeds one prepared read interval")
    end
    physical_read_duration <= definition.readout_duration ||
        _detector_acquisition_event_error(:readout_duration,
            "up-the-ramp readout duration must include the final physical read")
    return offsets
end

@inline _require_up_the_ramp_products(
    products::UpTheRampReadoutProducts) = products

function _require_up_the_ramp_products(::FrameReadoutProducts)
    _detector_acquisition_event_error(:prepared_products,
        "scheduled up-the-ramp acquisition requires prepared ramp products")
end

@inline function _prepare_detector_event_products!(
    ::_PostExposureDetectorReadout, ::Detector,
    ::Memory{PlantDuration}, ::Type{<:AbstractFloat})
    return nothing
end

function _prepare_detector_event_products!(::_ScheduledUpTheRampReadout,
    det::Detector, offsets::Memory{PlantDuration},
    ::Type{T}) where {T<:AbstractFloat}
    products = _require_up_the_ramp_products(
        ensure_up_the_ramp_products!(det, length(offsets)))
    @inbounds for read_index in eachindex(offsets)
        products.read_times[read_index] =
            plant_duration_seconds(offsets[read_index], T)
    end
    return nothing
end

@inline function _require_global_shutter_acquisition(det::Detector)
    is_global_shutter(det.params.timing_model) ||
        _detector_acquisition_event_error(:unsupported_timing,
            "explicit detector acquisition events currently require a global shutter")
    return nothing
end

@inline function _require_detector_event_idle(det::Detector)
    iszero(det.state.integrated_time) && det.state.readout_ready ||
        _detector_acquisition_event_error(:detector_busy,
            "detector must be idle before event acquisition preparation")
    return nothing
end

"""
    prepare_global_shutter_acquisition(detector, map, definition;
        normalized_to_photon_rate=nothing)

Prepare an exact virtual-time lifecycle over one existing detector intensity-map
contract. The detector's floating integration time must equal the definition's
nanosecond duration when converted to its numeric type. Scheduled up-the-ramp
read instants use endpoint-preserving integer-nanosecond floor quantization;
the detector's floating read duration is rounded to the nearest nanosecond for
spacing and readout validation.
"""
function prepare_global_shutter_acquisition(det::Detector, map::IntensityMap,
    definition::GlobalShutterAcquisitionDefinition;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    _require_global_shutter_acquisition(det)
    plan = prepare_detector_acquisition(det, map;
        normalized_to_photon_rate=normalized_to_photon_rate)
    return prepare_global_shutter_acquisition(det, map, plan, definition)
end

function prepare_global_shutter_acquisition(det::Detector, map::IntensityMap,
    plan::DetectorAcquisitionPlan,
    definition::GlobalShutterAcquisitionDefinition)
    _require_global_shutter_acquisition(det)
    _require_detector_event_idle(det)
    _require_prepared_acquisition(det, map, plan)
    T = eltype(det.state.frame)
    exposure_seconds = plant_duration_seconds(definition.exposure_duration, T)
    isequal(det.params.integration_time, exposure_seconds) ||
        _detector_acquisition_event_error(:exposure_duration,
            "detector integration time must exactly match the prepared plant duration")
    readout_style = _detector_event_readout_style(det.params.sensor)
    read_offsets = _prepare_detector_read_offsets(readout_style,
        det.params.sensor, det, definition)
    _prepare_detector_event_products!(readout_style, det, read_offsets, T)
    detection_efficiency = plan.rate_scale * plan.quantum_efficiency
    isfinite(detection_efficiency) && detection_efficiency >= zero(T) ||
        _detector_acquisition_event_error(:detection_efficiency,
            "prepared detector rate scaling and quantum efficiency are not representable")
    return PreparedGlobalShutterAcquisition(
        _PREPARED_GLOBAL_SHUTTER_ACQUISITION_TOKEN,
        _GlobalShutterAcquisitionBinding(), det, plan,
        det.state.readout_products, definition, readout_style,
        exposure_seconds, detection_efficiency, read_offsets)
end

@inline function _require_detector_event_binding(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    state.binding === prepared.binding ||
        _detector_acquisition_event_error(:foreign_state,
            "acquisition state belongs to another prepared detector acquisition")
    det = prepared.detector
    plan = prepared.plan
    det.params === plan.detector_params &&
        det.state === plan.detector_state &&
        det.state.frame === plan.detector_frame &&
        det.state.readout_products === prepared.readout_products &&
        typeof(backend(det)) === typeof(plan.detector_backend) ||
        _detector_acquisition_event_error(:prepared_binding,
            "detector storage changed after event acquisition preparation")
    return nothing
end

@inline function _expected_integrated_time(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    elapsed = state.integrated_through - state.exposure_start
    return plant_duration_seconds(elapsed,
        eltype(prepared.detector.state.frame))
end

@inline function _require_detector_event_progress(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    det = prepared.detector
    !det.state.readout_ready && isequal(det.state.integrated_time,
        _expected_integrated_time(prepared, state)) ||
        _detector_acquisition_event_error(:detector_state_changed,
            "detector integration state changed outside its prepared event owner")
    return nothing
end

@inline function _clear_detector_event_products!(
    ::_PostExposureDetectorReadout, ::PreparedGlobalShutterAcquisition)
    return nothing
end

function _clear_detector_event_products!(::_ScheduledUpTheRampReadout,
    prepared::PreparedGlobalShutterAcquisition)
    products = _require_up_the_ramp_products(prepared.readout_products)
    fill!(products.slope_frame, zero(eltype(products.slope_frame)))
    fill!(products.intercept_frame, zero(eltype(products.intercept_frame)))
    fill!(products.integrated_frame, zero(eltype(products.integrated_frame)))
    fill!(products.read_cube, zero(eltype(products.read_cube)))
    fill!(products.workspace_slope, zero(eltype(products.workspace_slope)))
    fill!(products.workspace_intercept,
        zero(eltype(products.workspace_intercept)))
    fill!(products.workspace_integrated,
        zero(eltype(products.workspace_integrated)))
    fill!(products.workspace_cube, zero(eltype(products.workspace_cube)))
    T = eltype(products.read_times)
    @inbounds for read_index in eachindex(prepared.read_offsets)
        products.read_times[read_index] = plant_duration_seconds(
            prepared.read_offsets[read_index], T)
    end
    return nothing
end

"""Begin one exposure at an exact plant timestamp."""
function begin_exposure!(prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorAcquisitionReady ||
        _detector_acquisition_event_error(:retrigger,
            "cannot begin an exposure while its detector acquisition is busy")
    det = prepared.detector
    iszero(det.state.integrated_time) && det.state.readout_ready ||
        _detector_acquisition_event_error(:detector_busy,
            "detector is not ready to begin an event-driven exposure")
    state.sequence != typemax(UInt64) ||
        _detector_acquisition_event_error(:sequence_overflow,
            "detector acquisition sequence is exhausted")
    iszero(state.sequence) || timestamp >= state.readiness ||
        _detector_acquisition_event_error(:time_regression,
            "exposure begin timestamp precedes the prior readiness transition")

    exposure_close = timestamp + prepared.definition.exposure_duration
    readout_complete = exposure_close + prepared.definition.readout_duration
    readiness = readout_complete + prepared.definition.readiness_delay
    next_sequence = state.sequence + UInt64(1)

    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    _clear_detector_event_products!(prepared.readout_style, prepared)
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = false
    state.status = DetectorExposureActive
    state.sequence = next_sequence
    state.exposure_start = timestamp
    state.exposure_close = exposure_close
    state.integrated_through = timestamp
    state.readout_complete = readout_complete
    state.readiness = readiness
    state.next_read_index = 1
    return nothing
end

@inline function _require_interval_before_next_read(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, stop::PlantTimestamp)
    state.next_read_index <= length(prepared.read_offsets) || return nothing
    next_read = state.exposure_start +
        @inbounds(prepared.read_offsets[state.next_read_index])
    stop <= next_read || _detector_acquisition_event_error(
        :missed_nondestructive_read,
        "integration interval crosses a pending nondestructive-read boundary")
    return nothing
end

@inline function _apply_detector_event_interval_statistics!(
    ::_PostExposureDetectorReadout, ::Detector, ::AbstractRNG)
    return nothing
end

@inline function _apply_detector_event_interval_statistics!(
    ::_ScheduledUpTheRampReadout, det::Detector, rng::AbstractRNG)
    sensor = det.params.sensor
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

"""
    accumulate_exposure_interval!(prepared, state, start, stop, rng)

Integrate the currently bound photon-arrival-rate map over exactly one positive
half-open interval `[start, stop)`. Intervals must be contiguous and cannot
cross exposure close or a pending nondestructive-read event.
"""
function accumulate_exposure_interval!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, start::PlantTimestamp,
    stop::PlantTimestamp, rng::AbstractRNG)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "integration requires an active detector exposure")
    _require_detector_event_progress(prepared, state)
    start == state.integrated_through ||
        _detector_acquisition_event_error(:noncontiguous_interval,
            "integration interval must begin at the previous interval boundary")
    start < stop || _detector_acquisition_event_error(:invalid_interval,
        "integration interval must have positive duration")
    stop <= state.exposure_close ||
        _detector_acquisition_event_error(:interval_after_close,
            "integration interval extends beyond exposure close")
    _require_interval_before_next_read(prepared, state, stop)

    det = prepared.detector
    T = eltype(det.state.frame)
    interval_seconds = plant_duration_seconds(stop - start, T)
    exposure_start = start == state.exposure_start
    capture_signal_pipeline!(det, prepared.plan.input_values, rng,
        interval_seconds, prepared.detection_efficiency, exposure_start,
        prepared.exposure_seconds)
    accumulate_incremental_charge_generation!(det, rng, interval_seconds)
    # Avalanche statistics belong to newly generated charge. Applying them
    # before accumulation preserves that realization in later nondestructive
    # reads instead of redrawing noise for charge already observed.
    _apply_detector_event_interval_statistics!(prepared.readout_style, det,
        rng)
    det.state.accum_buffer .+= det.state.frame
    advance_thermal!(det, interval_seconds)
    det.state.integrated_time = plant_duration_seconds(
        stop - state.exposure_start, T)
    state.integrated_through = stop
    return nothing
end

@inline function _take_nondestructive_read!(
    ::_PostExposureDetectorReadout,
    ::PreparedGlobalShutterAcquisition,
    ::GlobalShutterAcquisitionState, ::PlantTimestamp, ::AbstractRNG)
    _detector_acquisition_event_error(:nondestructive_read_unsupported,
        "prepared detector acquisition has no nondestructive-read events")
end

function _take_nondestructive_read!(::_ScheduledUpTheRampReadout,
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    active = state.status == DetectorExposureActive
    final_after_close = state.status == DetectorReadoutPending &&
        timestamp == state.exposure_close
    active || final_after_close ||
        _detector_acquisition_event_error(:invalid_transition,
            "nondestructive read requires an active exposure or its close boundary")
    _require_detector_event_progress(prepared, state)
    read_index = state.next_read_index
    read_index <= length(prepared.read_offsets) ||
        _detector_acquisition_event_error(:nondestructive_read_complete,
            "all prepared nondestructive reads have already been taken")
    expected_timestamp = state.exposure_start +
        @inbounds(prepared.read_offsets[read_index])
    timestamp == expected_timestamp ||
        _detector_acquisition_event_error(:nondestructive_read_not_due,
            "nondestructive read does not match the next prepared timestamp")
    state.integrated_through == timestamp ||
        _detector_acquisition_event_error(:incomplete_integration,
            "nondestructive read requires integration through its snapshot boundary")

    det = prepared.detector
    products = _require_up_the_ramp_products(prepared.readout_products)
    copyto!(det.state.frame, det.state.accum_buffer)
    finalize_charge_transport!(det, rng)
    target = @view products.workspace_cube[:, :, read_index]
    sample_frame_read!(det.params.sensor, det, target, det.state.frame,
        _raw_sampling_sigma(det), rng)
    if products.read_cube !== products.workspace_cube
        _copy_windowed_sampling_plane!(products.read_cube, read_index,
            products.workspace_cube, read_index, det)
    end
    state.next_read_index = read_index + 1
    return nothing
end

"""Take the next exact scheduled nondestructive charge snapshot."""
function take_nondestructive_read!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _require_detector_event_binding(prepared, state)
    return _take_nondestructive_read!(prepared.readout_style, prepared, state,
        timestamp, rng)
end

@inline function _require_no_missed_read_at_close(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    state.next_read_index <= length(prepared.read_offsets) || return nothing
    next_offset = @inbounds prepared.read_offsets[state.next_read_index]
    next_offset == prepared.definition.exposure_duration ||
        _detector_acquisition_event_error(:missing_nondestructive_read,
            "exposure cannot close after a missed nondestructive read")
    return nothing
end

"""Close an exposure at its exact prepared half-open boundary."""
function close_exposure!(prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorExposureActive ||
        _detector_acquisition_event_error(:invalid_transition,
            "close requires an active detector exposure")
    _require_detector_event_progress(prepared, state)
    timestamp == state.exposure_close ||
        _detector_acquisition_event_error(:close_not_due,
            "exposure close does not match its prepared timestamp")
    state.integrated_through == state.exposure_close ||
        _detector_acquisition_event_error(:incomplete_integration,
            "exposure cannot close before integration reaches its boundary")
    _require_no_missed_read_at_close(prepared, state)
    state.status = DetectorReadoutPending
    return nothing
end

@inline function _require_detector_reads_complete(
    ::_PostExposureDetectorReadout,
    ::PreparedGlobalShutterAcquisition,
    ::GlobalShutterAcquisitionState)
    return nothing
end

@inline function _require_detector_reads_complete(
    ::_ScheduledUpTheRampReadout,
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    state.next_read_index > length(prepared.read_offsets) ||
        _detector_acquisition_event_error(:missing_nondestructive_read,
            "readout cannot complete before every scheduled read is taken")
    return nothing
end

function _complete_detector_readout!(::_PostExposureDetectorReadout,
    prepared::PreparedGlobalShutterAcquisition, rng::AbstractRNG)
    det = prepared.detector
    copyto!(det.state.frame, det.state.accum_buffer)
    finalize_incremental_capture!(det, rng, prepared.exposure_seconds)
    return write_output!(det)
end

function _complete_detector_readout!(::_ScheduledUpTheRampReadout,
    prepared::PreparedGlobalShutterAcquisition, rng::AbstractRNG)
    det = prepared.detector
    products = _require_up_the_ramp_products(prepared.readout_products)
    finalize_scheduled_up_the_ramp_readout_products!(products, det,
        prepared.exposure_seconds)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(det.params.sensor, det,
        prepared.exposure_seconds)
    return write_output!(det)
end

"""Finalize detector readout at its exact prepared completion timestamp."""
function complete_readout!(prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorReadoutPending ||
        _detector_acquisition_event_error(:invalid_transition,
            "readout completion requires a closed detector exposure")
    _require_detector_event_progress(prepared, state)
    timestamp == state.readout_complete ||
        _detector_acquisition_event_error(:readout_not_due,
            "readout completion does not match its prepared timestamp")
    _require_detector_reads_complete(prepared.readout_style, prepared, state)

    output = _complete_detector_readout!(prepared.readout_style, prepared, rng)
    det = prepared.detector
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    det.state.integrated_time = zero(det.state.integrated_time)
    state.status = DetectorReadoutComplete
    return output
end

"""Publish acquisition readiness and permit the next exposure."""
function mark_acquisition_ready!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp)
    _require_detector_event_binding(prepared, state)
    state.status == DetectorReadoutComplete ||
        _detector_acquisition_event_error(:invalid_transition,
            "readiness requires completed detector readout")
    timestamp == state.readiness ||
        _detector_acquisition_event_error(:readiness_not_due,
            "acquisition readiness does not match its prepared timestamp")
    det = prepared.detector
    iszero(det.state.integrated_time) && !det.state.readout_ready ||
        _detector_acquisition_event_error(:detector_state_changed,
            "detector state changed before acquisition readiness")
    det.state.readout_ready = true
    state.status = DetectorAcquisitionReady
    return output_frame(det)
end
