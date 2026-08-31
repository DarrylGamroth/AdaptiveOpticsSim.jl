module AdaptiveOpticsSimPipeWireAOExt

import AdaptiveOpticsSim: AlgorithmGraphs
import PipeWireAO

const AG = AlgorithmGraphs
const PWAO = PipeWireAO

struct PipeWireAcquisitionTimeOrigin{Reference}
    nanoseconds::Int64
    timebase::PWAO.AcquisitionTimebase
    ptp_reference::Reference
end

struct PipeWireAcquisitionProvenance{Reference}
    identity::PWAO.AcquisitionIdentity
    exposure_start_nanoseconds::Int64
    timebase::PWAO.AcquisitionTimebase
    ptp_reference::Reference
    exposure_duration_nanoseconds::UInt64
end

@noinline function _capture_error(message::AbstractString)
    throw(AG.AlgorithmGraphError(String(message)))
end

function _acquisition_clock(metadata::PWAO.AcquisitionMetadata)
    start = PWAO.acquisition_exposure_start(metadata)
    start === nothing && _capture_error(
        "PipeWire acquisition metadata has no exposure start",
    )
    timebase = PWAO.acquisition_timebase(metadata)
    timebase === nothing && _capture_error(
        "PipeWire acquisition metadata has no exposure timebase",
    )
    reference = PWAO.acquisition_ptp_reference(metadata)
    if timebase == PWAO.ACQUISITION_TIMEBASE_TAI
        reference === nothing && _capture_error(
            "PipeWire TAI acquisition metadata has no PTP reference",
        )
    elseif reference !== nothing
        _capture_error(
            "PipeWire monotonic acquisition metadata carries a PTP reference",
        )
    end
    return start, timebase, reference
end

function AlgorithmGraphs.capture_model_time_origin(
    metadata::PWAO.AcquisitionMetadata,
)
    PWAO.acquisition_identity(metadata) === nothing && _capture_error(
        "PipeWire model-time capture requires an acquisition identity",
    )
    start, timebase, reference = _acquisition_clock(metadata)
    start.nanoseconds >= 0 || _capture_error(
        "PipeWire exposure start cannot map from a negative clock reading",
    )
    return PipeWireAcquisitionTimeOrigin(
        start.nanoseconds,
        timebase,
        reference,
    )
end

function AlgorithmGraphs.capture_model_timestamp(
    metadata::PWAO.AcquisitionMetadata,
    origin::PipeWireAcquisitionTimeOrigin,
)
    identity = PWAO.acquisition_identity(metadata)
    identity === nothing && _capture_error(
        "PipeWire model-time capture requires an acquisition identity",
    )
    duration = PWAO.acquisition_exposure_duration(metadata)
    duration === nothing && _capture_error(
        "PipeWire model-time capture requires an exposure duration",
    )
    start, timebase, reference = _acquisition_clock(metadata)
    timebase == origin.timebase || _capture_error(
        "PipeWire acquisition timebase differs from the capture origin",
    )
    isequal(reference, origin.ptp_reference) || _capture_error(
        "PipeWire acquisition PTP reference differs from the capture origin",
    )
    start.nanoseconds >= origin.nanoseconds || _capture_error(
        "PipeWire acquisition timestamp precedes the capture origin",
    )
    start.uncertainty_nanoseconds <= typemax(Int64) || _capture_error(
        "PipeWire acquisition uncertainty exceeds ModelDuration range",
    )

    timestamp = AG.ModelTimestamp(start.nanoseconds - origin.nanoseconds)
    uncertainty = AG.ModelDuration(start.uncertainty_nanoseconds)
    provenance = PipeWireAcquisitionProvenance(
        identity,
        start.nanoseconds,
        timebase,
        reference,
        duration,
    )
    return AG.CapturedModelTimestamp(timestamp, uncertainty, provenance)
end

end # module AdaptiveOpticsSimPipeWireAOExt
