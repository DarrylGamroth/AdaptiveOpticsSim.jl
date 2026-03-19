module AdaptiveOpticsSimCSVExt

using AdaptiveOpticsSim
import CSV

function AdaptiveOpticsSim.write_telemetry_csv(path::AbstractString, telemetry::AdaptiveOpticsSim.Telemetry; kwargs...)
    return CSV.write(path, telemetry; kwargs...)
end

function AdaptiveOpticsSim.write_telemetry_csv(path::AbstractString,
    telemetry::Union{
        AdaptiveOpticsSim.ClosedLoopTrace,
        AdaptiveOpticsSim.GSCClosedLoopTrace,
        AdaptiveOpticsSim.GSCAtmosphereReplayTrace,
    };
    kwargs...)
    return CSV.write(path, telemetry; kwargs...)
end

end
