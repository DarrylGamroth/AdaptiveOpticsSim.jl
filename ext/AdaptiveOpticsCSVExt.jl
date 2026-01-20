module AdaptiveOpticsCSVExt

using AdaptiveOptics
import CSV

function AdaptiveOptics.write_telemetry_csv(path::AbstractString, telemetry::AdaptiveOptics.Telemetry; kwargs...)
    return CSV.write(path, telemetry; kwargs...)
end

end
