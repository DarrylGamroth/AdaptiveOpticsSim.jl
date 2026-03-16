module AdaptiveOpticsSimTablesExt

using AdaptiveOpticsSim
import Tables

Tables.istable(::Type{<:AdaptiveOpticsSim.Telemetry}) = true
Tables.rowaccess(::Type{<:AdaptiveOpticsSim.Telemetry}) = true
Tables.rows(t::AdaptiveOpticsSim.Telemetry) = t.rows

end
