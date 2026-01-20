module AdaptiveOpticsTablesExt

using AdaptiveOptics
import Tables

Tables.istable(::Type{<:AdaptiveOptics.Telemetry}) = true
Tables.rowaccess(::Type{<:AdaptiveOptics.Telemetry}) = true
Tables.rows(t::AdaptiveOptics.Telemetry) = t.rows

end
