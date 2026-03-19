module AdaptiveOpticsSimTablesExt

using AdaptiveOpticsSim
import Tables

Tables.istable(::Type{<:AdaptiveOpticsSim.Telemetry}) = true
Tables.rowaccess(::Type{<:AdaptiveOpticsSim.Telemetry}) = true
Tables.rows(t::AdaptiveOpticsSim.Telemetry) = t.rows

Tables.istable(::Type{<:AdaptiveOpticsSim.ClosedLoopTrace}) = true
Tables.rowaccess(::Type{<:AdaptiveOpticsSim.ClosedLoopTrace}) = true
Tables.rows(t::AdaptiveOpticsSim.ClosedLoopTrace) = t.rows

Tables.istable(::Type{<:AdaptiveOpticsSim.GSCClosedLoopTrace}) = true
Tables.rowaccess(::Type{<:AdaptiveOpticsSim.GSCClosedLoopTrace}) = true
Tables.rows(t::AdaptiveOpticsSim.GSCClosedLoopTrace) = t.rows

Tables.istable(::Type{<:AdaptiveOpticsSim.GSCAtmosphereReplayTrace}) = true
Tables.rowaccess(::Type{<:AdaptiveOpticsSim.GSCAtmosphereReplayTrace}) = true
Tables.rows(t::AdaptiveOpticsSim.GSCAtmosphereReplayTrace) = t.rows

end
