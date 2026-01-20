module AdaptiveOpticsJSON3Ext

using AdaptiveOptics
import JSON3

function AdaptiveOptics.write_config_json(path::AbstractString, config; kwargs...)
    cfg = AdaptiveOptics.config_dict(config)
    open(path, "w") do io
        JSON3.write(io, cfg; kwargs...)
    end
    return path
end

end
