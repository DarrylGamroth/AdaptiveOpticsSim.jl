module AdaptiveOpticsSimJSON3Ext

using AdaptiveOpticsSim
import JSON3

function AdaptiveOpticsSim.write_config_json(path::AbstractString, config; kwargs...)
    cfg = AdaptiveOpticsSim.config_dict(config)
    open(path, "w") do io
        JSON3.write(io, cfg; kwargs...)
    end
    return path
end

end
