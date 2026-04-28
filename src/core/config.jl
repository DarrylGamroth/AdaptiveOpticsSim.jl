import TOML

config_value(x::Symbol) = String(x)
config_value(x::AbstractArray) = [config_value(v) for v in x]
config_value(x::Tuple) = [config_value(v) for v in x]
config_value(x::NamedTuple) = Dict(string(k) => config_value(v) for (k, v) in pairs(x))
config_value(x::AbstractDict) = Dict(string(k) => config_value(v) for (k, v) in pairs(x))
config_value(x::AbstractOpticalElement) = config_dict(x)
config_value(x::Nothing) = nothing
config_value(x) = x

function config_dict(x::AbstractDict)
    out = Dict{String, Any}()
    for (k, v) in pairs(x)
        val = config_value(v)
        if val !== nothing
            out[string(k)] = val
        end
    end
    return out
end

function config_dict(x::NamedTuple)
    out = Dict{String, Any}()
    for (k, v) in pairs(x)
        val = config_value(v)
        if val !== nothing
            out[string(k)] = val
        end
    end
    return out
end

function config_dict(x)
    T = typeof(x)
    if isstructtype(T)
        if hasproperty(x, :params)
            return config_dict(getfield(x, :params))
        end
        out = Dict{String, Any}()
        for name in fieldnames(T)
            val = config_value(getfield(x, name))
            if val !== nothing
                out[string(name)] = val
            end
        end
        return out
    end
    return config_value(x)
end

function snapshot_config(; kwargs...)
    return Dict(string(k) => config_dict(v) for (k, v) in pairs(kwargs))
end

function write_config_toml(path::AbstractString, config)
    cfg = config_dict(config)
    open(path, "w") do io
        TOML.print(io, cfg)
    end
    return path
end

function write_config_json(args...; kwargs...)
    throw(InvalidConfiguration("JSON3.jl not available; load JSON3 to enable write_config_json."))
end
