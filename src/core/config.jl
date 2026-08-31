config_value(x::Symbol) = String(x)
function config_value(x::AbstractArray)
    normalized_type = Union{}
    for value in x
        normalized_type = Union{
            normalized_type,
            typeof(config_value(value)),
        }
    end
    normalized = Array{normalized_type}(undef, size(x))
    index = 1
    for value in x
        normalized[index] = config_value(value)
        index += 1
    end
    return normalized
end

config_value(x::Tuple) = _concrete_union_vector(map(config_value, x))
config_value(x::NamedTuple) = config_dict(x)
config_value(x::AbstractDict) = config_dict(x)
config_value(x::Optics.AbstractOpticalElement) = config_dict(x)
config_value(x::Nothing) = nothing
config_value(x) = x

@inline function _config_entry(key, value)
    normalized = config_value(value)
    normalized === nothing && return nothing
    return string(key) => normalized
end

function _config_entries(values)
    candidates = (
        _config_entry(key, value)
        for (key, value) in pairs(values)
    )
    return Tuple(Iterators.filter(
        entry -> entry !== nothing,
        candidates,
    ))
end

function _concrete_config_dict(entries::Tuple)
    Value = Union{map(entry -> typeof(last(entry)), entries)...}
    out = Dict{String,Value}()
    for entry in entries
        out[first(entry)] = last(entry)
    end
    return out
end

config_dict(x::AbstractDict) = _concrete_config_dict(_config_entries(x))

function config_dict(x::NamedTuple)
    return _concrete_config_dict(_config_entries(x))
end

function config_dict(x)
    T = typeof(x)
    if isstructtype(T)
        if hasproperty(x, :params)
            return config_dict(getfield(x, :params))
        end
        candidates = (
            _config_entry(name, getfield(x, name))
            for name in fieldnames(T)
        )
        entries = Tuple(Iterators.filter(
            entry -> entry !== nothing,
            candidates,
        ))
        return _concrete_config_dict(entries)
    end
    return config_value(x)
end

function snapshot_config(; kwargs...)
    return config_dict((; kwargs...))
end

function write_config_json(args...; kwargs...)
    throw(InvalidConfiguration("JSON3.jl not available; load JSON3 to enable write_config_json."))
end
