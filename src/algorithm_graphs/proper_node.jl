struct ProperOpticalTrainNodeFactory{T<:Union{Float32,Float64}} end

"""
    proper_optical_train_node_factory(; T=Float64)

Return the cold factory for a fixed scalar or polychromatic Proper optical train.
Loading `Proper` activates the optional implementation. The built-in declarative graph
map admits this factory as `proper_optical_train_f64` and
`proper_optical_train_f32`; no executable source is read from configuration.
"""
proper_optical_train_node_factory(; T::Type{<:Union{Float32,Float64}}=Float64) =
    ProperOpticalTrainNodeFactory{T}()

function (factory::ProperOpticalTrainNodeFactory{T})(name::Symbol,
        config::NamedTuple, props::NamedTuple) where {T}
    return _proper_optical_train_node(T, name, config, props, NamedTuple())
end

function (factory::ProperOpticalTrainNodeFactory{T})(name::Symbol,
        config::NamedTuple, props::NamedTuple, bindings::NamedTuple) where {T}
    return _proper_optical_train_node(T, name, config, props, bindings)
end

function _proper_optical_train_node(::Type, name::Symbol,
        config::NamedTuple, props::NamedTuple, bindings::NamedTuple)
    throw(AlgorithmGraphError("Proper optical train requires Proper.jl; load it with `using Proper`"))
end

const _PROPER_CONFIG_FIELDS = (:resolution, :pupil_resolution, :diameter_m, :wavelength_um,
    :steps, :channels, :output_sampling_m, :pupil_grid, :pupil_opd_schema, :pupil_amplitude_schema, :output_schema)
const _PROPER_MAX_SAMPLED_STEPS = 256

function _proper_channel_configs(config::NamedTuple)
    hasproperty(config, :channels) || return ()
    channels = _file_tuple(config.channels, "channels")
    1 <= length(channels) <= 64 ||
        throw(AlgorithmGraphError("Proper optical train requires 1:64 channels"))
    for channel in channels
        channel isa NamedTuple || throw(AlgorithmGraphError("channel must be a mapping"))
        _require_named_fields(channel, (:wavelength_um, :weight, :transmission, :steps),
            (:wavelength_um, :weight, :steps), "Proper optical train channel")
    end
    return channels
end

function _proper_has_runtime_opd(config::NamedTuple)
    hasproperty(config, :steps) || throw(AlgorithmGraphError("Proper optical train config is missing steps"))
    count = 0
    for train in (config, _proper_channel_configs(config)...)
        for step in _file_tuple(train.steps, "steps")
            step isa NamedTuple || throw(AlgorithmGraphError("step must be a mapping"))
            hasproperty(step, :op) || throw(AlgorithmGraphError("optical step is missing op"))
            _file_string(step.op, "optical step op") == "add_opd" || continue
            count += 1
            count <= 64 || throw(AlgorithmGraphError("Proper optical train admits at most 64 runtime OPD inputs"))
        end
    end
    if count > 0 || hasproperty(config, :pupil_grid)
        length(config.steps) <= _PROPER_MAX_SAMPLED_STEPS ||
            throw(AlgorithmGraphError("runtime OPD steps and pupil_grid admit at most $_PROPER_MAX_SAMPLED_STEPS total optical steps"))
        (hasproperty(config, :channels) || hasproperty(config, :output_sampling_m)) &&
            throw(AlgorithmGraphError("runtime OPD steps and pupil_grid require a scalar train without channels or output_sampling_m"))
    end
    return count > 0
end

function _proper_mask_roles(config::NamedTuple)
    hasproperty(config, :steps) || throw(AlgorithmGraphError("Proper optical train config is missing steps"))
    roles = Symbol[]
    trains = (config, _proper_channel_configs(config)...)
    for train in trains, step in _file_tuple(train.steps, "steps")
        step isa NamedTuple || throw(AlgorithmGraphError("step must be a mapping"))
        hasproperty(step, :op) || throw(AlgorithmGraphError("optical step is missing op"))
        _file_string(step.op, "optical step op") == "sampled_mask" || continue
        hasproperty(step, :amplitude) || throw(AlgorithmGraphError("sampled mask is missing amplitude role"))
        for key in (:amplitude, :phase)
            hasproperty(step, key) || continue
            role = _file_identifier(getproperty(step, key), "sampled mask $key role")
            role in _PROPER_CONFIG_FIELDS && throw(AlgorithmGraphError("sampled mask role '$role' collides with optical train configuration"))
            role in roles || push!(roles, role)
        end
    end
    return Tuple(roles)
end

function graph_node_bindings(::ProperOpticalTrainNodeFactory, config::NamedTuple)
    _proper_has_runtime_opd(config)
    roles = _proper_mask_roles(config)
    _require_named_fields(config, (_PROPER_CONFIG_FIELDS..., roles...),
        (:resolution, :diameter_m, :wavelength_um, :steps, roles...), "Proper optical train config")
    names = map(roles) do role
        _file_identifier(getproperty(config, role), "sampled mask binding '$role'")
    end
    return Tuple(unique(names))
end

"""
    optical_sampling_m(owner)

Return the physical sample pitch in metres on both axes of a prepared Proper
optical train's raw output intensity. For a scalar train, read this after stepping
the graph: it reports the retained wavefront's current sampling. CUDA cold
preparation may already trace the output geometry before the first user step. A
mapped train always returns the declared fixed pitch or cold-prepared reference
output sampling.
Obtain the owner with `prepared_graph_node(graph, node_name)`.
"""
function optical_sampling_m end
