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
    :steps, :channels, :output_sampling_m, :output_quantity, :pupil_grid, :pupil_opd_schema, :pupil_amplitude_schema, :output_schema)
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
        for key in (:amplitude, :phase, :opd)
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
    optical_sampling_m(plan)

Return the physical sample pitch in metres on both axes of a Proper optical
train's raw output intensity or complex electric field. A prepared scalar
owner reports its frozen native output pitch, including before the first step.
A mapped train reports its declared fixed pitch or cold-prepared reference
output pitch. Obtain an owner with `prepared_graph_node(graph, node_name)`.
Pixel-integration owners and configurations instead report their declared output
pitch, and require equal pitches on both axes for this scalar interface.
Compiled-path owners report the artifact pitch unless an explicit
`output_sampling_calibration` declares detector coordinates for its existing
samples. That calibration changes coordinates only, without resampling values.

For a scalar train's run-immutable plan, obtained as `node.config` from the
Proper node factory, perform the geometry traversal on the host during cold
composition. This validates the selected output reference and returns its
native pitch without binding execution inputs, outputs or an accelerator.
"""
function optical_sampling_m end

# Model parsing is cold I/O, owned by the optional YAML/FITSIO extension.
function _read_compiled_optical_model(path,hash)
    throw(AlgorithmGraphError("compiled optical model reading requires YAML.jl and FITSIO.jl; load both before preparation"))
end

struct CompiledOpticalPathsNodeFactory{T<:Union{Float32,Float64}} end

"""
    compiled_optical_paths_node_factory(; T=Float64)

Return the cold factory for a hash-bound compiled multi-path Fresnel model.
Proper supplies execution; YAML and FITSIO supply optional model artifact I/O.
Optional `output_sampling_calibration=(sampling_m=..., provenance="...")`
declares the physical pitch of existing output samples in detector coordinates.
It requires a finite positive square-pixel pitch and nonempty application-owned
calibration provenance. It changes no intensities, dimensions, propagation,
radiometry or detector response, and supplies no absolute detector ROI origin.
"""
compiled_optical_paths_node_factory(; T::Type{<:Union{Float32,Float64}}=Float64) =
    CompiledOpticalPathsNodeFactory{T}()

const _COMPILED_PATHS_FIELDS = (:model_path, :model_sha256, :dm_grid, :pupil_grid,
    :dm_opd_schema, :pupil_opd_schema, :pupil_amplitude_schema, :output_schema,
    :summation, :entrance_amplitude_mode, :entrance_reference_amplitude,
    :output_sampling_calibration)

function _compiled_output_sampling_calibration(config, ::Type{T}) where {T}
    hasproperty(config, :output_sampling_calibration) || return nothing
    calibration = config.output_sampling_calibration
    context = "compiled output_sampling_calibration"
    calibration isa NamedTuple || throw(_file_error(context, "expected a table"))
    _require_named_fields(calibration, (:sampling_m, :provenance),
        (:sampling_m, :provenance), context)
    sampling = T(_file_real(calibration.sampling_m, "$context.sampling_m"))
    isfinite(sampling) && sampling > zero(T) || throw(_file_error(
        context, "sampling_m must be finite and positive in the node numeric type"))
    provenance = _file_string(calibration.provenance, "$context.provenance")
    isempty(strip(provenance)) && throw(_file_error(context, "provenance must not be empty"))
    return (sampling_m=sampling, provenance=provenance)
end

function graph_node_bindings(::CompiledOpticalPathsNodeFactory{T}, config::NamedTuple) where {T}
    _require_named_fields(config, _COMPILED_PATHS_FIELDS,
        (:model_path, :model_sha256, :dm_grid, :pupil_grid), "compiled optical paths config")
    _compiled_output_sampling_calibration(config, T)
    mode = hasproperty(config,:entrance_amplitude_mode) ? config.entrance_amplitude_mode : "reference_ratio"
    mode in ("reference_ratio","log_attenuation") || throw(AlgorithmGraphError("invalid entrance_amplitude_mode"))
    if mode == "reference_ratio"
        hasproperty(config,:entrance_reference_amplitude) || throw(AlgorithmGraphError("reference_ratio requires entrance_reference_amplitude binding"))
        return (_file_identifier(config.entrance_reference_amplitude,"entrance reference amplitude binding"),)
    end
    hasproperty(config,:entrance_reference_amplitude) && throw(AlgorithmGraphError("log_attenuation does not use an entrance reference amplitude"))
    return ()
end

function (factory::CompiledOpticalPathsNodeFactory{T})(name::Symbol,config::NamedTuple,
        props::NamedTuple,bindings::NamedTuple=NamedTuple()) where {T}
    return _compiled_optical_paths_node(T,name,config,props,bindings)
end

function _compiled_optical_paths_node(::Type,name,config,props,bindings)
    throw(AlgorithmGraphError("compiled optical paths require Proper.jl"))
end
