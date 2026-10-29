module AdaptiveOpticsSimProperExt

using AdaptiveOpticsSim
using Proper
using FixedSizeArrays: FixedSizeMatrixDefault, FixedSizeVectorDefault
import AdaptiveOpticsSim.AlgorithmGraphs as AOG
import AdaptiveOpticsSim.Backends

const DEFAULT_OPD_SCHEMA = "org.adaptiveopticssim.pupil-opd/1"
const DEFAULT_AMPLITUDE_SCHEMA = "org.adaptiveopticssim.pupil-amplitude/1"
const DEFAULT_INTENSITY_SCHEMA = "org.adaptiveopticssim.optical-intensity/1"
const DEFAULT_FIELD_SCHEMA = "org.adaptiveopticssim.electric-field/1"

function _proper_output_quantity(quantity::Symbol)
    quantity in (:intensity, :electric_field) ||
        throw(AOG.AlgorithmGraphError("Proper output_quantity must be intensity or electric_field"))
    return Val(quantity)
end
_proper_output_quantity(quantity::AbstractString) = _proper_output_quantity(Symbol(quantity))
_proper_output_quantity(quantity) =
    throw(AOG.AlgorithmGraphError("Proper output_quantity must be intensity or electric_field"))

_output_schema(::Val{:intensity}) = DEFAULT_INTENSITY_SCHEMA
_output_schema(::Val{:electric_field}) = DEFAULT_FIELD_SCHEMA
_output_type(::Type{T}, ::Val{:intensity}) where {T} = T
_output_type(::Type{T}, ::Val{:electric_field}) where {T} = Complex{T}

_require_output_reference(::Val{:intensity}, wavefront) = nothing
function _require_output_reference(::Val{:electric_field}, wavefront)
    wavefront.reference_surface === Proper.PLANAR || throw(AOG.AlgorithmGraphError(
        "Proper electric_field output requires a planar reference surface; spherical-reference fields need an explicit physical phase conversion"))
    return nothing
end

_finish_optical_train!(wavefront, output, ::Val{:intensity}) =
    Proper.prop_end!(output, wavefront)
_finish_optical_train!(wavefront, output, ::Val{:electric_field}) =
    Proper.prop_end!(output, wavefront; noabs=true)

function _train_step(step::NamedTuple, ::Type{T}, config, bindings) where {T}
    hasproperty(step, :op) || throw(AOG.AlgorithmGraphError("optical step is missing op"))
    op = AOG._file_string(step.op, "optical step op")
    op == "add_opd" && return _runtime_opd_step(step, T)
    fields = if op == "propagate"
        (:distance_m,)
    elseif op == "lens"
        (:focal_length_m,)
    elseif op in ("circular_aperture", "circular_obscuration")
        (:radius_m,)
    elseif op == "multiply"
        (:amplitude,)
    elseif op == "circular_mask"
        (:circles, :mode)
    elseif op == "sampled_mask"
        (:amplitude, :sampling_m, :center_pixel, :outside_amplitude, :outside_phase_rad)
    else
        throw(AOG.AlgorithmGraphError("unsupported optical operation '$op'"))
    end
    allowed = op == "sampled_mask" ? (:op, :name, fields..., :phase) : (:op, :name, fields...)
    required = op == "sampled_mask" ? (:name, fields...) : fields
    AOG._require_named_fields(step, allowed, required, "optical step")
    name = hasproperty(step, :name) ? AOG._file_string(step.name, "optical step name") : op
    if op == "sampled_mask"
        amplitude_role = AOG._file_identifier(step.amplitude, "sampled mask amplitude role")
        amplitude = _mask_snapshot(AOG._require_binding(bindings,
            getproperty(config, amplitude_role), "sampled mask amplitude"), T, :amplitude)
        phase = if hasproperty(step, :phase)
            phase_role = AOG._file_identifier(step.phase, "sampled mask phase role")
            _mask_snapshot(AOG._require_binding(bindings,
                getproperty(config, phase_role), "sampled mask phase"), T, :phase)
        else
            nothing
        end
        phase === nothing || size(phase) == size(amplitude) ||
            throw(AOG.AlgorithmGraphError("sampled mask amplitude and phase shapes must match"))
        sampling_m = _mask_pair(step.sampling_m, T, "sampling_m")
        all(>(zero(T)), sampling_m) || throw(AOG.AlgorithmGraphError("sampled mask sampling_m must be positive"))
        center_pixel = _mask_pair(step.center_pixel, T, "center_pixel")
        outside_amplitude = T(AOG._file_real(step.outside_amplitude, "outside_amplitude"))
        0 <= step.outside_amplitude <= 1 && isfinite(outside_amplitude) ||
            throw(AOG.AlgorithmGraphError("outside_amplitude must be in [0,1]"))
        outside_phase_rad = T(AOG._file_real(step.outside_phase_rad, "outside_phase_rad"))
        isfinite(outside_phase_rad) || throw(AOG.AlgorithmGraphError("outside_phase_rad is not representable"))
        return (; op=:sampled_mask, name, amplitude, phase, sampling_m, center_pixel,
            outside=outside_amplitude*cis(outside_phase_rad))
    end
    if op == "circular_mask"
        mode = AOG._file_string(step.mode, "circular mask mode")
        mode in ("transmit", "reflect") || throw(AOG.AlgorithmGraphError("invalid circular mask mode"))
        step.circles isa Tuple || throw(AOG.AlgorithmGraphError("circles must be a normalized tuple"))
        isempty(step.circles) && throw(AOG.AlgorithmGraphError("circles must not be empty"))
        circles = map(step.circles) do circle
            circle isa NamedTuple || throw(AOG.AlgorithmGraphError("circle must be a mapping"))
            AOG._require_named_fields(circle, (:radius_m, :x_m, :y_m), (:radius_m, :x_m, :y_m), "circle")
            values = map(k -> T(AOG._file_real(getproperty(circle, k), "circle $k")), (:radius_m, :x_m, :y_m))
            all(isfinite, values) && values[1] > 0 || throw(AOG.AlgorithmGraphError("invalid circle dimensions"))
            NamedTuple{(:radius_m, :x_m, :y_m)}(values)
        end
        return (; op=Symbol(op), name, circles, mode=Symbol(mode))
    end
    key = only(fields)
    value = T(AOG._file_real(getproperty(step, key), "optical step $key"))
    isfinite(value) || throw(AOG.AlgorithmGraphError("optical step $key is not representable"))
    if key == :radius_m
        value > 0 || throw(AOG.AlgorithmGraphError("radius_m must be positive"))
    elseif key == :focal_length_m
        !iszero(value) || throw(AOG.AlgorithmGraphError("focal_length_m must be nonzero"))
    elseif key == :amplitude
        0 <= value <= 1 || throw(AOG.AlgorithmGraphError("amplitude must be in [0,1]"))
    end
    return merge((; op=Symbol(op), name), NamedTuple{(key,)}((value,)))
end

function _mask_pair(value, ::Type{T}, name) where {T}
    pair = AOG._file_tuple(value, "sampled mask $name")
    length(pair) == 2 || throw(AOG.AlgorithmGraphError("sampled mask $name must contain two values"))
    converted = map(x -> T(AOG._file_real(x, "sampled mask $name")), pair)
    all(isfinite, converted) || throw(AOG.AlgorithmGraphError("sampled mask $name is not representable"))
    return converted
end

function _mask_snapshot(array, ::Type{T}, quantity) where {T}
    array isa AbstractMatrix && eltype(array) <: Real ||
        throw(AOG.AlgorithmGraphError("sampled mask $quantity must be a real matrix"))
    !Base.has_offset_axes(array) ||
        throw(AOG.AlgorithmGraphError("sampled mask $quantity requires one-based axes"))
    all(>(0), size(array)) || throw(AOG.AlgorithmGraphError("sampled mask $quantity must not be empty"))
    host = Backends.host_array(array)
    all(isfinite, host) || throw(AOG.AlgorithmGraphError("sampled mask $quantity must be finite"))
    quantity === :amplitude && !all(x -> 0 <= x <= 1, host) &&
        throw(AOG.AlgorithmGraphError("sampled mask amplitude must be in [0,1]"))
    converted = T.(host)
    all(isfinite, converted) || throw(AOG.AlgorithmGraphError("sampled mask $quantity is not representable"))
    return FixedSizeMatrixDefault{T}(converted)
end

function _train_optic!(wf, step, mask, context)
    if step.op === :propagate
        Proper.prop_propagate(wf, step.distance_m, context, step.name)
    elseif step.op === :lens
        Proper.prop_lens(wf, step.focal_length_m, context, step.name)
    elseif step.op === :circular_aperture
        Proper.prop_circular_aperture(wf, step.radius_m)
    elseif step.op === :circular_obscuration
        Proper.prop_circular_obscuration(wf, step.radius_m)
    elseif step.op === :multiply
        Proper.prop_multiply(wf, step.amplitude)
    elseif step.op === :circular_mask || step.op === :sampled_mask
        Proper.prop_multiply(wf, mask)
    end
    return nothing
end

function _train_optic!(wf, step, mask, context, shifted_mask)
    if step.op === :sampled_mask
        # Proper's planned FFT publishes its complex scratch as wf.field.
        # Its allocating centered-map fallback would then be required for a
        # complex mask. Retain a separate shift destination so neither the live
        # field nor the immutable centered mask is overwritten during capture.
        Proper.prop_shift_center!(shifted_mask, mask; inverse=true)
        wf.field .*= shifted_mask
    else
        _train_optic!(wf, step, mask, context)
    end
    return nothing
end


struct OpticalTrainNode{T} end

struct OpticalTrainPlan{T,S,G,Q}
    steps::S
    resolution::Int
    pupil_resolution::Int
    diameter_m::T
    wavelength_um::T
    schemas::NTuple{3,String}
    pupil_grid::G
    output_quantity::Q
end

function AOG._proper_optical_train_node(::Type{T}, name::Symbol,
        config::NamedTuple, props::NamedTuple, bindings::NamedTuple) where {T<:Union{Float32,Float64}}
    AOG.graph_node_bindings(AOG.ProperOpticalTrainNodeFactory{T}(), config)
    AOG._require_named_fields(props, (), (), "Proper optical train props")
    steps = map(AOG._file_tuple(config.steps, "steps")) do step
        step isa NamedTuple || throw(AOG.AlgorithmGraphError("step must be a mapping"))
        _train_step(step, T, config, bindings)
    end
    _validate_runtime_opd_inputs(steps)
    n = AOG._file_integer(config.resolution, "resolution")
    iseven(n) || throw(AOG.AlgorithmGraphError("Proper optical train resolution must be even for centered propagation"))
    p = hasproperty(config, :pupil_resolution) ?
        AOG._file_integer(config.pupil_resolution, "pupil_resolution") : n
    0 < p <= n || throw(AOG.AlgorithmGraphError("pupil_resolution must be in 1:resolution"))
    diameter = T(AOG._file_real(config.diameter_m, "diameter_m"))
    wavelength = T(AOG._file_real(config.wavelength_um, "wavelength_um"))
    all(x -> isfinite(x) && x > 0, (diameter, wavelength)) ||
        throw(AOG.AlgorithmGraphError("diameter and wavelength must be finite and positive"))
    output_quantity = hasproperty(config, :output_quantity) ?
        _proper_output_quantity(config.output_quantity) : Val(:intensity)
    if output_quantity === Val(:electric_field) &&
            (hasproperty(config, :channels) || hasproperty(config, :output_sampling_m))
        throw(AOG.AlgorithmGraphError(
            "Proper electric_field output requires a scalar train without channels or output_sampling_m"))
    end
    schemas = map((:pupil_opd_schema, :pupil_amplitude_schema, :output_schema),
        (DEFAULT_OPD_SCHEMA, DEFAULT_AMPLITUDE_SCHEMA, _output_schema(output_quantity))) do key, default
        value = hasproperty(config, key) ? AOG._file_string(getproperty(config, key), string(key)) : default
        isempty(value) && throw(AOG.AlgorithmGraphError("$key must not be empty"))
        value
    end
    pupil_grid = hasproperty(config, :pupil_grid) ? _entrance_pupil_grid(config.pupil_grid, T, p) : nothing
    plan = OpticalTrainPlan(steps, n, p, diameter, wavelength, schemas, pupil_grid, output_quantity)
    if hasproperty(config, :channels) || hasproperty(config, :output_sampling_m)
        return _polychromatic_optical_train_node(T, name, plan, config, bindings)
    end
    return AOG.algorithm_node(name, OpticalTrainNode{T}, plan; props=NamedTuple())
end

function AOG.graph_node_ports(::Type{OpticalTrainNode{T}}, plan::OpticalTrainPlan) where {T}
    pupil = (plan.pupil_resolution, plan.pupil_resolution)
    output = (plan.resolution, plan.resolution)
    ports = (AOG.graph_port_contract(:pupil_opd, :input, :data, T, pupil, plan.schemas[1], :column_major),
        AOG.graph_port_contract(:pupil_amplitude, :input, :data, T, pupil, plan.schemas[2], :column_major),
        AOG.graph_port_contract(:output, :output, :data, _output_type(T, plan.output_quantity), output, plan.schemas[3], :column_major))
    return (ports..., _runtime_opd_ports(T, plan.steps)...)
end

struct OpticalTrainWorkspace{W,C,A,M,S}
    wavefront::W
    context::C
    padded_opd::A
    padded_amplitude::A
    masks::M
    shifted_mask::S
end

struct OpticalTrainOwner{P,W,D,A,O,T}
    plan::P
    workspace::W
    pupil_opd::D
    pupil_amplitude::A
    output::O
    sampling_m::T
end

function _mask_sample(step, row, col)
    if 1 <= row <= size(step.amplitude, 1) && 1 <= col <= size(step.amplitude, 2)
        amplitude = step.amplitude[row, col]
        return step.phase === nothing ? complex(amplitude) : amplitude*cis(step.phase[row, col])
    end
    return step.outside
end

function _resample_mask(step, n, sampling::T) where {T}
    isfinite(sampling) && sampling > 0 ||
        throw(AOG.AlgorithmGraphError("sampled mask has invalid plane sampling"))
    mask = Matrix{Complex{T}}(undef, n, n)
    center = fld(n, 2)+1
    for col in 1:n, row in 1:n
        source_row = (row-center)*sampling/step.sampling_m[1]+step.center_pixel[1]
        source_col = (col-center)*sampling/step.sampling_m[2]+step.center_pixel[2]
        if !(0 <= source_row <= size(step.amplitude, 1)+1 &&
                0 <= source_col <= size(step.amplitude, 2)+1)
            mask[row, col] = step.outside
            continue
        end
        i, j = floor(Int, source_row), floor(Int, source_col)
        dy, dx = source_row-i, source_col-j
        # Interpolate the phasors, including each outside neighbor, to preserve
        # continuity across phase wraps rather than interpolating phase angles.
        mask[row, col] = (1-dy)*((1-dx)*_mask_sample(step, i, j)+dx*_mask_sample(step, i, j+1))+
            dy*((1-dx)*_mask_sample(step, i+1, j)+dx*_mask_sample(step, i+1, j+1))
    end
    return mask
end

function _prepare_masks(plan::OpticalTrainPlan{T}, target) where {T}
    n = plan.resolution
    context = Proper.RunContext(Matrix{Complex{T}})
    wf = Proper.prop_begin(plan.diameter_m, plan.wavelength_um*T(1e-6), n;
        beam_diam_fraction=T(plan.pupil_resolution/n), context)
    masks = map(plan.steps) do step
        mask = if step.op === :circular_mask
            sampling = Proper.prop_get_sampling(wf)
            isfinite(sampling) && sampling > 0 ||
                throw(AOG.AlgorithmGraphError("circular mask has invalid plane sampling"))
            center = fld(n,2)+1
            for (index, circle) in pairs(step.circles)
                j = round(Int, clamp(circle.x_m/sampling+center, one(T), T(n)))
                i = round(Int, clamp(circle.y_m/sampling+center, one(T), T(n)))
                hypot((j-center)*sampling-circle.x_m, (i-center)*sampling-circle.y_m) <= circle.radius_m ||
                    throw(AOG.AlgorithmGraphError(
                        "circular mask '$(step.name)' circle $index has no represented sample on the $(n)x$(n) grid at sampling $(sampling) m; increase sampling support or revise the declared ideal mask"))
            end
            centered = Matrix{T}(undef, n, n)
            for j in 1:n, i in 1:n
                x = (j-fld(n,2)-1)*sampling
                y = (i-fld(n,2)-1)*sampling
                inside = any(c -> hypot(x-c.x_m,y-c.y_m) <= c.radius_m, step.circles)
                centered[i,j] = T(step.mode === :transmit ? inside : !inside)
            end
            centered
        elseif step.op === :sampled_mask
            _resample_mask(step, n, T(Proper.prop_get_sampling(wf)))
        else
            nothing
        end
        _train_optic!(wf, step, mask, context)
        if mask === nothing
            nothing
        else
            device_mask = Backends.allocate_device_array(target, eltype(mask), n, n)
            copyto!(device_mask, mask)
            device_mask
        end
    end
    _require_output_reference(plan.output_quantity, wf)
    return (; masks, sampling_m=T(Proper.prop_get_sampling(wf)))
end

function AOG.prepare_graph_node(::Type{OpticalTrainNode{T}}, plan::OpticalTrainPlan{T},
        ::NamedTuple{()}, inputs::NamedTuple,
        outputs::NamedTuple{(:output,)}, ::NamedTuple{()}, target) where {T}
    (plan.pupil_grid !== nothing || any(step -> step.op === :add_opd, plan.steps)) &&
        return _prepare_runtime_opd_train(plan, inputs, outputs, target)
    return _prepare_optical_train(plan, inputs, outputs, target).owner
end

function _prepare_optical_train(plan::OpticalTrainPlan{T}, inputs, outputs, target) where {T}
    n = plan.resolution
    for (array, shape, numeric_type) in
            ((inputs.pupil_opd, (plan.pupil_resolution, plan.pupil_resolution), T),
             (inputs.pupil_amplitude, (plan.pupil_resolution, plan.pupil_resolution), T),
             (outputs.output, (n,n), _output_type(T, plan.output_quantity)))
        eltype(array) === numeric_type && size(array) == shape && Backends.compute_device(array) == target ||
            throw(AOG.AlgorithmGraphError("Proper optical train buffer contract mismatch"))
    end
    field = Backends.allocate_device_array(target, Complex{T}, n, n)
    context = Proper.RunContext(typeof(field))
    wavefront = Proper.prop_begin!(field, plan.diameter_m, plan.wavelength_um*T(1e-6);
        beam_diam_fraction=T(plan.pupil_resolution/n), context)
    padded_opd = Backends.allocate_device_array(target, T, n, n)
    padded_amplitude = similar(padded_opd)
    fill!(padded_opd, 0)
    fill!(padded_amplitude, 0)
    geometry = _prepare_masks(plan, target)
    shifted_mask = any(step -> step.op === :sampled_mask, plan.steps) ?
        Backends.allocate_device_array(target, Complex{T}, n, n) : nothing
    workspace = OpticalTrainWorkspace(wavefront, context, padded_opd, padded_amplitude,
        geometry.masks, shifted_mask)
    owner = OpticalTrainOwner(plan, workspace, inputs.pupil_opd, inputs.pupil_amplitude, outputs.output, geometry.sampling_m)
    _prepare_train_capture!(owner, target)
    return (; owner, geometry.sampling_m)
end

_prepare_train_capture!(owner, target) = nothing

function _prepare_train_capture!(owner,
        ::Backends.AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}})
    # All optical geometry depends only on the immutable plan. Exercise the
    # complete trajectory once while cold: this creates FFT plans, fills fixed
    # workspaces and settles the live field's alias to the FFT scratch before
    # any enclosing application capture records addresses. Reset clears only
    # scientific buffers and leaves those plans and addresses intact.
    AOG.step_graph_node!(owner)
    AOG.reset_graph_node!(owner)
    return nothing
end

function _begin_optical_train!(owner::OpticalTrainOwner)
    plan, workspace = owner.plan, owner.workspace
    opd, amplitude = workspace.padded_opd, workspace.padded_amplitude
    fill!(opd, 0)
    fill!(amplitude, 0)
    offset = fld(plan.resolution,2)-fld(plan.pupil_resolution,2)
    indices = (1:plan.pupil_resolution) .+ offset
    Backends.copyto_backend_async!(view(opd, indices, indices), owner.pupil_opd)
    Backends.copyto_backend_async!(view(amplitude, indices, indices), owner.pupil_amplitude)
    wf = workspace.wavefront
    Proper.prop_begin!(wf, plan.diameter_m,
        plan.wavelength_um*oftype(plan.wavelength_um, 1e-6);
        beam_diam_fraction=plan.pupil_resolution/plan.resolution)
    Proper.prop_multiply(wf, amplitude)
    Proper.prop_add_phase(wf, opd)
    return nothing
end

function AOG.step_graph_node!(owner::OpticalTrainOwner)
    _begin_optical_train!(owner)
    plan, workspace = owner.plan, owner.workspace
    wf = workspace.wavefront
    foreach(plan.steps, workspace.masks) do step, mask
        _train_optic!(wf, step, mask, workspace.context, workspace.shifted_mask)
    end
    _finish_optical_train!(wf, owner.output, plan.output_quantity)
    return nothing
end

AOG.optical_sampling_m(plan::OpticalTrainPlan) =
    _prepare_masks(plan, Backends.HostComputeDevice()).sampling_m

AOG.optical_sampling_m(owner::OpticalTrainOwner) = owner.sampling_m
_optical_train_capture_capability(::Backends.AbstractComputeDevice) = AOG.GraphNodeCaptureUnsupported()
_optical_train_capture_capability(::Backends.AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}}) = AOG.GraphNodeCaptureSafe()
AOG.graph_node_capture_capability(owner::OpticalTrainOwner) =
    _optical_train_capture_capability(Backends.compute_device(owner.output))

function AOG.reset_graph_node!(owner::OpticalTrainOwner)
    fill!(owner.workspace.padded_opd, 0)
    fill!(owner.workspace.padded_amplitude, 0)
    fill!(owner.output, 0)
    return nothing
end

include("proper_opd.jl")
include("proper_spectrum.jl")

end
