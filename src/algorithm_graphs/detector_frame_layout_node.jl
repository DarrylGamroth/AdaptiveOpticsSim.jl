# Detector placement and affine conversion for configured AOS graphs.
# Boundary semantics follow the former GPI2Sim FrameLayoutNode.
# Dimensions are [rows, columns]; the detector origin counts from zero.

const _MAX_DETECTOR_FRAME_LAYOUT_BYTES = 512 * 1024 * 1024

struct DetectorFrameLayoutNode end
struct DetectorFrameLayoutConfig{T,R,Q}
    input_size::Tuple{Int,Int}
    output_size::Tuple{Int,Int}
    origin::Tuple{Int,Int}
    scale::Float32
    bias::Float32
    range::R
    rounding::Q
    input_schema::String
    output_schema::String
end
struct DetectorFrameLayoutOwner{I,S,O,R,Q}
    detector_frame::I
    science::S
    frame::O
    scale::Float32
    bias::Float32
    range::R
    rounding::Q
end

function _detector_layout_integer(value::Integer, context)
    typemin(Int) <= value <= typemax(Int) ||
        throw(AlgorithmGraphError("$context exceeds machine integer bounds"))
    return Int(value)
end
_detector_layout_integer(value, context) = throw(AlgorithmGraphError("$context must be an integer"))
_detector_layout_integer(value::Bool, context) =
    throw(AlgorithmGraphError("$context must be an integer"))

function _detector_layout_string(value, context)
    text = _file_string(value, context)
    isempty(text) && throw(AlgorithmGraphError("$context must not be empty"))
    return text
end

function _detector_layout_tuple(value, rank, context)
    values = _file_tuple(value, context)
    length(values) == rank || throw(AlgorithmGraphError("$context must have $rank entries"))
    return values
end

function _detector_layout_shape(value, rank, context)
    dimensions = map(_detector_layout_tuple(value, rank, context)) do dimension
        count = _detector_layout_integer(dimension, context)
        count > 0 || throw(AlgorithmGraphError("$context dimensions must be positive"))
        return count
    end
    # Test division before multiplication, including the Float32 element width.
    bytes = sizeof(Float32)
    for count in dimensions
        count <= div(_MAX_DETECTOR_FRAME_LAYOUT_BYTES, bytes) ||
            throw(AlgorithmGraphError("$context exceeds the 512 MiB array limit"))
        bytes *= count
    end
    return dimensions
end

function _detector_layout_f32(value::Real, context)
    isfinite(value) && abs(value) <= floatmax(Float32) ||
        throw(AlgorithmGraphError("$context must be finite and within Float32 range"))
    return Float32(value)
end
_detector_layout_f32(value, context) = throw(AlgorithmGraphError("$context must be a real number"))
_detector_layout_f32(value::Bool, context) = throw(AlgorithmGraphError("$context must be a real number"))

function _detector_layout_range(value::Tuple, context)
    _detector_layout_tuple(value, 2, context)
    lower, upper = map(bound -> _detector_layout_f32(bound, context), value)
    value[1] <= value[2] ||
        throw(AlgorithmGraphError("$context lower bound exceeds upper bound"))
    return (lower, upper)
end
_detector_layout_range(::Nothing, context) = nothing
_detector_layout_range(value, context) =
    throw(AlgorithmGraphError("$context must contain two real bounds"))

function _detector_layout_rounding(value::AbstractString, context)
    value == "none" && return nothing
    value == "nearest_even" && return RoundNearest
    throw(AlgorithmGraphError("$context must be 'none' or 'nearest_even'"))
end
_detector_layout_rounding(value, context) =
    throw(AlgorithmGraphError("$context must be 'none' or 'nearest_even'"))

function _detector_layout_rounding_range(range::Tuple, ::typeof(RoundNearest), context)
    all(isinteger, range) || throw(
        AlgorithmGraphError("$context requires integral bounds with nearest_even rounding")
    )
    return nothing
end
_detector_layout_rounding_range(range, rounding, context) = nothing

_detector_layout_output_type(value::AbstractString) = _detector_layout_output_type(Val(Symbol(value)))
_detector_layout_output_type(::Val{:float32}) = Float32
_detector_layout_output_type(::Val{:uint16}) = UInt16
_detector_layout_output_type(value) =
    throw(AlgorithmGraphError("output_element_type must be 'float32' or 'uint16'"))
_detector_layout_output_contract(::Type{Float32}, range, rounding) = nothing
function _detector_layout_output_contract(::Type{UInt16}, range::Tuple, ::typeof(RoundNearest))
    0 <= range[1] <= range[2] <= typemax(UInt16) ||
        throw(AlgorithmGraphError("UInt16 camera output requires range within 0:65535"))
    return nothing
end
_detector_layout_output_contract(::Type{UInt16}, range, rounding) = throw(
    AlgorithmGraphError(
        "UInt16 camera output requires a bounded range and nearest_even rounding"
    ),
)

"""
Built-in detector frame layout graph factory:

- `detector_frame_layout_f32`: `input_size` and `output_size` in [rows, columns],
  zero-based `origin`, `scale`, `bias`, `input_schema`, and `output_schema`.
  Input `detector_frame` is placed in output `frame` with zero padding. Omit
  `range` for `input * scale + bias`; provide `[lower, upper]` to clamp it.
  Optional `rounding` is `"none"` (the default) or `"nearest_even"` to round
  the converted, clipped value to the nearest integer with ties to even while
  retaining Float32 storage by default. Set `output_element_type: uint16` for
  a typed ADC output; it requires nearest-even rounding and an integral range
  within 0:65535. Rounded ranges require integral bounds.

The node resets its output buffer to zero and requires one-based Float32
input arrays and Float32 or UInt16 output arrays on the declared compute
device. CUDA and AMDGPU owners support fixed-address command-graph capture;
host owners retain stream execution. Both declared shapes retain the 512 MiB
admission limit, calculated using Float32 element width.
"""

function _detector_frame_layout_f32_node(
    name::Symbol, config::NamedTuple, props::NamedTuple
)
    fields = (
        :input_size,
        :output_size,
        :origin,
        :scale,
        :bias,
        :range,
        :rounding,
        :input_schema,
        :output_schema,
        :output_element_type,
    )
    required = (
        :input_size, :output_size, :origin, :scale, :bias, :input_schema, :output_schema
    )
    context = "node '$name' config"
    _require_named_fields(config, fields, required, context)
    _require_named_fields(props, (), (), "node '$name' props")
    input_size = _detector_layout_shape(config.input_size, 2, "$context.input_size")
    output_size = _detector_layout_shape(config.output_size, 2, "$context.output_size")
    origin = map(_detector_layout_tuple(config.origin, 2, "$context.origin")) do index
        _detector_layout_integer(index, "$context.origin")
    end
    all(
        d ->
            0 <= origin[d] <= output_size[d] && input_size[d] <= output_size[d] - origin[d],
        1:2,
    ) || throw(AlgorithmGraphError("input rectangle does not fit inside the output at origin"))
    scale = _detector_layout_f32(config.scale, "$context.scale")
    bias = _detector_layout_f32(config.bias, "$context.bias")
    range = _detector_layout_range(get(config, :range, nothing), "$context.range")
    rounding = _detector_layout_rounding(get(config, :rounding, "none"), "$context.rounding")
    _detector_layout_rounding_range(get(config, :range, nothing), rounding, "$context.range")
    T = _detector_layout_output_type(get(config, :output_element_type, "float32"))
    _detector_layout_output_contract(T, range, rounding)
    input_schema = _detector_layout_string(config.input_schema, "$context.input_schema")
    output_schema = _detector_layout_string(config.output_schema, "$context.output_schema")
    return algorithm_node(
        name,
        DetectorFrameLayoutNode,
        DetectorFrameLayoutConfig{T,typeof(range),typeof(rounding)}(
            input_size,
            output_size,
            origin,
            scale,
            bias,
            range,
            rounding,
            input_schema,
            output_schema,
        ),
    )
end

function graph_node_ports(
    ::Type{DetectorFrameLayoutNode}, config::DetectorFrameLayoutConfig{T}
) where {T}
    return (
        graph_port_contract(
            :detector_frame,
            :input,
            :data,
            Float32,
            config.input_size,
            config.input_schema,
            :column_major,
        ),
        graph_port_contract(
            :frame,
            :output,
            :data,
            T,
            config.output_size,
            config.output_schema,
            :column_major,
        ),
    )
end

function _detector_layout_buffer(
    values::AbstractArray{T,N}, shape::NTuple{N,Int}, ::Type{T} = Float32
) where {T,N}
    Base.require_one_based_indexing(values)
    size(values) == shape ||
        throw(DimensionMismatch("graph adapter buffer shape differs from declaration"))
    return values
end
_detector_layout_buffer(values, shape, type = Float32) =
    throw(AlgorithmGraphError("graph adapter requires $type arrays of the declared rank"))

function prepare_graph_node(
    ::Type{DetectorFrameLayoutNode},
    config::DetectorFrameLayoutConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:detector_frame,)},
    outputs::NamedTuple{(:frame,)},
    ::NamedTuple{()},
    target::Backends.AbstractComputeDevice,
) where {T}
    detector_frame = _detector_layout_buffer(inputs.detector_frame, config.input_size)
    frame = _detector_layout_buffer(outputs.frame, config.output_size, T)
    Backends.compute_device(detector_frame) == target &&
    Backends.compute_device(frame) == target ||
        throw(AlgorithmGraphError("graph adapter buffers must use the declared compute device"))
    Base.mightalias(frame, detector_frame) &&
        throw(AlgorithmGraphError("frame output must not alias its detector input"))
    rows = (config.origin[1]+1):(config.origin[1]+config.input_size[1])
    columns = (config.origin[2]+1):(config.origin[2]+config.input_size[2])
    return DetectorFrameLayoutOwner(
        detector_frame,
        view(frame, rows, columns),
        frame,
        config.scale,
        config.bias,
        config.range,
        config.rounding,
    )
end

# Device fill and fused broadcast retain input/output addresses and enqueue on
# the graph's prepared context, without host scientific-state updates or waits.
_detector_layout_capture_capability(::Backends.AbstractComputeDevice) =
    GraphNodeCaptureUnsupported()
_detector_layout_capture_capability(
    ::Backends.AcceleratorComputeDevice{
        <:Union{Backends.CUDABackend,Backends.AMDGPUBackend}
    },
) = GraphNodeCaptureSafe()
graph_node_capture_capability(owner::DetectorFrameLayoutOwner) =
    _detector_layout_capture_capability(Backends.compute_device(owner.frame))

_layout_value(value, scale, bias, ::Nothing) = value * scale + bias
_layout_value(value, scale, bias, range::Tuple{Float32,Float32}) =
    clamp(value * scale + bias, range[1], range[2])
_layout_round(value, ::Nothing) = value
_layout_round(value, rounding::typeof(RoundNearest)) = round(value, rounding)
_layout_value(value, scale, bias, range, rounding) =
    _layout_round(_layout_value(value, scale, bias, range), rounding)

function step_graph_node!(owner::DetectorFrameLayoutOwner)
    fill!(owner.frame, zero(eltype(owner.frame)))
    owner.science .= _layout_value.(
        owner.detector_frame, owner.scale, owner.bias, Ref(owner.range), Ref(owner.rounding)
    )
    return nothing
end
function reset_graph_node!(owner::DetectorFrameLayoutOwner)
    fill!(owner.frame, zero(eltype(owner.frame)))
    return nothing
end
