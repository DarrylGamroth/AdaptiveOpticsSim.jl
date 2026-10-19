using FixedSizeArrays: FixedSizeMatrixDefault
import ..Backends
import ..Detectors

struct PixelMTFNode{T<:AbstractFloat} end

"""Cold graph contract and owned calibration snapshot for a pixel MTF filter."""
struct PixelMTFNodeConfig{A}
    input_shape::Tuple{Int,Int}
    padding::Int
    input_schema::String
    output_schema::String
    transfer::A
end

"""
    pixel_mtf_node(name; input_shape, padding, T, input_schema, output_schema,
                   transfer)

Declare a deterministic filter around `Detectors.prepare_pixel_mtf`. Input and
output ports are `:input` and `:output`, with the same shape and precision.
Caller-supplied schemas declare their scientific quantity and pixel coordinates.
The real corner-origin transfer has unit DC and the padded grid shape; it is
copied into bounded host configuration storage, then copied to the graph's exact
compute device during preparation. No detector response or unit defaults are
inferred, and the filter does not clamp its output.
"""
function pixel_mtf_node(name::Symbol; input_shape, padding, T,
        input_schema::AbstractString, output_schema::AbstractString,
        transfer::AbstractMatrix)
    T in (Float32, Float64) || throw(AlgorithmGraphError(
        "pixel MTF graph precision must be Float32 or Float64"))
    shape = Detectors._pixel_mtf_shape(input_shape)
    padding isa Integer && !(padding isa Bool) &&
        0 <= padding <= typemax(Int) ÷ 2 || throw(AlgorithmGraphError(
        "pixel MTF graph padding must be a nonnegative integer"))
    pad = Int(padding)
    all(n -> n <= typemax(Int) - 2pad, shape) || throw(AlgorithmGraphError(
        "pixel MTF graph padded dimensions overflow Int"))
    isempty(input_schema) && throw(AlgorithmGraphError(
        "pixel MTF graph input_schema must not be empty"))
    isempty(output_schema) && throw(AlgorithmGraphError(
        "pixel MTF graph output_schema must not be empty"))
    eltype(transfer) <: Real || throw(AlgorithmGraphError(
        "pixel MTF graph transfer must have real elements"))
    Base.require_one_based_indexing(transfer)
    size(transfer) == shape .+ 2pad || throw(AlgorithmGraphError(
        "pixel MTF graph transfer must match the padded input shape"))
    host_transfer = T.(Backends.host_array(transfer))
    all(isfinite, host_transfer) || throw(AlgorithmGraphError(
        "pixel MTF graph transfer must remain finite at graph precision"))
    isapprox(host_transfer[1, 1], one(T); rtol=8eps(T), atol=zero(T)) ||
        throw(AlgorithmGraphError("pixel MTF graph transfer must have unit DC"))
    config = PixelMTFNodeConfig(shape, pad, String(input_schema),
        String(output_schema), FixedSizeMatrixDefault{T}(host_transfer))
    return algorithm_node(name, PixelMTFNode{T}, config; props=NamedTuple())
end

function graph_node_ports(::Type{PixelMTFNode{T}}, config::PixelMTFNodeConfig) where {T}
    return (graph_port_contract(:input, :input, :data, T, config.input_shape,
            config.input_schema, :column_major),
        graph_port_contract(:output, :output, :data, T, config.input_shape,
            config.output_schema, :column_major))
end

struct _PixelMTFOwner{P}
    prepared::P
end

function prepare_graph_node(::Type{PixelMTFNode{T}}, config::PixelMTFNodeConfig,
        ::NamedTuple{()}, inputs::NamedTuple{(:input,)},
        outputs::NamedTuple{(:output,)}, ::NamedTuple{()}, target) where {T}
    # This calibration is known host storage; decode it before the detector
    # preparation boundary, which resolves array backend traits for its input.
    prepared = Detectors.prepare_pixel_mtf(inputs.input, Matrix(config.transfer);
        padding=config.padding, output=outputs.output)
    return _PixelMTFOwner(prepared)
end

@inline function step_graph_node!(owner::_PixelMTFOwner)
    Detectors.apply_pixel_mtf!(owner.prepared)
    return nothing
end

@inline _pixel_mtf_capture_capability(backend) = GraphNodeCaptureUnsupported()
@inline _pixel_mtf_capture_capability(::Backends.CUDABackend) = GraphNodeCaptureSafe()
@inline _pixel_mtf_capture_capability(::Backends.AMDGPUBackend) = GraphNodeCaptureSafe()

@inline graph_node_capture_capability(owner::_PixelMTFOwner) =
    _pixel_mtf_capture_capability(owner.prepared.plan.backend)

function reset_graph_node!(owner::_PixelMTFOwner)
    fill!(owner.prepared.products.output, zero(eltype(owner.prepared.products.output)))
    fill!(owner.prepared.workspace.buffer, zero(eltype(owner.prepared.workspace.buffer)))
    return nothing
end
