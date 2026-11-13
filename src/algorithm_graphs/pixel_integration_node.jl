struct PixelIntegrationNode{T<:Union{Float32,Float64}} end

"""Cold geometry and quantity contract for rectangular pixel integration."""
struct PixelIntegrationNodeConfig
    input_shape::NTuple{2,Int}
    output_shape::NTuple{2,Int}
    input_sampling_m::NTuple{2,Float64}
    output_sampling_m::NTuple{2,Float64}
    input_center_pixel::NTuple{2,Float64}
    output_center_pixel::NTuple{2,Float64}
    input_schema::String
    output_schema::String
end

function _pixel_integration_shape(value, context)
    pair = _file_tuple(value, context)
    length(pair) == 2 || throw(AlgorithmGraphError("$context requires two dimensions"))
    all(n -> n isa Integer && !(n isa Bool) && 0 < n <= min(typemax(Int)-1, 2^52-1), pair) ||
        throw(AlgorithmGraphError("$context requires positive dimensions in the exact host coordinate range"))
    shape = map(Int, pair)
    shape[1] <= typemax(Int) ÷ shape[2] ||
        throw(AlgorithmGraphError("$context exceeds the host index range"))
    return shape
end

function _pixel_integration_pair(value, context; positive=false)
    pair = _file_tuple(value, context)
    length(pair) == 2 || throw(AlgorithmGraphError("$context requires two coordinates"))
    values = map(x -> Float64(_file_real(x, context)), pair)
    all(x -> isfinite(x) && (!positive || x > 0), values) ||
        throw(AlgorithmGraphError("$context must be finite and representable"))
    return values
end

@inline function _pixel_integration_edges(index, ratio, source_center, destination_center)
    offset = Float64(index) - destination_center
    return (offset - 0.5) * ratio + source_center,
        (offset + 0.5) * ratio + source_center
end

"""
    pixel_integration_node(name; input_shape, output_shape, input_sampling_m,
        output_sampling_m, input_center_pixel, output_center_pixel,
        input_schema, output_schema, T=Float32)

Integrate cell-integrated flux onto axis-aligned rectangular output pixels.
All geometry is explicit in logical `(y,x)` array-axis order. Pitches are in
metres; centers are fractional one-based pixel coordinates at physical zero.
Input flux is uniform within each source cell. Each output receives the sum of
input flux times overlap area divided by input-cell area. Outside the source
grid the flux is zero; cropped flux is lost without renormalization. Quantity
units are preserved, including cell-integrated photon rate or relative power;
schemas must state the caller's quantity and coordinates. This operation adds
no exposure, optical magnification, detector response or noise.

Preparation builds compact separable overlap maps using Float64 geometry and
binds distinct one-based buffers on the selected CPU, CUDA or AMDGPU device.
Execution accepts finite nonnegative Float32 or Float64 flux and uses prepared
storage. Captured accelerator execution checks the full input before publishing
a completed graph step; invalid values require a graph reset before reuse.
`optical_sampling_m` reports the output pitch only for square output pixels.
"""
function pixel_integration_node(name::Symbol; input_shape, output_shape,
        input_sampling_m, output_sampling_m, input_center_pixel, output_center_pixel,
        input_schema::AbstractString, output_schema::AbstractString,
        T::Type{<:Union{Float32,Float64}}=Float32)
    source_shape = _pixel_integration_shape(input_shape, "pixel integration input_shape")
    destination_shape = _pixel_integration_shape(output_shape, "pixel integration output_shape")
    source_sampling = _pixel_integration_pair(input_sampling_m,
        "pixel integration input_sampling_m"; positive=true)
    destination_sampling = _pixel_integration_pair(output_sampling_m,
        "pixel integration output_sampling_m"; positive=true)
    source_center = _pixel_integration_pair(input_center_pixel, "pixel integration input_center_pixel")
    destination_center = _pixel_integration_pair(output_center_pixel, "pixel integration output_center_pixel")
    isempty(input_schema) && throw(AlgorithmGraphError("pixel integration input_schema must not be empty"))
    isempty(output_schema) && throw(AlgorithmGraphError("pixel integration output_schema must not be empty"))
    for axis in 1:2
        ratio = destination_sampling[axis] / source_sampling[axis]
        isfinite(ratio) && ratio > 0 || throw(AlgorithmGraphError(
            "pixel integration pitch ratio must be finite and positive"))
        for index in (1, destination_shape[axis])
            lo, hi = _pixel_integration_edges(index, ratio, source_center[axis], destination_center[axis])
            isfinite(lo) && isfinite(hi) && lo < hi || throw(AlgorithmGraphError(
                "pixel integration cell boundaries must be finite and distinguishable"))
        end
    end
    config = PixelIntegrationNodeConfig(source_shape, destination_shape,
        source_sampling, destination_sampling, source_center, destination_center,
        String(input_schema), String(output_schema))
    return algorithm_node(name, PixelIntegrationNode{T}, config; props=NamedTuple())
end

function _pixel_integration_node(::Type{T}, name::Symbol, config::NamedTuple,
        props::NamedTuple) where {T}
    fields = (:input_shape, :output_shape, :input_sampling_m, :output_sampling_m,
        :input_center_pixel, :output_center_pixel, :input_schema, :output_schema)
    _require_named_fields(config, fields, fields, "pixel integration config")
    _require_named_fields(props, (), (), "pixel integration props")
    return pixel_integration_node(name; config..., T)
end

_pixel_integration_f32_node(name, config, props) = _pixel_integration_node(Float32, name, config, props)
_pixel_integration_f64_node(name, config, props) = _pixel_integration_node(Float64, name, config, props)

function graph_node_ports(::Type{PixelIntegrationNode{T}}, config::PixelIntegrationNodeConfig) where {T}
    return (graph_port_contract(:input, :input, :data, T, config.input_shape,
            config.input_schema, :column_major),
        graph_port_contract(:output, :output, :data, T, config.output_shape,
            config.output_schema, :column_major))
end

struct _PixelIntegrationAxis{Offsets,Indices,Weights}
    offsets::Offsets
    indices::Indices
    weights::Weights
end

function _pixel_integration_axis(::Type{T}, config, axis) where {T}
    source_count = config.input_shape[axis]
    output_count = config.output_shape[axis]
    ratio = config.output_sampling_m[axis] / config.input_sampling_m[axis]
    offsets = Vector{Int}(undef, output_count + 1)
    indices, weights = Int[], T[]
    for index in 1:output_count
        offsets[index] = length(indices) + 1
        lo, hi = _pixel_integration_edges(index, ratio,
            config.input_center_pixel[axis], config.output_center_pixel[axis])
        isfinite(lo) && isfinite(hi) && lo < hi || throw(AlgorithmGraphError(
            "pixel integration cell boundaries must be finite and distinguishable"))
        # Clip before integer conversion so arbitrarily displaced finite grids
        # cannot overflow an index. Boundaries here are in source-pixel units.
        left, right = max(lo, 0.5), min(hi, source_count + 0.5)
        if left < right
            first_cell = max(1, floor(Int, left + 0.5))
            last_cell = min(source_count, ceil(Int, right - 0.5))
            for cell in first_cell:last_cell
                overlap = min(right, cell + 0.5) - max(left, cell - 0.5)
                if overlap > 0
                    push!(indices, cell)
                    push!(weights, T(overlap))
                end
            end
        end
    end
    offsets[end] = length(indices) + 1
    return _PixelIntegrationAxis(FixedSizeVectorDefault{Int}(offsets),
        FixedSizeVectorDefault{Int}(indices), FixedSizeVectorDefault{T}(weights))
end

struct _PixelIntegrationPlan{R,C}
    config::PixelIntegrationNodeConfig
    rows::R
    columns::C
end

struct _PixelIntegrationValidation{D,H,C}
    device::D
    host::H
    context::C
end

struct _PixelIntegrationOwner{P,I,O,V}
    plan::P
    input::I
    output::O
    validation::V
end

_pixel_integration_axis(target::HostComputeDevice, axis::_PixelIntegrationAxis) = axis
function _pixel_integration_axis(target::AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}},
        axis::_PixelIntegrationAxis)
    arrays = map((axis.offsets, axis.indices, axis.weights)) do values
        result = allocate_device_array(target, eltype(values), length(values))
        copyto!(result, collect(values))
        result
    end
    return _PixelIntegrationAxis(arrays...)
end
_pixel_integration_axis(target, axis::_PixelIntegrationAxis) =
    throw(AlgorithmGraphError("pixel integration requires CPU, CUDA or AMDGPU execution"))

_pixel_integration_validation(::HostComputeDevice) = nothing
function _pixel_integration_validation(target::AcceleratorComputeDevice)
    device = allocate_device_array(target, UInt32, 1)
    fill!(device, UInt32(0))
    # A fixed owner with Vector backing permits the backend's ordinary bulk
    # device-to-host copy without a new allocation or a scalar device read.
    host = FixedSizeVector{UInt32,Vector{UInt32}}(undef, 1)
    fill!(host, UInt32(0))
    return _PixelIntegrationValidation(device, host, Backends._current_device_execution_context(device))
end

function prepare_graph_node(::Type{PixelIntegrationNode{T}}, config::PixelIntegrationNodeConfig,
        ::NamedTuple{()}, inputs::NamedTuple{(:input,)}, outputs::NamedTuple{(:output,)},
        ::NamedTuple{()}, target) where {T}
    for (array, shape) in ((inputs.input, config.input_shape), (outputs.output, config.output_shape))
        size(array) == shape && eltype(array) === T && !Base.has_offset_axes(array) &&
            compute_device(array) == target ||
            throw(AlgorithmGraphError("pixel integration buffer contract mismatch"))
    end
    Base.mightalias(inputs.input, outputs.output) &&
        throw(AlgorithmGraphError("pixel integration requires distinct input and output storage"))
    plan = _PixelIntegrationPlan(config,
        _pixel_integration_axis(target, _pixel_integration_axis(T, config, 1)),
        _pixel_integration_axis(target, _pixel_integration_axis(T, config, 2)))
    return _PixelIntegrationOwner(plan, inputs.input, outputs.output,
        _pixel_integration_validation(target))
end

function _integrate_pixels!(::ScalarCPUStyle, owner::_PixelIntegrationOwner)
    input, output = owner.input, owner.output
    all(x -> isfinite(x) && x >= 0, input) || throw(AlgorithmGraphError(
        "pixel integration requires finite nonnegative cell-integrated flux"))
    rows, columns = owner.plan.rows, owner.plan.columns
    for column in axes(output, 2), row in axes(output, 1)
        flux = zero(eltype(output))
        for x in columns.offsets[column]:(columns.offsets[column+1]-1),
                y in rows.offsets[row]:(rows.offsets[row+1]-1)
            flux += input[rows.indices[y],columns.indices[x]] * rows.weights[y] * columns.weights[x]
        end
        output[row,column] = flux
    end
    return nothing
end

@kernel function _validate_pixel_integration_kernel!(failure, input)
    index = @index(Global, Linear)
    if index <= length(input)
        value = input[index]
        if !isfinite(value) || value < 0
            KernelAbstractions.@atomic failure[1] = UInt32(1)
        end
    end
end

@kernel function _pixel_integration_kernel!(output, input, row_offsets, row_indices,
        row_weights, column_offsets, column_indices, column_weights)
    row, column = @index(Global, NTuple)
    if row <= size(output,1) && column <= size(output,2)
        flux = zero(eltype(output))
        for x in column_offsets[column]:(column_offsets[column+1]-1),
                y in row_offsets[row]:(row_offsets[row+1]-1)
            flux += input[row_indices[y],column_indices[x]] * row_weights[y] * column_weights[x]
        end
        output[row,column] = flux
    end
end

function _enqueue_pixel_integration!(style::AcceleratorStyle, owner)
    fill!(owner.validation.device, UInt32(0))
    launch_kernel_async!(style, _validate_pixel_integration_kernel!,
        owner.validation.device, owner.input; ndrange=length(owner.input))
    rows, columns = owner.plan.rows, owner.plan.columns
    launch_kernel_async!(style, _pixel_integration_kernel!, owner.output, owner.input,
        rows.offsets, rows.indices, rows.weights,
        columns.offsets, columns.indices, columns.weights; ndrange=size(owner.output))
    return nothing
end

function _check_pixel_integration!(validation::_PixelIntegrationValidation)
    Backends._read_prepared_device_buffer!(validation.context, parent(validation.host), validation.device)
    iszero(only(validation.host)) || throw(AlgorithmGraphError(
        "pixel integration requires finite nonnegative cell-integrated flux"))
    return nothing
end
_check_pixel_integration!(::Nothing) = nothing

function _integrate_pixels!(style::AcceleratorStyle, owner::_PixelIntegrationOwner)
    _enqueue_pixel_integration!(style, owner)
    _check_pixel_integration!(owner.validation)
    return nothing
end

step_graph_node!(owner::_PixelIntegrationOwner) =
    _integrate_pixels!(execution_style(owner.output), owner)
enqueue_captured_graph_node!(owner::_PixelIntegrationOwner) =
    _enqueue_pixel_integration!(execution_style(owner.output), owner)
_complete_captured_graph_node!(owner::_PixelIntegrationOwner) =
    _check_pixel_integration!(owner.validation)

_pixel_integration_capture_capability(::AbstractComputeDevice) = GraphNodeCaptureUnsupported()
_pixel_integration_capture_capability(
    ::AcceleratorComputeDevice{<:Union{Backends.CUDABackend,Backends.AMDGPUBackend}}) = GraphNodeCaptureSafe()
graph_node_capture_capability(owner::_PixelIntegrationOwner) =
    _pixel_integration_capture_capability(compute_device(owner.output))

function reset_graph_node!(owner::_PixelIntegrationOwner)
    fill!(owner.output, zero(eltype(owner.output)))
    _reset_pixel_integration_validation!(owner.validation)
    return nothing
end
_reset_pixel_integration_validation!(::Nothing) = nothing
function _reset_pixel_integration_validation!(validation::_PixelIntegrationValidation)
    fill!(validation.device, UInt32(0))
    fill!(validation.host, UInt32(0))
    return nothing
end

function optical_sampling_m(config::PixelIntegrationNodeConfig)
    y, x = config.output_sampling_m
    y == x || throw(AlgorithmGraphError("scalar optical_sampling_m requires square output pixels"))
    return y
end

optical_sampling_m(owner::_PixelIntegrationOwner) = optical_sampling_m(owner.plan.config)
