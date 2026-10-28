struct IdealPupilRelayNode{T<:Union{Float32,Float64}} end

"""Run-immutable matched-grid ideal pupil relay contract."""
struct IdealPupilRelayPlan{T<:AbstractFloat}
    input_size::NTuple{2,Int}
    output_size::NTuple{2,Int}
    axes::NTuple{2,Int}
    offsets::NTuple{2,Int}
    input_sampling_m::NTuple{2,T}
    output_sampling_m::NTuple{2,T}
    wavelength_m::T
    input_schema::String
    output_schema::String
end

function _relay_size(value, context)
    pair = _file_tuple(value, context)
    length(pair) == 2 || throw(AlgorithmGraphError("$context requires two dimensions"))
    dimensions = map(x -> _file_integer(x, context), pair)
    all(>(0), dimensions) || throw(AlgorithmGraphError("$context must be positive"))
    dimensions[1] <= typemax(Int) ÷ dimensions[2] ||
        throw(AlgorithmGraphError("$context exceeds the host index range"))
    return dimensions
end

function _relay_sampling(value, ::Type{T}, context) where {T}
    pair = _file_tuple(value, context)
    length(pair) == 2 || throw(AlgorithmGraphError("$context requires two sample pitches"))
    sampling = map(x -> T(_file_real(x, context)), pair)
    all(x -> isfinite(x) && x > zero(T), sampling) ||
        throw(AlgorithmGraphError("$context must be finite, positive and representable"))
    return sampling
end

"""
    ideal_pupil_relay_node(name; input_size, output_size, axes=(2,1),
        input_sampling_m, output_sampling_m, wavelength_m,
        input_schema, output_schema, T=Float32)

Declare an ideal isotropic pupil relay between two explicitly sampled planes.
The `:input` and `:output` ports carry dimensionless `Complex{T}` field samples.
`axes[k]` selects the source array axis for output axis `k`. Both grids place
zero on sample `fld(size[k],2)+1`, including even dimensions. A centered window
is copied one sample per sample, without interpolation, phase conversion,
source scaling or renormalization of light outside the window. Sample pitches
are in array-axis order, and their ratios must declare one isotropic physical
magnification after the axis permutation. Output dimensions must fit inside the
permuted input dimensions. This contract declares an ideal optical relay; it
does not infer measured registration from an array shape.
"""
function ideal_pupil_relay_node(name::Symbol;
        input_size, output_size, axes=(2,1), input_sampling_m,
        output_sampling_m, wavelength_m::Real,
        input_schema::AbstractString, output_schema::AbstractString,
        T::Type{<:Union{Float32,Float64}}=Float32)
    source_size = _relay_size(input_size, "ideal pupil relay input_size")
    destination_size = _relay_size(output_size, "ideal pupil relay output_size")
    permutation = _file_tuple(axes, "ideal pupil relay axes")
    length(permutation) == 2 || throw(AlgorithmGraphError("ideal pupil relay axes requires two axes"))
    permutation = map(x -> _file_integer(x, "ideal pupil relay axis"), permutation)
    permutation in ((1,2), (2,1)) ||
        throw(AlgorithmGraphError("ideal pupil relay axes must be (1,2) or (2,1)"))
    all(k -> destination_size[k] <= source_size[permutation[k]], 1:2) ||
        throw(AlgorithmGraphError("ideal pupil relay output_size must fit inside the permuted input_size"))
    source_sampling = _relay_sampling(input_sampling_m, T, "ideal pupil relay input_sampling_m")
    destination_sampling = _relay_sampling(output_sampling_m, T, "ideal pupil relay output_sampling_m")
    magnification = ntuple(k -> destination_sampling[k] / source_sampling[permutation[k]], 2)
    all(x -> isfinite(x) && x > zero(T), magnification) &&
        isapprox(magnification[1], magnification[2]; rtol=8eps(T), atol=zero(T)) ||
        throw(AlgorithmGraphError("ideal pupil relay requires finite isotropic magnification"))
    wavelength = T(_file_real(wavelength_m, "ideal pupil relay wavelength_m"))
    isfinite(wavelength) && wavelength > zero(T) ||
        throw(AlgorithmGraphError("ideal pupil relay wavelength_m must be positive and representable"))
    isempty(input_schema) && throw(AlgorithmGraphError("ideal pupil relay input_schema must not be empty"))
    isempty(output_schema) && throw(AlgorithmGraphError("ideal pupil relay output_schema must not be empty"))
    offsets = ntuple(k -> fld(source_size[permutation[k]], 2) - fld(destination_size[k], 2), 2)
    plan = IdealPupilRelayPlan(source_size, destination_size, permutation, offsets,
        source_sampling, destination_sampling, wavelength, String(input_schema), String(output_schema))
    return algorithm_node(name, IdealPupilRelayNode{T}, plan; props=NamedTuple())
end

function _ideal_pupil_relay_node(::Type{T}, name::Symbol,
        config::NamedTuple, props::NamedTuple) where {T}
    fields = (:input_size, :output_size, :axes, :input_sampling_m,
        :output_sampling_m, :wavelength_m, :input_schema, :output_schema)
    required = (:input_size, :output_size, :input_sampling_m,
        :output_sampling_m, :wavelength_m, :input_schema, :output_schema)
    _require_named_fields(config, fields, required, "ideal pupil relay config")
    _require_named_fields(props, (), (), "ideal pupil relay props")
    return ideal_pupil_relay_node(name; config..., T)
end

_ideal_pupil_relay_f32_node(name, config, props) =
    _ideal_pupil_relay_node(Float32, name, config, props)
_ideal_pupil_relay_f64_node(name, config, props) =
    _ideal_pupil_relay_node(Float64, name, config, props)

function graph_node_ports(::Type{IdealPupilRelayNode{T}}, plan::IdealPupilRelayPlan{T}) where {T}
    return (graph_port_contract(:input, :input, :data, Complex{T}, plan.input_size,
            plan.input_schema, :column_major),
        graph_port_contract(:output, :output, :data, Complex{T}, plan.output_size,
            plan.output_schema, :column_major))
end

struct _IdealPupilRelayOwner{P,I,O}
    plan::P
    input::I
    output::O
end

function prepare_graph_node(::Type{IdealPupilRelayNode{T}}, plan::IdealPupilRelayPlan{T},
        ::NamedTuple{()}, inputs::NamedTuple{(:input,)}, outputs::NamedTuple{(:output,)},
        ::NamedTuple{()}, target) where {T}
    for (array, shape) in ((inputs.input, plan.input_size), (outputs.output, plan.output_size))
        size(array) == shape && eltype(array) === Complex{T} &&
            !Base.has_offset_axes(array) && compute_device(array) == target ||
            throw(AlgorithmGraphError("ideal pupil relay buffer contract mismatch"))
    end
    Base.mightalias(inputs.input, outputs.output) &&
        throw(AlgorithmGraphError("ideal pupil relay requires distinct input and output storage"))
    return _IdealPupilRelayOwner(plan, inputs.input, outputs.output)
end

@inline function _relay_sample(input, row, column, offsets, permutation)
    i, j = row + offsets[1], column + offsets[2]
    return permutation[1] == 1 ? input[i,j] : input[j,i]
end

function _copy_pupil_relay!(::ScalarCPUStyle, output, input, plan)
    for column in axes(output,2), row in axes(output,1)
        output[row,column] = _relay_sample(input, row, column, plan.offsets, plan.axes)
    end
    return nothing
end

@kernel function _pupil_relay_kernel!(output, input, offsets, permutation)
    row, column = @index(Global, NTuple)
    if row <= size(output,1) && column <= size(output,2)
        output[row,column] = _relay_sample(input, row, column, offsets, permutation)
    end
end

function _copy_pupil_relay!(style::AcceleratorStyle, output, input, plan)
    launch_kernel_async!(style, _pupil_relay_kernel!, output, input,
        plan.offsets, plan.axes; ndrange=size(output))
    return nothing
end

@inline function step_graph_node!(owner::_IdealPupilRelayOwner)
    _copy_pupil_relay!(execution_style(owner.output), owner.output, owner.input, owner.plan)
    return nothing
end

@inline _pupil_relay_capture_capability(::ScalarCPUStyle) = GraphNodeCaptureUnsupported()
@inline _pupil_relay_capture_capability(::AcceleratorStyle) = GraphNodeCaptureSafe()
@inline graph_node_capture_capability(owner::_IdealPupilRelayOwner) =
    _pupil_relay_capture_capability(execution_style(owner.output))

function reset_graph_node!(owner::_IdealPupilRelayOwner)
    fill!(owner.output, zero(eltype(owner.output)))
    return nothing
end
