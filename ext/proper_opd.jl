# Runtime sampled OPD steps use the adopted Proper.prop_add_phase convention.
function _runtime_opd_step(step, ::Type{T}) where {T}
    fields = (:op, :name, :input, :size, :axes, :sampling_m, :center_pixel,
        :schema, :outside_opd_m)
    AOG._require_named_fields(step, fields, fields, "runtime OPD step")
    name = AOG._file_string(step.name, "runtime OPD name")
    input = AOG._file_identifier(step.input, "runtime OPD input")
    input in (:pupil_opd, :pupil_amplitude, :output, AOG._PROPER_CONFIG_FIELDS...) &&
        throw(AOG.AlgorithmGraphError("runtime OPD input '$input' is reserved"))
    grid = _sampled_grid(step, T, "runtime OPD")
    schema = AOG._file_string(step.schema, "runtime OPD schema")
    isempty(schema) && throw(AOG.AlgorithmGraphError("runtime OPD schema must not be empty"))
    outside_opd_m = T(AOG._file_real(step.outside_opd_m, "outside_opd_m"))
    isfinite(outside_opd_m) || throw(AOG.AlgorithmGraphError("outside_opd_m must be finite and representable"))
    return merge((; op=:add_opd, name, input), grid, (; schema, outside_opd_m))
end

function _sampled_grid(grid, ::Type{T}, context) where {T}
    dimensions = AOG._file_tuple(grid.size, "$context size")
    length(dimensions) == 2 || throw(AOG.AlgorithmGraphError("$context size requires two dimensions"))
    shape = map(dimensions) do dimension
        value = AOG._file_integer(dimension, "$context dimension")
        value > 0 || throw(AOG.AlgorithmGraphError("$context dimensions must be positive"))
        value
    end
    shape[1] <= typemax(Int) ÷ shape[2] ||
        throw(AOG.AlgorithmGraphError("$context shape is too large"))
    permutation = AOG._file_tuple(grid.axes, "$context axes")
    length(permutation) == 2 || throw(AOG.AlgorithmGraphError("$context axes requires two axes"))
    axes = map(x -> AOG._file_integer(x, "$context axis"), permutation)
    axes in ((1, 2), (2, 1)) ||
        throw(AOG.AlgorithmGraphError("$context axes must be [1,2] or [2,1]"))
    sampling_m = _mask_pair(grid.sampling_m, T, "$context sampling_m")
    all(>(zero(T)), sampling_m) || throw(AOG.AlgorithmGraphError("$context sampling_m must be positive"))
    center_pixel = _mask_pair(grid.center_pixel, T, "$context center_pixel")
    return (; size=shape, axes, sampling_m, center_pixel)
end

function _entrance_pupil_grid(config::NamedTuple, ::Type{T}, p) where {T}
    fields = (:size, :axes, :sampling_m, :center_pixel)
    AOG._require_named_fields(config, fields, fields, "pupil_grid")
    grid = _sampled_grid(config, T, "pupil_grid")
    grid.size == (p,p) || throw(AOG.AlgorithmGraphError("pupil_grid size must match the square entrance pupil ports"))
    return grid
end

_entrance_pupil_grid(config, ::Type, p) =
    throw(AOG.AlgorithmGraphError("pupil_grid must be a mapping"))

function _validate_runtime_opd_inputs(steps)
    inputs = Symbol[]
    for step in steps
        step.op === :add_opd || continue
        step.input in inputs && throw(AOG.AlgorithmGraphError("duplicate runtime OPD input '$(step.input)'"))
        push!(inputs, step.input)
    end
    return nothing
end

function _runtime_opd_ports(::Type{T}, steps) where {T}
    return Tuple(AOG.graph_port_contract(step.input, :input, :data, T,
        step.size, step.schema, :column_major) for step in steps if step.op === :add_opd)
end

# Separable bilinear interpolation of real scalar samples. Index zero denotes an
# outside neighbor. Geometry contains no input/output or scratch bindings.
struct BilinearInterpolationPlan{I,W,A,T}
    row_indices::I
    column_indices::I
    row_weights::W
    column_weights::W
    axes::A
    outside::T
end

function _grid_axis_map(n, source_size, sampling::T, source_sampling, center) where {T}
    indices = Vector{NTuple{2,Int}}(undef, n)
    weights = Vector{T}(undef, n)
    for index in 1:n
        coordinate = T(index-fld(n,2)-1)*sampling/source_sampling+center
        isfinite(coordinate) || throw(AOG.AlgorithmGraphError("runtime OPD coordinate is not representable"))
        if 0 <= coordinate < source_size+1
            lower = floor(Int, coordinate)
            indices[index] = (1 <= lower <= source_size ? lower : 0,
                1 <= lower+1 <= source_size ? lower+1 : 0)
            weights[index] = coordinate-lower
        else
            indices[index] = (0, 0)
            weights[index] = zero(T)
        end
    end
    return FixedSizeVectorDefault{NTuple{2,Int}}(indices), FixedSizeVectorDefault{T}(weights)
end

function _bilinear_interpolation_plan(grid, n, sampling::T, outside::T) where {T}
    isfinite(sampling) && sampling > 0 || throw(AOG.AlgorithmGraphError("runtime OPD plane sampling must be finite and positive"))
    rows, dy = _grid_axis_map(n, grid.size[grid.axes[1]], sampling,
        grid.sampling_m[1], grid.center_pixel[1])
    columns, dx = _grid_axis_map(n, grid.size[grid.axes[2]], sampling,
        grid.sampling_m[2], grid.center_pixel[2])
    return BilinearInterpolationPlan(rows, columns, dy, dx, grid.axes, outside)
end

function _prepare_opd_interpolation(plan::OpticalTrainPlan{T}, masks) where {T}
    context = Proper.RunContext(Matrix{Complex{T}})
    wf = Proper.prop_begin(plan.diameter_m, plan.wavelength_um*T(1e-6), plan.resolution;
        beam_diam_fraction=T(plan.pupil_resolution/plan.resolution), context)
    return map(plan.steps, masks) do step, mask
        interpolation = step.op === :add_opd ?
            _bilinear_interpolation_plan(step, plan.resolution,
                T(Proper.prop_get_sampling(wf)), step.outside_opd_m) : nothing
        _train_optic!(wf, step, _interpolation_host_mask(mask), context)
        interpolation
    end
end

_interpolation_host_mask(::Nothing) = nothing
_interpolation_host_mask(mask) = Backends.host_array(mask)

_prepare_grid_interpolation(::Nothing, target) = nothing
_prepare_grid_interpolation(plan::BilinearInterpolationPlan, ::Backends.HostComputeDevice) = plan

function _interpolation_device_array(source, target)
    destination = Backends.allocate_device_array(target, eltype(source), length(source))
    copyto!(destination, collect(source))
    return destination
end

function _prepare_grid_interpolation(plan::BilinearInterpolationPlan,
        target::Backends.AcceleratorComputeDevice{Backends.CUDABackend})
    return BilinearInterpolationPlan(_interpolation_device_array(plan.row_indices, target),
        _interpolation_device_array(plan.column_indices, target),
        _interpolation_device_array(plan.row_weights, target),
        _interpolation_device_array(plan.column_weights, target), plan.axes, plan.outside)
end

struct RuntimeOPDTrainOwner{O,P,I,G}
    train::O
    interpolation::P
    inputs::I
    entrance::G
end

_require_sampled_train_target(::Backends.HostComputeDevice) = nothing
_require_sampled_train_target(::Backends.AcceleratorComputeDevice{Backends.CUDABackend}) = nothing
_require_sampled_train_target(target) =
    throw(AOG.AlgorithmGraphError("runtime OPD steps and pupil_grid currently support host and CUDA compute devices only"))

function _prepare_runtime_opd_train(plan::OpticalTrainPlan{T}, inputs, outputs, target) where {T}
    _require_sampled_train_target(target)
    length(plan.steps) <= AOG._PROPER_MAX_SAMPLED_STEPS ||
        throw(AOG.AlgorithmGraphError("sampled scalar train exceeds the total optical step limit"))
    roles = Tuple(step.input for step in plan.steps if step.op === :add_opd)
    keys(inputs) == (:pupil_opd, :pupil_amplitude, roles...) ||
        throw(AOG.AlgorithmGraphError("runtime OPD input roles do not match the train"))
    for step in plan.steps
        step.op === :add_opd || continue
        array = getproperty(inputs, step.input)
        eltype(array) === T && size(array) == step.size &&
            !Base.has_offset_axes(array) && Backends.compute_device(array) == target ||
            throw(AOG.AlgorithmGraphError("runtime OPD buffer contract mismatch for '$(step.input)'"))
    end
    train = _prepare_optical_train(plan, inputs, outputs, target).owner
    interpolation = map(_prepare_opd_interpolation(plan, train.workspace.masks)) do grid
        _prepare_grid_interpolation(grid, target)
    end
    opd_inputs = map(plan.steps) do step
        step.op === :add_opd ? getproperty(inputs, step.input) : nothing
    end
    entrance_geometry = plan.pupil_grid === nothing ? nothing :
        _bilinear_interpolation_plan(plan.pupil_grid, plan.resolution,
            plan.diameter_m/T(plan.pupil_resolution), zero(T))
    entrance = _prepare_grid_interpolation(entrance_geometry, target)
    owner = RuntimeOPDTrainOwner(train, interpolation, opd_inputs, entrance)
    _prepare_train_capture!(owner, target)
    return owner
end

@inline function _grid_sample(input, row, column, permutation, outside)
    (iszero(row) || iszero(column)) && return outside
    return permutation[1] == 1 ? input[row, column] : input[column, row]
end

@inline function _interpolated_sample(input, row, column, row_indices, column_indices,
        row_weights, column_weights, permutation, outside)
    i, j = row_indices[row], column_indices[column]
    dy, dx = row_weights[row], column_weights[column]
    return (1-dy)*((1-dx)*_grid_sample(input, i[1], j[1], permutation, outside)+
        dx*_grid_sample(input, i[1], j[2], permutation, outside))+
        dy*((1-dx)*_grid_sample(input, i[2], j[1], permutation, outside)+
        dx*_grid_sample(input, i[2], j[2], permutation, outside))
end

function _interpolate_grid!(output, input, plan::BilinearInterpolationPlan)
    _interpolate_grid!(Backends.execution_style(output), output, input, plan)
    return nothing
end

function _interpolate_grid!(::Backends.ScalarCPUStyle, output, input, plan::BilinearInterpolationPlan)
    for column in axes(output, 2), row in axes(output, 1)
        output[row,column] = _interpolated_sample(input, row, column,
            plan.row_indices, plan.column_indices, plan.row_weights,
            plan.column_weights, plan.axes, plan.outside)
    end
    return nothing
end

Backends.KernelAbstractions.@kernel function _interpolate_grid_kernel!(output, input,
        row_indices, column_indices, row_weights, column_weights, permutation, outside)
    row, column = Backends.KernelAbstractions.@index(Global, NTuple)
    if row <= size(output, 1) && column <= size(output, 2)
        output[row,column] = _interpolated_sample(input, row, column, row_indices,
            column_indices, row_weights, column_weights, permutation, outside)
    end
end

function _interpolate_grid!(style::Backends.AcceleratorStyle, output, input,
        plan::BilinearInterpolationPlan)
    Backends.launch_kernel_async!(style, _interpolate_grid_kernel!, output, input,
        plan.row_indices, plan.column_indices, plan.row_weights, plan.column_weights,
        plan.axes, plan.outside; ndrange=size(output))
    return nothing
end

# Literal tuple accesses avoid boxed tails on long heterogeneous trains. Cold
# admission bounds the generated body to at most 256 ordered optical steps.
@generated function _step_opd_train!(workspace, inputs::I, steps::S, masks::M,
        interpolation::P) where {I<:Tuple,S<:Tuple,M<:Tuple,P<:Tuple}
    calls = [:( _step_opd_optic!(workspace, inputs[$i], steps[$i], masks[$i], interpolation[$i]) )
        for i in 1:fieldcount(S)]
    return quote
        $(calls...)
        nothing
    end
end

function _step_opd_optic!(workspace, input, step, mask, ::Nothing)
    _train_optic!(workspace.wavefront, step, mask, workspace.context, workspace.shifted_mask)
    return nothing
end

function _step_opd_optic!(workspace, input, step, mask, interpolation::BilinearInterpolationPlan)
    _interpolate_grid!(workspace.padded_opd, input, interpolation)
    Proper.prop_add_phase(workspace.wavefront, workspace.padded_opd)
    return nothing
end

function AOG.step_graph_node!(owner::RuntimeOPDTrainOwner)
    train = owner.train
    _begin_sampled_pupil!(train, owner.entrance)
    workspace = train.workspace
    _step_opd_train!(workspace, owner.inputs, train.plan.steps, workspace.masks, owner.interpolation)
    Proper.prop_end(workspace.wavefront, train.output)
    return nothing
end

_begin_sampled_pupil!(train, ::Nothing) = _begin_optical_train!(train)

function _begin_sampled_pupil!(train, entrance::BilinearInterpolationPlan)
    plan, workspace = train.plan, train.workspace
    _interpolate_grid!(workspace.padded_opd, train.pupil_opd, entrance)
    _interpolate_grid!(workspace.padded_amplitude, train.pupil_amplitude, entrance)
    wf = workspace.wavefront
    Proper.prop_begin!(wf, plan.diameter_m, plan.wavelength_um*oftype(plan.wavelength_um,1e-6);
        beam_diam_fraction=plan.pupil_resolution/plan.resolution)
    Proper.prop_multiply(wf, workspace.padded_amplitude)
    Proper.prop_add_phase(wf, workspace.padded_opd)
    return nothing
end

AOG.optical_sampling_m(owner::RuntimeOPDTrainOwner) = AOG.optical_sampling_m(owner.train)
AOG.graph_node_capture_capability(owner::RuntimeOPDTrainOwner) =
    _optical_train_capture_capability(Backends.compute_device(owner.train.output))
AOG.reset_graph_node!(owner::RuntimeOPDTrainOwner) = AOG.reset_graph_node!(owner.train)
