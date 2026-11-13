# Coordinated scalar Proper propagation. This is an execution optimization of
# existing independent paths, not an optical authoring or continuation format.
struct SharedOpticalTrainNode{T} end

struct SharedOpticalTrainPlan{P}
    paths::P
    prefix_count::Int
end

struct SharedOpticalTrainOwner{P,B,G,I,V}
    prefix::P
    branches::B
    geometry::G
    inputs::I
    validation::V
end

_scalar_train(owner::OpticalTrainOwner) = owner
_scalar_train(owner::RuntimeOPDTrainOwner) = owner.train

_shared_port(index, role) = Symbol(:path_, index, :_, role)
_shared_output(index) = Symbol(:path_, index)
_shared_pixels(index) = Symbol(:pixels_, index)

_shared_pixel_configuration(owner, input) = nothing
function _shared_pixel_configuration(owner::AOG._PixelIntegrationOwner, input)
    return owner.input === input ? owner.plan.config : nothing
end

function _shared_path_configuration(graph, train)
    nodes = values(graph.nodes)
    all(node -> isempty(node.parameters), nodes) && length(graph.outputs) == 1 || return nothing
    if length(nodes) == 1
        only(values(graph.outputs)) === train.output || return nothing
        return (; pixels=nothing)
    end
    pixels = _shared_pixel_configuration(last(nodes).owner, train.output)
    pixels === nothing && return nothing
    only(values(graph.outputs)) === last(nodes).outputs.output || return nothing
    return (; pixels)
end

function _same_train_entrance(a, b)
    return typeof(a.diameter_m) === typeof(b.diameter_m) &&
        a.resolution == b.resolution && a.pupil_resolution == b.pupil_resolution &&
        a.diameter_m == b.diameter_m && a.wavelength_um == b.wavelength_um &&
        a.schemas == b.schemas && isequal(a.pupil_grid, b.pupil_grid) &&
        a.output_quantity === b.output_quantity
end

function AOG._prepare_external_pupil_propagation(
        owners::Tuple{Vararg{Union{OpticalTrainOwner,RuntimeOPDTrainOwner}}},
        graphs, observations)
    length(owners) >= 2 || return nothing
    trains = map(_scalar_train, owners)
    paths = map(train -> train.plan, trains)
    reference = first(paths)
    all(plan -> _same_train_entrance(reference, plan), paths) || return nothing
    policies = map(AOG.graph_execution_policy, graphs)
    policy = first(policies)
    policy isa Union{AOG.StreamGraphExecution,AOG.CapturedGraphExecution} || return nothing
    all(other -> typeof(other) === typeof(policy), policies) || return nothing
    target = AOG.compute_device(first(graphs))
    all(graph -> AOG.compute_device(graph) == target, graphs) || return nothing
    # Admit only a scalar root and optional ordinary terminal pixel integration.
    # Exact input identity proves the terminal reads this root's raw product.
    configurations = map(_shared_path_configuration, graphs, trains)
    any(isnothing, configurations) && return nothing
    # Independent graphs can encode an implicit ordering dependency through
    # caller-owned root inputs. Preserve their original chronology rather than
    # batching an input that a different path's raw or exported product changes.
    products = Tuple(output for (graph, train) in zip(graphs, trains)
        for output in (train.output, only(values(graph.outputs))))
    any(input -> any(output -> Base.mightalias(input, output), products),
        (input for graph in graphs for input in values(first(values(graph.nodes)).inputs))) && return nothing
    staged = Tuple(getproperty(first(values(graph.nodes)).inputs, role)
        for graph in graphs for role in (:pupil_opd,:pupil_amplitude))
    # Distinct overlapping native views can change an earlier path's pupil
    # while a later path is staged. Exact same-quantity sharing is safe; the
    # external-pupil contract already rejects cross-quantity aliases.
    for index in eachindex(staged), other in 1:index-1
        staged[index] === staged[other] && isodd(index) == isodd(other) && continue
        Base.mightalias(staged[index],staged[other]) && return nothing
    end
    runtime_inputs = (getproperty(first(values(graph.nodes)).inputs, role)
        for graph in graphs for role in keys(first(values(graph.nodes)).inputs)
        if role ∉ (:pupil_opd,:pupil_amplitude))
    any(input -> any(pupil -> Base.mightalias(input,pupil),staged),runtime_inputs) && return nothing
    prefix_count = 0
    for index in 1:minimum(plan -> length(plan.steps), paths)
        step = reference.steps[index]
        all(plan -> isequal(step, plan.steps[index]), paths) || break
        prefix_count = index
    end
    # Keep ordinary complete trajectories when there is no useful fork.
    0 < prefix_count < minimum(plan -> length(plan.steps), paths) || return nothing
    any(step -> step.op in (:propagate, :propagate_to_focus),
        reference.steps[1:prefix_count]) || return nothing
    inputs = Tuple(AOG.graph_input(_shared_port(index, role),
            :propagation => _shared_port(index, role), getproperty(node.inputs, role))
        for (index, graph) in enumerate(graphs)
        for node in (first(values(graph.nodes)),)
        for role in keys(node.inputs))
    outputs = Tuple(AOG.graph_output(_shared_output(index),
        configurations[index].pixels === nothing ? :propagation => _shared_output(index) :
            _shared_pixels(index) => :output) for index in eachindex(paths))
    plan = SharedOpticalTrainPlan(paths, prefix_count)
    T = typeof(reference.diameter_m)
    node = AOG.algorithm_node(:propagation, SharedOpticalTrainNode{T}, plan)
    pixels = Tuple(AOG.algorithm_node(_shared_pixels(index), AOG.PixelIntegrationNode{T},
        configurations[index].pixels) for index in eachindex(paths)
        if configurations[index].pixels !== nothing)
    links = Tuple(AOG.link(:propagation => _shared_output(index), _shared_pixels(index) => :input)
        for index in eachindex(paths) if configurations[index].pixels !== nothing)
    graph = AOG.prepare_algorithm_graph(AOG.algorithm_graph((node,pixels...);
        name=:shared_optical_prefix, inputs, outputs, links); target, execution=policy)
    bindings = Tuple(binding for index in eachindex(trains)
        for binding in (configurations[index].pixels === nothing ?
            ((AOG.graph_output(graph,_shared_output(index)),trains[index].output),) :
            ((getproperty(graph.nodes.propagation.outputs,_shared_output(index)),trains[index].output),
             (AOG.graph_output(graph,_shared_output(index)),only(values(graphs[index].outputs))))))
    publication = AOG._prepare_external_pupil_publication(graph, bindings)
    return (; graph, publication, state=AOG._ExternalPupilPropagationState(true))
end

function AOG.graph_node_ports(::Type{SharedOpticalTrainNode{T}}, plan::SharedOpticalTrainPlan) where {T}
    return Tuple(begin
            direction = AOG._port_direction(port)
            name = direction === :output ? _shared_output(index) :
                _shared_port(index, AOG._port_name(port))
            numeric_type = direction === :output ? _output_type(T, path.output_quantity) : T
            AOG.graph_port_contract(name, direction, AOG._port_role(port),
                numeric_type, port.format.shape, port.format.schema,
                port.format.layout)
        end for (index, path) in enumerate(plan.paths)
        for port in AOG.graph_node_ports(OpticalTrainNode{T}, path))
end

function _shared_train_inputs(inputs, index, plan, ::Type{T}) where {T}
    roles = Tuple(AOG._port_name(port) for port in AOG.graph_node_ports(OpticalTrainNode{T}, plan)
        if AOG._port_direction(port) === :input)
    return NamedTuple{roles}(map(role -> getproperty(inputs, _shared_port(index, role)), roles))
end

function _train_with_steps(plan::OpticalTrainPlan, steps)
    return OpticalTrainPlan(steps, plan.resolution, plan.pupil_resolution,
        plan.diameter_m, plan.wavelength_um, plan.schemas, plan.pupil_grid,
        plan.output_quantity)
end

# Exact adopted Proper.WaveFront scalar geometry. Field/workspace identities
# belong to separate owners and must never be restored from a fork snapshot.
function _wavefront_geometry(wf)
    return (; wf.wavelength_m, wf.sampling_m, wf.z_m, wf.beam_diameter_m,
        wf.z_w0_m, wf.w0_m, wf.z_rayleigh_m, wf.current_fratio,
        wf.reference_surface, wf.beam_type_old, wf.propagator_type, wf.rayleigh_factor)
end

function _restore_wavefront_geometry!(wf, geometry)
    for name in keys(geometry)
        setproperty!(wf, name, getproperty(geometry, name))
    end
    return nothing
end

function _enqueue_train_prefix!(owner::OpticalTrainOwner)
    _begin_optical_train!(owner)
    workspace = owner.workspace
    foreach(owner.plan.steps, workspace.masks) do step, mask
        _train_optic!(workspace.wavefront, step, mask, workspace.context, workspace.shifted_mask)
    end
    return nothing
end

function _enqueue_train_prefix!(owner::RuntimeOPDTrainOwner)
    train = owner.train
    _begin_sampled_pupil!(train, owner.entrance)
    _step_opd_train!(train.workspace, owner.inputs, train.plan.steps,
        train.workspace.masks, owner.interpolation)
    return nothing
end

function _shared_train_suffix(owner::OpticalTrainOwner, prefix_count)
    return (; owner, steps=owner.plan.steps[prefix_count+1:end],
        masks=owner.workspace.masks[prefix_count+1:end])
end

function _shared_train_suffix(owner::RuntimeOPDTrainOwner, prefix_count)
    return (; owner, steps=owner.train.plan.steps[prefix_count+1:end],
        masks=owner.train.workspace.masks[prefix_count+1:end],
        inputs=owner.inputs[prefix_count+1:end],
        interpolation=owner.interpolation[prefix_count+1:end])
end

function _enqueue_train_suffix!(branch, owner::OpticalTrainOwner)
    workspace = owner.workspace
    foreach(branch.steps, branch.masks) do step, mask
        _train_optic!(workspace.wavefront, step, mask, workspace.context, workspace.shifted_mask)
    end
    return nothing
end

function _enqueue_train_suffix!(branch, owner::RuntimeOPDTrainOwner)
    _step_opd_train!(owner.train.workspace, branch.inputs, branch.steps,
        branch.masks, branch.interpolation)
    return nothing
end

function AOG.prepare_graph_node(::Type{SharedOpticalTrainNode{T}}, plan::SharedOpticalTrainPlan,
        ::NamedTuple{()}, inputs::NamedTuple, outputs::NamedTuple,
        ::NamedTuple{()}, target) where {T}
    path_inputs = map(eachindex(plan.paths)) do index
        _shared_train_inputs(inputs, index, plan.paths[index], T)
    end
    reference = first(plan.paths)
    prefix_plan = _train_with_steps(reference, reference.steps[1:plan.prefix_count])
    prefix_output = Backends.allocate_device_array(target, _output_type(T, reference.output_quantity),
        reference.resolution, reference.resolution)
    prefix_inputs = Tuple(_shared_train_inputs(inputs, index, prefix_plan, T)
        for index in eachindex(plan.paths))
    prefix = AOG.prepare_graph_node(OpticalTrainNode{T}, prefix_plan, (;),
        first(prefix_inputs), (;output=prefix_output), (;), target)
    # Propagation geometry depends only on the frozen train. Establish the
    # split on the host before warming or recording any suffix trajectory.
    _enqueue_train_prefix!(prefix)
    geometry = _wavefront_geometry(_scalar_train(prefix).workspace.wavefront)
    branches = map(eachindex(plan.paths)) do index
        owner = AOG.prepare_graph_node(OpticalTrainNode{T}, plan.paths[index], (;),
            path_inputs[index], (;output=getproperty(outputs, _shared_output(index))), (;), target)
        _shared_train_suffix(owner, plan.prefix_count)
    end
    validation = _prepare_shared_train_validation(target, prefix_inputs)
    return SharedOpticalTrainOwner(prefix, Tuple(branches), geometry, prefix_inputs, validation)
end

struct SharedTrainValidation{C,D,H,A,P,F}
    context::C
    device::D
    host::H
    arrays::A
    runtime::P
    full::F
end

@kernel function _shared_train_equal!(flags, a, b, start)
    index = @index(Global, Linear)
    lane = @index(Local, Linear)
    block = @index(Group, Linear)
    mismatch = @localmem UInt32 (256,)
    mismatch[lane] = UInt32(index <= length(a) && a[index] != b[index])
    @synchronize
    for offset in (128,64,32,16,8,4,2,1)
        if lane <= offset
            mismatch[lane] |= mismatch[lane+offset]
        end
        @synchronize
    end
    if lane == 1
        flags[start+block] = UInt8(mismatch[1])
    end
end

_prepare_shared_train_validation(target, inputs) = Backends._prepare_device_execution_context(target)

function _prepare_shared_train_validation(
        target::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
        inputs)
    reference = first(inputs)
    roles = (:pupil_opd, :pupil_amplitude)
    pairs = Tuple((getproperty(reference, role), getproperty(branch, role))
        for role in keys(reference) for branch in Base.tail(inputs) if role ∉ roles)
    pairs = (pairs..., Tuple((getproperty(reference, role), getproperty(branch, role))
        for role in roles for branch in Base.tail(inputs))...)
    runtime_pairs = length(pairs) - length(roles) * length(Base.tail(inputs))
    counts = map(pair -> cld(length(first(pair)), 256), pairs)
    context = Backends._prepare_device_execution_context(target)
    # Leave unusual shapes on the portable comparison path.
    all(pair -> axes(first(pair)) == axes(last(pair)) && !Base.has_offset_axes(first(pair)), pairs) ||
        return context
    all(>(0), counts) || return context
    return Backends._with_prepared_device_execution_context(context) do
        device = Backends.allocate_device_array(target, UInt8, sum(counts))
        host = FixedSizeVector{UInt8,Vector{UInt8}}(undef, length(device))
        style = Backends.execution_style(device)
        function enqueue(selected)
            start = 0
            for ((a,b), blocks) in zip(selected, counts)
                Backends.launch_kernel_async!(style, _shared_train_equal!, device, a, b, start;
                    ndrange=256blocks)
                start += blocks
            end
            return nothing
        end
        runtime = pairs[1:runtime_pairs]
        enqueue(pairs)
        Backends._synchronize_prepared_device_execution_context!(context)
        full = Backends._capture_prepared_device_kernel_graph(() -> enqueue(pairs), context)
        prefix = isempty(runtime) ? nothing :
            Backends._capture_prepared_device_kernel_graph(() -> enqueue(runtime), context)
        arrays = (device, Tuple(array for pair in pairs for array in pair)...)
        return SharedTrainValidation(context, device, host, arrays,
            (prefix, sum(counts[1:runtime_pairs]; init=0)), full)
    end
end

function AOG._shared_prefix_inputs_match(owner::SharedOpticalTrainOwner, include_pupil::Bool=true)
    return _shared_prefix_inputs_match(owner.validation, owner.inputs, include_pupil)
end

function _shared_prefix_inputs_match(context, inputs, include_pupil)
    return Backends._with_prepared_device_execution_context(context) do
        _shared_train_inputs_match(inputs, include_pupil)
    end
end

function _shared_prefix_inputs_match(validation::SharedTrainValidation, inputs, include_pupil)
    capture, count = include_pupil ? (validation.full, length(validation.host)) : validation.runtime
    capture === nothing && return true
    context = validation.context
    Backends._with_prepared_device_execution_context(context) do
        Backends._launch_prepared_device_graph!(capture, context, validation.arrays)
        Backends._copy_prepared_device_buffer!(context, parent(validation.host), validation.device)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
    end
    for index in 1:count
        iszero(validation.host[index]) || return false
    end
    return true
end

function _shared_train_inputs_match(inputs, include_pupil)
    reference = first(inputs)
    for branch in Base.tail(inputs), role in keys(reference)
        !include_pupil && role in (:pupil_opd,:pupil_amplitude) && continue
        a, b = getproperty(reference, role), getproperty(branch, role)
        # Bulk reduction stays outside capture and never permits unequal
        # branch-bound DM/phase/amplitude values to use the first path's field.
        mapreduce(==, &, a, b) || return false
    end
    return true
end

function _validate_shared_train_inputs!(owner)
    AOG._shared_prefix_inputs_match(owner) || throw(AOG.AlgorithmGraphError(
        "shared optical prefix requires equal runtime inputs in every branch"))
    return nothing
end

function _enqueue_shared_train!(owner)
    _enqueue_train_prefix!(owner.prefix)
    source = _scalar_train(owner.prefix).workspace.wavefront.field
    for branch in owner.branches
        train = _scalar_train(branch.owner)
        _restore_wavefront_geometry!(train.workspace.wavefront, owner.geometry)
        Backends.copyto_backend_async!(train.workspace.wavefront.field, source)
        _enqueue_train_suffix!(branch, branch.owner)
        _finish_optical_train!(train.workspace.wavefront, train.output, train.plan.output_quantity)
    end
    return nothing
end

function AOG.step_graph_node!(owner::SharedOpticalTrainOwner)
    _validate_shared_train_inputs!(owner)
    return _enqueue_shared_train!(owner)
end

AOG.enqueue_captured_graph_node!(owner::SharedOpticalTrainOwner) = _enqueue_shared_train!(owner)
AOG._preflight_captured_graph_node!(owner::SharedOpticalTrainOwner) = _validate_shared_train_inputs!(owner)
AOG.graph_node_capture_capability(owner::SharedOpticalTrainOwner) =
    _optical_train_capture_capability(Backends.compute_device(_scalar_train(owner.prefix).output))

function AOG.reset_graph_node!(owner::SharedOpticalTrainOwner)
    AOG.reset_graph_node!(owner.prefix)
    for branch in owner.branches
        AOG.reset_graph_node!(branch.owner)
    end
    return nothing
end
