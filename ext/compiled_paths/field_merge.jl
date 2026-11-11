# The fitted SCC's two masked fields are independent through the second mask.
# The admitted suffix is common and linear, so it can propagate one combined
# field. Build a prefix-only checkpoint schedule so no selected traversal
# consumes a checkpoint from an omitted closed or pinhole-only suffix.
struct SCCFieldMerge{S,P,R,I,D}
    steps::S
    prefixes::P
    resampling::R
    image_resample::I
    defocus::D
end

function _prepare_scc_field_merge(owner; share_defocus=false)
    plan, w = owner.plan, owner.workspace
    w.branch_batch === nothing && return nothing
    w.image_resample === nothing && return nothing
    plan.image_size < plan.grid || return nothing
    offset = 0
    for path in plan.paths
        closed, pinhole = path.closed, path.only
        count = length(closed)
        count == length(pinhole) && count >= 7 && (count - 1) % 3 == 0 || return nothing
        offset += length(path.prefix)
        for index in 1:count
            op, other = closed[index], pinhole[index]
            if index == count
                op.kind == PathResample || return nothing
            elseif (index - 1) % 3 == 0
                op.kind == other.kind == PathMultiply || return nothing
            else
                op.kind in (PathFourierTransform, PathTransferFunction) || return nothing
                w.propagation[offset + index] isa Proper.PreparedPropagation || return nothing
            end
            if index == 1 || index == 4
                _compiled_paths_same_operation(op, other) && return nothing
            else
                _compiled_paths_same_operation(op, other) || return nothing
            end
        end
        last_transform = w.propagation[offset + count - 1]
        closed[count - 1].kind == PathTransferFunction && last_transform.direction == 2 || return nothing
        offset += 2count
    end

    empty = _compiled_paths_seal(eltype(first(plan.paths).prefix)[])
    prefix_paths = map(plan.paths) do path
        OpticalPathPlan(path.path_index, path.anchor_index, path.entrance, path.prefix, empty, empty)
    end
    sharing, cache_count = _compiled_paths_shared_prefixes(prefix_paths)
    cache_count <= length(w.shared_prefixes) || return nothing
    steps, prefixes = PathExecutionStep[], Tuple{PathSequenceExecution,Bool}[]
    offset = 0
    outputs = eltype(w.closed_sources)[]
    for (path, shared) in zip(plan.paths, sharing)
        schedule = _compiled_paths_prepare_sequence!(steps, path.prefix, w.propagation, offset, shared.prefix)
        push!(prefixes, (schedule, shared.prefix.count == 0))
        offset += length(path.prefix) + length(path.closed) + length(path.only)
        slot = _compiled_paths_slot(path.path_index, path.anchor_index, plan.n_anchor)
        push!(outputs, w.closed_sources[slot])
    end
    # Reuse the existing image FFT and GEMM workspace with one output per path.
    batch = w.image_resample
    resample = PathResampleBatch(batch.batch_size, batch.field, batch.temporary,
        batch.kernel, batch.forward, _compiled_paths_seal(outputs), batch.coefficients)
    crop = _compiled_paths_prepare_branch_resampling(plan, w.field, plan.paths, resample)
    defocus = share_defocus ? _prepare_shared_scc_defocus(owner) : nothing
    return SCCFieldMerge(_compiled_paths_seal(steps), _compiled_paths_seal(prefixes),
        crop, resample, defocus)
end

struct SharedSCCDefocus{A,F,S,G}
    field::A
    forward::F
    source_slots::S
    group_slots::G
end

function _prepare_shared_scc_defocus(owner)
    plan, fast = owner.plan, owner.workspace.fast
    batch = fast.defocus
    batch === nothing && return nothing
    source_slots = Int[]
    group_slots = zeros(Int32, size(batch.field, 3))
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        plane = fast.defocus_slots[_compiled_paths_slot(k, a, plan.n_anchor)]
        iszero(plane) && continue
        source = _compiled_paths_slot(path.source_index, a, plan.n_anchor)
        group = findfirst(==(source), source_slots)
        if group === nothing
            push!(source_slots, source)
            group = length(source_slots)
        end
        group_slots[plane] = group
    end
    length(source_slots) < length(group_slots) || return nothing
    field = similar(batch.field, size(batch.field, 1), size(batch.field, 2), length(source_slots))
    fill!(field, zero(eltype(field)))
    groups = _compiled_paths_resident(plan.target, group_slots, Int32)
    return SharedSCCDefocus(field, Backends.plan_fft_backend!(field, (1, 2)),
        _compiled_paths_seal(source_slots), groups)
end

KernelAbstractions.@kernel function _apply_shared_scc_defocus!(out, input, transfer, groups)
    y, x, plane = @index(Global, NTuple)
    if y <= size(out, 1) && x <= size(out, 2) && plane <= size(out, 3)
        @inbounds out[y, x, plane] = input[y, x, groups[plane]] * transfer[y, x, plane]
    end
end

function _compiled_paths_defocus!(sources, owner, shared::SharedSCCDefocus)
    plan, batch = owner.plan, owner.workspace.fast.defocus
    for (group, source) in enumerate(shared.source_slots)
        _compiled_paths_pack_filter!(shared.field, sources[source], group, plan.phantom_pad)
    end
    Backends.enqueue_fft_plan!(shared.field, shared.forward)
    backend = KernelAbstractions.get_backend(batch.field)
    _apply_shared_scc_defocus!(backend, (16, 16, 1))(
        batch.field, shared.field, batch.transfer, shared.group_slots; ndrange=size(batch.field))
    Backends.enqueue_fft_plan!(batch.field, batch.inverse)
    return nothing
end

KernelAbstractions.@kernel function _merge_scc_fields!(out, fields, after, scale, gain, closed_mask, pinhole_mask,
        pinhole_gain, edge, before, shift)
    y, x = @index(Global, NTuple)
    n = size(out, 1)
    if y <= n && x <= n
        sy, sx = y + shift, x + shift
        sy > n && (sy -= n)
        sx > n && (sx -= n)
        @inbounds begin
            closed = fields[y, x, 1] * scale
            pinhole = fields[y, x, 2] * scale
            if !isempty(after)
                phase = Proper._propagation_phase(after, y, x, n, n)
                closed *= phase
                pinhole *= phase
            end
            closed = Proper._propagation_input(closed * gain, closed_mask, sy, sx)
            pinhole = Proper._propagation_input(pinhole * gain, pinhole_mask, sy, sx)
            value = Proper._propagation_input(closed + pinhole_gain * pinhole, edge, sy, sx)
            out[y, x] = isempty(before) ? value : value * Proper._propagation_phase(before, y, x, n, n)
        end
    end
end

function _pack_scc_field!(owner, path, offset, ::Val{:fringed}, backend, shape)
    w, plan = owner.workspace, owner.plan
    batch = w.branch_batch
    closed, pinhole = path.closed, path.only
    initial, second, third = (w.propagation[offset + j] for j in (2, 3, 5))
    shift = fld(first(shape), 2)
    Proper._propagation_batch_pack_kernel!(backend, (16, 16))(
        batch.field, (w.field, w.field), (closed[1].array, pinhole[1].array),
        _compiled_paths_edge(w.edge_multiplier, closed[2]), initial.before, shift; ndrange=shape)
    _compiled_paths_branch_fft!(batch, initial)
    Proper._propagation_batch_join_kernel!(backend, (16, 16))(
        batch.field, initial.after, initial.scale, _compiled_paths_branch_gain(closed[2]),
        (nothing, nothing), _compiled_paths_edge(w.edge_multiplier, closed[3]),
        second.before, shift; ndrange=shape)
    _compiled_paths_branch_fft!(batch, second)
    _merge_scc_fields!(backend, (16, 16))(
        third.field, batch.field, second.after, second.scale, _compiled_paths_branch_gain(closed[3]),
        closed[4].array, pinhole[4].array, plan.pinhole_gain[path.anchor_index],
        _compiled_paths_edge(w.edge_multiplier, closed[5]), third.before, shift; ndrange=shape)
    return 5
end

function _pack_scc_field!(owner, path, offset, ::Val{:unfringed}, backend, shape)
    w, closed = owner.workspace, path.closed
    initial = w.propagation[offset + 2]
    Proper._propagation_pack_kernel!(backend, (16, 16))(
        initial.field, w.field, (closed[1].array, _compiled_paths_edge(w.edge_multiplier, closed[2])),
        initial.before, fld(first(shape), 2); ndrange=shape)
    return 2
end

function _propagate_scc_field!(owner, path, offset, index, merge, product)
    plan, w = owner.plan, owner.workspace
    closed = path.closed
    backend = KernelAbstractions.get_backend(w.field)
    shape = size(w.field)
    n = first(shape)
    start = _pack_scc_field!(owner, path, offset, product, backend, shape)
    terminal = length(closed) - 1
    previous = w.propagation[offset + start]
    Proper._propagation_fft!(previous)
    for j in start+1:terminal
        closed[j].kind == PathMultiply && continue
        next = w.propagation[offset + j]
        mask = closed[j-1].kind == PathMultiply ? closed[j-1].array : nothing
        edge = _compiled_paths_edge(w.edge_multiplier, closed[j])
        prior = mask === nothing ? closed[j-1] : closed[j-2]
        Proper._propagation_join_kernel!(backend, (16, 16))(
            previous.field, previous.after, previous.scale, _compiled_paths_branch_gain(prior),
            (mask, edge), next.before, fld(n, 2); ndrange=shape)
        # Keep the terminal transfer in Fourier order for the cropped GEMMs.
        if j == terminal
            Backends.enqueue_fft_plan!(next.field, next.forward)
        else
            Proper._propagation_fft!(next)
        end
        previous = next
    end
    Proper._propagation_transfer_unpack_kernel!(backend, (16, 16))(
        w.physical_fft.scratch, previous.field, previous.transfer, previous.after,
        n^2 * previous.scale, one(closed[terminal].gouy), cld(n, 2); ndrange=shape)
    slot = _compiled_paths_slot(path.path_index, path.anchor_index, plan.n_anchor)
    work = w.operations[offset + length(closed)]
    _compiled_paths_stage_branch_image!(merge.resampling, owner, last(closed), work,
        index, index, w.closed_sources[slot], merge.image_resample)
    return nothing
end

function _compiled_paths_merged!(owner, product, merge::SCCFieldMerge)
    plan, w = owner.plan, owner.workspace
    offset = 0
    for (index, path) in enumerate(plan.paths)
        prefix, from_entrance = merge.prefixes[index]
        from_entrance && Backends.copyto_backend_async!(w.field, path.entrance)
        _compiled_paths_execute_schedule!(owner, path.prefix, offset, prefix, merge.steps)
        offset += length(path.prefix)
        _propagate_scc_field!(owner, path, offset, index, merge, product)
        offset += length(path.closed) + length(path.only)
    end
    for a in 1:plan.n_anchor
        w.closed_sources[a] .*= plan.main_phase
    end
    return _compiled_paths_selected_intensity!(owner, product, merge.defocus)
end
