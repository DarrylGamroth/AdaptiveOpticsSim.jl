# Derived from SPIDERSSim.jl b479871e017e8f1df2e56880bba6f47f80d1273a
# src/optics/physical15.jl. Preserves the lowered spectre Fresnel and spectral
# arithmetic; instrument-specific selection and export remain cold data.
using FixedSizeArrays: FixedSizeVectorDefault, FixedSizeMatrixDefault
import KernelAbstractions
using KernelAbstractions: @index, @localmem, @synchronize

# Representation and field primitives precede cold preparation. The executor
# below keeps ordinary numerical routines beside their prepared lowerings.
include("compiled_paths/types.jl")
include("compiled_paths/fields.jl")
include("compiled_paths/scheduling.jl")
include("compiled_paths/preparation.jl")
include("compiled_paths/execution.jl")

# Execution: prepared and ordinary propagation, checkpoint continuations,
# independent field batches, then spectral synthesis and pixel response.
@inline _compiled_paths_edge(::Nothing, op) = isempty(op.auxiliary) ? nothing : op.auxiliary
@inline _compiled_paths_edge(cache::Tuple, op) =
    op.auxiliary === first(cache) ? last(cache) : _compiled_paths_edge(nothing, op)
@inline _compiled_paths_mul_coefficients(w, prepared::Bool) =
    prepared ? w.prepared_mul_coefficients : w.mul_coefficients

_compiled_paths_proper!(workspace, op, ::Nothing) = _compiled_paths_proper!(workspace, op)

function _compiled_paths_proper!(workspace, op, prepared::Proper.PreparedPropagation)
    edge = _compiled_paths_edge(workspace.edge_multiplier, op)
    Proper.propagate!(workspace.field, prepared, workspace.field; centered=true,
        input_multiplier=edge,
        output_multiplier=op.kind == PathFourierTransform ? op.gouy : one(op.gouy))
    return nothing
end

function _compiled_paths_proper!(workspace, op)
    wf, ctx = workspace.wavefront, workspace.context
    if !isempty(op.auxiliary)
        workspace.field .*= op.auxiliary
    end
    _compiled_paths_shift!(wf.field, workspace.field, true)
    wf.wavelength_m = op.wavelength_m
    wf.sampling_m = op.dx
    wf.z_m = zero(op.dx)
    if op.kind == PathTransferFunction
        wf.reference_surface = Proper.PLANAR
        Proper.prop_ptp(wf, op.dz, ctx)
    elseif op.from_waist
        wf.reference_surface = Proper.PLANAR
        Proper.prop_wts(wf, op.dz, ctx)
    else
        wf.reference_surface = Proper.SPHERICAL
        Proper.prop_stw(wf, op.dz, ctx)
    end
    _compiled_paths_shift!(workspace.field, wf.field)
    if op.kind == PathFourierTransform
        workspace.field .*= op.gouy
    end
    return nothing
end

function _compiled_paths_proper_pair!(workspace, op, next, first::Proper.PreparedPropagation, second::Proper.PreparedPropagation)
    edge(value) = _compiled_paths_edge(workspace.edge_multiplier, value)
    gain(value) = value.kind == PathFourierTransform ? value.gouy : one(value.gouy)
    Proper.propagate!(workspace.field, (first,second), workspace.field; centered=true,
        input_multipliers=(edge(op),edge(next)), output_multipliers=(gain(op),gain(next)))
    return nothing
end
function _compiled_paths_proper_pair!(workspace, op, next, ::Nothing, ::Nothing)
    _compiled_paths_proper!(workspace, op)
    _compiled_paths_proper!(workspace, next)
    return nothing
end

function _compiled_paths_proper_group!(workspace, operations, masks,
    prepared::NTuple{N,P},input=workspace.field,initial_mask=nothing) where {N,P<:Proper.PreparedPropagation}
    edge(op)=_compiled_paths_edge(workspace.edge_multiplier, op)
    gain(op)=op.kind==PathFourierTransform ? op.gouy : one(op.gouy)
    inputs=ntuple(Val(N)) do i
        i==1 ? (initial_mask,edge(operations[i])) :
            isodd(i) ? (masks[fld(i-1,2)].array,edge(operations[i])) : edge(operations[i])
    end
    Proper.propagate!(workspace.field,prepared,input;centered=true,
        input_multipliers=inputs,output_multipliers=map(gain,operations))
    return nothing
end
function _compiled_paths_proper_group!(workspace,operations,masks,::NTuple{N,Nothing},
    input=workspace.field,initial_mask=nothing) where {N}
    input === workspace.field || Backends.copyto_backend_async!(workspace.field,input)
    initial_mask === nothing || (workspace.field .*=initial_mask)
    for i in 1:N
        if i>1 && isodd(i)
            workspace.field .*=masks[fld(i-1,2)].array
        end
        _compiled_paths_proper!(workspace,operations[i])
    end
    return nothing
end

function _compiled_paths_apply_resampling_roll!(frequency, op, work, prepared)
    if !isempty(op.auxiliary) && !(prepared && work.identity_roll)
        frequency .*= op.auxiliary
    end
    return nothing
end

function _compiled_paths_resample_frequency!(w,op,work,prepared::Bool=true,finalize::Bool=true)
    _compiled_paths_apply_resampling_roll!(w.physical_fft.scratch, op, work, prepared)
    Backends._enqueue_mul!(work.temporary,op.array,w.physical_fft.scratch,_compiled_paths_mul_coefficients(w,prepared))
    Backends._enqueue_mul!(w.field,work.temporary,transpose(op.array),_compiled_paths_mul_coefficients(w,prepared))
    if finalize
        normalization = size(op.array,2)^2
        if prepared && !isempty(op.post)
            @. w.field = (w.field / normalization) * op.post
        else
            w.field ./=normalization
            isempty(op.post) || (w.field .*=op.post)
        end
    end
    return nothing
end
function _compiled_paths_resample!(w,op,work,prepared::Bool=true)
    _compiled_paths_fft!(w.physical_fft.scratch,w.field,w.physical_fft)
    return _compiled_paths_resample_frequency!(w,op,work,prepared)
end

function _compiled_paths_proper_resample!(w,op,tf,resample,work,
    first::Proper.PreparedPropagation,last::Proper.PreparedPropagation,::Val{true})
    # PTP ends with BFFT followed by scale; resampling immediately applies FFT.
    # FFT(BFFT(G)*s)=length(G)*s*G. Preserve phases, edge/gain positions and the
    # centered Fourier coordinates; only the redundant FFT pair is removed.
    edge(value)=_compiled_paths_edge(w.edge_multiplier, value)
    gain=op.kind==PathFourierTransform ? op.gouy : one(op.gouy)
    backend=KernelAbstractions.get_backend(first.field)
    shift=fld(size(first.field,1),2)
    Proper._propagation_pack_kernel!(backend,(16,16))(
        first.field,w.field,edge(op),first.before,shift;ndrange=size(w.field))
    Proper._propagation_fft!(first)
    Proper._propagation_join_kernel!(backend,(16,16))(
        first.field,first.after,first.scale,gain,edge(tf),last.before,shift;ndrange=size(w.field))
    last.forward*last.field
    Proper._propagation_transfer_kernel!(backend,(16,16))(
        last.field,last.transfer;ndrange=size(w.field))
    Proper._propagation_unpack_kernel!(backend,(16,16))(
        w.physical_fft.scratch,last.field,last.after,length(last.field)*last.scale,
        one(tf.gouy),cld(size(last.field,1),2);ndrange=size(w.field))
    return _compiled_paths_resample_frequency!(w,resample,work)
end
function _compiled_paths_proper_resample!(w,op,tf,resample,work,first,last,::Val{false})
    _compiled_paths_proper_pair!(w,op,tf,first,last)
    return _compiled_paths_resample!(w,resample,work)
end
function _compiled_paths_execute_group!(owner,sequence,offset,index,::Val{N},
    input=owner.workspace.field,initial_mask=nothing) where {N}
    return _compiled_paths_execute_group!(owner,sequence,offset,index,Val(N),input,initial_mask,
        owner.workspace.propagation[offset+index])
end
function _compiled_paths_execute_group!(owner,sequence,offset,index,::Val{N},
    input,initial_mask,::P) where {N,P}
    positions=ntuple(i->index+3fld(i-1,2)+mod(i-1,2),Val(N))
    masks=ntuple(i->sequence[index+3i-1],Val(N÷2-1))
    # Preparation seals one propagation type (or all Nothing). Narrow it before
    # constructing the tuple, rather than boxing Union-valued tuple elements.
    prepared=ntuple(i->owner.workspace.propagation[offset+positions[i]]::P,Val(N))
    _compiled_paths_proper_group!(owner.workspace,map(i->sequence[i],positions),masks,
        prepared,input,initial_mask)
    return nothing
end

function _compiled_paths_interpolate!(work, kernel, source, coefficients)
    work.input .= source
    Backends._enqueue_mul!(work.temporary, kernel, work.input, coefficients)
    Backends._enqueue_mul!(work.mapped, work.temporary, transpose(kernel), coefficients)
    work.opd .= real.(work.mapped)
    return work.opd
end

# Independent ordinary traversal: no schedule, grouping or prepared propagation.
# Keep the original operation arithmetic as the numerical oracle.
function _compiled_paths_execute_reference!(owner, sequence, offset, first, last)
    w = owner.workspace
    for index in first:last
        op = sequence[index]
        work = w.operations[offset + index]
        if op.kind == PathMultiply
            w.field .*= op.array
        elseif op.kind == PathNormalize
            _compiled_paths_normalize!(w.field, w.normalization)
        elseif op.kind in (PathTransferFunction, PathFourierTransform)
            _compiled_paths_proper!(w, op)
        elseif op.kind == PathSpectralMultiplier
            _compiled_paths_filter!(w.field, w.field, op.array, w.physical_fft)
        elseif op.kind == PathResample
            _compiled_paths_resample!(w,op,work,false)
        elseif op.kind == PathLiveScreen
            input = op.source == 1 ? owner.inputs.dm_opd : owner.inputs.entrance_opd
            _compiled_paths_interpolate!(work, op.array, input, _compiled_paths_mul_coefficients(w,false))
            AdaptiveOpticsSim.Optics.apply_phase_async!(work.electric_field, work.opd)
            if op.source == 2
                _compiled_paths_interpolate!(work, op.array, owner.inputs.entrance_attenuation,
                    _compiled_paths_mul_coefficients(w,false))
                @. w.field *= exp(-work.opd)
            end
        end
    end
    return nothing
end

@inline _compiled_paths_slot(path, anchor, n_anchor) = (path - 1) * n_anchor + anchor

function _compiled_paths_sequence!(owner,sequence,offset,sharing,share_prefix,prepared,stage)
    if prepared
        schedule = share_prefix ? owner.execution.shared[stage] : owner.execution.unshared[stage]
        _compiled_paths_execute_schedule!(owner,sequence,offset,schedule)
    elseif !share_prefix
        _compiled_paths_execute_reference!(owner,sequence,offset,1,length(sequence))
    else
        first=sharing.count+1
        sharing.count>0 && Backends.copyto_backend_async!(owner.workspace.field,
            owner.workspace.shared_prefixes[sharing.cache_index])
        for checkpoint in sharing.checkpoints
            _compiled_paths_execute_reference!(owner,sequence,offset,first,checkpoint.index)
            Backends.copyto_backend_async!(owner.workspace.shared_prefixes[checkpoint.cache_index],owner.workspace.field)
            first=checkpoint.index+1
        end
        _compiled_paths_execute_reference!(owner,sequence,offset,first,length(sequence))
    end
    return offset + length(sequence)
end

# Inputs remain independent through the resampling FFT and both matrix
# products. Batch outputs are consumed only after all slow paths complete.
function _compiled_paths_prepare_image_resample(plan, closed, only, coefficients, ::Val{N}) where {N}
    N > 0 || throw(ArgumentError("compiled image resampling batch must be positive"))
    n = plan.image_size
    field = _compiled_paths_array(plan.target, eltype(plan.resample_kernel), n, n, N)
    temporary = similar(field)
    # fftshift(F) = P*F*transpose(P). Copy K*P by permuting K's columns,
    # so both products accept raw FFT order. Coefficients are not reevaluated;
    # the changed dot-product summation order can change roundoff.
    columns = mod1.((1:n) .+ fld(n, 2), n)
    values = reshape(Array(plan.resample_kernel)[:, columns], n, n, 1)
    kernel = _compiled_paths_resident(plan.target, values, eltype(field))
    outputs = eltype(closed)[]
    for path in plan.paths
        slot = _compiled_paths_slot(path.path_index, path.anchor_index, plan.n_anchor)
        push!(outputs, closed[slot], only[slot])
    end
    return PathResampleBatch(Val(N), field, temporary, kernel,
        Backends.plan_fft_backend!(field, (1, 2)), _compiled_paths_seal(outputs), coefficients)
end

# Preparation checks crop bounds. The wrapped coordinates stay in 1:n.
# A final partial batch clears unused planes to prevent stale values growing
# on replay; transforms and products never combine separate batch planes.
@inline _compiled_paths_resample_pixel(value, y, x, ::Nothing) = value

@inline function _compiled_paths_resample_pixel(value, y, x, finalization::Tuple)
    normalization, post = finalization
    value = value / normalization
    return isempty(post) ? value : value * post[y, x]
end

KernelAbstractions.@kernel function _compiled_paths_pack_resample_kernel!(out, input, start, slot, valid, finalization)
    y, x = @index(Global, NTuple)
    n = size(out, 1)
    if y <= n && x <= n
        sy, sx = y + fld(n, 2), x + fld(n, 2)
        sy > n && (sy -= n)
        sx > n && (sx -= n)
        @inbounds out[y, x, slot] = _compiled_paths_resample_pixel(
            input[sy + start, sx + start], sy + start, sx + start, finalization)
        for plane in valid + 1:size(out, 3)
            @inbounds out[y, x, plane] = zero(eltype(out))
        end
    end
end

KernelAbstractions.@kernel function _compiled_paths_unpack_resample_kernel!(outputs, field, valid, normalization)
    y, x, plane = @index(Global, NTuple)
    if y <= size(field, 1) && x <= size(field, 2) && plane <= valid
        @inbounds outputs[plane][y, x] = field[y, x, plane] / normalization
    end
end

function _compiled_paths_stage_image!(::Nothing, owner, number, out, prepared, finalization=nothing)
    w, plan = owner.workspace, owner.plan
    _compiled_paths_crop!(w.crop, w.field, plan.crop_start)
    _compiled_paths_resample_image!(out, w.crop, owner, prepared)
    return nothing
end

function _compiled_paths_stage_image!(batch::PathResampleBatch{N}, owner, number, out, prepared,
        finalization=nothing, input=owner.workspace.field, start=owner.plan.crop_start) where {N}
    prepared || return _compiled_paths_stage_image!(nothing, owner, number, out, prepared)
    n = owner.plan.image_size
    slot = mod1(number, N)
    valid = number == length(batch.outputs) ? slot : N
    backend = KernelAbstractions.get_backend(batch.field)
    _compiled_paths_pack_resample_kernel!(backend, (16, 16))(
        batch.field, input, start, slot, valid, finalization; ndrange = (n, n))
    if slot == N || number == length(batch.outputs)
        Backends.enqueue_fft_plan!(batch.field, batch.forward)
        AdaptiveOpticsSim.Optics._compiled_paths_resample_mul!(batch.field, batch.temporary,
            batch.kernel, batch.coefficients)
        first = number - slot
        outputs = ntuple(i -> batch.outputs[min(first + i, length(batch.outputs))], Val(N))
        _compiled_paths_unpack_resample_kernel!(backend, (16, 16, 1))(
            outputs, batch.field, slot, n^2; ndrange = size(batch.field))
    end
    return nothing
end

_compiled_paths_terminal_finalization(::Nothing, op) = nothing

_compiled_paths_terminal_finalization(::PathResampleBatch, op) = (size(op.array, 2)^2, op.post)

# Match independent suffixes without changing the existing prefix trie or its
# cache allocation. Restrict cuts to common multiplier/transform-pair regions.
# Other sequences retain their ordinary executor.
function _compiled_paths_branch_pair(path, sharing, propagation, offset)
    a, b = path.closed, path.only
    length(a) == length(b) && !isempty(a) && sharing.closed.count == sharing.only.count || return nothing
    a[end].kind == b[end].kind == PathResample || return nothing
    _compiled_paths_same_operation(a[end], b[end]) || return nothing
    left, right = sharing.closed.checkpoints, sharing.only.checkpoints
    length(left) == length(right) || return nothing
    all(i -> left[i].index == right[i].index, eachindex(left)) || return nothing
    # Only may depend on a checkpoint produced by this closed suffix. It must
    # then wait for that producer in the original serial order.
    sharing.only.count > 0 && any(p -> p.cache_index == sharing.only.cache_index, left) && return nothing
    any(x -> any(y -> x.cache_index == y.cache_index, right), left) && return nothing
    segments = [(last=p.index, closed_cache=p.cache_index, only_cache=right[i].cache_index)
        for (i, p) in enumerate(left)]
    push!(segments, (last=length(a)-1, closed_cache=0, only_cache=0))
    prepared_segments = PathBranchSegment[]
    first = sharing.closed.count + 1
    for segment in segments
        last = segment.last
        last >= first && (last - first + 1) % 3 == 0 || return nothing
        for i in first:3:last
            a[i].kind == b[i].kind == PathMultiply || return nothing
            for j in i + 1:i + 2
                a[j].kind in (PathTransferFunction, PathFourierTransform) &&
                    _compiled_paths_same_operation(a[j], b[j]) || return nothing
                p, q = propagation[offset + j], propagation[offset + length(a) + j]
                p !== nothing && p === q || return nothing
            end
        end
        transforms = Union{PathBranchTransform{false},PathBranchTransform{true}}[]
        for j in first+2:last
            a[j].kind == PathMultiply && continue
            mask = a[j-1].kind == PathMultiply ? j-1 : 0
            push!(transforms,PathBranchTransform{j==length(a)-1}(j,mask))
        end
        push!(prepared_segments,PathBranchSegment(last,segment.closed_cache,
            segment.only_cache,_compiled_paths_seal(transforms)))
        first = last + 1
    end
    a[end - 1].kind == PathTransferFunction || return nothing
    propagation[offset + length(a) - 1].direction == 2 || return nothing
    return PathBranchPair(sharing.closed.count + 1, _compiled_paths_seal(prepared_segments))
end

function _compiled_paths_branch_pairs(plan, propagation)
    S = typeof(_compiled_paths_seal(PathBranchSegment[]))
    pairs = Union{Nothing,PathBranchPair{S}}[]
    offset = 0
    for (i, path) in enumerate(plan.paths)
        offset += length(path.prefix)
        push!(pairs, _compiled_paths_branch_pair(path, plan.shared_prefixes[i], propagation, offset))
        offset += length(path.closed) + length(path.only)
    end
    return _compiled_paths_seal(pairs)
end
function _compiled_paths_branch_batch(field, pairs, resampling=nothing)
    n = size(field, 1)
    planes = similar(field, n, n, 2)
    fill!(planes, zero(eltype(planes)))
    return PathBranchBatch(planes, (view(planes, :, :, 1), view(planes, :, :, 2)),
        Backends.plan_fft_backend!(planes, (1, 2)),
        Backends.plan_bfft_backend!(planes, (1, 2)), pairs, resampling)
end

function _compiled_paths_prepare_branch_resampling(plan, field::A, pairs, image_resample) where {A}
    image_resample === nothing && return nothing
    plan.image_size < plan.grid || return nothing
    rows = plan.crop_start .+ (1:plan.image_size)
    kernels, posts = IdDict{A,A}(), IdDict{A,A}()
    operators = Union{Nothing,Tuple{A,A}}[]
    for (path, pair) in zip(plan.paths, pairs)
        if pair === nothing
            push!(operators, nothing)
            continue
        end
        # Pair admission proves both terminal operations share these arrays.
        op = last(path.closed)
        kernel = get!(kernels, op.array) do
            op.array[rows, :]
        end
        post = if isempty(op.post)
            op.post
        else
            get!(posts, op.post) do
                op.post[rows, rows]
            end
        end
        push!(operators, (kernel, post))
    end
    return PathBranchResampling(_compiled_paths_seal(operators),
        similar(field, plan.image_size, plan.grid))
end

function _compiled_paths_prepare_branch_batch(plan, field, propagation, image_resample)
    AdaptiveOpticsSim.Optics._compiled_paths_branch_batch_enabled(plan.target, field) || return nothing
    pairs = _compiled_paths_branch_pairs(plan, propagation)
    all(isnothing, pairs) && return nothing
    resampling = _compiled_paths_prepare_branch_resampling(plan, field, pairs, image_resample)
    return _compiled_paths_branch_batch(field, pairs, resampling)
end

function _compiled_paths_stage_branch_image!(::Nothing, owner, op, work, index, number, output)
    w = owner.workspace
    finalization = _compiled_paths_terminal_finalization(w.image_resample, op)
    _compiled_paths_resample_frequency!(w, op, work, true, finalization === nothing)
    _compiled_paths_stage_image!(w.image_resample, owner, number, output, true, finalization)
    return nothing
end

function _compiled_paths_stage_branch_image!(resampling::PathBranchResampling{O,A},
        owner, op, work, index, number, output, image_resample=owner.workspace.image_resample) where {O,A}
    w = owner.workspace
    kernel, post = resampling.operators[index]::Tuple{A,A}
    frequency = w.physical_fft.scratch
    _compiled_paths_apply_resampling_roll!(frequency, op, work, true)
    coefficients = w.prepared_mul_coefficients
    # The camera consumes only the crop, so select the output rows on both
    # sides of the Fourier re-evaluation. Keep the full frequency axes.
    Backends._enqueue_mul!(resampling.temporary, kernel, frequency, coefficients)
    Backends._enqueue_mul!(w.crop, resampling.temporary, transpose(kernel), coefficients)
    _compiled_paths_stage_image!(image_resample, owner, number, output, true,
        (owner.plan.grid^2, post), w.crop, 0)
    return nothing
end
@inline _compiled_paths_branch_pair(::Nothing, index) = nothing

@inline _compiled_paths_branch_pair(batch::PathBranchBatch, index) = batch.pairs[index]

@inline _compiled_paths_branch_gain(op) = op.kind == PathFourierTransform ? op.gouy : one(op.gouy)

function _compiled_paths_branch_fft!(batch, propagation, ::Val{Terminal}=Val(false)) where {Terminal}
    if propagation.direction == 2
        Backends.enqueue_fft_plan!(batch.field, batch.forward)
        n = size(batch.field, 1)
        if !Terminal
            Proper._propagation_batch_transfer_kernel!(KernelAbstractions.get_backend(batch.field), (16, 16))(
                batch.field, propagation.transfer; ndrange = (n, n))
            Backends.enqueue_fft_plan!(batch.field, batch.backward)
        end
    elseif propagation.direction == 1
        Backends.enqueue_fft_plan!(batch.field, batch.forward)
    elseif propagation.direction == -1
        Backends.enqueue_fft_plan!(batch.field, batch.backward)
    end
    return nothing
end

_compiled_paths_branch_fft!(batch, propagation, ::PathBranchTransform{Terminal}) where {Terminal} =
    _compiled_paths_branch_fft!(batch,propagation,Val(Terminal))

function _compiled_paths_branches!(owner, path, sharing, offset, index, pair)
    w = owner.workspace
    batch = w.branch_batch
    n = size(batch.field, 1)
    backend = KernelAbstractions.get_backend(batch.field)
    a, b = path.closed, path.only
    other = offset + length(a)
    # Both initial cache reads precede either suffix's new stores. Stage-level
    # cache lifetimes therefore remain conservative; paired cache outputs are
    # disjoint and no current closed producer feeds the paired only input.
    inputs = (sharing.closed.count > 0 ? w.shared_prefixes[sharing.closed.cache_index] : w.field,
        sharing.only.count > 0 ? w.shared_prefixes[sharing.only.cache_index] : w.field)
    first = pair.first
    for segment in pair.segments
        initial = first + 1
        previous = w.propagation[offset + initial]
        previous_op = a[initial]
        Proper._propagation_batch_pack_kernel!(backend, (16, 16))(
            batch.field, inputs, (a[first].array, b[first].array),
            _compiled_paths_edge(w.edge_multiplier, previous_op), previous.before, fld(n, 2); ndrange = (n, n))
        _compiled_paths_branch_fft!(batch, previous)
        for step in segment.transforms
            j = step.index
            next = w.propagation[offset + j]
            masks = step.mask == 0 ? (nothing,nothing) : (a[step.mask].array,b[step.mask].array)
            Proper._propagation_batch_join_kernel!(backend, (16, 16))(
                batch.field, previous.after, previous.scale, _compiled_paths_branch_gain(previous_op),
                masks, _compiled_paths_edge(w.edge_multiplier, a[j]), next.before, fld(n, 2); ndrange = (n, n))
            _compiled_paths_branch_fft!(batch,next,step)
            previous, previous_op = next, a[j]
        end
        if segment.closed_cache == 0
            # Preserve FFT(BFFT(G)*s)=length(G)*s*G at the terminal PTP.
            # Consume one plane at a time through the existing scalar mapping
            # scratch; neither resampling GEMM writes the paired raw field.
            slot = _compiled_paths_slot(path.path_index, path.anchor_index, owner.plan.n_anchor)
            for plane in 1:2
                Proper._propagation_transfer_unpack_kernel!(backend, (16, 16))(
                    w.physical_fft.scratch, batch.planes[plane], previous.transfer, previous.after,
                    n^2 * previous.scale, one(previous_op.gouy), cld(n, 2); ndrange = (n, n))
                operation = plane == 1 ? a[end] : b[end]
                work = w.operations[(plane == 1 ? offset : other) + length(a)]
                output = plane == 1 ? w.closed_sources[slot] : w.only_sources[slot]
                _compiled_paths_stage_branch_image!(batch.resampling, owner, operation, work,
                    index, 2index - 2 + plane, output)
            end
        else
            outputs = (w.shared_prefixes[segment.closed_cache], w.shared_prefixes[segment.only_cache])
            Proper._propagation_batch_unpack_kernel!(backend, (16, 16))(
                outputs, batch.field, previous.after, previous.scale,
                _compiled_paths_branch_gain(previous_op), cld(n, 2); ndrange = (n, n))
            inputs = outputs
        end
        first = segment.last + 1
    end
    return other+length(b)
end

function _compiled_paths_slow!(owner,share_prefix::Bool=true,prepared::Bool=true)
    plan,w=owner.plan,owner.workspace
    offset=0
    for (index,path) in enumerate(plan.paths)
        sharing=plan.shared_prefixes[index]
        # Each reused prefix is produced earlier in this invocation, including
        # during capture. Scratch offsets remain those of the original program;
        # cached fields do not restore wavefront geometry or FFT bindings.
        if !share_prefix || sharing.prefix.count==0
            Backends.copyto_backend_async!(w.field,path.entrance)
        end
        offset=_compiled_paths_sequence!(owner,path.prefix,offset,sharing.prefix,share_prefix,prepared,3index-2)
        pair = prepared && share_prefix ? _compiled_paths_branch_pair(w.branch_batch,index) : nothing
        if pair !== nothing
            offset = _compiled_paths_branches!(owner,path,sharing,offset,index,pair)
            continue
        end
        if !share_prefix || sharing.only.count==0
            Backends.copyto_backend_async!(w.prefix,w.field)
        end
        offset=_compiled_paths_sequence!(owner,path.closed,offset,sharing.closed,share_prefix,prepared,3index-1)
        slot=_compiled_paths_slot(path.path_index,path.anchor_index,plan.n_anchor)
        _compiled_paths_stage_image!(w.image_resample,owner,2index-1,w.closed_sources[slot],prepared)
        if !share_prefix || sharing.only.count==0
            Backends.copyto_backend_async!(w.field,w.prefix)
        end
        offset=_compiled_paths_sequence!(owner,path.only,offset,sharing.only,share_prefix,prepared,3index)
        _compiled_paths_stage_image!(w.image_resample,owner,2index,w.only_sources[slot],prepared)
    end
    for a in 1:plan.n_anchor
        w.closed_sources[a] .*=plan.main_phase
        w.only_sources[a] .*=plan.main_phase
    end
    return nothing
end

function _compiled_paths_resample_image!(out, input, owner, prepared::Bool=true)
    plan, w = owner.plan, owner.workspace
    _compiled_paths_fft!(w.image_fft.scratch, input, w.image_fft)
    Backends._enqueue_mul!(w.resample_temporary, plan.resample_kernel, w.image_fft.scratch,
        _compiled_paths_mul_coefficients(w,prepared))
    Backends._enqueue_mul!(out, w.resample_temporary, transpose(plan.resample_kernel),
        _compiled_paths_mul_coefficients(w,prepared))
    out ./= plan.image_size^2
    return out
end

function _compiled_paths_fast_scalar!(envelopes, sources, owner, free_ghost)
    plan, w = owner.plan, owner.workspace
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        out = envelopes[_compiled_paths_slot(k, a, plan.n_anchor)]
        if path.source_index == 0
            if free_ghost
                Backends.copyto_backend_async!(out, path.ghost[a])
            else
                fill!(out, zero(eltype(out)))
            end
        else
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            if isempty(path.defocus_transfer)
                Backends.copyto_backend_async!(out, source)
            else
                _compiled_paths_pad!(w.phantom_fft.scratch, source, plan.phantom_pad)
                _compiled_paths_filter!(w.phantom_fft.scratch, w.phantom_fft.scratch,
                    path.defocus_transfer[a], w.phantom_fft)
                _compiled_paths_crop!(out, w.phantom_fft.scratch, plan.phantom_pad)
            end
        end
        _compiled_paths_filter!(out, out, path.shift_transfer, w.image_fft)
        out .*= path.phase_multiplier[a]
    end
    return nothing
end

function _compiled_paths_fast!(envelopes, sources, owner, free_ghost)
    plan, w = owner.plan, owner.workspace
    fast = w.fast
    # All bindings and batch slots are frozen at preparation. Closed and only
    # branches reuse the same scratch in stream order, including capture replay.
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        defocus_slot = fast.defocus_slots[slot]
        if defocus_slot != 0
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            _compiled_paths_pack_filter!(fast.defocus.field, source, defocus_slot, plan.phantom_pad)
        end
    end
    _compiled_paths_filter_batch!(fast.defocus)
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        defocus_slot = fast.defocus_slots[slot]
        if defocus_slot != 0
            _compiled_paths_pack_defocus!(fast.image.field, fast.defocus.field,
                slot, defocus_slot, plan.phantom_pad)
        elseif path.source_index == 0
            _compiled_paths_pack_filter!(fast.image.field, path.ghost[a], slot, 0, free_ghost)
        else
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            _compiled_paths_pack_filter!(fast.image.field, source, slot, 0)
        end
    end
    _compiled_paths_filter_batch!(fast.image)
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        _compiled_paths_unpack_filter!(envelopes[slot], fast.image.field, path.phase_multiplier[a], slot)
    end
    return nothing
end

function _compiled_paths_prepare_pixel_mtf(plan, products)
    plan.apply_pixel_mtf || return nothing
    centered = Array(plan.pixel_mtf)
    all(iszero, imag.(centered)) || throw(ArgumentError(
        "compiled optical paths pixel MTF transfer must have real values"))
    # The artifact stores fftshift(fftfreq) ordering; AOS accepts corner-origin
    # FFT ordering. Conversion belongs to cold calibration preparation.
    corner = circshift(real.(centered),
        (-fld(size(centered, 1), 2), -fld(size(centered, 2), 2)))
    return map((products.fringed, products.unfringed, products.reference)) do input
        output = _compiled_paths_array(plan.target, eltype(input), size(input)...)
        AdaptiveOpticsSim.Detectors.prepare_pixel_mtf(input, corner;
            padding=plan.pixel_pad, output)
    end
end

_compiled_paths_apply_pixel_mtf!(::Nothing) = nothing

function _compiled_paths_apply_pixel_mtf!(owners::Tuple)
    for owner in owners
        output = AdaptiveOpticsSim.Detectors.apply_pixel_mtf!(owner)
        Backends.copyto_backend_async!(owner.input, output)
    end
    return nothing
end

function _enqueue_compiled_paths!(owner::PreparedCompiledPaths, share_prefix::Bool=true, batch_fast::Bool=true, prepared::Bool=true)
    plan, w, products = owner.plan, owner.workspace, owner.products
    _compiled_paths_slow!(owner, share_prefix, prepared)
    if batch_fast
        _compiled_paths_fast!(w.closed_envelopes, w.closed_sources, owner, true)
        _compiled_paths_fast!(w.only_envelopes, w.only_sources, owner, false)
    else
        _compiled_paths_fast_scalar!(w.closed_envelopes, w.closed_sources, owner, true)
        _compiled_paths_fast_scalar!(w.only_envelopes, w.only_sources, owner, false)
    end
    for k in 1:length(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        gain = plan.pinhole_gain[a]
        open = w.open_envelopes[slot]
        closed = w.closed_envelopes[slot]
        only = w.only_envelopes[slot]
        @. open = closed + gain * only
    end
    for synthesis in owner.synthesis
        AdaptiveOpticsSim.Optics.synthesize_spectrum!(synthesis)
    end
    _compiled_paths_apply_pixel_mtf!(w.pixel_mtf)
    return products
end

function _compiled_paths_selected!(owner, product::Val{:unfringed})
    _compiled_paths_slow!(owner)
    return _compiled_paths_selected_intensity!(owner, product)
end

function _compiled_paths_selected_intensity!(owner, ::Val{:unfringed})
    w = owner.workspace
    _compiled_paths_fast!(w.closed_envelopes, w.closed_sources, owner, true)
    AdaptiveOpticsSim.Optics.synthesize_spectrum!(owner.synthesis[1])
    return owner.products.unfringed
end

function _compiled_paths_selected!(owner, product::Val{:fringed})
    plan, w = owner.plan, owner.workspace
    _compiled_paths_slow!(owner)
    # Each subsequent shift, defocus and phase operation is linear and uses
    # the same operator for both fields. The free ghosts enter only once.
    for slot in eachindex(w.closed_sources)
        gain = plan.pinhole_gain[mod1(slot, plan.n_anchor)]
        closed = w.closed_sources[slot]
        only = w.only_sources[slot]
        @. closed = closed + gain * only
    end
    return _compiled_paths_selected_intensity!(owner, product)
end

function _compiled_paths_selected_intensity!(owner, ::Val{:fringed})
    w = owner.workspace
    _compiled_paths_fast!(w.open_envelopes, w.closed_sources, owner, true)
    AdaptiveOpticsSim.Optics.synthesize_spectrum!(owner.synthesis[2])
    return owner.products.fringed
end

include("compiled_paths/field_merge.jl")

"""Update the prepared canonical `(y,x)` fringed, unfringed and reference images."""
function AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner::PreparedCompiledPaths)
    owner.execution_context === nothing && throw(ArgumentError(
        "graph-owned compiled_paths propagation requires its caller's prepared execution context"))
    return Backends._with_prepared_device_execution_context(owner.execution_context) do
        products = _enqueue_compiled_paths!(owner)
        Backends._synchronize_prepared_device_execution_context!(owner.execution_context)
        products
    end
end


function _validate_compiled_paths_buffer(target,values,::Type{T},shape,role) where {T}
    eltype(values)===T || throw(ArgumentError("$role element type must be $T"))
    size(values)==shape || throw(DimensionMismatch("$role shape must be $shape"))
    Backends.compute_device(values)==target || throw(ArgumentError("$role occupies another compute device"))
    return values
end
