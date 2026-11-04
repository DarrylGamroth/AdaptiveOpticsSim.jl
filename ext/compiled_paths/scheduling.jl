# Resolve topology once, after backend-qualified propagation has been prepared.
# The closed instruction set bounds specialization independently of path length,
# spectral anchor count, wavelength values and checkpoint indices.
const _COMPILED_PATH_GROUP_SIZES = (10, 8, 6, 4, 2)

function _compiled_paths_group_matches(sequence,index,last,::Val{N}) where {N}
    stop=index+3(N÷2)-2
    stop<=last || return false
    # Leave a final PTP and its Fourier resampling to the transform pair
    # lowering below. A checkpoint beyond `last` must still receive a field.
    stop<last && sequence[stop].kind==PathTransferFunction &&
        sequence[stop+1].kind==PathResample && return false
    for pair in 0:N÷2-1
        first=index+3pair
        sequence[first].kind in (PathTransferFunction,PathFourierTransform) &&
            sequence[first+1].kind in (PathTransferFunction,PathFourierTransform) || return false
        pair==N÷2-1 || sequence[first+2].kind==PathMultiply || return false
    end
    return true
end

function _compiled_paths_group_size(sequence, index, last; include_pair=true)
    for size in _COMPILED_PATH_GROUP_SIZES
        size == 2 && !include_pair && break
        _compiled_paths_group_matches(sequence,index,last,Val(size)) && return size
    end
    return 0
end

function _compiled_paths_group_kind(size)
    size == 10 && return ExecuteGroup10
    size == 8 && return ExecuteGroup8
    size == 6 && return ExecuteGroup6
    size == 4 && return ExecuteGroup4
    size == 2 && return ExecuteGroup2
    throw(ArgumentError("unsupported compiled propagation group size"))
end

function _compiled_paths_prepare_multiplied_group!(steps, sequence, propagation,
        offset, first, last, cache)
    first+2 <= last && sequence[first].kind == PathMultiply &&
        sequence[first+1].kind in (PathTransferFunction,PathFourierTransform) &&
        propagation[offset+first+1] !== nothing || return first-1
    size = _compiled_paths_group_size(sequence,first+1,last)
    size == 0 && return first-1
    push!(steps,PathExecutionStep(_compiled_paths_group_kind(size),first+1,first,cache))
    return first+3(size÷2)-1
end

function _compiled_paths_prepare_range!(steps, sequence, propagation, offset, first, last)
    index = first
    while index <= last
        op = sequence[index]
        if op.kind == PathMultiply
            consumed = _compiled_paths_prepare_multiplied_group!(steps,sequence,
                propagation,offset,index,last,0)
            if consumed >= index
                index = consumed+1
                continue
            end
            kind = ExecuteMultiply
        elseif op.kind == PathNormalize
            kind = ExecuteNormalize
        elseif op.kind in (PathTransferFunction,PathFourierTransform)
            # A terminal transform/PTP/resample has priority over a longer group.
            # `last` is a checkpoint boundary: cancellation cannot cross it.
            if index+2 <= last && sequence[index+1].kind == PathTransferFunction &&
                    sequence[index+2].kind == PathResample
                tail = propagation[offset+index+1]
                kind = tail !== nothing && tail.direction == 2 ?
                    ExecuteCancelledResample : ExecutePropagationResample
                push!(steps,PathExecutionStep(kind,index,0,0))
                index += 3
                continue
            end
            size = _compiled_paths_group_size(sequence,index,last;include_pair=false)
            if size > 0
                push!(steps,PathExecutionStep(_compiled_paths_group_kind(size),index,0,0))
                index += 3(size÷2)-1
                continue
            elseif index < last && sequence[index+1].kind in (PathTransferFunction,PathFourierTransform)
                push!(steps,PathExecutionStep(ExecutePropagationPair,index,0,0))
                index += 2
                continue
            end
            kind = ExecutePropagation
        elseif op.kind == PathSpectralMultiplier
            kind = ExecuteSpectralMultiplier
        elseif op.kind == PathResample
            kind = ExecuteResample
        elseif op.kind == PathLiveScreen
            kind = op.source == 1 ? ExecuteLiveDM : ExecuteLiveEntrance
        else
            throw(ArgumentError("unknown compiled optical operation"))
        end
        push!(steps,PathExecutionStep(kind,index,0,0))
        index += 1
    end
    return nothing
end

function _compiled_paths_prepare_sequence!(steps, sequence, propagation, offset, sharing)
    start = length(steps)+1
    first = sharing.count+1
    if sharing.count > 0
        last = isempty(sharing.checkpoints) ? length(sequence) : sharing.checkpoints[1].index
        consumed = _compiled_paths_prepare_multiplied_group!(steps,sequence,
            propagation,offset,first,last,sharing.cache_index)
        if consumed < first
            push!(steps,PathExecutionStep(LoadCheckpoint,0,0,sharing.cache_index))
        else
            first = consumed+1
        end
    end
    for checkpoint in sharing.checkpoints
        _compiled_paths_prepare_range!(steps,sequence,propagation,offset,first,checkpoint.index)
        push!(steps,PathExecutionStep(StoreCheckpoint,0,0,checkpoint.cache_index))
        first = checkpoint.index+1
    end
    _compiled_paths_prepare_range!(steps,sequence,propagation,offset,first,length(sequence))
    return PathSequenceExecution(start,length(steps))
end

function _compiled_paths_prepare_execution(plan, propagation)
    steps = PathExecutionStep[]
    shared = PathSequenceExecution[]
    unshared = PathSequenceExecution[]
    offset = 0
    for (index,path) in enumerate(plan.paths)
        sharing = plan.shared_prefixes[index]
        for (sequence,reuse) in zip((path.prefix,path.closed,path.only),
                (sharing.prefix,sharing.closed,sharing.only))
            start = length(steps)+1
            _compiled_paths_prepare_range!(steps,sequence,propagation,offset,1,length(sequence))
            push!(unshared,PathSequenceExecution(start,length(steps)))
            push!(shared,_compiled_paths_prepare_sequence!(steps,sequence,propagation,offset,reuse))
            offset += length(sequence)
        end
    end
    return PathExecutionPlan(_compiled_paths_seal(steps),
        _compiled_paths_seal(shared),_compiled_paths_seal(unshared))
end
