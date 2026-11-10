function _compiled_paths_execute_scheduled_group!(owner, sequence, offset, step, size::Val)
    input = step.cache == 0 ? owner.workspace.field : owner.workspace.shared_prefixes[step.cache]
    mask = step.mask == 0 ? nothing : sequence[step.mask].array
    _compiled_paths_execute_group!(owner,sequence,offset,step.index,size,input,mask)
    return nothing
end

# One bounded opcode selection per prepared action. Group arities enter the
# numerical executor as literal Val parameters; there is no topology search,
# dictionary, arbitrary callable registry or model-specific generated method.
function _compiled_paths_execute_schedule!(owner, sequence, offset, schedule, steps=owner.execution.steps)
    w = owner.workspace
    for position in schedule.first:schedule.last
        step = steps[position]
        kind = step.kind
        if kind == LoadCheckpoint
            Backends.copyto_backend_async!(w.field,w.shared_prefixes[step.cache])
            continue
        elseif kind == StoreCheckpoint
            Backends.copyto_backend_async!(w.shared_prefixes[step.cache],w.field)
            continue
        end
        index = step.index
        op = sequence[index]
        work = w.operations[offset+index]
        if kind == ExecuteMultiply
            w.field .*= op.array
        elseif kind == ExecuteNormalize
            _compiled_paths_normalize!(w.field,w.normalization)
        elseif kind == ExecutePropagation
            _compiled_paths_proper!(w,op,w.propagation[offset+index])
        elseif kind == ExecutePropagationPair
            _compiled_paths_proper_pair!(w,op,sequence[index+1],
                w.propagation[offset+index],w.propagation[offset+index+1])
        elseif kind == ExecuteGroup10
            _compiled_paths_execute_scheduled_group!(owner,sequence,offset,step,Val(10))
        elseif kind == ExecuteGroup8
            _compiled_paths_execute_scheduled_group!(owner,sequence,offset,step,Val(8))
        elseif kind == ExecuteGroup6
            _compiled_paths_execute_scheduled_group!(owner,sequence,offset,step,Val(6))
        elseif kind == ExecuteGroup4
            _compiled_paths_execute_scheduled_group!(owner,sequence,offset,step,Val(4))
        elseif kind == ExecuteGroup2
            _compiled_paths_execute_scheduled_group!(owner,sequence,offset,step,Val(2))
        elseif kind == ExecutePropagationResample
            _compiled_paths_proper_resample!(w,op,sequence[index+1],sequence[index+2],
                w.operations[offset+index+2],w.propagation[offset+index],
                w.propagation[offset+index+1],Val(false))
        elseif kind == ExecuteCancelledResample
            _compiled_paths_proper_resample!(w,op,sequence[index+1],sequence[index+2],
                w.operations[offset+index+2],w.propagation[offset+index],
                w.propagation[offset+index+1],Val(true))
        elseif kind == ExecuteSpectralMultiplier
            _compiled_paths_filter!(w.field,w.field,op.array,w.physical_fft)
        elseif kind == ExecuteResample
            _compiled_paths_resample!(w,op,work)
        elseif kind == ExecuteLiveDM
            _compiled_paths_interpolate!(work,op.array,owner.inputs.dm_opd,w.prepared_mul_coefficients)
            AdaptiveOpticsSim.Optics.apply_phase_async!(work.electric_field,work.opd)
        elseif kind == ExecuteLiveEntrance
            _compiled_paths_interpolate!(work,op.array,owner.inputs.entrance_opd,w.prepared_mul_coefficients)
            AdaptiveOpticsSim.Optics.apply_phase_async!(work.electric_field,work.opd)
            _compiled_paths_interpolate!(work,op.array,owner.inputs.entrance_attenuation,w.prepared_mul_coefficients)
            @. w.field *= exp(-work.opd)
        else
            throw(ArgumentError("unknown compiled execution step"))
        end
    end
    return nothing
end
