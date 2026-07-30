module AdaptiveOpticsSimAcceleratedKernelsExt

using AdaptiveOpticsSim
import AdaptiveOpticsSim: Ensembles
import AcceleratedKernels

struct AcceleratedKernelsSchedulerState{TP}
    partitioner::TP
end

function Ensembles.init_ensemble_scheduler(
    policy::Ensembles.AcceleratedKernelsExecution,
    members::Tuple,
)
    partitioner = AcceleratedKernels.TaskPartitioner(
        length(members),
        policy.max_tasks,
        policy.min_members_per_task,
    )
    return AcceleratedKernelsSchedulerState(partitioner)
end

function Ensembles.execute_ensemble!(
    ::Ensembles.AcceleratedKernelsExecution,
    state::AcceleratedKernelsSchedulerState,
    f!,
    members::Tuple,
)
    if length(state.partitioner) == 1
        @inbounds for i in eachindex(members)
            Ensembles._run_ensemble_member!(f!, members[i])
        end
        return members
    end
    AcceleratedKernels.task_partition(state.partitioner) do member_indices
        @inbounds for i in member_indices
            Ensembles._run_ensemble_member!(f!, members[i])
        end
    end
    return members
end

end
