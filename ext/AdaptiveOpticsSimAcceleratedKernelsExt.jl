module AdaptiveOpticsSimAcceleratedKernelsExt

using AdaptiveOpticsSim
import AcceleratedKernels

struct AcceleratedKernelsSchedulerState{TP}
    partitioner::TP
end

function AdaptiveOpticsSim.init_ensemble_scheduler(
    policy::AdaptiveOpticsSim.AcceleratedKernelsExecution,
    members::Tuple,
)
    partitioner = AcceleratedKernels.TaskPartitioner(
        length(members),
        policy.max_tasks,
        policy.min_members_per_task,
    )
    return AcceleratedKernelsSchedulerState(partitioner)
end

function AdaptiveOpticsSim.execute_ensemble!(
    ::AdaptiveOpticsSim.AcceleratedKernelsExecution,
    state::AcceleratedKernelsSchedulerState,
    f!,
    members::Tuple,
)
    if length(state.partitioner) == 1
        @inbounds for i in eachindex(members)
            AdaptiveOpticsSim._run_ensemble_member!(f!, members[i])
        end
        return members
    end
    AcceleratedKernels.task_partition(state.partitioner) do member_indices
        @inbounds for i in member_indices
            AdaptiveOpticsSim._run_ensemble_member!(f!, members[i])
        end
    end
    return members
end

end
