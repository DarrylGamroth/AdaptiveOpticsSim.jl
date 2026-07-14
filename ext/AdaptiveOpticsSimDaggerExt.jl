module AdaptiveOpticsSimDaggerExt

using AdaptiveOpticsSim
import Dagger

struct DaggerSchedulerState end

AdaptiveOpticsSim.init_ensemble_scheduler(
    ::AdaptiveOpticsSim.DaggerExecution,
    ::Tuple,
) = DaggerSchedulerState()

@inline function _spawn_member(
    ::AdaptiveOpticsSim.DaggerExecution{Nothing},
    f!,
    member,
)
    return Dagger.spawn(AdaptiveOpticsSim._run_ensemble_member!, f!, member)
end

@inline function _spawn_member(
    policy::AdaptiveOpticsSim.DaggerExecution,
    f!,
    member,
)
    options = Dagger.Options(; scope=policy.scope)
    return Dagger.spawn(
        AdaptiveOpticsSim._run_ensemble_member!,
        options,
        f!,
        member,
    )
end

function AdaptiveOpticsSim.execute_ensemble!(
    policy::AdaptiveOpticsSim.DaggerExecution,
    ::DaggerSchedulerState,
    f!,
    members::Tuple,
)
    tasks = map(member -> _spawn_member(policy, f!, member), members)
    return map(fetch, tasks)
end

end
