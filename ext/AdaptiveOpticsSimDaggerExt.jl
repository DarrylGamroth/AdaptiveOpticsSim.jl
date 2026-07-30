module AdaptiveOpticsSimDaggerExt

using AdaptiveOpticsSim
import AdaptiveOpticsSim: Ensembles
import Dagger

struct DaggerSchedulerState end

Ensembles.init_ensemble_scheduler(
    ::Ensembles.DaggerExecution,
    ::Tuple,
) = DaggerSchedulerState()

@inline function _spawn_member(
    ::Ensembles.DaggerExecution{Nothing},
    f!,
    member,
)
    return Dagger.spawn(Ensembles._run_ensemble_member!, f!, member)
end

@inline function _spawn_member(
    policy::Ensembles.DaggerExecution,
    f!,
    member,
)
    options = Dagger.Options(; scope=policy.scope)
    return Dagger.spawn(
        Ensembles._run_ensemble_member!,
        options,
        f!,
        member,
    )
end

function Ensembles.execute_ensemble!(
    policy::Ensembles.DaggerExecution,
    ::DaggerSchedulerState,
    f!,
    members::Tuple,
)
    tasks = map(member -> _spawn_member(policy, f!, member), members)
    return map(fetch, tasks)
end

end
