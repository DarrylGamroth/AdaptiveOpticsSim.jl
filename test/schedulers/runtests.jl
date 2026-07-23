using AcceleratedKernels
using AdaptiveOpticsSim
using Dagger
using Test

mutable struct SchedulerCounter
    count::Int
end

increment!(counter::SchedulerCounter) =
    (counter.count += 1; counter)
increment_ten!(counter::SchedulerCounter) =
    (counter.count += 10; counter)

@testset "AcceleratedKernels ensemble extension" begin
    @test !isnothing(Base.get_extension(
        AdaptiveOpticsSim,
        :AdaptiveOpticsSimAcceleratedKernelsExt,
    ))
    counters = ntuple(_ -> SchedulerCounter(0), 8)
    ensemble = SimulationEnsemble(counters;
        policy=AcceleratedKernelsExecution(
            max_tasks=min(4, Threads.nthreads()),
            min_members_per_task=1,
        ))
    @test AdaptiveOpticsSim.run_ensemble!(increment!, ensemble) === ensemble
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(ensemble))
    @test AdaptiveOpticsSim.run_ensemble!(increment_ten!, ensemble) ===
        ensemble
    @test all(counter -> counter.count == 11,
        AdaptiveOpticsSim.ensemble_members(ensemble))

    single_task = SimulationEnsemble(
        ntuple(_ -> SchedulerCounter(0), 2);
        policy=AcceleratedKernelsExecution(
            max_tasks=1,
            min_members_per_task=1,
        ),
    )
    @test AdaptiveOpticsSim.run_ensemble!(increment!, single_task) ===
        single_task
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(single_task))
end

@testset "Dagger ensemble extension" begin
    @test !isnothing(Base.get_extension(
        AdaptiveOpticsSim,
        :AdaptiveOpticsSimDaggerExt,
    ))
    counters = ntuple(_ -> SchedulerCounter(0), 4)
    ensemble = SimulationEnsemble(counters; policy=DaggerExecution())
    @test AdaptiveOpticsSim.run_ensemble!(increment!, ensemble) === ensemble
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(ensemble))
    @test AdaptiveOpticsSim.run_ensemble!(increment_ten!, ensemble) ===
        ensemble
    @test all(counter -> counter.count == 11,
        AdaptiveOpticsSim.ensemble_members(ensemble))

    local_scope = Dagger.scope(worker=1)
    scoped = SimulationEnsemble(
        ntuple(_ -> SchedulerCounter(0), 2);
        policy=DaggerExecution(scope=local_scope),
    )
    @test AdaptiveOpticsSim.run_ensemble!(increment!, scoped) === scoped
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(scoped))
end
