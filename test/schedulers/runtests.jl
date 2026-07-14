using AcceleratedKernels
using AdaptiveOpticsSim
using Dagger
using Test

mutable struct SchedulerCounter <: AbstractControlSimulation
    count::Int
end

AdaptiveOpticsSim.step!(counter::SchedulerCounter) =
    (counter.count += 1; counter)
AdaptiveOpticsSim.sense!(counter::SchedulerCounter) =
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
    @test step!(ensemble) === ensemble
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(ensemble))
    @test sense!(ensemble) === ensemble
    @test all(counter -> counter.count == 11,
        AdaptiveOpticsSim.ensemble_members(ensemble))
end

@testset "Dagger ensemble extension" begin
    @test !isnothing(Base.get_extension(
        AdaptiveOpticsSim,
        :AdaptiveOpticsSimDaggerExt,
    ))
    counters = ntuple(_ -> SchedulerCounter(0), 4)
    ensemble = SimulationEnsemble(counters; policy=DaggerExecution())
    @test step!(ensemble) === ensemble
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(ensemble))
    @test sense!(ensemble) === ensemble
    @test all(counter -> counter.count == 11,
        AdaptiveOpticsSim.ensemble_members(ensemble))

    local_scope = Dagger.scope(worker=1)
    scoped = SimulationEnsemble(
        ntuple(_ -> SchedulerCounter(0), 2);
        policy=DaggerExecution(scope=local_scope),
    )
    @test step!(scoped) === scoped
    @test all(counter -> counter.count == 1,
        AdaptiveOpticsSim.ensemble_members(scoped))
end
