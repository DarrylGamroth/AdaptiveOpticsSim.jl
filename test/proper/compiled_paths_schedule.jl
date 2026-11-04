@testset "Prepared path schedules retain fixed topology and bounded storage" begin
    fixture=compiled_paths_synthetic_fixture()
    declaration=only(fixture.metadata["propagations"])
    original=copy(declaration["prefix"])
    declaration["prefix"]=copy(original[1:3])
    for _ in 1:6
        append!(declaration["prefix"],[Dict("kind"=>"mult","array"=>"mask"),
            deepcopy(original[4]),deepcopy(original[5])])
    end
    owner=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
        key->fixture.arrays[key])
    w=owner.workspace
    propagation=CP._compiled_paths_prepare_propagation(AdaptiveOpticsSim.Backends.CUDABackend(),
        owner.plan.paths,w.wavefront,w.context)
    execution=CP._compiled_paths_prepare_execution(owner.plan,propagation)
    range=execution.unshared[1]
    actions=execution.steps[range.first:range.last]
    @test [step.kind for step in actions]==[CP.ExecuteLiveEntrance,CP.ExecuteNormalize,
        CP.ExecuteLiveDM,CP.ExecuteGroup10,CP.ExecuteGroup2]
    @test [step.mask for step in actions[4:5]]==[4,19]
    @test isconcretetype(typeof(execution))
    @test isbitstype(eltype(execution.steps))
    @test all(isconcretetype,fieldtypes(eltype(execution.steps)))

    # A producer checkpoint after the first pair forbids grouping across it.
    # A consumer can fuse its initial cache read only with prepared propagation.
    sequence=only(owner.plan.paths).prefix
    producer=CP.PathPrefixPlan(0,0,CP._compiled_paths_seal([CP.PathPrefixCheckpoint(6,1)]))
    steps=CP.PathExecutionStep[]
    CP._compiled_paths_prepare_sequence!(steps,sequence,propagation,0,producer)
    @test [step.kind for step in steps[4:6]]==[CP.ExecuteGroup2,CP.StoreCheckpoint,CP.ExecuteGroup10]
    consumer=CP.PathPrefixPlan(1,3,CP._compiled_paths_seal(CP.PathPrefixCheckpoint[]))
    empty!(steps)
    CP._compiled_paths_prepare_sequence!(steps,sequence,propagation,0,consumer)
    @test first(steps).kind==CP.ExecuteGroup10 && first(steps).cache==1
    empty!(steps)
    CP._compiled_paths_prepare_sequence!(steps,sequence,w.propagation,0,consumer)
    @test first(steps).kind==CP.LoadCheckpoint
    @test steps[2].kind==CP.ExecuteMultiply

    # Instruction storage types are independent of geometry, anchor values
    # and program length. Long user prescriptions grow data, not tuple types.
    for n in (7,16), count in (1,40)
        other=compiled_paths_synthetic_fixture(;n)
        prefix=only(other.metadata["propagations"])["prefix"]
        append!(prefix,[deepcopy(prefix[2]) for _ in 1:count])
        candidate=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(other.metadata,
            key->other.arrays[key])
        @test typeof(candidate.execution)===typeof(execution)
    end
end

@testset "Terminal resampling decisions respect checkpoint materialization" begin
    for distance in (-.004,0.,.004)
        fixture=compiled_paths_resample_fixture(;distance)
        owner=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
            key->fixture.arrays[key])
        w=owner.workspace
        propagation=CP._compiled_paths_prepare_propagation(AdaptiveOpticsSim.Backends.CUDABackend(),
            owner.plan.paths,w.wavefront,w.context)
        sequence=only(owner.plan.paths).prefix
        steps=CP.PathExecutionStep[]
        CP._compiled_paths_prepare_range!(steps,sequence,propagation,0,5,7)
        @test only(steps).kind==(distance==0 ? CP.ExecutePropagationResample : CP.ExecuteCancelledResample)
        empty!(steps)
        CP._compiled_paths_prepare_range!(steps,sequence,propagation,0,5,6)
        CP._compiled_paths_prepare_range!(steps,sequence,propagation,0,7,7)
        @test [step.kind for step in steps]==[CP.ExecutePropagationPair,CP.ExecuteResample]
    end
end
