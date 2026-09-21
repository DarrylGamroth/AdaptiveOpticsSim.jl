import KernelAbstractions

const PE02_REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

mutable struct EnsembleCounter{O}
    value::Int
    owner::O
end

struct ImmutableEnsembleMember{O}
    owner::O
end

struct EnsembleMemberError <: Exception end

EnsembleCounter(value::Int=0) = EnsembleCounter(value, Ref(value))

Ensembles.ensemble_ownership_roots(counter::EnsembleCounter) =
    (counter.owner,)
Ensembles.ensemble_ownership_roots(member::ImmutableEnsembleMember) =
    (member.owner,)

function increment_counter!(counter::EnsembleCounter)
    counter.value += 1
    return counter
end

function pe02_maintained_text()
    paths = String[]
    for directory in ("src", "ext", "test", "docs")
        root = joinpath(PE02_REPOSITORY_ROOT, directory)
        for (parent, _, files) in walkdir(root)
            for file in files
                extension = splitext(file)[2]
                extension in (".jl", ".toml", ".md") || continue
                push!(paths, joinpath(parent, file))
            end
        end
    end
    sort!(paths)
    return join((read(path, String) for path in paths), '\n')
end

@testset "PE-02 execution-strategy terminology" begin
    family_types = (
        (Atmospheres, :AbstractAtmosphericFieldExecutionStrategy),
        (Backends, :AbstractReductionExecutionStrategy),
        (Detectors, :AbstractDetectorExecutionStrategy),
        (WavefrontSensors, :AbstractGroupedAccumulationStrategy),
        (WavefrontSensors,
            :AbstractPyramidModulationPropagationStrategy),
        (WavefrontSensors, :AbstractShackHartmannWFSSensingStrategy),
    )
    selector_types = (
        (Atmospheres, :GeometricFieldSynchronousStrategy),
        (Atmospheres, :GeometricFieldAsyncStrategy),
        (Atmospheres, :LayeredFresnelFieldSynchronousStrategy),
        (Atmospheres, :LayeredFresnelFieldAsyncStrategy),
        (Backends, :DirectReductionStrategy),
        (Backends, :HostMirrorReductionStrategy),
        (Detectors, :DetectorDirectStrategy),
        (Detectors, :DetectorHostMirrorStrategy),
        (WavefrontSensors, :GroupedDirectAccumulateStrategy),
        (WavefrontSensors, :GroupedStackReduceStrategy),
        (WavefrontSensors, :GroupedStaged2DStrategy),
        (WavefrontSensors, :PyramidPupilTiltStrategy),
        (WavefrontSensors, :PyramidShiftedMaskStrategy),
        (WavefrontSensors, :ShackHartmannWFSScalarStrategy),
        (WavefrontSensors, :ShackHartmannWFSBatchedStrategy),
        (WavefrontSensors, :ShackHartmannWFSDeviceStatsStrategy),
        (WavefrontSensors, :ShackHartmannWFSROCmSafeStrategy),
        (WavefrontSensors, :ShackHartmannWFSROCmHostStatsStrategy),
    )
    selection_functions = (
        :atmospheric_field_execution_strategy,
        :reduction_execution_strategy,
        :detector_execution_strategy,
        :counting_output_execution_strategy,
        :grouped_accumulation_strategy,
        :sh_sensing_execution_strategy,
    )
    @test length(family_types) == 6
    @test length(selector_types) == 18
    @test length(selection_functions) == 6
    @test all(endswith("Strategy"), String(name) for (_, name) in family_types)
    @test all(endswith("Strategy"), String(name) for (_, name) in selector_types)
    @test all(name -> !endswith(String(name), "_plan"), selection_functions)
    @test all(isabstracttype(getfield(owner, name))
        for (owner, name) in family_types)
    @test all(isconcretetype(getfield(owner, name))
        for (owner, name) in selector_types)
    @test all(fieldcount(getfield(owner, name)) == 0
        for (owner, name) in selector_types)

    previous_suffix = string("Pl", "an")
    previous_type_stems = (
        "AbstractAtmosphericFieldExecution",
        "GeometricFieldSynchronous",
        "GeometricFieldAsync",
        "LayeredFresnelFieldSynchronous",
        "LayeredFresnelFieldAsync",
        "AbstractReductionExecution",
        "DirectReduction",
        "HostMirrorReduction",
        "AbstractDetectorExecution",
        "DetectorDirect",
        "DetectorHostMirror",
        "AbstractGroupedAccumulation",
        "GroupedDirectAccumulate",
        "GroupedStackReduce",
        "GroupedStaged2D",
        "AbstractShackHartmannWFSSensing",
        "ShackHartmannWFSScalar",
        "ShackHartmannWFSBatched",
        "ShackHartmannWFSDeviceStats",
        "ShackHartmannWFSRocmSafe",
        "ShackHartmannWFSRocmHostStats",
    )
    previous_function_stems = (
        "atmospheric_field_execution",
        "reduction_execution",
        "detector_execution",
        "counting_output_execution",
        "grouped_accumulation",
        "sh_sensing_execution",
        "sh_uses_rocm_safe_sensing",
        "sh_uses_host_stats_sensing",
        "sh_uses_batched_sensing",
        "sh_uses_device_stats_sensing",
        "_detector_value",
    )
    previous_names = (
        (string(stem, previous_suffix) for stem in previous_type_stems)...,
        (string(stem, '_', lowercase(previous_suffix))
            for stem in previous_function_stems)...,
    )
    maintained_text = pe02_maintained_text()
    @test all(name -> !occursin(name, maintained_text), previous_names)
end

@testset "PE-02 execution-strategy inference" begin
    scalar = Backends.ScalarCPUStyle()
    accelerated = Backends.AcceleratorStyle(KernelAbstractions.CPU())

    @test @inferred(Atmospheres.atmospheric_field_execution_strategy(
        scalar,
        Atmospheres.GeometricAtmosphericPropagation(),
    )) isa Atmospheres.GeometricFieldSynchronousStrategy
    @test @inferred(Atmospheres.atmospheric_field_execution_strategy(
        accelerated,
        Atmospheres.LayeredFresnelAtmosphericPropagation(),
    )) isa Atmospheres.LayeredFresnelFieldAsyncStrategy

    values = zeros(2, 2)
    @test @inferred(Backends.reduction_execution_strategy(
        scalar, values)) isa Backends.DirectReductionStrategy
    @test @inferred(Backends.reduction_execution_strategy(
        accelerated, values)) isa Backends.HostMirrorReductionStrategy

    detector = Detector(noise=NoiseNone())
    @test @inferred(Detectors.detector_execution_strategy(
        scalar, detector)) isa Detectors.DetectorDirectStrategy
    @test @inferred(Detectors.counting_output_execution_strategy(
        typeof(scalar),
        Detectors.SPADArrayDetector,
        Matrix{Float64},
    )) isa Detectors.DetectorDirectStrategy

    @test @inferred(WavefrontSensors.grouped_accumulation_strategy(
        typeof(scalar),
        WavefrontSensors.ShackHartmannWFS,
    )) isa WavefrontSensors.GroupedStackReduceStrategy
    @test @inferred(WavefrontSensors.sh_sensing_execution_strategy(
        typeof(scalar),
        WavefrontSensors.ShackHartmannWFS,
    )) isa WavefrontSensors.ShackHartmannWFSScalarStrategy
    @test @inferred(WavefrontSensors.sh_sensing_execution_strategy(
        typeof(accelerated),
        WavefrontSensors.ShackHartmannWFS,
    )) isa WavefrontSensors.ShackHartmannWFSBatchedStrategy
end

@testset "Generic coarse ensembles" begin
    first_counter = EnsembleCounter()
    second_counter = EnsembleCounter(10)
    ensemble = @inferred SimulationEnsemble(first_counter, second_counter)
    @test @inferred(ensemble_members(ensemble)) ===
        (first_counter, second_counter)
    @test @inferred(execution_policy(ensemble)) isa SequentialExecution
    @test @inferred(Ensembles.run_ensemble!(
        increment_counter!,
        ensemble,
    )) === ensemble
    @test (first_counter.value, second_counter.value) == (1, 11)
    if !coverage_instrumented()
        @test @allocated(Ensembles.run_ensemble!(
            increment_counter!, ensemble)) == 0
    end

    visit_order = Int[]
    ordered = SimulationEnsemble(
        EnsembleCounter(1),
        EnsembleCounter(2),
        EnsembleCounter(3),
    )
    Ensembles.run_ensemble!(
        member -> (push!(visit_order, member.value); member),
        ordered,
    )
    @test visit_order == [1, 2, 3]

    failing = SimulationEnsemble(
        EnsembleCounter(1),
        EnsembleCounter(2),
        EnsembleCounter(3),
    )
    @test_throws EnsembleMemberError Ensembles.run_ensemble!(
        member -> (member.value == 2 && throw(EnsembleMemberError()); member),
        failing,
    )
    @test map(member -> member.value, ensemble_members(failing)) == (1, 2, 3)

    immutable_values = SimulationEnsemble(1, 1)
    @test ensemble_members(immutable_values) === (1, 1)

    threaded = SimulationEnsemble(
        EnsembleCounter(),
        EnsembleCounter();
        policy=ThreadedExecution(),
    )
    Ensembles.run_ensemble!(increment_counter!, threaded)
    @test all(counter -> counter.value == 1, ensemble_members(threaded))

    shared_owner = Ref(0)
    @test_throws InvalidConfiguration SimulationEnsemble(
        EnsembleCounter(0, shared_owner),
        EnsembleCounter(0, shared_owner);
        policy=ThreadedExecution(),
    )
    @test_throws InvalidConfiguration SimulationEnsemble(
        ImmutableEnsembleMember(shared_owner),
        ImmutableEnsembleMember(shared_owner);
        policy=ThreadedExecution(),
    )
    @test_throws InvalidConfiguration SimulationEnsemble(())

    for policy in (
        BackendStreamExecution(),
        AcceleratedKernelsExecution(),
        DaggerExecution(),
    )
        unsupported = SimulationEnsemble(
            EnsembleCounter();
            policy=policy,
        )
        @test_throws UnsupportedAlgorithm Ensembles.run_ensemble!(
            increment_counter!,
            unsupported,
        )
    end

    if Threads.nthreads() == 1
        deterministic = SimulationEnsemble(
            EnsembleCounter();
            policy=DeterministicExecution(),
        )
        Ensembles.run_ensemble!(
            increment_counter!,
            deterministic,
        )
        @test first(ensemble_members(deterministic)).value == 1
        @test BLAS.get_num_threads() == 1
    else
        @test_throws InvalidConfiguration SimulationEnsemble(
            EnsembleCounter();
            policy=DeterministicExecution(),
        )
    end
end
