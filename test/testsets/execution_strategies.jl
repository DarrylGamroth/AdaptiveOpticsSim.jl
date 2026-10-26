import KernelAbstractions

const PE02_REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function pe02_maintained_text()
    paths = String[]
    for directory in ("src", "ext", "test", "docs")
        root = joinpath(PE02_REPOSITORY_ROOT, directory)
        for (parent, _, files) in walkdir(root)
            for file in files
                extension = splitext(file)[2]
                extension in (".jl", ".toml", ".yaml", ".yml", ".md") || continue
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
    )
    selection_functions = (
        :atmospheric_field_execution_strategy,
        :reduction_execution_strategy,
        :detector_execution_strategy,
        :counting_output_execution_strategy,
        :grouped_accumulation_strategy,
    )
    @test length(family_types) == 5
    @test length(selector_types) == 12
    @test length(selection_functions) == 5
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
    )
    previous_function_stems = (
        "atmospheric_field_execution",
        "reduction_execution",
        "detector_execution",
        "counting_output_execution",
        "grouped_accumulation",
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
end
