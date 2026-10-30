# Explicit external SPIDERS fitted-program qualification, not a production simulator fallback.
# HSA_OVERRIDE_GFX_VERSION=<local override> JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
# julia --project=<source environment> test/amdgpu_compiled_paths_fitted.jl /path/to/instrument
using Test, AdaptiveOpticsSim, Proper, AMDGPU, FITSIO, YAML, SHA, LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends
AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const FITTED_TARGET=Backends.AcceleratorComputeDevice(Backends.AMDGPUBackend(),1)

function fitted_input_graph(config,reference,input;target=Backends.HostComputeDevice(),
        execution=StreamGraphExecution())
    node=AlgorithmGraphs.compiled_optical_paths_node_factory(T=Float32)(:optics,config,(;),
        NamedTuple{(Symbol(config.entrance_reference_amplitude),)}((reference,)))
    definition=algorithm_graph((node,);name=:frozen_fitted_paths,
        inputs=(graph_input(:dm,:optics=>:dm_opd,input.dm),
            graph_input(:opd,:optics=>:pupil_opd,input.opd),
            graph_input(:amplitude,:optics=>:pupil_amplitude,input.amplitude)),
        outputs=(graph_output(:fringed,:optics=>:fringed),
            graph_output(:unfringed,:optics=>:unfringed),
            graph_output(:reference,:optics=>:reference)))
    return prepare_algorithm_graph(definition;target,execution)
end

function fitted_capture_check(instrument)
    graph_path=joinpath(instrument,"optics/generated/fitted/scc.yaml")
    graph_document=YAML.load_file(graph_path)
    entry=only(graph_document["nodes"])
    config=deepcopy(entry["config"])
    config["model_path"]=joinpath(instrument,config["model_path"])
    factory=AlgorithmGraphs.compiled_optical_paths_node_factory(T=Float32)
    config=AlgorithmGraphs.graph_node_configuration(factory,
        AlgorithmGraphs._normalize_graph_value(config),
        AlgorithmGraphs._normalize_graph_value(entry["bindings"]))
    descriptor=YAML.load_file(joinpath(instrument,"optics/generated/bindings.yaml"))[config.entrance_reference_amplitude]
    reference_path=joinpath(instrument,descriptor["path"])
    @test bytes2hex(sha256(read(reference_path)))==descriptor["sha256"]
    reference=FITS(reference_path,"r") do file
        permutedims(read(file[descriptor["hdu"]]),Tuple(descriptor["axes"]))
    end
    @test eltype(reference)===Float32
    @test size(reference)==Tuple(descriptor["shape"])==config.pupil_grid.size
    host=(dm=zeros(Float32,config.dm_grid.size),opd=zeros(Float32,config.pupil_grid.size),
        amplitude=copy(reference))
    device=map(ROCArray,host)
    println("Preparing CPU numerical reference from hash-verified frozen program");flush(stdout)
    cpu=fitted_input_graph(config,reference,host)
    println("Preparing one captured AMDGPU owner");flush(stdout)
    gpu=fitted_input_graph(config,reference,device;target=FITTED_TARGET,execution=CapturedGraphExecution())
    @test captured_graph_node_count(gpu)==1
    @test graph_execution_policy(gpu) isa CapturedGraphExecution
    owner=AlgorithmGraphs.prepared_graph_node(gpu,:optics)
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    @test Backends.compute_device(owner.propagation.workspace.field)==FITTED_TARGET
    retained=(field=owner.propagation.workspace.field,dm=owner.propagation.inputs.dm_opd,
        interpolation=owner.mapping,reference=owner.workspace.reference,flags=owner.workspace.validation.device)
    initial=nothing
    previous=nothing
    for frame in 1:4
        frame==2 && (host.dm .= Float32[3e-8sin(.15i+.09j) for i in axes(host.dm,1),j in axes(host.dm,2)])
        frame==3 && (host.opd .= Float32[2e-8cos(.07i-.11j) for i in axes(host.opd,1),j in axes(host.opd,2)])
        frame==4 && (host.amplitude .= reference .* Float32[.7+.2cos(.03i+.05j)^2
            for i in axes(reference,1),j in axes(reference,2)])
        for (destination,source) in zip(device,host);copyto!(destination,source);end
        AMDGPU.synchronize()
        step_graph!(cpu);step_graph!(gpu)
        actual=map(name->Array(graph_output(gpu,name)),(:fringed,:unfringed,:reference))
        for (name,values) in zip((:fringed,:unfringed,:reference),actual)
            expected=graph_output(cpu,name)
            relative=norm(values-expected)/norm(expected)
            peak=maximum(abs,values-expected)/maximum(abs,expected)
            @test all(isfinite,values) && sum(values)>0
            @test relative<=1e-3
            @test peak<=2e-3
            println("frame=",frame," product=",name," relative_l2=",relative," peak_relative_error=",peak)
        end
        if frame==1
            initial=map(copy,actual)
        else
            @test norm(first(actual)-previous)/norm(previous)>1e-6
        end
        previous=copy(first(actual))
        @test owner.propagation.workspace.field===retained.field
        @test owner.propagation.inputs.dm_opd===retained.dm
        @test owner.mapping===retained.interpolation
        @test owner.workspace.reference===retained.reference
        @test owner.workspace.validation.device===retained.flags
        flush(stdout)
    end
    # Measure only repeated retained-graph steps, after all initialization/JIT,
    # reference work and native capture. This is functional evidence, not an SLO.
    elapsed=[@elapsed(step_graph!(gpu)) for _ in 1:3]
    println("warm captured seconds=",elapsed);flush(stdout)
    fill!(host.dm,0);fill!(host.opd,0);copyto!(host.amplitude,reference)
    for (destination,source) in zip(device,host);copyto!(destination,source);end
    AMDGPU.synchronize()
    reset_graph!(gpu)
    for name in (:fringed,:unfringed,:reference)
        @test all(iszero,Array(graph_output(gpu,name)))
    end
    step_graph!(gpu)
    for (name,expected) in zip((:fringed,:unfringed,:reference),initial)
        actual=Array(graph_output(gpu,name))
        @test norm(actual-expected)/norm(expected)<=1e-6
    end
    # Validate the complete input grid, including a corner outside pupil support.
    host.opd[1,1]=NaN32
    copyto!(device.opd,host.opd);AMDGPU.synchronize()
    sequence=graph_step_sequence(gpu)
    ticket=step_graph_async!(gpu)
    @test_throws r"entrance OPD" wait_graph_step!(ticket)
    @test graph_step_sequence(gpu)==sequence
    @test graph_failed(gpu)
    host.opd[1,1]=0
    copyto!(device.opd,host.opd);AMDGPU.synchronize()
    reset_graph!(gpu);step_graph!(gpu)
    @test graph_step_sequence(gpu)==1
    @test all(isfinite,Array(graph_output(gpu,:fringed)))
    @test bytes2hex(sha256(read(config.model_path)))==config.model_sha256
    @test bytes2hex(sha256(read(reference_path)))==descriptor["sha256"]
    return nothing
end

length(ARGS)==1 || error("pass one frozen fitted instrument directory")
@testset "real frozen fitted compiled program on captured AMDGPU" begin
    fitted_capture_check(abspath(only(ARGS)))
end
@test !any(id->id.name=="CUDA",keys(Base.loaded_modules))
