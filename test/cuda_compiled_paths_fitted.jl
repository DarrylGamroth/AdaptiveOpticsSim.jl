# Offline frozen fitted SCC gate; no core, simulator service or GUI.
# julia --startup-file=no -O2 --threads=1 --project=ENV test/cuda_compiled_paths_fitted.jl APPLICATION [--benchmark|--profile]
using Test, AdaptiveOpticsSim, Proper, CUDA, FITSIO, YAML, SHA, LinearAlgebra, Statistics
using AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const CP = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimProperExt)
CUDA.allowscalar(false)
BLAS.set_num_threads(1)
B.set_fft_provider_threads!(1)

function fitted_graph(application, name)
    path = joinpath(application, "optics/generated/fitted/$name.yaml")
    entry = only(YAML.load_file(path)["nodes"])
    raw = deepcopy(entry["config"])
    raw["model_path"] = joinpath(application, raw["model_path"])
    factory = AlgorithmGraphs.compiled_optical_paths_node_factory(T=Float32)
    config = AlgorithmGraphs.graph_node_configuration(factory,
        AlgorithmGraphs._normalize_graph_value(raw), AlgorithmGraphs._normalize_graph_value(entry["bindings"]))
    descriptor = YAML.load_file(joinpath(application, "optics/generated/bindings.yaml"))[config.entrance_reference_amplitude]
    reference_path = joinpath(application, descriptor["path"])
    @test bytes2hex(sha256(read(reference_path))) == descriptor["sha256"]
    reference = FITS(reference_path, "r") do file
        permutedims(read(file[descriptor["hdu"]]), Tuple(descriptor["axes"]))
    end
    host = (dm=zeros(Float32, config.dm_grid.size), opd=zeros(Float32, config.pupil_grid.size), amplitude=copy(reference))
    device = map(CuArray, host)
    node = factory(:optics, config, (;), NamedTuple{(Symbol(config.entrance_reference_amplitude),)}((reference,)))
    definition = algorithm_graph((node,); name=:fitted_scc,
        inputs=(graph_input(:dm, :optics=>:dm_opd, device.dm),
            graph_input(:opd, :optics=>:pupil_opd, device.opd),
            graph_input(:amplitude, :optics=>:pupil_amplitude, device.amplitude)),
        outputs=(graph_output(:fringed, :optics=>:fringed),
            graph_output(:unfringed, :optics=>:unfringed), graph_output(:reference, :optics=>:reference)))
    graph = prepare_algorithm_graph(definition; target=B.AcceleratorComputeDevice(B.CUDABackend(),0),
        execution=CapturedGraphExecution())
    return (; graph, host, device, reference, config)
end

function fitted_images(owner)
    return map(name -> Array(getproperty(owner.products, name)), (:fringed,:unfringed,:reference))
end
function capture_paths(graph, operation)
    return B._with_prepared_device_execution_context(graph.context) do
        Proper.prepare_cuda_graph(operation; warmup=2)
    end
end
function replay_paths(graph, captured)
    B._with_prepared_device_execution_context(graph.context) do
        Proper.launch_cuda_graph!(captured)
        CUDA.synchronize()
    end
end
function node_count(graph)
    count=Ref{Csize_t}(0)
    CUDA.CUDACore.cuGraphGetNodes(graph,C_NULL,count)
    return count[]
end
function compare_fitted(application, name; benchmark=false, profile=false)
    preparation_start=time_ns()
    fixture = fitted_graph(application, name)
    preparation_seconds=(time_ns()-preparation_start)*1e-9
    (;graph,host,device,reference,config) = fixture
    owner = AlgorithmGraphs.prepared_graph_node(graph,:optics).propagation
    retained = (owner.workspace.field, owner.inputs, owner.workspace.propagation)
    phase_arrays=Base.IdSet{typeof(owner.workspace.field)}()
    for propagation in owner.workspace.propagation
        propagation === nothing && continue
        for phase in (propagation.before,propagation.transfer,propagation.after)
            stored=Proper._propagation_phase_storage(phase)
            stored === nothing || isempty(stored) || push!(phase_arrays,stored)
        end
    end
    phase_bytes=sum(length(phase)*sizeof(eltype(phase)) for phase in phase_arrays)
    println((;name,phase_arrays=length(phase_arrays),phase_bytes,
        shared_fields=length(owner.workspace.shared_prefixes),
        shared_field_bytes=sum(length(field)*sizeof(eltype(field)) for field in owner.workspace.shared_prefixes)));flush(stdout)
    # Both graphs reuse the exact same owner, FFT plans, buffers and inputs.
    capture_start=time_ns()
    ordinary = capture_paths(graph, () -> CP._enqueue_compiled_paths!(owner,true,true,false))
    original = capture_paths(graph, () -> begin
        CP._enqueue_compiled_paths!(owner,false,false,false)
        for synthesis in owner.synthesis
            AdaptiveOpticsSim.Optics._synthesize_spectrum_reference!(synthesis)
        end
        CP._compiled_paths_apply_pixel_mtf!(owner.workspace.pixel_mtf)
    end)
    prepared = capture_paths(graph, () -> CP._enqueue_compiled_paths!(owner))
    capture_seconds=(time_ns()-capture_start)*1e-9
    first_images = nothing
    previous = nothing
    for state in 1:3
        host.dm .= Float32[3e-8sin(state*(.15i+.09j)) for i in axes(host.dm,1), j in axes(host.dm,2)]
        host.opd .= Float32[2e-8cos(state*(.07i-.11j)) for i in axes(host.opd,1), j in axes(host.opd,2)]
        host.amplitude .= reference .* Float32[.7+.2cos(state*(.03i+.05j))^2
            for i in axes(reference,1), j in axes(reference,2)]
        foreach(pair -> copyto!(pair...), zip(device,host))
        CUDA.synchronize()
        step_graph!(graph)
        actual = fitted_images(owner)
        replay_paths(graph, original)
        expected = fitted_images(owner)
        for (product,a,b) in zip((:fringed,:unfringed,:reference),actual,expected)
            relative=norm(a-b)/norm(b)
            peak=maximum(abs,a-b)/maximum(abs,b)
            println((;name,state,product,relative,peak,exact=a==b));flush(stdout)
            @test all(isfinite,a) && sum(a)>0
            @test relative <= 1e-5 && peak <= 1e-5
        end
        if state==1
            first_images=map(copy,actual)
        else
            @test norm(first(actual)-previous)/norm(previous)>1e-6
        end
        previous=copy(first(actual))
        @test owner.workspace.field === retained[1] && owner.inputs === retained[2] && owner.workspace.propagation === retained[3]
    end
    # Reset and captured replay must reproduce the first changed-input state.
    host.dm .= Float32[3e-8sin(.15i+.09j) for i in axes(host.dm,1), j in axes(host.dm,2)]
    host.opd .= Float32[2e-8cos(.07i-.11j) for i in axes(host.opd,1), j in axes(host.opd,2)]
    host.amplitude .= reference .* Float32[.7+.2cos(.03i+.05j)^2 for i in axes(reference,1), j in axes(reference,2)]
    foreach(pair -> copyto!(pair...), zip(device,host))
    CUDA.synchronize()
    reset_graph!(graph)
    @test all(name -> all(iszero,Array(graph_output(graph,name))),(:fringed,:unfringed,:reference))
    step_graph!(graph)
    @test all(pair -> pair[1]==pair[2],zip(fitted_images(owner),first_images))
    @test bytes2hex(sha256(read(config.model_path)))==config.model_sha256
    if benchmark
        times = (Float64[], Float64[])
        # Alternate blocks in one process; exclude capture/JIT/preparation.
        for _ in 1:8, (index,captured) in enumerate((ordinary,prepared))
            for _ in 1:4; replay_paths(graph,captured); end
            for _ in 1:8; push!(times[index], @elapsed replay_paths(graph,captured)); end
        end
        println((;name,preparation_seconds,capture_seconds,julia=string(VERSION),
            device=CUDA.name(CUDA.device()),baseline_nodes=node_count(ordinary.graph),prepared_nodes=node_count(prepared.graph),
            baseline_ms=1000median(times[1]),prepared_ms=1000median(times[2]),
            baseline_p90_ms=1000quantile(times[1],.9),prepared_p90_ms=1000quantile(times[2],.9),
            baseline_max_ms=1000maximum(times[1]),prepared_max_ms=1000maximum(times[2]),
            samples_per_case=length(times[1])));flush(stdout)
    end
    if profile
        CUDA.@profile external=true begin
            for _ in 1:4; replay_paths(graph,prepared); end
        end
    end
end

isempty(ARGS) && error("pass frozen fitted application directory")
@testset "frozen fitted SCC captured propagation" begin
    for name in ("scc","scc_unfringed")
        compare_fitted(abspath(first(ARGS)),name;benchmark="--benchmark" in ARGS,
            profile="--profile" in ARGS && name=="scc")
        GC.gc();CUDA.reclaim()
    end
end
