# Offline frozen LLOWFS model gate. No installed simulator, GUI or private core.
# julia --startup-file=no -O2 --threads=1 --project=ENV test/cuda_proper_train_fitted.jl APP [--benchmark|--profile]
using Test,AdaptiveOpticsSim,Proper,CUDA,FITSIO,YAML,SHA,LinearAlgebra,Statistics
using AdaptiveOpticsSim.AlgorithmGraphs
const TRAIN_BACKEND=AdaptiveOpticsSim.Backends
const TRAIN_EXTENSION=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
CUDA.allowscalar(false);BLAS.set_num_threads(1);TRAIN_BACKEND.set_fft_provider_threads!(1)

function train_resolve(value,document)
    if value isa AbstractDict
        if haskey(value,"\$ref")
            return foldl((entry,key)->entry[key],split(value["\$ref"],'/';keepempty=false);init=document)
        end
        return Dict(key=>train_resolve(item,document) for (key,item) in value)
    elseif value isa AbstractVector
        return map(item->train_resolve(item,document),value)
    end
    return value
end
function fitted_train_graph(application)
    descriptions=YAML.load_file(joinpath(application,"optics/generated/bindings.yaml"))
    function read_binding(name)
        descriptor=descriptions[name];path=joinpath(application,descriptor["path"])
        @test bytes2hex(sha256(read(path)))==descriptor["sha256"]
        FITS(path,"r") do file
            permutedims(read(file[descriptor["hdu"]]),Tuple(descriptor["axes"]))
        end
    end
    reference=read_binding("pupil_amplitude")
    static=(fpm_amplitude=read_binding("fpm_amplitude"),fpm_opd=read_binding("fpm_opd"),
        lyot_reflected=read_binding("lyot_reflected"),apodizer=read_binding("apodizer"))
    host=(pupil_opd=zeros(Float32,size(reference)),pupil_amplitude=copy(reference),dm_opd=zeros(Float32,128,128))
    device=map(CuArray,host)
    instrument=Dict("instrument"=>YAML.load_file(joinpath(application,"instrument.yaml")))
    path=joinpath(application,"optics/generated/llowfs-camera.yaml")
    document=train_resolve(YAML.load_file(path),instrument)
    definition=AlgorithmGraphs.load_algorithm_graph(document;bindings=merge(static,device))
    graph=prepare_algorithm_graph(definition;target=TRAIN_BACKEND.AcceleratorComputeDevice(
        TRAIN_BACKEND.CUDABackend(),0),execution=CapturedGraphExecution())
    return (;graph,host,device,reference,path,identity=bytes2hex(sha256(read(path))))
end
function train_capture(graph,operation)
    TRAIN_BACKEND._with_prepared_device_execution_context(graph.context) do
        Proper.prepare_cuda_graph(operation;warmup=2)
    end
end
function train_replay(graph,captured)
    TRAIN_BACKEND._with_prepared_device_execution_context(graph.context) do
        Proper.launch_cuda_graph!(captured);CUDA.synchronize()
        AlgorithmGraphs._complete_captured_graph_node!(AlgorithmGraphs.prepared_graph_node(graph,:camera_pixels))
    end
end
function train_node_count(graph)
    count=Ref{Csize_t}(0);CUDA.CUDACore.cuGraphGetNodes(graph,C_NULL,count);count[]
end
function train_storage(owner)
    execution=owner.train.workspace.execution;field=owner.train.workspace.wavefront.field
    storage=Base.IdSet{typeof(field)}()
    for segment in execution.segments
        segment isa TRAIN_EXTENSION.TrainPropagationSegment || continue
        for p in segment.propagation,phase in (p.before,p.transfer,p.after)
            values=Proper._propagation_phase_storage(phase)
            values===nothing || isempty(values) || push!(storage,values)
        end
        for masks in (segment.multipliers...,segment.tail),mask in masks
            push!(storage,parent(mask))
        end
    end
    return (;arrays=length(storage),cache_bytes=sum(sizeof,storage)+size(field,1)*sizeof(Int32),
        extra_field_bytes=0,field_is_fft_scratch=field===Proper.fft_workspace(owner.train.workspace.context).scratch)
end
function compare_fitted_train(application;benchmark=false,profile=false,profile_original=false)
    preparation=@timed fitted_train_graph(application)
    f=preparation.value
    owner=AlgorithmGraphs.prepared_graph_node(f.graph,:llowfs)
    pixels=AlgorithmGraphs.prepared_graph_node(f.graph,:camera_pixels)
    @test owner.train.workspace.execution!==nothing
    original()=begin
        TRAIN_EXTENSION._step_runtime_opd_reference!(owner)
        AlgorithmGraphs.enqueue_captured_graph_node!(pixels)
    end
    prepared()=begin
        AlgorithmGraphs.step_graph_node!(owner)
        AlgorithmGraphs.enqueue_captured_graph_node!(pixels)
    end
    start=time_ns();before=train_capture(f.graph,original);after=train_capture(f.graph,prepared)
    capture_seconds=(time_ns()-start)*1e-9
    println((;preparation_seconds=preparation.time,preparation_compile_seconds=preparation.compile_time,
        preparation_recompile_seconds=preparation.recompile_time,capture_seconds,before_nodes=train_node_count(before.graph),
        after_nodes=train_node_count(after.graph),storage=train_storage(owner)));flush(stdout)
    retained=(owner.train.workspace.wavefront.field,owner.train.workspace.execution)
    first_image=nothing
    for state in 1:3
        f.host.pupil_opd .= Float32[2e-8cos(state*(.07i-.11j)) for i in axes(f.host.pupil_opd,1),j in axes(f.host.pupil_opd,2)]
        f.host.dm_opd .= Float32[3e-8sin(state*(.15i+.09j)) for i in axes(f.host.dm_opd,1),j in axes(f.host.dm_opd,2)]
        f.host.pupil_amplitude .= f.reference .* Float32[.7+.2cos(state*(.03i+.05j))^2 for i in axes(f.reference,1),j in axes(f.reference,2)]
        foreach(pair->copyto!(pair...),zip(f.device,f.host));CUDA.synchronize()
        train_replay(f.graph,before)
        reference=Array(owner.train.output);reference_field=Array(owner.train.workspace.wavefront.field)
        reference_pixels=Array(graph_output(f.graph,:intensity))
        geometry=TRAIN_EXTENSION._wavefront_geometry(owner.train.workspace.wavefront)
        step_graph!(f.graph)
        actual=Array(owner.train.output);actual_pixels=Array(graph_output(f.graph,:intensity))
        relative=norm(actual-reference)/norm(reference);peak=maximum(abs,actual-reference)/maximum(abs,reference)
        pixel_relative=norm(actual_pixels-reference_pixels)/norm(reference_pixels)
        pixel_peak=maximum(abs,actual_pixels-reference_pixels)/maximum(abs,reference_pixels)
        field_relative=norm(Array(owner.train.workspace.wavefront.field)-reference_field)/norm(reference_field)
        println((;state,relative,peak,pixel_relative,pixel_peak,field_relative));flush(stdout)
        @test relative<=1e-5 && peak<=1e-5
        @test pixel_relative<=1e-5 && pixel_peak<=1e-5
        @test field_relative<=1e-5
        @test all(isfinite,actual_pixels) && sum(actual_pixels)>0
        @test TRAIN_EXTENSION._wavefront_geometry(owner.train.workspace.wavefront)==geometry
        @test owner.train.workspace.wavefront.field===retained[1] && owner.train.workspace.execution===retained[2]
        state==1 && (first_image=copy(actual_pixels))
    end
    f.host.pupil_opd .= Float32[2e-8cos(.07i-.11j) for i in axes(f.host.pupil_opd,1),j in axes(f.host.pupil_opd,2)]
    f.host.dm_opd .= Float32[3e-8sin(.15i+.09j) for i in axes(f.host.dm_opd,1),j in axes(f.host.dm_opd,2)]
    f.host.pupil_amplitude .= f.reference .* Float32[.7+.2cos(.03i+.05j)^2 for i in axes(f.reference,1),j in axes(f.reference,2)]
    foreach(pair->copyto!(pair...),zip(f.device,f.host));CUDA.synchronize()
    reset_graph!(f.graph);@test all(iszero,Array(graph_output(f.graph,:intensity)))
    step_graph!(f.graph);@test Array(graph_output(f.graph,:intensity))==first_image
    @test bytes2hex(sha256(read(f.path)))==f.identity
    if benchmark
        times=(Float64[],Float64[])
        for block in 1:8
            for index in (isodd(block) ? (1,2) : (2,1))
                captured=index==1 ? before : after
                for _ in 1:4;train_replay(f.graph,captured);end
                for _ in 1:8;push!(times[index],@elapsed train_replay(f.graph,captured));end
            end
        end
        for (index,t) in enumerate(times)
            println((;prepared=index==2,samples=length(t),median_ms=1e3median(t),p90_ms=1e3quantile(t,.9),maximum_ms=1e3maximum(t)));flush(stdout)
        end
    end
    if profile
        selected=profile_original ? before : after
        for _ in 1:4;train_replay(f.graph,selected);end
        CUDA.@profile external=true begin train_replay(f.graph,selected);end
    end
end
isempty(ARGS) && error("pass frozen fitted application directory")
@testset "Frozen LLOWFS prepared train" begin
    compare_fitted_train(abspath(first(ARGS));benchmark="--benchmark" in ARGS,
        profile="--profile" in ARGS,profile_original="--profile-original" in ARGS)
end
