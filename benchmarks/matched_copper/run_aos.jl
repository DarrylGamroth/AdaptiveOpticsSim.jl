# julia --project=ENV run_aos.jl FIXTURE OUTPUT cpu|cuda FFT_THREADS preflight|validate|benchmark
using AdaptiveOpticsSim, AdaptiveOpticsSim.Backends, AdaptiveOpticsSim.AlgorithmGraphs
using LinearAlgebra, Statistics, TOML, SHA
import CUDA
const MatchedWFS=AdaptiveOpticsSim.WavefrontSensors
function read_array(root,meta)
    data=read(joinpath(root,meta["file"]))
    bytes2hex(sha256(data))==meta["sha256"] || error("fixture digest mismatch")
    return reshape(copy(reinterpret(Float32,data)),Tuple(meta["shape"]))
end
function store_array(root,name,values)
    path=joinpath(root,name*".f32")
    host=Float32.(Array(values));open(io->write(io,vec(host)),path,"w")
    Dict("file"=>basename(path),"shape"=>collect(size(host)),"dtype"=>"float32",
        "order"=>"F","sha256"=>bytes2hex(sha256(read(path))))
end
function persist_record(root,record)
    p=joinpath(root,"record.toml");open(io->TOML.print(io,record;sorted=true),p*".tmp","w");mv(p*".tmp",p;force=true)
end
function make_graph(c,gpu)
    opd=gpu ? CUDA.zeros(Float32,240,240) : zeros(Float32,240,240)
    definition=algorithm_graph((pyramid_rate_node(:pwfs;resolution=240,telescope_diameter_m=c["diameter_m"],
        central_obstruction_ratio=0.,pupil_reflectivity=1.,T=Float32,
        pupil_samples=30,modulation=5.,modulation_points=96,
        modulation_propagation_strategy=MatchedWFS.PyramidPupilTiltStrategy(),
        diffraction_padding=2,psf_centering=true,n_pix_separation=2,n_pix_edge=1,binning=1,
        source_wavelength_m=c["wavelength_m"],source_photon_irradiance_m2_s=c["photon_irradiance"],
        opd_schema="benchmark.matched-copper.opd-m/1",photon_rate_schema="benchmark.matched-copper.rate/1"),);
        inputs=(graph_input(:opd,:pwfs=>:opd,opd),),outputs=(graph_output(:rate,:pwfs=>:photon_rate),),name=:matched_copper)
    return prepare_algorithm_graph(definition;target=compute_device(opd),
        execution=gpu ? CapturedGraphExecution() : StreamGraphExecution())
end
function run_main(args)
    length(args)==6 || error("FIXTURE OUTPUT cpu|cuda FFT_THREADS preflight|validate|benchmark ROUND")
    fixture,output,backend,threadstr,mode,roundstr=args
    backend in ("cpu","cuda") ||error("backend")
    mode in ("preflight","validate","benchmark") ||error("mode")
    isdir(output) && !isempty(readdir(output)) &&error("output must be empty")
    mkpath(output)
    nt=parse(Int,threadstr);round=parse(Int,roundstr);gpu=backend=="cuda"
    Threads.nthreads()==1 ||error("use one Julia thread")
    BLAS.set_num_threads(1);Backends.set_fft_provider_threads!(nt);CUDA.allowscalar(false)
    meta=TOML.parsefile(joinpath(fixture,"fixture.toml"));c=meta["contract"]
    probes=read_array(fixture,meta["probes"]);amplitude=read_array(fixture,meta["amplitude"])
    protocol=TOML.parsefile(joinpath(@__DIR__,"protocol.toml"))
    prepared_start=time_ns();graph=make_graph(c,gpu);preparation_s=(time_ns()-prepared_start)/1e9
    owner=AlgorithmGraphs.prepared_graph_node(graph,:pwfs).prepared
    Array(owner.input.amplitude)==amplitude ||error("shared pupil differs")
    pad=size(owner.workspace.field,1)
    pad==512 ||error("FFT side mismatch")
    rate=graph_output(graph,:rate);size(rate)==(64,64)||error("rate shape")
    eltype(rate)===Float32 ||error("rate precision")
    input=graph_input(graph,:opd)
    staging=zeros(Float32,240,240)
    function stage!(host)
        copyto!(staging,host)
        copyto!(input,staging)
        gpu &&CUDA.synchronize(;blocking=true)
        return nothing
    end
    rates=Array{Float32}(undef,64,64,size(probes,3))
    for i in axes(probes,3)
        stage!(view(probes,:,:,i));step_graph!(graph);rates[:,:,i].=Array(rate)
    end
    all(isfinite,rates) &&minimum(rates)>=0 ||error("invalid rates")
    rates[:,:,1]!=rates[:,:,8] ||error("stale captured input")
    batch=owner.workspace.modulation_batch
    actual_batch=hasproperty(batch,:batch_size) ? batch.batch_size : 1
    record=Dict{String,Any}("engine"=>"aos","backend"=>backend,"cpu_workers"=>nt,
        "mode"=>mode,"round"=>round,"complete"=>false,"prescription"=>"common",
        "revision"=>readchomp(`git -C $(dirname(dirname(pathof(AdaptiveOpticsSim)))) rev-parse HEAD`),
        "source_sha256"=>Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in (@__FILE__,joinpath(@__DIR__,"protocol.toml"))),
        "fixture_sha256"=>bytes2hex(sha256(read(joinpath(fixture,"fixture.toml")))),
        "julia"=>string(VERSION),"fft_side"=>pad,"output_dtype"=>string(eltype(rate)),
        "fft_dtype"=>string(eltype(owner.workspace.field)),"actual_batch"=>actual_batch,
        "cuda_graph"=>gpu,"preparation_seconds"=>preparation_s,
        "cpu"=>Sys.cpu_info()[1].model,"rates"=>store_array(output,"rates",rates),
        "boundary"=>protocol["scope"],"arrival"=>protocol["arrival"])
    if gpu
        record["cuda_runtime"]=sprint(io->CUDA.versioninfo(io))
        record["gpu"]=CUDA.name(CUDA.device())
    end
    persist_record(output,record)
    if mode=="validate"
        basis=read_array(fixture,meta["calibration_basis"])
        positive=Array{Float32}(undef,64,64,277);negative=similar(positive)
        scratch=zeros(Float32,240,240);poke=Float32(c["calibration_poke_m"])
        for k in 1:277
            @views @. scratch=poke*basis[:,:,k]
            stage!(scratch);step_graph!(graph);positive[:,:,k].=Array(rate)
            @views @. scratch=-poke*basis[:,:,k]
            stage!(scratch);step_graph!(graph);negative[:,:,k].=Array(rate)
            if k%32==0;println("calibration $k/277");flush(stdout);end
        end
        record["calibration_positive"]=store_array(output,"calibration_positive",positive)
        record["calibration_negative"]=store_array(output,"calibration_negative",negative)
        persist_record(output,record)
    elseif mode=="benchmark"
        p=protocol["performance"];warmups=p["warmups"];samples=p["samples"]
        for k in 1:warmups
            i=mod1(k,size(probes,3));stage!(view(probes,:,:,i));step_graph!(graph)
        end
        stage!(view(probes,:,:,13));step_graph!(graph)
        record["host_allocation_bytes"]=[@allocated(step_graph!(graph)) for _ in 1:8]
        gpu &&!all(iszero,record["host_allocation_bytes"]) &&error("AOS CUDA host allocation")
        memory_before=gpu ? CUDA.used_memory() : 0
        times=zeros(Float64,samples)
        for k in 1:samples
            i=mod1(k,size(probes,3));stage!(view(probes,:,:,i))
            start=time_ns();step_graph!(graph);times[k]=(time_ns()-start)/1e6
        end
        Array(rate)==rates[:,:,mod1(samples,size(probes,3))] ||error("replay changed rates")
        gpu &&CUDA.used_memory()!=memory_before &&error("CUDA pool changed")
        record["service"]=Dict("samples_ms"=>times,"median_ms"=>median(times),"minimum_ms"=>minimum(times),
            "maximum_ms"=>maximum(times),"p95_ms"=>quantile(times,.95),"warmups"=>warmups,"samples"=>samples)
        record["observed_cuda_pool_bytes"]=memory_before
    end
    record["complete"]=true;persist_record(output,record)
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_main(ARGS)
end
