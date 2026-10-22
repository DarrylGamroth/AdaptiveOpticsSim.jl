# Fixed numerical sampling; preparation-only batch override, captured AOS rate node.
# julia --project=CUDA_ENV sweep_modulation_batches_cuda.jl OUTPUT.toml [CASE,...]
# AOS_BATCH_COMPARE_POLICY=1 measures the production choice as well as baseline.
include("modulation_batch_override.jl")
using AdaptiveOpticsSim.Backends, AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.WavefrontSensors
using LinearAlgebra, Statistics, Random, TOML, SHA
import CUDA

const BATCH_CASES = (
    (name="small", n=64, pupils=16, points=8, separation=2, edge=1),
    (name="medium", n=128, pupils=32, points=32, separation=2, edge=1),
    (name="wide", n=256, pupils=64, points=64, separation=2, edge=1),
    (name="gpi", n=240, pupils=60, points=24, separation=60, edge=30),
    (name="copper", n=480, pupils=30, points=32, separation=2, edge=1),
    (name="many", n=240, pupils=30, points=96, separation=2, edge=1),
    (name="large", n=384, pupils=48, points=128, separation=2, edge=1),
    (name="fine", n=720, pupils=60, points=192, separation=60, edge=30),
    (name="pow256", n=120, pupils=30, points=32, separation=2, edge=1),
    (name="pow512", n=240, pupils=30, points=32, separation=2, edge=1),
    (name="pow512short", n=240, pupils=30, points=12, separation=2, edge=1),
    (name="pow768", n=360, pupils=30, points=32, separation=2, edge=1),
    (name="nonpow", n=192, pupils=32, points=32, separation=2, edge=1),
    (name="odd24", n=128, pupils=32, points=24, separation=2, edge=1),
    (name="odd30", n=128, pupils=32, points=30, separation=2, edge=1),
    (name="prime31", n=128, pupils=32, points=31, separation=2, edge=1),
    (name="edge560", n=272, pupils=68, points=32, separation=2, edge=1),
    (name="edge576", n=280, pupils=70, points=32, separation=2, edge=1),
    (name="pow512six", n=240, pupils=30, points=6, separation=2, edge=1),
    (name="pow512eight", n=240, pupils=30, points=8, separation=2, edge=1),
    (name="pow512sixteen", n=240, pupils=30, points=16, separation=2, edge=1),
    (name="pow512twentyfour", n=240, pupils=30, points=24, separation=2, edge=1),
    (name="pow512sixtyfour", n=240, pupils=30, points=64, separation=2, edge=1),
)

function batch_graph(c, strategy; T=Float32)
    opd=CUDA.zeros(T,c.n,c.n)
    host=T[40e-9*(sinpi(2i/c.n)+0.4cospi(2j/c.n)+0.2sinpi((i+j)/c.n)) for i in 1:c.n,j in 1:c.n]
    copyto!(opd,host)
    definition=algorithm_graph((pyramid_rate_node(:pwfs; resolution=c.n,
        telescope_diameter_m=8.,pupil_samples=c.pupils,modulation=3.,
        modulation_points=c.points,modulation_propagation_strategy=strategy,
        n_pix_separation=c.separation,n_pix_edge=c.edge,
        source_wavelength_m=0.75e-6,source_photon_irradiance_m2_s=2e8,
        opd_schema="benchmark.batch.opd/1",photon_rate_schema="benchmark.batch.rate/1"),);
        inputs=(graph_input(:opd,:pwfs=>:opd,opd),),
        outputs=(graph_output(:rate,:pwfs=>:photon_rate),),name=:modulation_batch_sweep)
    prepare_algorithm_graph(definition;target=compute_device(opd),execution=CapturedGraphExecution())
end

function measure_batch(c,strategy,batch,reference; warmups=16,samples=48)
    MODULATION_BATCH[]=batch
    graph=batch_graph(c,strategy)
    workspace=AlgorithmGraphs.prepared_graph_node(graph,:pwfs).prepared.workspace.modulation_batch
    actual=workspace.batch_size
    for _ in 1:warmups;step_graph!(graph);end
    output=Array(graph_output(graph,:rate))
    relative=reference===nothing ? 0. : norm(Float64.(output).-reference)/norm(reference)
    reference===nothing || @assert relative < 1e-5
    bytes=[@allocated(step_graph!(graph)) for _ in 1:4]
    @assert all(iszero,bytes)
    times=zeros(Float64,samples)
    for i in eachindex(times)
        t=time_ns();step_graph!(graph);times[i]=(time_ns()-t)/1e6
    end
    record=Dict("batch"=>actual,"samples_ms"=>times,"median_ms"=>median(times),
        "p95_ms"=>quantile(times,.95),"maximum_ms"=>maximum(times),
        "host_bytes"=>bytes,"relative_l2"=>relative,"flux"=>sum(Float64,output),
        "field_stack_bytes"=>sizeof(eltype(workspace.field_stack))*length(workspace.field_stack),
        "padded_side"=>size(workspace.field_stack,1),"algorithm_nodes"=>captured_graph_node_count(graph))
    println("$(c.name) $(nameof(typeof(strategy))) batch=$actual median_ms=$(record["median_ms"]) relative=$relative");flush(stdout)
    return record,output
end

function sweep_main(args)
    length(args) in (1,2) || error("OUTPUT.toml [CASE,...]")
    CUDA.allowscalar(false);BLAS.set_num_threads(1);Backends.set_fft_provider_threads!(1)
    selected=length(args)==2 ? split(args[2],',') : [c.name for c in BATCH_CASES]
    strategies=get(ENV,"AOS_BATCH_STRATEGIES","pupil")
    strategies == "pupil" || error("only physical pupil-tilt propagation is supported")
    rounds=parse(Int,get(ENV,"AOS_BATCH_ROUNDS","1"))
    samples=parse(Int,get(ENV,"AOS_BATCH_SAMPLES","48"))
    output=abspath(args[1]);mkpath(dirname(output))
    result=Dict{String,Any}("revision"=>readchomp(`git -C $(@__DIR__) rev-parse HEAD`),
        "status"=>readchomp(`git -C $(@__DIR__) status --porcelain`),
        "julia"=>string(VERSION),"cuda"=>string(pkgversion(CUDA)),
        "gpu"=>CUDA.name(CUDA.device()),"runtime"=>sprint(io->CUDA.versioninfo(io)),
        "l2_bytes"=>CUDA.attribute(CUDA.device(),CUDA.DEVICE_ATTRIBUTE_L2_CACHE_SIZE),
        "harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "override_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"modulation_batch_override.jl")))),
        "setup_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","wfs","pyramid","setup.jl")))),
        "cuda_extension_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"..","ext","AdaptiveOpticsSimCUDAExt.jl")))),
        "boundary"=>"Captured AOS photon-rate node through completed step; no detector/RTC/host copies; closed-loop unpaced service samples",
        "warmups"=>16,"samples"=>samples,"rounds"=>rounds,"cases"=>Any[])
    persist()=open(io->TOML.print(io,result;sorted=true),output,"w")
    persist()
    for c in BATCH_CASES
        c.name in selected || continue
        for (name,strategy) in (("pupil",WavefrontSensors.PyramidPupilTiltStrategy()),)
            name in split(strategies,',') || continue
            base,reference=measure_batch(c,strategy,0,nothing;samples)
            pad=base["padded_side"]
            bound=max(1,div(256*1024^2,pad^2*sizeof(ComplexF32)))
            candidates=filter(b->c.points%b==0 && b<=bound,1:64)
            case=Dict{String,Any}("name"=>c.name,"configuration"=>Dict(string(k)=>v for (k,v) in pairs(c)),"strategy"=>name,
                "baseline"=>base,"trials"=>Any[])
            push!(result["cases"],case);persist()
            if get(ENV,"AOS_BATCH_COMPARE_POLICY","0")=="1"
                selected_record,_=measure_batch(c,strategy,-1,reference;samples)
                case["selected"]=selected_record;persist()
            end
            for round in 1:rounds
                policy_only=get(ENV,"AOS_BATCH_POLICY_ONLY","0")=="1"
                order=shuffle(MersenneTwister(812+round),policy_only ? [0,-1] : collect(candidates))
                for batch in order
                    GC.gc(true);CUDA.reclaim()
                    record,_=measure_batch(c,strategy,batch,reference;samples)
                    record["round"]=round
                    policy_only && (record["policy"]=batch==0 ? "baseline" : "selected")
                    push!(case["trials"],record);persist()
                end
            end
        end
    end
    result["complete"]=true;persist()
end
if abspath(PROGRAM_FILE) == @__FILE__
    sweep_main(ARGS)
end
