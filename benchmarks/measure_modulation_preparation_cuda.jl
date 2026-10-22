# Fresh-process preparation measurement; environment is prepared by the caller.
# julia --project=ENV --startup-file=no HARNESS COPPER_ROOT OUTPUT_TOML N M
# Warm a disposable n60/m4 captured calibration owner, release it, then measure
# one target preparation while retaining its owner through all measurements.
using AdaptiveOpticsSim, AdaptiveOpticsSim.AlgorithmGraphs, REVOLTCopperSim
using LinearAlgebra, TOML, SHA
import CUDA

# These files have guarded command-line entry points; including them here only
# defines helpers. The Copper source is selected explicitly by the caller.
include(joinpath(@__DIR__,"validate_modulation_storage_cuda.jl"))

function preparation_main(args)
    length(args)==4 || error("usage: HARNESS COPPER_ROOT OUTPUT_TOML N M")
    root=realpath(args[1]); output=abspath(args[2])
    n=parse(Int,args[3]); points=parse(Int,args[4])
    n>0 && points>0 || error("N and M must be positive")
    isfile(output) && error("OUTPUT_TOML already exists")
    realpath(dirname(dirname(pathof(REVOLTCopperSim))))==root || error("loaded Copper source differs from COPPER_ROOT")
    CUDA.functional() || error("functional CUDA required")
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    AOS.Backends.set_fft_provider_threads!(1)
    target=AOS.Backends.compute_device(CUDA.zeros(Float32,1))
    warm=Main.CopperConvergence.prepare_case(60,4;target)
    step_graph!(warm.graph)
    CUDA.synchronize()
    warm=nothing
    GC.gc(true)
    CUDA.reclaim()
    CUDA.synchronize()
    pool_before=pool_record()
    rss_before=Sys.maxrss()
    # @timed measures preparation's host allocations and wall/GC time. The
    # final sync includes all device preparation in the measured interval.
    system=nothing
    measurement=@timed begin
        system=Main.CopperConvergence.prepare_case(n,points;target)
        CUDA.synchronize()
        nothing
    end
    pool_after=pool_record()
    rss_after=Sys.maxrss()
    modulation=modulation_record(system.graph,n,points)
    direct_host=sum(record["bytes"] for record in values(modulation["arrays"]) if record["backend"]=="host";init=0)
    direct_device=sum(record["bytes"] for record in values(modulation["arrays"]) if record["backend"]=="CUDA";init=0)
    project=something(Base.active_project(),"")
    helper=joinpath(root,"benchmarks","copper_convergence.jl")
    result=Dict{String,Any}(
        "complete"=>true,"resolution"=>n,"modulation_points"=>points,
        "host_allocated_bytes"=>measurement.bytes,
        "preparation_seconds"=>measurement.time,"gc_seconds"=>measurement.gctime,
        "maxrss_before_bytes"=>rss_before,"maxrss_after_bytes"=>rss_after,
        "maxrss_highwater_increase_bytes"=>rss_after-rss_before,
        "maxrss_scope"=>"process lifetime high-water resident bytes, including startup, compilation and warm graph; difference is not isolated target peak; compare fresh processes",
        "pool_before"=>pool_before,"pool_after"=>pool_after,
        "modulation"=>modulation,"direct_modulation_host_bytes"=>direct_host,
        "direct_modulation_device_bytes"=>direct_device,
        "warmup"=>"one disposable captured calibration graph n60/m4 and one completed graph frame, then GC/reclaim",
        "measurement_boundary"=>"one CopperConvergence.prepare_case including captured graph and HIL boundary preparation plus final CUDA synchronization; target owner remains live through memory inspection",
        "aos"=>source_record(AOS),"copper"=>source_record(REVOLTCopperSim),
        "julia"=>string(VERSION),"cuda"=>string(pkgversion(CUDA)),
        "device"=>CUDA.name(CUDA.device()),"julia_threads"=>Threads.nthreads(),
        "active_project"=>project,"arguments"=>copy(args),"harness"=>abspath(@__FILE__),
        "harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "validation_harness_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"validate_modulation_storage_cuda.jl")))),
        "copper_helper_sha256"=>bytes2hex(sha256(read(helper))),
        "command"=>join(Base.julia_cmd().exec," ")*" --project="*project*" --startup-file=no "*join(vcat([abspath(@__FILE__)],args)," "))
    for path in (project,joinpath(dirname(project),"Manifest.toml"))
        isfile(path) && (result[basename(path)*"_sha256"]=bytes2hex(sha256(read(path))))
    end
    # Preserve liveness until after the pool and direct storage observations.
    GC.@preserve system begin
        mkpath(dirname(output))
        open(io->TOML.print(io,result;sorted=true),output*".tmp","w")
        mv(output*".tmp",output;force=true)
    end
    println("Wrote ",output)
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==4 || error("usage: HARNESS COPPER_ROOT OUTPUT_TOML N M")
    include(joinpath(realpath(ARGS[1]),"benchmarks","copper_convergence.jl"))
    preparation_main(ARGS)
end
