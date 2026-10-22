# Fixed physical sampling: captured completed photon-rate nodes, no detector/RTC.
include("fft_layout_override.jl")
include("fft_pair_layout_override.jl")
include("profile_physical_pyramid_cuda.jl")
_benchmark_fft_storage_bytes(plan::FFTLayoutExperiment.TransposeFFTPlan) =
    FFTLayoutExperiment.scratch_bytes(plan)
_benchmark_fft_storage_bytes(plan::FFTPairLayoutExperiment.PairFFTPlan) =
    FFTLayoutExperiment.scratch_bytes(plan)

# Deliberately enlarge only the experimental stack budget; production stays fixed.
@eval WavefrontSensors function _pyramid_modulation_batch_size(points::Int, bytes::Int,
    limit::Int=_PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT,
    budget_bytes::Int=_PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES)
    points >= 1 && bytes >= 1 || throw(InvalidConfiguration("invalid benchmark dimensions"))
    requested = Main.MODULATION_BATCH[]
    if requested > 0
        requested <= max(1, div(1536*1024^2, bytes)) && points % requested == 0 ||
            throw(InvalidConfiguration("experimental batch must divide points and fit 1.5 GiB stack budget"))
        return requested
    end
    selected_budget = requested == 0 ? _PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES : budget_bytes
    selected_limit = requested == 0 ? _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT : limit
    bound = max(1, div(selected_budget, bytes))
    for candidate in min(points, selected_limit, bound):-1:1
        points % candidate == 0 && return candidate
    end
    return 1
end

function layout_sweep(args)
    length(args)==2 || error("OUTPUT.toml CASE")
    output=abspath(args[1]); c=only(filter(c->c.name==args[2],PHYSICAL_PROFILE_CASES))
    CUDA.allowscalar(false);BLAS.set_num_threads(1);Backends.set_fft_provider_threads!(1)
    warmups=parse(Int,get(ENV,"AOS_FFT_WARMUPS","16"))
    samples=parse(Int,get(ENV,"AOS_FFT_SAMPLES","24"))
    rounds=parse(Int,get(ENV,"AOS_FFT_ROUNDS","2"))
    batches=parse.(Int,split(get(ENV,"AOS_FFT_BATCHES",c.name=="gpi_fine" ? "4,8,16" : "2,4,8,16"),','))
    layout=get(ENV,"AOS_FFT_LAYOUT","native")
    result=Dict{String,Any}("case"=>c.name,"layout"=>layout,"pair_layout"=>FFTPairLayoutExperiment.ENABLED,"transpose_kind"=>FFTTransposeExperiment.TRANSPOSE,"julia"=>string(VERSION),
        "source_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f))))
            for f in ("sweep_fft_layout_cuda.jl","fft_layout_override.jl",
                "fft_pair_layout_override.jl","fft_transpose_cuda.jl","profile_physical_pyramid_cuda.jl",
                "sweep_modulation_batches_cuda.jl","modulation_batch_override.jl")),
        "revision"=>readchomp(`git -C $(@__DIR__) rev-parse HEAD`),
        "dirty_status"=>readchomp(`git -C $(@__DIR__) status --porcelain`),
        "runtime"=>sprint(io->CUDA.versioninfo(io)),"warmups"=>warmups,"samples"=>samples,
        "rounds"=>rounds,"stack_budget_bytes"=>1536*1024^2,
        "boundary"=>"Unpaced completed captured photon-rate node; fixed sampling, single outstanding call; setup and transfers excluded",
        "trials"=>Any[])
    persist()=open(io->TOML.print(io,result;sorted=true),output,"w")
    mkpath(dirname(output));persist()
    reference=nothing
    for round in 1:rounds, batch in shuffle(MersenneTwister(1967+round),batches)
        GC.gc(true);CUDA.reclaim()
        trial,rate=measure_batch(c,WavefrontSensors.PyramidPupilTiltStrategy(),batch,reference;warmups,samples)
        reference===nothing && (reference=Float64.(rate))
        trial["round"]=round;trial["layout"]=layout
        path=output*".r$(round).b$(batch).f32";open(io->write(io,vec(rate)),path,"w")
        trial["rate_file"]=basename(path);trial["rate_sha256"]=bytes2hex(sha256(read(path)))
        push!(result["trials"],trial);persist()
    end
    result["complete"]=true;persist()
end
if abspath(PROGRAM_FILE) == @__FILE__
    layout_sweep(ARGS)
end
