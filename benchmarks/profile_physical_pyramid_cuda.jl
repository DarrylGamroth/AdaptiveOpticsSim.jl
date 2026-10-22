# Completed-node baseline and a warmed CUDA-profiler range; fixed physical sampling.
# Usage: OUTPUT.toml CASE [BATCH], where BATCH=-1 selects production policy.
include("sweep_modulation_batches_cuda.jl")

const PHYSICAL_PROFILE_CASES = (
    (name="gpi_default", n=240, pupils=60, points=24, separation=60, edge=30),
    (name="gpi_fine", n=720, pupils=60, points=192, separation=60, edge=30),
    (name="copper_default", n=480, pupils=30, points=32, separation=2, edge=1),
    (name="copper_fast", n=240, pupils=30, points=96, separation=2, edge=1),
    (name="copper_fine", n=1440, pupils=30, points=256, separation=2, edge=1),
)

function physical_profile_main(args)
    length(args) in (2,3) || error("OUTPUT.toml CASE [BATCH]")
    matches = filter(c -> c.name == args[2], PHYSICAL_PROFILE_CASES)
    length(matches) == 1 || error("unknown physical profile case: $(args[2])")
    c = only(matches)
    MODULATION_BATCH[] = length(args)==3 ? parse(Int,args[3]) : -1
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    Threads.nthreads() == 1 || error("use one Julia thread")
    graph = batch_graph(c,WavefrontSensors.PyramidPupilTiltStrategy())
    prepared = AlgorithmGraphs.prepared_graph_node(graph,:pwfs).prepared
    workspace = prepared.workspace.modulation_batch
    phases = prepared.front_end.modulation.phases
    warmups = parse(Int,get(ENV,"AOS_PROFILE_WARMUPS","16"))
    samples = parse(Int,get(ENV,"AOS_PROFILE_SAMPLES","32"))
    repetitions = parse(Int,get(ENV,"AOS_PROFILE_REPLAYS","4"))
    for _ in 1:warmups; step_graph!(graph); end
    reference = Array(graph_output(graph,:rate))
    all(isfinite,reference) && sum(reference)>0 || error("invalid rate map")
    allocations = [@allocated(step_graph!(graph)) for _ in 1:16]
    all(iszero,allocations) || error("steady host allocation")
    memory_before = CUDA.used_memory()
    times = zeros(Float64,samples)
    for i in eachindex(times)
        start = time_ns(); step_graph!(graph); times[i]=(time_ns()-start)/1e6
    end
    CUDA.used_memory() == memory_before || error("device memory changed")
    if get(ENV,"AOS_PROFILE_EXTERNAL","0") == "1"
        CUDA.@profile external=true begin
            for _ in 1:repetitions; step_graph!(graph); end
        end
    end
    after = Array(graph_output(graph,:rate))
    after == reference || error("deterministic replay changed rate map")
    record=Dict("case"=>Dict(string(k)=>v for (k,v) in pairs(c)),
        "revision"=>readchomp(`git -C $(@__DIR__) rev-parse HEAD`),
        "harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "julia"=>string(VERSION),"runtime"=>sprint(io->CUDA.versioninfo(io)),
        "gpu"=>CUDA.name(CUDA.device()),"batch"=>workspace.batch_size,
        "padded_side"=>size(workspace.field_stack,1),
        "field_stack_bytes"=>sizeof(eltype(workspace.field_stack))*length(workspace.field_stack),
        "phase_bytes"=>sizeof(eltype(phases))*length(phases),
        "device_used_bytes"=>memory_before,"host_allocated_bytes"=>allocations,
        "warmups"=>warmups,"samples_ms"=>times,"median_ms"=>median(times),
        "p95_ms"=>quantile(times,.95),"maximum_ms"=>maximum(times),
        "minimum_ms"=>minimum(times),"profile_replays"=>repetitions,
        "output_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(after)))),
        "output_sum"=>sum(Float64,after),"output_shape"=>collect(size(after)),
        "profile_enabled"=>get(ENV,"AOS_PROFILE_EXTERNAL","0")=="1",
        "boundary"=>"One captured physical AOS photon-rate node through completion; excludes detector, RTC, transfers, setup; unpaced single outstanding service")
    mkpath(dirname(abspath(args[1])))
    open(io->TOML.print(io,record;sorted=true),args[1],"w")
    println("$(c.name): pad=$(record["padded_side"]) batch=$(record["batch"]) median_ms=$(record["median_ms"]) phase_bytes=$(record["phase_bytes"])")
end

if abspath(PROGRAM_FILE) == @__FILE__
    physical_profile_main(ARGS)
end
