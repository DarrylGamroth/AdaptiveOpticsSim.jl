# Run the identical file against each checkout in separate, prepared environments:
# julia --project=ISOLATED_ENV --startup-file=no HARNESS COPPER_ROOT OUTPUT_DIR [480:32,240:96]
# Optional AOS_MODULATION_MEMORY_ONLY=1440:256 prepares that extra case without stepping.
# Dependencies: AdaptiveOpticsSim, REVOLTCopperSim, CUDA; stdlibs only otherwise.
# This script never activates, installs, develops, or changes an environment.
using AdaptiveOpticsSim, AdaptiveOpticsSim.AlgorithmGraphs, REVOLTCopperSim
using LinearAlgebra, Statistics, TOML, SHA
import CUDA

const AOS = AdaptiveOpticsSim
const AG = AOS.AlgorithmGraphs

function source_record(mod)
    root = dirname(dirname(pathof(mod)))
    Dict("root"=>root, "version"=>string(pkgversion(mod)),
        "revision"=>readchomp(`git -C $root rev-parse HEAD`),
        "status"=>readchomp(`git -C $root status --porcelain`),
        "diff_sha256"=>bytes2hex(sha256(read(`git -C $root diff HEAD`))))
end
function write_array(directory, name, values::AbstractArray{T}) where {T<:Union{Float32,Float64,ComplexF32,ComplexF64}}
    ENDIAN_BOM == 0x04030201 || error("little-endian host required")
    extension = T == Float32 ? ".f32" : T == Float64 ? ".f64" : T == ComplexF32 ? ".c64" : ".c128"
    path = joinpath(directory, name * extension)
    open(io->write(io, vec(Array(values))), path, "w")
    Dict("file"=>basename(path), "shape"=>collect(size(values)), "eltype"=>string(T),
        "order"=>"Julia column-major little-endian", "bytes"=>filesize(path),
        "sha256"=>bytes2hex(sha256(read(path))))
end
function persist(directory, result)
    path = joinpath(directory, "validation.toml")
    open(io->TOML.print(io, result; sorted=true), path*".tmp", "w")
    mv(path*".tmp", path; force=true)
end
function pool_record()
    CUDA.synchronize()
    # CUDA.used_memory is documented as application-live pool memory. It returns
    # missing on unsupported allocators, which must not be reported as zero.
    used = isdefined(CUDA, :used_memory) ? CUDA.used_memory() : missing
    cached = isdefined(CUDA, :cached_memory) ? CUDA.cached_memory() : missing
    Dict("used_bytes"=>ismissing(used) ? "unavailable" : used,
        "reserved_bytes"=>ismissing(cached) ? "unavailable" : cached)
end
function modulation_record(graph, n, points)
    # Diagnostic backend-owner access stays in this benchmark. Account each
    # numeric storage array once, without traversing CUDA handles or FFT plans.
    prepared = AG.prepared_graph_node(graph, :pwfs).prepared
    modulation = prepared.front_end.modulation
    arrays = Dict{String,Any}()
    seen = IdDict()
    for name in fieldnames(typeof(modulation))
        value = getfield(modulation, name)
        value isa AbstractArray || continue
        haskey(seen, value) && continue
        seen[value] = true
        arrays[string(name)] = Dict("bytes"=>sizeof(eltype(value))*length(value),
            "shape"=>collect(size(value)), "eltype"=>string(eltype(value)),
            "backend"=>value isa CUDA.CuArray ? "CUDA" : "host")
    end
    Dict("owner_type"=>string(typeof(modulation)), "arrays"=>arrays,
        "retained_array_bytes"=>sum(record["bytes"] for record in values(arrays)),
        "baseline_dense_phase_bytes"=>n*n*points*sizeof(ComplexF32),
        "scope"=>"direct retained modulation arrays; excludes policy, shared propagation, masks, FFT plans, capture and allocator reserve")
end
function prepare_case(n, points, target; calibration, execution)
    definition = REVOLTCopperSim._graph_definition(:grid_gaussian, target;
        calibration, pupil_resolution=n, modulation_points=points)
    prepare_algorithm_graph(definition; target, execution)
end
function timings!(graph; warmups=32, samples=128, require_zero=true)
    for _ in 1:warmups
        step_graph!(graph)
    end
    bytes = [@allocated(step_graph!(graph)) for _ in 1:16]
    require_zero && !all(iszero, bytes) && error("graph step allocated host bytes: $bytes")
    times = Vector{Float64}(undef, samples)
    for i in eachindex(times)
        start = time_ns()
        step_graph!(graph) # prepared graph completion is included
        times[i] = (time_ns()-start)/1e6
    end
    Dict("samples_ms"=>times, "minimum_ms"=>minimum(times),
        "median_ms"=>median(times), "p95_ms"=>quantile(times,0.95),
        "maximum_ms"=>maximum(times), "host_bytes"=>bytes,
        "warmups"=>warmups, "samples"=>samples, "zero_host_allocation_required"=>require_zero)
end
function replay!(graph; frames=32)
    reset_graph!(graph)
    fill!(graph_input(graph,:pdm_command),0f0)
    CUDA.synchronize()
    rates = Array{Float32}(undef,64,64,frames)
    images = similar(rates)
    for index in 1:frames
        step_graph!(graph)
        rates[:,:,index] .= Array(graph_output(graph,:pwfs_photon_rate))
        images[:,:,index] .= Array(graph_output(graph,:pwfs_frame))
    end
    (; rates, images)
end
function compare(a,b)
    delta = Float64.(a).-Float64.(b)
    scale = norm(Float64.(b))
    Dict("identical"=>a==b, "max_abs"=>maximum(abs,delta),
        "relative_l2"=>iszero(scale) ? norm(delta) : norm(delta)/scale)
end
function interaction64!(system, divisor, convergence)
    command = graph_input(system.graph,:pdm_command)
    host = zeros(Float32,277)
    matrix = Matrix{Float64}(undef,2length(convergence.OFFSETS),277)
    for actuator in 1:277
        host[actuator]=Float32(convergence.POKE)
        copyto!(command,host); CUDA.synchronize(); step_graph!(system.graph)
        positive = convergence.raw_slopes(Array(graph_output(system.graph,:pwfs_photon_rate)),divisor)
        host[actuator]=-Float32(convergence.POKE)
        copyto!(command,host); CUDA.synchronize(); step_graph!(system.graph)
        negative = convergence.raw_slopes(Array(graph_output(system.graph,:pwfs_photon_rate)),divisor)
        matrix[:,actuator] .= (positive.-negative)./(2convergence.POKE)
        host[actuator]=0f0
        actuator % 32 == 0 && println("interaction ",actuator,"/277")
    end
    fill!(command,0f0); CUDA.synchronize()
    all(isfinite,matrix) || error("nonfinite interaction")
    all(column->norm(column)>0,eachcol(matrix)) || error("zero interaction column")
    matrix
end
function dense_preparation!(directory)
    result = Dict{String,Any}()
    O = AOS.Optics
    for T in (Float32,Float64), n in (7,8)
        policies = (O.NoModulation(), O.CircularModulation(T(2.3);samples=9,phase_offset_rad=T(0.37),T),
            O.SampledModulation([(0.2,-0.7),(1.1,0.4),(-0.3,0.8)];weights=[1,3,0],T))
        for (index,policy) in enumerate(policies)
            host = O.prepare_focal_plane_modulation(policy,n,zeros(T,1),T)
            gpu = O.prepare_focal_plane_modulation(policy,n,CUDA.zeros(T,1),T)
            CUDA.synchronize()
            phases = Array(gpu.phases)
            error = maximum(abs, phases.-host.phases)
            isapprox(phases, host.phases; rtol=8eps(T), atol=8eps(T)) || Base.error("dense modulation CPU/CUDA mismatch: $T n=$n policy=$index error=$error")
            key = "dense-$(T)-n$(n)-policy$(index)"
            result[key] = Dict("cpu"=>write_array(directory,key*"-cpu",host.phases),
                "cuda"=>write_array(directory,key*"-cuda",phases), "max_abs"=>Float64(error),
                "rtol"=>Float64(8eps(T)), "atol"=>Float64(8eps(T)), "amplitude_weights"=>Float64.(host.amplitude_weights))
        end
    end
    result
end
function run_case!(directory,n,points,target,convergence,progress!)
    prefix="n$(n)-m$(points)"
    GC.gc(true); CUDA.reclaim()
    before=pool_record()
    system=convergence.prepare_case(n,points;target)
    result=Dict{String,Any}("resolution"=>n,"points"=>points,
        "pool_before"=>before,"pool_after"=>pool_record(),
        "modulation"=>modulation_record(system.graph,n,points))
    progress!(result)
    rates=convergence.probe_rates(system,n)
    divisor=convergence.flat_divisor(rates[:,:,1])
    slopes=hcat((convergence.raw_slopes(rates[:,:,i],divisor) for i in axes(rates,3))...)
    size(slopes)==(1432,13) || error("unexpected Copper slope/probe shape")
    result["rates"]=write_array(directory,prefix*"-rates",rates)
    result["slopes"]=write_array(directory,prefix*"-slopes",slopes)
    result["divisor"]=divisor
    result["probes"]=collect(convergence.PROBES)
    progress!(result)
    result["interaction"]=write_array(directory,prefix*"-interaction",interaction64!(system,divisor,convergence))
    result["interaction_poke_m"]=convergence.POKE
    progress!(result)
    result["calibration_captured"]=timings!(system.graph)
    progress!(result)
    system=nothing; GC.gc(true); CUDA.reclaim()
    # Prepare each noisy owner independently: equal seed, physical inputs and
    # initialization. Save full rates and noisy frames for cross-revision parity.
    captured=prepare_case(n,points,target;calibration=false,execution=CapturedGraphExecution())
    first=replay!(captured)
    second=replay!(captured)
    first.rates==second.rates && first.images==second.images || error("captured reset is not deterministic")
    result["captured_reset_rates"]=compare(first.rates,second.rates)
    result["captured_reset_frames"]=compare(first.images,second.images)
    result["atmosphere_rates"]=write_array(directory,prefix*"-atmosphere-rates",first.rates)
    result["atmosphere_frames"]=write_array(directory,prefix*"-atmosphere-frames",first.images)
    progress!(result)
    result["atmosphere_captured"]=timings!(captured)
    progress!(result)
    captured=nothing; GC.gc(true); CUDA.reclaim()
    stream=prepare_case(n,points,target;calibration=false,execution=StreamGraphExecution())
    streamed=replay!(stream)
    result["stream_rates"]=write_array(directory,prefix*"-stream-rates",streamed.rates)
    result["stream_frames"]=write_array(directory,prefix*"-stream-frames",streamed.images)
    result["stream_captured_rates"]=compare(streamed.rates,first.rates)
    result["stream_captured_frames"]=compare(streamed.images,first.images)
    progress!(result)
    result["atmosphere_stream"]=timings!(stream;require_zero=false)
    result["complete"]=true
    progress!(result)
    stream=nothing; GC.gc(true); CUDA.reclaim()
    result
end
function main(args)
    2<=length(args)<=3 || error("usage: HARNESS COPPER_ROOT OUTPUT_DIR [480:32,240:96]")
    copper=realpath(args[1]); directory=abspath(args[2]); mkpath(directory)
    isempty(readdir(directory)) || error("OUTPUT_DIR must be empty")
    realpath(dirname(dirname(pathof(REVOLTCopperSim))))==copper || error("loaded Copper source differs from argument")
    helper=joinpath(copper,"benchmarks","copper_convergence.jl")
    convergence=Main.CopperConvergence
    CUDA.functional() || error("functional CUDA required")
    CUDA.allowscalar(false); BLAS.set_num_threads(1); AOS.Backends.set_fft_provider_threads!(1)
    target=AOS.Backends.compute_device(CUDA.zeros(Float32,1))
    project=something(Base.active_project(),"")
    result=Dict{String,Any}("aos"=>source_record(AOS),"copper"=>source_record(REVOLTCopperSim),
        "julia"=>string(VERSION),"cuda"=>string(pkgversion(CUDA)),"device"=>CUDA.name(CUDA.device()),
        "harness"=>abspath(@__FILE__),"harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "helper_sha256"=>bytes2hex(sha256(read(helper))),"arguments"=>copy(args),
        "command"=>join(Base.julia_cmd().exec," ")*" --project="*project*" --startup-file=no "*join(vcat([abspath(@__FILE__)],args)," "),
        "active_project"=>project,"julia_threads"=>Threads.nthreads(),
        "timing_boundary"=>"host wall elapsed around step_graph!, includes prepared graph completion; warmups32 samples128; evolving noisy atmosphere",
        "cases"=>Dict{String,Any}())
    for filename in (project,joinpath(dirname(project),"Manifest.toml"))
        isfile(filename) && (result[basename(filename)*"_sha256"]=bytes2hex(sha256(read(filename))))
    end
    persist(directory,result)
    result["dense_preparation"]=dense_preparation!(directory)
    persist(directory,result)
    for pair in split(length(args)==3 ? args[3] : "480:32,240:96",',')
        n,points=parse.(Int,split(pair,':'))
        println("Validating ",n,"/",points); flush(stdout)
        progress! = function (partial)
            result["cases"]["n$(n)-m$(points)"] = partial
            persist(directory,result)
        end
        run_case!(directory,n,points,target,convergence,progress!)
        persist(directory,result)
    end
    if haskey(ENV,"AOS_MODULATION_MEMORY_ONLY")
        n,points=parse.(Int,split(ENV["AOS_MODULATION_MEMORY_ONLY"],':'))
        GC.gc(true); CUDA.reclaim(); before=pool_record()
        system=convergence.prepare_case(n,points;target)
        result["memory_only"]=Dict("resolution"=>n,"points"=>points,
            "pool_before"=>before,"pool_after"=>pool_record(),
            "modulation"=>modulation_record(system.graph,n,points))
    end
    result["complete"]=true
    persist(directory,result)
    println("Wrote ",joinpath(directory,"validation.toml"))
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=2 || error("usage: HARNESS COPPER_ROOT OUTPUT_DIR [480:32,240:96]")
    include(joinpath(realpath(ARGS[1]),"benchmarks","copper_convergence.jl"))
    main(ARGS)
end
