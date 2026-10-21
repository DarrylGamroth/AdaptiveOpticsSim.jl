# Run with one unchanged application environment from adaptive-optics-simulations.
# This is a review harness: no transport, RTC service, or dependency mutation.
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using LinearAlgebra, Statistics, TOML
import CUDA
const NVTX = CUDA.CUDATools.Profile.NVTX
const INSTRUMENT = first(ARGS)
if INSTRUMENT == "classic"
    import REVOLTClassicSim
    const Instrument = REVOLTClassicSim
    const RATE_OUTPUT = :shwfs_photon_rate
elseif INSTRUMENT == "copper"
    import REVOLTCopperSim
    const Instrument = REVOLTCopperSim
    const RATE_OUTPUT = :pwfs_photon_rate
elseif INSTRUMENT == "gpi2"
    import GPI2Sim
    const Instrument = GPI2Sim
    const RATE_OUTPUT = :photon_rate
else
    error("expected classic, copper, or gpi2")
end
const PROFILE = get(ENV, "AOS_INSTRUMENT_PROFILE", "0") == "1"
const SAMPLES = parse(Int, get(ENV, "AOS_INSTRUMENT_SAMPLES", "32"))
const OUTPUT = ENV["AOS_INSTRUMENT_RESULTS"]
function fields(io, values::NamedTuple)
    for (key, value) in pairs(values)
        TOML.print(io, Dict{String,typeof(value)}(String(key)=>value))
    end
end
prepare_system(target, execution) = Instrument.prepare_hil_system(; target, execution)
function prepare_calibration(target, execution=StreamGraphExecution())
    if INSTRUMENT == "gpi2"
        return Instrument.prepare_calibration_system(; target, execution,
            parameters=Instrument.GPI2Parameters(noise=false))
    end
    return Instrument.prepare_calibration_system(; target, execution, profile=:grid_gaussian)
end
function cycle!(boundary)
    sequence = step_hil_frame!(boundary)
    adopt_hil_command!(boundary, sequence)
    return nothing
end
function probe!(system)
    opd = system.uncompensated_opd
    n = size(opd, 1)
    host_opd = [2f-8 * (sin(Float32(2pi*i/n)) + cos(Float32(2pi*j/n)))
        for i in 1:n, j in 1:n]
    copyto!(opd, host_opd)
    sequence = step_hil_frame!(system.boundary)
    command = hil_command_buffer(system.boundary)
    fill!(command, 0f0)
    command[1], command[cld(length(command),2)], command[end] = 1f-8, -2f-8, 3f-9
    adopt_hil_command!(system.boundary, sequence)
    cycle!(system.boundary)
    return nothing
end
function numerical_check(target)
    # CPU reference executes before accelerator preparation. These are matching
    # full-resolution deterministic models, not an independent physics oracle.
    cpu = prepare_calibration(AdaptiveOpticsSim.Backends.HostComputeDevice())
    probe!(cpu)
    @assert all(isfinite, hil_frame_buffer(cpu.boundary))
    gpu = prepare_calibration(target)
    probe!(gpu)
    expected = hil_frame_buffer(cpu.boundary)
    actual = hil_frame_buffer(gpu.boundary)
    frame_error = sqrt(sum(abs2, Float64.(actual) .- expected) / sum(abs2, expected))
    expected_rate = Array(graph_output(cpu.graph, RATE_OUTPUT))
    actual_rate = Array(graph_output(gpu.graph, RATE_OUTPUT))
    rate_error = sqrt(sum(abs2, Float64.(actual_rate) .- expected_rate) / sum(abs2, expected_rate))
    @assert frame_error < 1e-3
    @assert rate_error < 1e-3
    @assert all(isfinite, actual)
    @assert all(isfinite, actual_rate)
    captured_error = ""
    captured = try
        prepare_calibration(target, CapturedGraphExecution())
    catch exception
        captured_error = sprint(showerror, exception)
        nothing
    end
    captured_frame_error, captured_rate_error = -1.0, -1.0
    if captured !== nothing
        probe!(captured)
        captured_frame = hil_frame_buffer(captured.boundary)
        captured_rate = Array(graph_output(captured.graph, RATE_OUTPUT))
        captured_frame_error = sqrt(sum(abs2, Float64.(captured_frame) .- actual) / sum(abs2, actual))
        captured_rate_error = sqrt(sum(abs2, Float64.(captured_rate) .- actual_rate) / sum(abs2, actual_rate))
        @assert captured_frame_error < 1e-3
        @assert captured_rate_error < 1e-3
    end
    return (; frame_relative_l2_error=frame_error, rate_relative_l2_error=rate_error,
        pupil_resolution=size(cpu.uncompensated_opd,1), frame_shape=collect(size(actual)),
        captured_calibration_admitted=captured !== nothing,
        captured_frame_relative_l2_error=captured_frame_error,
        captured_rate_relative_l2_error=captured_rate_error,
        capture_preparation_error=captured_error)
end
function measure(f, value)
    for _ in 1:32
        f(value)
    end
    bytes = @allocated f(value)
    memory_before = CUDA.used_memory()
    times = zeros(UInt64, SAMPLES)
    gc_count = 0
    for index in eachindex(times)
        gc_before = Base.gc_num()
        start = time_ns()
        f(value)
        times[index] = time_ns() - start
        gc_count += Base.GC_Diff(Base.gc_num(), gc_before).pause
    end
    memory_after = CUDA.used_memory()
    return (; samples_ns=times, p50_ns=median(times), host_bytes=bytes,
        timed_gc_collections=gc_count, device_used_bytes_before=memory_before,
        device_used_bytes_after=memory_after)
end
function inspect_system(io, system, policy)
    graph = system.graph
    println(io, "\n[[execution]]")
    fields(io, (; policy, actual_policy=string(typeof(graph_execution_policy(graph))),
        captured_node_count=captured_graph_node_count(graph)))
    if !PROFILE
        graph_stats = measure(step_graph!, graph)
        println(io, "\n[execution.graph]")
        fields(io, graph_stats)
        reset_hil_boundary!(system.boundary)
        cycle_stats = measure(cycle!, system.boundary)
        println(io, "\n[execution.hil_cycle]")
        fields(io, cycle_stats)
        reset_hil_boundary!(system.boundary)
        cycle!(system.boundary)
        first_frame = copy(hil_frame_buffer(system.boundary))
        reset_hil_boundary!(system.boundary)
        cycle!(system.boundary)
        @assert hil_frame_buffer(system.boundary) == first_frame
        fields(io, (; seeded_reset_reproducible=true,
            final_frame_sum=Float64(sum(first_frame))))
    else
        # Compile every range body before starting the external profiler.
        for _ in 1:32
            step_graph!(graph)
        end
        reset_hil_boundary!(system.boundary)
        for _ in 1:32
            cycle!(system.boundary)
        end
        CUDA.synchronize()
        CUDA.@profile external=true begin
            NVTX.@range string(INSTRUMENT, "_", policy, "_graph") begin
                for _ in 1:64
                    step_graph!(graph)
                end
            end
            reset_hil_boundary!(system.boundary)
            NVTX.@range string(INSTRUMENT, "_", policy, "_hil_cycle") begin
                for _ in 1:64
                    cycle!(system.boundary)
                end
            end
            CUDA.synchronize()
        end
    end
    flush(io)
    println("completed ", INSTRUMENT, " ", policy)
end
function main()
    CUDA.functional() || error("CUDA hardware required")
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.Backends.set_fft_provider_threads!(1)
    target = AdaptiveOpticsSim.Backends.compute_device(CUDA.zeros(Float32, 1))
    open(OUTPUT, "w") do io
        fields(io, (; instrument=INSTRUMENT, julia=string(VERSION),
            aos_version=string(Base.pkgversion(AdaptiveOpticsSim)),
            aos_source=dirname(dirname(pathof(AdaptiveOpticsSim))),
            instrument_source=dirname(dirname(pathof(Instrument))),
            device=CUDA.name(CUDA.device()), profile=PROFILE,
            endpoint="Completed graph GPU outputs; separate complete HIL cycle stages host frame and adopts same-sequence zero command. No transport or RTC.",
            oracle="Full-resolution noiseless CPU/GPU self parity, deterministic OPD and command probes; not an independent scientific oracle."))
        if !PROFILE
            println(io, "\n[numerical]")
            fields(io, numerical_check(target))
            println("completed full-resolution CPU/GPU check ", INSTRUMENT)
        end
        stream = prepare_system(target, StreamGraphExecution())
        inspect_system(io, stream, "stream")
        captured = try
            prepare_system(target, CapturedGraphExecution())
        catch exception
            # Report preparation failure without disguising a stream fallback
            # as capture. Execution and oracle assertions remain uncaught.
            println(io, "\n[capture_failure]")
            fields(io, (; stage="preparation", type=string(typeof(exception)),
                message=sprint(showerror,exception)))
            println("capture preparation failed ", INSTRUMENT, ": ", sprint(showerror,exception))
            nothing
        end
        if captured !== nothing
            inspect_system(io, captured, "captured")
        end
    end
end
main()
