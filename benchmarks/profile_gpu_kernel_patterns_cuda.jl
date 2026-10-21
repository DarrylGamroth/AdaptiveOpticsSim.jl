#=
Synthetic diagnostics for existing AOS KernelAbstractions kernels on CUDA.

    julia --project=test/cuda --startup-file=no benchmarks/profile_gpu_kernel_patterns_cuda.jl [output.toml]
    AOS_GPU_PATTERN_PROFILE=1 ncu --target-processes all --profile-from-start off \\
      julia --project=test/cuda --startup-file=no benchmarks/profile_gpu_kernel_patterns_cuda.jl

Production Backends.preferred_workgroupsize selection is preserved on CUDA. Each measured operation is a zero
fill followed by one existing kernel, including kernels that overwrite output.
Timing uses a native CUDA graph containing 32 operations and 64 event samples;
reported seconds are per complete operation. Compilation, inputs, parity checks,
capture and host transfers are excluded. Graph scheduling and the fill are included.
Profile mode emits ONE warmed fill+kernel operation per case, with profiler API
start/stop boundaries, and collects no latency measurements. Case order is printed
and recorded; each profiler region contains only that case's fill and kernel.

CPU execution of the same KA kernel checks backend parity, not independent model
validation. The reduction utility additionally checks a Float64 mathematical sum;
it has no production callers in this source revision. Inputs are synthetic and
these measurements do not predict application throughput or physical validity.
No package launch policy or package implementation is changed.
=#
using AdaptiveOpticsSim
using KernelAbstractions
using CUDA
using TOML
using Dates
using Statistics

const AOS = AdaptiveOpticsSim
const REPEATS = 32
const SAMPLES = 64
const PARITY_ATOL = 2f-5
const PARITY_RTOL = 3f-4

struct KernelFixture{K,O,A,R,D,W}
    label::String
    kernel_name::String
    source_file::String
    kernel::K
    output::O
    inputs::A
    ndrange::R
    dimensions::D
    workgroupsize::W
    production::Bool
end

function fixture(label, kernel, source, output, inputs, ndrange, dimensions;
    production=true)
    return KernelFixture(label, string(nameof(kernel)), source,
        kernel(KernelAbstractions.get_backend(output)), output, inputs,
        ndrange, dimensions, nothing, production)
end

# This is deliberately the whole operation in every timing and profiler region.
function operation!(f::KernelFixture)
    fill!(f.output, zero(eltype(f.output)))
    launch_fixture!(f, f.workgroupsize)
    return nothing
end

launch_fixture!(f::KernelFixture, ::Nothing) =
    f.kernel(f.output, f.inputs...; ndrange=f.ndrange)
launch_fixture!(f::KernelFixture, workgroupsize::Tuple) =
    f.kernel(f.output, f.inputs...; ndrange=f.ndrange, workgroupsize)

hostcopy(x::AbstractArray) = copy(Array(x))
hostcopy(x) = x
readback(x::AbstractArray) = Array(x)
readback(x) = x
devicecopy(x::AbstractArray) = CuArray(x)
devicecopy(x) = x
function device_fixture(f::KernelFixture, kernel)
    output = CUDA.zeros(eltype(f.output), size(f.output)...)
    style = AOS.Backends.execution_style(output)
    workgroupsize = AOS.Backends.preferred_workgroupsize(style, f.ndrange)
    return KernelFixture(f.label, f.kernel_name, f.source_file,
        kernel(KernelAbstractions.get_backend(output)),
        output, map(devicecopy, f.inputs),
        f.ndrange, f.dimensions, workgroupsize, f.production)
end

function field(n, m)
    return [sin(0.017f0 * Float32(i)) + cos(0.013f0 * Float32(j))
        for i in 1:n, j in 1:m]
end
pupil(n) = [(i - (n + 1) / 2)^2 + (j - (n + 1) / 2)^2 <= (0.45n)^2
    for i in 1:n, j in 1:n]

function atmosphere_case(n)
    k = AOS.Atmospheres.moving_layer_replay_kernel!
    f = fixture("atmosphere_n$(n)", k, "src/atmosphere/multilayer.jl",
        zeros(Float32, n, n),
        (field(2n, 2n), Float32[1.25], Float32[2.5], 1, 1f0, 0.7f0,
            false, n, 2n), (n, n), (; n, screen_n=2n))
    return (f, device_fixture(f, k))
end

function dm_cases(n, actuators)
    x = field(n, actuators) ./ Float32(actuators)
    c = field(actuators, actuators)
    y = permutedims(field(n, actuators))
    left = AOS.Optics.dm_separable_left_kernel!
    right = AOS.Optics.dm_separable_right_pupil_kernel!
    dims = (; n, actuators)
    f = fixture("dm_left_n$(n)_k$(actuators)", left,
        "src/optics/deformable_mirror.jl", zeros(Float32, n, actuators),
        (x, c), (n, actuators), dims)
    # Right pass receives a fixed deterministic left product, rather than a
    # changing dependency; this measures each existing kernel independently.
    g = fixture("dm_right_n$(n)_k$(actuators)", right,
        "src/optics/deformable_mirror.jl", zeros(Float32, n, n),
        (x * c, y, pupil(n)), (n, n), dims)
    return ((f, device_fixture(f, left)), (g, device_fixture(g, right)))
end

function slopes_case(n, n_sub)
    sub = n ÷ n_sub
    k = AOS.WavefrontSensors.geometric_slopes_kernel!
    f = fixture("slopes_n$(n)_sub$(n_sub)", k, "src/wfs/geometric_slopes.jl",
        zeros(Float32, 2n_sub^2),
        (field(n, n), fill(true, n_sub, n_sub), sub, n_sub, n_sub^2, 1f0, 1f0),
        (n_sub, n_sub), (; n, n_sub, pixels_per_subaperture=sub))
    return (f, device_fixture(f, k))
end

function detector_cases(n, taps)
    weights = fill(1f0 / Float32(taps), taps)
    radius = taps ÷ 2
    rows = AOS.Detectors.separable_response_rows_kernel!
    cols = AOS.Detectors.separable_response_cols_kernel!
    inputs = (field(n, n), weights, radius, n, n, taps)
    f = fixture("detector_rows_n$(n)_taps$(taps)", rows,
        "src/detectors/interface.jl", zeros(Float32, n, n), inputs,
        (n, n), (; n, taps))
    g = fixture("detector_cols_n$(n)_taps$(taps)", cols,
        "src/detectors/interface.jl", zeros(Float32, n, n), map(hostcopy, inputs),
        (n, n), (; n, taps))
    return ((f, device_fixture(f, rows)), (g, device_fixture(g, cols)))
end

function pyramid_case(n, batches)
    stack = [ComplexF32(sin(0.01f0 * i), cos(0.02f0 * (j + b)))
        for i in 1:n, j in 1:n, b in 1:batches]
    k = AOS.WavefrontSensors.pyramid_modulation_batch_intensity_kernel!
    f = fixture("pyramid_n$(n)_batches$(batches)", k,
        "src/wfs/pyramid/setup.jl", zeros(Float32, n, n),
        (stack, 1f0 / Float32(batches), n, batches), (n, n), (; n, batches))
    return (f, device_fixture(f, k))
end

function reduction_case(n)
    # Binary fractions keep this diagnostic independent of summation order.
    values = [1f0 + Float32(mod(i + j, 8)) / 16f0 for i in 1:n, j in 1:n]
    k = AOS.Backends.masked_sum2d_kernel!
    f = fixture("utility_masked_sum_n$(n)", k, "src/backends/reductions.jl",
        zeros(Float32, 1), (values, pupil(n), 0, 0, n, n), (n, n), (; n);
        production=false)
    return (f, device_fixture(f, k))
end

function check_inputs(f, originals)
    for (input, original) in zip(f.inputs, originals)
        isequal(readback(input), original) || error("$(f.label): input mutated")
    end
end

# Keep serialization cold and concretely typed, including heterogeneous fields.
function write_fields(io, fields::NamedTuple)
    for (key, value) in pairs(fields)
        TOML.print(io, Dict(String(key) => value))
    end
end

function assess!(io, pair, profile)
    cpu, gpu = pair
    originals = map(hostcopy, cpu.inputs)
    operation!(cpu)
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(cpu.output))
    expected = copy(cpu.output)
    for _ in 1:32
        operation!(gpu)
    end
    CUDA.synchronize()
    actual = Array(gpu.output)
    all(isapprox.(actual, expected; atol=PARITY_ATOL, rtol=PARITY_RTOL)) ||
        error("$(cpu.label): CPU/CUDA parity failed")
    parity_error = maximum(abs.(Float64.(actual) .- Float64.(expected)))
    if !cpu.production
        values, mask = cpu.inputs[1:2]
        oracle = sum(Float64(values[i]) for i in eachindex(values) if mask[i])
        isapprox(only(actual), oracle; atol=PARITY_ATOL, rtol=PARITY_RTOL) ||
            error("$(cpu.label): mathematical reduction check failed")
    end
    println("case=$(cpu.label) source=$(cpu.kernel_name)")
    times = Float64[]
    if profile
        CUDA.@profile external=true begin
            operation!(gpu)
            CUDA.synchronize()
        end
    else
        graph = CUDA.capture() do
            for _ in 1:REPEATS
                operation!(gpu)
            end
        end
        executable = CUDA.instantiate(graph)
        CUDA.launch(executable)
        CUDA.synchronize()
        start_event, stop_event = CUDA.CuEvent(), CUDA.CuEvent()
        for _ in 1:SAMPLES
            CUDA.record(start_event)
            CUDA.launch(executable)
            CUDA.record(stop_event)
            CUDA.synchronize(stop_event)
            push!(times, CUDA.elapsed(start_event, stop_event) / REPEATS)
        end
    end
    check_inputs(cpu, originals)
    check_inputs(gpu, originals)
    isequal(cpu.output, expected) || error("CPU reference mutated")
    all(isapprox.(Array(gpu.output), expected; atol=PARITY_ATOL, rtol=PARITY_RTOL)) ||
        error("$(cpu.label): repeated execution is not idempotent")
    println(io, "\n[[cases]]")
    write_fields(io, (; label=cpu.label, kernel_name=cpu.kernel_name,
        source_file=cpu.source_file, production_caller=cpu.production,
        operation="fill_zero_then_one_kernel", workgroupsize=collect(gpu.workgroupsize),
        parity_max_abs=parity_error,
        parity_atol=PARITY_ATOL, parity_rtol=PARITY_RTOL,
        event_seconds_per_operation=times))
    isempty(times) || write_fields(io, (; median_seconds=median(times),
        min_seconds=minimum(times), max_seconds=maximum(times)))
    println(io, "[cases.dimensions]")
    write_fields(io, cpu.dimensions)
    flush(io)
    return nothing
end

function main()
    CUDA.functional() || error("A functional CUDA device is required")
    CUDA.allowscalar(false)
    profile = get(ENV, "AOS_GPU_PATTERN_PROFILE", "0") == "1"
    output = isempty(ARGS) ? joinpath(pwd(), "gpu_kernel_patterns_cuda.toml") : only(ARGS)
    cases = (atmosphere_case(512), dm_cases(256, 32)..., dm_cases(512, 64)...,
        slopes_case(512, 32), slopes_case(1024, 32), detector_cases(512, 7)...,
        detector_cases(512, 31)..., pyramid_case(256, 8), pyramid_case(256, 32),
        reduction_case(128))
    CUDA.synchronize() # All cold input uploads complete before warmup/capture.
    open(output, "w") do io
        write_fields(io, (; schema_version=1, timestamp_utc=string(now(UTC)),
            julia_version=string(VERSION), cuda_version=string(pkgversion(CUDA)),
            device_name=CUDA.name(CUDA.device()),
            compute_capability=string(CUDA.capability(CUDA.device())),
            driver_version=string(CUDA.driver_version()),
            runtime_version=string(CUDA.runtime_version()),
            julia_threads=Threads.nthreads(), profiler_mode=profile,
            repeats_per_graph=REPEATS, event_samples=profile ? 0 : SAMPLES,
            launch_policy="Production Backends.preferred_workgroupsize; selected before warmup",
            timing_boundary="CUDA events around native graph replay; per reset+kernel operation",
            parity_scope="CPU execution of same KA kernels; backend parity, not independent physical validation"))
        for pair in cases
            assess!(io, pair, profile)
        end
    end
    println("results=$output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
