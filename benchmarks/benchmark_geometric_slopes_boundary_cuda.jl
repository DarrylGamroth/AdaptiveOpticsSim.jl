#=
Standalone benchmark experiment, not a qualified numerical policy or package API.

    julia --project=test/cuda --startup-file=no benchmarks/benchmark_geometric_slopes_boundary_cuda.jl [output.toml]

For unmasked pixels within a valid subaperture, the ordinary geometric finite
 differences telescope in exact real arithmetic:
    sx = Σy (opd[xe,y] - opd[xs,y]) / max(sub*(sub-1), 1)
    sy = Σx (opd[x,ye] - opd[x,ys]) / max(sub*(sub-1), 1)
The experimental kernel below uses this identity without fastmath. Float32
reassociation can change results. This experiment uses the original CPU KA
kernel as a backend/numerical comparison, not independent model validation.
The validity mask selects entire subapertures, as in the original kernel.

Both variants use the production-selected 16×16 workgroups, identical inputs,
zero-fill+kernel timing boundaries, 32-operation native graphs, and 64 CUDA event
samples per block. Four ABBA rounds retain all raw per-operation event times.
Each of three shapes checks a bounded smooth and a bounded nonsmooth OPD field
with mixed valid subapertures. Compilation, capture, warmup, parity and transfers
are outside timing. No production source, dependency, or tests are changed.
=#
include(joinpath(@__DIR__, "profile_gpu_kernel_patterns_cuda.jl"))

@kernel function boundary_geometric_slopes_kernel!(slopes, opd, valid_mask,
    sub::Int, n_sub::Int, offset::Int, scale_x, scale_y)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = i + (j - 1) * n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        sx = zero(eltype(slopes))
        sy = zero(eltype(slopes))
        if @inbounds valid_mask[i, j]
            @inbounds for y in ys:ye
                sx += opd[xe, y] - opd[xs, y]
            end
            @inbounds for x in xs:xe
                sy += opd[x, ye] - opd[x, ys]
            end
            count = max(sub * (sub - 1), 1)
            @inbounds slopes[idx] = scale_x * sx / count
            @inbounds slopes[idx + offset] = scale_y * sy / count
        else
            @inbounds slopes[idx] = zero(eltype(slopes))
            @inbounds slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

smooth_opd(n) = [0.4f0 * sin(0.017f0 * i) + 0.5f0 * cos(0.013f0 * j)
    for i in 1:n, j in 1:n]
nonsmooth_opd(n) = [Float32(mod(37i + 19j + i*j, 257)) / 128f0 - 1f0
    for i in 1:n, j in 1:n]

function boundary_fixtures(n, n_sub, dataset, opd)
    sub = n ÷ n_sub
    mask = [mod(i + 2j, 5) != 0 for i in 1:n_sub, j in 1:n_sub]
    inputs = (opd, mask, sub, n_sub, n_sub^2, 1f0, 1f0)
    dims = (; n, n_sub, pixels_per_subaperture=sub)
    original = AOS.WavefrontSensors.geometric_slopes_kernel!
    cpu = fixture("$(dataset)_n$(n)_l$(n_sub)", original,
        "src/wfs/geometric_slopes.jl", zeros(Float32, 2n_sub^2),
        inputs, (n_sub, n_sub), dims)
    alternative = fixture(cpu.label, boundary_geometric_slopes_kernel!,
        "benchmarks/benchmark_geometric_slopes_boundary_cuda.jl",
        zeros(Float32, 2n_sub^2), map(hostcopy, inputs), (n_sub, n_sub), dims;
        production=false)
    gpu_original = device_fixture(cpu, original)
    gpu_boundary = device_fixture(alternative, boundary_geometric_slopes_kernel!)
    return (; cpu, gpu_original, gpu_boundary, dataset)
end

function captured_operation(f)
    graph = CUDA.capture() do
        for _ in 1:REPEATS
            operation!(f)
        end
    end
    executable = CUDA.instantiate(graph)
    CUDA.launch(executable)
    CUDA.synchronize()
    return (; graph, executable)
end

function event_block(captured, start_event, stop_event)
    samples = Vector{Float64}(undef, SAMPLES)
    for i in eachindex(samples)
        CUDA.record(start_event)
        CUDA.launch(captured.executable)
        CUDA.record(stop_event)
        CUDA.synchronize(stop_event)
        samples[i] = CUDA.elapsed(start_event, stop_event) / REPEATS
    end
    return samples
end

function compare_boundary!(io, f)
    originals = map(hostcopy, f.cpu.inputs)
    operation!(f.cpu)
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(f.cpu.output))
    expected = copy(f.cpu.output)
    for _ in 1:32
        operation!(f.gpu_original)
        operation!(f.gpu_boundary)
    end
    CUDA.synchronize()
    actual_original = Array(f.gpu_original.output)
    actual_boundary = Array(f.gpu_boundary.output)
    for actual in (actual_original, actual_boundary)
        all(isapprox.(actual, expected; atol=PARITY_ATOL, rtol=PARITY_RTOL)) ||
            error("$(f.cpu.label): numerical comparison failed")
    end
    error_original = maximum(abs.(Float64.(actual_original) .- Float64.(expected)))
    error_boundary = maximum(abs.(Float64.(actual_boundary) .- Float64.(expected)))
    difference = maximum(abs.(Float64.(actual_boundary) .- Float64.(actual_original)))
    original_graph = captured_operation(f.gpu_original)
    boundary_graph = captured_operation(f.gpu_boundary)
    start_event, stop_event = CUDA.CuEvent(), CUDA.CuEvent()
    blocks = ntuple(16) do block
        # Each round is original, boundary, boundary, original (ABBA).
        original = mod(block - 1, 4) in (0, 3)
        samples = event_block(original ? original_graph : boundary_graph,
            start_event, stop_event)
        return (; round=(block - 1) ÷ 4 + 1, position=mod(block - 1, 4) + 1,
            variant=original ? "original" : "boundary", samples)
    end
    check_inputs(f.cpu, originals)
    check_inputs(f.gpu_original, originals)
    check_inputs(f.gpu_boundary, originals)
    isequal(f.cpu.output, expected) || error("CPU reference mutated")
    for gpu in (f.gpu_original, f.gpu_boundary)
        all(isapprox.(Array(gpu.output), expected;
            atol=PARITY_ATOL, rtol=PARITY_RTOL)) || error("Repeated output drift")
    end
    original_times = reduce(vcat, map(b -> b.samples, filter(b -> b.variant == "original", blocks)))
    boundary_times = reduce(vcat, map(b -> b.samples, filter(b -> b.variant == "boundary", blocks)))
    println(io, "\n[[cases]]")
    write_fields(io, (; label=f.cpu.label, dataset=f.dataset,
        original_kernel=f.cpu.kernel_name, original_source=f.cpu.source_file,
        boundary_kernel=f.gpu_boundary.kernel_name,
        boundary_source=f.gpu_boundary.source_file,
        original_vs_cpu_max_abs=error_original,
        boundary_vs_cpu_max_abs=error_boundary, boundary_vs_original_max_abs=difference,
        comparison_atol=PARITY_ATOL, comparison_rtol=PARITY_RTOL,
        original_median_seconds=median(original_times),
        boundary_median_seconds=median(boundary_times),
        median_ratio_original_over_boundary=median(original_times)/median(boundary_times),
        workgroupsize=collect(f.gpu_original.workgroupsize)))
    println(io, "[cases.dimensions]")
    write_fields(io, f.cpu.dimensions)
    for block in blocks
        println(io, "[[cases.blocks]]")
        write_fields(io, (; round=block.round, position=block.position,
            variant=block.variant, event_seconds_per_operation=block.samples))
    end
    flush(io)
    println("$(f.cpu.label): original=$(median(original_times))s boundary=$(median(boundary_times))s max_difference=$difference")
    return nothing
end

function boundary_main()
    CUDA.functional() || error("A functional CUDA device is required")
    CUDA.allowscalar(false)
    output = isempty(ARGS) ? "geometric_slopes_boundary_cuda.toml" : only(ARGS)
    shapes = ((64, 64), (512, 32), (1024, 32))
    cases = map(shapes) do (n, n_sub)
        (boundary_fixtures(n, n_sub, "smooth", smooth_opd(n)),
            boundary_fixtures(n, n_sub, "nonsmooth", nonsmooth_opd(n)))
    end
    CUDA.synchronize() # Complete all input writes before any measurements.
    open(output, "w") do io
        write_fields(io, (; schema_version=1, timestamp_utc=string(now(UTC)),
            julia_version=string(VERSION), cuda_version=string(pkgversion(CUDA)),
            device_name=CUDA.name(CUDA.device()),
            compute_capability=string(CUDA.capability(CUDA.device())),
            driver_version=string(CUDA.driver_version()),
            runtime_version=string(CUDA.runtime_version()),
            julia_threads=Threads.nthreads(), repeats_per_graph=REPEATS,
            samples_per_block=SAMPLES, paired_order="ABBA", rounds=4,
            launch_policy="Production Backends.preferred_workgroupsize; unchanged",
            timing_boundary="CUDA events around native graph replay; per zero-fill+kernel operation",
            qualification="Benchmark-only algebraic reassociation experiment; numerical policy not qualified",
            comparison_scope="Original CPU KA kernel; backend/numerical comparison, not independent physical validation"))
        for shape_cases in cases, f in shape_cases
            compare_boundary!(io, f)
        end
    end
    println("results=$output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    boundary_main()
end
