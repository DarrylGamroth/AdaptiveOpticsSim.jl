# Benchmark-only mask locality probe. Usage: julia --project=ENV FILE OUTPUT.toml
using AdaptiveOpticsSim, KernelAbstractions
using Random, Statistics, LinearAlgebra, TOML, SHA
import CUDA
const MaskProbeWFS = AdaptiveOpticsSim.WavefrontSensors

@kernel function shared_transposed_mask_kernel!(stack, mask, side::Int, batch::Int)
    x, y = @index(Local, NTuple)
    bx, by = @index(Group, NTuple)
    tile = @localmem eltype(mask) (33, 32)
    original_x = (by - 1) * 32 + x
    original_y = (bx - 1) * 32 + y
    for offset in 0:8:24
        if original_x <= side && original_y + offset <= side
            @inbounds tile[x, y + offset] = mask[original_x, original_y + offset]
        end
    end
    @synchronize
    i = (bx - 1) * 32 + x
    j = (by - 1) * 32 + y
    for offset in 0:8:24
        if i <= side && j + offset <= side
            mask_value = @inbounds tile[y + offset, x]
            @inbounds for b in 1:batch
                stack[i, j + offset, b] *= mask_value
            end
        end
    end
end

function launch_mask_probe!(::Val{:current}, stack, mask)
    n = size(stack, 1)
    MaskProbeWFS.pyramid_modulation_transposed_mask_kernel!(get_backend(stack))(
        stack, mask, n, size(stack, 3); ndrange=(n, n), workgroupsize=(16, 16))
    return nothing
end

function launch_mask_probe!(::Val{:shared}, stack, mask)
    n = size(stack, 1)
    shared_transposed_mask_kernel!(get_backend(stack), (32, 8))(
        stack, mask, n, size(stack, 3); ndrange=(cld(n, 32) * 32, cld(n, 32) * 8))
    return nothing
end

function capture_mask_probe(variant, stack, mask)
    launch_mask_probe!(variant, stack, mask)
    CUDA.synchronize(; blocking=true)
    graph = CUDA.capture() do
        launch_mask_probe!(variant, stack, mask)
    end
    return CUDA.instantiate(graph)
end

function small_mask_parity()
    rng = Xoshiro(0x72193)
    records = Dict{String,Any}[]
    for T in (ComplexF32, ComplexF64), n in (17, 33), batch in (1, 3)
        host = rand(rng, T, n, n, batch)
        host_mask = rand(rng, T, n, n)
        expected = host .* reshape(permutedims(host_mask, (2, 1)), n, n, 1)
        original = CUDA.CuArray(host)
        mask = CUDA.CuArray(host_mask)
        current = copy(original)
        shared = copy(original)
        executables = (capture_mask_probe(Val(:current), current, mask),
            capture_mask_probe(Val(:shared), shared, mask))
        # Capture warmup has mutated the outputs; restore before the replay.
        copyto!(current, original)
        copyto!(shared, original)
        CUDA.synchronize(; blocking=true)
        for executable in executables
            CUDA.launch(executable)
        end
        CUDA.synchronize(; blocking=true)
        old = Array(current)
        new = Array(shared)
        bitwise = reinterpret(UInt8, vec(old)) == reinterpret(UInt8, vec(new))
        bitwise || error("mask kernels differ for $T/$n/$batch")
        relative = norm(new - expected) / norm(expected)
        tolerance = 2eps(real(T))
        relative <= tolerance || error("CPU mask parity failed: $relative > $tolerance")
        push!(records, Dict("type"=>string(T), "side"=>n, "batch"=>batch,
            "old_new_bitwise_equal"=>bitwise, "cpu_relative_l2"=>relative,
            "cpu_relative_l2_bound"=>tolerance))
    end
    return records
end

function prepare_large_mask_probe()
    n, batch = 2880, 16
    rng = Xoshiro(0x834271)
    host = rand(rng, ComplexF32, n, n, batch)
    input_hash = bytes2hex(sha256(reinterpret(UInt8, vec(host))))
    original = CUDA.CuArray(host)
    mask = CUDA.CuArray(rand(rng, ComplexF32, n, n))
    stack = similar(original)
    reference = similar(original)
    return (; original, stack, reference, mask, input_hash, n, batch)
end

function mask_probe_main(args)
    length(args) == 1 || error("Usage: OUTPUT.toml")
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Threads.nthreads() == 1 || error("use one Julia thread")
    small = small_mask_parity()
    GC.gc(true)
    CUDA.reclaim()
    f = prepare_large_mask_probe()
    GC.gc(true)
    copyto!(f.stack, f.original)
    CUDA.synchronize(; blocking=true)
    variants = ((name="current", executable=capture_mask_probe(Val(:current), f.stack, f.mask)),
        (name="shared", executable=capture_mask_probe(Val(:shared), f.stack, f.mask)))
    copyto!(f.reference, f.original)
    launch_mask_probe!(Val(:current), f.reference, f.mask)
    CUDA.synchronize(; blocking=true)
    parity = Dict{String,Any}[]
    for variant in variants
        copyto!(f.stack, f.original)
        CUDA.synchronize(; blocking=true)
        CUDA.launch(variant.executable)
        CUDA.synchronize(; blocking=true)
        equal = all(isequal.(f.stack, f.reference))
        equal || error("large mask parity failed for $(variant.name)")
        push!(parity, Dict("variant"=>variant.name, "all_elements_equal"=>equal))
    end
    GC.gc(true)
    CUDA.reclaim()
    rng = Xoshiro(0x921346)
    trials = Dict{String,Any}[]
    for round in 1:2
        order = randperm(rng, length(variants))
        for (position, index) in enumerate(order)
            variant = variants[index]
            for _ in 1:8
                copyto!(f.stack, f.original)
                CUDA.synchronize(; blocking=true)
                CUDA.launch(variant.executable)
                CUDA.synchronize(; blocking=true)
            end
            samples = zeros(Float64, 32)
            for i in eachindex(samples)
                # All variants see identical data. This D2D restore and its
                # completion boundary are deliberately outside the timer.
                copyto!(f.stack, f.original)
                CUDA.synchronize(; blocking=true)
                start = time_ns()
                CUDA.launch(variant.executable)
                CUDA.synchronize(; blocking=true)
                samples[i] = (time_ns() - start) / 1e6
            end
            push!(trials, Dict("round"=>round, "position"=>position,
                "variant"=>variant.name, "samples_ms"=>samples,
                "median_ms"=>median(samples), "minimum_ms"=>minimum(samples),
                "p95_ms"=>quantile(samples, .95)))
        end
    end
    array_bytes = sizeof(eltype(f.original)) * length(f.original)
    record = Dict("revision"=>readchomp(`git -C $(dirname(@__DIR__)) rev-parse HEAD`),
        "harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "julia"=>string(VERSION), "runtime"=>sprint(io->CUDA.versioninfo(io)),
        "gpu"=>CUDA.name(CUDA.device()), "side"=>f.n, "batch"=>f.batch,
        "original_bytes"=>array_bytes, "stack_bytes"=>array_bytes,
        "reference_bytes"=>array_bytes, "mask_bytes"=>sizeof(eltype(f.mask))*length(f.mask),
        "device_used_bytes"=>CUDA.used_memory(), "input_sha256"=>f.input_hash,
        "small_seed"=>0x72193, "large_seed"=>0x834271, "order_seed"=>0x921346,
        "small_parity"=>small, "large_parity"=>parity, "trials"=>trials,
        "warmups_per_trial"=>8, "samples_per_trial"=>32,
        "boundary"=>"Captured mask kernel launch through blocking completion; D2D input restore and wait outside every sample; unpaced single outstanding call",
        "limitations"=>"Kernel microbenchmark; restore changes cache state; full FFT pipeline validation required before production selection")
    mkpath(dirname(abspath(args[1])))
    open(io->TOML.print(io, record; sorted=true), args[1], "w")
    for trial in trials
        println("round=$(trial["round"]) $(trial["variant"]) median_ms=$(trial["median_ms"])")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mask_probe_main(ARGS)
end
