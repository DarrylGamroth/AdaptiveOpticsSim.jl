using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

AdaptiveOpticsSim.set_fft_provider_threads!(1)
BLAS.set_num_threads(parse(Int, get(ENV, "AOS_CONTROL_BLAS_THREADS", "1")))

function measure_operation!(f!; warmup::Int, samples::Int, trials::Int)
    for _ in 1:warmup
        f!()
    end
    trial_medians = Vector{Float64}(undef, trials)
    trial_p95 = similar(trial_medians)
    trial_p99 = similar(trial_medians)
    for trial in 1:trials
        GC.gc()
        latencies = Vector{UInt64}(undef, samples)
        @inbounds for i in eachindex(latencies)
            start_ns = time_ns()
            f!()
            latencies[i] = time_ns() - start_ns
        end
        values = Float64.(latencies)
        trial_medians[trial] = median(values)
        trial_p95[trial] = quantile(values, 0.95)
        trial_p99[trial] = quantile(values, 0.99)
    end
    allocated_bytes = @allocated f!()
    return (
        trial_medians_ns=trial_medians,
        median_ns=median(trial_medians),
        median_p95_ns=median(trial_p95),
        median_p99_ns=median(trial_p99),
        allocated_bytes=allocated_bytes,
    )
end

operator_storage_bytes(operator) = sum(
    Base.summarysize,
    AdaptiveOpticsSim.runtime_reconstructor_storage(operator),
)

function print_result(label, result, storage_bytes)
    println(
        label,
        ": trial_medians_ns=", round.(Int, result.trial_medians_ns),
        " median=", round(Int, result.median_ns), " ns",
        " p95=", round(Int, result.median_p95_ns), " ns",
        " p99=", round(Int, result.median_p99_ns), " ns",
        " alloc=", result.allocated_bytes, " B",
        " storage=", storage_bytes, " B",
    )
end

function main()
    n_slopes = parse(Int, get(ENV, "AOS_CONTROL_SLOPES", "512"))
    n_commands = parse(Int, get(ENV, "AOS_CONTROL_COMMANDS", "256"))
    rank = parse(Int, get(ENV, "AOS_CONTROL_RANK", "32"))
    warmup = parse(Int, get(ENV, "AOS_CONTROL_WARMUP", "10"))
    samples = parse(Int, get(ENV, "AOS_CONTROL_SAMPLES", "200"))
    trials = parse(Int, get(ENV, "AOS_CONTROL_TRIALS", "3"))
    samples >= 100 || error("AOS_CONTROL_SAMPLES must be at least 100 for p99")
    0 < rank <= min(n_slopes, n_commands) ||
        error("AOS_CONTROL_RANK must be within the matrix rank bound")

    rng = MersenneTwister(20260713)
    D = randn(rng, n_slopes, n_commands)
    input = randn(rng, n_slopes)
    interaction = InteractionMatrix(D, 0.1)
    policy = TSVDInverse(rtol=0.0,
        n_trunc=min(n_slopes, n_commands) - rank)

    dense_build_start = time_ns()
    dense = ModalReconstructor(interaction; policy=policy)
    dense_build_ns = time_ns() - dense_build_start
    factorized_build_start = time_ns()
    factorized = FactorizedReconstructor(interaction; policy=policy)
    factorized_build_ns = time_ns() - factorized_build_start
    controlled = ControlledReconstructor(
        factorized,
        DiscreteIntegratorController(n_commands;
            gain=0.5, tau=0.01);
        dt=1e-3,
    )

    dense_output = zeros(n_commands)
    factorized_output = zeros(n_commands)
    controlled_output = zeros(n_commands)
    reconstruct!(dense_output, dense, input)
    reconstruct!(factorized_output, factorized, input)
    relative_error = norm(dense_output - factorized_output) /
        max(norm(dense_output), eps())

    dense_result = measure_operation!(
        () -> reconstruct!(dense_output, dense, input);
        warmup=warmup,
        samples=samples,
        trials=trials,
    )
    factorized_result = measure_operation!(
        () -> reconstruct!(factorized_output, factorized, input);
        warmup=warmup,
        samples=samples,
        trials=trials,
    )
    controlled_result = measure_operation!(
        () -> reconstruct!(controlled_output, controlled, input);
        warmup=warmup,
        samples=samples,
        trials=trials,
    )

    println("AdaptiveOpticsSim scalable control-operator benchmark")
    println("julia=", VERSION,
        " cpu_threads=", Sys.CPU_THREADS,
        " julia_threads=", Threads.nthreads(),
        " blas_threads=", BLAS.get_num_threads(),
        " slopes=", n_slopes,
        " commands=", n_commands,
        " retained_rank=", rank,
        " warmup=", warmup,
        " samples=", samples,
        " trials=", trials)
    println("dense_build_ns=", dense_build_ns,
        " factorized_build_ns=", factorized_build_ns,
        " relative_output_error=", relative_error)
    print_result("ModalReconstructor", dense_result,
        operator_storage_bytes(dense))
    print_result("FactorizedReconstructor", factorized_result,
        operator_storage_bytes(factorized))
    print_result("ControlledReconstructor", controlled_result,
        operator_storage_bytes(controlled))
end

main()
