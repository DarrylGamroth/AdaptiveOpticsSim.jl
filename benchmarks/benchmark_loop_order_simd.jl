using AdaptiveOpticsSim
using Random
using Statistics

function measure_operation!(f!; warmup::Int=5, samples::Int=100,
    trials::Int=3)
    for _ in 1:warmup
        f!()
    end
    trial_medians = Vector{Float64}(undef, trials)
    trial_p95 = similar(trial_medians)
    trial_p99 = similar(trial_medians)
    for trial in 1:trials
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
    return (
        trial_medians_ns=trial_medians,
        median_ns=median(trial_medians),
        p95_ns=median(trial_p95),
        p99_ns=median(trial_p99),
        allocated_bytes=@allocated(f!()),
    )
end

function legacy_conv2d_same!(dest, src, kernel)
    n, m = size(src)
    kh, kw = size(kernel)
    cx = div(kh, 2)
    cy = div(kw, 2)
    norm = sum(kernel)
    T = eltype(src)
    inv_norm = norm == 0 ? one(T) : inv(T(norm))
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for ki in 1:kh, kj in 1:kw
            ii = AdaptiveOpticsSim.symm_index(i + ki - cx - 1, n)
            jj = AdaptiveOpticsSim.symm_index(j + kj - cy - 1, m)
            acc += src[ii, jj] * kernel[ki, kj]
        end
        dest[i, j] = acc * inv_norm
    end
    return dest
end

function legacy_conv2d_same_separable!(dest, tmp, src, row_kernel,
    col_kernel)
    n, m = size(src)
    kr = length(row_kernel)
    kc = length(col_kernel)
    cx = div(kr, 2)
    cy = div(kc, 2)
    T = eltype(src)
    inv_norm = inv(sum(row_kernel) * sum(col_kernel))
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for ki in 1:kr
            ii = AdaptiveOpticsSim.symm_index(i + ki - cx - 1, n)
            acc += src[ii, j] * row_kernel[ki]
        end
        tmp[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for kj in 1:kc
            jj = AdaptiveOpticsSim.symm_index(j + kj - cy - 1, m)
            acc += tmp[i, jj] * col_kernel[kj]
        end
        dest[i, j] = acc * inv_norm
    end
    return dest
end

function legacy_elongation!(dest, src, kernel)
    n1, n2 = size(src)
    half = fld(length(kernel), 2)
    @inbounds for i in 1:n1, j in 1:n2
        acc = zero(eltype(dest))
        for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += src[i, jj] * kernel[k + half + 1]
        end
        dest[i, j] = acc
    end
    return dest
end

function legacy_elongation_stack!(dest, src, kernel)
    n1, n2, n3 = size(src)
    half = fld(length(kernel), 2)
    @inbounds for plane in 1:n3, i in 1:n1, j in 1:n2
        acc = zero(eltype(dest))
        for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += src[i, jj, plane] * kernel[k + half + 1]
        end
        dest[i, j, plane] = acc
    end
    return dest
end

function show_result(name, legacy, candidate, exact)
    speedup = legacy.median_ns / candidate.median_ns
    println(name,
        " legacy_trial_medians_ns=", round.(Int, legacy.trial_medians_ns),
        " current_trial_medians_ns=", round.(Int,
            candidate.trial_medians_ns),
        " legacy_median_ns=", round(Int, legacy.median_ns),
        " current_median_ns=", round(Int, candidate.median_ns),
        " current_p95_ns=", round(Int, candidate.p95_ns),
        " current_p99_ns=", round(Int, candidate.p99_ns),
        " speedup=", round(speedup; digits=2),
        " legacy_alloc=", legacy.allocated_bytes,
        " current_alloc=", candidate.allocated_bytes,
        " bitwise_equal=", exact)
end

function main()
    rng = MersenneTwister(20260713)
    T = Float32
    kernel2d = rand(rng, T, 5, 5)
    row_kernel = rand(rng, T, 9)
    col_kernel = rand(rng, T, 9)

    src_conv = rand(rng, T, 96, 96)
    legacy_conv = similar(src_conv)
    candidate_conv = similar(src_conv)
    legacy_conv2d_same!(legacy_conv, src_conv, kernel2d)
    AdaptiveOpticsSim.conv2d_same!(candidate_conv, src_conv, kernel2d)
    conv_exact = legacy_conv == candidate_conv
    legacy_conv_result = measure_operation!(
        () -> legacy_conv2d_same!(legacy_conv, src_conv, kernel2d);
        samples=40)
    candidate_conv_result = measure_operation!(
        () -> AdaptiveOpticsSim.conv2d_same!(candidate_conv, src_conv,
            kernel2d);
        samples=40)

    src_sep = rand(rng, T, 128, 128)
    legacy_sep = similar(src_sep)
    candidate_sep = similar(src_sep)
    legacy_tmp = similar(src_sep)
    candidate_tmp = similar(src_sep)
    legacy_conv2d_same_separable!(legacy_sep, legacy_tmp, src_sep,
        row_kernel, col_kernel)
    AdaptiveOpticsSim.conv2d_same_separable!(candidate_sep, candidate_tmp,
        src_sep, row_kernel, col_kernel)
    sep_exact = legacy_sep == candidate_sep
    legacy_sep_result = measure_operation!(
        () -> legacy_conv2d_same_separable!(legacy_sep, legacy_tmp,
            src_sep, row_kernel, col_kernel))
    candidate_sep_result = measure_operation!(
        () -> AdaptiveOpticsSim.conv2d_same_separable!(candidate_sep,
            candidate_tmp, src_sep, row_kernel, col_kernel))

    src_elongation = rand(rng, T, 256, 256)
    legacy_elongation = similar(src_elongation)
    candidate_elongation = similar(src_elongation)
    legacy_elongation!(legacy_elongation, src_elongation, row_kernel)
    AdaptiveOpticsSim._apply_elongation!(AdaptiveOpticsSim.ScalarCPUStyle(), src_elongation,
        candidate_elongation, row_kernel, 4, 256, 256)
    elongation_exact = legacy_elongation == candidate_elongation
    legacy_elongation_result = measure_operation!(
        () -> legacy_elongation!(legacy_elongation, src_elongation,
            row_kernel))
    candidate_elongation_result = measure_operation!(
        () -> AdaptiveOpticsSim._apply_elongation!(AdaptiveOpticsSim.ScalarCPUStyle(),
            src_elongation, candidate_elongation, row_kernel, 4, 256, 256))

    src_stack = rand(rng, T, 32, 32, 256)
    legacy_stack = similar(src_stack)
    candidate_stack = similar(src_stack)
    legacy_elongation_stack!(legacy_stack, src_stack, row_kernel)
    AdaptiveOpticsSim._apply_elongation_stack!(AdaptiveOpticsSim.ScalarCPUStyle(), src_stack,
        candidate_stack, row_kernel, 4, 32, 32, 256)
    stack_exact = legacy_stack == candidate_stack
    legacy_stack_result = measure_operation!(
        () -> legacy_elongation_stack!(legacy_stack, src_stack, row_kernel);
        samples=30)
    candidate_stack_result = measure_operation!(
        () -> AdaptiveOpticsSim._apply_elongation_stack!(AdaptiveOpticsSim.ScalarCPUStyle(),
            src_stack, candidate_stack, row_kernel, 4, 32, 32, 256);
        samples=30)

    println("AdaptiveOpticsSim loop-order/native-SIMD comparison")
    println("julia=", VERSION,
        " cpu=", Sys.CPU_NAME,
        " threads=", Threads.nthreads(),
        " element_type=", T)
    show_result("conv2d_same", legacy_conv_result, candidate_conv_result,
        conv_exact)
    show_result("conv2d_same_separable", legacy_sep_result,
        candidate_sep_result, sep_exact)
    show_result("elongation_frame", legacy_elongation_result,
        candidate_elongation_result, elongation_exact)
    show_result("elongation_stack", legacy_stack_result,
        candidate_stack_result, stack_exact)
end

main()
