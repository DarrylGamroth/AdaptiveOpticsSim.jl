module AdaptiveOpticsSimAMDGPUExt

using AdaptiveOpticsSim
using AMDGPU
using AbstractFFTs
using LinearAlgebra
using Random

#
# AMDGPU backend extension
#
# This extension supplies backend-native dense linear algebra and FFT plumbing
# for the maintained ROCArray execution paths. The main mathematical surfaces
# implemented here are:
#
# - pseudoinverse construction from SVD
# - stable Hermitian right division used by tomography/calibration
# - normal-equation solves and SVD fallback for LiFT
#
# The key rule is that the algorithms match the core implementation, while the
# execution is specialized to rocBLAS / rocSOLVER / rocFFT where that improves
# performance or avoids host fallback.
#
AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.ROCArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = :amdgpu
AdaptiveOpticsSim.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)
AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AMDGPU.rocFFT.ROCFFTPlan) = (plan * buffer; buffer)
AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AbstractFFTs.ScaledPlan) = (plan * buffer; buffer)
AdaptiveOpticsSim.default_build_backend(::AMDGPU.ROCArray) = AdaptiveOpticsSim.GPUArrayBuildBackend(AdaptiveOpticsSim.AMDGPUBackendTag)
AdaptiveOpticsSim.prepare_build_matrix(::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag}, A::AbstractMatrix) = Matrix(A)
AdaptiveOpticsSim.grouped_accumulation_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.PyramidWFS},
) = AdaptiveOpticsSim.GroupedStaged2DPlan()
AdaptiveOpticsSim.grouped_accumulation_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.BioEdgeWFS},
) = AdaptiveOpticsSim.GroupedStaged2DPlan()
AdaptiveOpticsSim.sh_sensing_execution_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.ShackHartmann},
) = AdaptiveOpticsSim.ShackHartmannRocmSafePlan()
AdaptiveOpticsSim.detector_execution_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detector},
) = AdaptiveOpticsSim.DetectorHostMirrorPlan()
AdaptiveOpticsSim.reduction_execution_plan(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AMDGPU.ROCArray,
) = AdaptiveOpticsSim.HostMirrorReductionPlan()
AdaptiveOpticsSim.randn_backend_async!(::AdaptiveOpticsSim.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
AdaptiveOpticsSim._randn_backend!(::AdaptiveOpticsSim.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
function AdaptiveOpticsSim._randn_frame_noise!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    out::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    AdaptiveOpticsSim.randn_backend!(rng, out)
    return out
end
function AdaptiveOpticsSim._randn_frame_noise!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    AdaptiveOpticsSim.randn_backend!(rng, cube)
    return cube
end
function AdaptiveOpticsSim.randn_phase_noise!(rng::AbstractRNG, out::AMDGPU.ROCArray{T,2}, host::Matrix{T}) where {T<:AbstractFloat}
    if size(host) != size(out)
        host = Matrix{T}(undef, size(out)...)
    end
    randn!(rng, host)
    copyto!(out, host)
    return host
end
function roc_sh_spot_cutoff_stats_kernel!(stats, spot_cube, valid_mask, cutoff, n_sub::Int, n1::Int, n2::Int)
    lane = Int(AMDGPU.workitemIdx().x)
    idx = Int(AMDGPU.workgroupIdx().x)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        if @inbounds valid_mask[i, j]
            T = eltype(stats)
            total_sh = AMDGPU.@ROCStaticLocalArray(T, 256, false)
            sx_sh = AMDGPU.@ROCStaticLocalArray(T, 256, false)
            sy_sh = AMDGPU.@ROCStaticLocalArray(T, 256, false)
            total = zero(T)
            sx = zero(T)
            sy = zero(T)
            npix = n1 * n2
            p = lane
            while p <= npix
                x = (p - 1) ÷ n2 + 1
                y = p - (x - 1) * n2
                val = @inbounds spot_cube[idx, x, y]
                if val >= cutoff
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
                p += 256
            end
            @inbounds begin
                total_sh[lane] = total
                sx_sh[lane] = sx
                sy_sh[lane] = sy
            end
            AMDGPU.sync_workgroup()
            stride = 128
            while stride >= 1
                if lane <= stride
                    @inbounds begin
                        total_sh[lane] += total_sh[lane + stride]
                        sx_sh[lane] += sx_sh[lane + stride]
                        sy_sh[lane] += sy_sh[lane + stride]
                    end
                end
                AMDGPU.sync_workgroup()
                stride ÷= 2
            end
            if lane == 1
                base = 3 * (idx - 1)
                @inbounds begin
                    stats[base + 1] = total_sh[1]
                    stats[base + 2] = sx_sh[1]
                    stats[base + 3] = sy_sh[1]
                end
            end
        end
    end
    return
end
function AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    wfs::AdaptiveOpticsSim.ShackHartmann,
    cutoff::T,
) where {T<:AbstractFloat}
    n_sub = wfs.params.n_subap
    offset = n_sub * n_sub
    n1 = size(wfs.state.spot_cube, 2)
    n2 = size(wfs.state.spot_cube, 3)
    fill!(wfs.state.spot_stats, zero(eltype(wfs.state.spot_stats)))
    AMDGPU.@roc groupsize=256 gridsize=offset roc_sh_spot_cutoff_stats_kernel!(
        wfs.state.spot_stats,
        wfs.state.spot_cube,
        wfs.state.valid_mask,
        cutoff,
        n_sub,
        n1,
        n2,
    )
    AMDGPU.synchronize()
    style = AdaptiveOpticsSim.execution_style(wfs.state.slopes)
    AdaptiveOpticsSim.launch_kernel!(style, AdaptiveOpticsSim.sh_finalize_spot_slopes_kernel!, wfs.state.slopes,
        wfs.state.spot_stats, wfs.state.valid_mask, n_sub, offset; ndrange=offset)
    return wfs.state.slopes
end
AdaptiveOpticsSim.backend_matmul(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'N', A, B)
AdaptiveOpticsSim.backend_matmul_transpose_right(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'T', A, B)

function dense_copy_to_roc(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    out = AMDGPU.ROCArray{T}(undef, size(A)...)
    copyto!(out, A)
    return out
end

function roc_svd(A::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    F = copy(A)
    U, S, Vt = AMDGPU.rocSOLVER.gesvd!('S', 'S', F)
    return (; U, S, Vt, s_host=AdaptiveOpticsSim.singular_values_host(S))
end

function roc_lu_solve!(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    n = size(A, 1)
    ipiv = AMDGPU.ROCArray{Cint}(undef, n)
    AMDGPU.rocSOLVER.getrf!(A, ipiv)
    AMDGPU.rocSOLVER.getrs!('N', A, ipiv, B)
    return B
end

function roc_cholesky_solve!(
    A::AMDGPU.ROCArray{T,2},
    rhs::AMDGPU.ROCArray{T,N};
    check::Bool=false,
) where {T<:AbstractFloat,N}
    rhs_mat = if N == 1
        out = AMDGPU.ROCArray{T}(undef, length(rhs), 1)
        copyto!(vec(out), rhs)
        out
    else
        rhs
    end
    chol = cholesky!(Hermitian(A), check=check)
    if issuccess(chol)
        ldiv!(chol, rhs_mat)
    end
    if N == 1
        copyto!(rhs, vec(rhs_mat))
        return rhs, chol
    end
    return rhs_mat, chol
end

"""
    pseudoinverse_from_roc_svd(backend, U, S, Vt, inv_s_host)

Assemble the pseudoinverse `V * Diagonal(inv_s) * U'` from the compact rocSOLVER
SVD factors.

`rocSOLVER.gesvd!` returns `Vt`, so the final matrix product is expressed with
transpose flags rather than materialized transposes.
"""
function pseudoinverse_from_roc_svd(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    U::AMDGPU.ROCArray{T,2}, S::AMDGPU.ROCArray{T,1}, Vt::AMDGPU.ROCArray{T,2},
    inv_s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s = AdaptiveOpticsSim.materialize_build(backend, S, inv_s_host)
    U_scaled = copy(U)
    U_scaled .*= reshape(inv_s, 1, :)
    # Mathematically this is V * U_scaled', but rocSOLVER returns Vt and the
    # arrays are already column-major, so transpose flags are enough here.
    return AMDGPU.rocBLAS.gemm('T', 'T', Vt, U_scaled)
end

function inverse_scaling_and_stats(::AdaptiveOpticsSim.ExactPseudoInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s_host = similar(s_host)
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = iszero(s_host[i]) ? zero(T) : inv(s_host[i])
    end
    effective_rank = count(!iszero, s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function inverse_scaling_and_stats(policy::AdaptiveOpticsSim.TSVDInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.n_trunc >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("TSVD n_trunc must be >= 0"))
    isempty(s_host) && return similar(s_host), AdaptiveOpticsSim.InverseStats(s_host, T(Inf), 0, 0)
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s_host)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s_host = similar(s_host)
    fill!(inv_s_host, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s_host[i] = inv(s_host[i])
    end
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, length(s_host) - effective_rank)
end

function inverse_scaling_and_stats(policy::AdaptiveOpticsSim.TikhonovInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.lambda >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("Tikhonov lambda must be >= 0"))
    inv_s_host = similar(s_host)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = s_host[i] / (s_host[i]^2 + λ2)
    end
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s_host)
    denom = max(isempty(s_host) ? zero(T) : s_host[end], T(policy.lambda))
    cond = (isempty(s_host) || iszero(denom)) ? T(Inf) : s_host[begin] / denom
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function roc_inverse_operator(
    backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2},
    policy::AdaptiveOpticsSim.InversePolicy,
) where {T<:AbstractFloat}
    svd_parts = roc_svd(A)
    inv_s_host, stats = inverse_scaling_and_stats(policy, svd_parts.s_host)
    if isempty(svd_parts.s_host)
        empty_inverse = AdaptiveOpticsSim.materialize_build(backend, similar(A, T, size(A, 2), size(A, 1)))
        return empty_inverse, stats
    end
    M = pseudoinverse_from_roc_svd(backend, svd_parts.U, svd_parts.S, svd_parts.Vt, inv_s_host)
    return M, stats
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, ::AdaptiveOpticsSim.ExactPseudoInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, AdaptiveOpticsSim.ExactPseudoInverse())
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TSVDInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, policy)
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TikhonovInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, policy)
end

"""
    stable_hermitian_right_division(_, rhs, gram)

Solve the right-division `rhs / gram` through a left solve on the transposed
system.

The preferred path is Cholesky on the Hermitian Gram matrix. If that fails, the
implementation falls back to LU so the higher-level algorithm remains robust on
ill-conditioned runtime/calibration cases.
"""
function AdaptiveOpticsSim.stable_hermitian_right_division(
    _backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    rhs::AMDGPU.ROCArray{T,2},
    gram::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    gram_factor = copy(gram)
    rhs_t = permutedims(rhs, (2, 1))
    _, fact = roc_cholesky_solve!(gram_factor, rhs_t; check=false)
    if issuccess(fact)
        return permutedims(rhs_t, (2, 1))
    end
    gram_lu = copy(gram)
    roc_lu_solve!(gram_lu, rhs_t)
    return permutedims(rhs_t, (2, 1))
end

function AdaptiveOpticsSim.solve_lift_fallback!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T},
    rhs::AMDGPU.ROCArray{T,1}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::AdaptiveOpticsSim.LiFTDampingMode) where {T<:AbstractFloat}
    H_mat = dense_copy_to_roc(H)
    svd_parts = roc_svd(H_mat)
    λ = damping isa AdaptiveOpticsSim.LiFTLevenbergMarquardt ?
        AdaptiveOpticsSim.damping_lambda(damping, adjoint(H_mat) * H_mat) : zero(T)
    work = AMDGPU.ROCArray{T}(undef, length(svd_parts.S))
    mul!(work, transpose(svd_parts.U), residual)
    @. work = ifelse(iszero(svd_parts.S^2 + λ^2), zero(T), (svd_parts.S * work) / (svd_parts.S^2 + λ^2))
    mul!(rhs, adjoint(svd_parts.Vt), work)
    return rhs
end

"""
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)

Solve the LiFT normal equations on ROCArray inputs.

The main path uses Cholesky on the normal matrix, optionally with diagonal
loading from the damping policy. If repeated factorization attempts fail, the
code falls back to the SVD-based Levenberg-Marquardt solve to preserve the same
robustness guarantees as the CPU implementation.
"""
function AdaptiveOpticsSim.solve_normal_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    ::AdaptiveOpticsSim.LiFTDampingNone) where {T<:AbstractFloat}
    factor_mat = factor isa AMDGPU.ROCArray{T,2} ? factor : dense_copy_to_roc(factor)
    copyto!(factor_mat, normal)
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = AdaptiveOpticsSim.regularization_load(normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if !issuccess(chol)
            λ *= T(10)
            @views factor_mat[diagind(factor_mat)] .+= λ
            _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
            if !issuccess(chol)
                return AdaptiveOpticsSim.solve_lift_fallback!(diag, rhs, H, residual,
                    AdaptiveOpticsSim.LiFTLevenbergMarquardt(lambda0=λ))
            end
        end
    end
    diag.regularization = λ
    return rhs
end

function AdaptiveOpticsSim.solve_normal_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::AdaptiveOpticsSim.LiFTLevenbergMarquardt) where {T<:AbstractFloat}
    factor_mat = factor isa AMDGPU.ROCArray{T,2} ? factor : dense_copy_to_roc(factor)
    copyto!(factor_mat, normal)
    λ = AdaptiveOpticsSim.damping_lambda(damping, normal)
    if λ > zero(T)
        @views factor_mat[diagind(factor_mat)] .+= λ
    end
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), AdaptiveOpticsSim.regularization_load(normal))
        copyto!(factor_mat, normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if λ > T(1e12)
            return AdaptiveOpticsSim.solve_lift_fallback!(diag, rhs, H, residual, damping)
        end
    end
    diag.regularization = λ
    return rhs
end

end
