module AdaptiveOpticsSimAMDGPUExt

using AdaptiveOpticsSim
using AMDGPU
using AbstractFFTs
using LinearAlgebra

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

function _amdgpu_dense_copy(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    out = AMDGPU.ROCArray{T}(undef, size(A)...)
    copyto!(out, A)
    return out
end

function _amdgpu_svd(A::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    F = copy(A)
    U, S, Vt = AMDGPU.rocSOLVER.gesvd!('S', 'S', F)
    return U, S, Vt, AdaptiveOpticsSim.singular_values_host(S)
end

function _amdgpu_lu_solve!(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    n = size(A, 1)
    ipiv = AMDGPU.ROCArray{Cint}(undef, n)
    AMDGPU.rocSOLVER.getrf!(A, ipiv)
    AMDGPU.rocSOLVER.getrs!('N', A, ipiv, B)
    return B
end

function _amdgpu_cholesky_rhs(rhs::AMDGPU.ROCArray{T,1}) where {T<:AbstractFloat}
    rhs_mat = AMDGPU.ROCArray{T}(undef, length(rhs), 1)
    copyto!(vec(rhs_mat), rhs)
    return rhs_mat
end

_amdgpu_cholesky_rhs(rhs::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} = rhs

function _amdgpu_finalize_rhs!(dest::AMDGPU.ROCArray{T,1}, rhs_mat::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    copyto!(dest, vec(rhs_mat))
    return dest
end

_amdgpu_finalize_rhs!(::AMDGPU.ROCArray{T,2}, rhs_mat::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} = rhs_mat

function _amdgpu_cholesky_solve!(
    A::AMDGPU.ROCArray{T,2},
    rhs::AMDGPU.ROCArray{T,N};
    check::Bool=false,
) where {T<:AbstractFloat,N}
    rhs_mat = _amdgpu_cholesky_rhs(rhs)
    chol = cholesky!(Hermitian(A), check=check)
    if issuccess(chol)
        ldiv!(chol, rhs_mat)
    end
    return _amdgpu_finalize_rhs!(rhs, rhs_mat), chol
end

function _amdgpu_pseudoinverse(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    U::AMDGPU.ROCArray{T,2}, S::AMDGPU.ROCArray{T,1}, Vt::AMDGPU.ROCArray{T,2},
    inv_s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s = AdaptiveOpticsSim.materialize_build(backend, S, inv_s_host)
    U_scaled = similar(U)
    copyto!(U_scaled, U)
    U_scaled .*= reshape(inv_s, 1, :)
    M = adjoint(Vt) * adjoint(U_scaled)
    return M
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, ::AdaptiveOpticsSim.ExactPseudoInverse) where {T<:AbstractFloat}
    U, S, Vt, s_host = _amdgpu_svd(A)
    inv_s_host = similar(s_host)
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = iszero(s_host[i]) ? zero(T) : inv(s_host[i])
    end
    M = _amdgpu_pseudoinverse(backend, U, S, Vt, inv_s_host)
    effective_rank = count(!iszero, s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return M, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TSVDInverse) where {T<:AbstractFloat}
    U, S, Vt, s_host = _amdgpu_svd(A)
    isempty(s_host) && return AdaptiveOpticsSim.materialize_build(backend, similar(A, T, size(A, 2), size(A, 1))),
        AdaptiveOpticsSim.InverseStats(s_host, T(Inf), 0, 0)
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s_host)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s_host = similar(s_host)
    fill!(inv_s_host, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s_host[i] = inv(s_host[i])
    end
    M = _amdgpu_pseudoinverse(backend, U, S, Vt, inv_s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return M, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, length(s_host) - effective_rank)
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TikhonovInverse) where {T<:AbstractFloat}
    U, S, Vt, s_host = _amdgpu_svd(A)
    inv_s_host = similar(s_host)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = s_host[i] / (s_host[i]^2 + λ2)
    end
    M = _amdgpu_pseudoinverse(backend, U, S, Vt, inv_s_host)
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s_host)
    denom = max(isempty(s_host) ? zero(T) : s_host[end], T(policy.lambda))
    cond = (isempty(s_host) || iszero(denom)) ? T(Inf) : s_host[begin] / denom
    return M, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function AdaptiveOpticsSim.stable_hermitian_right_division(
    backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    rhs::AMDGPU.ROCArray{T,2},
    css::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    css_mat = copy(css)
    rhs_t = permutedims(rhs, (2, 1))
    _, fact = _amdgpu_cholesky_solve!(css_mat, rhs_t; check=false)
    if issuccess(fact)
        return permutedims(rhs_t, (2, 1))
    end
    css_lu = copy(css)
    _amdgpu_lu_solve!(css_lu, rhs_t)
    return permutedims(rhs_t, (2, 1))
end

function AdaptiveOpticsSim.solve_lift_fallback!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T},
    rhs::AMDGPU.ROCArray{T,1}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::AdaptiveOpticsSim.LiFTDampingMode) where {T<:AbstractFloat}
    H_mat = _amdgpu_dense_copy(H)
    U, S, Vt, _ = _amdgpu_svd(H_mat)
    λ = damping isa AdaptiveOpticsSim.LiFTLevenbergMarquardt ?
        AdaptiveOpticsSim.damping_lambda(damping, adjoint(H_mat) * H_mat) : zero(T)
    work = AMDGPU.ROCArray{T}(undef, length(S))
    mul!(work, transpose(U), residual)
    @. work = ifelse(iszero(S^2 + λ^2), zero(T), (S * work) / (S^2 + λ^2))
    mul!(rhs, adjoint(Vt), work)
    return rhs
end

function AdaptiveOpticsSim.solve_normal_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    ::AdaptiveOpticsSim.LiFTDampingNone) where {T<:AbstractFloat}
    factor_mat = _amdgpu_dense_copy(normal)
    _, chol = _amdgpu_cholesky_solve!(factor_mat, rhs; check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = AdaptiveOpticsSim.regularization_load(normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = _amdgpu_cholesky_solve!(factor_mat, rhs; check=false)
        if !issuccess(chol)
            λ *= T(10)
            @views factor_mat[diagind(factor_mat)] .+= λ
            _, chol = _amdgpu_cholesky_solve!(factor_mat, rhs; check=false)
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
    factor_mat = _amdgpu_dense_copy(normal)
    λ = AdaptiveOpticsSim.damping_lambda(damping, normal)
    if λ > zero(T)
        @views factor_mat[diagind(factor_mat)] .+= λ
    end
    _, chol = _amdgpu_cholesky_solve!(factor_mat, rhs; check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), AdaptiveOpticsSim.regularization_load(normal))
        copyto!(factor_mat, normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = _amdgpu_cholesky_solve!(factor_mat, rhs; check=false)
        if λ > T(1e12)
            return AdaptiveOpticsSim.solve_lift_fallback!(diag, rhs, H, residual, damping)
        end
    end
    diag.regularization = λ
    return rhs
end

end
