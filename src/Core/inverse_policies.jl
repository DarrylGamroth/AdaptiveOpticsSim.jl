using LinearAlgebra

abstract type InversePolicy end

struct ExactPseudoInverse <: InversePolicy end

struct TSVDInverse{T<:AbstractFloat} <: InversePolicy
    rtol::T
    atol::T
    n_trunc::Int
end

struct TikhonovInverse{T<:AbstractFloat} <: InversePolicy
    lambda::T
    rtol::T
    atol::T
end

TSVDInverse(; rtol::Real=eps(Float64), atol::Real=0.0, n_trunc::Integer=0) =
    TSVDInverse(float(rtol), float(atol), Int(n_trunc))

TikhonovInverse(lambda::Real; rtol::Real=eps(Float64), atol::Real=0.0) =
    TikhonovInverse(float(lambda), float(rtol), float(atol))

struct InverseStats{T<:AbstractFloat,V<:AbstractVector{T}}
    singular_values::V
    cond::T
    effective_rank::Int
    n_trunc::Int
end

function _inverse_cutoff(s::AbstractVector{T}, rtol::T, atol::T) where {T<:AbstractFloat}
    isempty(s) && return zero(T)
    return max(atol, rtol * s[begin])
end

function inverse_operator(A::AbstractMatrix{T}, ::ExactPseudoInverse) where {T<:AbstractFloat}
    F = svd(A; full=false)
    s = F.S
    inv_s = similar(s)
    @inbounds for i in eachindex(s)
        inv_s[i] = iszero(s[i]) ? zero(T) : inv(s[i])
    end
    M = Matrix(F.V * Diagonal(inv_s) * F.U')
    effective_rank = count(!iszero, s)
    cond = effective_rank == 0 ? T(Inf) : s[begin] / s[effective_rank]
    return M, InverseStats(s, cond, effective_rank, 0)
end

function inverse_operator(A::AbstractMatrix{T}, policy::TSVDInverse) where {T<:AbstractFloat}
    policy.n_trunc >= 0 || throw(InvalidConfiguration("TSVD n_trunc must be >= 0"))
    F = svd(A; full=false)
    s = F.S
    isempty(s) && return Matrix{T}(undef, size(A, 2), size(A, 1)), InverseStats(s, T(Inf), 0, 0)
    cutoff = _inverse_cutoff(s, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s = similar(s)
    fill!(inv_s, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s[i] = inv(s[i])
    end
    M = Matrix(F.V * Diagonal(inv_s) * F.U')
    cond = effective_rank == 0 ? T(Inf) : s[begin] / s[effective_rank]
    return M, InverseStats(s, cond, effective_rank, length(s) - effective_rank)
end

function inverse_operator(A::AbstractMatrix{T}, policy::TikhonovInverse) where {T<:AbstractFloat}
    policy.lambda >= 0 || throw(InvalidConfiguration("Tikhonov lambda must be >= 0"))
    F = svd(A; full=false)
    s = F.S
    inv_s = similar(s)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s)
        inv_s[i] = s[i] / (s[i]^2 + λ2)
    end
    cutoff = _inverse_cutoff(s, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s)
    denom = max(isempty(s) ? zero(T) : s[end], T(policy.lambda))
    cond = (isempty(s) || iszero(denom)) ? T(Inf) : s[begin] / denom
    M = Matrix(F.V * Diagonal(inv_s) * F.U')
    return M, InverseStats(s, cond, effective_rank, 0)
end
