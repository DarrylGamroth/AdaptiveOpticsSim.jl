using LinearAlgebra

struct CalibrationVault{T<:AbstractFloat,
    D<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T},
    P<:InversePolicy}
    D::D
    M::Union{Nothing,M}
    Mtrunc::Union{Nothing,M}
    singular_values::V
    cond::T
    effective_rank::Int
    n_trunc::Int
    policy::P
end

function CalibrationVault(D::AbstractMatrix{T}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(T)) where {T<:AbstractFloat}
    if !invert
        empty_vals = similar(D, T, 0)
        return CalibrationVault{T, typeof(D), Matrix{T}, typeof(empty_vals), typeof(policy)}(
            D,
            nothing,
            nothing,
            empty_vals,
            T(NaN),
            0,
            0,
            policy,
        )
    end
    if n_trunc < 0
        throw(InvalidConfiguration("n_trunc must be >= 0"))
    end
    effective_policy = policy isa TSVDInverse ? TSVDInverse(
        rtol=policy.rtol,
        atol=policy.atol,
        n_trunc=n_trunc,
    ) : policy
    Minv, stats = inverse_operator(D, effective_policy)
    Mtrunc = effective_policy isa TSVDInverse ? Minv : nothing
    return CalibrationVault{T, typeof(D), typeof(Minv), typeof(stats.singular_values), typeof(effective_policy)}(
        D,
        Minv,
        Mtrunc,
        stats.singular_values,
        stats.cond,
        stats.effective_rank,
        stats.n_trunc,
        effective_policy,
    )
end

function CalibrationVault(D::AbstractMatrix{S}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(float(S))) where {S<:Real}
    return CalibrationVault(float.(D); n_trunc=n_trunc, invert=invert, policy=policy)
end

function with_truncation(vault::CalibrationVault, n_trunc::Int)
    return CalibrationVault(vault.D; n_trunc=n_trunc, invert=true, policy=vault.policy)
end
