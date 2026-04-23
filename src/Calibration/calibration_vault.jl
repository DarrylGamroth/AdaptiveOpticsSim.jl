using LinearAlgebra

#
# Calibration inverse storage
#
# `CalibrationVault` packages a measured forward operator together with its
# chosen inverse representation.
#
# The typical flow is:
# 1. measure or assemble the direct interaction matrix `D`
# 2. choose an inverse policy (pseudo-inverse, TSVD, damping, ...)
# 3. compute the inverse operator and retain inversion diagnostics
#
# Runtime reconstructors then use the stored inverse matrix while the vault
# keeps enough metadata to report condition number, effective rank, and any
# applied truncation.
#
struct CalibrationVault{T<:AbstractFloat,
    D<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T},
    P<:InversePolicy,
    B<:BuildBackend}
    D::D
    M::Union{Nothing,M}
    Mtrunc::Union{Nothing,M}
    singular_values::V
    cond::T
    effective_rank::Int
    n_trunc::Int
    policy::P
    build_backend::B
end

@inline forward_operator(vault::CalibrationVault) = vault.D
@inline inverse_operator_matrix(vault::CalibrationVault) = vault.M
@inline inverse_policy(vault::CalibrationVault) = vault.policy
@inline singular_values(vault::CalibrationVault) = vault.singular_values
@inline condition_number(vault::CalibrationVault) = vault.cond
@inline effective_rank(vault::CalibrationVault) = vault.effective_rank
@inline truncation_count(vault::CalibrationVault) = vault.n_trunc

"""
    CalibrationVault(D; n_trunc=0, invert=true, policy=..., build_backend=...)

Construct the stored inverse representation for a calibration operator.

If `invert=true`, the vault computes the inverse operator selected by `policy`
and retains diagnostics such as singular values, condition number, effective
rank, and any explicit truncation count.
"""
function CalibrationVault(D::AbstractMatrix{T}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(T),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D)) where {T<:AbstractFloat}
    if !invert
        empty_vals = similar(D, T, 0)
        return CalibrationVault{T, typeof(D), Matrix{T}, typeof(empty_vals), typeof(policy), typeof(build_backend)}(
            D,
            nothing,
            nothing,
            empty_vals,
            T(NaN),
            0,
            0,
            policy,
            build_backend,
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
    Minv, stats = inverse_operator(build_backend, D, effective_policy)
    Minv_native = build_backend isa CPUBuildBackend ?
        materialize_build(NativeBuildBackend(), D, Minv) : Minv
    Mtrunc = effective_policy isa TSVDInverse ? Minv_native : nothing
    singular_values = build_backend isa CPUBuildBackend ?
        materialize_build(NativeBuildBackend(), similar(D, T, 0), stats.singular_values) :
        stats.singular_values
    return CalibrationVault{T, typeof(D), typeof(Minv_native), typeof(singular_values), typeof(effective_policy), typeof(build_backend)}(
        D,
        Minv_native,
        Mtrunc,
        singular_values,
        stats.cond,
        stats.effective_rank,
        stats.n_trunc,
        effective_policy,
        build_backend,
    )
end

function CalibrationVault(D::AbstractMatrix{S}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(float(S)),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D)) where {S<:Real}
    return CalibrationVault(float.(D); n_trunc=n_trunc, invert=invert, policy=policy, build_backend=build_backend)
end

"""
    with_truncation(vault, n_trunc)

Rebuild a calibration vault with the same forward operator and inverse policy
but a different explicit TSVD truncation count.
"""
function with_truncation(vault::CalibrationVault, n_trunc::Int)
    return CalibrationVault(vault.D; n_trunc=n_trunc, invert=true, policy=vault.policy,
        build_backend=vault.build_backend)
end
