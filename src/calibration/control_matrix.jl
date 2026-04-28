using LinearAlgebra

#
# Control matrix storage
#
# `ControlMatrix` packages a measured forward operator together with its
# chosen inverse representation.
#
# The typical flow is:
# 1. measure or assemble the direct interaction matrix
# 2. choose an inverse policy (pseudo-inverse, TSVD, damping, ...)
# 3. compute the inverse operator and retain inversion diagnostics
#
# Runtime reconstructors use the stored control matrix while the object keeps
# enough metadata to report condition number, effective rank, and any applied
# truncation.
#
struct ControlMatrix{T<:AbstractFloat,
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

@inline forward_operator(control_matrix::ControlMatrix) = control_matrix.D
@inline inverse_operator_matrix(control_matrix::ControlMatrix) = control_matrix.M
@inline inverse_policy(control_matrix::ControlMatrix) = control_matrix.policy
@inline singular_values(control_matrix::ControlMatrix) = control_matrix.singular_values
@inline condition_number(control_matrix::ControlMatrix) = control_matrix.cond
@inline effective_rank(control_matrix::ControlMatrix) = control_matrix.effective_rank
@inline truncation_count(control_matrix::ControlMatrix) = control_matrix.n_trunc

"""
    ControlMatrix(D; n_trunc=0, invert=true, policy=..., build_backend=...)

Construct the stored inverse representation for a calibration operator.

If `invert=true`, the constructor computes the control matrix selected by
`policy` and retains diagnostics such as singular values, condition number,
effective rank, and any explicit truncation count.
"""
function ControlMatrix(D::AbstractMatrix{T}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(T),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D)) where {T<:AbstractFloat}
    if !invert
        empty_vals = similar(D, T, 0)
        return ControlMatrix{T, typeof(D), Matrix{T}, typeof(empty_vals), typeof(policy), typeof(build_backend)}(
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
    effective_policy = calibration_effective_policy(policy, n_trunc)
    Minv, stats = inverse_operator(build_backend, D, effective_policy)
    Minv_native = materialize_runtime_build_result(build_backend, D, Minv)
    Mtrunc = calibration_truncated_matrix(effective_policy, Minv_native)
    singular_values = materialize_runtime_build_result(build_backend, similar(D, T, 0), stats.singular_values)
    return ControlMatrix{T, typeof(D), typeof(Minv_native), typeof(singular_values), typeof(effective_policy), typeof(build_backend)}(
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

calibration_effective_policy(policy::InversePolicy, n_trunc::Int) = policy
calibration_effective_policy(policy::TSVDInverse, n_trunc::Int) = TSVDInverse(
    rtol=policy.rtol,
    atol=policy.atol,
    n_trunc=n_trunc,
)

materialize_runtime_build_result(::CPUBuildBackend, ref, data) =
    materialize_build(NativeBuildBackend(), ref, data)
materialize_runtime_build_result(::BuildBackend, ref, data) = data

calibration_truncated_matrix(::InversePolicy, Minv) = nothing
calibration_truncated_matrix(::TSVDInverse, Minv) = Minv

function ControlMatrix(D::AbstractMatrix{S}; n_trunc::Int=0, invert::Bool=true,
    policy::InversePolicy=default_calibration_inverse_policy(float(S)),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D)) where {S<:Real}
    return ControlMatrix(float.(D); n_trunc=n_trunc, invert=invert, policy=policy, build_backend=build_backend)
end

"""
    with_truncation(control_matrix, n_trunc)

Rebuild a control matrix with the same forward operator and inverse policy
but a different explicit TSVD truncation count.
"""
function with_truncation(control_matrix::ControlMatrix, n_trunc::Int)
    return ControlMatrix(control_matrix.D; n_trunc=n_trunc, invert=true, policy=control_matrix.policy,
        build_backend=control_matrix.build_backend)
end
