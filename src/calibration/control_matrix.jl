#
# Control matrix storage
#
# `ControlMatrix` retains a measured forward operator together with a
# Calibration-authored inverse product materialized on the requested runtime
# backend. AdaptiveOpticsCalibration owns the compact-SVD equations and
# diagnostics; this AOS owner retains only plant-facing storage.
#

function _with_truncation(
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse,
    n_trunc::Nothing,
)
    return method
end

function _with_truncation(
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse,
    n_trunc::Integer,
)
    n_trunc >= 0 || throw(InvalidConfiguration("n_trunc must be >= 0"))
    return method
end

function _with_truncation(
    method::_AOC_RECONSTRUCTORS.TSVDInverse,
    n_trunc::Integer,
)
    n_trunc >= 0 || throw(InvalidConfiguration("n_trunc must be >= 0"))
    return _AOC_RECONSTRUCTORS.TSVDInverse(
        rtol=method.rtol,
        atol=method.atol,
        n_trunc=n_trunc,
    )
end

struct ControlMatrix{
    T<:AbstractFloat,
    D<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T},
    C<:_AOC_RECONSTRUCTORS.AbstractSVDInverse,
    B<:BuildBackend,
}
    D::D
    M::Union{Nothing,M}
    singular_values::V
    cond::T
    effective_rank::Int
    n_trunc::Int
    method::C
    build_backend::B
end

@inline forward_operator(control_matrix::ControlMatrix) = control_matrix.D
@inline inverse_operator_matrix(control_matrix::ControlMatrix) = control_matrix.M
@inline calibration_method(control_matrix::ControlMatrix) = control_matrix.method
@inline singular_values(control_matrix::ControlMatrix) = control_matrix.singular_values
@inline condition_number(control_matrix::ControlMatrix) = control_matrix.cond
@inline effective_rank(control_matrix::ControlMatrix) = control_matrix.effective_rank
@inline truncation_count(control_matrix::ControlMatrix) = control_matrix.n_trunc

"""
    ControlMatrix(D; n_trunc=nothing, invert=true, method=..., build_backend=...)

Materialize a Calibration-authored compact-SVD inverse for a measured plant
response. Calibration preparation is host-side and may allocate. The accepted
dense inverse and singular values are then copied to `build_backend` for later
runtime use.
"""
function ControlMatrix(
    D::AbstractMatrix{T};
    n_trunc::Union{Nothing,Integer}=nothing,
    invert::Bool=true,
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse=
        _default_svd_inverse_method(T),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D),
) where {T<:AbstractFloat}
    if !invert
        isnothing(n_trunc) || n_trunc >= 0 ||
            throw(InvalidConfiguration("n_trunc must be >= 0"))
        empty_values = similar(D, T, 0)
        return ControlMatrix{
            T,
            typeof(D),
            Matrix{T},
            typeof(empty_values),
            typeof(method),
            typeof(build_backend),
        }(
            D,
            nothing,
            empty_values,
            T(NaN),
            0,
            something(n_trunc, 0),
            method,
            build_backend,
        )
    end

    effective_method = _with_truncation(method, n_trunc)
    product = _prepare_svd_reconstructor(D, effective_method)
    inverse = materialize_runtime_build_result(
        build_backend,
        D,
        _AOC_RECONSTRUCTORS.reconstructor(product),
    )
    values = materialize_runtime_build_result(
        build_backend,
        similar(D, T, 0),
        _AOC_RECONSTRUCTORS.singular_values(product),
    )
    return ControlMatrix{
        T,
        typeof(D),
        typeof(inverse),
        typeof(values),
        typeof(effective_method),
        typeof(build_backend),
    }(
        D,
        inverse,
        values,
        _AOC_RECONSTRUCTORS.condition_number(product),
        _AOC_RECONSTRUCTORS.effective_rank(product),
        _AOC_RECONSTRUCTORS.truncation_count(product),
        effective_method,
        build_backend,
    )
end

function ControlMatrix(
    D::AbstractMatrix{S};
    n_trunc::Union{Nothing,Integer}=nothing,
    invert::Bool=true,
    method::_AOC_RECONSTRUCTORS.AbstractSVDInverse=
        _default_svd_inverse_method(float(S)),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(D),
) where {S<:Real}
    return ControlMatrix(
        float.(D);
        n_trunc=n_trunc,
        invert=invert,
        method=method,
        build_backend=build_backend,
    )
end

"""
    with_truncation(control_matrix, n_trunc)

Rebuild a control matrix from the same measured response and calibration method
with a different explicit TSVD truncation count.
"""
function with_truncation(control_matrix::ControlMatrix, n_trunc::Integer)
    return ControlMatrix(
        control_matrix.D;
        n_trunc=n_trunc,
        invert=true,
        method=control_matrix.method,
        build_backend=control_matrix.build_backend,
    )
end
