using LinearAlgebra

struct CalibrationVault{T<:AbstractFloat,
    D<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T}}
    D::D
    M::Union{Nothing,M}
    Mtrunc::Union{Nothing,M}
    singular_values::V
    cond::T
    n_trunc::Int
end

function CalibrationVault(D::AbstractMatrix{T}; n_trunc::Int=0, invert::Bool=true) where {T<:AbstractFloat}
    if !invert
        empty_vals = similar(D, T, 0)
        return CalibrationVault{T, typeof(D), Matrix{T}, typeof(empty_vals)}(
            D,
            nothing,
            nothing,
            empty_vals,
            T(NaN),
            0,
        )
    end
    if n_trunc < 0
        throw(InvalidConfiguration("n_trunc must be >= 0"))
    end

    F = svd(D; full=false)
    s = F.S
    if n_trunc >= length(s)
        throw(InvalidConfiguration("n_trunc must be < $(length(s))"))
    end

    inv_s = similar(s)
    @inbounds for i in eachindex(s)
        inv_s[i] = s[i] == 0 ? zero(T) : one(T) / s[i]
    end

    Minv = Matrix(F.V * Diagonal(inv_s) * F.U')
    n_keep = length(s) - n_trunc
    Vtrunc = F.V[:, 1:n_keep]
    Utrunc = F.U[:, 1:n_keep]
    inv_s_trunc = inv_s[1:n_keep]
    Mtrunc = Matrix(Vtrunc * Diagonal(inv_s_trunc) * Utrunc')
    cond = s[1] / s[end - n_trunc]
    return CalibrationVault{T, typeof(D), typeof(Minv), typeof(s)}(D, Minv, Mtrunc, s, cond, n_trunc)
end

function CalibrationVault(D::AbstractMatrix{S}; n_trunc::Int=0, invert::Bool=true) where {S<:Real}
    return CalibrationVault(float.(D); n_trunc=n_trunc, invert=invert)
end

function with_truncation(vault::CalibrationVault, n_trunc::Int)
    return CalibrationVault(vault.D; n_trunc=n_trunc, invert=true)
end
