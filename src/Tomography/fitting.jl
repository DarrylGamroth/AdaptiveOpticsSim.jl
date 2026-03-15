using LinearAlgebra

struct TomographyFitting{T<:AbstractFloat, M<:AbstractMatrix{T}}
    influence_functions::M
    fitting_matrix::M
    resolution::Int
end

function TomographyFitting(
    influence_functions::AbstractMatrix{T};
    regularization::Real=sqrt(eps(T)),
    resolution::Integer=round(Int, sqrt(size(influence_functions, 1))),
) where {T<:AbstractFloat}
    size(influence_functions, 1) > 0 ||
        throw(InvalidConfiguration("influence_functions must have at least one row"))
    size(influence_functions, 2) > 0 ||
        throw(InvalidConfiguration("influence_functions must have at least one column"))
    regularization >= 0 || throw(InvalidConfiguration("regularization must be non-negative"))

    gram = Matrix{T}(undef, size(influence_functions, 2), size(influence_functions, 2))
    mul!(gram, adjoint(influence_functions), influence_functions)
    @inbounds for k in axes(gram, 1)
        gram[k, k] += T(regularization)
    end
    rhs = Matrix(adjoint(influence_functions))
    fitting_matrix = gram \ rhs
    return TomographyFitting{T, typeof(fitting_matrix)}(
        Matrix(influence_functions),
        fitting_matrix,
        Int(resolution),
    )
end

function TomographyFitting(dm::DeformableMirror; regularization::Real=sqrt(eps(eltype(dm.state.modes))))
    resolution = round(Int, sqrt(size(dm.state.modes, 1)))
    return TomographyFitting(dm.state.modes; regularization=regularization, resolution=resolution)
end

function fit_commands!(
    out::AbstractVector{T},
    fitting::TomographyFitting{T},
    opd::AbstractVector{T},
) where {T<:AbstractFloat}
    size(fitting.fitting_matrix, 1) == length(out) ||
        throw(DimensionMismatchError("output length must match fitting matrix row count"))
    size(fitting.fitting_matrix, 2) == length(opd) ||
        throw(DimensionMismatchError("opd length must match fitting matrix column count"))
    mul!(out, fitting.fitting_matrix, opd)
    return out
end

function fit_commands(fitting::TomographyFitting{T}, opd::AbstractVector{T}) where {T<:AbstractFloat}
    out = Vector{T}(undef, size(fitting.fitting_matrix, 1))
    return fit_commands!(out, fitting, opd)
end
