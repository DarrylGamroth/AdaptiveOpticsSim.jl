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

    # Match pyTomoAO/NumPy fitting behavior with a pseudoinverse rather than
    # a regularized normal-equations solve.
    pinv_rtol = regularization == sqrt(eps(T)) ? T(1e-15) : T(regularization)
    fitting_matrix = pinv(Matrix(influence_functions); rtol=pinv_rtol)
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

function _extract_actuator_coordinates(valid_actuators::AbstractMatrix{Bool})
    return findall(valid_actuators)
end

function _map_actuators_to_grid(
    actuator_coords::AbstractVector{<:CartesianIndex{2}},
    original_shape::Tuple{Int,Int},
    new_shape::Tuple{Int,Int};
    stretch_factor::T,
) where {T<:AbstractFloat}
    orig_height, orig_width = original_shape
    new_height, new_width = new_shape
    mapped = Vector{NTuple{2,T}}(undef, length(actuator_coords))
    @inbounds for (k, idx) in pairs(actuator_coords)
        y = idx[1] - 1
        x = idx[2] - 1
        normalized_y = T(y) / T(orig_height - 1)
        normalized_x = T(x) / T(orig_width - 1)
        mapped_y = (T(2) * normalized_y - one(T)) * stretch_factor
        mapped_x = (T(2) * normalized_x - one(T)) * stretch_factor
        mapped_y = mapped_y * T(new_height - 1) / T(2) + T(new_height - 1) / T(2)
        mapped_x = mapped_x * T(new_width - 1) / T(2) + T(new_width - 1) / T(2)
        mapped[k] = (mapped_y, mapped_x)
    end
    return mapped
end

function influence_functions(
    dm::TomographyDMParams{T};
    resolution::Integer=size(dm_valid_support(dm), 1),
    w1::Real=2,
    w2::Real=-1,
    sigma1::Real=1.0,
    sigma2::Real=1.7,
    stretch_factor::Real=1.03,
) where {T<:AbstractFloat}
    resolution > 0 || throw(InvalidConfiguration("resolution must be positive"))
    support = dm_valid_support(dm)
    actuator_coords = _extract_actuator_coordinates(support)
    mapped = _map_actuators_to_grid(
        actuator_coords,
        size(support),
        (Int(resolution), Int(resolution));
        stretch_factor=T(stretch_factor),
    )

    n_modes = count(dm.valid_actuators)
    modes = Matrix{T}(undef, Int(resolution)^2, n_modes)
    two_sigma1_sq = T(2) * T(sigma1)^2
    two_sigma2_sq = T(2) * T(sigma2)^2
    norm1 = T(w1) / (T(2π) * T(sigma1)^2)
    norm2 = T(w2) / (T(2π) * T(sigma2)^2)

    @inbounds for (mode_idx, (cy, cx)) in pairs(mapped)
        col = 1
        for y in 0:Int(resolution)-1
            y_t = T(y)
            for x in 0:Int(resolution)-1
                x_t = T(x)
                r2 = (x_t - cx)^2 + (y_t - cy)^2
                modes[col, mode_idx] = norm1 * exp(-r2 / two_sigma1_sq) +
                    norm2 * exp(-r2 / two_sigma2_sq)
                col += 1
            end
        end
    end
    return modes
end

function TomographyFitting(
    dm::TomographyDMParams{T};
    regularization::Real=sqrt(eps(T)),
    resolution::Integer=size(dm_valid_support(dm), 1),
    w1::Real=2,
    w2::Real=-1,
    sigma1::Real=1.0,
    sigma2::Real=1.7,
    stretch_factor::Real=1.03,
) where {T<:AbstractFloat}
    modes = influence_functions(
        dm;
        resolution=resolution,
        w1=w1,
        w2=w2,
        sigma1=sigma1,
        sigma2=sigma2,
        stretch_factor=stretch_factor,
    )
    return TomographyFitting(modes; regularization=regularization, resolution=resolution)
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
