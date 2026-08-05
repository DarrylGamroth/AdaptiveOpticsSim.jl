#
# LiFT phase retrieval
#
# LiFT fits modal coefficients by matching a separately prepared focal-plane
# forward model to a caller-owned observation. Focal-plane propagation never owns or
# triggers a detector; acquisition timing, QE, and stochastic readout remain at
# the detector boundary.
#
# Forward model:
# 1. combine modal coefficients into an OPD map
# 2. add the configured diversity term
# 3. propagate to the focal plane and form intensity
# 4. optionally convolve with an object kernel
#
# Inverse model:
# - `LiFTAnalyticJacobian` builds the Jacobian from focal-plane field derivatives
# - `LiFTNumericalJacobian` builds the Jacobian by centered finite differences
#
# The update step is then solved by QR or normal equations with explicit
# damping/fallback logic.
#
@kernel function lift_scatter_update_kernel!(coeffs, delta, mode_ids, n_modes::Int)
    i = @index(Global, Linear)
    if i <= n_modes
        @inbounds coeffs[mode_ids[i]] += delta[i]
    end
end
@kernel function lift_gather_kernel!(out, coeffs, mode_ids, n_modes::Int)
    i = @index(Global, Linear)
    if i <= n_modes
        @inbounds out[i] = coeffs[mode_ids[i]]
    end
end

@kernel function lift_sqrt_weights_kernel!(weights, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds weights[i] = sqrt(weights[i])
    end
end

@kernel function lift_affine_basis_mode_kernel!(dest, base, basis,
    scale, mode_offset::Int, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = base[i] + scale * basis[i + mode_offset]
    end
end

@kernel function lift_scaled_basis_mode_kernel!(dest, amplitude, basis,
    scale, mode_offset::Int, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = amplitude[i] * scale * basis[i + mode_offset]
    end
end

@kernel function lift_copy_column_kernel!(dest, column::Int, src, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i, column] = src[i]
    end
end

@kernel function lift_residual_kernel!(dest, observation, model, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = observation[i] - model[i]
    end
end

@kernel function lift_row_weights_kernel!(matrix, weights,
    n_rows::Int, n_cols::Int)
    i, j = @index(Global, NTuple)
    if i <= n_rows && j <= n_cols
        @inbounds matrix[i, j] *= weights[i]
    end
end

@kernel function lift_inverse_variance_kernel!(dest, values,
    scale, offset, floor_value, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = inv(max(values[i] * scale + offset,
            floor_value))
    end
end

@kernel function lift_add_diagonal_kernel!(matrix, value, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds matrix[i, i] += value
    end
end

@kernel function lift_dense_convolution_kernel!(dest, src, kernel,
    inv_norm, n::Int, m::Int, kh::Int, kw::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        cx = div(kh, 2)
        cy = div(kw, 2)
        @inbounds for ki in 1:kh, kj in 1:kw
            ii = symm_index(i + ki - cx - 1, n)
            jj = symm_index(j + kj - cy - 1, m)
            acc += src[ii, jj] * kernel[ki, kj]
        end
        @inbounds dest[i, j] = acc * inv_norm
    end
end

@kernel function lift_row_convolution_kernel!(dest, src, kernel,
    n::Int, m::Int, nk::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        center = div(nk, 2)
        @inbounds for k in 1:nk
            ii = symm_index(i + k - center - 1, n)
            acc += src[ii, j] * kernel[k]
        end
        @inbounds dest[i, j] = acc
    end
end

@kernel function lift_column_convolution_kernel!(dest, src, kernel,
    inv_norm, n::Int, m::Int, nk::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        center = div(nk, 2)
        @inbounds for k in 1:nk
            jj = symm_index(j + k - center - 1, m)
            acc += src[i, jj] * kernel[k]
        end
        @inbounds dest[i, j] = acc * inv_norm
    end
end
