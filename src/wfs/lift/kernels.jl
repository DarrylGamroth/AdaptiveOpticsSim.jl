#
# LiFT physical forward model
#
# Focal-plane propagation never owns or triggers a detector; acquisition
# timing, QE, and stochastic readout remain at the detector boundary.
# AdaptiveOpticsCalibration owns the iterative inverse estimator.
#
# Forward model:
# 1. combine modal coefficients into an OPD map
# 2. add the configured diversity term
# 3. propagate to the focal plane and form intensity
# 4. optionally convolve with an object kernel
#
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
