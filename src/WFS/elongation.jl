@inline function ensure_kernel(kernel::Vector{T}, needed::Int) where {T}
    if length(kernel) < needed
        resize!(kernel, needed)
    end
    return kernel
end

@inline function ensure_kernel(kernel::AbstractVector{T}, needed::Int) where {T}
    if length(kernel) < needed
        kernel = similar(kernel, needed)
    end
    return kernel
end

function apply_elongation!(intensity::AbstractMatrix{T}, factor::Real, tmp::AbstractMatrix{T},
    kernel::AbstractVector{T}) where {T<:AbstractFloat}
    if factor <= 1
        return kernel
    end
    sigma = T(0.5) * (T(factor) - one(T))
    if sigma <= 0
        return kernel
    end
    half = max(1, ceil(Int, 2 * sigma))
    needed = 2 * half + 1
    kernel = ensure_kernel(kernel, needed)
    @inbounds for k in -half:half
        kernel[k + half + 1] = exp(-T(0.5) * (T(k) / sigma)^2)
    end
    view_kernel = @view kernel[1:needed]
    view_kernel ./= sum(view_kernel)

    n1, n2 = size(intensity)
    @inbounds for i in 1:n1, j in 1:n2
        acc = zero(T)
        for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += intensity[i, jj] * kernel[k + half + 1]
        end
        tmp[i, j] = acc
    end
    copyto!(intensity, tmp)
    return kernel
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Int, wavelength::Real)
    return 206265 * float(wavelength) / float(diameter) / diffraction_padding
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Real, wavelength::Real)
    return 206265 * float(wavelength) / float(diameter) / float(diffraction_padding)
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Int, src::LGSSource)
    return lgs_pixel_scale(diameter, diffraction_padding, wavelength(src))
end

function lgs_reference_vector(tel::Telescope, x0::Real, y0::Real, altitude::Real)
    x = (tel.params.diameter / 4) * (x0 / altitude)
    y = (tel.params.diameter / 4) * (y0 / altitude)
    z = (tel.params.diameter^2) / (8 * altitude) / (2 * sqrt(3))
    return (x, y, z)
end

function lgs_spot_shift(vec::Tuple{T,T,T}, tel::Telescope, x_subap::Real, y_subap::Real,
    pixel_scale::Real) where {T<:Real}
    dx0 = vec[1] * (4 / tel.params.diameter)
    dy0 = vec[2] * (4 / tel.params.diameter)
    dx1 = vec[3] * (sqrt(3) * (4 / tel.params.diameter)^2) * x_subap
    dy1 = vec[3] * (sqrt(3) * (4 / tel.params.diameter)^2) * y_subap
    shift_x = 206265 * (dx0 + dx1) / pixel_scale
    shift_y = 206265 * (dy0 + dy1) / pixel_scale
    return shift_x, shift_y
end

function gaussian_shift!(kernel::AbstractMatrix{T}, sigma::Real, cx::Real, cy::Real,
    weight::Real) where {T<:AbstractFloat}
    if weight == 0
        return kernel
    end
    n = size(kernel, 1)
    if sigma <= 0
        ix = clamp(round(Int, cx), 1, n)
        iy = clamp(round(Int, cy), 1, n)
        kernel[ix, iy] += T(weight)
        return kernel
    end
    inv_2sigma2 = T(0.5) / (T(sigma) * T(sigma))
    @inbounds for i in 1:n, j in 1:n
        dx = T(i) - T(cx)
        dy = T(j) - T(cy)
        kernel[i, j] += T(weight) * exp(-(dx * dx + dy * dy) * inv_2sigma2)
    end
    return kernel
end

function lgs_spot_kernel!(kernel::AbstractMatrix{T}, tel::Telescope, src::LGSSource,
    altitudes::AbstractVector, weights::AbstractVector, pixel_scale::Real, sigma_px::Real,
    center::Real, x_subap::Real, y_subap::Real,
    ref_vec::Tuple{<:Real,<:Real,<:Real}) where {T<:AbstractFloat}
    fill!(kernel, zero(T))
    x0 = src.params.laser_coordinates[2]
    y0 = -src.params.laser_coordinates[1]
    @inbounds for k in eachindex(altitudes)
        vec = lgs_reference_vector(tel, x0, y0, altitudes[k])
        vec = (vec[1] - ref_vec[1], vec[2] - ref_vec[2], vec[3] - ref_vec[3])
        shift_x, shift_y = lgs_spot_shift(vec, tel, x_subap, y_subap, pixel_scale)
        gaussian_shift!(kernel, sigma_px, center + shift_x, center + shift_y, weights[k])
    end
    total = sum(kernel)
    if total > 0
        kernel ./= total
    end
    return kernel
end

function apply_lgs_convolution!(intensity::AbstractMatrix{T}, kernel_fft::AbstractMatrix{Complex{T}},
    fft_buffer::AbstractMatrix{Complex{T}}, fft_plan, ifft_plan) where {T<:AbstractFloat}
    if size(kernel_fft, 1) == 0
        return intensity
    end
    n = size(intensity, 1)
    @inbounds for i in 1:n, j in 1:n
        fft_buffer[i, j] = complex(intensity[i, j], zero(T))
    end
    mul!(fft_buffer, fft_plan, fft_buffer)
    @. fft_buffer *= kernel_fft
    mul!(fft_buffer, ifft_plan, fft_buffer)
    scale = T(1) / (n * n)
    @inbounds for i in 1:n, j in 1:n
        intensity[i, j] = real(fft_buffer[i, j]) * scale
    end
    return intensity
end

function apply_lgs_convolution!(intensity::AbstractMatrix{T}, kernel_fft::AbstractMatrix{Complex{T}},
    fft_buffer::AbstractMatrix{Complex{T}}, fft_plan,
    ifft_buffer::AbstractMatrix{Complex{T}}, ifft_plan) where {T<:AbstractFloat}
    if size(kernel_fft, 1) == 0
        return intensity
    end
    n = size(intensity, 1)
    @inbounds for i in 1:n, j in 1:n
        fft_buffer[i, j] = complex(intensity[i, j], zero(T))
    end
    mul!(fft_buffer, fft_plan, fft_buffer)
    @. fft_buffer *= kernel_fft
    copyto!(ifft_buffer, fft_buffer)
    mul!(ifft_buffer, ifft_plan, ifft_buffer)
    scale = T(1) / (n * n)
    @inbounds for i in 1:n, j in 1:n
        intensity[i, j] = real(ifft_buffer[i, j]) * scale
    end
    return intensity
end

function lgs_average_kernel_fft(tel::Telescope, src::LGSSource, pad::Int, n_subap::Int,
    pixel_scale::Real, fft_buffer::AbstractMatrix{Complex{T}}, fft_plan) where {T<:AbstractFloat}
    na_profile = src.params.na_profile
    if na_profile === nothing
        return similar(fft_buffer, Complex{T}, 0, 0)
    end
    altitudes = na_profile[1, :]
    weights = na_profile[2, :]
    if length(altitudes) == 0
        return similar(fft_buffer, Complex{T}, 0, 0)
    end

    fwhm_px = src.params.fwhm_spot_up / pixel_scale
    sigma_px = fwhm_px / (2 * sqrt(2 * log(T(2))))
    center = (pad + 1) / 2

    x_subap = range(-tel.params.diameter / 2, tel.params.diameter / 2; length=n_subap)
    y_subap = range(-tel.params.diameter / 2, tel.params.diameter / 2; length=n_subap)
    kernel = similar(fft_buffer, T, pad, pad)
    temp = similar(kernel)
    fill!(kernel, zero(T))

    x0 = src.params.laser_coordinates[2]
    y0 = -src.params.laser_coordinates[1]
    ref_idx = Int(cld(length(altitudes), 2))
    ref_vec = lgs_reference_vector(tel, x0, y0, altitudes[ref_idx])

    @inbounds for iy in 1:n_subap, ix in 1:n_subap
        lgs_spot_kernel!(temp, tel, src, altitudes, weights, pixel_scale, sigma_px, center,
            x_subap[ix], y_subap[iy], ref_vec)
        kernel .+= temp
    end
    kernel ./= max(n_subap * n_subap, 1)

    @inbounds for i in 1:pad, j in 1:pad
        fft_buffer[i, j] = complex(kernel[i, j], zero(T))
    end
    mul!(fft_buffer, fft_plan, fft_buffer)
    kernel_fft = similar(fft_buffer, Complex{T}, pad, pad)
    copyto!(kernel_fft, fft_buffer)
    return kernel_fft
end
