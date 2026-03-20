const ARCSEC_PER_RAD = 180 * 3600 / π

@kernel function elongation_apply_kernel!(tmp, intensity, kernel, half::Int, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n1 && j <= n2
        acc = zero(eltype(tmp))
        @inbounds for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += intensity[i, jj] * kernel[k + half + 1]
        end
        @inbounds tmp[i, j] = acc
    end
end

@kernel function elongation_apply_stack_kernel!(tmp, intensity_stack, kernel, half::Int, n1::Int, n2::Int, n3::Int)
    i, j, k = @index(Global, NTuple)
    if i <= n1 && j <= n2 && k <= n3
        acc = zero(eltype(tmp))
        @inbounds for dk in -half:half
            jj = clamp(j + dk, 1, n2)
            acc += intensity_stack[i, jj, k] * kernel[dk + half + 1]
        end
        @inbounds tmp[i, j, k] = acc
    end
end

@kernel function real_to_complex_stack_kernel!(dest, src, n1::Int, n2::Int, n3::Int)
    i, j, k = @index(Global, NTuple)
    if i <= n1 && j <= n2 && k <= n3
        T = eltype(src)
        @inbounds dest[i, j, k] = Complex{T}(src[i, j, k], zero(T))
    end
end

@kernel function multiply_kernel_fft_stack_kernel!(fft_stack, kernels_fft, n1::Int, n2::Int, n3::Int)
    i, j, k = @index(Global, NTuple)
    if i <= n1 && j <= n2 && k <= n3
        @inbounds fft_stack[i, j, k] *= kernels_fft[i, j, k]
    end
end

@kernel function complex_to_real_scaled_stack_kernel!(dest, src, scale, n1::Int, n2::Int, n3::Int)
    i, j, k = @index(Global, NTuple)
    if i <= n1 && j <= n2 && k <= n3
        @inbounds dest[i, j, k] = real(src[i, j, k]) * scale
    end
end

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
    host_kernel = Vector{T}(undef, needed)
    @inbounds for k in -half:half
        host_kernel[k + half + 1] = exp(-T(0.5) * (T(k) / sigma)^2)
    end
    host_kernel ./= sum(host_kernel)
    @views copyto!(kernel[1:needed], host_kernel)

    n1, n2 = size(intensity)
    _apply_elongation!(execution_style(intensity), intensity, tmp, kernel, half, n1, n2)
    copyto!(intensity, tmp)
    return kernel
end

function apply_elongation_stack!(intensity_stack::AbstractArray{T,3}, factor::Real, tmp::AbstractArray{T,3},
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
    host_kernel = Vector{T}(undef, needed)
    @inbounds for k in -half:half
        host_kernel[k + half + 1] = exp(-T(0.5) * (T(k) / sigma)^2)
    end
    host_kernel ./= sum(host_kernel)
    @views copyto!(kernel[1:needed], host_kernel)

    n1, n2, n3 = size(intensity_stack)
    _apply_elongation_stack!(execution_style(intensity_stack), intensity_stack, tmp, kernel, half, n1, n2, n3)
    copyto!(intensity_stack, tmp)
    return kernel
end

function _apply_elongation!(::ScalarCPUStyle, intensity::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    kernel::AbstractVector{T}, half::Int, n1::Int, n2::Int) where {T<:AbstractFloat}
    @inbounds for i in 1:n1, j in 1:n2
        acc = zero(T)
        for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += intensity[i, jj] * kernel[k + half + 1]
        end
        tmp[i, j] = acc
    end
    return tmp
end

function _apply_elongation!(style::AcceleratorStyle, intensity::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    kernel::AbstractVector{T}, half::Int, n1::Int, n2::Int) where {T<:AbstractFloat}
    launch_kernel!(style, elongation_apply_kernel!, tmp, intensity, kernel, half, n1, n2; ndrange=size(intensity))
    return tmp
end

function _apply_elongation_stack!(::ScalarCPUStyle, intensity_stack::AbstractArray{T,3}, tmp::AbstractArray{T,3},
    kernel::AbstractVector{T}, half::Int, n1::Int, n2::Int, n3::Int) where {T<:AbstractFloat}
    @inbounds for k in 1:n3, i in 1:n1, j in 1:n2
        acc = zero(T)
        for dk in -half:half
            jj = clamp(j + dk, 1, n2)
            acc += intensity_stack[i, jj, k] * kernel[dk + half + 1]
        end
        tmp[i, j, k] = acc
    end
    return tmp
end

function _apply_elongation_stack!(style::AcceleratorStyle, intensity_stack::AbstractArray{T,3}, tmp::AbstractArray{T,3},
    kernel::AbstractVector{T}, half::Int, n1::Int, n2::Int, n3::Int) where {T<:AbstractFloat}
    launch_kernel!(style, elongation_apply_stack_kernel!, tmp, intensity_stack, kernel, half, n1, n2, n3;
        ndrange=size(intensity_stack))
    return tmp
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Int, wavelength::Real)
    return ARCSEC_PER_RAD * float(wavelength) / float(diameter) / diffraction_padding
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Real, wavelength::Real)
    return ARCSEC_PER_RAD * float(wavelength) / float(diameter) / float(diffraction_padding)
end

@inline function lgs_pixel_scale(diameter::Real, diffraction_padding::Int, src::LGSSource)
    return lgs_pixel_scale(diameter, diffraction_padding, wavelength(src))
end

function lgs_reference_vector(tel::Telescope, x0::Real, y0::Real, altitude::Real)
    T = promote_type(typeof(tel.params.diameter), typeof(x0), typeof(y0), typeof(altitude))
    diameter = T(tel.params.diameter)
    x0_t = T(x0)
    y0_t = T(y0)
    altitude_t = T(altitude)
    x = (diameter / T(4)) * (x0_t / altitude_t)
    y = (diameter / T(4)) * (y0_t / altitude_t)
    z = (diameter^2) / (T(8) * altitude_t) / (T(2) * sqrt(T(3)))
    return (x, y, z)
end

function lgs_spot_shift(vec::Tuple{T,T,T}, tel::Telescope, x_subap::Real, y_subap::Real,
    pixel_scale::Real) where {T<:Real}
    dx0 = vec[1] * (4 / tel.params.diameter)
    dy0 = vec[2] * (4 / tel.params.diameter)
    dx1 = vec[3] * (sqrt(3) * (4 / tel.params.diameter)^2) * x_subap
    dy1 = vec[3] * (sqrt(3) * (4 / tel.params.diameter)^2) * y_subap
    shift_x = ARCSEC_PER_RAD * (dx0 + dx1) / pixel_scale
    shift_y = ARCSEC_PER_RAD * (dy0 + dy1) / pixel_scale
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
    @. fft_buffer = complex(intensity, zero(T))
    execute_fft_plan!(fft_buffer, fft_plan)
    @. fft_buffer *= kernel_fft
    execute_fft_plan!(fft_buffer, ifft_plan)
    scale = T(1) / (n * n)
    @. intensity = real(fft_buffer) * scale
    return intensity
end

function apply_lgs_convolution!(intensity::AbstractMatrix{T}, kernel_fft::AbstractMatrix{Complex{T}},
    fft_buffer::AbstractMatrix{Complex{T}}, fft_plan,
    ifft_buffer::AbstractMatrix{Complex{T}}, ifft_plan) where {T<:AbstractFloat}
    if size(kernel_fft, 1) == 0
        return intensity
    end
    n = size(intensity, 1)
    @. fft_buffer = complex(intensity, zero(T))
    execute_fft_plan!(fft_buffer, fft_plan)
    @. fft_buffer *= kernel_fft
    copyto!(ifft_buffer, fft_buffer)
    execute_fft_plan!(ifft_buffer, ifft_plan)
    scale = T(1) / (n * n)
    @. intensity = real(ifft_buffer) * scale
    return intensity
end

function apply_lgs_convolution_stack!(intensity_stack::AbstractArray{T,3}, kernels_fft::AbstractArray{Complex{T},3},
    fft_stack::AbstractArray{Complex{T},3}, fft_plan, ifft_plan) where {T<:AbstractFloat}
    size(kernels_fft, 3) == 0 && return intensity_stack
    _apply_lgs_convolution_stack!(execution_style(intensity_stack), intensity_stack, kernels_fft, fft_stack, fft_plan, ifft_plan)
    return intensity_stack
end

function _apply_lgs_convolution_stack!(::ScalarCPUStyle, intensity_stack::AbstractArray{T,3}, kernels_fft::AbstractArray{Complex{T},3},
    fft_stack::AbstractArray{Complex{T},3}, fft_plan, ifft_plan) where {T<:AbstractFloat}
    n1, n2, n3 = size(intensity_stack)
    scale = T(1) / (n1 * n2)
    @inbounds for k in 1:n3
        @views fft_stack[:, :, k] .= complex.(intensity_stack[:, :, k], zero(T))
    end
    execute_fft_plan!(fft_stack, fft_plan)
    @. fft_stack *= kernels_fft
    execute_fft_plan!(fft_stack, ifft_plan)
    @inbounds for k in 1:n3
        @views intensity_stack[:, :, k] .= real.(fft_stack[:, :, k]) .* scale
    end
    return intensity_stack
end

function _apply_lgs_convolution_stack!(style::AcceleratorStyle, intensity_stack::AbstractArray{T,3}, kernels_fft::AbstractArray{Complex{T},3},
    fft_stack::AbstractArray{Complex{T},3}, fft_plan, ifft_plan) where {T<:AbstractFloat}
    n1, n2, n3 = size(intensity_stack)
    launch_kernel_async!(style, real_to_complex_stack_kernel!, fft_stack, intensity_stack, n1, n2, n3; ndrange=size(fft_stack))
    synchronize_backend!(style)
    execute_fft_plan!(fft_stack, fft_plan)
    launch_kernel_async!(style, multiply_kernel_fft_stack_kernel!, fft_stack, kernels_fft, n1, n2, n3; ndrange=size(fft_stack))
    synchronize_backend!(style)
    execute_fft_plan!(fft_stack, ifft_plan)
    scale = T(1) / (n1 * n2)
    launch_kernel!(style, complex_to_real_scaled_stack_kernel!, intensity_stack, fft_stack, scale, n1, n2, n3; ndrange=size(fft_stack))
    return intensity_stack
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
    kernel = Matrix{T}(undef, pad, pad)
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

    @. fft_buffer = complex(kernel, zero(T))
    execute_fft_plan!(fft_buffer, fft_plan)
    kernel_fft = similar(fft_buffer, Complex{T}, pad, pad)
    copyto!(kernel_fft, fft_buffer)
    return kernel_fft
end
