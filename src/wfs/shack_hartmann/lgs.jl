function apply_lgs_elongation!(::LGSProfileNone, wfs::ShackHartmannWFS, ::Telescope, src::LGSSource, ::Int)
    wfs.front_end.propagation.elongation_kernel = apply_elongation!(
        wfs.front_end.propagation.intensity,
        lgs_elongation_factor(src),
        wfs.front_end.propagation.temp,
        wfs.front_end.propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource, idx::Int)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution!(
        wfs.front_end.propagation.intensity,
        wfs.front_end.propagation.lgs_kernel_fft,
        wfs.front_end.propagation.fft_buffer,
        wfs.front_end.propagation.fft_plan,
        wfs.front_end.propagation.ifft_plan,
        idx,
    )
    return wfs
end

function ensure_lgs_kernels!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    dimensions = (tel.params.resolution, tel.params.resolution)
    ensure_lgs_kernels!(wfs.front_end, src, dimensions,
        tel.params.diameter,
        tel.aperture.sampling_m, tel.aperture.origin_m, wavelength(src))
    return wfs
end

function ensure_lgs_kernels!(front_end::ShackHartmannOpticalFrontEnd,
    src::LGSSource,
    pupil_dimensions::NTuple{2,Int}, pupil_diameter::Real,
    pupil_sampling::NTuple{2,<:Real}, pupil_origin::NTuple{2,<:Real},
    wavelength_m::Real)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return front_end
    end
    propagation = front_end.propagation
    pad = size(propagation.intensity, 1)
    n_sub = n_lenslets(front_end)
    pixel_scale = lgs_pixel_scale(
        pupil_diameter / n_sub,
        propagation.effective_padding,
        wavelength_m,
    )
    tag = lgs_kernel_signature(
        src,
        pad,
        n_sub,
        pixel_scale,
        eltype(propagation.intensity),
        pupil_dimensions,
        pupil_diameter,
        pupil_sampling,
        pupil_origin;
        wavelength_m,
        model=:per_subaperture,
        threshold=sh_threshold_convolution(front_end),
    )
    if size(propagation.lgs_kernel_fft, 1) == pad &&
        size(propagation.lgs_kernel_fft, 3) == n_sub * n_sub &&
        propagation.lgs_kernel_tag == tag
        return front_end
    end
    propagation.lgs_kernel_fft = lgs_spot_kernels_fft(
        pupil_diameter, front_end, src, pad, wavelength_m)
    propagation.lgs_kernel_tag = tag
    return front_end
end

function apply_lgs_convolution!(intensity::AbstractMatrix{T}, kernels_fft::AbstractArray{Complex{T},3},
    fft_buffer::AbstractMatrix{Complex{T}}, fft_plan, ifft_plan, idx::Int) where {T<:AbstractFloat}
    if size(kernels_fft, 3) < idx
        return intensity
    end
    kernel = @view kernels_fft[:, :, idx]
    apply_lgs_convolution!(intensity, kernel, fft_buffer, fft_plan, ifft_plan)
    return intensity
end

function lgs_spot_kernels_fft(tel::Telescope, wfs::ShackHartmannWFS, src::LGSSource, pad::Int)
    return lgs_spot_kernels_fft(tel.params.diameter, wfs.front_end, src,
        pad,
        wavelength(src))
end

function lgs_spot_kernels_fft(pupil_diameter::Real,
    front_end::ShackHartmannOpticalFrontEnd, src::LGSSource, pad::Int,
    wavelength_m::Real)
    propagation = front_end.propagation
    T = eltype(propagation.intensity)
    n_sub = n_lenslets(front_end)
    na_profile = src.params.na_profile
    altitudes = na_profile[1, :]
    weights = na_profile[2, :]
    if length(altitudes) == 0
        return similar(propagation.fft_buffer, Complex{T}, 0, 0, 0)
    end

    pixel_scale = lgs_pixel_scale(
        pupil_diameter / n_lenslets(front_end),
        propagation.effective_padding,
        wavelength_m,
    )
    fwhm_px = src.params.fwhm_spot_up / pixel_scale
    sigma_px = fwhm_px / (2 * sqrt(2 * log(T(2))))
    center = (pad + 1) / 2

    x_subap = range(-pupil_diameter / 2, pupil_diameter / 2; length=n_sub)
    y_subap = range(-pupil_diameter / 2, pupil_diameter / 2; length=n_sub)
    kernels_fft = similar(propagation.fft_buffer, Complex{T}, pad, pad,
        n_sub * n_sub)

    x0 = src.params.laser_coordinates[2]
    y0 = -src.params.laser_coordinates[1]
    ref_idx = Int(cld(length(altitudes), 2))
    ref_vec = lgs_reference_vector(pupil_diameter, x0, y0,
        altitudes[ref_idx])

    idx = 1
    kernel = Matrix{T}(undef, pad, pad)
    for iy in 1:n_sub, ix in 1:n_sub
        lgs_spot_kernel!(kernel, pupil_diameter, src, altitudes, weights,
            pixel_scale, sigma_px, center,
            x_subap[ix], y_subap[iy], ref_vec)
        peak = maximum(kernel)
        if peak > 0
            cutoff = sh_threshold_convolution(front_end) * peak
            @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
                if kernel[i, j] < cutoff
                    kernel[i, j] = zero(T)
                end
            end
            total = sum(kernel)
            if total > 0
                kernel ./= total
            end
        end
        fft_buffer = propagation.fft_buffer
        _copy_host_real_kernel_to_complex!(execution_style(fft_buffer),
            fft_buffer, kernel)
        execute_fft_plan!(fft_buffer, propagation.fft_plan)
        @views kernels_fft[:, :, idx] .= fft_buffer
        idx += 1
    end

    return kernels_fft
end
