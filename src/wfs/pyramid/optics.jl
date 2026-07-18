function build_pyramid_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_pyramid_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_pyramid_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_pyramid_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, pyramid_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_pyramid_mask!(wfs::PyramidWFS, tel::Telescope)
    mask = wfs.front_end.propagation.pyramid_mask
    copyto!(mask, host_pyramid_mask(wfs, tel))
    return mask
end

function host_pyramid_mask(wfs::PyramidWFS, tel::Telescope)
    n = size(wfs.front_end.propagation.pyramid_mask, 1)
    T = eltype(wfs.estimator.state.slopes)
    host = Matrix{Complex{T}}(undef, n, n)
    if wfs.front_end.phase_mask.old_mask
        build_pyramid_mask_old_host!(host, wfs, tel)
    else
        build_pyramid_mask_new_host!(host, wfs, tel)
    end
    return host
end

function build_pyramid_mask_new_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    rooftop_pixels = wfs.front_end.phase_mask.rooftop * wfs.front_end.phase_mask.diffraction_padding / sqrt(T(2))
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(π)
    if wfs.front_end.phase_mask.psf_centering
        xvals = range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        xvals = range(-lim, lim; length=n + 1)[1:n]
    end
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    sx = wfs.estimator.state.shift_x
    sy = wfs.estimator.state.shift_y
    θ = wfs.front_end.phase_mask.theta_rotation
    cθ = cos(θ)
    sθ = sin(θ)
    @inbounds for i in 1:n, j in 1:n
        x_ = xvals[i]
        y_ = xvals[j]
        x = x_ * cθ - y_ * sθ
        y = y_ * cθ + x_ * sθ
        p1 = x * r + x_ * sx[1] + y * r - y_ * sy[1] + rooftop_pixels
        p2 = -x * r + x_ * sx[2] + y * r - y_ * sy[2]
        p3 = -x * r + x_ * sx[3] - y * r - y_ * sy[3] + rooftop_pixels
        p4 = x * r + x_ * sx[4] - y * r - y_ * sy[4]
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function build_pyramid_mask_old_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n_tot = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    sx = wfs.estimator.state.shift_x
    sy = wfs.estimator.state.shift_y
    norma = (T(tel.params.resolution) / T(n_sub)) / 4
    fill!(mask, complex(zero(T), zero(T)))
    if wfs.front_end.phase_mask.psf_centering
        tip = centered_grid(T, n_tot ÷ 2, true)
        tilt = centered_grid(T, n_tot ÷ 2, true)
        Tip = repeat(tip', n_tot ÷ 2, 1)
        Tilt = repeat(tilt, 1, n_tot ÷ 2)
        Tip .-= mean(Tip)
        Tilt .-= mean(Tilt)
        q = T(n_sub + sep) * wfs.front_end.phase_mask.mask_scale
        @views begin
            mask[1:n_tot÷2, 1:n_tot÷2] .= cis.((Tip .* (q + sx[1]) .+ Tilt .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[2]) .+ Tilt .* (q - sy[2])) .* norma)
            mask[n_tot÷2+1:end, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[3]) .- Tilt .* (q + sy[3])) .* norma)
            mask[n_tot÷2+1:end, 1:n_tot÷2] .= cis.((Tip .* (q + sx[4]) .- Tilt .* (q + sy[4])) .* norma)
        end
    else
        d_pix = T(π) / T(n_tot)
        lim_p = T(π)
        lim_m = T(π) - 2 * d_pix
        tip1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip2 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tip4 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tilt2 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt4 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        Tip1 = repeat(tip1', length(tilt1), 1); Tilt1 = repeat(tilt1, 1, length(tip1))
        Tip2 = repeat(tip2', length(tilt2), 1); Tilt2 = repeat(tilt2, 1, length(tip2))
        Tip3 = repeat(tip3', length(tilt3), 1); Tilt3 = repeat(tilt3, 1, length(tip3))
        Tip4 = repeat(tip4', length(tilt4), 1); Tilt4 = repeat(tilt4, 1, length(tip4))
        Tip1 .-= mean(Tip1); Tilt1 .-= mean(Tilt1)
        Tip2 .-= mean(Tip2); Tilt2 .-= mean(Tilt2)
        Tip3 .-= mean(Tip3); Tilt3 .-= mean(Tilt3)
        Tip4 .-= mean(Tip4); Tilt4 .-= mean(Tilt4)
        q = T(n_sub + sep) * wfs.front_end.phase_mask.mask_scale
        @views begin
            mask[1:n_tot÷2+1, 1:n_tot÷2+1] .= cis.((Tip1 .* (q + sx[1]) .+ Tilt1 .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2+1, n_tot÷2:end] .= cis.((-Tip4 .* (q - sx[2]) .+ Tilt4 .* (q - sy[2])) .* norma)
            mask[n_tot÷2:end, n_tot÷2:end] .= cis.((-Tip3 .* (q - sx[3]) .- Tilt3 .* (q + sy[3])) .* norma)
            mask[n_tot÷2:end, 1:n_tot÷2+1] .= cis.((Tip2 .* (q + sx[4]) .- Tilt2 .* (q + sy[4])) .* norma)
        end
    end
    return mask
end

centered_grid(::Type{T}, n::Int, endpoint::Bool) where {T<:AbstractFloat} =
    axis_values(T, n, -T(π), T(π); endpoint=endpoint)

function axis_values(::Type{T}, n::Int, lo::T, hi::T; endpoint::Bool) where {T<:AbstractFloat}
    if endpoint
        return reshape(range(lo, hi; length=n), n, 1)
    end
    return reshape(range(lo; step=(hi - lo) / n, length=n), n, 1)
end

function _build_pyramid_mask!(::ScalarCPUStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    pix_per_subap = T(tel.params.resolution) / T(n_sub)
    norma = pix_per_subap
    lim = T(pi)
    x_vals = if wfs.front_end.phase_mask.psf_centering
        range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        range(-lim, lim; length=n, endpoint=false)
    end
    @inbounds for i in 1:n, j in 1:n
        x = x_vals[i]
        y = x_vals[j]
        p1 = x * r + y * r
        p2 = -x * r + y * r
        p3 = -x * r - y * r
        p4 = x * r - y * r
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function _build_pyramid_mask!(style::AcceleratorStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(pi)
    start = wfs.front_end.phase_mask.psf_centering ? -lim * (one(T) - one(T) / T(n)) : -lim
    step = T(2) * lim / T(n)
    launch_kernel!(style, pyramid_mask_kernel!, mask, r, norma, start, step, n; ndrange=size(mask))
    return mask
end

function accumulate_pyramid_focal_intensity!(out::AbstractMatrix,
    front_end::PyramidOpticalFrontEnd)
    propagation = front_end.propagation
    if front_end.phase_mask.psf_centering
        @. propagation.focal_field = propagation.focal_field * propagation.phasor
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        @. propagation.focal_field = propagation.focal_field * propagation.pyramid_mask
        copyto!(propagation.pupil_field, propagation.focal_field)
        execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    else
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        fftshift2d!(propagation.pupil_field, propagation.focal_field)
        @. propagation.pupil_field = propagation.pupil_field * propagation.pyramid_mask
        execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    end
    @. propagation.temp = abs2(propagation.pupil_field)
    out .+= propagation.temp
    return out
end

@inline accumulate_pyramid_focal_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS{<:Diffractive}) =
    accumulate_pyramid_focal_intensity!(out, wfs.front_end)

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    tel::Telescope, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, tel, src,
        pyramid_operating_modulation(wfs))
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    tel::Telescope, src::AbstractSource,
    modulation::PreparedFocalPlaneModulation) where {T<:AbstractFloat}
    return pyramid_intensity_core!(execution_style(out), out, wfs, tel, src,
        modulation)
end

function pyramid_intensity_core!(::ScalarCPUStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource, modulation::PreparedFocalPlaneModulation) where {
    T<:AbstractFloat,
}
    prepare_pyramid_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    reflectivity = pupil_reflectivity(tel)

    fill!(out, zero(T))
    fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
    @views @. wfs.front_end.propagation.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * sqrt(reflectivity) *
        cispi(opd_to_cycles * tel.state.opd)
    for p in 1:modulation_point_count(modulation)
        copyto!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.field)
        @views @. wfs.front_end.propagation.focal_field[ox+1:ox+n, oy+1:oy+n] *=
            modulation.amplitude_weights[p] * modulation.phases[:, :, p]
        accumulate_pyramid_focal_intensity!(out, wfs)
    end
    return out
end

function pyramid_intensity_core!(::AcceleratorStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource, modulation::PreparedFocalPlaneModulation) where {
    T<:AbstractFloat,
}
    prepare_pyramid_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    reflectivity = pupil_reflectivity(tel)

    fill!(out, zero(T))
    for p in 1:modulation_point_count(modulation)
        fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. wfs.front_end.propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * sqrt(reflectivity) *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * tel.state.opd)
        copyto!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.field)
        accumulate_pyramid_focal_intensity!(out, wfs)
    end
    return out
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, tel, src)
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::ExtendedSource) where {T<:AbstractFloat}
    return accumulate_pyramid_extended_intensity!(execution_style(out), out, wfs, tel, src)
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    pyramid_intensity_core!(out, wfs, tel, src)
    apply_lgs_elongation!(lgs_profile(src), out, wfs, tel, src)
    return out
end

function pyramid_modulation_frame!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource) where {T<:AbstractFloat}
    prepare_pyramid_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.front_end.propagation.field, 1)
    fft_scale2 = inv(T(pad) * T(pad))
    if size(out) != (pad, pad)
        throw(DimensionMismatchError("modulation frame size must match pyramid sampling"))
    end
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    modulation = pyramid_operating_modulation(wfs)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    reflectivity = pupil_reflectivity(tel)

    fill!(out, zero(T))
    for p in 1:modulation_point_count(modulation)
        fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. wfs.front_end.propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * sqrt(reflectivity) *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * tel.state.opd)
        copyto!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.field)
        if wfs.front_end.phase_mask.psf_centering
            @. wfs.front_end.propagation.focal_field = wfs.front_end.propagation.focal_field * wfs.front_end.propagation.phasor
            execute_fft_plan!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.fft_plan)
        else
            execute_fft_plan!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.fft_plan)
            @. wfs.front_end.propagation.intensity = abs2(wfs.front_end.propagation.focal_field) * fft_scale2
            fftshift2d!(wfs.front_end.propagation.temp, wfs.front_end.propagation.intensity)
            out .+= wfs.front_end.propagation.temp
            continue
        end
        @. wfs.front_end.propagation.temp = abs2(wfs.front_end.propagation.focal_field) * fft_scale2
        out .+= wfs.front_end.propagation.temp
    end
    return out
end

"""
    pyramid_modulation_frame(wfs, telescope, source)

Allocate and form the focal-plane modulation-cycle intensity used by a
gain-sensing camera. The returned array has the prepared Pyramid front end's
native focal-plane sampling.
"""
function pyramid_modulation_frame(wfs::PyramidWFS{<:Diffractive},
    tel::Telescope, src::AbstractSource)
    out = similar(wfs.front_end.propagation.intensity)
    return pyramid_modulation_frame!(out, wfs, tel, src)
end

function apply_lgs_elongation!(::LGSProfileNone, out::AbstractMatrix{T}, wfs::PyramidWFS,
    ::Telescope, src::LGSSource) where {T<:AbstractFloat}
    wfs.front_end.propagation.elongation_kernel = apply_elongation!(
        out,
        lgs_elongation_factor(src),
        wfs.front_end.propagation.scratch,
        wfs.front_end.propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, out::AbstractMatrix{T}, wfs::PyramidWFS,
    tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, tel, src)
    apply_lgs_convolution!(
        out,
        wfs.front_end.propagation.lgs_kernel_fft,
        wfs.front_end.propagation.focal_field,
        wfs.front_end.propagation.fft_plan,
        wfs.front_end.propagation.pupil_field,
        wfs.front_end.propagation.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    pad = size(wfs.front_end.propagation.intensity, 1)
    padding = wfs.front_end.propagation.effective_resolution / tel.params.resolution
    pixel_scale = lgs_pixel_scale(tel.params.diameter, padding,
        wavelength(src))
    tag = lgs_kernel_signature(
        tel,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        eltype(wfs.front_end.propagation.intensity);
        model=:subaperture_average,
    )
    if size(wfs.front_end.propagation.lgs_kernel_fft, 1) == pad && wfs.front_end.propagation.lgs_kernel_tag == tag
        return wfs
    end
    wfs.front_end.propagation.lgs_kernel_fft = lgs_average_kernel_fft(
        tel,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        wfs.front_end.propagation.focal_field,
        wfs.front_end.propagation.fft_plan,
    )
    wfs.front_end.propagation.lgs_kernel_tag = tag
    return wfs
end
