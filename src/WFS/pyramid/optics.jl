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
    mask = wfs.state.pyramid_mask
    copyto!(mask, host_pyramid_mask(wfs, tel))
    return mask
end

function host_pyramid_mask(wfs::PyramidWFS, tel::Telescope)
    n = size(wfs.state.pyramid_mask, 1)
    T = eltype(wfs.state.slopes)
    host = Matrix{Complex{T}}(undef, n, n)
    if wfs.params.old_mask
        build_pyramid_mask_old_host!(host, wfs, tel)
    else
        build_pyramid_mask_new_host!(host, wfs, tel)
    end
    return host
end

function build_pyramid_mask_new_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    rooftop_pixels = wfs.params.rooftop * wfs.params.diffraction_padding / sqrt(T(2))
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(π)
    if wfs.params.psf_centering
        xvals = range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        xvals = range(-lim, lim; length=n + 1)[1:n]
    end
    r = (T(n_sub) + T(sep)) / 2
    sx = wfs.state.shift_x
    sy = wfs.state.shift_y
    θ = wfs.params.theta_rotation
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
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    sx = wfs.state.shift_x
    sy = wfs.state.shift_y
    norma = (T(tel.params.resolution) / T(n_sub)) / 4
    fill!(mask, complex(zero(T), zero(T)))
    if wfs.params.psf_centering
        tip = centered_grid(T, n_tot ÷ 2, true)
        tilt = centered_grid(T, n_tot ÷ 2, true)
        Tip = repeat(tip', n_tot ÷ 2, 1)
        Tilt = repeat(tilt, 1, n_tot ÷ 2)
        Tip .-= mean(Tip)
        Tilt .-= mean(Tilt)
        q = T(n_sub + sep)
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
        q = T(n_sub + sep)
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
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.params.mask_scale / 2
    pix_per_subap = T(tel.params.resolution) / T(n_sub)
    norma = pix_per_subap
    lim = T(pi)
    x_vals = if wfs.params.psf_centering
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
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.params.mask_scale / 2
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(pi)
    start = wfs.params.psf_centering ? -lim * (one(T) - one(T) / T(n)) : -lim
    step = T(2) * lim / T(n)
    launch_kernel!(style, pyramid_mask_kernel!, mask, r, norma, start, step, n; ndrange=size(mask))
    return mask
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::AbstractSource) where {T<:AbstractFloat}
    prepare_pyramid_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2 /
        wfs.params.modulation_points))
    fft_scale2 = inv(T(pad) * T(pad))

    fill!(out, zero(T))
    for p in 1:wfs.params.modulation_points
        fill!(wfs.state.field, zero(eltype(wfs.state.field)))
        @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
            wfs.state.modulation_phases[:, :, p] * cispi(opd_to_cycles * tel.state.opd)
        copyto!(wfs.state.focal_field, wfs.state.field)
        if wfs.params.psf_centering
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.pyramid_mask
            copyto!(wfs.state.pupil_field, wfs.state.focal_field)
            execute_fft_plan!(wfs.state.pupil_field, wfs.state.ifft_plan)
        else
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
            fftshift2d!(wfs.state.pupil_field, wfs.state.focal_field)
            @. wfs.state.pupil_field = wfs.state.pupil_field * wfs.state.pyramid_mask
            execute_fft_plan!(wfs.state.pupil_field, wfs.state.ifft_plan)
        end
        @. wfs.state.temp = abs2(wfs.state.pupil_field)
        out .+= wfs.state.temp
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
    pad = size(wfs.state.field, 1)
    fft_scale2 = inv(T(pad) * T(pad))
    if size(out) != (pad, pad)
        throw(DimensionMismatchError("modulation frame size must match pyramid sampling"))
    end
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2 /
        wfs.params.modulation_points))

    fill!(out, zero(T))
    for p in 1:wfs.params.modulation_points
        fill!(wfs.state.field, zero(eltype(wfs.state.field)))
        @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
            wfs.state.modulation_phases[:, :, p] * cispi(opd_to_cycles * tel.state.opd)
        copyto!(wfs.state.focal_field, wfs.state.field)
        if wfs.params.psf_centering
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
        else
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
            @. wfs.state.intensity = abs2(wfs.state.focal_field) * fft_scale2
            fftshift2d!(wfs.state.temp, wfs.state.intensity)
            out .+= wfs.state.temp
            continue
        end
        @. wfs.state.temp = abs2(wfs.state.focal_field) * fft_scale2
        out .+= wfs.state.temp
    end
    return out
end

function apply_lgs_elongation!(::LGSProfileNone, out::AbstractMatrix{T}, wfs::PyramidWFS,
    ::Telescope, src::LGSSource) where {T<:AbstractFloat}
    wfs.state.elongation_kernel = apply_elongation!(
        out,
        lgs_elongation_factor(src),
        wfs.state.scratch,
        wfs.state.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, out::AbstractMatrix{T}, wfs::PyramidWFS,
    tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, tel, src)
    apply_lgs_convolution!(
        out,
        wfs.state.lgs_kernel_fft,
        wfs.state.focal_field,
        wfs.state.fft_plan,
        wfs.state.pupil_field,
        wfs.state.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    pad = size(wfs.state.intensity, 1)
    tag = objectid(na_profile) ⊻ hash(src.params.laser_coordinates) ⊻ hash(src.params.fwhm_spot_up) ⊻
        hash(pad) ⊻ hash(wfs.params.n_subap)
    if size(wfs.state.lgs_kernel_fft, 1) == pad && wfs.state.lgs_kernel_tag == tag
        return wfs
    end
    padding = wfs.state.effective_resolution / tel.params.resolution
    pixel_scale = lgs_pixel_scale(tel.params.diameter, padding, wavelength(src))
    wfs.state.lgs_kernel_fft = lgs_average_kernel_fft(
        tel,
        src,
        pad,
        wfs.params.n_subap,
        pixel_scale,
        wfs.state.focal_field,
        wfs.state.fft_plan,
    )
    wfs.state.lgs_kernel_tag = tag
    return wfs
end

function build_modulation_phases!(wfs::PyramidWFS, tel::Telescope)
    phases = wfs.state.modulation_phases
    copyto!(phases, host_modulation_phases(eltype(wfs.state.slopes), tel, wfs.params.modulation,
        wfs.params.modulation_points, wfs.params.delta_theta, wfs.params.user_modulation_path))
    return phases
end

function host_modulation_phases(::Type{T}, tel::Telescope, modulation::T, n_pts::Int,
    delta_theta::T, user_modulation_path) where {T<:AbstractFloat}
    n = tel.params.resolution
    phases = Array{Complex{T}}(undef, n, n, n_pts)
    if n_pts == 1 && modulation == 0 && user_modulation_path === nothing
        fill!(phases, complex(one(T), zero(T)))
        return phases
    end
    tip_vals = range(-T(π), T(π); length=n)
    Tilt = reshape(tip_vals, n, 1)
    Tip = reshape(tip_vals, 1, n)
    pupil = Array(tel.state.pupil)
    if user_modulation_path === nothing
        @inbounds for p in 1:n_pts
            θ = delta_theta + T(2 * π * (p - 1) / n_pts)
            sθ, cθ = sincos(θ)
            mx = modulation * cθ
            my = modulation * sθ
            for i in 1:n, j in 1:n
                phase = pupil[i, j] ? (mx * Tip[j] + my * Tilt[i]) : zero(T)
                phases[i, j, p] = cis(phase)
            end
        end
    else
        @inbounds for p in 1:n_pts
            mx = T(user_modulation_path[p][1])
            my = T(user_modulation_path[p][2])
            for i in 1:n, j in 1:n
                phase = pupil[i, j] ? (mx * Tip[j] + my * Tilt[i]) : zero(T)
                phases[i, j, p] = cis(phase)
            end
        end
    end
    return phases
end

function _build_modulation_phases!(::ScalarCPUStyle, phases::AbstractArray{Complex{T},3}, mod, center, n_pts::Int, n::Int) where {T<:AbstractFloat}
    x_coords = Vector{T}(undef, n)
    for i in 1:n
        x_coords[i] = (i - center) / n
    end
    @inbounds for p in 1:n_pts
        angle = 2 * pi * (p - 1) / n_pts
        s, c = sincos(angle)
        for i in 1:n, j in 1:n
            phase = 2 * pi * mod * (x_coords[i] * c + x_coords[j] * s)
            phases[i, j, p] = cis(phase)
        end
    end
    return phases
end

function _build_modulation_phases!(style::AcceleratorStyle, phases::AbstractArray{Complex{T},3}, mod, center, n_pts::Int, n::Int) where {T<:AbstractFloat}
    launch_kernel!(style, modulation_phases_kernel!, phases, T(mod), T(center), n_pts, n; ndrange=size(phases))
    return phases
end
