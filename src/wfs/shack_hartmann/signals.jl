
@inline function centroid_from_intensity!(intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    return centroid_from_intensity!(execution_style(intensity), intensity, threshold)
end

@inline function centroid_from_intensity!(::ScalarCPUStyle, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    peak = maximum(intensity)
    if peak <= 0
        return zero(T), zero(T), zero(T)
    end
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), intensity, threshold * peak)
end

@inline function centroid_from_intensity_cutoff!(intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_intensity_cutoff!(execution_style(intensity), intensity, cutoff)
end

@inline function centroid_from_intensity_cutoff!(::ScalarCPUStyle, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    total = zero(T)
    sx = zero(T)
    sy = zero(T)
    n1 = size(intensity, 1)
    n2 = size(intensity, 2)
    @inbounds for x in 1:n1, y in 1:n2
        val = intensity[x, y]
        if val < cutoff
            intensity[x, y] = zero(T)
        else
            total += val
            sx += T(x - 1) * val
            sy += T(y - 1) * val
        end
    end
    if total <= 0
        return zero(T), zero(T), zero(T)
    end
    return total, sx / total, sy / total
end

@inline function centroid_from_intensity!(::AcceleratorStyle, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    host_intensity = Array(intensity)
    return centroid_from_intensity!(ScalarCPUStyle(), host_intensity, threshold)
end

@inline function centroid_from_intensity_cutoff!(::AcceleratorStyle, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    host_intensity = Array(intensity)
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), host_intensity, cutoff)
end

@inline function centroid_from_spot!(wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    return centroid_from_spot!(execution_style(intensity), wfs, intensity, threshold)
end

@inline function centroid_from_spot!(wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    return centroid_from_spot!(wfs, intensity, centroid_threshold(wfs))
end

@inline function centroid_from_spot!(::ScalarCPUStyle, ::ShackHartmannWFS, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    return centroid_from_intensity!(ScalarCPUStyle(), intensity, threshold)
end

@inline function centroid_from_spot!(::AcceleratorStyle, wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    if size(wfs.state.centroid_host) != size(intensity)
        wfs.state.centroid_host = Matrix{T}(undef, size(intensity)...)
    end
    copyto!(wfs.state.centroid_host, intensity)
    return centroid_from_intensity!(ScalarCPUStyle(), wfs.state.centroid_host, threshold)
end

@inline function centroid_from_spot_cutoff!(wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_spot_cutoff!(execution_style(intensity), wfs, intensity, cutoff)
end

@inline function centroid_from_spot_cutoff!(::ScalarCPUStyle, ::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), intensity, cutoff)
end

@inline function centroid_from_spot_cutoff!(::AcceleratorStyle, wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    if size(wfs.state.centroid_host) != size(intensity)
        wfs.state.centroid_host = Matrix{T}(undef, size(intensity)...)
    end
    copyto!(wfs.state.centroid_host, intensity)
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), wfs.state.centroid_host, cutoff)
end

function sh_signal_from_spots_device_stats!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = wfs.params.n_lenslets
    offset = n_sub * n_sub
    fill!(wfs.state.spot_stats, zero(eltype(wfs.state.spot_stats)))
    launch_kernel!(style, sh_spot_cutoff_stats_kernel!, wfs.state.spot_stats, wfs.state.spot_cube,
        wfs.state.valid_mask, cutoff, n_sub, size(wfs.state.spot_cube, 2),
        size(wfs.state.spot_cube, 3); ndrange=offset)
    launch_kernel!(style, sh_finalize_spot_slopes_kernel!, wfs.state.slopes, wfs.state.spot_stats,
        wfs.state.valid_mask, n_sub, offset; ndrange=offset)
    return wfs.state.slopes
end

function sh_signal_from_spots_calibrated_device_stats!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = wfs.params.n_lenslets
    offset = n_sub * n_sub
    reference = vec(wfs.state.reference_signal_2d)
    fill!(wfs.state.spot_stats, zero(eltype(wfs.state.spot_stats)))
    launch_kernel!(style, sh_spot_cutoff_stats_kernel!, wfs.state.spot_stats, wfs.state.spot_cube,
        wfs.state.valid_mask, cutoff, n_sub, size(wfs.state.spot_cube, 2),
        size(wfs.state.spot_cube, 3); ndrange=offset)
    launch_kernel!(style, sh_finalize_spot_slopes_reference_scale_kernel!, wfs.state.slopes, wfs.state.spot_stats,
        reference, wfs.state.valid_mask, wfs.state.slopes_units, n_sub, offset; ndrange=offset)
    return wfs.state.slopes
end

@inline function sh_spot_view(wfs::ShackHartmannWFS, idx::Int)
    return @view wfs.state.spot_cube[idx, :, :]
end

@inline function sync_sh_staged_spot!(style::AcceleratorStyle, spot::AbstractMatrix)
    synchronize_backend!(style)
    synchronize_backend!(execution_style(spot))
    return spot
end

@inline function sync_sh_staged_view!(style::AcceleratorStyle, spot_view::AbstractMatrix)
    synchronize_backend!(style)
    synchronize_backend!(execution_style(spot_view))
    return spot_view
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    prepare_sampling!(wfs, tel, src)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource)
    ref = spectral_reference_source(src)
    prepare_sampling!(wfs, tel, ref)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource)
    ast = extended_source_asterism(src)
    prepare_sampling!(wfs, tel, ast.sources[1])
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_lenslets
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            sample_spot!(wfs, wfs.state.intensity)
            copyto!(spot_view, wfs.state.spot)
            peak = max(peak, sh_safe_peak_value(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n = tel.params.resolution
        n_sub = wfs.params.n_lenslets
        sub = div(n, n_sub)
        pad = size(wfs.state.field, 1)
        ox = div(pad - sub, 2)
        oy = div(pad - sub, 2)
        peak = zero(eltype(wfs.state.slopes))
        idx = 1
        @inbounds for i in 1:n_sub, j in 1:n_sub
            spot_view = sh_spot_view(wfs, idx)
            if wfs.state.valid_mask_host[i, j]
                xs = (i - 1) * sub + 1
                ys = (j - 1) * sub + 1
                xe = min(i * sub, n)
                ye = min(j * sub, n)
                compute_intensity_safe!(style, wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
                sample_spot!(wfs, wfs.state.intensity)
                sync_sh_staged_spot!(style, wfs.state.spot)
                copyto!(spot_view, wfs.state.spot)
                sync_sh_staged_view!(style, spot_view)
                peak = max(peak, sh_safe_peak_value(spot_view))
            else
                fill!(spot_view, zero(eltype(spot_view)))
            end
            idx += 1
        end
        return peak
    end
    compute_intensity_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource)
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    peak = zero(eltype(wfs.state.slopes))
    total_flux = photon_flux(src)
    @inbounds for sample in src.bundle.samples
        variant = source_with_wavelength_and_flux(src, sample.wavelength,
            eltype(wfs.state.slopes)(total_flux * sample.weight))
        peak = max(peak, sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, variant))
        wfs.state.spot_cube_accum .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    return max(peak, sh_safe_peak_value(wfs.state.spot_cube))
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
        peak = zero(eltype(wfs.state.slopes))
        total_flux = photon_flux(src)
        @inbounds for sample in src.bundle.samples
            variant = source_with_wavelength_and_flux(src, sample.wavelength,
                eltype(wfs.state.slopes)(total_flux * sample.weight))
            peak = max(peak, sampled_spots_peak!(style, wfs, tel, variant))
            @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
        end
        copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
        return max(peak, sh_safe_peak_value(wfs.state.spot_cube))
    end
    if is_lgs_source(src)
        fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
        peak = zero(eltype(wfs.state.slopes))
        total_flux = photon_flux(src)
        @inbounds for sample in src.bundle.samples
            variant = source_with_wavelength_and_flux(src, sample.wavelength,
                eltype(wfs.state.slopes)(total_flux * sample.weight))
            peak = max(peak, sampled_spots_peak!(style, wfs, tel, variant))
            @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
        end
        copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
        return max(peak, sh_safe_peak_value(wfs.state.spot_cube))
    end
    compute_intensity_spectral_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, ast.sources[1])
    end
    return sampled_spots_peak_asterism_stacked!(ScalarCPUStyle(), wfs, tel, ast)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(style, wfs, tel, ast.sources[1])
    end
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        return sampled_spots_peak_asterism_stacked!(style, wfs, tel, ast)
    end
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    peak = zero(eltype(wfs.state.slopes))
    @inbounds for sample in ast.sources
        peak = max(peak, sampled_spots_peak!(style, wfs, tel, sample))
        @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    return max(peak, sh_safe_peak_value(wfs.state.spot_cube))
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    prepare_sampling!(wfs, tel, src)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_lenslets
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
            sample_spot!(wfs, wfs.state.intensity)
            copyto!(spot_view, wfs.state.spot)
            peak = max(peak, sh_safe_peak_value(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, tel, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    prepare_sampling!(wfs, tel, src)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    ref = spectral_reference_source(src)
    prepare_sampling!(wfs, tel, ref)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    prepare_sampling!(wfs, tel, ast.sources[1])
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(ScalarCPUStyle(), wfs, tel, src)
    sample_spot_stack!(ScalarCPUStyle(), wfs)
    capture_sampled_spot_stack!(wfs, det, rng)
    zero_invalid_sh_spot_cube!(ScalarCPUStyle(), wfs.state.spot_cube, wfs.state.valid_mask)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n = tel.params.resolution
        n_sub = wfs.params.n_lenslets
        sub = div(n, n_sub)
        pad = size(wfs.state.field, 1)
        ox = div(pad - sub, 2)
        oy = div(pad - sub, 2)
        peak = zero(eltype(wfs.state.slopes))
        idx = 1
        @inbounds for i in 1:n_sub, j in 1:n_sub
            spot_view = sh_spot_view(wfs, idx)
            if wfs.state.valid_mask_host[i, j]
                xs = (i - 1) * sub + 1
                ys = (j - 1) * sub + 1
                xe = min(i * sub, n)
                ye = min(j * sub, n)
                compute_intensity_safe!(style, wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
                sample_spot!(wfs, wfs.state.intensity)
                sync_sh_staged_spot!(style, wfs.state.spot)
                frame = capture!(det, wfs.state.spot; rng=rng)
                sync_sh_staged_spot!(style, frame)
                copyto!(spot_view, frame)
                sync_sh_staged_view!(style, spot_view)
                peak = max(peak, sh_safe_peak_value(spot_view))
            else
                fill!(spot_view, zero(eltype(spot_view)))
            end
            idx += 1
        end
        return peak
    end
    compute_intensity_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_lenslets
    capture_sampled_spot_stack!(wfs, det, rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    total_flux = photon_flux(src)
    @inbounds for sample in src.bundle.samples
        variant = source_with_wavelength_and_flux(src, sample.wavelength,
            eltype(wfs.state.slopes)(total_flux * sample.weight))
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, variant)
        wfs.state.spot_cube_accum .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    capture_stack!(det, wfs.state.spot_cube, wfs.state.spot_cube_accum; rng=rng)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
        total_flux = photon_flux(src)
        @inbounds for sample in src.bundle.samples
            variant = source_with_wavelength_and_flux(src, sample.wavelength,
                eltype(wfs.state.slopes)(total_flux * sample.weight))
            sampled_spots_peak!(style, wfs, tel, variant, det, rng)
            @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
        end
        copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
        capture_stack!(det, wfs.state.spot_cube, wfs.state.spot_cube_accum; rng=rng)
        return sh_safe_peak_value(wfs.state.spot_cube)
    end
    if is_lgs_source(src)
        fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
        total_flux = photon_flux(src)
        @inbounds for sample in src.bundle.samples
            variant = source_with_wavelength_and_flux(src, sample.wavelength,
                eltype(wfs.state.slopes)(total_flux * sample.weight))
            sampled_spots_peak!(style, wfs, tel, variant)
            @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
        end
        copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
        n_sub = wfs.params.n_lenslets
        capture_stack!(det, wfs.state.spot_cube, wfs.state.spot_cube_accum; rng=rng)
        launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
            n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
        return sh_safe_peak_value(wfs.state.spot_cube)
    end
    compute_intensity_spectral_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_lenslets
    capture_sampled_spot_stack!(wfs, det, rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, ast.sources[1], det, rng)
    end
    return sampled_spots_peak_asterism_stacked!(ScalarCPUStyle(), wfs, tel, ast, det, rng)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(style, wfs, tel, ast.sources[1], det, rng)
    end
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        return sampled_spots_peak_asterism_stacked!(style, wfs, tel, ast, det, rng)
    end
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    peak = zero(eltype(wfs.state.slopes))
    @inbounds for sample in ast.sources
        peak = max(peak, sampled_spots_peak!(style, wfs, tel, sample, det, rng))
        @. wfs.state.spot_cube_accum = wfs.state.spot_cube_accum + wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    return max(peak, sh_safe_peak_value(wfs.state.spot_cube))
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    prepare_sampling!(wfs, tel, src)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = tel.params.resolution
    n_sub = wfs.params.n_lenslets
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
            sample_spot!(wfs, wfs.state.intensity)
            frame = capture!(det, wfs.state.spot; rng=rng)
            copyto!(spot_view, frame)
            peak = max(peak, sh_safe_peak_value(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, tel, src, det, rng)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    compute_intensity_stack!(style, wfs, tel, src)
    tmp_view = @view wfs.state.intensity_tmp_stack[:, :, 1:size(wfs.state.intensity_stack, 3)]
    apply_elongation_stack!(wfs.state.intensity_stack, lgs_elongation_factor(src),
        tmp_view, wfs.state.elongation_kernel)
    sample_spot_stack!(style, wfs)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(style, wfs, tel, src)
    tmp_view = @view wfs.state.intensity_tmp_stack[:, :, 1:size(wfs.state.intensity_stack, 3)]
    apply_elongation_stack!(wfs.state.intensity_stack, lgs_elongation_factor(src),
        tmp_view, wfs.state.elongation_kernel)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_lenslets
    capture_sampled_spot_stack!(wfs, det, rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    n_sub = wfs.params.n_lenslets
    compute_intensity_stack!(style, wfs, tel, src)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution_stack!(wfs.state.intensity_stack, wfs.state.lgs_kernel_fft,
        wfs.state.fft_stack, wfs.state.fft_stack_plan, wfs.state.ifft_stack_plan)
    sample_spot_stack!(style, wfs)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n_sub = wfs.params.n_lenslets
    compute_intensity_stack!(style, wfs, tel, src)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution_stack!(wfs.state.intensity_stack, wfs.state.lgs_kernel_fft,
        wfs.state.fft_stack, wfs.state.fft_stack_plan, wfs.state.ifft_stack_plan)
    sample_spot_stack!(style, wfs)
    capture_sampled_spot_stack!(wfs, det, rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sh_signal_from_spots!(wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    return sh_signal_from_spots!(execution_style(wfs.state.slopes), wfs, cutoff)
end

function sh_signal_from_spots!(::ScalarCPUStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = wfs.params.n_lenslets
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if wfs.state.valid_mask[i, j]
            total, sx, sy = centroid_from_spot_cube_cutoff!(wfs.state.spot_cube, idx, cutoff)
            if total <= 0
                wfs.state.slopes[idx] = zero(T)
                wfs.state.slopes[idx + n_sub * n_sub] = zero(T)
            else
                wfs.state.slopes[idx] = sy
                wfs.state.slopes[idx + n_sub * n_sub] = sx
            end
        else
            wfs.state.slopes[idx] = zero(T)
            wfs.state.slopes[idx + n_sub * n_sub] = zero(T)
        end
        idx += 1
    end
    return wfs.state.slopes
end

function sh_signal_from_spots!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n_sub = wfs.params.n_lenslets
        offset = n_sub * n_sub
        host_slopes = wfs.state.slopes_host
        idx = 1
        @inbounds for i in 1:n_sub, j in 1:n_sub
            if wfs.state.valid_mask_host[i, j]
                total, sx, sy = centroid_from_spot_cutoff!(wfs, sh_spot_view(wfs, idx), cutoff)
                if total <= 0
                    host_slopes[idx] = zero(T)
                    host_slopes[idx + offset] = zero(T)
                else
                    host_slopes[idx] = sy
                    host_slopes[idx + offset] = sx
                end
            else
                host_slopes[idx] = zero(T)
                host_slopes[idx + offset] = zero(T)
            end
            idx += 1
        end
        copyto!(wfs.state.slopes, host_slopes)
        return wfs.state.slopes
    end
    if sh_uses_device_stats_sensing_plan(wfs)
        return sh_signal_from_spots_device_stats!(style, wfs, cutoff)
    end
    n_sub = wfs.params.n_lenslets
    offset = n_sub * n_sub
    launch_kernel!(style, sh_spot_centroid_kernel!, wfs.state.slopes, wfs.state.spot_cube,
        wfs.state.valid_mask, cutoff, n_sub, offset, size(wfs.state.spot_cube, 2),
        size(wfs.state.spot_cube, 3); ndrange=offset)
    return wfs.state.slopes
end

function sh_signal_from_spots_calibrated!(wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    return sh_signal_from_spots_calibrated!(execution_style(wfs.state.slopes), wfs, cutoff)
end

function sh_signal_from_spots_calibrated!(::ScalarCPUStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    sh_signal_from_spots!(ScalarCPUStyle(), wfs, cutoff)
    subtract_reference_and_scale!(ScalarCPUStyle(), wfs)
    return wfs.state.slopes
end

function sh_signal_from_spots_calibrated!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n_sub = wfs.params.n_lenslets
        offset = n_sub * n_sub
        host_slopes = wfs.state.slopes_host
        reference = wfs.state.reference_signal_host
        inv_units = inv(wfs.state.slopes_units)
        idx = 1
        @inbounds for i in 1:n_sub, j in 1:n_sub
            if wfs.state.valid_mask_host[i, j]
                total, sx, sy = centroid_from_spot_cutoff!(wfs, sh_spot_view(wfs, idx), cutoff)
                if total <= 0
                    host_slopes[idx] = -reference[idx] * inv_units
                    host_slopes[idx + offset] = -reference[idx + offset] * inv_units
                else
                    host_slopes[idx] = (sy - reference[idx]) * inv_units
                    host_slopes[idx + offset] = (sx - reference[idx + offset]) * inv_units
                end
            else
                host_slopes[idx] = zero(T)
                host_slopes[idx + offset] = zero(T)
            end
            idx += 1
        end
        copyto!(wfs.state.slopes, host_slopes)
        return wfs.state.slopes
    end
    if sh_uses_device_stats_sensing_plan(wfs)
        return sh_signal_from_spots_calibrated_device_stats!(style, wfs, cutoff)
    end
    n_sub = wfs.params.n_lenslets
    offset = n_sub * n_sub
    reference = vec(wfs.state.reference_signal_2d)
    launch_kernel!(style, sh_spot_centroid_reference_scale_kernel!, wfs.state.slopes, wfs.state.spot_cube,
        reference, wfs.state.valid_mask, cutoff, wfs.state.slopes_units, n_sub, offset,
        size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=offset)
    return wfs.state.slopes
end

@inline function zero_invalid_sh_slopes!(::ScalarCPUStyle, slopes::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if !valid_mask[i, j]
            slopes[idx] = zero(T)
            slopes[idx + offset] = zero(T)
        end
        idx += 1
    end
    return slopes
end

@inline function zero_invalid_sh_slopes!(style::AcceleratorStyle, slopes::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    launch_kernel!(style, zero_invalid_sh_slopes_kernel!, slopes, valid_mask, n_sub, offset; ndrange=offset)
    return slopes
end

function sh_signal_from_spots!(wfs::ShackHartmannWFS, peak::T, threshold::T) where {T<:AbstractFloat}
    cutoff = peak <= 0 ? zero(T) : threshold * peak
    return sh_signal_from_spots!(wfs, cutoff)
end

function sh_signal_from_spots!(wfs::ShackHartmannWFS, peak::T, extraction::CenterOfGravityExtraction{T}) where {T<:AbstractFloat}
    return sh_signal_from_spots!(wfs, centroid_cutoff(extraction, peak))
end

function sh_signal_from_spots_calibrated!(wfs::ShackHartmannWFS, peak::T, threshold::T) where {T<:AbstractFloat}
    cutoff = peak <= 0 ? zero(T) : threshold * peak
    return sh_signal_from_spots_calibrated!(wfs, cutoff)
end

function sh_signal_from_spots_calibrated!(wfs::ShackHartmannWFS, peak::T, extraction::CenterOfGravityExtraction{T}) where {T<:AbstractFloat}
    return sh_signal_from_spots_calibrated!(wfs, centroid_cutoff(extraction, peak))
end

function mean_valid_signal(signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return mean_valid_signal(execution_style(signal), signal, valid_mask)
end

function mean_valid_signal(::ScalarCPUStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return packed_valid_pair_mean(ScalarCPUStyle(), signal, valid_mask)
end

function mean_valid_signal(::AcceleratorStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return packed_valid_pair_mean(execution_style(signal), signal, valid_mask)
end

@inline function centroid_from_spot_cube_cutoff!(spot_cube::AbstractArray{T,3}, idx::Int, cutoff::T) where {T<:AbstractFloat}
    total = zero(T)
    sx = zero(T)
    sy = zero(T)
    n1 = size(spot_cube, 2)
    n2 = size(spot_cube, 3)
    @inbounds for x in 1:n1, y in 1:n2
        val = spot_cube[idx, x, y]
        if val < cutoff
            spot_cube[idx, x, y] = zero(T)
        else
            total += val
            sx += T(x - 1) * val
            sy += T(y - 1) * val
        end
    end
    if total <= 0
        return zero(T), zero(T), zero(T)
    end
    return total, sx / total, sy / total
end

@inline function zero_invalid_sh_spot_cube!(::ScalarCPUStyle, spot_cube::AbstractArray{T,3}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if !valid_mask[i, j]
            for x in axes(spot_cube, 2), y in axes(spot_cube, 3)
                spot_cube[idx, x, y] = zero(T)
            end
        end
        idx += 1
    end
    return spot_cube
end

@inline function zero_invalid_sh_spot_cube!(style::AcceleratorStyle, spot_cube::AbstractArray{T,3}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    launch_kernel!(style, zero_invalid_spots_kernel!, spot_cube, valid_mask,
        n_sub, size(spot_cube, 2), size(spot_cube, 3); ndrange=size(spot_cube))
    return spot_cube
end

function subtract_reference_and_scale!(wfs::ShackHartmannWFS)
    return subtract_reference_and_scale!(execution_style(wfs.state.slopes), wfs)
end

function subtract_reference_and_scale!(::ScalarCPUStyle, wfs::ShackHartmannWFS)
    inv_units = inv(wfs.state.slopes_units)
    ref = vec(wfs.state.reference_signal_2d)
    @inbounds for idx in eachindex(wfs.state.slopes, ref)
        wfs.state.slopes[idx] = (wfs.state.slopes[idx] - ref[idx]) * inv_units
    end
    return wfs.state.slopes
end

function subtract_reference_and_scale!(::AcceleratorStyle, wfs::ShackHartmannWFS)
    inv_units = inv(wfs.state.slopes_units)
    ref = vec(wfs.state.reference_signal_2d)
    @. wfs.state.slopes = (wfs.state.slopes - ref) * inv_units
    return wfs.state.slopes
end

function subtract_reference!(wfs::ShackHartmannWFS)
    ref = vec(wfs.state.reference_signal_2d)
    @. wfs.state.slopes = wfs.state.slopes - ref
    return wfs.state.slopes
end

@kernel function calibration_ramp_kernel!(opd, scale, step, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        xi = -π + (i - 1) * step
        xj = -π + (j - 1) * step
        @inbounds opd[i, j] = (xj + xi) * scale
    end
end

function fill_calibration_ramp!(::ScalarCPUStyle, opd::AbstractMatrix{T}, scale::T, n::Int) where {T<:AbstractFloat}
    step = T(2π) / n
    xvals = range(-T(π); step=step, length=n)
    @inbounds for j in 1:n, i in 1:n
        opd[i, j] = (xvals[j] + xvals[i]) * scale
    end
    return opd
end

function fill_calibration_ramp!(style::AcceleratorStyle, opd::AbstractMatrix{T}, scale::T, n::Int) where {T<:AbstractFloat}
    step = T(2π) / n
    launch_kernel!(style, calibration_ramp_kernel!, opd, scale, step, n; ndrange=size(opd))
    return opd
end

function centroid_sums!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs, wfs.state.intensity)
    copyto!(sh_spot_view(wfs, idx), spot)
    return centroid_from_spot!(wfs, spot)
end

function centroid_sums!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
    spot = sample_spot!(wfs, wfs.state.intensity)
    copyto!(sh_spot_view(wfs, idx), spot)
    return centroid_from_spot!(wfs, spot)
end

function centroid_sums!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, idx::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs, wfs.state.intensity)
    frame = capture!(det, spot; rng=rng)
    copyto!(sh_spot_view(wfs, idx), frame)
    return centroid_from_spot!(wfs, frame)
end

function centroid_sums!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, idx::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
    spot = sample_spot!(wfs, wfs.state.intensity)
    frame = capture!(det, spot; rng=rng)
    copyto!(sh_spot_view(wfs, idx), frame)
    return centroid_from_spot!(wfs, frame)
end

function sh_reference_signal!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    peak = sampled_spots_peak!(wfs, tel, src)
    sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
    set_reference_signal!(subaperture_calibration(wfs), reshape(wfs.state.slopes, size(wfs.state.reference_signal_2d)))
    return wfs
end

function ensure_sh_calibration!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    λ = calibration_wavelength(src, eltype(wfs.state.slopes))
    sig = calibration_signature(src)
    if calibration_matches(wfs.state.calibrated, wfs.state.calibration_wavelength, λ,
        wfs.state.calibration_signature, sig)
        return wfs
    end
    opd_saved = save_zero_opd!(tel)
    sh_reference_signal!(wfs, tel, src)
    restore_opd!(tel, opd_saved)

    T = eltype(wfs.state.slopes)
    n = tel.params.resolution
    pixel_scale_init = sh_pixel_scale_init(tel.params.diameter / wfs.params.n_lenslets, wfs.state.effective_padding, src)
    pixel_scale = T(wfs.state.binning_pixel_scale) * pixel_scale_init
    rad2arcsec = T(180 * 3600 / π)
    scale = T(T(tel.params.diameter) * pixel_scale / (T(2π) * rad2arcsec))
    fill_calibration_ramp!(execution_style(tel.state.opd), tel.state.opd, scale, n)
    peak = sampled_spots_peak!(wfs, tel, src)
    sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
    subtract_reference!(wfs)
    wfs.state.slopes_units = mean_valid_signal(wfs.state.slopes, wfs.state.valid_mask)
    if !isfinite(wfs.state.slopes_units) || wfs.state.slopes_units == zero(T)
        wfs.state.slopes_units = one(T)
    end
    restore_opd!(tel, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    wfs.state.calibration_signature = sig
    set_calibration_state!(subaperture_calibration(wfs);
        slopes_units=wfs.state.slopes_units,
        calibrated=wfs.state.calibrated,
        wavelength=wfs.state.calibration_wavelength,
        signature=wfs.state.calibration_signature)
    return wfs
end
