
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
    @inbounds for y in 1:n2, x in 1:n1
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
    result = centroid_from_intensity!(ScalarCPUStyle(), host_intensity, threshold)
    copyto!(intensity, host_intensity)
    return result
end

@inline function centroid_from_intensity_cutoff!(::AcceleratorStyle, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    host_intensity = Array(intensity)
    result = centroid_from_intensity_cutoff!(ScalarCPUStyle(), host_intensity, cutoff)
    copyto!(intensity, host_intensity)
    return result
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
    if size(wfs.estimator.centroid_host) != size(intensity)
        wfs.estimator.centroid_host = Matrix{T}(undef, size(intensity)...)
    end
    copyto!(wfs.estimator.centroid_host, intensity)
    result = centroid_from_intensity!(ScalarCPUStyle(), wfs.estimator.centroid_host, threshold)
    copyto!(intensity, wfs.estimator.centroid_host)
    return result
end

@inline function centroid_from_spot_cutoff!(wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_spot_cutoff!(execution_style(intensity), wfs, intensity, cutoff)
end

@inline function centroid_from_spot_cutoff!(::ScalarCPUStyle, ::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), intensity, cutoff)
end

@inline function centroid_from_spot_cutoff!(::AcceleratorStyle, wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    if size(wfs.estimator.centroid_host) != size(intensity)
        wfs.estimator.centroid_host = Matrix{T}(undef, size(intensity)...)
    end
    copyto!(wfs.estimator.centroid_host, intensity)
    result = centroid_from_intensity_cutoff!(ScalarCPUStyle(), wfs.estimator.centroid_host, cutoff)
    copyto!(intensity, wfs.estimator.centroid_host)
    return result
end

function sh_signal_from_spots_device_stats!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = n_lenslets(wfs)
    offset = n_sub * n_sub
    fill!(wfs.estimator.spot_stats, zero(eltype(wfs.estimator.spot_stats)))
    launch_kernel!(style, sh_spot_cutoff_stats_kernel!, wfs.estimator.spot_stats, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask, cutoff, n_sub, size(wfs.acquisition.spot_cube, 2),
        size(wfs.acquisition.spot_cube, 3); ndrange=(n_sub, n_sub))
    launch_kernel!(style, sh_finalize_spot_slopes_kernel!, wfs.estimator.slopes, wfs.estimator.spot_stats,
        wfs.front_end.layout.valid_mask, n_sub, offset; ndrange=(n_sub, n_sub))
    return wfs.estimator.slopes
end

function sh_signal_from_spots_calibrated_device_stats!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = n_lenslets(wfs)
    offset = n_sub * n_sub
    reference = vec(wfs.calibration.reference_signal_2d)
    fill!(wfs.estimator.spot_stats, zero(eltype(wfs.estimator.spot_stats)))
    launch_kernel!(style, sh_spot_cutoff_stats_kernel!, wfs.estimator.spot_stats, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask, cutoff, n_sub, size(wfs.acquisition.spot_cube, 2),
        size(wfs.acquisition.spot_cube, 3); ndrange=(n_sub, n_sub))
    launch_kernel!(style, sh_finalize_spot_slopes_reference_scale_kernel!, wfs.estimator.slopes, wfs.estimator.spot_stats,
        reference, wfs.front_end.layout.valid_mask, wfs.calibration.centroid_response, n_sub, offset;
        ndrange=(n_sub, n_sub))
    return wfs.estimator.slopes
end

@inline function sh_spot_view(wfs::ShackHartmannWFS, idx::Int)
    return @view wfs.acquisition.spot_cube[idx, :, :]
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

@inline function sh_spectral_source_variant(wfs::ShackHartmannWFS,
    src::SpectralSource, sample::SpectralSample, radiometric_value::Real)
    T = eltype(wfs.estimator.slopes)
    return source_with_wavelength_and_radiometric_value(src,
        T(sample.wavelength), T(radiometric_value))
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource)
    ast = extended_source_asterism(src)
    prepare_sampling!(wfs, pupil, ast.sources[1])
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    compute_intensity_stack!(ScalarCPUStyle(), wfs, pupil, src)
    sample_spot_stack!(ScalarCPUStyle(), wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n = _pupil_resolution(pupil)
        n_sub = n_lenslets(wfs)
        sub = div(n, n_sub)
        pad = size(wfs.front_end.propagation.field, 1)
        ox = div(pad - sub, 2)
        oy = div(pad - sub, 2)
        peak = zero(eltype(wfs.estimator.slopes))
        idx = 1
        @inbounds for j in 1:n_sub, i in 1:n_sub
            spot_view = sh_spot_view(wfs, idx)
            if wfs.front_end.layout.valid_mask_host[i, j]
                xs = (i - 1) * sub + 1
                ys = (j - 1) * sub + 1
                xe = min(i * sub, n)
                ye = min(j * sub, n)
                compute_intensity_safe!(style, wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
                sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
                sync_sh_staged_spot!(style, wfs.front_end.propagation.spot)
                copyto!(spot_view, wfs.front_end.propagation.spot)
                sync_sh_staged_view!(style, spot_view)
                peak = max(peak, sh_safe_peak_value(spot_view))
            else
                fill!(spot_view, zero(eltype(spot_view)))
            end
            idx += 1
        end
        return peak
    end
    compute_intensity_stack!(style, wfs, pupil, src)
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource)
    require_sh_common_spectral_grid(wfs, src)
    fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
    peak = zero(eltype(wfs.estimator.slopes))
    total_irradiance = photon_irradiance(src)
    @inbounds for sample in src.bundle.samples
        variant = sh_spectral_source_variant(wfs, src, sample,
            total_irradiance * sample.weight)
        peak = max(peak, sampled_spots_peak!(ScalarCPUStyle(), wfs, pupil, variant))
        wfs.front_end.propagation.spot_cube_accum .+= wfs.acquisition.spot_cube
    end
    copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
    return max(peak, sh_safe_peak_value(wfs.acquisition.spot_cube))
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource)
    require_sh_common_spectral_grid(wfs, src)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
        peak = zero(eltype(wfs.estimator.slopes))
        total_irradiance = photon_irradiance(src)
        @inbounds for sample in src.bundle.samples
            variant = sh_spectral_source_variant(wfs, src, sample,
                total_irradiance * sample.weight)
            peak = max(peak, sampled_spots_peak!(style, wfs, pupil, variant))
            @. wfs.front_end.propagation.spot_cube_accum = wfs.front_end.propagation.spot_cube_accum + wfs.acquisition.spot_cube
        end
        copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
        return max(peak, sh_safe_peak_value(wfs.acquisition.spot_cube))
    end
    if is_lgs_source(src)
        fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
        peak = zero(eltype(wfs.estimator.slopes))
        total_irradiance = photon_irradiance(src)
        @inbounds for sample in src.bundle.samples
            variant = sh_spectral_source_variant(wfs, src, sample,
                total_irradiance * sample.weight)
            peak = max(peak, sampled_spots_peak!(style, wfs, pupil, variant))
            @. wfs.front_end.propagation.spot_cube_accum = wfs.front_end.propagation.spot_cube_accum + wfs.acquisition.spot_cube
        end
        copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
        return max(peak, sh_safe_peak_value(wfs.acquisition.spot_cube))
    end
    compute_intensity_spectral_stack!(style, wfs, pupil, src)
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(ScalarCPUStyle(), wfs, pupil, ast.sources[1])
    end
    return sampled_spots_peak_asterism_stacked!(ScalarCPUStyle(), wfs, pupil, ast)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(style, wfs, pupil, ast.sources[1])
    end
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        return sampled_spots_peak_asterism_stacked!(style, wfs, pupil, ast)
    end
    fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
    peak = zero(eltype(wfs.estimator.slopes))
    @inbounds for sample in ast.sources
        peak = max(peak, sampled_spots_peak!(style, wfs, pupil, sample))
        @. wfs.front_end.propagation.spot_cube_accum = wfs.front_end.propagation.spot_cube_accum + wfs.acquisition.spot_cube
    end
    copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
    return max(peak, sh_safe_peak_value(wfs.acquisition.spot_cube))
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    n = _pupil_resolution(pupil)
    n_sub = n_lenslets(wfs)
    sub = div(n, n_sub)
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.estimator.slopes))
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.front_end.layout.valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, pupil, src, idx)
            sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
            copyto!(spot_view, wfs.front_end.propagation.spot)
            peak = max(peak, sh_safe_peak_value(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, pupil, src)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src, det, rng)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src, det, rng)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    prepare_sampling!(wfs, pupil, ast.sources[1])
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(ScalarCPUStyle(), wfs, pupil, src)
    sample_spot_stack!(ScalarCPUStyle(), wfs.front_end)
    capture_sampled_spot_stack!(wfs, src, det, rng)
    zero_invalid_sh_spot_cube!(ScalarCPUStyle(), wfs.acquisition.spot_cube, wfs.front_end.layout.valid_mask)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n = _pupil_resolution(pupil)
        n_sub = n_lenslets(wfs)
        sub = div(n, n_sub)
        pad = size(wfs.front_end.propagation.field, 1)
        ox = div(pad - sub, 2)
        oy = div(pad - sub, 2)
        peak = zero(eltype(wfs.estimator.slopes))
        idx = 1
        @inbounds for j in 1:n_sub, i in 1:n_sub
            spot_view = sh_spot_view(wfs, idx)
            if wfs.front_end.layout.valid_mask_host[i, j]
                xs = (i - 1) * sub + 1
                ys = (j - 1) * sub + 1
                xe = min(i * sub, n)
                ye = min(j * sub, n)
                compute_intensity_safe!(style, wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
                sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
                sync_sh_staged_spot!(style, wfs.front_end.propagation.spot)
                frame = capture!(det, wfs.front_end.propagation.spot, src; rng=rng)
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
    compute_intensity_stack!(style, wfs, pupil, src)
    sample_spot_stack!(style, wfs.front_end)
    n_sub = n_lenslets(wfs)
    capture_sampled_spot_stack!(wfs, src, det, rng)
    n1, n2 = size(wfs.acquisition.spot_cube, 2), size(wfs.acquisition.spot_cube, 3)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.acquisition.spot_cube, wfs.front_end.layout.valid_mask,
        n_sub, n1, n2; ndrange=(n_sub, n_sub, n1, n2))
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    require_sh_common_spectral_grid(wfs, src)
    fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
    total_irradiance = photon_irradiance(src)
    qe_model = quantum_efficiency_model(det)
    @inbounds for sample in src.bundle.samples
        channel_qe = qe_at(qe_model, sample.wavelength)
        variant = sh_spectral_source_variant(wfs, src, sample,
            total_irradiance * sample.weight * channel_qe)
        sampled_spots_peak!(ScalarCPUStyle(), wfs, pupil, variant)
        wfs.front_end.propagation.spot_cube_accum .+= wfs.acquisition.spot_cube
    end
    copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
    return capture_sh_qe_weighted_spots!(ScalarCPUStyle(), wfs, det, rng)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector, rng::AbstractRNG)
    require_sh_common_spectral_grid(wfs, src)
    if sh_uses_rocm_safe_sensing_plan(wfs)
        fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
        total_irradiance = photon_irradiance(src)
        qe_model = quantum_efficiency_model(det)
        @inbounds for sample in src.bundle.samples
            channel_qe = qe_at(qe_model, sample.wavelength)
            variant = sh_spectral_source_variant(wfs, src, sample,
                total_irradiance * sample.weight * channel_qe)
            sampled_spots_peak!(style, wfs, pupil, variant)
            @. wfs.front_end.propagation.spot_cube_accum = wfs.front_end.propagation.spot_cube_accum + wfs.acquisition.spot_cube
        end
        copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
        return capture_sh_qe_weighted_spots!(style, wfs, det, rng)
    end
    if is_lgs_source(src)
        fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
        total_irradiance = photon_irradiance(src)
        qe_model = quantum_efficiency_model(det)
        @inbounds for sample in src.bundle.samples
            channel_qe = qe_at(qe_model, sample.wavelength)
            variant = sh_spectral_source_variant(wfs, src, sample,
                total_irradiance * sample.weight * channel_qe)
            sampled_spots_peak!(style, wfs, pupil, variant)
            @. wfs.front_end.propagation.spot_cube_accum = wfs.front_end.propagation.spot_cube_accum + wfs.acquisition.spot_cube
        end
        copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
        return capture_sh_qe_weighted_spots!(style, wfs, det, rng)
    end
    compute_intensity_spectral_stack!(style, wfs, pupil, src,
        quantum_efficiency_model(det))
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return capture_sh_qe_weighted_spots!(style, wfs, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(ScalarCPUStyle(), wfs, pupil, ast.sources[1], det, rng)
    end
    return sampled_spots_peak_asterism_stacked!(ScalarCPUStyle(), wfs, pupil, ast, det, rng)
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::ExtendedSource,
    det::AbstractDetector, rng::AbstractRNG)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return sampled_spots_peak!(style, wfs, pupil, ast.sources[1], det, rng)
    end
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        return sampled_spots_peak_asterism_stacked!(style, wfs, pupil, ast, det, rng)
    end
    accumulate_sh_asterism_spots!(style, wfs, pupil, ast)
    return capture_sh_asterism_spots!(style, wfs, ast, det, rng)
end

function sampled_spots_peak!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    prepare_sampling!(wfs, pupil, src)
    return sampled_spots_peak!(execution_style(wfs.front_end.layout.valid_mask), wfs, pupil, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = _pupil_resolution(pupil)
    n_sub = n_lenslets(wfs)
    sub = div(n, n_sub)
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.estimator.slopes))
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.front_end.layout.valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, pupil, src, idx)
            sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
            frame = capture!(det, wfs.front_end.propagation.spot, src; rng=rng)
            copyto!(spot_view, frame)
            peak = max(peak, sh_safe_peak_value(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, pupil, src, det, rng)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    compute_intensity_stack!(style, wfs, pupil, src)
    tmp_view = @view wfs.front_end.propagation.intensity_tmp_stack[:, :, 1:size(wfs.front_end.propagation.intensity_stack, 3)]
    apply_elongation_stack!(wfs.front_end.propagation.intensity_stack, lgs_elongation_factor(src),
        tmp_view, wfs.front_end.propagation.elongation_kernel)
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(style, wfs, pupil, src)
    tmp_view = @view wfs.front_end.propagation.intensity_tmp_stack[:, :, 1:size(wfs.front_end.propagation.intensity_stack, 3)]
    apply_elongation_stack!(wfs.front_end.propagation.intensity_stack, lgs_elongation_factor(src),
        tmp_view, wfs.front_end.propagation.elongation_kernel)
    sample_spot_stack!(style, wfs.front_end)
    n_sub = n_lenslets(wfs)
    capture_sampled_spot_stack!(wfs, src, det, rng)
    n1, n2 = size(wfs.acquisition.spot_cube, 2), size(wfs.acquisition.spot_cube, 3)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.acquisition.spot_cube, wfs.front_end.layout.valid_mask,
        n_sub, n1, n2; ndrange=(n_sub, n_sub, n1, n2))
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    n_sub = n_lenslets(wfs)
    compute_intensity_stack!(style, wfs, pupil, src)
    ensure_lgs_kernels!(wfs, pupil, src)
    apply_lgs_convolution_stack!(wfs.front_end.propagation.intensity_stack, wfs.front_end.propagation.lgs_kernel_fft,
        wfs.front_end.propagation.fft_stack, wfs.front_end.propagation.fft_stack_plan, wfs.front_end.propagation.ifft_stack_plan)
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n_sub = n_lenslets(wfs)
    compute_intensity_stack!(style, wfs, pupil, src)
    ensure_lgs_kernels!(wfs, pupil, src)
    apply_lgs_convolution_stack!(wfs.front_end.propagation.intensity_stack, wfs.front_end.propagation.lgs_kernel_fft,
        wfs.front_end.propagation.fft_stack, wfs.front_end.propagation.fft_stack_plan, wfs.front_end.propagation.ifft_stack_plan)
    sample_spot_stack!(style, wfs.front_end)
    capture_sampled_spot_stack!(wfs, src, det, rng)
    n1, n2 = size(wfs.acquisition.spot_cube, 2), size(wfs.acquisition.spot_cube, 3)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.acquisition.spot_cube, wfs.front_end.layout.valid_mask,
        n_sub, n1, n2; ndrange=(n_sub, n_sub, n1, n2))
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sh_signal_from_spots!(wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    return sh_signal_from_spots!(execution_style(wfs.estimator.slopes), wfs, cutoff)
end

function sh_signal_from_spots!(::ScalarCPUStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    n_sub = n_lenslets(wfs)
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if wfs.front_end.layout.valid_mask[i, j]
            total, sx, sy = centroid_from_spot_cube_cutoff!(wfs.acquisition.spot_cube, idx, cutoff)
            if total <= 0
                wfs.estimator.slopes[idx] = zero(T)
                wfs.estimator.slopes[idx + n_sub * n_sub] = zero(T)
            else
                wfs.estimator.slopes[idx] = sx
                wfs.estimator.slopes[idx + n_sub * n_sub] = sy
            end
        else
            wfs.estimator.slopes[idx] = zero(T)
            wfs.estimator.slopes[idx + n_sub * n_sub] = zero(T)
        end
        idx += 1
    end
    return wfs.estimator.slopes
end

function sh_signal_from_spots!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    if sh_uses_host_stats_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n_sub = n_lenslets(wfs)
        offset = n_sub * n_sub
        host_slopes = wfs.estimator.slopes_host
        idx = 1
        @inbounds for j in 1:n_sub, i in 1:n_sub
            if wfs.front_end.layout.valid_mask_host[i, j]
                total, sx, sy = centroid_from_spot_cutoff!(wfs, sh_spot_view(wfs, idx), cutoff)
                if total <= 0
                    host_slopes[idx] = zero(T)
                    host_slopes[idx + offset] = zero(T)
                else
                    host_slopes[idx] = sx
                    host_slopes[idx + offset] = sy
                end
            else
                host_slopes[idx] = zero(T)
                host_slopes[idx + offset] = zero(T)
            end
            idx += 1
        end
        copyto!(wfs.estimator.slopes, host_slopes)
        return wfs.estimator.slopes
    end
    if sh_uses_device_stats_sensing_plan(wfs)
        return sh_signal_from_spots_device_stats!(style, wfs, cutoff)
    end
    n_sub = n_lenslets(wfs)
    offset = n_sub * n_sub
    launch_kernel!(style, sh_spot_centroid_kernel!, wfs.estimator.slopes, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask, cutoff, n_sub, offset, size(wfs.acquisition.spot_cube, 2),
        size(wfs.acquisition.spot_cube, 3); ndrange=(n_sub, n_sub))
    return wfs.estimator.slopes
end

function sh_signal_from_spots_calibrated!(wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    return sh_signal_from_spots_calibrated!(execution_style(wfs.estimator.slopes), wfs, cutoff)
end

function sh_signal_from_spots_calibrated!(::ScalarCPUStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    sh_signal_from_spots!(ScalarCPUStyle(), wfs, cutoff)
    subtract_reference_and_scale!(ScalarCPUStyle(), wfs)
    return wfs.estimator.slopes
end

function sh_signal_from_spots_calibrated!(style::AcceleratorStyle, wfs::ShackHartmannWFS, cutoff::T) where {T<:AbstractFloat}
    if sh_uses_host_stats_sensing_plan(wfs)
        sh_refresh_valid_mask_host!(wfs)
        n_sub = n_lenslets(wfs)
        offset = n_sub * n_sub
        host_slopes = wfs.estimator.slopes_host
        reference = wfs.calibration.reference_signal_host
        inv_units = inv(wfs.calibration.centroid_response)
        idx = 1
        @inbounds for j in 1:n_sub, i in 1:n_sub
            if wfs.front_end.layout.valid_mask_host[i, j]
                total, sx, sy = centroid_from_spot_cutoff!(wfs, sh_spot_view(wfs, idx), cutoff)
                if total <= 0
                    host_slopes[idx] = -reference[idx] * inv_units
                    host_slopes[idx + offset] = -reference[idx + offset] * inv_units
                else
                    host_slopes[idx] = (sx - reference[idx]) * inv_units
                    host_slopes[idx + offset] = (sy - reference[idx + offset]) * inv_units
                end
            else
                host_slopes[idx] = zero(T)
                host_slopes[idx + offset] = zero(T)
            end
            idx += 1
        end
        copyto!(wfs.estimator.slopes, host_slopes)
        return wfs.estimator.slopes
    end
    if sh_uses_device_stats_sensing_plan(wfs)
        return sh_signal_from_spots_calibrated_device_stats!(style, wfs, cutoff)
    end
    n_sub = n_lenslets(wfs)
    offset = n_sub * n_sub
    reference = vec(wfs.calibration.reference_signal_2d)
    launch_kernel!(style, sh_spot_centroid_reference_scale_kernel!, wfs.estimator.slopes, wfs.acquisition.spot_cube,
        reference, wfs.front_end.layout.valid_mask, cutoff, wfs.calibration.centroid_response, n_sub, offset,
        size(wfs.acquisition.spot_cube, 2), size(wfs.acquisition.spot_cube, 3);
        ndrange=(n_sub, n_sub))
    return wfs.estimator.slopes
end

@inline function zero_invalid_sh_slopes!(::ScalarCPUStyle, slopes::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
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
    launch_kernel!(style, zero_invalid_sh_slopes_kernel!, slopes, valid_mask, n_sub, offset;
        ndrange=(n_sub, n_sub))
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
    @inbounds for y in 1:n2, x in 1:n1
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
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if !valid_mask[i, j]
            for y in axes(spot_cube, 3), x in axes(spot_cube, 2)
                spot_cube[idx, x, y] = zero(T)
            end
        end
        idx += 1
    end
    return spot_cube
end

@inline function zero_invalid_sh_spot_cube!(style::AcceleratorStyle, spot_cube::AbstractArray{T,3}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    n1, n2 = size(spot_cube, 2), size(spot_cube, 3)
    launch_kernel!(style, zero_invalid_spots_kernel!, spot_cube, valid_mask,
        n_sub, n1, n2; ndrange=(n_sub, n_sub, n1, n2))
    return spot_cube
end

function subtract_reference_and_scale!(wfs::ShackHartmannWFS)
    return subtract_reference_and_scale!(execution_style(wfs.estimator.slopes), wfs)
end

function subtract_reference_and_scale!(::ScalarCPUStyle, wfs::ShackHartmannWFS)
    inv_units = inv(wfs.calibration.centroid_response)
    ref = vec(wfs.calibration.reference_signal_2d)
    @inbounds for idx in eachindex(wfs.estimator.slopes, ref)
        wfs.estimator.slopes[idx] = (wfs.estimator.slopes[idx] - ref[idx]) * inv_units
    end
    return wfs.estimator.slopes
end

function subtract_reference_and_scale!(::AcceleratorStyle, wfs::ShackHartmannWFS)
    inv_units = inv(wfs.calibration.centroid_response)
    ref = vec(wfs.calibration.reference_signal_2d)
    @. wfs.estimator.slopes = (wfs.estimator.slopes - ref) * inv_units
    return wfs.estimator.slopes
end

function subtract_reference!(wfs::ShackHartmannWFS)
    ref = vec(wfs.calibration.reference_signal_2d)
    @. wfs.estimator.slopes = wfs.estimator.slopes - ref
    return wfs.estimator.slopes
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

function centroid_sums!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
    copyto!(sh_spot_view(wfs, idx), spot)
    return centroid_from_spot!(wfs, spot)
end

function centroid_sums!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, pupil, src, idx)
    spot = sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
    copyto!(sh_spot_view(wfs, idx), spot)
    return centroid_from_spot!(wfs, spot)
end

function centroid_sums!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, idx::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
    frame = capture!(det, spot, src; rng=rng)
    copyto!(sh_spot_view(wfs, idx), frame)
    return centroid_from_spot!(wfs, frame)
end

function centroid_sums!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, idx::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, pupil, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, pupil, src, idx)
    spot = sample_spot!(wfs.front_end, wfs.front_end.propagation.intensity)
    frame = capture!(det, spot, src; rng=rng)
    copyto!(sh_spot_view(wfs, idx), frame)
    return centroid_from_spot!(wfs, frame)
end

function sh_reference_signal!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    peak = sampled_spots_peak!(wfs, pupil, src)
    sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
    set_reference_signal!(subaperture_calibration(wfs), reshape(wfs.estimator.slopes, size(wfs.calibration.reference_signal_2d)))
    return wfs
end

function ensure_sh_calibration!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    λ = sh_calibration_wavelength(wfs, src)
    sig = pupil_aperture_calibration_signature(pupil,
        calibration_signature(src))
    if calibration_matches(wfs.calibration.calibrated, wfs.calibration.wavelength, λ,
        wfs.calibration.signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    T = eltype(wfs.estimator.slopes)
    centroid_response = one(T)
    try
        sh_reference_signal!(wfs, pupil, src)
        restore_opd!(pupil, opd_saved)

        n = _pupil_resolution(pupil)
        pixel_scale_init = sh_pixel_scale_init(
            _pupil_diameter_m(pupil) / n_lenslets(wfs),
            wfs.front_end.propagation.effective_padding, λ)
        pixel_scale = T(wfs.front_end.propagation.binning_pixel_scale) * pixel_scale_init
        rad2arcsec = T(180 * 3600 / π)
        scale = T(T(_pupil_diameter_m(pupil)) * pixel_scale /
            (T(2π) * rad2arcsec))
        fill_calibration_ramp!(execution_style(pupil.opd), pupil.opd,
            scale, n)
        peak = sampled_spots_peak!(wfs, pupil, src)
        sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
        subtract_reference!(wfs)
        centroid_response = mean_valid_signal(wfs.estimator.slopes,
            wfs.front_end.layout.valid_mask)
        if !isfinite(centroid_response) || iszero(centroid_response)
            centroid_response = one(T)
        end
    finally
        restore_opd!(pupil, opd_saved)
    end
    set_subaperture_calibration!(subaperture_calibration(wfs),
        wfs.calibration.reference_signal_2d;
        centroid_response, output_units=:pixel, wavelength=λ,
        signature=sig)
    return wfs
end

@inline sh_calibration_wavelength(wfs::ShackHartmannWFS,
    src::AbstractSource) = calibration_wavelength(src,
    eltype(wfs.estimator.slopes))

@inline sh_calibration_wavelength(wfs::ShackHartmannWFS,
    src::SpectralSource) = require_sh_common_spectral_grid(wfs, src)
