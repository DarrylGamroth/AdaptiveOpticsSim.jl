function sample_spot!(wfs::ShackHartmannWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.state.binning_pixel_scale
    spot_in = intensity
    if binning > 1
        pad = size(intensity, 1)
        if pad % binning != 0
            throw(InvalidConfiguration("lenslet sampling is not divisible by binning_pixel_scale"))
        end
        n_binned = div(pad, binning)
        if size(wfs.state.bin_buffer) != (n_binned, n_binned)
            wfs.state.bin_buffer = similar(wfs.state.bin_buffer, n_binned, n_binned)
        end
        bin2d!(wfs.state.bin_buffer, intensity, binning)
        spot_in = wfs.state.bin_buffer
    end
    center_resize2d!(wfs.state.spot, spot_in)
    return wfs.state.spot
end

function measure!(mode::Geometric, wfs::ShackHartmannWFS, tel::Telescope)
    geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.n_lenslets
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive ShackHartmannWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    Base.require_one_based_indexing(tel.state.opd)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource)
    Base.require_one_based_indexing(tel.state.opd)
    ref = spectral_reference_source(src)
    prepare_sampling!(wfs, tel, ref)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.state.slopes
end

"""
    measure!(Diffractive(), wfs, tel, src)
    measure!(Diffractive(), wfs, tel, src, det; rng)

Compute diffractive Shack-Hartmann slopes from focal-plane lenslet spots.

The algorithm is:

1. propagate each subaperture field with FFT-based Fraunhofer propagation
2. sample/crop the resulting lenslet spots
3. optionally apply detector physics
4. convert spot centroids to x/y slopes
5. subtract the calibrated reference signal and apply the stored slope scaling
"""
function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src, det, rng)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    ref = spectral_reference_source(src)
    prepare_sampling!(wfs, tel, ref)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src, det, rng)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    n = tel.params.resolution
    n_sub = wfs.params.n_lenslets
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.state.slopes), wfs, tel, ast)
        sync_exported_spots!(wfs)
        sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
        subtract_reference_and_scale!(wfs)
        return wfs.state.slopes
    end
    style = execution_style(wfs.state.slopes)
    if sh_uses_accelerator_batched_sensing(style, wfs)
        return measure_sh_asterism_batched!(style, wfs, tel, ast)
    end
    return measure_sh_asterism_diffractive!(style, wfs, tel, ast, n, n_sub, sub, pad, ox, oy)
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    n = tel.params.resolution
    n_sub = wfs.params.n_lenslets
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.state.slopes), wfs, tel, ast, det, rng)
        sync_exported_spots!(wfs)
        sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
        subtract_reference_and_scale!(wfs)
        return wfs.state.slopes
    end
    style = execution_style(wfs.state.slopes)
    if sh_uses_accelerator_batched_sensing(style, wfs)
        return measure_sh_asterism_batched!(style, wfs, tel, ast, det, rng)
    end
    return measure_sh_asterism_diffractive!(style, wfs, tel, ast, det, rng, n, n_sub, sub, pad, ox, oy)
end

@inline sh_uses_accelerator_batched_sensing(::ScalarCPUStyle, ::ShackHartmannWFS) = false
@inline sh_uses_accelerator_batched_sensing(::AcceleratorStyle, wfs::ShackHartmannWFS) =
    sh_uses_batched_sensing_plan(wfs)

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            total_sum = zero(eltype(wfs.state.slopes))
            sx_sum = zero(eltype(wfs.state.slopes))
            sy_sum = zero(eltype(wfs.state.slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
                @views wfs.state.spot_cube_accum[idx, :, :] .+= wfs.state.spot_cube[idx, :, :]
            end
            if total_sum <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                wfs.state.slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_2d[idx]) / wfs.state.slopes_units
                wfs.state.slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_2d[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    zero_invalid_sh_slopes!(ScalarCPUStyle(), wfs.state.slopes, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    sh_refresh_valid_mask_host!(wfs)
    host_slopes = wfs.state.slopes_host
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask_host[i, j]
            total_sum = zero(eltype(host_slopes))
            sx_sum = zero(eltype(host_slopes))
            sy_sum = zero(eltype(host_slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
                @views wfs.state.spot_cube_accum[idx, :, :] .+= wfs.state.spot_cube[idx, :, :]
            end
            if total_sum <= 0
                host_slopes[idx] = zero(eltype(host_slopes))
                host_slopes[idx + n_sub * n_sub] = zero(eltype(host_slopes))
            else
                nsrc = length(ast.sources)
                host_slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_host[idx]) / wfs.state.slopes_units
                host_slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_host[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            host_slopes[idx] = zero(eltype(host_slopes))
            host_slopes[idx + n_sub * n_sub] = zero(eltype(host_slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    copyto!(wfs.state.slopes, host_slopes)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            total_sum = zero(eltype(wfs.state.slopes))
            sx_sum = zero(eltype(wfs.state.slopes))
            sy_sum = zero(eltype(wfs.state.slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx, det, rng)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
                @views wfs.state.spot_cube_accum[idx, :, :] .+= wfs.state.spot_cube[idx, :, :]
            end
            if total_sum <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                wfs.state.slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_2d[idx]) / wfs.state.slopes_units
                wfs.state.slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_2d[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    zero_invalid_sh_slopes!(ScalarCPUStyle(), wfs.state.slopes, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    sh_refresh_valid_mask_host!(wfs)
    host_slopes = wfs.state.slopes_host
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask_host[i, j]
            total_sum = zero(eltype(host_slopes))
            sx_sum = zero(eltype(host_slopes))
            sy_sum = zero(eltype(host_slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx, det, rng)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
                @views wfs.state.spot_cube_accum[idx, :, :] .+= wfs.state.spot_cube[idx, :, :]
            end
            if total_sum <= 0
                host_slopes[idx] = zero(eltype(host_slopes))
                host_slopes[idx + n_sub * n_sub] = zero(eltype(host_slopes))
            else
                nsrc = length(ast.sources)
                host_slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_host[idx]) / wfs.state.slopes_units
                host_slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_host[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            host_slopes[idx] = zero(eltype(host_slopes))
            host_slopes[idx + n_sub * n_sub] = zero(eltype(host_slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    copyto!(wfs.state.slopes, host_slopes)
    return wfs.state.slopes
end
