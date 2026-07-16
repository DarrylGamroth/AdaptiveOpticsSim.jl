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
    common_source = common_wfs_calibration_source(ast, "ShackHartmannWFS")
    prepare_sampling!(wfs, tel, common_source)
    ensure_sh_calibration!(wfs, tel, common_source)
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
    common_source = common_wfs_calibration_source(ast, "ShackHartmannWFS")
    prepare_sampling!(wfs, tel, common_source)
    ensure_sh_calibration!(wfs, tel, common_source)
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
    peak = accumulate_sh_asterism_spots!(ScalarCPUStyle(), wfs, tel, ast)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    peak = accumulate_sh_asterism_spots!(style, wfs, tel, ast)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    accumulate_sh_asterism_spots!(ScalarCPUStyle(), wfs, tel, ast)
    peak = capture_sh_asterism_spots!(ScalarCPUStyle(), wfs, ast, det, rng)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    accumulate_sh_asterism_spots!(style, wfs, tel, ast)
    peak = capture_sh_asterism_spots!(style, wfs, ast, det, rng)
    return finalize_sh_asterism_slopes!(wfs, peak)
end
