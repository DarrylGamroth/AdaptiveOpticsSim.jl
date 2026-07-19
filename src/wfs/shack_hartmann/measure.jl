function sample_spot!(front_end::ShackHartmannOpticalFrontEnd,
    intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    propagation = front_end.propagation
    binning = propagation.binning_pixel_scale
    spot_in = intensity
    if binning > 1
        pad = size(intensity, 1)
        if pad % binning != 0
            throw(InvalidConfiguration("lenslet sampling is not divisible by binning_pixel_scale"))
        end
        n_binned = div(pad, binning)
        if size(propagation.bin_buffer) != (n_binned, n_binned)
            propagation.bin_buffer = similar(propagation.bin_buffer,
                n_binned, n_binned)
        end
        bin2d!(propagation.bin_buffer, intensity, binning)
        spot_in = propagation.bin_buffer
    end
    center_resize2d!(propagation.spot, spot_in)
    return propagation.spot
end

function measure!(mode::Geometric, wfs::ShackHartmannWFS, pupil::PupilFunction)
    geometric_slopes!(wfs.estimator.slopes, pupil.opd,
        wfs.front_end.layout.valid_mask)
    return wfs.estimator.slopes
end

function measure!(::Geometric, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(Geometric(), wfs, pupil)
end

function measure!(::Geometric, wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    slopes = measure!(Geometric(), wfs, pupil)
    n_sub = n_lenslets(wfs)
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction)
    throw(InvalidConfiguration("Diffractive ShackHartmannWFS requires a source; call measure!(wfs, pupil, src)."))
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction)
    return measure!(sensing_mode(wfs), wfs, pupil)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::ShackHartmannWFS, pupil::PupilFunction, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource)
    Base.require_one_based_indexing(pupil.opd)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    peak = sampled_spots_peak!(wfs, pupil, src)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.estimator.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource)
    Base.require_one_based_indexing(pupil.opd)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    peak = sampled_spots_peak!(wfs, pupil, src)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.estimator.slopes
end

"""
    measure!(Diffractive(), wfs, pupil, src)
    measure!(Diffractive(), wfs, pupil, src, det; rng)

Compute diffractive Shack-Hartmann slopes from focal-plane lenslet spots.

The algorithm is:

1. propagate each subaperture field with FFT-based Fraunhofer propagation
2. sample/crop the resulting lenslet spots
3. optionally apply detector physics
4. convert spot centroids to x/y slopes
5. subtract the calibrated reference signal and apply the stored slope scaling
"""
function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    peak = sampled_spots_peak!(wfs, pupil, src, det, rng)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.estimator.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    prepare_sampling!(wfs, pupil, src)
    ensure_sh_calibration!(wfs, pupil, src)
    peak = sampled_spots_peak!(wfs, pupil, src, det, rng)
    sync_exported_spots!(wfs)
    sh_signal_from_spots_calibrated!(wfs, peak, slope_extraction_model(wfs))
    return wfs.estimator.slopes
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, ast::Asterism)
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "ShackHartmannWFS")
    prepare_sampling!(wfs, pupil, common_source)
    ensure_sh_calibration!(wfs, pupil, common_source)
    n = _pupil_resolution(pupil)
    n_sub = n_lenslets(wfs)
    sub = div(n, n_sub)
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.estimator.slopes), wfs, pupil, ast)
        sync_exported_spots!(wfs)
        sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
        subtract_reference_and_scale!(wfs)
        return wfs.estimator.slopes
    end
    style = execution_style(wfs.estimator.slopes)
    if sh_uses_accelerator_batched_sensing(style, wfs)
        return measure_sh_asterism_batched!(style, wfs, pupil, ast)
    end
    return measure_sh_asterism_diffractive!(style, wfs, pupil, ast, n, n_sub, sub, pad, ox, oy)
end

function measure!(::Diffractive, wfs::ShackHartmannWFS, pupil::PupilFunction, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "ShackHartmannWFS")
    prepare_sampling!(wfs, pupil, common_source)
    ensure_sh_calibration!(wfs, pupil, common_source)
    n = _pupil_resolution(pupil)
    n_sub = n_lenslets(wfs)
    sub = div(n, n_sub)
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    if sh_stacked_asterism_compatible(ast) && sh_uses_batched_sensing_plan(wfs)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.estimator.slopes), wfs, pupil, ast, det, rng)
        sync_exported_spots!(wfs)
        sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
        subtract_reference_and_scale!(wfs)
        return wfs.estimator.slopes
    end
    style = execution_style(wfs.estimator.slopes)
    if sh_uses_accelerator_batched_sensing(style, wfs)
        return measure_sh_asterism_batched!(style, wfs, pupil, ast, det, rng)
    end
    return measure_sh_asterism_diffractive!(style, wfs, pupil, ast, det, rng, n, n_sub, sub, pad, ox, oy)
end

@inline sh_uses_accelerator_batched_sensing(::ScalarCPUStyle, ::ShackHartmannWFS) = false
@inline sh_uses_accelerator_batched_sensing(::AcceleratorStyle, wfs::ShackHartmannWFS) =
    sh_uses_batched_sensing_plan(wfs)

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    peak = accumulate_sh_asterism_spots!(ScalarCPUStyle(), wfs, pupil, ast)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    peak = accumulate_sh_asterism_spots!(style, wfs, pupil, ast)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmannWFS, pupil::PupilFunction,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    accumulate_sh_asterism_spots!(ScalarCPUStyle(), wfs, pupil, ast)
    peak = capture_sh_asterism_spots!(ScalarCPUStyle(), wfs, ast, det, rng)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_diffractive!(style::AcceleratorStyle, wfs::ShackHartmannWFS, pupil::PupilFunction,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    accumulate_sh_asterism_spots!(style, wfs, pupil, ast)
    peak = capture_sh_asterism_spots!(style, wfs, ast, det, rng)
    return finalize_sh_asterism_slopes!(wfs, peak)
end
