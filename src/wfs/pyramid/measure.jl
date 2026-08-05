function measure!(mode::Geometric, wfs::PyramidWFS, pupil::PupilFunction)
    slopes = pyramid_estimator_products(wfs).slopes
    geometric_slopes!(slopes, pupil.opd, wfs.estimator.state.valid_mask)
    gain = inv(1 + wfs.estimator.params.geometric_modulation_radius)
    @. slopes = gain * slopes * wfs.estimator.state.optical_gain
    return slopes
end

@inline function apply_pyramid_optical_gain!(wfs::PyramidWFS)
    slopes = pyramid_estimator_products(wfs).slopes
    @. slopes *= wfs.estimator.state.optical_gain
    return slopes
end

function measure!(::Geometric, wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(Geometric(), wfs, pupil)
end

function measure!(::Geometric, wfs::PyramidWFS, pupil::PupilFunction, src::LGSSource)
    slopes = measure!(Geometric(), wfs, pupil)
    n_sub = wfs.estimator.params.pupil_samples
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction)
    throw(InvalidConfiguration("Diffractive PyramidWFS requires a source; call measure!(wfs, pupil, src)."))
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction)
    return measure!(sensing_mode(wfs), wfs, pupil)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    ensure_pyramid_calibration!(wfs, pupil, src)
    propagation = pyramid_propagation_workspace(wfs)
    pyramid_intensity!(propagation.intensity, wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, src)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource)
    ensure_pyramid_calibration!(wfs, pupil, src)
    propagation = pyramid_propagation_workspace(wfs)
    accumulate_pyramid_spectral_intensity!(
        execution_style(propagation.intensity), wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, src)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    propagation = pyramid_propagation_workspace(wfs)
    pyramid_intensity!(propagation.intensity, wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    propagation = pyramid_propagation_workspace(wfs)
    accumulate_pyramid_spectral_intensity!(execution_style(propagation.intensity),
        wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction,
    src::SpectralSource, det::Detector;
    rng::AbstractRNG=Random.default_rng())
    qe_model = quantum_efficiency_model(det)
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    propagation = pyramid_propagation_workspace(wfs)
    accumulate_pyramid_spectral_intensity!(execution_style(propagation.intensity),
        wfs, pupil, src, qe_model)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    frame = capture_with_quantum_efficiency!(det, intensity,
        one(eltype(intensity)), rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    ensure_pyramid_calibration!(wfs, pupil, common_source)
    propagation = pyramid_propagation_workspace(wfs)
    accumulate_pyramid_asterism_intensity!(
        execution_style(propagation.intensity), wfs, pupil, ast)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, ast)
    return apply_pyramid_optical_gain!(wfs)
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    ensure_pyramid_calibration!(wfs, pupil, common_source, det)
    propagation = pyramid_propagation_workspace(wfs)
    accumulate_pyramid_asterism_intensity!(
        execution_style(propagation.intensity), wfs, pupil, ast)
    intensity = sample_pyramid_intensity!(wfs, pupil, propagation.intensity)
    frame = capture!(det, intensity, common_source; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, common_source,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, ast, normalization_scale)
    return apply_pyramid_optical_gain!(wfs)
end
