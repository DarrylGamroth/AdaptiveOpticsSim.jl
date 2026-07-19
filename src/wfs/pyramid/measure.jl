function measure!(mode::Geometric, wfs::PyramidWFS, pupil::PupilFunction)
    geometric_slopes!(wfs.estimator.state.slopes, pupil.opd, wfs.estimator.state.valid_mask)
    gain = inv(1 + wfs.estimator.params.geometric_modulation_radius)
    @. wfs.estimator.state.slopes = gain * wfs.estimator.state.slopes * wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
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
    pyramid_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, src)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource)
    ensure_pyramid_calibration!(wfs, pupil, src)
    accumulate_pyramid_spectral_intensity!(execution_style(wfs.front_end.propagation.intensity), wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, src)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    pyramid_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    accumulate_pyramid_spectral_intensity!(execution_style(wfs.front_end.propagation.intensity),
        wfs, pupil, src)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction,
    src::SpectralSource, det::Detector;
    rng::AbstractRNG=Random.default_rng())
    qe_model = quantum_efficiency_model(det)
    ensure_pyramid_calibration!(wfs, pupil, src, det)
    accumulate_pyramid_spectral_intensity!(execution_style(wfs.front_end.propagation.intensity),
        wfs, pupil, src, qe_model)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture_with_quantum_efficiency!(det, intensity,
        one(eltype(intensity)), rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, src, normalization_scale)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    ensure_pyramid_calibration!(wfs, pupil, common_source)
    accumulate_pyramid_asterism_intensity!(execution_style(wfs.front_end.propagation.intensity), wfs, pupil, ast)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    pyramid_signal!(wfs, pupil, intensity, ast)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "PyramidWFS")
    ensure_pyramid_calibration!(wfs, pupil, common_source, det)
    accumulate_pyramid_asterism_intensity!(execution_style(wfs.front_end.propagation.intensity), wfs, pupil, ast)
    intensity = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, common_source; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    normalization_scale = wfs_detector_incidence_scale(det, common_source,
        eltype(frame))
    pyramid_signal!(wfs, pupil, frame, ast, normalization_scale)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end
