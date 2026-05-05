function measure!(mode::Geometric, wfs::PyramidWFS, tel::Telescope)
    geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask)
    gain = inv(1 + wfs.params.modulation)
    @. wfs.state.slopes = gain * wfs.state.slopes * wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.pupil_samples
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive PyramidWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::PyramidWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::SpectralSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::PyramidWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    ensure_pyramid_calibration!(wfs, tel, src)
    pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    pyramid_signal!(wfs, tel, intensity, src)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    ensure_pyramid_calibration!(wfs, tel, src)
    accumulate_pyramid_spectral_intensity!(execution_style(wfs.state.intensity), wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    pyramid_signal!(wfs, tel, intensity, spectral_reference_source(src))
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, tel, src)
    pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    pyramid_signal!(wfs, tel, frame, src)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::SpectralSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_pyramid_calibration!(wfs, tel, src)
    accumulate_pyramid_spectral_intensity!(execution_style(wfs.state.intensity), wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    pyramid_signal!(wfs, tel, frame, spectral_reference_source(src))
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    ensure_pyramid_calibration!(wfs, tel, ast.sources[1])
    accumulate_pyramid_asterism_intensity!(execution_style(wfs.state.intensity), wfs, tel, ast)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    pyramid_signal!(wfs, tel, intensity)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    ensure_pyramid_calibration!(wfs, tel, ast.sources[1])
    accumulate_pyramid_asterism_intensity!(execution_style(wfs.state.intensity), wfs, tel, ast)
    intensity = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    pyramid_signal!(wfs, tel, frame)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end
