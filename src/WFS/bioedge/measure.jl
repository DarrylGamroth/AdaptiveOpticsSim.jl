function bioedge_intensity_core!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource; apply_lgs::Bool=false) where {T<:AbstractFloat}
    prepare_bioedge_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2 /
        wfs.params.modulation_points))

    fill!(out, zero(T))
    profile = apply_lgs ? lgs_profile(src) : LGSProfileNone()
    apply_lgs && ensure_bioedge_lgs_kernel!(profile, wfs, tel, src)

    @inbounds for p in 1:wfs.params.modulation_points
        fill!(wfs.state.field, zero(eltype(wfs.state.field)))
        @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
            wfs.state.modulation_phases[:, :, p] * cispi(opd_to_cycles * tel.state.opd)
        copyto!(wfs.state.focal_field, wfs.state.field)
        if wfs.params.psf_centering
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
        else
            execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
            fftshift2d!(wfs.state.fft_buffer, wfs.state.focal_field)
        end

        focal_source = wfs.params.psf_centering ? wfs.state.focal_field : wfs.state.fft_buffer
        lgs_fft_buffer = wfs.params.psf_centering ? wfs.state.fft_buffer : wfs.state.focal_field
        lgs_ifft_buffer = wfs.state.pupil_field
        @inbounds for k in 1:4
            @views @. wfs.state.pupil_field = focal_source * wfs.state.bioedge_masks[:, :, k]
            execute_fft_plan!(wfs.state.pupil_field, wfs.state.ifft_plan)
            @. wfs.state.temp = abs2(wfs.state.pupil_field)
            if apply_lgs
                apply_bioedge_lgs_profile!(profile, wfs, src, lgs_fft_buffer, lgs_ifft_buffer)
            end
            oxq = k in (3, 4) ? pad : 0
            oyq = k in (2, 4) ? pad : 0
            @views out[oxq+1:oxq+pad, oyq+1:oyq+pad] .+= wfs.state.temp
        end
    end
    return out
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, tel, src; apply_lgs=false)
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::LGSSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, tel, src; apply_lgs=true)
end

function measure!(mode::Geometric, wfs::BioEdgeWFS, tel::Telescope)
    edge_geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask, wfs.state.edge_mask)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.pupil_samples
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive BioEdgeWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::BioEdgeWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, ast::Asterism)
    return measure!(sensing_mode(wfs), wfs, tel, ast)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    return bioedge_signal!(wfs, tel, intensity, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    return bioedge_signal!(wfs, tel, intensity, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1))
    return bioedge_signal!(wfs, tel, frame, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1))
    return bioedge_signal!(wfs, tel, frame, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    ensure_bioedge_calibration!(wfs, tel, ast.sources[1])
    accumulate_bioedge_asterism_intensity!(execution_style(wfs.state.intensity), wfs, tel, ast)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    return bioedge_signal!(wfs, tel, intensity)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    ensure_bioedge_calibration!(wfs, tel, ast.sources[1])
    accumulate_bioedge_asterism_intensity!(execution_style(wfs.state.intensity), wfs, tel, ast)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1))
    return bioedge_signal!(wfs, tel, frame)
end

