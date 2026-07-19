function accumulate_bioedge_masked_pupils!(out::AbstractMatrix,
    front_end::BioEdgeOpticalFrontEnd)
    return accumulate_bioedge_masked_pupils!(out, front_end,
        NoPreparedFourPupilLGS())
end

function accumulate_bioedge_masked_pupils!(out::AbstractMatrix,
    front_end::BioEdgeOpticalFrontEnd,
    lgs_model::AbstractPreparedFourPupilLGS)
    propagation = front_end.propagation
    pad = size(propagation.field, 1)
    if front_end.amplitude_mask.psf_centering
        @. propagation.focal_field = propagation.focal_field *
            propagation.phasor
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        focal_source = propagation.focal_field
    else
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        fftshift2d!(propagation.fft_buffer, propagation.focal_field)
        focal_source = propagation.fft_buffer
    end
    @inbounds for branch in 1:4
        @views @. propagation.pupil_field = focal_source *
            propagation.bioedge_masks[:, :, branch]
        execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
        @. propagation.temp = abs2(propagation.pupil_field)
        lgs_fft_buffer = front_end.amplitude_mask.psf_centering ?
            propagation.fft_buffer : propagation.focal_field
        apply_prepared_four_pupil_lgs!(lgs_model, propagation.temp,
            propagation.scratch, lgs_fft_buffer, propagation.fft_plan,
            propagation.pupil_field, propagation.ifft_plan)
        axis_1_offset = branch in (3, 4) ? pad : 0
        axis_2_offset = branch in (2, 4) ? pad : 0
        @views out[axis_1_offset+1:axis_1_offset+pad,
            axis_2_offset+1:axis_2_offset+pad] .+= propagation.temp
    end
    return out
end

function bioedge_intensity_core!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, pupil::PupilFunction,
    src::AbstractSource; apply_lgs::Bool=false) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, pupil, src,
        bioedge_operating_modulation(wfs); apply_lgs)
end

function bioedge_intensity_core!(out::AbstractMatrix{T}, wfs::BioEdgeWFS,
    pupil::PupilFunction, src::AbstractSource,
    modulation::PreparedFocalPlaneModulation;
    apply_lgs::Bool=false) where {T<:AbstractFloat}
    require_leaf_source(src, "BioEdgeWFS")
    prepare_bioedge_sampling!(wfs, pupil)
    n = _pupil_resolution(pupil)
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2))
    amplitude = pupil.amplitude

    fill!(out, zero(T))
    profile = apply_lgs ? lgs_profile(src) : LGSProfileNone()
    apply_lgs && ensure_bioedge_lgs_kernel!(profile, wfs, pupil, src)

    @inbounds for p in 1:modulation_point_count(modulation)
        fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. wfs.front_end.propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * amplitude *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * pupil.opd)
        copyto!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.field)
        if !apply_lgs
            accumulate_bioedge_masked_pupils!(out, wfs.front_end)
            continue
        end
        if wfs.front_end.amplitude_mask.psf_centering
            @. wfs.front_end.propagation.focal_field = wfs.front_end.propagation.focal_field * wfs.front_end.propagation.phasor
            execute_fft_plan!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.fft_plan)
        else
            execute_fft_plan!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.fft_plan)
            fftshift2d!(wfs.front_end.propagation.fft_buffer, wfs.front_end.propagation.focal_field)
        end

        focal_source = wfs.front_end.amplitude_mask.psf_centering ? wfs.front_end.propagation.focal_field : wfs.front_end.propagation.fft_buffer
        lgs_fft_buffer = wfs.front_end.amplitude_mask.psf_centering ? wfs.front_end.propagation.fft_buffer : wfs.front_end.propagation.focal_field
        lgs_ifft_buffer = wfs.front_end.propagation.pupil_field
        @inbounds for k in 1:4
            @views @. wfs.front_end.propagation.pupil_field = focal_source * wfs.front_end.propagation.bioedge_masks[:, :, k]
            execute_fft_plan!(wfs.front_end.propagation.pupil_field, wfs.front_end.propagation.ifft_plan)
            @. wfs.front_end.propagation.temp = abs2(wfs.front_end.propagation.pupil_field)
            if apply_lgs
                apply_bioedge_lgs_profile!(profile, wfs, src, lgs_fft_buffer, lgs_ifft_buffer)
            end
            oxq = k in (3, 4) ? pad : 0
            oyq = k in (2, 4) ? pad : 0
            @views out[oxq+1:oxq+pad, oyq+1:oyq+pad] .+= wfs.front_end.propagation.temp
        end
    end
    return out
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, pupil::PupilFunction,
    src::AbstractSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, pupil, src; apply_lgs=false)
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, pupil::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, pupil, src; apply_lgs=true)
end

function measure!(mode::Geometric, wfs::BioEdgeWFS, pupil::PupilFunction)
    edge_geometric_slopes!(wfs.estimator.state.slopes, pupil.opd, wfs.estimator.state.valid_mask, wfs.estimator.state.edge_mask)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function measure!(::Geometric, wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "geometric BioEdgeWFS")
    return measure!(Geometric(), wfs, pupil)
end

function measure!(::Geometric, wfs::BioEdgeWFS, pupil::PupilFunction, src::LGSSource)
    slopes = measure!(Geometric(), wfs, pupil)
    n_sub = wfs.estimator.params.pupil_samples
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction)
    throw(InvalidConfiguration("Diffractive BioEdgeWFS requires a source; call measure!(wfs, pupil, src)."))
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction)
    return measure!(sensing_mode(wfs), wfs, pupil)
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction, ast::Asterism)
    return measure!(sensing_mode(wfs), wfs, pupil, ast)
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::BioEdgeWFS, pupil::PupilFunction, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    ensure_bioedge_calibration!(wfs, pupil, src)
    bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    return bioedge_signal!(wfs, pupil, intensity, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, src::LGSSource)
    ensure_bioedge_calibration!(wfs, pupil, src)
    bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    return bioedge_signal!(wfs, pupil, intensity, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_bioedge_calibration!(wfs, pupil, src, det)
    bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1), det)
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    return bioedge_signal!(wfs, pupil, frame, src, normalization_scale)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, src::LGSSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_bioedge_calibration!(wfs, pupil, src, det)
    bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, src; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1), det)
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    return bioedge_signal!(wfs, pupil, frame, src, normalization_scale)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, ast::Asterism)
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "BioEdgeWFS")
    ensure_bioedge_calibration!(wfs, pupil, common_source)
    accumulate_bioedge_asterism_intensity!(execution_style(wfs.front_end.propagation.intensity), wfs, pupil, ast)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    return bioedge_signal!(wfs, pupil, intensity, ast)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, pupil::PupilFunction, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(pupil.opd)
    common_source = common_wfs_calibration_source(ast, "BioEdgeWFS")
    ensure_bioedge_calibration!(wfs, pupil, common_source, det)
    accumulate_bioedge_asterism_intensity!(execution_style(wfs.front_end.propagation.intensity), wfs, pupil, ast)
    intensity = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = capture!(det, intensity, common_source; rng=rng)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1), det)
    normalization_scale = wfs_detector_incidence_scale(det, common_source,
        eltype(frame))
    return bioedge_signal!(wfs, pupil, frame, ast, normalization_scale)
end
