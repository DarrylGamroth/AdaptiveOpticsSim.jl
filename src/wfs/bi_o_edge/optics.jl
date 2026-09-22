function accumulate_bi_o_edge_masked_pupils!(out::AbstractMatrix,
    front_end::BiOEdgeOpticalFrontEnd)
    return accumulate_bi_o_edge_masked_pupils!(out, front_end,
        NoPreparedFourPupilLGS())
end

function accumulate_bi_o_edge_masked_pupils!(out::AbstractMatrix,
    front_end::BiOEdgeOpticalFrontEnd,
    lgs_model::AbstractPreparedFourPupilLGS)
    propagation = bi_o_edge_propagation_workspace(front_end)
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
            propagation.bi_o_edge_masks[:, :, branch]
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

function bi_o_edge_intensity_core!(out::AbstractMatrix{T}, wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::AbstractSource; apply_lgs::Bool=false) where {T<:AbstractFloat}
    return bi_o_edge_intensity_core!(out, wfs, pupil, src,
        bi_o_edge_operating_modulation(wfs); apply_lgs)
end

function bi_o_edge_intensity_core!(out::AbstractMatrix{T}, wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource,
    modulation::PreparedFocalPlaneModulation;
    apply_lgs::Bool=false) where {T<:AbstractFloat}
    require_leaf_source(src, "Bi-O-edge WFS")
    prepare_bi_o_edge_sampling!(wfs, pupil)
    propagation = bi_o_edge_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2))
    amplitude = pupil.amplitude

    fill!(out, zero(T))
    profile_style = apply_lgs ? sodium_layer_profile_style(src) :
        NoSodiumLayerProfileStyle()
    apply_lgs && ensure_bi_o_edge_lgs_kernel!(profile_style, wfs, pupil, src)

    @inbounds for p in 1:modulation_point_count(modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * amplitude *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * pupil.opd)
        copyto!(propagation.focal_field, propagation.field)
        if !apply_lgs
            accumulate_bi_o_edge_masked_pupils!(out, wfs.front_end)
            continue
        end
        if wfs.front_end.amplitude_mask.psf_centering
            @. propagation.focal_field = propagation.focal_field * propagation.phasor
            execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        else
            execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
            fftshift2d!(propagation.fft_buffer, propagation.focal_field)
        end

        focal_source = wfs.front_end.amplitude_mask.psf_centering ?
            propagation.focal_field : propagation.fft_buffer
        lgs_fft_buffer = wfs.front_end.amplitude_mask.psf_centering ?
            propagation.fft_buffer : propagation.focal_field
        lgs_ifft_buffer = propagation.pupil_field
        @inbounds for k in 1:4
            @views @. propagation.pupil_field = focal_source *
                propagation.bi_o_edge_masks[:, :, k]
            execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
            @. propagation.temp = abs2(propagation.pupil_field)
            if apply_lgs
                apply_bi_o_edge_sodium_layer_profile!(profile_style, wfs, src,
                    lgs_fft_buffer, lgs_ifft_buffer)
            end
            oxq = k in (3, 4) ? pad : 0
            oyq = k in (2, 4) ? pad : 0
            @views out[oxq+1:oxq+pad, oyq+1:oyq+pad] .+= propagation.temp
        end
    end
    return out
end

function bi_o_edge_intensity!(out::AbstractMatrix{T}, wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::AbstractSource) where {T<:AbstractFloat}
    return bi_o_edge_intensity_core!(out, wfs, pupil, src; apply_lgs=false)
end

function bi_o_edge_intensity!(out::AbstractMatrix{T}, wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    return bi_o_edge_intensity_core!(out, wfs, pupil, src; apply_lgs=true)
end

function apply_lgs_elongation!(::NoSodiumLayerProfileStyle,
    intensity::AbstractMatrix{T}, wfs::BiOEdgeWFS, ::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    propagation = bi_o_edge_propagation_workspace(wfs)
    propagation.elongation_kernel = apply_elongation!(
        intensity,
        lgs_elongation_factor(src),
        propagation.scratch,
        propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::SampledSodiumLayerProfileStyle,
    intensity::AbstractMatrix{T}, wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, pupil, src)
    propagation = bi_o_edge_propagation_workspace(wfs)
    apply_lgs_convolution!(
        intensity,
        propagation.lgs_kernel_fft,
        propagation.fft_buffer,
        propagation.fft_plan,
        propagation.pupil_field,
        propagation.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::BiOEdgeWFS, pupil::PupilFunction, src::LGSSource)
    profile = src.params.sodium_layer_profile
    if profile === nothing
        return wfs
    end
    propagation = bi_o_edge_propagation_workspace(wfs)
    pad = size(propagation.fft_buffer, 1)
    padding = propagation.effective_resolution / _pupil_resolution(pupil)
    pixel_scale = lgs_pixel_scale(_pupil_diameter_m(pupil), padding,
        wavelength(src))
    tag = lgs_kernel_signature(
        pupil,
        src,
        pad,
        wfs.front_end.pupil_samples,
        pixel_scale,
        eltype(propagation.intensity);
        model=:subaperture_average,
    )
    if size(propagation.lgs_kernel_fft, 1) == pad &&
        propagation.lgs_kernel_tag == tag
        return wfs
    end
    propagation.lgs_kernel_fft = lgs_average_kernel_fft(
        pupil,
        src,
        pad,
        wfs.front_end.pupil_samples,
        pixel_scale,
        propagation.fft_buffer,
        propagation.fft_plan,
    )
    propagation.lgs_kernel_tag = tag
    return wfs
end
