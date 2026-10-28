@inline function _pyramid_modulated_pupil_field!(out, amplitude, opd,
    modulation::PreparedFocalPlaneModulation, point, amplitude_scale, opd_to_cycles)
    weight = modulation.amplitude_weights[point]
    @views @. out = amplitude_scale * weight * amplitude *
        modulation.phases[:, :, point] * cispi(opd_to_cycles * opd)
    return out
end

@inline function _pyramid_modulated_pupil_field!(out, amplitude, opd,
    modulation::Optics.PreparedSeparableFocalPlaneModulation, point,
    amplitude_scale, opd_to_cycles)
    weight = modulation.amplitude_weights[point]
    axis_1 = view(modulation.axis_1_factors, :, point)
    axis_2 = transpose(view(modulation.axis_2_factors, :, point))
    @. out = amplitude_scale * weight * amplitude * (axis_1 * axis_2) *
        cispi(opd_to_cycles * opd)
    return out
end

@inline function _pyramid_modulated_electric_field!(out, field,
    modulation::PreparedFocalPlaneModulation, point)
    weight = modulation.amplitude_weights[point]
    @views @. out = weight * field * modulation.phases[:, :, point]
    return out
end

@inline function _pyramid_modulated_electric_field!(out, field,
    modulation::Optics.PreparedSeparableFocalPlaneModulation, point)
    weight = modulation.amplitude_weights[point]
    axis_1 = view(modulation.axis_1_factors, :, point)
    axis_2 = transpose(view(modulation.axis_2_factors, :, point))
    @. out = weight * field * (axis_1 * axis_2)
    return out
end

function build_pyramid_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_pyramid_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_pyramid_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    # This is cold setup. Form the phase in at least Float64 before rounding
    # to runtime storage: Float32 phase error grows with the padded FFT side.
    W = promote_type(T, Float64)
    scale = -W(n + 1) / W(n)
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cispi(phase)
    end
    return phasor
end

function _build_pyramid_phasor!(::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    host = Matrix{Complex{T}}(undef, size(phasor))
    _build_pyramid_phasor!(ScalarCPUStyle(), host)
    copyto!(phasor, host)
    return phasor
end

function build_pyramid_mask!(wfs::PyramidWFS, pupil::PupilFunction)
    mask = pyramid_propagation_workspace(wfs).pyramid_mask
    _build_pyramid_mask!(execution_style(mask), mask, wfs, pupil)
    return mask
end

@inline function _build_pyramid_mask!(::ScalarCPUStyle, mask,
    wfs::PyramidWFS, pupil::PupilFunction)
    build_pyramid_mask_host!(mask, wfs, pupil)
    return mask
end

function _build_pyramid_mask!(::AcceleratorStyle, mask,
    wfs::PyramidWFS, pupil::PupilFunction)
    copyto!(mask, host_pyramid_mask(wfs, pupil))
    return mask
end

function host_pyramid_mask(wfs::PyramidWFS, pupil::PupilFunction)
    n = size(pyramid_propagation_workspace(wfs).pyramid_mask, 1)
    T = real(eltype(pyramid_propagation_workspace(wfs).pyramid_mask))
    host = Matrix{Complex{T}}(undef, n, n)
    build_pyramid_mask_host!(host, wfs, pupil)
    return host
end

function build_pyramid_mask_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, pupil::PupilFunction) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.front_end.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    rooftop_pixels = wfs.front_end.phase_mask.rooftop * wfs.front_end.phase_mask.diffraction_padding / sqrt(T(2))
    norma = T(_pupil_resolution(pupil)) / T(n_sub)
    lim = T(π)
    if wfs.front_end.phase_mask.psf_centering
        xvals = range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        xvals = range(-lim, lim; length=n + 1)[1:n]
    end
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    sx = wfs.front_end.phase_mask.pupil_shift_x_pixels
    sy = wfs.front_end.phase_mask.pupil_shift_y_pixels
    θ = wfs.front_end.phase_mask.rotation_rad
    cθ = cos(θ)
    sθ = sin(θ)
    @inbounds for i in 1:n, j in 1:n
        x_ = xvals[i]
        y_ = xvals[j]
        x = x_ * cθ - y_ * sθ
        y = y_ * cθ + x_ * sθ
        p1 = x * r + x_ * sx[1] + y * r - y_ * sy[1] + rooftop_pixels
        p2 = -x * r + x_ * sx[2] + y * r - y_ * sy[2]
        p3 = -x * r + x_ * sx[3] - y * r - y_ * sy[3] + rooftop_pixels
        p4 = x * r + x_ * sx[4] - y * r - y_ * sy[4]
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function accumulate_pyramid_focal_intensity!(out::AbstractMatrix,
    front_end::PyramidOpticalFrontEnd)
    propagation = pyramid_propagation_workspace(front_end)
    if front_end.phase_mask.psf_centering
        @. propagation.focal_field = propagation.focal_field * propagation.phasor
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        @. propagation.focal_field = propagation.focal_field * propagation.pyramid_mask
        copyto!(propagation.pupil_field, propagation.focal_field)
        execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    else
        execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        fftshift2d!(propagation.pupil_field, propagation.focal_field)
        @. propagation.pupil_field = propagation.pupil_field * propagation.pyramid_mask
        execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    end
    @. propagation.temp = abs2(propagation.pupil_field)
    out .+= propagation.temp
    return out
end

@inline accumulate_pyramid_focal_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS) =
    accumulate_pyramid_focal_intensity!(out, wfs.front_end)

@inline function _pyramid_modulation_batch_weights(batch,
    front_end::PyramidOpticalFrontEnd, modulation)
    return modulation === front_end.modulation ? batch.operating_weights : nothing
end

@inline function pyramid_bfft_intensity_scale(
    ::Type{T}, pad::Integer,
) where {T<:AbstractFloat}
    inverse_scale = inv(T(pad) * T(pad))
    return inverse_scale * inverse_scale
end

@inline _finish_pyramid_transpose!(::typeof(launch_kernel_async!), style) = nothing
@inline _finish_pyramid_transpose!(::typeof(launch_kernel!), style) =
    synchronize_backend!(style)

@inline function _pyramid_transpose!(style::AcceleratorStyle, destination, source,
    launch::L) where {L}
    pyramid_modulation_transpose_kernel!(style.backend, (32, 8, 1))(
        destination, source; ndrange=(cld(size(source, 1), 32) * 32,
            cld(size(source, 2), 32) * 8, size(source, 3)))
    _finish_pyramid_transpose!(launch, style)
    return destination
end

@inline function _pyramid_fft_mask_bfft!(::NoPyramidModulationFFTWorkspace,
    batch, mask, style, launch::L, execute_fft!::F) where {L,F}
    execute_fft!(batch.field_stack, batch.fft_plan)
    launch(style, pyramid_modulation_batch_mask_kernel!, batch.field_stack, mask,
        size(mask, 1), batch.batch_size; ndrange=size(mask))
    execute_fft!(batch.field_stack, batch.bfft_plan)
    return nothing
end

@inline function _pyramid_fft_mask_bfft!(workspace::PyramidTransposeFFTWorkspace,
    batch, mask, style, launch::L, execute_fft!::F) where {L,F}
    stack = batch.field_stack
    transposed = workspace.field_stack
    execute_fft!(stack, batch.fft_plan)
    _pyramid_transpose!(style, transposed, stack, launch)
    execute_fft!(transposed, workspace.fft_plan)
    launch(style, pyramid_modulation_transposed_mask_kernel!, transposed, mask,
        size(mask, 1), batch.batch_size; ndrange=size(mask))
    execute_fft!(transposed, workspace.bfft_plan)
    _pyramid_transpose!(style, stack, transposed, launch)
    execute_fft!(stack, batch.bfft_plan)
    return nothing
end

@inline function _propagate_pyramid_modulation_batch!(out,
    front_end::PyramidOpticalFrontEnd,
    batch::PyramidModulationBatchWorkspace, style::AcceleratorStyle)
    stack = batch.field_stack
    pad = size(out, 1)
    intensity_scale = pyramid_bfft_intensity_scale(eltype(out), pad)
    _pyramid_fft_mask_bfft!(batch.fft_workspace, batch,
        pyramid_propagation_workspace(front_end).pyramid_mask,
        style, launch_kernel!, execute_fft_plan!)
    launch_kernel!(style, pyramid_modulation_batch_intensity_kernel!,
        out, stack, intensity_scale, pad, batch.batch_size;
        ndrange=size(out))
    return nothing
end

@inline function _pyramid_pupil_modulation_batch!(
    ::NoPyramidModulationBatchWorkspace, out,
    front_end::PyramidOpticalFrontEnd, amplitude, opd, modulation,
    amplitude_scale, opd_to_cycles, offset::Int, resolution::Int)
    return false
end

function _pyramid_pupil_modulation_batch!(
    batch::PyramidModulationBatchWorkspace, out,
    front_end::PyramidOpticalFrontEnd, amplitude, opd, modulation,
    amplitude_scale, opd_to_cycles, offset::Int, resolution::Int)
    weights = _pyramid_modulation_batch_weights(batch, front_end, modulation)
    isnothing(weights) && return false
    style = execution_style(out)
    pad = size(out, 1)
    point_count = modulation_point_count(modulation)
    launch_kernel!(style, pyramid_pupil_phase_kernel!, batch.pupil_phase,
        opd, opd_to_cycles, resolution; ndrange=size(batch.pupil_phase))
    phases = Optics._modulation_phase_storage(modulation)
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel!(style, pyramid_pupil_modulation_batch_kernel!,
            batch.field_stack, amplitude, batch.pupil_phase, phases, weights,
            pyramid_propagation_workspace(front_end).phasor,
            amplitude_scale, offset, resolution, first_point,
            pad, batch.batch_size; ndrange=size(batch.field_stack))
        _propagate_pyramid_modulation_batch!(out, front_end, batch, style)
    end
    return true
end

@inline function _pyramid_electric_field_modulation_batch!(
    ::NoPyramidModulationBatchWorkspace, out,
    front_end::PyramidOpticalFrontEnd, field, modulation,
    offset::Int, resolution::Int)
    return false
end

function _pyramid_electric_field_modulation_batch!(
    batch::PyramidModulationBatchWorkspace, out,
    front_end::PyramidOpticalFrontEnd, field, modulation,
    offset::Int, resolution::Int)
    weights = _pyramid_modulation_batch_weights(batch, front_end, modulation)
    isnothing(weights) && return false
    style = execution_style(out)
    pad = size(out, 1)
    point_count = modulation_point_count(modulation)
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel!(style,
            pyramid_electric_field_modulation_batch_kernel!,
            batch.field_stack, field, Optics._modulation_phase_storage(modulation), weights,
            pyramid_propagation_workspace(front_end).phasor,
            offset, resolution, first_point, pad, batch.batch_size;
            ndrange=size(batch.field_stack))
        _propagate_pyramid_modulation_batch!(out, front_end, batch, style)
    end
    return true
end

@inline function _enqueue_pyramid_rate_batch!(
    style::AcceleratorStyle,
    out,
    front_end::PyramidOpticalFrontEnd,
    input::PupilFunction,
    source::Source,
    batch::PyramidModulationBatchWorkspace,
)
    propagation = pyramid_propagation_workspace(front_end)
    resolution = input.metadata.dimensions[1]
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    amplitude_scale = sqrt(T(photon_irradiance(source)) *
        T(input.metadata.sampling[1] * input.metadata.sampling[2]))
    opd_to_cycles = T(2) / T(wavelength(source))
    point_count = modulation_point_count(front_end.modulation)
    intensity_scale = pyramid_bfft_intensity_scale(T, pad)

    fill!(out, zero(T))
    launch_kernel_async!(style, pyramid_pupil_phase_kernel!, batch.pupil_phase,
        input.opd, opd_to_cycles, resolution; ndrange=size(batch.pupil_phase))
    phases = Optics._modulation_phase_storage(front_end.modulation)
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel_async!(style, pyramid_pupil_modulation_batch_kernel!,
            batch.field_stack, input.amplitude, batch.pupil_phase,
            phases, batch.operating_weights,
            propagation.phasor, amplitude_scale, offset,
            resolution, first_point, pad, batch.batch_size;
            ndrange=size(batch.field_stack))
        _pyramid_fft_mask_bfft!(batch.fft_workspace, batch,
            propagation.pyramid_mask, style, launch_kernel_async!, enqueue_fft_plan!)
        launch_kernel_async!(style,
            pyramid_modulation_batch_intensity_kernel!, out,
            batch.field_stack, intensity_scale, pad, batch.batch_size;
            ndrange=size(out))
    end
    return out
end

@inline function _enqueue_pyramid_rate_batch!(
    style::AcceleratorStyle,
    out,
    front_end::PyramidOpticalFrontEnd,
    input::ElectricField,
    ::Nothing,
    batch::PyramidModulationBatchWorkspace,
)
    propagation = pyramid_propagation_workspace(front_end)
    resolution = input.metadata.dimensions[1]
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    point_count = modulation_point_count(front_end.modulation)
    intensity_scale = pyramid_bfft_intensity_scale(T, pad)
    phases = Optics._modulation_phase_storage(front_end.modulation)

    fill!(out, zero(T))
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel_async!(style, pyramid_electric_field_modulation_batch_kernel!,
            batch.field_stack, input.values, phases, batch.operating_weights,
            propagation.phasor, offset, resolution, first_point, pad,
            batch.batch_size; ndrange=size(batch.field_stack))
        _pyramid_fft_mask_bfft!(batch.fft_workspace, batch,
            propagation.pyramid_mask, style, launch_kernel_async!, enqueue_fft_plan!)
        launch_kernel_async!(style,
            pyramid_modulation_batch_intensity_kernel!, out,
            batch.field_stack, intensity_scale, pad, batch.batch_size;
            ndrange=size(out))
    end
    return out
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, pupil, src,
        pyramid_operating_modulation(wfs))
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource,
    modulation::Union{PreparedFocalPlaneModulation,
        Optics.PreparedSeparableFocalPlaneModulation}) where {T<:AbstractFloat}
    return pyramid_intensity_core!(execution_style(out), out, wfs, pupil, src,
        modulation)
end

function pyramid_intensity_core!(::ScalarCPUStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource, modulation::PreparedFocalPlaneModulation) where {
    T<:AbstractFloat,
}
    prepare_pyramid_sampling!(wfs, pupil)
    propagation = pyramid_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2))
    amplitude = pupil.amplitude

    fill!(out, zero(T))
    _pyramid_pupil_modulation_batch!(
        propagation.modulation_batch,
        out,
        wfs.front_end,
        amplitude,
        pupil.opd,
        modulation,
        amp_scale,
        opd_to_cycles,
        ox,
        n,
    ) && return out
    fill!(propagation.field, zero(eltype(propagation.field)))
    for p in 1:modulation_point_count(modulation)
        _pyramid_modulated_pupil_field!(
            view(propagation.field, ox+1:ox+n, oy+1:oy+n), amplitude,
            pupil.opd, modulation, p, amp_scale, opd_to_cycles)
        copyto!(propagation.focal_field, propagation.field)
        accumulate_pyramid_focal_intensity!(out, wfs)
    end
    return out
end

function pyramid_intensity_core!(::AcceleratorStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource, modulation::Union{PreparedFocalPlaneModulation,
        Optics.PreparedSeparableFocalPlaneModulation}) where {
    T<:AbstractFloat,
}
    prepare_pyramid_sampling!(wfs, pupil)
    propagation = pyramid_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2))
    amplitude = pupil.amplitude

    fill!(out, zero(T))
    _pyramid_pupil_modulation_batch!(
        propagation.modulation_batch,
        out,
        wfs.front_end,
        amplitude,
        pupil.opd,
        modulation,
        amp_scale,
        opd_to_cycles,
        ox,
        n,
    ) && return out
    for p in 1:modulation_point_count(modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        _pyramid_modulated_pupil_field!(
            view(propagation.field, ox+1:ox+n, oy+1:oy+n), amplitude,
            pupil.opd, modulation, p, amp_scale, opd_to_cycles)
        copyto!(propagation.focal_field, propagation.field)
        accumulate_pyramid_focal_intensity!(out, wfs)
    end
    return out
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, pupil, src)
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction, src::ExtendedSource) where {T<:AbstractFloat}
    return accumulate_pyramid_extended_intensity!(execution_style(out), out, wfs, pupil, src)
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction, src::LGSSource) where {T<:AbstractFloat}
    pyramid_intensity_core!(out, wfs, pupil, src)
    apply_lgs_elongation!(sodium_layer_profile_style(src), out, wfs,
        pupil, src)
    return out
end

function pyramid_modulation_frame!(out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource) where {T<:AbstractFloat}
    prepare_pyramid_sampling!(wfs, pupil)
    propagation = pyramid_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    fft_scale2 = inv(T(pad) * T(pad))
    if size(out) != (pad, pad)
        throw(DimensionMismatchError("modulation frame size must match pyramid sampling"))
    end
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    modulation = pyramid_operating_modulation(wfs)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2))
    amplitude = pupil.amplitude

    fill!(out, zero(T))
    for p in 1:modulation_point_count(modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        _pyramid_modulated_pupil_field!(
            view(propagation.field, ox+1:ox+n, oy+1:oy+n), amplitude,
            pupil.opd, modulation, p, amp_scale, opd_to_cycles)
        copyto!(propagation.focal_field, propagation.field)
        if wfs.front_end.phase_mask.psf_centering
            @. propagation.focal_field = propagation.focal_field * propagation.phasor
            execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
        else
            execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
            @. propagation.intensity = abs2(propagation.focal_field) * fft_scale2
            fftshift2d!(propagation.temp, propagation.intensity)
            out .+= propagation.temp
            continue
        end
        @. propagation.temp = abs2(propagation.focal_field) * fft_scale2
        out .+= propagation.temp
    end
    return out
end

"""
    pyramid_modulation_frame(wfs, pupil, source)

Allocate and form the focal-plane modulation-cycle intensity used by a
gain-sensing camera. The returned array has the prepared Pyramid front end's
native focal-plane sampling.
"""
function pyramid_modulation_frame(wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource)
    out = similar(pyramid_propagation_workspace(wfs).intensity)
    return pyramid_modulation_frame!(out, wfs, pupil, src)
end

function apply_lgs_elongation!(::NoSodiumLayerProfileStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, ::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    propagation = pyramid_propagation_workspace(wfs)
    propagation.elongation_kernel = apply_elongation!(
        out,
        lgs_elongation_factor(src),
        propagation.scratch,
        propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::SampledSodiumLayerProfileStyle,
    out::AbstractMatrix{T}, wfs::PyramidWFS, pupil::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, pupil, src)
    propagation = pyramid_propagation_workspace(wfs)
    apply_lgs_convolution!(
        out,
        propagation.lgs_kernel_fft,
        propagation.focal_field,
        propagation.fft_plan,
        propagation.pupil_field,
        propagation.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::PyramidWFS, pupil::PupilFunction, src::LGSSource)
    profile = src.params.sodium_layer_profile
    if profile === nothing
        return wfs
    end
    propagation = pyramid_propagation_workspace(wfs)
    pad = size(propagation.intensity, 1)
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
        propagation.focal_field,
        propagation.fft_plan,
    )
    propagation.lgs_kernel_tag = tag
    return wfs
end

"""
    Optics.update_cycle_averaged_circular_modulation!(wfs::PyramidWFS, radius; enabled=true)

Update prepared circular pupil-tilt phases in place between frames, preserving
captured buffer bindings.
The sample count, weights, precision and phase origin remain fixed. This is a
synchronous control operation, not a captured waveform node.
"""
function Optics.update_cycle_averaged_circular_modulation!(
    wfs::PyramidWFS, radius::T; enabled::Bool=true) where {T<:AbstractFloat}
    modulation = pyramid_operating_modulation(wfs)
    Optics.update_cycle_averaged_circular_modulation!(modulation, radius; enabled)
    return wfs
end
