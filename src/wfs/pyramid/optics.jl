function build_pyramid_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_pyramid_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_pyramid_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_pyramid_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, pyramid_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_pyramid_mask!(wfs::PyramidWFS, pupil::PupilFunction)
    mask = pyramid_propagation_workspace(wfs).pyramid_mask
    copyto!(mask, host_pyramid_mask(wfs, pupil))
    return mask
end

@inline build_pyramid_shifted_masks!(::PyramidWFS{<:Geometric},
    ::PupilFunction) = nothing

function build_pyramid_shifted_masks!(wfs::PyramidWFS{<:Diffractive},
    pupil::PupilFunction)
    batch = pyramid_propagation_workspace(wfs).modulation_batch
    return _build_pyramid_shifted_masks!(
        execution_style(pyramid_propagation_workspace(wfs).pyramid_mask),
        batch,
        wfs,
        pupil,
    )
end

@inline _build_pyramid_shifted_masks!(::ExecutionStyle,
    ::NoPyramidModulationBatchWorkspace, ::PyramidWFS,
    ::PupilFunction) = nothing

@inline _build_pyramid_shifted_masks!(::ExecutionStyle,
    ::PyramidModulationBatchWorkspace, ::PyramidWFS,
    ::PupilFunction) = nothing

@inline function _pyramid_shifted_mask_parameters(
    batch::PyramidShiftedMaskModulationWorkspace,
    wfs::PyramidWFS,
    pupil::PupilFunction,
)
    masks = batch.shifted_masks
    T = real(eltype(masks))
    pad = size(masks, 1)
    n_sub = wfs.estimator.params.pupil_samples
    phase_mask = wfs.front_end.phase_mask
    separation = something(phase_mask.n_pix_separation, 0)
    r = (T(n_sub) + T(separation)) * phase_mask.mask_scale / T(2)
    norma = T(_pupil_resolution(pupil)) / T(n_sub)
    rooftop_pixels = phase_mask.rooftop * phase_mask.diffraction_padding /
        sqrt(T(2))
    coordinate_start = -T(pi) * (one(T) - inv(T(pad)))
    coordinate_step = T(2pi) / T(pad)
    rotation_sin, rotation_cos = sincos(phase_mask.rotation_rad)
    shift_x = wfs.estimator.state.shift_x
    shift_y = wfs.estimator.state.shift_y
    return (
        masks,
        batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad,
        r,
        norma,
        rooftop_pixels,
        coordinate_start,
        coordinate_step,
        rotation_cos,
        rotation_sin,
        shift_x,
        shift_y,
        pad,
        size(masks, 3),
    )
end

function _build_pyramid_shifted_masks!(::ScalarCPUStyle,
    batch::PyramidShiftedMaskModulationWorkspace,
    wfs::PyramidWFS,
    pupil::PupilFunction)
    masks, axis_1_shifts_rad, axis_2_shifts_rad, r, norma,
        rooftop_pixels, coordinate_start, coordinate_step, rotation_cos,
        rotation_sin, shift_x, shift_y, pad, point_count =
        _pyramid_shifted_mask_parameters(batch, wfs, pupil)
    @inbounds for point in 1:point_count, axis_2 in 1:pad, axis_1 in 1:pad
        unrotated_x = coordinate_start + (axis_1 - 1) * coordinate_step +
            axis_1_shifts_rad[point]
        unrotated_y = coordinate_start + (axis_2 - 1) * coordinate_step +
            axis_2_shifts_rad[point]
        x = unrotated_x * rotation_cos - unrotated_y * rotation_sin
        y = unrotated_y * rotation_cos + unrotated_x * rotation_sin
        phase_1 = x * r + unrotated_x * shift_x[1] + y * r -
            unrotated_y * shift_y[1] + rooftop_pixels
        phase_2 = -x * r + unrotated_x * shift_x[2] + y * r -
            unrotated_y * shift_y[2]
        phase_3 = -x * r + unrotated_x * shift_x[3] - y * r -
            unrotated_y * shift_y[3] + rooftop_pixels
        phase_4 = x * r + unrotated_x * shift_x[4] - y * r -
            unrotated_y * shift_y[4]
        phase = -max(max(phase_1, phase_2), max(phase_3, phase_4)) * norma
        masks[axis_1, axis_2, point] = cis(phase)
    end
    return masks
end

function _build_pyramid_shifted_masks!(style::AcceleratorStyle,
    batch::PyramidShiftedMaskModulationWorkspace,
    wfs::PyramidWFS,
    pupil::PupilFunction)
    masks, axis_1_shifts_rad, axis_2_shifts_rad, r, norma,
        rooftop_pixels, coordinate_start, coordinate_step, rotation_cos,
        rotation_sin, shift_x, shift_y, pad, point_count =
        _pyramid_shifted_mask_parameters(batch, wfs, pupil)
    launch_kernel!(style, pyramid_shifted_mask_stack_kernel!, masks,
        axis_1_shifts_rad, axis_2_shifts_rad, r, norma, rooftop_pixels,
        coordinate_start, coordinate_step, rotation_cos, rotation_sin,
        shift_x[1], shift_y[1], shift_x[2], shift_y[2],
        shift_x[3], shift_y[3], shift_x[4], shift_y[4],
        pad, point_count; ndrange=size(masks))
    return masks
end

function host_pyramid_mask(wfs::PyramidWFS, pupil::PupilFunction)
    n = size(pyramid_propagation_workspace(wfs).pyramid_mask, 1)
    T = eltype(pyramid_estimator_products(wfs).slopes)
    host = Matrix{Complex{T}}(undef, n, n)
    if wfs.front_end.phase_mask.old_mask
        build_pyramid_mask_old_host!(host, wfs, pupil)
    else
        build_pyramid_mask_new_host!(host, wfs, pupil)
    end
    return host
end

function build_pyramid_mask_new_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, pupil::PupilFunction) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
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
    sx = wfs.estimator.state.shift_x
    sy = wfs.estimator.state.shift_y
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

function build_pyramid_mask_old_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, pupil::PupilFunction) where {T<:AbstractFloat}
    n_tot = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    sx = wfs.estimator.state.shift_x
    sy = wfs.estimator.state.shift_y
    norma = (T(_pupil_resolution(pupil)) / T(n_sub)) / 4
    fill!(mask, complex(zero(T), zero(T)))
    if wfs.front_end.phase_mask.psf_centering
        tip = centered_grid(T, n_tot ÷ 2, true)
        tilt = centered_grid(T, n_tot ÷ 2, true)
        Tip = repeat(tip', n_tot ÷ 2, 1)
        Tilt = repeat(tilt, 1, n_tot ÷ 2)
        Tip .-= mean(Tip)
        Tilt .-= mean(Tilt)
        q = T(n_sub + sep) * wfs.front_end.phase_mask.mask_scale
        @views begin
            mask[1:n_tot÷2, 1:n_tot÷2] .= cis.((Tip .* (q + sx[1]) .+ Tilt .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[2]) .+ Tilt .* (q - sy[2])) .* norma)
            mask[n_tot÷2+1:end, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[3]) .- Tilt .* (q + sy[3])) .* norma)
            mask[n_tot÷2+1:end, 1:n_tot÷2] .= cis.((Tip .* (q + sx[4]) .- Tilt .* (q + sy[4])) .* norma)
        end
    else
        d_pix = T(π) / T(n_tot)
        lim_p = T(π)
        lim_m = T(π) - 2 * d_pix
        tip1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip2 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tip4 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tilt2 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt4 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        Tip1 = repeat(tip1', length(tilt1), 1); Tilt1 = repeat(tilt1, 1, length(tip1))
        Tip2 = repeat(tip2', length(tilt2), 1); Tilt2 = repeat(tilt2, 1, length(tip2))
        Tip3 = repeat(tip3', length(tilt3), 1); Tilt3 = repeat(tilt3, 1, length(tip3))
        Tip4 = repeat(tip4', length(tilt4), 1); Tilt4 = repeat(tilt4, 1, length(tip4))
        Tip1 .-= mean(Tip1); Tilt1 .-= mean(Tilt1)
        Tip2 .-= mean(Tip2); Tilt2 .-= mean(Tilt2)
        Tip3 .-= mean(Tip3); Tilt3 .-= mean(Tilt3)
        Tip4 .-= mean(Tip4); Tilt4 .-= mean(Tilt4)
        q = T(n_sub + sep) * wfs.front_end.phase_mask.mask_scale
        @views begin
            mask[1:n_tot÷2+1, 1:n_tot÷2+1] .= cis.((Tip1 .* (q + sx[1]) .+ Tilt1 .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2+1, n_tot÷2:end] .= cis.((-Tip4 .* (q - sx[2]) .+ Tilt4 .* (q - sy[2])) .* norma)
            mask[n_tot÷2:end, n_tot÷2:end] .= cis.((-Tip3 .* (q - sx[3]) .- Tilt3 .* (q + sy[3])) .* norma)
            mask[n_tot÷2:end, 1:n_tot÷2+1] .= cis.((Tip2 .* (q + sx[4]) .- Tilt2 .* (q + sy[4])) .* norma)
        end
    end
    return mask
end

centered_grid(::Type{T}, n::Int, endpoint::Bool) where {T<:AbstractFloat} =
    axis_values(T, n, -T(π), T(π); endpoint=endpoint)

function axis_values(::Type{T}, n::Int, lo::T, hi::T; endpoint::Bool) where {T<:AbstractFloat}
    if endpoint
        return reshape(range(lo, hi; length=n), n, 1)
    end
    return reshape(range(lo; step=(hi - lo) / n, length=n), n, 1)
end

function _build_pyramid_mask!(::ScalarCPUStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, pupil::PupilFunction) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    pix_per_subap = T(_pupil_resolution(pupil)) / T(n_sub)
    norma = pix_per_subap
    lim = T(pi)
    x_vals = if wfs.front_end.phase_mask.psf_centering
        range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        range(-lim, lim; length=n, endpoint=false)
    end
    @inbounds for i in 1:n, j in 1:n
        x = x_vals[i]
        y = x_vals[j]
        p1 = x * r + y * r
        p2 = -x * r + y * r
        p3 = -x * r - y * r
        p4 = x * r - y * r
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function _build_pyramid_mask!(style::AcceleratorStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, pupil::PupilFunction) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.estimator.params.pupil_samples
    sep = wfs.front_end.phase_mask.n_pix_separation === nothing ? 0 : wfs.front_end.phase_mask.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.front_end.phase_mask.mask_scale / 2
    norma = T(_pupil_resolution(pupil)) / T(n_sub)
    lim = T(pi)
    start = wfs.front_end.phase_mask.psf_centering ? -lim * (one(T) - one(T) / T(n)) : -lim
    step = T(2) * lim / T(n)
    launch_kernel!(style, pyramid_mask_kernel!, mask, r, norma, start, step, n; ndrange=size(mask))
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
    wfs::PyramidWFS{<:Diffractive}) =
    accumulate_pyramid_focal_intensity!(out, wfs.front_end)

@inline function _pyramid_modulation_batch_weights(batch,
    front_end::PyramidOpticalFrontEnd, modulation)
    modulation === front_end.modulation && return batch.operating_weights
    modulation === front_end.calibration_modulation &&
        return batch.calibration_weights
    return nothing
end

function _prepare_pyramid_shifted_focal_field!(::ScalarCPUStyle,
    front_end::PyramidOpticalFrontEnd, amplitude, opd, amplitude_scale,
    opd_to_cycles, offset::Int, resolution::Int)
    propagation = pyramid_propagation_workspace(front_end)
    focal_field = propagation.focal_field
    fill!(focal_field, zero(eltype(focal_field)))
    @views @. focal_field[offset+1:offset+resolution,
        offset+1:offset+resolution] = amplitude_scale * amplitude *
        cispi(opd_to_cycles * opd) *
        propagation.phasor[offset+1:offset+resolution,
            offset+1:offset+resolution]
    execute_fft_plan!(focal_field, propagation.fft_plan)
    return focal_field
end

function _prepare_pyramid_shifted_focal_field!(style::AcceleratorStyle,
    front_end::PyramidOpticalFrontEnd, amplitude, opd, amplitude_scale,
    opd_to_cycles, offset::Int, resolution::Int)
    propagation = pyramid_propagation_workspace(front_end)
    focal_field = propagation.focal_field
    pad = size(focal_field, 1)
    launch_kernel!(style, pyramid_unmodulated_pupil_field_kernel!,
        focal_field, amplitude, opd, propagation.phasor, amplitude_scale,
        opd_to_cycles, offset, resolution, pad; ndrange=size(focal_field))
    execute_fft_plan!(focal_field, propagation.fft_plan)
    return focal_field
end

function _prepare_pyramid_shifted_focal_field!(::ScalarCPUStyle,
    front_end::PyramidOpticalFrontEnd, field, offset::Int, resolution::Int)
    propagation = pyramid_propagation_workspace(front_end)
    focal_field = propagation.focal_field
    fill!(focal_field, zero(eltype(focal_field)))
    @views @. focal_field[offset+1:offset+resolution,
        offset+1:offset+resolution] = field *
        propagation.phasor[offset+1:offset+resolution,
            offset+1:offset+resolution]
    execute_fft_plan!(focal_field, propagation.fft_plan)
    return focal_field
end

function _prepare_pyramid_shifted_focal_field!(style::AcceleratorStyle,
    front_end::PyramidOpticalFrontEnd, field, offset::Int, resolution::Int)
    propagation = pyramid_propagation_workspace(front_end)
    focal_field = propagation.focal_field
    pad = size(focal_field, 1)
    launch_kernel!(style, pyramid_unmodulated_electric_field_kernel!,
        focal_field, field, propagation.phasor, offset, resolution, pad;
        ndrange=size(focal_field))
    execute_fft_plan!(focal_field, propagation.fft_plan)
    return focal_field
end

function _propagate_pyramid_shifted_masks!(::ScalarCPUStyle, out,
    front_end::PyramidOpticalFrontEnd,
    batch::PyramidShiftedMaskModulationWorkspace)
    stack = batch.field_stack
    focal_field = pyramid_propagation_workspace(front_end).focal_field
    pad = size(out, 1)
    point_count = size(batch.shifted_masks, 3)
    @inbounds for first_point in 1:batch.batch_size:point_count
        for batch_index in 1:batch.batch_size
            point = first_point + batch_index - 1
            @views @. stack[:, :, batch_index] =
                batch.operating_weights[point] * focal_field *
                batch.shifted_masks[:, :, point]
        end
        execute_fft_plan!(stack, batch.ifft_plan)
        for j in 1:pad, i in 1:pad
            value = out[i, j]
            for batch_index in 1:batch.batch_size
                value += abs2(stack[i, j, batch_index])
            end
            out[i, j] = value
        end
    end
    return nothing
end

function _propagate_pyramid_shifted_masks!(style::AcceleratorStyle, out,
    front_end::PyramidOpticalFrontEnd,
    batch::PyramidShiftedMaskModulationWorkspace)
    stack = batch.field_stack
    focal_field = pyramid_propagation_workspace(front_end).focal_field
    pad = size(out, 1)
    point_count = size(batch.shifted_masks, 3)
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel!(style, pyramid_shifted_mask_batch_kernel!, stack,
            focal_field, batch.shifted_masks, batch.operating_weights,
            first_point, pad, batch.batch_size; ndrange=size(stack))
        execute_fft_plan!(stack, batch.ifft_plan)
        launch_kernel!(style, pyramid_modulation_batch_intensity_kernel!,
            out, stack, pad, batch.batch_size; ndrange=size(out))
    end
    return nothing
end

@inline function _propagate_pyramid_modulation_batch!(out,
    front_end::PyramidOpticalFrontEnd,
    batch::PyramidModulationBatchWorkspace, style::AcceleratorStyle)
    stack = batch.field_stack
    pad = size(out, 1)
    execute_fft_plan!(stack, batch.fft_plan)
    launch_kernel!(style, pyramid_modulation_batch_mask_kernel!,
        stack, pyramid_propagation_workspace(front_end).pyramid_mask,
        pad, batch.batch_size; ndrange=size(stack))
    execute_fft_plan!(stack, batch.ifft_plan)
    launch_kernel!(style, pyramid_modulation_batch_intensity_kernel!,
        out, stack, pad, batch.batch_size; ndrange=size(out))
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
    @inbounds for first_point in 1:batch.batch_size:point_count
        launch_kernel!(style, pyramid_pupil_modulation_batch_kernel!,
            batch.field_stack, amplitude, opd, modulation.phases, weights,
            pyramid_propagation_workspace(front_end).phasor,
            amplitude_scale, opd_to_cycles, offset, resolution, first_point,
            pad, batch.batch_size; ndrange=size(batch.field_stack))
        _propagate_pyramid_modulation_batch!(out, front_end, batch, style)
    end
    return true
end

function _pyramid_pupil_modulation_batch!(
    batch::PyramidShiftedMaskModulationWorkspace, out,
    front_end::PyramidOpticalFrontEnd, amplitude, opd, modulation,
    amplitude_scale, opd_to_cycles, offset::Int, resolution::Int)
    modulation === front_end.modulation || return false
    style = execution_style(out)
    _prepare_pyramid_shifted_focal_field!(style, front_end, amplitude, opd,
        amplitude_scale, opd_to_cycles, offset, resolution)
    _propagate_pyramid_shifted_masks!(style, out, front_end, batch)
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
            batch.field_stack, field, modulation.phases, weights,
            pyramid_propagation_workspace(front_end).phasor,
            offset, resolution, first_point, pad, batch.batch_size;
            ndrange=size(batch.field_stack))
        _propagate_pyramid_modulation_batch!(out, front_end, batch, style)
    end
    return true
end

function _pyramid_electric_field_modulation_batch!(
    batch::PyramidShiftedMaskModulationWorkspace, out,
    front_end::PyramidOpticalFrontEnd, field, modulation,
    offset::Int, resolution::Int)
    modulation === front_end.modulation || return false
    style = execution_style(out)
    _prepare_pyramid_shifted_focal_field!(style, front_end, field, offset,
        resolution)
    _propagate_pyramid_shifted_masks!(style, out, front_end, batch)
    return true
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, pupil, src,
        pyramid_operating_modulation(wfs))
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource,
    modulation::PreparedFocalPlaneModulation) where {T<:AbstractFloat}
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
    @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * amplitude *
        cispi(opd_to_cycles * pupil.opd)
    for p in 1:modulation_point_count(modulation)
        copyto!(propagation.focal_field, propagation.field)
        @views @. propagation.focal_field[ox+1:ox+n, oy+1:oy+n] *=
            modulation.amplitude_weights[p] * modulation.phases[:, :, p]
        accumulate_pyramid_focal_intensity!(out, wfs)
    end
    return out
end

function pyramid_intensity_core!(::AcceleratorStyle,
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
    for p in 1:modulation_point_count(modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * amplitude *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * pupil.opd)
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
        amplitude_weight = modulation.amplitude_weights[p]
        @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] =
            amp_scale * amplitude_weight * amplitude *
            modulation.phases[:, :, p] * cispi(opd_to_cycles * pupil.opd)
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
function pyramid_modulation_frame(wfs::PyramidWFS{<:Diffractive},
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
        wfs.estimator.params.pupil_samples,
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
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        propagation.focal_field,
        propagation.fft_plan,
    )
    propagation.lgs_kernel_tag = tag
    return wfs
end
