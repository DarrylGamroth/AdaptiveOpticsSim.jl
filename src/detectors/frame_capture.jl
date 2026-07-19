function ensure_latent_buffer!(det::Detector)
    if size(det.state.latent_buffer) != size(det.state.frame)
        det.state.latent_buffer = similar(det.state.frame, size(det.state.frame)...)
        fill!(det.state.latent_buffer, zero(eltype(det.state.latent_buffer)))
    end
    return det.state.latent_buffer
end

apply_signal_defects!(::NullDetectorDefectModel, det::Detector, exposure_time::Real) = det.state.frame
apply_dark_defects!(::NullDetectorDefectModel, det::Detector, exposure_time::Real) = det.state.frame

function apply_signal_defects!(model::PixelResponseNonuniformity, det::Detector, exposure_time::Real)
    _require_detector_defect_shape(model, size(det.state.frame))
    det.state.frame .*= model.gain_map
    return det.state.frame
end

apply_dark_defects!(::PixelResponseNonuniformity, det::Detector, exposure_time::Real) = det.state.frame
apply_signal_defects!(::DarkSignalNonuniformity, det::Detector, exposure_time::Real) = det.state.frame

function apply_dark_defects!(model::DarkSignalNonuniformity, det::Detector, exposure_time::Real)
    _require_detector_defect_shape(model, size(det.state.frame))
    det.state.frame .+= model.dark_map .* exposure_time
    return det.state.frame
end

function apply_signal_defects!(model::BadPixelMask, det::Detector, exposure_time::Real)
    _require_detector_defect_shape(model, size(det.state.frame))
    throughput = model.throughput
    throughput == one(throughput) && return det.state.frame
    det.state.frame .= ifelse.(model.mask, throughput .* det.state.frame, det.state.frame)
    return det.state.frame
end

apply_dark_defects!(::BadPixelMask, det::Detector, exposure_time::Real) = det.state.frame

function apply_signal_defects!(model::CompositeDetectorDefectModel, det::Detector, exposure_time::Real)
    foreach(stage -> apply_signal_defects!(stage, det, exposure_time), model.stages)
    return det.state.frame
end

function apply_dark_defects!(model::CompositeDetectorDefectModel, det::Detector, exposure_time::Real)
    foreach(stage -> apply_dark_defects!(stage, det, exposure_time), model.stages)
    return det.state.frame
end

apply_frame_nonlinearity!(::NullFrameNonlinearity, det::Detector) = det.state.frame

function apply_frame_nonlinearity!(model::SaturatingFrameNonlinearity, det::Detector)
    coeff = model.coefficient
    coeff <= zero(coeff) && return det.state.frame
    @. det.state.frame = det.state.frame / (1 + coeff * det.state.frame)
    return det.state.frame
end

apply_sensor_persistence!(::FrameSensorType, det::Detector, exposure_time::Real) = det.state.frame
update_sensor_persistence!(::FrameSensorType, det::Detector, exposure_time::Real) = det.state.frame

function detector_host_buffer!(det::Detector, ::Type{T}, dims::Tuple{Int,Int}) where {T<:AbstractFloat}
    host = det.state.noise_buffer_host
    if size(host) != dims
        host = Matrix{T}(undef, dims...)
        det.state.noise_buffer_host = host
    end
    return host
end

function detector_host_frame!(det::Detector, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_buffer!(det, T, size(frame))
    copyto!(host, frame)
    return host
end

function _poisson_noise_frame!(::DetectorDirectPlan, det::Detector, rng::AbstractRNG, img::AbstractMatrix{T}) where {T<:AbstractFloat}
    poisson_noise!(rng, img)
    return img
end

function _poisson_noise_frame!(::DetectorDirectPlan, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    poisson_noise!(rng, cube)
    return cube
end

function _poisson_noise_frame!(::DetectorHostMirrorPlan, det::Detector, rng::AbstractRNG, img::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_frame!(det, img)
    _poisson_noise!(ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end

function _poisson_noise_frame!(plan::DetectorHostMirrorPlan, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _poisson_noise_frame!(plan, det, rng, @view(cube[b, :, :]))
    end
    return cube
end

function poisson_noise_frame!(det::Detector, rng::AbstractRNG, img::AbstractArray{T}) where {T<:AbstractFloat}
    plan = detector_execution_plan(typeof(execution_style(img)), typeof(det))
    return _poisson_noise_frame!(plan, det, rng, img)
end

function _randn_frame_noise!(::DetectorDirectPlan, det::Detector, rng::AbstractRNG, out::AbstractMatrix{T}) where {T<:AbstractFloat}
    randn_backend!(rng, out)
    return out
end

function _randn_frame_noise!(::DetectorDirectPlan, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    randn_backend!(rng, cube)
    return cube
end

function _randn_frame_noise!(::DetectorHostMirrorPlan, det::Detector, rng::AbstractRNG, out::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_buffer!(det, T, size(out))
    randn!(rng, host)
    copyto!(out, host)
    return out
end

function _randn_frame_noise!(plan::DetectorHostMirrorPlan, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _randn_frame_noise!(plan, det, rng, @view(cube[b, :, :]))
    end
    return cube
end

function randn_frame_noise!(det::Detector, rng::AbstractRNG, out::AbstractArray{T}) where {T<:AbstractFloat}
    plan = detector_execution_plan(typeof(execution_style(out)), typeof(det))
    return _randn_frame_noise!(plan, det, rng, out)
end

function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time)
    return nothing
end
function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time, qe)
    return nothing
end

function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time)
    return nothing
end
function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time, qe)
    return nothing
end

function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time)
    return nothing
end
function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time, qe)
    return nothing
end

function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time)
    return nothing
end
function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_time, qe)
    return nothing
end

apply_background_flux!(::NoBackground, det::Detector, rng::AbstractRNG, exposure_time::Real) = det.state.frame

function apply_background_flux!(background::ScalarBackground, det::Detector, rng::AbstractRNG, exposure_time::Real)
    return add_poisson_rate!(det.state.frame, det, rng, background.level * exposure_time)
end

function apply_background_flux!(background::BackgroundFrame, det::Detector, rng::AbstractRNG, exposure_time::Real)
    _require_background_flux_shape(background, size(det.state.frame))
    copyto!(det.state.noise_buffer, background.map)
    det.state.noise_buffer .*= exposure_time
    poisson_noise_frame!(det, rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_dark_current!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    dark_signal = effective_dark_current(det) * effective_dark_current_time(det.params.sensor, exposure_time)
    return add_poisson_rate!(det.state.frame, det, rng, dark_signal)
end

sensor_saturation_limit(det::Detector) = sensor_saturation_limit(det.params.sensor, det)
sensor_saturation_limit(::FrameSensorType, det::Detector) = det.params.full_well

function _apply_saturation!(::DetectorDirectPlan, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    clamp_array!(det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

function _apply_saturation!(::DetectorHostMirrorPlan, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    host = detector_host_frame!(det, det.state.frame)
    clamp!(host, zero(eltype(host)), full_well)
    copyto!(det.state.frame, host)
    return det.state.frame
end

function _apply_saturation!(::ScalarCPUStyle, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    clamp!(det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

function _apply_saturation!(style::AcceleratorStyle, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    _clamp_array!(style, det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

@inline _detector_value_plan(plan::DetectorDirectPlan, ::ExecutionStyle) = plan
@inline _detector_value_plan(plan::DetectorHostMirrorPlan, ::ScalarCPUStyle) = plan
@inline _detector_value_plan(::DetectorHostMirrorPlan, style::AcceleratorStyle) = style

function apply_saturation!(det::Detector)
    style = execution_style(det.state.frame)
    plan = detector_execution_plan(typeof(style), typeof(det))
    return _apply_saturation!(_detector_value_plan(plan, style), det)
end

apply_sensor_statistics!(sensor::FrameSensorType, det::Detector,
    rng::AbstractRNG, exposure_time::Real) = det.state.frame

function apply_avalanche_excess_noise!(factor, det::Detector, rng::AbstractRNG)
    factor <= one(factor) && return det.state.frame
    randn_backend!(rng, det.state.noise_buffer)
    scale2 = factor * factor - one(factor)
    zero_t = zero(eltype(det.state.frame))
    @. det.state.frame += sqrt(max(scale2 * det.state.frame, zero_t)) * det.state.noise_buffer
    return det.state.frame
end

apply_pre_readout_gain!(::FrameSensorType, det::Detector, rng::AbstractRNG) = det.state.frame
apply_post_readout_gain!(::FrameSensorType, det::Detector) = det.state.frame
apply_detection_output!(::FrameSensorType, det::Detector,
    rng::AbstractRNG) = det.state.frame
apply_charge_transfer!(::FrameSensorType, det::Detector) = det.state.frame
reset_readout_products!(det::Detector) = (det.state.readout_products = NoFrameReadoutProducts(); det)

apply_readout_noise!(det::Detector{NoiseNone}, rng::AbstractRNG) = det.state.frame
apply_readout_noise!(det::Detector{NoisePhoton}, rng::AbstractRNG) = det.state.frame

function apply_readout_noise!(det::Detector{<:NoiseReadout}, rng::AbstractRNG)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    return add_gaussian_noise!(det.state.frame, det, rng, sigma)
end

function apply_readout_noise!(det::Detector{<:NoisePhotonReadout}, rng::AbstractRNG)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    return add_gaussian_noise!(det.state.frame, det, rng, sigma)
end

apply_sensor_readout_noise!(::FrameSensorType, det::Detector,
    rng::AbstractRNG) = det.state.frame

function _apply_quantization!(::DetectorDirectPlan, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    levels = exp2(eltype(det.state.frame)(bits))
    full_well = something(det.params.full_well)
    det.state.frame .*= (levels - one(levels)) / full_well
    clamp_array!(det.state.frame, zero(eltype(det.state.frame)), levels - one(levels))
    return det.state.frame
end

function _apply_quantization!(::DetectorHostMirrorPlan, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    host = detector_host_frame!(det, det.state.frame)
    levels = exp2(eltype(host)(bits))
    full_well = something(det.params.full_well)
    host .*= (levels - one(levels)) / full_well
    clamp!(host, zero(eltype(host)), levels - one(levels))
    copyto!(det.state.frame, host)
    return det.state.frame
end

function _apply_quantization!(::ScalarCPUStyle, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    levels = exp2(eltype(det.state.frame)(bits))
    full_well = something(det.params.full_well)
    det.state.frame .*= (levels - one(levels)) / full_well
    clamp!(det.state.frame, zero(eltype(det.state.frame)), levels - one(levels))
    return det.state.frame
end

function _apply_quantization!(style::AcceleratorStyle, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    levels = exp2(eltype(det.state.frame)(bits))
    full_well = something(det.params.full_well)
    det.state.frame .*= (levels - one(levels)) / full_well
    _clamp_array!(style, det.state.frame, zero(eltype(det.state.frame)),
        levels - one(levels))
    return det.state.frame
end

function apply_quantization!(det::Detector)
    style = execution_style(det.state.frame)
    plan = detector_execution_plan(typeof(style), typeof(det))
    return _apply_quantization!(_detector_value_plan(plan, style), det)
end

subtract_background_map!(::NoBackground, det::Detector) = det.state.frame

function subtract_background_map!(background::ScalarBackground, det::Detector)
    det.state.frame .-= background.level
    return det.state.frame
end

function subtract_background_map!(background::BackgroundFrame, det::Detector)
    _require_background_map_shape(background, size(det.state.frame))
    det.state.frame .-= background.map
    return det.state.frame
end

readout_product_shape(det::Detector) = det.params.readout_window === nothing ?
    size(det.state.frame) :
    (length(det.params.readout_window.rows), length(det.params.readout_window.cols))

function _copy_windowed_frame(frame::AbstractMatrix, det::Detector)
    plan = detector_execution_plan(typeof(execution_style(frame)), typeof(det))
    return _copy_windowed_frame(plan, frame, det)
end

function _copy_windowed_frame(::DetectorDirectPlan, frame::AbstractMatrix,
    det::Detector)
    window = det.params.readout_window
    if window === nothing
        return copy(frame)
    end
    return copy(@view(frame[window.rows, window.cols]))
end

function _copy_windowed_cube(cube::AbstractArray{T,3}, det::Detector) where {T}
    plan = detector_execution_plan(typeof(execution_style(cube)), typeof(det))
    return _copy_windowed_cube(plan, cube, det)
end

function _copy_windowed_cube(::DetectorDirectPlan, cube::AbstractArray{T,3},
    det::Detector) where {T}
    window = det.params.readout_window
    if window === nothing
        return copy(cube)
    end
    return copy(@view(cube[window.rows, window.cols, :]))
end

function _copy_windowed_frame(::DetectorHostMirrorPlan, frame::AbstractMatrix,
    det::Detector)
    frame_host = Array(frame)
    window = det.params.readout_window
    selected = window === nothing ? frame_host :
        copy(@view(frame_host[window.rows, window.cols]))
    output = similar(frame, size(selected)...)
    copyto!(output, selected)
    return output
end

function _copy_windowed_cube(::DetectorHostMirrorPlan, cube::AbstractArray{T,3},
    det::Detector) where {T}
    cube_host = Array(cube)
    window = det.params.readout_window
    selected = window === nothing ? cube_host :
        copy(@view(cube_host[window.rows, window.cols, :]))
    output = similar(cube, size(selected)...)
    copyto!(output, selected)
    return output
end

apply_readout_correction!(::NullFrameReadoutCorrection, frame::AbstractMatrix) = frame

function _reference_pixel_bias(model::ReferencePixelCommonModeCorrection, frame::AbstractMatrix{T}) where {T}
    n, m = size(frame)
    n_edge_rows = min(model.edge_rows, n)
    n_edge_cols = min(model.edge_cols, m)
    total = zero(T)
    count = 0
    if n_edge_rows > 0
        total += sum(@view(frame[1:n_edge_rows, :]))
        count += n_edge_rows * m
        if n > n_edge_rows
            total += sum(@view(frame[n - n_edge_rows + 1:n, :]))
            count += n_edge_rows * m
        end
    end
    row_lo = n_edge_rows + 1
    row_hi = n - n_edge_rows
    if n_edge_cols > 0 && row_lo <= row_hi
        total += sum(@view(frame[row_lo:row_hi, 1:n_edge_cols]))
        count += (row_hi - row_lo + 1) * n_edge_cols
        if m > n_edge_cols
            total += sum(@view(frame[row_lo:row_hi, m - n_edge_cols + 1:m]))
            count += (row_hi - row_lo + 1) * n_edge_cols
        end
    end
    count > 0 || return zero(T)
    return total / T(count)
end

function _row_reference_bias(edge_cols::Int, row::AbstractVector{T}) where {T}
    n = length(row)
    n_edge = min(edge_cols, n)
    total = zero(T)
    count = 0
    if n_edge > 0
        total += sum(@view(row[1:n_edge]))
        count += n_edge
        if n > n_edge
            total += sum(@view(row[n - n_edge + 1:n]))
            count += n_edge
        end
    end
    count > 0 || return zero(T)
    return total / T(count)
end

function _column_reference_bias(edge_rows::Int, col::AbstractVector{T}) where {T}
    n = length(col)
    n_edge = min(edge_rows, n)
    total = zero(T)
    count = 0
    if n_edge > 0
        total += sum(@view(col[1:n_edge]))
        count += n_edge
        if n > n_edge
            total += sum(@view(col[n - n_edge + 1:n]))
            count += n_edge
        end
    end
    count > 0 || return zero(T)
    return total / T(count)
end

function _output_reference_bias(model::ReferenceOutputCommonModeCorrection, block::AbstractMatrix{T}) where {T}
    n, m = size(block)
    n_edge_rows = min(model.edge_rows, n)
    n_edge_cols = min(model.edge_cols, m)
    total = zero(T)
    count = 0
    if n_edge_rows > 0
        total += sum(@view(block[1:n_edge_rows, :]))
        count += n_edge_rows * m
        if n > n_edge_rows
            total += sum(@view(block[n - n_edge_rows + 1:n, :]))
            count += n_edge_rows * m
        end
    end
    row_lo = n_edge_rows + 1
    row_hi = n - n_edge_rows
    if n_edge_cols > 0 && row_lo <= row_hi
        total += sum(@view(block[row_lo:row_hi, 1:n_edge_cols]))
        count += (row_hi - row_lo + 1) * n_edge_cols
        if m > n_edge_cols
            total += sum(@view(block[row_lo:row_hi, m - n_edge_cols + 1:m]))
            count += (row_hi - row_lo + 1) * n_edge_cols
        end
    end
    count > 0 || return zero(T)
    return total / T(count)
end

function apply_readout_correction!(model::ReferencePixelCommonModeCorrection, frame::AbstractMatrix)
    frame .-= _reference_pixel_bias(model, frame)
    return frame
end

function apply_readout_correction!(model::ReferenceRowCommonModeCorrection, frame::AbstractMatrix)
    for row_idx in axes(frame, 1)
        row = @view(frame[row_idx, :])
        row .-= _row_reference_bias(model.edge_cols, row)
    end
    return frame
end

function apply_readout_correction!(model::ReferenceColumnCommonModeCorrection, frame::AbstractMatrix)
    for col_idx in axes(frame, 2)
        col = @view(frame[:, col_idx])
        col .-= _column_reference_bias(model.edge_rows, col)
    end
    return frame
end

function apply_readout_correction!(model::ReferenceOutputCommonModeCorrection, frame::AbstractMatrix)
    n_cols = size(frame, 2)
    for col_lo in 1:model.output_cols:n_cols
        col_hi = min(col_lo + model.output_cols - 1, n_cols)
        block = @view(frame[:, col_lo:col_hi])
        block .-= _output_reference_bias(model, block)
    end
    return frame
end

function apply_readout_correction!(model::CompositeFrameReadoutCorrection, frame::AbstractMatrix)
    for stage in model.stages
        apply_readout_correction!(stage, frame)
    end
    return frame
end

function apply_readout_correction!(model::FrameReadoutCorrectionModel, cube::AbstractArray{T,3}) where {T}
    for k in axes(cube, 3)
        apply_readout_correction!(model, @view(cube[:, :, k]))
    end
    return cube
end

function _apply_readout_correction!(::DetectorDirectPlan, model::FrameReadoutCorrectionModel,
    frame::AbstractMatrix{T}, det::Detector) where {T<:AbstractFloat}
    return _apply_readout_correction!(execution_style(frame), model, frame, det)
end

function _apply_readout_correction!(::DetectorDirectPlan, model::FrameReadoutCorrectionModel,
    cube::AbstractArray{T,3}, det::Detector) where {T<:AbstractFloat}
    return apply_readout_correction!(model, cube)
end

function _detector_readout_scratch!(det::Detector, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    if size(det.state.noise_buffer) != size(frame)
        det.state.noise_buffer = similar(det.state.noise_buffer, size(frame)...)
        fill!(det.state.noise_buffer, zero(eltype(det.state.noise_buffer)))
    end
    return det.state.noise_buffer
end

function _apply_readout_correction!(::ScalarCPUStyle, model::FrameReadoutCorrectionModel,
    frame::AbstractMatrix{T}, det::Detector) where {T<:AbstractFloat}
    return apply_readout_correction!(model, frame)
end

function _apply_readout_correction!(style::AcceleratorStyle, model::FrameReadoutCorrectionModel,
    frame::AbstractMatrix{T}, det::Detector) where {T<:AbstractFloat}
    if supports_batched_readout_correction(model)
        scratch = _detector_readout_scratch!(det, frame)
        frame_cube = reshape(frame, 1, size(frame, 1), size(frame, 2))
        scratch_cube = reshape(scratch, 1, size(scratch, 1), size(scratch, 2))
        _batched_apply_readout_correction!(style, model, frame_cube, scratch_cube)
        return frame
    end
    return apply_readout_correction!(model, frame)
end

function _apply_readout_correction!(::DetectorHostMirrorPlan, model::FrameReadoutCorrectionModel,
    frame::AbstractMatrix{T}, det::Detector) where {T<:AbstractFloat}
    style = execution_style(frame)
    if can_apply_device_readout_correction(style, model)
        return _apply_readout_correction!(style, model, frame, det)
    end
    host = detector_host_frame!(det, frame)
    apply_readout_correction!(model, host)
    copyto!(frame, host)
    return frame
end

can_apply_device_readout_correction(::ScalarCPUStyle, ::FrameReadoutCorrectionModel) = false
can_apply_device_readout_correction(::AcceleratorStyle, model::FrameReadoutCorrectionModel) =
    supports_batched_readout_correction(model)

function _apply_readout_correction!(plan::DetectorHostMirrorPlan, model::FrameReadoutCorrectionModel,
    cube::AbstractArray{T,3}, det::Detector) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _apply_readout_correction!(plan, model, @view(cube[b, :, :]), det)
    end
    return cube
end

function apply_readout_correction!(model::FrameReadoutCorrectionModel, frame::AbstractArray{T}, det::Detector) where {T<:AbstractFloat}
    plan = detector_execution_plan(typeof(execution_style(frame)), typeof(det))
    return _apply_readout_correction!(plan, model, frame, det)
end

finalize_readout_products!(::FrameSensorType, det::Detector, rng::AbstractRNG, exposure_time::Real) =
    reset_readout_products!(det)

function _finalize_capture!(::FrameSensorType, det::Detector, rng::AbstractRNG,
    exposure_time::Real)
    return finalize_readout_pipeline!(det, rng, exposure_time)
end

function finalize_capture!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    return _finalize_capture!(det.params.sensor, det, rng, exposure_time)
end

function _finalize_incremental_capture!(::FrameSensorType, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_readout_pipeline!(det, rng, exposure_time,
        zero(exposure_time))
end

function finalize_incremental_capture!(det::Detector, rng::AbstractRNG,
    exposure_time::Real)
    return _finalize_incremental_capture!(det.params.sensor, det, rng,
        exposure_time)
end

validated_temporal_frame(frame::AbstractMatrix) = frame
validated_temporal_frame(frame) =
    throw(InvalidConfiguration("FunctionFrameSource must return an AbstractMatrix"))

function initial_temporal_frame(source::FunctionFrameSource, det::Detector,
    time, exposure_time)
    return validated_temporal_frame(source.f(time))
end

function initial_temporal_frame(source::InPlaceFrameSource, det::Detector,
    time, exposure_time)
    frame = ensure_temporal_buffer!(det, source.frame_size)
    sample_frame!(frame, source, time)
    return frame
end

function initial_temporal_frame(source::FunctionExposureFrameSource,
    det::Detector, time, exposure_time)
    return validated_temporal_frame(source.f(time, exposure_time))
end

function initial_temporal_frame(source::InPlaceExposureFrameSource,
    det::Detector, time, exposure_time)
    frame = ensure_temporal_buffer!(det, source.frame_size)
    sample_exposure_frame!(frame, source, time, exposure_time)
    return frame
end

function ensure_temporal_buffer!(det::Detector, dims::Tuple{Int,Int})
    if size(det.state.temporal_buffer) != dims
        det.state.temporal_buffer = similar(det.state.temporal_buffer, dims...)
    end
    return det.state.temporal_buffer
end

function capture_temporal_signal!(det::Detector, source::AbstractTemporalFrameSource, first_frame::AbstractMatrix,
    rng::AbstractRNG, exposure_time::Real, ::GlobalShutter)
    capture_signal_pipeline!(det, first_frame, rng, exposure_time)
    return det.state.frame
end

rolling_exposure_start(::RollingExposure, line_index, line_time, exposure_time, ::Type{T}) where {T<:AbstractFloat} =
    T(line_index) * T(line_time)
rolling_exposure_duration(::RollingExposure, line_index, line_time, exposure_time, ::Type{T}) where {T<:AbstractFloat} =
    T(exposure_time)
rolling_exposure_start(::GlobalResetExposure, line_index, line_time, exposure_time, ::Type{T}) where {T<:AbstractFloat} =
    zero(T)
rolling_exposure_duration(::GlobalResetExposure, line_index, line_time, exposure_time, ::Type{T}) where {T<:AbstractFloat} =
    T(exposure_time) + T(line_index) * T(line_time)

function capture_temporal_signal!(det::Detector, source::AbstractTemporalFrameSource, first_frame::AbstractMatrix,
    rng::AbstractRNG, exposure_time::Real, timing::RollingShutter)
    fill_frame!(det, first_frame, exposure_time)
    det.state.accum_buffer .= det.state.frame

    scratch = ensure_temporal_buffer!(det, size(first_frame))
    n_rows = size(det.state.frame, 1)
    group_size = timing.row_group_size
    value_type = eltype(det.state.frame)
    for row_lo in (firstindex(det.state.frame, 1) + group_size):group_size:n_rows
        row_hi = min(row_lo + group_size - 1, n_rows)
        line_index = div(row_lo - 1, group_size)
        sample_time = rolling_exposure_start(timing.exposure_mode, line_index, timing.line_time, exposure_time, value_type)
        group_exposure = rolling_exposure_duration(timing.exposure_mode, line_index, timing.line_time, exposure_time, value_type)
        sample_exposure_frame!(scratch, source, sample_time, group_exposure)
        fill_frame!(det, scratch, group_exposure)
        @views det.state.accum_buffer[row_lo:row_hi, :] .= det.state.frame[row_lo:row_hi, :]
    end

    det.state.frame .= det.state.accum_buffer
    apply_signal_defects!(det.params.defect_model, det, exposure_time)
    apply_sensor_persistence!(det.params.sensor, det, exposure_time)
    photon_noise_enabled(det) && poisson_noise_frame!(det, rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return det.state.frame
end

function _write_output!(::DetectorDirectPlan, det::Detector, output::AbstractMatrix,
    source::AbstractMatrix)
    if eltype(output) <: Integer
        write_integer_output!(output, source)
    else
        copyto!(output, source)
    end
    return output
end

function _write_output!(::DetectorHostMirrorPlan, det::Detector,
    output::AbstractMatrix, source::AbstractMatrix)
    frame_host = detector_host_frame!(det, det.state.frame)
    window = det.params.readout_window
    source_host = window === nothing ? frame_host :
        @view(frame_host[window.rows, window.cols])
    output_host = det.state.output_buffer_host
    output_host === nothing && throw(InvalidConfiguration(
        "Detector host-mirror output requires a host output buffer"))
    if eltype(output_host) <: Integer
        _write_integer_output!(ScalarCPUStyle(), output_host, source_host)
    else
        copyto!(output_host, source_host)
    end
    copyto!(output, output_host)
    return output
end

function write_output!(det::Detector)
    window = det.params.readout_window
    output = det.state.output_buffer
    if output === nothing
        window === nothing && return det.state.frame
        throw(InvalidConfiguration("Detector readout_window requires an allocated output buffer"))
    end
    source = window === nothing ? det.state.frame : @view(det.state.frame[window.rows, window.cols])
    plan = detector_execution_plan(typeof(execution_style(output)), typeof(det))
    return _write_output!(plan, det, output, source)
end

@inline function require_whole_capture_idle(det::Detector)
    iszero(det.state.integrated_time) && det.state.readout_ready ||
        throw(InvalidConfiguration(
            "cannot start a whole detector exposure while an incremental " *
            "exposure is pending; complete it with sample_duration or call " *
            "reset_integration!"))
    return nothing
end

@inline function completed_whole_capture_output!(det::Detector)
    output = write_output!(det)
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = true
    return output
end

function capture_with_quantum_efficiency!(det::Detector,
    photon_rate::AbstractMatrix{T}, quantum_efficiency::Real,
    rng::AbstractRNG) where {T}
    require_whole_capture_idle(det)
    exposure_time = det.params.integration_time
    capture_signal!(det, photon_rate, rng, exposure_time,
        quantum_efficiency)
    finalize_capture!(det, rng, exposure_time)
    advance_thermal!(det, exposure_time)
    return completed_whole_capture_output!(det)
end

"""
    capture!(detector, photon_rate, rng)

Legacy matrix acquisition path. Each matrix value is interpreted as a
cell-integrated photon-arrival rate on the input optical grid. The detector
applies its configured exposure duration exactly once. Use
`prepare_detector_acquisition` with an `IntensityMap` when geometry,
radiometry, backend, and device contracts must be checked explicitly.
"""
function capture!(det::Detector, psf::AbstractMatrix{T},
    rng::AbstractRNG) where {T}
    return capture_with_quantum_efficiency!(det, psf, det.params.qe, rng)
end

function capture!(det::Detector, psf::AbstractMatrix{T}, src::AbstractSource, rng::AbstractRNG) where {T}
    require_whole_capture_idle(det)
    return capture_with_quantum_efficiency!(det, psf,
        effective_qe(det, src, eltype(det.state.frame)), rng)
end

function capture!(det::Detector, psf::AbstractMatrix{T}, src::AbstractSource; rng::AbstractRNG=Random.default_rng()) where {T}
    return capture!(det, psf, src, rng)
end

function capture!(det::Detector, source::AbstractTemporalFrameSource, rng::AbstractRNG)
    require_whole_capture_idle(det)
    exposure_time = det.params.integration_time
    first_frame = initial_temporal_frame(source, det,
        zero(eltype(det.state.frame)), exposure_time)
    capture_temporal_signal!(det, source, first_frame, rng, exposure_time,
        det.params.timing_model)
    finalize_capture!(det, rng, exposure_time)
    advance_thermal!(det, exposure_time)
    return completed_whole_capture_output!(det)
end

function capture!(det::Detector, source::AbstractTemporalFrameSource; rng::AbstractRNG=Random.default_rng())
    return capture!(det, source, rng)
end

"""
    capture_incremental!(detector, photon_rate, rng, sample_duration,
        quantum_efficiency=detector.params.qe)

Accumulate one positive `sample_duration` in seconds from a cell-integrated
photon-arrival-rate matrix. This frame-step convenience finalizes automatically
when the configured integration duration is reached. `sample_duration` is not
an absolute timestamp; scheduled detector events own their timestamps and
completion semantics separately.
"""
function capture_incremental!(det::Detector, photon_rate::AbstractMatrix,
    rng::AbstractRNG, sample_duration::Real, qe=det.params.qe)
    if !iszero(det.state.integrated_time) || !det.state.readout_ready
        size(photon_rate) == size(det.state.presampling_buffer) ||
            throw(DimensionMismatchError(
                "incremental detector input dimensions cannot change while " *
                "an exposure is pending"))
    end
    prepare_detector_buffers!(det, size(photon_rate))
    T = eltype(det.state.frame)
    dt = T(sample_duration)
    isfinite(dt) && dt > zero(T) || throw(InvalidConfiguration(
        "sample_duration must be finite and > 0"))
    remaining = det.params.integration_time - det.state.integrated_time
    tolerance = T(8) * eps(det.params.integration_time) *
        max(one(T), abs(det.params.integration_time))
    dt <= remaining + tolerance || throw(InvalidConfiguration(
        "sample_duration exceeds the remaining detector integration duration"))
    dt = min(dt, remaining)

    exposure_start = iszero(det.state.integrated_time)
    exposure_start && fill!(det.state.accum_buffer,
        zero(eltype(det.state.accum_buffer)))
    capture_signal_pipeline!(det, photon_rate, rng, dt, qe, exposure_start,
        det.params.integration_time)
    accumulate_incremental_charge_generation!(det, rng, dt)
    det.state.accum_buffer .+= det.state.frame
    det.state.integrated_time += dt
    advance_thermal!(det, dt)
    det.state.readout_ready = false
    if det.state.integrated_time + tolerance >= det.params.integration_time
        det.state.frame .= det.state.accum_buffer
        finalize_incremental_capture!(det, rng, det.state.integrated_time)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    return write_output!(det)
end

function capture!(det::Detector, photon_rate::AbstractMatrix{T};
    rng::AbstractRNG=Random.default_rng(),
    sample_duration::Union{Nothing,Real}=nothing) where {T}
    if sample_duration === nothing
        return capture!(det, photon_rate, rng)
    end
    return capture_incremental!(det, photon_rate, rng, sample_duration)
end
