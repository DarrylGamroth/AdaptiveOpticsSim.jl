function ensure_latent_buffer!(det::Detector)
    if size(det.state.latent_buffer) != size(det.products.frame)
        det.state.latent_buffer = similar(det.products.frame, size(det.products.frame)...)
        fill!(det.state.latent_buffer, zero(eltype(det.state.latent_buffer)))
    end
    return det.state.latent_buffer
end

apply_signal_defects!(::NullDetectorDefectModel, det::Detector, exposure_duration::Real) = det.products.frame
apply_dark_defects!(::NullDetectorDefectModel, det::Detector, exposure_duration::Real) = det.products.frame

function apply_signal_defects!(model::PixelResponseNonuniformity, det::Detector, exposure_duration::Real)
    _require_detector_defect_shape(model, size(det.products.frame))
    det.products.frame .*= model.gain_map
    return det.products.frame
end

apply_dark_defects!(::PixelResponseNonuniformity, det::Detector, exposure_duration::Real) = det.products.frame
apply_signal_defects!(::DarkSignalNonuniformity, det::Detector, exposure_duration::Real) = det.products.frame

function apply_dark_defects!(model::DarkSignalNonuniformity, det::Detector, exposure_duration::Real)
    _require_detector_defect_shape(model, size(det.products.frame))
    det.products.frame .+= model.dark_map .* exposure_duration
    return det.products.frame
end

function apply_signal_defects!(model::BadPixelMask, det::Detector, exposure_duration::Real)
    _require_detector_defect_shape(model, size(det.products.frame))
    throughput = model.throughput
    throughput == one(throughput) && return det.products.frame
    det.products.frame .= ifelse.(model.mask, throughput .* det.products.frame, det.products.frame)
    return det.products.frame
end

apply_dark_defects!(::BadPixelMask, det::Detector, exposure_duration::Real) = det.products.frame

function apply_signal_defects!(model::CompositeDetectorDefectModel, det::Detector, exposure_duration::Real)
    foreach(stage -> apply_signal_defects!(stage, det, exposure_duration), model.stages)
    return det.products.frame
end

function apply_dark_defects!(model::CompositeDetectorDefectModel, det::Detector, exposure_duration::Real)
    foreach(stage -> apply_dark_defects!(stage, det, exposure_duration), model.stages)
    return det.products.frame
end

apply_frame_nonlinearity!(::NullFrameNonlinearity, det::Detector) = det.products.frame

function apply_frame_nonlinearity!(model::SaturatingFrameNonlinearity, det::Detector)
    coeff = model.coefficient
    coeff <= zero(coeff) && return det.products.frame
    @. det.products.frame = det.products.frame / (1 + coeff * det.products.frame)
    return det.products.frame
end

apply_sensor_persistence!(::AbstractFrameSensor, det::Detector, exposure_duration::Real) = det.products.frame
update_sensor_persistence!(::AbstractFrameSensor, det::Detector, exposure_duration::Real) = det.products.frame

function detector_host_buffer!(det::Detector, ::Type{T}, dims::Tuple{Int,Int}) where {T<:AbstractFloat}
    host = det.workspace.noise_buffer_host
    if size(host) != dims
        host = Matrix{T}(undef, dims...)
        det.workspace.noise_buffer_host = host
    end
    return host
end

function detector_host_frame!(det::Detector, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_buffer!(det, T, size(frame))
    copyto!(host, frame)
    return host
end

function _poisson_noise_frame!(::DetectorDirectStrategy, det::Detector, rng::AbstractRNG, img::AbstractMatrix{T}) where {T<:AbstractFloat}
    poisson_noise!(rng, img)
    return img
end

function _poisson_noise_frame!(::DetectorDirectStrategy, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    poisson_noise!(rng, cube)
    return cube
end

function _poisson_noise_frame!(::DetectorHostMirrorStrategy, det::Detector, rng::AbstractRNG, img::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_frame!(det, img)
    _poisson_noise!(ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end

function _poisson_noise_frame!(strategy::DetectorHostMirrorStrategy, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _poisson_noise_frame!(strategy, det, rng, @view(cube[b, :, :]))
    end
    return cube
end

function poisson_noise_frame!(det::Detector, rng::AbstractRNG, img::AbstractArray{T}) where {T<:AbstractFloat}
    strategy = detector_execution_strategy(typeof(execution_style(img)), typeof(det))
    return _poisson_noise_frame!(strategy, det, rng, img)
end

function _randn_frame_noise!(::DetectorDirectStrategy, det::Detector, rng::AbstractRNG, out::AbstractMatrix{T}) where {T<:AbstractFloat}
    randn_backend!(rng, out)
    return out
end

function _randn_frame_noise!(::DetectorDirectStrategy, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    randn_backend!(rng, cube)
    return cube
end

function _randn_frame_noise!(::DetectorHostMirrorStrategy, det::Detector, rng::AbstractRNG, out::AbstractMatrix{T}) where {T<:AbstractFloat}
    host = detector_host_buffer!(det, T, size(out))
    randn!(rng, host)
    copyto!(out, host)
    return out
end

function _randn_frame_noise!(strategy::DetectorHostMirrorStrategy, det::Detector, rng::AbstractRNG, cube::AbstractArray{T,3}) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _randn_frame_noise!(strategy, det, rng, @view(cube[b, :, :]))
    end
    return cube
end

function randn_frame_noise!(det::Detector, rng::AbstractRNG, out::AbstractArray{T}) where {T<:AbstractFloat}
    strategy = detector_execution_strategy(typeof(execution_style(out)), typeof(det))
    return _randn_frame_noise!(strategy, det, rng, out)
end

function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration)
    return nothing
end
function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration, qe)
    return nothing
end

function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration)
    return nothing
end
function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration, qe)
    return nothing
end

function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration)
    return nothing
end
function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration, qe)
    return nothing
end

function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration)
    return nothing
end
function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_duration::Real, qe) where {T}
    capture_signal_pipeline!(det, psf, rng, exposure_duration, qe)
    return nothing
end

apply_background_flux!(::NoBackground, det::Detector, rng::AbstractRNG, exposure_duration::Real) = det.products.frame

function apply_background_flux!(background::ScalarBackground, det::Detector, rng::AbstractRNG, exposure_duration::Real)
    return add_poisson_rate!(det.products.frame, det, rng, background.level * exposure_duration)
end

function apply_background_flux!(background::BackgroundFrame, det::Detector, rng::AbstractRNG, exposure_duration::Real)
    _require_background_flux_shape(background, size(det.products.frame))
    copyto!(det.workspace.noise_buffer, background.map)
    det.workspace.noise_buffer .*= exposure_duration
    poisson_noise_frame!(det, rng, det.workspace.noise_buffer)
    det.products.frame .+= det.workspace.noise_buffer
    return det.products.frame
end

function apply_dark_current!(det::Detector, rng::AbstractRNG, exposure_duration::Real)
    dark_signal = effective_dark_current(det) * effective_dark_current_duration(det.params.sensor, exposure_duration)
    return add_poisson_rate!(det.products.frame, det, rng, dark_signal)
end

sensor_saturation_limit(det::Detector) = sensor_saturation_limit(det.params.sensor, det)
sensor_saturation_limit(::AbstractFrameSensor, det::Detector) = det.params.full_well

function _apply_saturation!(::DetectorDirectStrategy, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.products.frame
    clamp_array!(det.products.frame, zero(eltype(det.products.frame)), full_well)
    return det.products.frame
end

function _apply_saturation!(::DetectorHostMirrorStrategy, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.products.frame
    host = detector_host_frame!(det, det.products.frame)
    clamp!(host, zero(eltype(host)), full_well)
    copyto!(det.products.frame, host)
    return det.products.frame
end

function _apply_saturation!(::ScalarCPUStyle, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.products.frame
    clamp!(det.products.frame, zero(eltype(det.products.frame)), full_well)
    return det.products.frame
end

function _apply_saturation!(style::AcceleratorStyle, det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.products.frame
    _clamp_array!(style, det.products.frame, zero(eltype(det.products.frame)), full_well)
    return det.products.frame
end

@inline _detector_value_strategy(strategy::DetectorDirectStrategy, ::ExecutionStyle) = strategy
@inline _detector_value_strategy(strategy::DetectorHostMirrorStrategy, ::ScalarCPUStyle) = strategy
@inline _detector_value_strategy(::DetectorHostMirrorStrategy, style::AcceleratorStyle) = style

function apply_saturation!(det::Detector)
    style = execution_style(det.products.frame)
    strategy = detector_execution_strategy(typeof(style), typeof(det))
    return _apply_saturation!(_detector_value_strategy(strategy, style), det)
end

apply_sensor_statistics!(sensor::AbstractFrameSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real) = det.products.frame

apply_pre_readout_gain!(::AbstractFrameSensor, det::Detector, rng::AbstractRNG) = det.products.frame
apply_post_readout_gain!(::AbstractFrameSensor, det::Detector) = det.products.frame
apply_detection_output!(::AbstractFrameSensor, det::Detector,
    rng::AbstractRNG) = det.products.frame
apply_charge_transfer!(::AbstractFrameSensor, det::Detector) = det.products.frame
reset_readout_products!(det::Detector) = (det.products.readout = NoFrameReadoutProducts(); det)

apply_readout_noise!(det::Detector{NoiseNone}, rng::AbstractRNG) = det.products.frame
apply_readout_noise!(det::Detector{NoisePhoton}, rng::AbstractRNG) = det.products.frame

function apply_readout_noise!(det::Detector{<:NoiseReadout}, rng::AbstractRNG)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    return add_gaussian_noise!(det.products.frame, det, rng, sigma)
end

function apply_readout_noise!(det::Detector{<:NoisePhotonReadout}, rng::AbstractRNG)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    return add_gaussian_noise!(det.products.frame, det, rng, sigma)
end

apply_sensor_readout_noise!(::AbstractFrameSensor, det::Detector,
    rng::AbstractRNG) = det.products.frame

function _apply_quantization!(::DetectorDirectStrategy, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.products.frame
    levels = exp2(eltype(det.products.frame)(bits))
    full_well = something(det.params.full_well)
    det.products.frame .*= (levels - one(levels)) / full_well
    clamp_array!(det.products.frame, zero(eltype(det.products.frame)), levels - one(levels))
    return det.products.frame
end

function _apply_quantization!(::DetectorHostMirrorStrategy, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.products.frame
    host = detector_host_frame!(det, det.products.frame)
    levels = exp2(eltype(host)(bits))
    full_well = something(det.params.full_well)
    host .*= (levels - one(levels)) / full_well
    clamp!(host, zero(eltype(host)), levels - one(levels))
    copyto!(det.products.frame, host)
    return det.products.frame
end

function _apply_quantization!(::ScalarCPUStyle, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.products.frame
    levels = exp2(eltype(det.products.frame)(bits))
    full_well = something(det.params.full_well)
    det.products.frame .*= (levels - one(levels)) / full_well
    clamp!(det.products.frame, zero(eltype(det.products.frame)), levels - one(levels))
    return det.products.frame
end

function _apply_quantization!(style::AcceleratorStyle, det::Detector)
    bits = det.params.bits
    bits === nothing && return det.products.frame
    levels = exp2(eltype(det.products.frame)(bits))
    full_well = something(det.params.full_well)
    det.products.frame .*= (levels - one(levels)) / full_well
    _clamp_array!(style, det.products.frame, zero(eltype(det.products.frame)),
        levels - one(levels))
    return det.products.frame
end

function apply_quantization!(det::Detector)
    style = execution_style(det.products.frame)
    strategy = detector_execution_strategy(typeof(style), typeof(det))
    return _apply_quantization!(_detector_value_strategy(strategy, style), det)
end

subtract_background_map!(::NoBackground, det::Detector) = det.products.frame

function subtract_background_map!(background::ScalarBackground, det::Detector)
    det.products.frame .-= background.level
    return det.products.frame
end

function subtract_background_map!(background::BackgroundFrame, det::Detector)
    _require_background_map_shape(background, size(det.products.frame))
    det.products.frame .-= background.map
    return det.products.frame
end

readout_product_shape(det::Detector) = det.params.readout_window === nothing ?
    size(det.products.frame) :
    (length(det.params.readout_window.rows), length(det.params.readout_window.cols))

function _copy_windowed_frame(frame::AbstractMatrix, det::Detector)
    strategy = detector_execution_strategy(typeof(execution_style(frame)), typeof(det))
    return _copy_windowed_frame(strategy, frame, det)
end

function _copy_windowed_frame(::DetectorDirectStrategy, frame::AbstractMatrix,
    det::Detector)
    window = det.params.readout_window
    if window === nothing
        return copy(frame)
    end
    return copy(@view(frame[window.rows, window.cols]))
end

function _copy_windowed_cube(cube::AbstractArray{T,3}, det::Detector) where {T}
    strategy = detector_execution_strategy(typeof(execution_style(cube)), typeof(det))
    return _copy_windowed_cube(strategy, cube, det)
end

function _copy_windowed_cube(::DetectorDirectStrategy, cube::AbstractArray{T,3},
    det::Detector) where {T}
    window = det.params.readout_window
    if window === nothing
        return copy(cube)
    end
    return copy(@view(cube[window.rows, window.cols, :]))
end

function _copy_windowed_frame(::DetectorHostMirrorStrategy, frame::AbstractMatrix,
    det::Detector)
    frame_host = Array(frame)
    window = det.params.readout_window
    selected = window === nothing ? frame_host :
        copy(@view(frame_host[window.rows, window.cols]))
    output = similar(frame, size(selected)...)
    copyto!(output, selected)
    return output
end

function _copy_windowed_cube(::DetectorHostMirrorStrategy, cube::AbstractArray{T,3},
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

function _apply_readout_correction!(::DetectorDirectStrategy, model::FrameReadoutCorrectionModel,
    frame::AbstractMatrix{T}, det::Detector) where {T<:AbstractFloat}
    return _apply_readout_correction!(execution_style(frame), model, frame, det)
end

function _apply_readout_correction!(::DetectorDirectStrategy, model::FrameReadoutCorrectionModel,
    cube::AbstractArray{T,3}, det::Detector) where {T<:AbstractFloat}
    return apply_readout_correction!(model, cube)
end

function _detector_readout_scratch!(det::Detector, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    if size(det.workspace.noise_buffer) != size(frame)
        det.workspace.noise_buffer = similar(det.workspace.noise_buffer, size(frame)...)
        fill!(det.workspace.noise_buffer, zero(eltype(det.workspace.noise_buffer)))
    end
    return det.workspace.noise_buffer
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

function _apply_readout_correction!(::DetectorHostMirrorStrategy, model::FrameReadoutCorrectionModel,
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

function _apply_readout_correction!(strategy::DetectorHostMirrorStrategy, model::FrameReadoutCorrectionModel,
    cube::AbstractArray{T,3}, det::Detector) where {T<:AbstractFloat}
    for b in axes(cube, 1)
        _apply_readout_correction!(strategy, model, @view(cube[b, :, :]), det)
    end
    return cube
end

function apply_readout_correction!(model::FrameReadoutCorrectionModel, frame::AbstractArray{T}, det::Detector) where {T<:AbstractFloat}
    strategy = detector_execution_strategy(typeof(execution_style(frame)), typeof(det))
    return _apply_readout_correction!(strategy, model, frame, det)
end

finalize_readout_products!(::AbstractFrameSensor, det::Detector, rng::AbstractRNG, exposure_duration::Real) =
    reset_readout_products!(det)

function _finalize_capture!(::AbstractFrameSensor, det::Detector, rng::AbstractRNG,
    exposure_duration::Real)
    return finalize_readout_pipeline!(det, rng, exposure_duration)
end

function finalize_capture!(det::Detector, rng::AbstractRNG, exposure_duration::Real)
    return _finalize_capture!(det.params.sensor, det, rng, exposure_duration)
end

function _finalize_incremental_capture!(::AbstractFrameSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    return finalize_readout_pipeline!(det, rng, exposure_duration,
        zero(exposure_duration))
end

function finalize_incremental_capture!(det::Detector, rng::AbstractRNG,
    exposure_duration::Real)
    return _finalize_incremental_capture!(det.params.sensor, det, rng,
        exposure_duration)
end

validated_temporal_frame(frame::AbstractMatrix) = frame
validated_temporal_frame(frame) =
    throw(InvalidConfiguration("FunctionFrameSource must return an AbstractMatrix"))

function initial_temporal_frame(source::FunctionFrameSource, det::Detector,
    sample_offset_s, exposure_duration)
    return validated_temporal_frame(source.f(sample_offset_s))
end

function initial_temporal_frame(source::InPlaceFrameSource, det::Detector,
    sample_offset_s, exposure_duration)
    frame = ensure_temporal_buffer!(det, source.frame_size)
    sample_frame!(frame, source, sample_offset_s)
    return frame
end

function initial_temporal_frame(source::FunctionExposureFrameSource,
    det::Detector, start_offset_s, exposure_duration)
    return validated_temporal_frame(source.f(start_offset_s,
        exposure_duration))
end

function initial_temporal_frame(source::InPlaceExposureFrameSource,
    det::Detector, start_offset_s, exposure_duration)
    frame = ensure_temporal_buffer!(det, source.frame_size)
    sample_exposure_frame!(frame, source, start_offset_s,
        exposure_duration)
    return frame
end

function ensure_temporal_buffer!(det::Detector, dims::Tuple{Int,Int})
    if size(det.workspace.temporal_buffer) != dims
        det.workspace.temporal_buffer = similar(det.workspace.temporal_buffer, dims...)
    end
    return det.workspace.temporal_buffer
end

function capture_temporal_signal!(det::Detector, source::AbstractTemporalFrameSource, first_frame::AbstractMatrix,
    rng::AbstractRNG, exposure_duration::Real, ::GlobalShutter)
    capture_signal_pipeline!(det, first_frame, rng, exposure_duration)
    return det.products.frame
end

rolling_exposure_start_offset_s(::RollingExposure, line_index, line_duration, exposure_duration, ::Type{T}) where {T<:AbstractFloat} =
    T(line_index) * T(line_duration)
rolling_exposure_duration(::RollingExposure, line_index, line_duration, exposure_duration, ::Type{T}) where {T<:AbstractFloat} =
    T(exposure_duration)
rolling_exposure_start_offset_s(::GlobalResetExposure, line_index, line_duration, exposure_duration, ::Type{T}) where {T<:AbstractFloat} =
    zero(T)
rolling_exposure_duration(::GlobalResetExposure, line_index, line_duration, exposure_duration, ::Type{T}) where {T<:AbstractFloat} =
    T(exposure_duration) + T(line_index) * T(line_duration)

function capture_temporal_signal!(det::Detector, source::AbstractTemporalFrameSource, first_frame::AbstractMatrix,
    rng::AbstractRNG, exposure_duration::Real, timing::RollingShutter)
    fill_frame!(det, first_frame, exposure_duration)
    det.state.accum_buffer .= det.products.frame

    scratch = ensure_temporal_buffer!(det, size(first_frame))
    n_rows = size(det.products.frame, 1)
    group_size = timing.row_group_size
    value_type = eltype(det.products.frame)
    for row_lo in (firstindex(det.products.frame, 1) + group_size):group_size:n_rows
        row_hi = min(row_lo + group_size - 1, n_rows)
        line_index = div(row_lo - 1, group_size)
        sample_offset_s = rolling_exposure_start_offset_s(timing.exposure_mode, line_index, timing.line_duration, exposure_duration, value_type)
        group_exposure = rolling_exposure_duration(timing.exposure_mode, line_index, timing.line_duration, exposure_duration, value_type)
        sample_exposure_frame!(scratch, source, sample_offset_s, group_exposure)
        fill_frame!(det, scratch, group_exposure)
        @views det.state.accum_buffer[row_lo:row_hi, :] .= det.products.frame[row_lo:row_hi, :]
    end

    det.products.frame .= det.state.accum_buffer
    apply_signal_defects!(det.params.defect_model, det, exposure_duration)
    apply_sensor_persistence!(det.params.sensor, det, exposure_duration)
    photon_noise_enabled(det) && poisson_noise_frame!(det, rng, det.products.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_duration)
    return det.products.frame
end

function _write_output!(::DetectorDirectStrategy, det::Detector, output::AbstractMatrix,
    source::AbstractMatrix)
    if eltype(output) <: Integer
        write_integer_output!(output, source)
    else
        copyto!(output, source)
    end
    return output
end

function _write_output!(::DetectorHostMirrorStrategy, det::Detector,
    output::AbstractMatrix, source::AbstractMatrix)
    frame_host = detector_host_frame!(det, det.products.frame)
    window = det.params.readout_window
    source_host = window === nothing ? frame_host :
        @view(frame_host[window.rows, window.cols])
    output_host = det.workspace.output_buffer_host
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
    output = det.products.output_buffer
    if output === nothing
        window === nothing && return det.products.frame
        throw(InvalidConfiguration("Detector readout_window requires an allocated output buffer"))
    end
    source = window === nothing ? det.products.frame : @view(det.products.frame[window.rows, window.cols])
    strategy = detector_execution_strategy(typeof(execution_style(output)), typeof(det))
    return _write_output!(strategy, det, output, source)
end

@inline function require_whole_capture_idle(det::Detector)
    iszero(det.state.integrated_time) && det.state.readout_ready ||
        throw(InvalidConfiguration(
            "cannot start a whole detector exposure while an incremental " *
            "exposure is pending; complete it with integration_duration or call " *
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
    exposure_duration = det.params.exposure_duration
    capture_signal!(det, photon_rate, rng, exposure_duration,
        quantum_efficiency)
    finalize_capture!(det, rng, exposure_duration)
    advance_thermal!(det, exposure_duration)
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
        effective_qe(det, src, eltype(det.products.frame)), rng)
end

function capture!(det::Detector, psf::AbstractMatrix{T}, src::AbstractSource;
    rng::AbstractRNG=runtime_rng()) where {T}
    return capture!(det, psf, src, rng)
end

function capture!(det::Detector, source::AbstractTemporalFrameSource, rng::AbstractRNG)
    require_whole_capture_idle(det)
    exposure_duration = det.params.exposure_duration
    first_frame = initial_temporal_frame(source, det,
        zero(eltype(det.products.frame)), exposure_duration)
    capture_temporal_signal!(det, source, first_frame, rng, exposure_duration,
        det.params.timing_model)
    finalize_capture!(det, rng, exposure_duration)
    advance_thermal!(det, exposure_duration)
    return completed_whole_capture_output!(det)
end

function capture!(det::Detector, source::AbstractTemporalFrameSource;
    rng::AbstractRNG=runtime_rng())
    return capture!(det, source, rng)
end

"""
    capture_incremental!(detector, photon_rate, rng, integration_duration,
        quantum_efficiency=detector.params.qe)

Accumulate one positive `integration_duration` in seconds from a cell-integrated
photon-arrival-rate matrix. This frame-step convenience finalizes automatically
when the configured exposure duration is reached. `integration_duration`
is neither an absolute timestamp nor the period between consecutive samples;
scheduled detector events own their timestamps and completion semantics.
"""
function capture_incremental!(det::Detector, photon_rate::AbstractMatrix,
    rng::AbstractRNG, integration_duration::Real, qe=det.params.qe)
    if !iszero(det.state.integrated_time) || !det.state.readout_ready
        size(photon_rate) == size(det.workspace.presampling_buffer) ||
            throw(DimensionMismatchError(
                "incremental detector input dimensions cannot change while " *
                "an exposure is pending"))
    end
    prepare_detector_buffers!(det, size(photon_rate))
    T = eltype(det.products.frame)
    dt = T(integration_duration)
    isfinite(dt) && dt > zero(T) || throw(InvalidConfiguration(
        "integration_duration must be finite and > 0"))
    remaining = det.params.exposure_duration - det.state.integrated_time
    tolerance = T(8) * eps(det.params.exposure_duration) *
        max(one(T), abs(det.params.exposure_duration))
    dt <= remaining + tolerance || throw(InvalidConfiguration(
        "integration_duration exceeds the remaining detector exposure " *
        "duration"))
    dt = min(dt, remaining)

    exposure_start = iszero(det.state.integrated_time)
    exposure_start && fill!(det.state.accum_buffer,
        zero(eltype(det.state.accum_buffer)))
    capture_signal_pipeline!(det, photon_rate, rng, dt, qe, exposure_start,
        det.params.exposure_duration)
    accumulate_incremental_charge_generation!(det, rng, dt)
    det.state.accum_buffer .+= det.products.frame
    det.state.integrated_time += dt
    advance_thermal!(det, dt)
    det.state.readout_ready = false
    if det.state.integrated_time + tolerance >= det.params.exposure_duration
        det.products.frame .= det.state.accum_buffer
        finalize_incremental_capture!(det, rng, det.state.integrated_time)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    return write_output!(det)
end

function capture!(det::Detector, photon_rate::AbstractMatrix{T};
    rng::AbstractRNG=runtime_rng(),
    integration_duration::Union{Nothing,Real}=nothing) where {T}
    if integration_duration === nothing
        return capture!(det, photon_rate, rng)
    end
    return capture_incremental!(det, photon_rate, rng, integration_duration)
end
