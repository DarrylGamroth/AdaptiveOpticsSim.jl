function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    poisson_noise!(rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    poisson_noise!(rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

apply_background_flux!(::NoBackground, det::Detector, rng::AbstractRNG, exposure_time::Real) = det.state.frame

function apply_background_flux!(background::ScalarBackground, det::Detector, rng::AbstractRNG, exposure_time::Real)
    fill!(det.state.noise_buffer, background.level * exposure_time)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_background_flux!(background::BackgroundFrame, det::Detector, rng::AbstractRNG, exposure_time::Real)
    if size(background.map) != size(det.state.frame)
        throw(DimensionMismatchError("background_flux size must match detector frame size"))
    end
    copyto!(det.state.noise_buffer, background.map)
    det.state.noise_buffer .*= exposure_time
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_dark_current!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    dark_signal = det.params.dark_current * effective_dark_current_time(det.params.sensor, exposure_time)
    if dark_signal <= 0
        return det.state.frame
    end
    fill!(det.state.noise_buffer, dark_signal)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

sensor_saturation_limit(det::Detector) = sensor_saturation_limit(det.params.sensor, det)
sensor_saturation_limit(::FrameSensorType, det::Detector) = det.params.full_well

function apply_saturation!(det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    clamp!(det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

apply_sensor_statistics!(sensor::FrameSensorType, det::Detector, rng::AbstractRNG) = det.state.frame

function apply_avalanche_excess_noise!(factor, det::Detector, rng::AbstractRNG)
    factor <= one(factor) && return det.state.frame
    randn_backend!(rng, det.state.noise_buffer)
    scale2 = factor * factor - one(factor)
    zero_t = zero(eltype(det.state.frame))
    @. det.state.frame += sqrt(max(scale2 * det.state.frame, zero_t)) * det.state.noise_buffer
    return det.state.frame
end

apply_pre_readout_gain!(::FrameSensorType, det::Detector) = det.state.frame
apply_post_readout_gain!(::FrameSensorType, det::Detector) = det.state.frame
reset_readout_products!(det::Detector) = (det.state.readout_products = NoFrameReadoutProducts(); det)

apply_readout_noise!(det::Detector{NoiseNone}, rng::AbstractRNG) = det.state.frame
apply_readout_noise!(det::Detector{NoisePhoton}, rng::AbstractRNG) = det.state.frame

function apply_readout_noise!(det::Detector{<:NoiseReadout}, rng::AbstractRNG)
    randn_backend!(rng, det.state.noise_buffer)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    det.state.frame .+= sigma .* det.state.noise_buffer
    return det.state.frame
end

function apply_readout_noise!(det::Detector{<:NoisePhotonReadout}, rng::AbstractRNG)
    randn_backend!(rng, det.state.noise_buffer)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    det.state.frame .+= sigma .* det.state.noise_buffer
    return det.state.frame
end

function apply_quantization!(det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    levels = exp2(eltype(det.state.frame)(bits))
    full_well = det.params.full_well
    if full_well === nothing
        peak = maximum(det.state.frame)
        if peak > 0
            det.state.frame .*= levels / peak
        end
    else
        det.state.frame .*= (levels - one(levels)) / full_well
        clamp!(det.state.frame, zero(eltype(det.state.frame)), levels - one(levels))
    end
    return det.state.frame
end

subtract_background_map!(::NoBackground, det::Detector) = det.state.frame

function subtract_background_map!(background::ScalarBackground, det::Detector)
    det.state.frame .-= background.level
    return det.state.frame
end

function subtract_background_map!(background::BackgroundFrame, det::Detector)
    if size(background.map) != size(det.state.frame)
        throw(DimensionMismatchError("background_map size must match detector frame size"))
    end
    det.state.frame .-= background.map
    return det.state.frame
end

readout_product_shape(det::Detector) = det.params.readout_window === nothing ?
    size(det.state.frame) :
    (length(det.params.readout_window.rows), length(det.params.readout_window.cols))

function _copy_windowed_frame(frame::AbstractMatrix, det::Detector)
    window = det.params.readout_window
    if window === nothing
        return copy(frame)
    end
    return copy(@view(frame[window.rows, window.cols]))
end

function _copy_windowed_cube(cube::AbstractArray{T,3}, det::Detector) where {T}
    window = det.params.readout_window
    if window === nothing
        return copy(cube)
    end
    return copy(@view(cube[window.rows, window.cols, :]))
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

function apply_readout_correction!(model::ReferencePixelCommonModeCorrection, frame::AbstractMatrix)
    frame .-= _reference_pixel_bias(model, frame)
    return frame
end

function apply_readout_correction!(model::FrameReadoutCorrectionModel, cube::AbstractArray{T,3}) where {T}
    for k in axes(cube, 3)
        apply_readout_correction!(model, @view(cube[:, :, k]))
    end
    return cube
end

finalize_readout_products!(::FrameSensorType, det::Detector, rng::AbstractRNG, exposure_time::Real) =
    reset_readout_products!(det)

function finalize_capture!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    apply_dark_current!(det, rng, exposure_time)
    apply_saturation!(det)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_pre_readout_gain!(det.params.sensor, det)
    apply_readout_noise!(det, rng)
    apply_post_readout_gain!(det.params.sensor, det)
    finalize_readout_products!(det.params.sensor, det, rng, exposure_time)
    apply_readout_correction!(det.params.correction_model, det.state.frame)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    return det.state.frame
end

function write_output!(det::Detector)
    window = det.params.readout_window
    output = det.state.output_buffer
    if output === nothing
        window === nothing && return det.state.frame
        throw(InvalidConfiguration("Detector readout_window requires an allocated output buffer"))
    end
    output_eltype = eltype(output)
    source = window === nothing ? det.state.frame : @view(det.state.frame[window.rows, window.cols])
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(source), lo, hi))
    else
        output .= source
    end
    return output
end

function capture!(det::Detector, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    capture_signal!(det, psf, rng, det.params.integration_time)
    finalize_capture!(det, rng, det.params.integration_time)
    return write_output!(det)
end

function capture!(det::Detector, psf::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng(),
    sample_time::Union{Nothing,Real}=nothing) where {T}
    if sample_time === nothing
        return capture!(det, psf, rng)
    end
    dt = eltype(det.state.frame)(sample_time)
    if dt <= 0
        throw(InvalidConfiguration("sample_time must be > 0"))
    end
    if det.params.integration_time < dt
        throw(InvalidConfiguration("sample_time must be <= detector integration_time"))
    end
    capture_signal!(det, psf, rng, dt)
    if size(det.state.accum_buffer) != size(det.state.frame)
        det.state.accum_buffer = similar(det.state.frame, size(det.state.frame)...)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    end
    det.state.accum_buffer .+= det.state.frame
    det.state.integrated_time += dt
    det.state.readout_ready = false
    if det.state.integrated_time + eps(dt) >= det.params.integration_time
        det.state.frame .= det.state.accum_buffer
        finalize_capture!(det, rng, det.state.integrated_time)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    return write_output!(det)
end
