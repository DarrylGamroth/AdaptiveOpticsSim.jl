apply_frame_response!(::NullFrameResponse, det::Detector) = det.state.frame

function apply_frame_response!(model::SeparableGaussianPixelResponse, det::Detector)
    return apply_frame_response!(execution_style(det.state.frame), model, det)
end

function apply_frame_response!(::ScalarCPUStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * frame[i, jj]
        end
        scratch[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * scratch[ii, j]
        end
        frame[i, j] = acc
    end
    return frame
end

function apply_frame_response!(style::AcceleratorStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    launch_kernel!(style, separable_response_rows_kernel!, scratch, frame, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    launch_kernel!(style, separable_response_cols_kernel!, frame, scratch, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    return frame
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}, exposure_time::Real) where {T}
    n_in, m_in = size(psf)
    sampling = det.params.psf_sampling
    binning = det.params.binning
    if sampling < 1 || binning < 1
        throw(InvalidConfiguration("psf_sampling and binning must be >= 1"))
    end
    if n_in % sampling != 0 || m_in % sampling != 0
        throw(DimensionMismatchError("psf_sampling must evenly divide input dimensions"))
    end
    n_mid = div(n_in, sampling)
    m_mid = div(m_in, sampling)
    if n_mid % binning != 0 || m_mid % binning != 0
        throw(DimensionMismatchError("binning must evenly divide sampled dimensions"))
    end
    n_out = div(n_mid, binning)
    m_out = div(m_mid, binning)
    ensure_buffers!(det, n_mid, m_mid, n_out, m_out)

    if sampling > 1
        if binning > 1
            bin2d!(det.state.frame, psf, sampling * binning)
        else
            bin2d!(det.state.bin_buffer, psf, sampling)
            det.state.frame .= det.state.bin_buffer
        end
    else
        if binning > 1
            bin2d!(det.state.frame, psf, binning)
        else
            det.state.frame .= psf
        end
    end
    @. det.state.frame *= det.params.qe * exposure_time
    apply_frame_response!(det.params.response_model, det)
    return det.state.frame
end

fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T} = fill_frame!(det, psf, det.params.integration_time)

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

function finalize_capture!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    apply_dark_current!(det, rng, exposure_time)
    apply_saturation!(det)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_pre_readout_gain!(det.params.sensor, det)
    apply_readout_noise!(det, rng)
    apply_post_readout_gain!(det.params.sensor, det)
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

function _require_batched_detector_compat(det::Detector, cube::AbstractArray, scratch::AbstractArray)
    size(cube) == size(scratch) ||
        throw(DimensionMismatchError("batched detector scratch must match cube size"))
    ndims(cube) == 3 || throw(DimensionMismatchError("batched detector input must be 3D"))
    det.params.psf_sampling == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires psf_sampling == 1"))
    det.params.binning == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires binning == 1"))
    det.params.output_precision === nothing ||
        throw(InvalidConfiguration("batched detector capture currently requires output_precision === nothing"))
    det.params.readout_window === nothing ||
        throw(InvalidConfiguration("batched detector capture currently requires full-frame readout"))
    supports_detector_mtf(det) &&
        throw(InvalidConfiguration("batched detector capture currently requires NullFrameResponse()"))
    _require_batched_sensor_compat(det.params.sensor)
    return nothing
end

_require_batched_sensor_compat(::FrameSensorType) = nothing

_batched_background_map!(::NoBackground, cube::AbstractArray, scratch::AbstractArray) = cube
_batched_background_map!(background::ScalarBackground, cube::AbstractArray, scratch::AbstractArray) = (cube .-= background.level; cube)

function _batched_background_map!(background::BackgroundFrame, cube::AbstractArray, scratch::AbstractArray)
    size(background.map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("background_map size must match detector frame size"))
    cube .-= reshape(background.map, size(background.map)..., 1)
    return cube
end

_batched_background_flux!(::NoBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real) = cube

function _batched_background_flux!(background::ScalarBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    fill!(scratch, background.level * exposure_time)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

function _batched_background_flux!(background::BackgroundFrame, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    size(background.map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("background_flux size must match detector frame size"))
    scratch .= reshape(background.map, size(background.map)..., 1)
    scratch .*= exposure_time
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

function _batched_dark_current!(det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real)
    dark_signal = det.params.dark_current * effective_dark_current_time(det.params.sensor, exposure_time)
    if dark_signal <= 0
        return cube
    end
    fill!(scratch, dark_signal)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_pre_readout_gain!(::FrameSensorType, det::Detector, cube::AbstractArray) = cube
_batched_sensor_statistics!(sensor::FrameSensorType, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube

function _batched_avalanche_excess_noise!(factor, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    factor <= one(factor) && return cube
    randn_backend!(rng, scratch)
    scale2 = factor * factor - one(factor)
    zero_t = zero(eltype(cube))
    @. cube += sqrt(max(scale2 * cube, zero_t)) * scratch
    return cube
end

_batched_post_readout_gain!(::FrameSensorType, det::Detector, cube::AbstractArray) = cube
_batched_readout_noise!(det::Detector{NoiseNone}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube
_batched_readout_noise!(det::Detector{NoisePhoton}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube

function _batched_readout_noise!(det::Detector{<:NoiseReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    randn_backend!(rng, scratch)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    cube .+= sigma .* scratch
    return cube
end

function _batched_readout_noise!(det::Detector{<:NoisePhotonReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    randn_backend!(rng, scratch)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    cube .+= sigma .* scratch
    return cube
end

function _batched_quantization!(det::Detector, cube::AbstractArray)
    bits = det.params.bits
    bits === nothing && return cube
    levels = exp2(eltype(cube)(bits))
    full_well = det.params.full_well
    if full_well === nothing
        peak = maximum(cube)
        if peak > 0
            cube .*= levels / peak
        end
    else
        cube .*= (levels - one(levels)) / full_well
        clamp!(cube, zero(eltype(cube)), levels - one(levels))
    end
    return cube
end

capture_stack_poisson_noise!(det::Detector, cube::AbstractArray, rng::AbstractRNG) = cube

function capture_stack_poisson_noise!(det::Detector{NoisePhoton}, cube::AbstractArray, rng::AbstractRNG)
    poisson_noise!(rng, cube)
    return cube
end

function capture_stack_poisson_noise!(det::Detector{<:NoisePhotonReadout}, cube::AbstractArray, rng::AbstractRNG)
    poisson_noise!(rng, cube)
    return cube
end

function capture_stack!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    _require_batched_detector_compat(det, cube, scratch)
    exposure_time = det.params.integration_time
    cube .*= det.params.qe * exposure_time
    capture_stack_poisson_noise!(det, cube, rng)
    _batched_background_flux!(det.background_flux, det, cube, scratch, rng, exposure_time)
    _batched_dark_current!(det, cube, scratch, rng, exposure_time)
    apply_saturation!(det, cube)
    _batched_sensor_statistics!(det.params.sensor, det, cube, scratch, rng)
    _batched_pre_readout_gain!(det.params.sensor, det, cube)
    _batched_readout_noise!(det, cube, scratch, rng)
    _batched_post_readout_gain!(det.params.sensor, det, cube)
    _batched_quantization!(det, cube)
    _batched_background_map!(det.background_map, cube, scratch)
    return cube
end

function apply_saturation!(det::Detector, cube::AbstractArray)
    full_well = det.params.full_well
    full_well === nothing && return cube
    clamp!(cube, zero(eltype(cube)), full_well)
    return cube
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

function ensure_buffers!(det::Detector, n_mid::Int, m_mid::Int, n_out::Int, m_out::Int)
    if size(det.state.frame) != (n_out, m_out)
        det.state.frame = similar(det.state.frame, n_out, m_out)
    end
    if size(det.state.response_buffer) != (n_out, m_out)
        det.state.response_buffer = similar(det.state.response_buffer, n_out, m_out)
    end
    if size(det.state.bin_buffer) != (n_mid, m_mid)
        det.state.bin_buffer = similar(det.state.bin_buffer, n_mid, m_mid)
    end
    if size(det.state.noise_buffer) != (n_out, m_out)
        det.state.noise_buffer = similar(det.state.noise_buffer, n_out, m_out)
    end
    if size(det.state.accum_buffer) != (n_out, m_out)
        det.state.accum_buffer = similar(det.state.accum_buffer, n_out, m_out)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    window = det.params.readout_window === nothing ? nothing : validate_readout_window(det.params.readout_window, n_out, m_out)
    out_rows = window === nothing ? n_out : length(window.rows)
    out_cols = window === nothing ? m_out : length(window.cols)
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != (out_rows, out_cols)
        det.state.output_buffer = similar(det.state.output_buffer, out_rows, out_cols)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end
