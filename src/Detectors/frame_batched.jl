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
    is_null_readout_correction(det.params.correction_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires null readout correction"))
    supports_separable_application(det.params.response_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires a maintained separable frame response"))
    is_global_shutter(det.params.timing_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires global-shutter timing semantics"))
    is_null_persistence(persistence_model(det.params.sensor)) ||
        throw(InvalidConfiguration("batched detector capture currently requires null detector persistence"))
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
    dark_signal = effective_dark_current(det) * effective_dark_current_time(det.params.sensor, exposure_time)
    if dark_signal <= 0
        return cube
    end
    fill!(scratch, dark_signal)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_signal_defects!(::NullDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube
_batched_dark_defects!(::NullDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_signal_defects!(model::PixelResponseNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.gain_map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("PixelResponseNonuniformity gain_map size must match detector frame size"))
    cube .*= reshape(model.gain_map, size(model.gain_map)..., 1)
    return cube
end

_batched_dark_defects!(::PixelResponseNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube
_batched_signal_defects!(::DarkSignalNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_dark_defects!(model::DarkSignalNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.dark_map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("DarkSignalNonuniformity dark_map size must match detector frame size"))
    cube .+= reshape(model.dark_map .* exposure_time, size(model.dark_map)..., 1)
    return cube
end

function _batched_signal_defects!(model::BadPixelMask, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.mask) == size(cube)[1:2] ||
        throw(DimensionMismatchError("BadPixelMask mask size must match detector frame size"))
    throughput = model.throughput
    throughput == one(throughput) && return cube
    cube .= ifelse.(reshape(model.mask, size(model.mask)..., 1), throughput .* cube, cube)
    return cube
end

_batched_dark_defects!(::BadPixelMask, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_signal_defects!(model::CompositeDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    foreach(stage -> _batched_signal_defects!(stage, cube, scratch, exposure_time), model.stages)
    return cube
end

function _batched_dark_defects!(model::CompositeDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    foreach(stage -> _batched_dark_defects!(stage, cube, scratch, exposure_time), model.stages)
    return cube
end

_batched_frame_nonlinearity!(::NullFrameNonlinearity, cube::AbstractArray) = cube

function _batched_frame_nonlinearity!(model::SaturatingFrameNonlinearity, cube::AbstractArray)
    coeff = model.coefficient
    coeff <= zero(coeff) && return cube
    @. cube = cube / (1 + coeff * cube)
    return cube
end

_batched_pre_readout_gain!(::FrameSensorType, det::Detector, cube::AbstractArray, rng::AbstractRNG) = cube
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

function _batched_apply_response!(::ExecutionStyle, ::NullFrameResponse, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    return cube
end

function _batched_apply_response!(::ScalarCPUStyle, ::NullFrameResponse, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    return cube
end

function _batched_apply_response!(::ScalarCPUStyle, model::AbstractFrameResponse, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    for b in axes(cube, 1)
        apply_response!(ScalarCPUStyle(), model, @view(cube[b, :, :]), @view(scratch[b, :, :]))
    end
    return cube
end

function _batched_apply_response!(::ScalarCPUStyle, model::GaussianPixelResponse,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(ScalarCPUStyle(), cube, scratch, model.kernel, model.kernel)
    return cube
end

function _batched_apply_response!(::ScalarCPUStyle, model::RectangularPixelAperture,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(ScalarCPUStyle(), cube, scratch, model.kernel_y, model.kernel_x)
    return cube
end

function _batched_apply_response!(::ScalarCPUStyle, model::SeparablePixelMTF,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(ScalarCPUStyle(), cube, scratch, model.kernel_y, model.kernel_x)
    return cube
end

function _batched_apply_response!(style::AcceleratorStyle, model::GaussianPixelResponse,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(style, cube, scratch, model.kernel, model.kernel)
    return cube
end

function _batched_apply_response!(style::AcceleratorStyle, model::RectangularPixelAperture,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(style, cube, scratch, model.kernel_y, model.kernel_x)
    return cube
end

function _batched_apply_response!(style::AcceleratorStyle, model::SeparablePixelMTF,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    _batched_apply_separable_response!(style, cube, scratch, model.kernel_y, model.kernel_x)
    return cube
end

function _batched_apply_separable_response!(style::AcceleratorStyle, cube::AbstractArray{T,3},
    scratch::AbstractArray{T,3}, kernel_y::AbstractVector, kernel_x::AbstractVector) where {T}
    n_batch, n, m = size(cube)
    radius_x = fld(length(kernel_x), 2)
    radius_y = fld(length(kernel_y), 2)
    launch_kernel_async!(style, separable_response_stack_rows_kernel!, scratch, cube,
        kernel_x, radius_x, n_batch, n, m, length(kernel_x); ndrange=size(cube))
    launch_kernel_async!(style, separable_response_stack_cols_kernel!, cube, scratch,
        kernel_y, radius_y, n_batch, n, m, length(kernel_y); ndrange=size(cube))
    return cube
end

function _batched_apply_separable_response!(::ScalarCPUStyle, cube::AbstractArray{T,3},
    scratch::AbstractArray{T,3}, kernel_y::AbstractVector, kernel_x::AbstractVector) where {T}
    n_batch, n, m = size(cube)
    radius_x = fld(length(kernel_x), 2)
    radius_y = fld(length(kernel_y), 2)
    @inbounds for j in 1:m, i in 1:n, b in 1:n_batch
        acc = zero(T)
        for kk in eachindex(kernel_x)
            jj = clamp(j + kk - radius_x - 1, 1, m)
            acc += kernel_x[kk] * cube[b, i, jj]
        end
        scratch[b, i, j] = acc
    end
    @inbounds for j in 1:m, i in 1:n, b in 1:n_batch
        acc = zero(T)
        for kk in eachindex(kernel_y)
            ii = clamp(i + kk - radius_y - 1, 1, n)
            acc += kernel_y[kk] * scratch[b, ii, j]
        end
        cube[b, i, j] = acc
    end
    return cube
end

function capture_stack!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    _require_batched_detector_compat(det, cube, scratch)
    exposure_time = det.params.integration_time
    cube .*= det.params.qe * exposure_time
    _batched_signal_defects!(det.params.defect_model, cube, scratch, exposure_time)
    _batched_apply_response!(execution_style(cube), det.params.response_model, cube, scratch)
    capture_stack_poisson_noise!(det, cube, rng)
    _batched_background_flux!(det.background_flux, det, cube, scratch, rng, exposure_time)
    _batched_dark_current!(det, cube, scratch, rng, exposure_time)
    _batched_dark_defects!(det.params.defect_model, cube, scratch, exposure_time)
    _batched_sensor_statistics!(det.params.sensor, det, cube, scratch, rng)
    _batched_frame_nonlinearity!(det.params.nonlinearity_model, cube)
    apply_saturation!(det, cube)
    _batched_pre_readout_gain!(det.params.sensor, det, cube, rng)
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
