function _require_batched_detector_compat(det::Detector, cube::AbstractArray, scratch::AbstractArray)
    style = execution_style(cube)
    size(cube) == size(scratch) ||
        throw(DimensionMismatchError("batched detector scratch must match cube size"))
    ndims(cube) == 3 || throw(DimensionMismatchError("batched detector input must be 3D"))
    det.params.psf_sampling == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires psf_sampling == 1"))
    det.params.binning == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires binning == 1"))
    det.params.readout_window === nothing ||
        throw(InvalidConfiguration("batched detector capture currently requires full-frame readout"))
    supports_batched_readout_correction(det.params.correction_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires a maintained batched readout correction"))
    supports_batched_response_application(style, det.params.response_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires a maintained batched frame response"))
    is_global_shutter(det.params.timing_model) ||
        throw(InvalidConfiguration("batched detector capture currently requires global-shutter timing semantics"))
    is_null_persistence(persistence_model(det.params.sensor)) ||
        throw(InvalidConfiguration("batched detector capture currently requires null detector persistence"))
    _require_batched_sensor_compat(det.params.sensor)
    return nothing
end

_require_batched_sensor_compat(::FrameSensorType) = nothing

_batched_frame_shape(cube::AbstractArray{T,3}) where {T} = (size(cube, 2), size(cube, 3))
_batched_frame_map(map::AbstractMatrix) = reshape(map, 1, size(map, 1), size(map, 2))

_batched_poisson_noise_async!(::ScalarCPUStyle, det::Detector, cube::AbstractArray, rng::AbstractRNG) =
    poisson_noise_frame!(det, rng, cube)
_batched_poisson_noise_async!(style::AcceleratorStyle, det::Detector, cube::AbstractArray{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    poisson_noise_async!(style, rng, cube)

_batched_randn_noise_async!(::ScalarCPUStyle, det::Detector, cube::AbstractArray, rng::AbstractRNG) =
    randn_frame_noise!(det, rng, cube)
_batched_randn_noise_async!(style::AcceleratorStyle, det::Detector, cube::AbstractArray{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    randn_backend_async!(style, rng, cube)

_batched_background_map!(::NoBackground, cube::AbstractArray, scratch::AbstractArray) = cube
_batched_background_map!(background::ScalarBackground, cube::AbstractArray, scratch::AbstractArray) = (cube .-= background.level; cube)

function _batched_background_map!(background::BackgroundFrame, cube::AbstractArray, scratch::AbstractArray)
    size(background.map) == _batched_frame_shape(cube) ||
        throw(DimensionMismatchError("background_map size must match detector frame size"))
    cube .-= _batched_frame_map(background.map)
    return cube
end

_batched_background_flux!(::NoBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real) = cube

function _batched_background_flux!(background::ScalarBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    fill!(scratch, background.level * exposure_time)
    _batched_poisson_noise_async!(execution_style(cube), det, scratch, rng)
    cube .+= scratch
    return cube
end

function _batched_background_flux!(background::BackgroundFrame, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    size(background.map) == _batched_frame_shape(cube) ||
        throw(DimensionMismatchError("background_flux size must match detector frame size"))
    scratch .= _batched_frame_map(background.map)
    scratch .*= exposure_time
    _batched_poisson_noise_async!(execution_style(cube), det, scratch, rng)
    cube .+= scratch
    return cube
end

function _batched_dark_current!(det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real)
    dark_signal = effective_dark_current(det) * effective_dark_current_time(det.params.sensor, exposure_time)
    if dark_signal <= 0
        return cube
    end
    fill!(scratch, dark_signal)
    _batched_poisson_noise_async!(execution_style(cube), det, scratch, rng)
    cube .+= scratch
    return cube
end

_batched_signal_defects!(::NullDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube
_batched_dark_defects!(::NullDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_signal_defects!(model::PixelResponseNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.gain_map) == _batched_frame_shape(cube) ||
        throw(DimensionMismatchError("PixelResponseNonuniformity gain_map size must match detector frame size"))
    cube .*= _batched_frame_map(model.gain_map)
    return cube
end

_batched_dark_defects!(::PixelResponseNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube
_batched_signal_defects!(::DarkSignalNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_dark_defects!(model::DarkSignalNonuniformity, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.dark_map) == _batched_frame_shape(cube) ||
        throw(DimensionMismatchError("DarkSignalNonuniformity dark_map size must match detector frame size"))
    cube .+= _batched_frame_map(model.dark_map .* exposure_time)
    return cube
end

function _batched_signal_defects!(model::BadPixelMask, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    size(model.mask) == _batched_frame_shape(cube) ||
        throw(DimensionMismatchError("BadPixelMask mask size must match detector frame size"))
    throughput = model.throughput
    throughput == one(throughput) && return cube
    cube .= ifelse.(_batched_frame_map(model.mask), throughput .* cube, cube)
    return cube
end

_batched_dark_defects!(::BadPixelMask, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real) = cube

function _batched_signal_defects!(model::CompositeDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    for stage in model.stages
        _batched_signal_defects!(stage, cube, scratch, exposure_time)
    end
    return cube
end

function _batched_dark_defects!(model::CompositeDetectorDefectModel, cube::AbstractArray, scratch::AbstractArray, exposure_time::Real)
    for stage in model.stages
        _batched_dark_defects!(stage, cube, scratch, exposure_time)
    end
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
_batched_apply_readout_correction!(::NullFrameReadoutCorrection, cube::AbstractArray{T,3}) where {T} = cube
_batched_apply_readout_correction!(::NullFrameReadoutCorrection, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T} = cube

function _batched_apply_readout_correction!(model::FrameReadoutCorrectionModel, cube::AbstractArray{T,3}) where {T}
    for b in axes(cube, 1)
        apply_readout_correction!(model, @view(cube[b, :, :]))
    end
    return cube
end

function _batched_apply_readout_correction!(model::FrameReadoutCorrectionModel, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    return _batched_apply_readout_correction!(execution_style(cube), model, cube, scratch)
end

function _batched_apply_readout_correction!(::ExecutionStyle, model::FrameReadoutCorrectionModel,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    return _batched_apply_readout_correction!(model, cube)
end

@kernel function reference_pixel_bias_kernel!(bias, cube, edge_rows::Int, edge_cols::Int, n::Int, m::Int)
    idx = @index(Global, Linear)
    if idx <= size(cube, 1)
        T = eltype(cube)
        n_edge_rows = min(edge_rows, n)
        n_edge_cols = min(edge_cols, m)
        total = zero(T)
        count = 0
        @inbounds begin
            if n_edge_rows > 0
                for i in 1:n_edge_rows, j in 1:m
                    total += cube[idx, i, j]
                end
                count += n_edge_rows * m
                if n > n_edge_rows
                    for i in (n - n_edge_rows + 1):n, j in 1:m
                        total += cube[idx, i, j]
                    end
                    count += n_edge_rows * m
                end
            end
            row_lo = n_edge_rows + 1
            row_hi = n - n_edge_rows
            if n_edge_cols > 0 && row_lo <= row_hi
                for i in row_lo:row_hi, j in 1:n_edge_cols
                    total += cube[idx, i, j]
                end
                count += (row_hi - row_lo + 1) * n_edge_cols
                if m > n_edge_cols
                    for i in row_lo:row_hi, j in (m - n_edge_cols + 1):m
                        total += cube[idx, i, j]
                    end
                    count += (row_hi - row_lo + 1) * n_edge_cols
                end
            end
        end
        @inbounds bias[idx] = count > 0 ? total / T(count) : zero(T)
    end
end

function _batched_apply_readout_correction!(style::AcceleratorStyle, model::ReferencePixelCommonModeCorrection,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    bias = @view scratch[:, 1, 1]
    launch_kernel!(style, reference_pixel_bias_kernel!, bias, cube, model.edge_rows, model.edge_cols,
        size(cube, 2), size(cube, 3); ndrange=size(cube, 1))
    cube .-= reshape(bias, :, 1, 1)
    return cube
end

@kernel function reference_row_bias_kernel!(scratch, cube, edge_cols::Int, n::Int, m::Int)
    batch_idx, row_idx = @index(Global, NTuple)
    if batch_idx <= size(cube, 1) && row_idx <= n
        T = eltype(cube)
        n_edge = min(edge_cols, m)
        total = zero(T)
        count = 0
        @inbounds begin
            if n_edge > 0
                for j in 1:n_edge
                    total += cube[batch_idx, row_idx, j]
                end
                count += n_edge
                if m > n_edge
                    for j in (m - n_edge + 1):m
                        total += cube[batch_idx, row_idx, j]
                    end
                    count += n_edge
                end
            end
            scratch[batch_idx, row_idx, 1] = count > 0 ? total / T(count) : zero(T)
        end
    end
end

@kernel function subtract_row_bias_kernel!(cube, scratch, n::Int, m::Int)
    batch_idx, row_idx, col_idx = @index(Global, NTuple)
    if batch_idx <= size(cube, 1) && row_idx <= n && col_idx <= m
        @inbounds cube[batch_idx, row_idx, col_idx] -= scratch[batch_idx, row_idx, 1]
    end
end

function _batched_apply_readout_correction!(style::AcceleratorStyle, model::ReferenceRowCommonModeCorrection,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    n = size(cube, 2)
    m = size(cube, 3)
    launch_kernel!(style, reference_row_bias_kernel!, scratch, cube, model.edge_cols, n, m;
        ndrange=(size(cube, 1), n))
    launch_kernel!(style, subtract_row_bias_kernel!, cube, scratch, n, m; ndrange=size(cube))
    return cube
end

@kernel function reference_column_bias_kernel!(scratch, cube, edge_rows::Int, n::Int, m::Int)
    batch_idx, col_idx = @index(Global, NTuple)
    if batch_idx <= size(cube, 1) && col_idx <= m
        T = eltype(cube)
        n_edge = min(edge_rows, n)
        total = zero(T)
        count = 0
        @inbounds begin
            if n_edge > 0
                for i in 1:n_edge
                    total += cube[batch_idx, i, col_idx]
                end
                count += n_edge
                if n > n_edge
                    for i in (n - n_edge + 1):n
                        total += cube[batch_idx, i, col_idx]
                    end
                    count += n_edge
                end
            end
            scratch[batch_idx, 1, col_idx] = count > 0 ? total / T(count) : zero(T)
        end
    end
end

@kernel function subtract_column_bias_kernel!(cube, scratch, n::Int, m::Int)
    batch_idx, row_idx, col_idx = @index(Global, NTuple)
    if batch_idx <= size(cube, 1) && row_idx <= n && col_idx <= m
        @inbounds cube[batch_idx, row_idx, col_idx] -= scratch[batch_idx, 1, col_idx]
    end
end

function _batched_apply_readout_correction!(style::AcceleratorStyle, model::ReferenceColumnCommonModeCorrection,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    n = size(cube, 2)
    m = size(cube, 3)
    launch_kernel!(style, reference_column_bias_kernel!, scratch, cube, model.edge_rows, n, m;
        ndrange=(size(cube, 1), m))
    launch_kernel!(style, subtract_column_bias_kernel!, cube, scratch, n, m; ndrange=size(cube))
    return cube
end

@kernel function reference_output_bias_kernel!(scratch, cube, output_cols::Int, edge_rows::Int, edge_cols::Int, n::Int, m::Int)
    batch_idx, block_idx = @index(Global, NTuple)
    n_blocks = cld(m, output_cols)
    if batch_idx <= size(cube, 1) && block_idx <= n_blocks
        T = eltype(cube)
        col_lo = (block_idx - 1) * output_cols + 1
        col_hi = min(col_lo + output_cols - 1, m)
        block_m = col_hi - col_lo + 1
        n_edge_rows = min(edge_rows, n)
        n_edge_cols = min(edge_cols, block_m)
        total = zero(T)
        count = 0
        @inbounds begin
            if n_edge_rows > 0
                for i in 1:n_edge_rows, j in col_lo:col_hi
                    total += cube[batch_idx, i, j]
                end
                count += n_edge_rows * block_m
                if n > n_edge_rows
                    for i in (n - n_edge_rows + 1):n, j in col_lo:col_hi
                        total += cube[batch_idx, i, j]
                    end
                    count += n_edge_rows * block_m
                end
            end
            row_lo = n_edge_rows + 1
            row_hi = n - n_edge_rows
            if n_edge_cols > 0 && row_lo <= row_hi
                for i in row_lo:row_hi, j in col_lo:(col_lo + n_edge_cols - 1)
                    total += cube[batch_idx, i, j]
                end
                count += (row_hi - row_lo + 1) * n_edge_cols
                if block_m > n_edge_cols
                    for i in row_lo:row_hi, j in (col_hi - n_edge_cols + 1):col_hi
                        total += cube[batch_idx, i, j]
                    end
                    count += (row_hi - row_lo + 1) * n_edge_cols
                end
            end
            scratch[batch_idx, block_idx, 1] = count > 0 ? total / T(count) : zero(T)
        end
    end
end

@kernel function subtract_output_bias_kernel!(cube, scratch, output_cols::Int, n::Int, m::Int)
    batch_idx, row_idx, col_idx = @index(Global, NTuple)
    if batch_idx <= size(cube, 1) && row_idx <= n && col_idx <= m
        block_idx = (col_idx - 1) ÷ output_cols + 1
        @inbounds cube[batch_idx, row_idx, col_idx] -= scratch[batch_idx, block_idx, 1]
    end
end

function _batched_apply_readout_correction!(style::AcceleratorStyle, model::ReferenceOutputCommonModeCorrection,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    n = size(cube, 2)
    m = size(cube, 3)
    n_blocks = cld(m, model.output_cols)
    launch_kernel!(style, reference_output_bias_kernel!, scratch, cube, model.output_cols, model.edge_rows,
        model.edge_cols, n, m; ndrange=(size(cube, 1), n_blocks))
    launch_kernel!(style, subtract_output_bias_kernel!, cube, scratch, model.output_cols, n, m; ndrange=size(cube))
    return cube
end

function _batched_apply_readout_correction!(style::AcceleratorStyle, model::CompositeFrameReadoutCorrection,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    for stage in model.stages
        _batched_apply_readout_correction!(style, stage, cube, scratch)
    end
    return cube
end

function _batched_readout_noise!(det::Detector{<:NoiseReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    _batched_randn_noise_async!(execution_style(cube), det, scratch, rng)
    sigma = effective_readout_sigma(det.params.sensor, det.noise.sigma)
    cube .+= sigma .* scratch
    return cube
end

function _batched_readout_noise!(det::Detector{<:NoisePhotonReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    _batched_randn_noise_async!(execution_style(cube), det, scratch, rng)
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
    _batched_poisson_noise_async!(execution_style(cube), det, cube, rng)
    return cube
end

function capture_stack_poisson_noise!(det::Detector{<:NoisePhotonReadout}, cube::AbstractArray, rng::AbstractRNG)
    _batched_poisson_noise_async!(execution_style(cube), det, cube, rng)
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

function _batched_apply_response!(::ScalarCPUStyle, model::SampledFrameResponse,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    n_batch, n, m = size(cube)
    kn, km = size(model.kernel)
    radius_i = fld(kn, 2)
    radius_j = fld(km, 2)
    @inbounds for j in 1:m, i in 1:n, b in 1:n_batch
        acc = zero(T)
        for ki in 1:kn
            ii = clamp(i + ki - radius_i - 1, 1, n)
            for kj in 1:km
                jj = clamp(j + kj - radius_j - 1, 1, m)
                acc += model.kernel[ki, kj] * cube[b, ii, jj]
            end
        end
        scratch[b, i, j] = acc
    end
    copyto!(cube, scratch)
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

function _batched_apply_response!(style::AcceleratorStyle, model::SampledFrameResponse,
    cube::AbstractArray{T,3}, scratch::AbstractArray{T,3}) where {T}
    n_batch, n, m = size(cube)
    kn, km = size(model.kernel)
    radius_i = fld(kn, 2)
    radius_j = fld(km, 2)
    launch_kernel_async!(style, sampled_response_stack_kernel!, scratch, cube, model.kernel,
        radius_i, radius_j, n_batch, n, m, kn, km; ndrange=size(cube))
    copyto!(cube, scratch)
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

function _capture_stack_fixed!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3},
    rng::AbstractRNG) where {T<:AbstractFloat}
    _require_batched_detector_compat(det, cube, scratch)
    return _apply_batched_detector_pipeline!(det, cube, scratch, rng)
end

function _capture_stack_fixed!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    return _capture_stack_fixed!(det, cube, scratch, rng)
end

function _apply_batched_detector_pipeline!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3},
    rng::AbstractRNG) where {T<:AbstractFloat}
    exposure_time = det.params.integration_time
    qe_scale = det.params.qe * exposure_time
    isone(qe_scale) || (cube .*= qe_scale)
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
    _batched_apply_readout_correction!(det.params.correction_model, cube, scratch)
    _batched_quantization!(det, cube)
    _batched_background_map!(det.background_map, cube, scratch)
    synchronize_backend!(execution_style(cube))
    return cube
end

function _apply_batched_detector_pipeline!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    return _apply_batched_detector_pipeline!(det, cube, scratch, rng)
end

function _require_generalized_batched_detector_compat(det::Detector, out_cube::AbstractArray, in_cube::AbstractArray)
    ndims(in_cube) == 3 || throw(DimensionMismatchError("generalized batched detector input must be 3D"))
    ndims(out_cube) == 3 || throw(DimensionMismatchError("generalized batched detector output must be 3D"))
    size(in_cube, 1) == size(out_cube, 1) ||
        throw(DimensionMismatchError("generalized batched detector batch count must match"))
    is_global_shutter(det.params.timing_model) ||
        throw(InvalidConfiguration("generalized batched detector capture currently requires global-shutter timing semantics"))
    is_null_persistence(persistence_model(det.params.sensor)) ||
        throw(InvalidConfiguration("generalized batched detector capture currently requires null detector persistence"))
    expected_shape = detector_output_shape(det, (size(in_cube, 2), size(in_cube, 3)))
    size(out_cube, 2) == expected_shape[1] && size(out_cube, 3) == expected_shape[2] ||
        throw(DimensionMismatchError("generalized batched detector output stack shape must match detector output shape"))
    return nothing
end

function _capture_stack_generalized!(det::Detector, out_cube::AbstractArray{TO,3}, in_cube::AbstractArray{TI,3};
    rng::AbstractRNG=Random.default_rng()) where {TO,TI}
    _require_generalized_batched_detector_compat(det, out_cube, in_cube)
    for b in axes(in_cube, 1)
        input_frame = @view(in_cube[b, :, :])
        output_frame = @view(out_cube[b, :, :])
        capture_signal!(det, input_frame, rng, det.params.integration_time)
        finalize_capture!(det, rng, det.params.integration_time)
        copyto!(output_frame, write_output!(det))
    end
    return out_cube
end

function capture_stack!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{S,3},
    rng::AbstractRNG) where {T<:AbstractFloat,S<:AbstractFloat}
    if size(cube) == size(scratch)
        return _capture_stack_fixed!(det, cube, scratch; rng=rng)
    end
    return _capture_stack_generalized!(det, cube, scratch; rng=rng)
end

function capture_stack!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{S,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat,S<:AbstractFloat}
    return capture_stack!(det, cube, scratch, rng)
end

function capture_stack!(det::Detector, out_cube::AbstractArray{TO,3}, in_cube::AbstractArray{TI,3},
    rng::AbstractRNG) where {TO,TI<:AbstractFloat}
    return _capture_stack_generalized!(det, out_cube, in_cube; rng=rng)
end

function capture_stack!(det::Detector, out_cube::AbstractArray{TO,3}, in_cube::AbstractArray{TI,3};
    rng::AbstractRNG=Random.default_rng()) where {TO,TI<:AbstractFloat}
    return capture_stack!(det, out_cube, in_cube, rng)
end

function apply_saturation!(det::Detector, cube::AbstractArray)
    full_well = det.params.full_well
    full_well === nothing && return cube
    clamp!(cube, zero(eltype(cube)), full_well)
    return cube
end
