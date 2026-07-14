@inline multi_read_sampling_mode(sensor::FrameSensorType) = SingleRead()

_raw_sampling_sigma(det::Detector{<:NoiseReadout}) = det.noise.sigma
_raw_sampling_sigma(det::Detector{<:NoisePhotonReadout}) = det.noise.sigma
_raw_sampling_sigma(det::Detector) = zero(eltype(det.state.frame))

function _sampling_average_cube!(out::AbstractMatrix{T}, cube::AbstractArray{T,3}) where {T}
    n_reads = size(cube, 3)
    fill!(out, zero(T))
    for read_idx in 1:n_reads
        @views out .+= cube[:, :, read_idx]
    end
    out ./= T(n_reads)
    return out
end

function _sampling_average_cube(cube::AbstractArray{T,3}) where {T}
    out = similar(cube, size(cube, 1), size(cube, 2))
    return _sampling_average_cube!(out, cube)
end

function _ensure_windowed_frame_buffer(current, det::Detector, full_frame::AbstractMatrix{T}) where {T}
    target_shape = readout_product_shape(det)
    if current === nothing || size(current) != target_shape
        return similar(full_frame, target_shape...)
    end
    return current
end

function _ensure_windowed_cube_buffer(current, det::Detector, full_frame::AbstractMatrix{T}, n_reads::Int) where {T}
    target_shape = (readout_product_shape(det)..., n_reads)
    if current === nothing || size(current) != target_shape
        return similar(full_frame, target_shape...)
    end
    return current
end

function _ensure_read_times_buffer(current, ::Type{T}, n_reads::Int) where {T<:AbstractFloat}
    if current === nothing || length(current) != n_reads || eltype(current) !== T
        return Vector{T}(undef, n_reads)
    end
    return current
end

function _copy_windowed_frame!(target::AbstractMatrix, full_frame::AbstractMatrix, det::Detector)
    window = det.params.readout_window
    if window === nothing
        copyto!(target, full_frame)
    else
        copyto!(target, @view(full_frame[window.rows, window.cols]))
    end
    return target
end

function _copy_windowed_cube!(target::AbstractArray{T,3}, full_cube::AbstractArray{T,3}, det::Detector) where {T}
    window = det.params.readout_window
    if window === nothing
        copyto!(target, full_cube)
    else
        copyto!(target, @view(full_cube[window.rows, window.cols, :]))
    end
    return target
end

_multi_read_products_ready(::FrameReadoutProducts, ::Detector,
    ::Int, ::Int, ::Int) = false

function _multi_read_products_ready(products::MultiReadFrameReadoutProducts,
    det::Detector, n_ref::Int, n_sig::Int, n_reads::Int)
    frame_shape = readout_product_shape(det)
    reference_frame_ready = n_ref <= 0 ? products.reference_frame === nothing :
        products.reference_frame !== nothing && size(products.reference_frame) == frame_shape
    reference_cube_ready = n_ref <= 0 ? products.reference_cube === nothing :
        products.reference_cube !== nothing && size(products.reference_cube) == (frame_shape..., n_ref)
    read_cube_ready = n_reads <= 1 ? products.read_cube === nothing :
        products.read_cube !== nothing && size(products.read_cube) == (frame_shape..., n_reads)
    read_times_ready = n_reads <= 1 ? products.read_times === nothing :
        products.read_times !== nothing && length(products.read_times) == n_reads &&
        eltype(products.read_times) === eltype(det.state.frame)
    return reference_frame_ready &&
        size(products.signal_frame) == frame_shape &&
        size(products.combined_frame) == frame_shape &&
        reference_cube_ready &&
        products.signal_cube !== nothing &&
        size(products.signal_cube) == (frame_shape..., n_sig) &&
        read_cube_ready && read_times_ready
end

function _ensure_multi_read_products!(sensor::FrameSensorType, det::Detector)
    current = readout_products(det)
    frame = det.state.frame
    n_ref = frame_sampling_reference_reads(sensor)
    n_sig = frame_sampling_signal_reads(sensor)
    n_reads = frame_sampling_reads(sensor)
    _multi_read_products_ready(current, det, n_ref, n_sig, n_reads) && return current

    reference_frame = detector_reference_frame(current)
    signal_frame = detector_signal_frame(current)
    combined_frame = detector_combined_frame(current)
    reference_cube = detector_reference_cube(current)
    signal_cube = detector_signal_cube(current)
    read_cube = detector_read_cube(current)
    read_times = detector_read_times(current)

    reference_frame = n_ref <= 0 ? nothing : _ensure_windowed_frame_buffer(reference_frame, det, frame)
    signal_frame = _ensure_windowed_frame_buffer(signal_frame, det, frame)
    combined_frame = _ensure_windowed_frame_buffer(combined_frame, det, frame)
    reference_cube = n_ref <= 0 ? nothing : _ensure_windowed_cube_buffer(reference_cube, det, frame, n_ref)
    signal_cube = _ensure_windowed_cube_buffer(signal_cube, det, frame, n_sig)
    read_cube = n_reads <= 1 ? nothing : _ensure_windowed_cube_buffer(read_cube, det, frame, n_reads)
    T = eltype(frame)
    read_times = n_reads <= 1 ? nothing : _ensure_read_times_buffer(read_times, T, n_reads)

    products = MultiReadFrameReadoutProducts(reference_frame, signal_frame, combined_frame,
        reference_cube, signal_cube, read_cube, read_times)
    det.state.readout_products = products
    return products
end

function sample_frame_read!(::FrameSensorType, det::Detector, target::AbstractMatrix, baseline::AbstractMatrix,
    sigma, rng::AbstractRNG)
    copyto!(target, baseline)
    add_gaussian_noise!(target, det, rng, sigma)
    target .*= det.params.gain
    apply_readout_correction!(det.params.correction_model, target, det)
    return target
end

function _sampling_reference_cube!(cube::AbstractArray{T,3}, sensor::FrameSensorType,
    det::Detector, sigma, rng::AbstractRNG, baseline::AbstractMatrix{T}) where {T}
    n_ref = frame_sampling_reference_reads(sensor)
    for read_idx in 1:n_ref
        fill!(baseline, zero(eltype(baseline)))
        sample_frame_read!(sensor, det, baseline, baseline, sigma, rng)
        @views copyto!(cube[:, :, read_idx], baseline)
    end
    return cube
end

function _sampling_reference_cube(sensor::FrameSensorType, det::Detector, sigma, rng::AbstractRNG)
    n_ref = frame_sampling_reference_reads(sensor)
    n_ref <= 0 && return nothing
    cube = similar(det.state.frame, size(det.state.frame)..., n_ref)
    baseline = similar(det.state.frame, size(det.state.frame)...)
    return _sampling_reference_cube!(cube, sensor, det, sigma, rng, baseline)
end

function _sampling_signal_cube!(cube::AbstractArray{T,3}, sensor::FrameSensorType,
    det::Detector, sigma, rng::AbstractRNG, baseline::AbstractMatrix{T}) where {T}
    n_sig = frame_sampling_signal_reads(sensor)
    for read_idx in 1:n_sig
        sample_frame_read!(sensor, det, baseline, det.state.frame, sigma, rng)
        @views copyto!(cube[:, :, read_idx], baseline)
    end
    return cube
end

function _sampling_signal_cube(sensor::FrameSensorType, det::Detector, sigma, rng::AbstractRNG)
    n_sig = frame_sampling_signal_reads(sensor)
    cube = similar(det.state.frame, size(det.state.frame)..., n_sig)
    baseline = similar(det.state.frame, size(det.state.frame)...)
    return _sampling_signal_cube!(cube, sensor, det, sigma, rng, baseline)
end

function _copy_sampling_plane!(target::AbstractArray{T,3}, target_index::Int,
    source::AbstractArray{T,3}, source_index::Int) where {T}
    return _copy_sampling_plane!(execution_style(target), target, target_index,
        source, source_index)
end

function _copy_sampling_plane!(::ScalarCPUStyle, target::AbstractArray{T,3},
    target_index::Int, source::AbstractArray{T,3}, source_index::Int) where {T}
    @inbounds for j in axes(target, 2), i in axes(target, 1)
        target[i, j, target_index] = source[i, j, source_index]
    end
    return target
end

function _copy_sampling_plane!(::AcceleratorStyle, target::AbstractArray{T,3},
    target_index::Int, source::AbstractArray{T,3}, source_index::Int) where {T}
    @views copyto!(target[:, :, target_index], source[:, :, source_index])
    return target
end

function _sampling_read_cube!(cube::AbstractArray{T,3}, ::FrameSensorType,
    ::Nothing, signal_cube::AbstractArray{T,3}) where {T}
    for read_idx in axes(signal_cube, 3)
        _copy_sampling_plane!(cube, read_idx, signal_cube, read_idx)
    end
    return cube
end

function _sampling_read_cube!(cube::AbstractArray{T,3}, ::FrameSensorType,
    reference_cube::AbstractArray{T,3}, signal_cube::AbstractArray{T,3}) where {T}
    for read_idx in axes(reference_cube, 3)
        _copy_sampling_plane!(cube, read_idx, reference_cube, read_idx)
    end
    offset = size(reference_cube, 3)
    for read_idx in axes(signal_cube, 3)
        _copy_sampling_plane!(cube, offset + read_idx, signal_cube, read_idx)
    end
    return cube
end

function _sampling_read_cube(sensor::FrameSensorType, reference_cube::Union{Nothing,AbstractArray{T,3}},
    signal_cube::AbstractArray{T,3}) where {T}
    n_reads = frame_sampling_reads(sensor)
    n_reads <= 1 && return nothing
    cube = similar(signal_cube, size(signal_cube, 1), size(signal_cube, 2), n_reads)
    return _sampling_read_cube!(cube, sensor, reference_cube, signal_cube)
end

function _sampling_read_times!(times::AbstractVector{T}, sensor::FrameSensorType,
    det::Detector, n_reads::Int) where {T}
    read_dt = sampling_read_time(sensor, size(det.state.frame), det.params.readout_window, T)
    for read_idx in 1:n_reads
        times[read_idx] = T(read_idx) * read_dt
    end
    return times
end

function _sampling_read_times(sensor::FrameSensorType, det::Detector, n_reads::Int)
    n_reads <= 0 && return nothing
    T = eltype(det.state.frame)
    times = Vector{T}(undef, n_reads)
    return _sampling_read_times!(times, sensor, det, n_reads)
end

function finalize_multi_read_readout_products!(sensor::FrameSensorType, det::Detector, rng::AbstractRNG)
    sigma = _raw_sampling_sigma(det)
    if det.params.readout_window !== nothing
        reference_cube_full = _sampling_reference_cube(sensor, det, sigma, rng)
        signal_cube_full = _sampling_signal_cube(sensor, det, sigma, rng)
        cube_full = _sampling_read_cube(sensor, reference_cube_full, signal_cube_full)
        reference_full = isnothing(reference_cube_full) ? nothing : _sampling_average_cube(reference_cube_full)
        signal_full = _sampling_average_cube(signal_cube_full)
        combined_full = isnothing(reference_full) ? signal_full : signal_full .- reference_full
        reference = isnothing(reference_full) ? nothing : _copy_windowed_frame(reference_full, det)
        signal = _copy_windowed_frame(signal_full, det)
        combined = _copy_windowed_frame(combined_full, det)
        reference_cube = isnothing(reference_cube_full) ? nothing : _copy_windowed_cube(reference_cube_full, det)
        signal_cube = _copy_windowed_cube(signal_cube_full, det)
        cube = isnothing(cube_full) ? nothing : _copy_windowed_cube(cube_full, det)
        read_times = isnothing(cube_full) ? nothing : _sampling_read_times(sensor, det, size(cube_full, 3))
        det.state.frame .= combined_full
        det.state.readout_products = MultiReadFrameReadoutProducts(reference, signal, combined, reference_cube, signal_cube, cube, read_times)
        return det.state.readout_products
    end

    products = _ensure_multi_read_products!(sensor, det)
    baseline = det.state.response_buffer
    reference_average = det.state.accum_buffer
    signal_average = det.state.bin_buffer

    if !isnothing(products.reference_cube)
        _sampling_reference_cube!(products.reference_cube, sensor, det, sigma, rng, baseline)
        _sampling_average_cube!(reference_average, products.reference_cube)
        _copy_windowed_frame!(products.reference_frame, reference_average, det)
    end

    _sampling_signal_cube!(products.signal_cube, sensor, det, sigma, rng, baseline)
    _sampling_average_cube!(signal_average, products.signal_cube)
    _copy_windowed_frame!(products.signal_frame, signal_average, det)

    if !isnothing(products.reference_frame)
        products.combined_frame .= products.signal_frame .- products.reference_frame
        det.state.frame .= signal_average .- reference_average
    else
        copyto!(products.combined_frame, products.signal_frame)
        det.state.frame .= signal_average
    end

    if !isnothing(products.read_cube)
        _sampling_read_cube!(products.read_cube, sensor, products.reference_cube, products.signal_cube)
        _sampling_read_times!(products.read_times, sensor, det, size(products.read_cube, 3))
    end
    return products
end

function _up_the_ramp_products_ready(products::UpTheRampReadoutProducts,
    det::Detector, n_reads::Int)
    output_shape = readout_product_shape(det)
    frame_shape = size(det.state.frame)
    return size(products.slope_frame) == output_shape &&
        size(products.intercept_frame) == output_shape &&
        size(products.integrated_frame) == output_shape &&
        size(products.read_cube) == (output_shape..., n_reads) &&
        length(products.read_times) == n_reads &&
        size(products.workspace_slope) == frame_shape &&
        size(products.workspace_intercept) == frame_shape &&
        size(products.workspace_integrated) == frame_shape &&
        size(products.workspace_cube) == (frame_shape..., n_reads)
end

_up_the_ramp_products_ready(::FrameReadoutProducts, det::Detector,
    n_reads::Int) = false

function ensure_up_the_ramp_products!(det::Detector, n_reads::Int)
    current = readout_products(det)
    _up_the_ramp_products_ready(current, det, n_reads) && return current

    frame = det.state.frame
    output_shape = readout_product_shape(det)
    frame_shape = size(frame)
    slope_frame = similar(frame, output_shape...)
    intercept_frame = similar(frame, output_shape...)
    integrated_frame = similar(frame, output_shape...)
    read_cube = similar(frame, output_shape..., n_reads)
    read_times = Vector{eltype(frame)}(undef, n_reads)

    if det.params.readout_window === nothing
        workspace_slope = slope_frame
        workspace_intercept = intercept_frame
        workspace_integrated = integrated_frame
        workspace_cube = read_cube
    else
        workspace_slope = similar(frame, frame_shape...)
        workspace_intercept = similar(frame, frame_shape...)
        workspace_integrated = similar(frame, frame_shape...)
        workspace_cube = similar(frame, frame_shape..., n_reads)
    end

    products = UpTheRampReadoutProducts(slope_frame, intercept_frame,
        integrated_frame, read_cube, read_times, workspace_slope,
        workspace_intercept, workspace_integrated, workspace_cube)
    det.state.readout_products = products
    return products
end

function validate_up_the_ramp_schedule(sensor::FrameSensorType, det::Detector,
    mode::UpTheRampSampling, exposure_time::Real)
    T = eltype(det.state.frame)
    read_time = sampling_read_time(sensor, size(det.state.frame),
        det.params.readout_window, T)
    read_spacing = T(exposure_time) / T(mode.n_reads - 1)
    read_time <= read_spacing + eps(read_spacing) || throw(InvalidConfiguration(
        "up-the-ramp read_time must not exceed the spacing between reads"))
    return read_spacing
end

function _fill_up_the_ramp_times!(times::AbstractVector{T},
    exposure_time::Real) where {T<:AbstractFloat}
    n_reads = length(times)
    dt = T(exposure_time) / T(n_reads - 1)
    for read_idx in eachindex(times)
        times[read_idx] = T(read_idx - firstindex(times)) * dt
    end
    times[end] = T(exposure_time)
    return times
end

function _sample_up_the_ramp_cube!(cube::AbstractArray{T,3},
    sensor::FrameSensorType, det::Detector, sigma, rng::AbstractRNG) where {T}
    n_reads = size(cube, 3)
    denominator = T(n_reads - 1)
    for read_idx in 1:n_reads
        fraction = T(read_idx - 1) / denominator
        target = @view cube[:, :, read_idx]
        @. target = fraction * det.state.frame
        add_gaussian_noise!(target, det, rng, sigma)
        target .*= det.params.gain
        apply_readout_correction!(det.params.correction_model, target, det)
    end
    return cube
end

@kernel function fit_up_the_ramp_kernel!(integrated, slope, intercept, cube,
    dt, mean_time, inv_denominator, inv_n_reads, exposure_time,
    n_reads::Int, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        sum_y = zero(eltype(slope))
        centered_dot = zero(eltype(slope))
        @inbounds for read_idx in 1:n_reads
            centered_time = (read_idx - 1) * dt - mean_time
            value = cube[i, j, read_idx]
            sum_y += value
            centered_dot = muladd(centered_time, value, centered_dot)
        end
        fitted_slope = centered_dot * inv_denominator
        fitted_intercept = sum_y * inv_n_reads - fitted_slope * mean_time
        @inbounds begin
            slope[i, j] = fitted_slope
            intercept[i, j] = fitted_intercept
            integrated[i, j] = fitted_slope * exposure_time
        end
    end
end

function _up_the_ramp_fit_coefficients(::Type{T}, n_reads::Int,
    exposure_time::Real) where {T<:AbstractFloat}
    exposure = T(exposure_time)
    n = T(n_reads)
    dt = exposure / T(n_reads - 1)
    mean_time = exposure / T(2)
    denominator = exposure * exposure * n * (n + one(T)) /
        (T(12) * (n - one(T)))
    return dt, mean_time, inv(denominator), inv(n), exposure
end

function _fit_up_the_ramp!(::ScalarCPUStyle, integrated::AbstractMatrix{T},
    slope::AbstractMatrix{T}, intercept::AbstractMatrix{T},
    cube::AbstractArray{T,3}, exposure_time::Real) where {T<:AbstractFloat}
    n_reads = size(cube, 3)
    dt, mean_time, inv_denominator, inv_n_reads, exposure =
        _up_the_ramp_fit_coefficients(T, n_reads, exposure_time)
    @inbounds for j in axes(slope, 2), i in axes(slope, 1)
        sum_y = zero(T)
        centered_dot = zero(T)
        for read_idx in 1:n_reads
            centered_time = T(read_idx - 1) * dt - mean_time
            value = cube[i, j, read_idx]
            sum_y += value
            centered_dot = muladd(centered_time, value, centered_dot)
        end
        fitted_slope = centered_dot * inv_denominator
        slope[i, j] = fitted_slope
        intercept[i, j] = sum_y * inv_n_reads - fitted_slope * mean_time
        integrated[i, j] = fitted_slope * exposure
    end
    return integrated
end

function _fit_up_the_ramp!(style::AcceleratorStyle,
    integrated::AbstractMatrix{T}, slope::AbstractMatrix{T},
    intercept::AbstractMatrix{T}, cube::AbstractArray{T,3},
    exposure_time::Real) where {T<:AbstractFloat}
    n_reads = size(cube, 3)
    dt, mean_time, inv_denominator, inv_n_reads, exposure =
        _up_the_ramp_fit_coefficients(T, n_reads, exposure_time)
    n, m = size(slope)
    launch_kernel!(style, fit_up_the_ramp_kernel!, integrated, slope, intercept,
        cube, dt, mean_time, inv_denominator, inv_n_reads, exposure, n_reads,
        n, m; ndrange=(n, m))
    return integrated
end

function finalize_up_the_ramp_readout_products!(mode::UpTheRampSampling,
    sensor::FrameSensorType, det::Detector, rng::AbstractRNG,
    exposure_time::Real)
    validate_up_the_ramp_schedule(sensor, det, mode, exposure_time)
    products = ensure_up_the_ramp_products!(det, mode.n_reads)
    _fill_up_the_ramp_times!(products.read_times, exposure_time)
    _sample_up_the_ramp_cube!(products.workspace_cube, sensor, det,
        _raw_sampling_sigma(det), rng)
    _fit_up_the_ramp!(execution_style(products.workspace_integrated),
        products.workspace_integrated, products.workspace_slope,
        products.workspace_intercept, products.workspace_cube, exposure_time)

    if det.params.readout_window !== nothing
        _copy_windowed_frame!(products.slope_frame, products.workspace_slope, det)
        _copy_windowed_frame!(products.intercept_frame,
            products.workspace_intercept, det)
        _copy_windowed_frame!(products.integrated_frame,
            products.workspace_integrated, det)
        _copy_windowed_cube!(products.read_cube, products.workspace_cube, det)
    end

    copyto!(det.state.frame, products.workspace_integrated)
    return products
end
