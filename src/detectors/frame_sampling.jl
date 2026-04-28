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

function _ensure_multi_read_products!(sensor::FrameSensorType, det::Detector)
    current = readout_products(det)
    frame = det.state.frame
    n_ref = frame_sampling_reference_reads(sensor)
    n_sig = frame_sampling_signal_reads(sensor)
    n_reads = frame_sampling_reads(sensor)

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
    fill!(baseline, zero(eltype(baseline)))
    for read_idx in 1:n_ref
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

function _sampling_read_cube!(cube::AbstractArray{T,3}, sensor::FrameSensorType,
    reference_cube::Union{Nothing,AbstractArray{T,3}}, signal_cube::AbstractArray{T,3}) where {T}
    offset = 0
    if !isnothing(reference_cube)
        for read_idx in axes(reference_cube, 3)
            @views copyto!(cube[:, :, read_idx], reference_cube[:, :, read_idx])
        end
        offset = size(reference_cube, 3)
    end
    for read_idx in axes(signal_cube, 3)
        @views copyto!(cube[:, :, offset + read_idx], signal_cube[:, :, read_idx])
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
