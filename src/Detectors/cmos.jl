struct CMOSSensor{T<:AbstractFloat} <: FrameSensorType
    column_readout_sigma::T
end

function CMOSSensor(; column_readout_sigma::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    column_readout_sigma >= 0 || throw(InvalidConfiguration("CMOSSensor column_readout_sigma must be >= 0"))
    return CMOSSensor{T}(T(column_readout_sigma))
end

detector_sensor_symbol(::CMOSSensor) = :cmos
supports_column_readout_noise(::CMOSSensor) = true

function apply_sensor_statistics!(sensor::CMOSSensor, det::Detector, rng::AbstractRNG)
    sigma = sensor.column_readout_sigma
    sigma <= zero(sigma) && return det.state.frame
    randn_backend!(rng, det.state.noise_buffer)
    return apply_column_noise!(execution_style(det.state.frame), det.state.frame, det.state.noise_buffer, sigma)
end

function apply_column_noise!(::ScalarCPUStyle, frame, noise, sigma)
    n, m = size(frame)
    @inbounds for j in 1:m
        offset = sigma * noise[1, j]
        for i in 1:n
            frame[i, j] += offset
        end
    end
    return frame
end

function apply_column_noise!(style::AcceleratorStyle, frame, noise, sigma)
    n, m = size(frame)
    launch_kernel!(style, add_column_noise_kernel!, frame, noise, sigma, n, m; ndrange=(n, m))
    return frame
end

function apply_post_readout_gain!(::CMOSSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function _require_batched_sensor_compat(sensor::CMOSSensor)
    sensor.column_readout_sigma <= zero(sensor.column_readout_sigma) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor column_readout_sigma"))
    return nothing
end

_batched_post_readout_gain!(::CMOSSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
