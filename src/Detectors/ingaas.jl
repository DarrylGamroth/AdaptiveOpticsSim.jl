struct InGaAsSensor{T<:AbstractFloat} <: FrameSensorType
    glow_rate::T
end

function InGaAsSensor(; glow_rate::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    glow_rate >= 0 || throw(InvalidConfiguration("InGaAsSensor glow_rate must be >= 0"))
    return InGaAsSensor{T}(T(glow_rate))
end

detector_sensor_symbol(::InGaAsSensor) = :ingaas
supports_sensor_glow(::InGaAsSensor) = true

function apply_sensor_statistics!(sensor::InGaAsSensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_post_readout_gain!(::InGaAsSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::InGaAsSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_post_readout_gain!(::InGaAsSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
