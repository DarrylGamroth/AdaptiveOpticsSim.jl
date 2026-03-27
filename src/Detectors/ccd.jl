struct CCDSensor{T<:AbstractFloat} <: FrameSensorType
    clock_induced_charge_rate::T
end

function CCDSensor(; clock_induced_charge_rate::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    clock_induced_charge_rate >= 0 || throw(InvalidConfiguration("CCDSensor clock_induced_charge_rate must be >= 0"))
    return CCDSensor{T}(T(clock_induced_charge_rate))
end

detector_sensor_symbol(::CCDSensor) = :ccd
supports_clock_induced_charge(::CCDSensor) = true

function apply_sensor_statistics!(sensor::CCDSensor, det::Detector, rng::AbstractRNG)
    rate = sensor.clock_induced_charge_rate * det.params.integration_time
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_post_readout_gain!(::CCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::CCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.clock_induced_charge_rate * det.params.integration_time
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_post_readout_gain!(::CCDSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
