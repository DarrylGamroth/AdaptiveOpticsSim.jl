struct EMCCDSensor{T<:AbstractFloat} <: AvalancheFrameSensorType
    excess_noise_factor::T
end

function EMCCDSensor(; excess_noise_factor::Real=1.0, T::Type{<:AbstractFloat}=Float64)
    excess_noise_factor >= 1 || throw(InvalidConfiguration("EMCCDSensor excess_noise_factor must be >= 1"))
    return EMCCDSensor{T}(T(excess_noise_factor))
end

detector_sensor_symbol(::EMCCDSensor) = :emccd

apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG) =
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)

function apply_pre_readout_gain!(::EMCCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(::EMCCDSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) =
    _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
