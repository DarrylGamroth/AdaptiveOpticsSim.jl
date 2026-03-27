apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG) =
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)

function apply_pre_readout_gain!(::EMCCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(::EMCCDSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) =
    _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
