function sensor_saturation_limit(sensor::SAPHIRASensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    if rate > zero(rate)
        fill!(det.state.noise_buffer, rate)
        poisson_noise!(rng, det.state.noise_buffer)
        det.state.frame .+= det.state.noise_buffer
    end
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

function apply_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

function apply_post_readout_gain!(::SAPHIRASensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise!(rng, scratch)
        cube .+= scratch
    end
    return _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
end

_batched_post_readout_gain!(::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
