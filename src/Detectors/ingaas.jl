struct InGaAsSensor{T<:AbstractFloat,P<:AbstractPersistenceModel} <: FrameSensorType
    glow_rate::T
    persistence_model::P
end

function InGaAsSensor(; glow_rate::Real=0.0, persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64)
    glow_rate >= 0 || throw(InvalidConfiguration("InGaAsSensor glow_rate must be >= 0"))
    converted_persistence = convert_persistence_model(persistence_model, T)
    validated_persistence = validate_persistence_model(converted_persistence)
    return InGaAsSensor{T,typeof(validated_persistence)}(T(glow_rate), validated_persistence)
end

detector_sensor_symbol(::InGaAsSensor) = :ingaas
supports_sensor_glow(::InGaAsSensor) = true
supports_detector_defect_maps(::InGaAsSensor) = true
supports_detector_persistence(::InGaAsSensor) = true
supports_detector_nonlinearity(::InGaAsSensor) = true
persistence_model(sensor::InGaAsSensor) = sensor.persistence_model

convert_persistence_model(::NullPersistence, ::Type{T}) where {T<:AbstractFloat} = NullPersistence()
convert_persistence_model(model::ExponentialPersistence, ::Type{T}) where {T<:AbstractFloat} =
    ExponentialPersistence{T}(T(model.coupling), T(model.decay))
validate_persistence_model(::NullPersistence) = NullPersistence()

function validate_persistence_model(model::ExponentialPersistence)
    zero(model.coupling) <= model.coupling <= one(model.coupling) ||
        throw(InvalidConfiguration("ExponentialPersistence coupling must lie in [0, 1]"))
    zero(model.decay) <= model.decay <= one(model.decay) ||
        throw(InvalidConfiguration("ExponentialPersistence decay must lie in [0, 1]"))
    return model
end

function apply_sensor_statistics!(sensor::InGaAsSensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

apply_sensor_persistence!(::InGaAsSensor{T,NullPersistence}, det::Detector, exposure_time::Real) where {T} = det.state.frame

function apply_sensor_persistence!(sensor::InGaAsSensor{T,<:ExponentialPersistence}, det::Detector, exposure_time::Real) where {T}
    ensure_latent_buffer!(det)
    det.state.frame .+= det.state.latent_buffer
    return det.state.frame
end

update_sensor_persistence!(::InGaAsSensor{T,NullPersistence}, det::Detector, exposure_time::Real) where {T} = det.state.frame

function update_sensor_persistence!(sensor::InGaAsSensor{T,<:ExponentialPersistence}, det::Detector, exposure_time::Real) where {T}
    ensure_latent_buffer!(det)
    model = sensor.persistence_model
    det.state.latent_buffer .= model.decay .* det.state.latent_buffer .+ model.coupling .* det.state.frame
    return det.state.latent_buffer
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
