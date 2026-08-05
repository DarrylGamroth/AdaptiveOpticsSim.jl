"""
    InGaAsSensor(; glow_rate=0, persistence_model=NullPersistence(), T=Float64)

Product-neutral InGaAs area-sensor model. `glow_rate` is a configured glow
charge rate in electrons per pixel per second of integration. It does not model
read-cadence-dependent glow. A presampling response is never inferred from the
detector material; configure one on the owning `Detector` when required.

`ExponentialPersistence` is an optional charge-domain, frame-to-frame reduced
model. It is updated after charge nonlinearity, saturation, and coupling and
before read noise, conversion gain, correction, quantization, and background
subtraction.
"""
struct InGaAsSensor{T<:AbstractFloat,P<:AbstractPersistenceModel} <: AbstractFrameSensor
    glow_rate::T
    persistence_model::P
end

function InGaAsSensor(; glow_rate::Real=0.0, persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64)
    glow = T(glow_rate)
    isfinite(glow) && glow >= zero(T) || throw(InvalidConfiguration(
        "InGaAsSensor glow_rate must be finite and >= 0"))
    converted_persistence = convert_persistence_model(persistence_model, T)
    validated_persistence = validate_persistence_model(converted_persistence)
    return InGaAsSensor{T,typeof(validated_persistence)}(
        glow, validated_persistence)
end

detector_sensor_symbol(::InGaAsSensor) = :ingaas
supports_sensor_glow(::InGaAsSensor) = true
supports_detector_defect_maps(::InGaAsSensor) = true
supports_detector_persistence(::InGaAsSensor) = true
supports_detector_nonlinearity(::InGaAsSensor) = true
default_response_model(::InGaAsSensor;
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend()) = NullFrameResponse()
persistence_model(sensor::InGaAsSensor) = sensor.persistence_model
configured_glow_rate(sensor::InGaAsSensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.glow_rate)

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

function apply_sensor_statistics!(sensor::InGaAsSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_duration(sensor, exposure_duration)
    rate <= zero(rate) && return det.products.frame
    fill!(det.workspace.noise_buffer, rate)
    poisson_noise!(rng, det.workspace.noise_buffer)
    det.products.frame .+= det.workspace.noise_buffer
    return det.products.frame
end

function apply_incremental_sensor_statistics!(sensor::InGaAsSensor,
    det::Detector, rng::AbstractRNG, exposure_duration::Real)
    rate = effective_glow_rate(det) * exposure_duration
    return add_poisson_rate!(det.products.frame, det, rng, rate)
end

apply_sensor_persistence!(::InGaAsSensor{T,NullPersistence}, det::Detector, exposure_duration::Real) where {T} = det.products.frame

function apply_sensor_persistence!(sensor::InGaAsSensor{T,<:ExponentialPersistence}, det::Detector, exposure_duration::Real) where {T}
    ensure_latent_buffer!(det)
    det.products.frame .+= det.state.latent_buffer
    return det.products.frame
end

update_sensor_persistence!(::InGaAsSensor{T,NullPersistence}, det::Detector, exposure_duration::Real) where {T} = det.products.frame

function update_sensor_persistence!(sensor::InGaAsSensor{T,<:ExponentialPersistence}, det::Detector, exposure_duration::Real) where {T}
    ensure_latent_buffer!(det)
    model = sensor.persistence_model
    det.state.latent_buffer .= model.decay .* det.state.latent_buffer .+ model.coupling .* det.products.frame
    return det.state.latent_buffer
end

function _finalize_capture!(::InGaAsSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    return finalize_ingaas_capture!(det, rng, exposure_duration, exposure_duration)
end

function _finalize_incremental_capture!(::InGaAsSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    return finalize_ingaas_capture!(det, rng, exposure_duration,
        zero(exposure_duration))
end

function finalize_ingaas_capture!(det::Detector, rng::AbstractRNG,
    exposure_duration::Real, charge_exposure_duration::Real)
    finalize_charge_generation!(det, rng, charge_exposure_duration)
    finalize_charge_transport!(det, rng)
    update_sensor_persistence!(det.params.sensor, det, exposure_duration)
    return finalize_electronics_without_persistence!(det, rng, exposure_duration)
end

function apply_post_readout_gain!(::InGaAsSensor, det::Detector)
    det.products.frame .*= det.params.gain
    return det.products.frame
end

function _batched_sensor_statistics!(sensor::InGaAsSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_duration::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_duration(sensor, exposure_duration)
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_post_readout_gain!(::InGaAsSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
