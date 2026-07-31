"""
    HgCdTeAvalancheArraySensor(; avalanche_gain=1,
        excess_noise_factor=1, kwargs...)

HgCdTe linear-avalanche photodiode array. Nondestructive-read configuration is
composed through `HgCdTeReadout`; conventional HgCdTe operation should use
`HgCdTeSensor` and does not require an avalanche-named type.
"""
struct HgCdTeAvalancheArraySensor{
    T<:AbstractFloat,
    R<:HgCdTeReadout,
    P<:AbstractPersistenceModel,
} <: HgCdTeAvalancheArraySensorType
    avalanche_gain::T
    excess_noise_factor::T
    readout::R
    persistence_model::P
end

function HgCdTeAvalancheArraySensor(
    readout::HgCdTeReadout{T};
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    persistence_model::AbstractPersistenceModel=NullPersistence(),
) where {T<:AbstractFloat}
    gain = T(avalanche_gain)
    noise_factor = T(excess_noise_factor)
    isfinite(gain) && gain >= one(T) || throw(InvalidConfiguration(
        "HgCdTeAvalancheArraySensor avalanche_gain must be finite and >= 1"))
    isfinite(noise_factor) && noise_factor >= one(T) ||
        throw(InvalidConfiguration(
            "HgCdTeAvalancheArraySensor excess_noise_factor must be finite and >= 1"))
    persistence = prepare_hgcdte_persistence_model(persistence_model, T)
    return HgCdTeAvalancheArraySensor{
        T,typeof(readout),typeof(persistence)}(
        gain, noise_factor, readout, persistence)
end

function HgCdTeAvalancheArraySensor(;
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    glow_rate::Real=0.0,
    read_time::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64,
)
    readout = HgCdTeReadout(
        glow_rate=glow_rate,
        read_time=read_time,
        sampling_mode=sampling_mode,
        T=T,
    )
    return HgCdTeAvalancheArraySensor(readout;
        avalanche_gain=avalanche_gain,
        excess_noise_factor=excess_noise_factor,
        persistence_model=persistence_model)
end

@inline hgcdte_readout(sensor::HgCdTeAvalancheArraySensor) = sensor.readout

detector_sensor_symbol(::HgCdTeAvalancheArraySensor) =
    :hgcdte_linear_avalanche_array

function sensor_saturation_limit(
    sensor::HgCdTeAvalancheArraySensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    _apply_hgcdte_glow!(sensor, det, rng, exposure_time)
    return apply_avalanche_excess_noise!(
        sensor.excess_noise_factor, det, rng)
end

function apply_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, ::AbstractRNG)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

_batched_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor,
    ::Detector, cube::AbstractArray, ::AbstractArray, ::AbstractRNG) =
    (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(
    sensor::HgCdTeAvalancheArraySensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_time::Real)
    _batched_hgcdte_glow!(
        sensor, det, cube, scratch, rng, exposure_time)
    return _batched_avalanche_excess_noise!(
        sensor.excess_noise_factor, cube, scratch, rng)
end
