abstract type AbstractHgCdTeAvalancheMultiplication end

"""
    ConditionalGammaAvalancheMultiplication()

Nonnegative CPU reference model for independent HgCdTe linear-avalanche
multiplication. Conditional on input charge `q`, the input-referred multiplied
charge is Gamma distributed with shape `q / (F - 1)` and scale `F - 1`, where
`F` is the detector-literature excess-noise factor. `F == 1` is deterministic.
"""
struct ConditionalGammaAvalancheMultiplication <:
    AbstractHgCdTeAvalancheMultiplication end

"""
    ClippedGaussianAvalancheMultiplicationApproximation()

Backend-portable conditional-moment approximation for HgCdTe linear-avalanche
multiplication. Its qualified regime is `q / (F - 1) >= 25`; lower-charge
output is defined but is not claimed as a faithful multiplication
distribution.
"""
struct ClippedGaussianAvalancheMultiplicationApproximation <:
    AbstractHgCdTeAvalancheMultiplication end

"""
    HgCdTeAvalancheArraySensor(; avalanche_gain=1,
        excess_noise_factor=1,
        multiplication_model=
            ClippedGaussianAvalancheMultiplicationApproximation(),
        kwargs...)

HgCdTe linear-avalanche photodiode array. `excess_noise_factor` uses
`F = E[M^2] / E[M]^2`, so independent input charge `q` has input-referred
multiplication variance `(F - 1)q`. Nondestructive-read configuration is
composed through `HgCdTeReadout`; conventional HgCdTe operation should use
`HgCdTeSensor`.
"""
struct HgCdTeAvalancheArraySensor{
    T<:AbstractFloat,
    M<:AbstractHgCdTeAvalancheMultiplication,
    R<:HgCdTeReadout,
    P<:AbstractPersistenceModel,
} <: AbstractHgCdTeAvalancheArraySensor
    avalanche_gain::T
    excess_noise_factor::T
    multiplication_model::M
    readout::R
    persistence_model::P
end

function HgCdTeAvalancheArraySensor(
    readout::HgCdTeReadout{T};
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    multiplication_model::AbstractHgCdTeAvalancheMultiplication=
        ClippedGaussianAvalancheMultiplicationApproximation(),
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
        T,typeof(multiplication_model),typeof(readout),typeof(persistence)}(
        gain, noise_factor, multiplication_model, readout, persistence)
end

function HgCdTeAvalancheArraySensor(;
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    multiplication_model::AbstractHgCdTeAvalancheMultiplication=
        ClippedGaussianAvalancheMultiplicationApproximation(),
    glow_rate::Real=0.0,
    read_duration::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64,
)
    readout = HgCdTeReadout(
        glow_rate=glow_rate,
        read_duration=read_duration,
        sampling_mode=sampling_mode,
        T=T,
    )
    return HgCdTeAvalancheArraySensor(readout;
        avalanche_gain=avalanche_gain,
        excess_noise_factor=excess_noise_factor,
        multiplication_model=multiplication_model,
        persistence_model=persistence_model)
end

function owned_frame_sensor(sensor::HgCdTeAvalancheArraySensor, ::Type{T},
    backend::AbstractArrayBackend) where {T<:AbstractFloat}
    readout = HgCdTeReadout(
        glow_rate=sensor.readout.glow_rate,
        read_duration=sensor.readout.read_duration,
        sampling_mode=sensor.readout.sampling_mode,
        T=T,
    )
    owned = HgCdTeAvalancheArraySensor(readout;
        avalanche_gain=sensor.avalanche_gain,
        excess_noise_factor=sensor.excess_noise_factor,
        multiplication_model=sensor.multiplication_model,
        persistence_model=sensor.persistence_model)
    validate_hgcdte_avalanche_backend(owned.multiplication_model, backend)
    return owned
end

validate_hgcdte_avalanche_backend(
    ::AbstractHgCdTeAvalancheMultiplication, ::AbstractArrayBackend) = nothing

validate_hgcdte_avalanche_backend(
    ::ConditionalGammaAvalancheMultiplication, ::CPUBackend) = nothing

function validate_hgcdte_avalanche_backend(
    ::ConditionalGammaAvalancheMultiplication,
    backend::AbstractArrayBackend)
    throw(InvalidConfiguration(
        "ConditionalGammaAvalancheMultiplication is CPU-only; use " *
        "ClippedGaussianAvalancheMultiplicationApproximation for backend " *
        "$(typeof(backend))"))
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

function _apply_hgcdte_avalanche_statistics!(
    ::ClippedGaussianAvalancheMultiplicationApproximation,
    sensor::HgCdTeAvalancheArraySensor, frame::AbstractArray{T},
    scratch::AbstractArray, rng::AbstractRNG) where {T<:AbstractFloat}
    factor_minus_one = T(sensor.excess_noise_factor - one(T))
    factor_minus_one <= zero(T) && return frame
    randn_backend!(rng, scratch)
    zero_t = zero(T)
    @. frame = max(
        frame + sqrt(max(factor_minus_one * frame, zero_t)) * scratch,
        zero_t)
    return frame
end

function _apply_hgcdte_avalanche_statistics!(
    ::ScalarCPUStyle, ::ConditionalGammaAvalancheMultiplication,
    sensor::HgCdTeAvalancheArraySensor, frame::AbstractArray{T},
    ::AbstractArray, rng::AbstractRNG) where {T<:AbstractFloat}
    factor_minus_one = T(sensor.excess_noise_factor - one(T))
    factor_minus_one <= zero(T) && return frame
    @inbounds for index in eachindex(frame)
        charge = max(frame[index], zero(T))
        frame[index] = charge <= zero(T) ? zero(T) :
            factor_minus_one *
            _gamma_unit_scale(rng, charge / factor_minus_one)
    end
    return frame
end

function _apply_hgcdte_avalanche_statistics!(
    ::AcceleratorStyle, ::ConditionalGammaAvalancheMultiplication,
    ::HgCdTeAvalancheArraySensor, ::AbstractArray, ::AbstractArray,
    ::AbstractRNG)
    throw(InvalidConfiguration(
        "ConditionalGammaAvalancheMultiplication is CPU-only"))
end

function _apply_hgcdte_avalanche_statistics!(
    model::ConditionalGammaAvalancheMultiplication,
    sensor::HgCdTeAvalancheArraySensor, frame::AbstractArray,
    scratch::AbstractArray, rng::AbstractRNG)
    return _apply_hgcdte_avalanche_statistics!(
        execution_style(frame), model, sensor, frame, scratch, rng)
end

function apply_hgcdte_avalanche_statistics!(
    sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG)
    return _apply_hgcdte_avalanche_statistics!(
        sensor.multiplication_model, sensor, det.products.frame,
        det.workspace.noise_buffer, rng)
end

function apply_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, rng::AbstractRNG, exposure_duration::Real)
    _apply_hgcdte_glow!(sensor, det, rng, exposure_duration)
    return apply_hgcdte_avalanche_statistics!(sensor, det, rng)
end

function apply_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, ::AbstractRNG)
    det.products.frame .*= sensor.avalanche_gain
    return det.products.frame
end

_batched_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor,
    ::Detector, cube::AbstractArray, ::AbstractArray, ::AbstractRNG) =
    (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(
    sensor::HgCdTeAvalancheArraySensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_duration::Real)
    _batched_hgcdte_glow!(
        sensor, det, cube, scratch, rng, exposure_duration)
    return _apply_hgcdte_avalanche_statistics!(
        sensor.multiplication_model, sensor, cube, scratch, rng)
end
