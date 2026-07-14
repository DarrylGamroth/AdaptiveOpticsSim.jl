abstract type AbstractLinearAPDTopology end

"""A single linear-mode APD with scalar optical input and vector storage."""
struct SingleElementAPD <: AbstractLinearAPDTopology end

"""A fixed-size bank of independent linear-mode APD channels."""
struct APDChannelBank <: AbstractLinearAPDTopology
    n_channels::Int
    function APDChannelBank(n_channels::Integer)
        n_channels >= 2 || throw(InvalidConfiguration(
            "APDChannelBank n_channels must be >= 2; use SingleElementAPD for one channel"))
        return new(Int(n_channels))
    end
end

linear_apd_channel_count(::SingleElementAPD) = 1
linear_apd_channel_count(topology::APDChannelBank) = topology.n_channels
linear_apd_topology_symbol(::SingleElementAPD) = :single_element
linear_apd_topology_symbol(::APDChannelBank) = :channel_bank

struct LinearAPDDetectorParams{T<:AbstractFloat,TP<:AbstractLinearAPDTopology}
    integration_time::T
    qe::T
    avalanche_gain::T
    excess_noise_factor::T
    dark_current::T
    conversion_gain::T
    topology::TP
end

mutable struct LinearAPDDetectorState{T<:AbstractFloat,A<:AbstractVector{T}}
    channels::A
    noise_buffer::A
end

"""
    LinearAPDDetector(; topology=SingleElementAPD(), ...)

Linear-mode single-element APD or fixed channel bank. Inputs are photon fluxes
in photons/channel/second and outputs are one-dimensional channel values. This
is intentionally separate from both area-detector frames and the Geiger-mode
counting `APDDetector`.
"""
struct LinearAPDDetector{N<:NoiseModel,P<:LinearAPDDetectorParams,
    S<:LinearAPDDetectorState,B<:AbstractArrayBackend} <: AbstractDetector
    noise::N
    params::P
    state::S
end

@inline backend(::LinearAPDDetector{<:Any,<:Any,<:Any,B}) where {B} = B()

function LinearAPDDetector(;
    topology::AbstractLinearAPDTopology=SingleElementAPD(),
    integration_time::Real=1.0,
    qe::Real=1.0,
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    dark_current::Real=0.0,
    conversion_gain::Real=1.0,
    noise::NoiseModel=NoisePhoton(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
)
    integration_time > 0 || throw(InvalidConfiguration(
        "LinearAPDDetector integration_time must be > 0"))
    zero(T) <= T(qe) <= one(T) || throw(InvalidConfiguration(
        "LinearAPDDetector qe must lie in [0, 1]"))
    avalanche_gain >= 1 || throw(InvalidConfiguration(
        "LinearAPDDetector avalanche_gain must be >= 1"))
    excess_noise_factor >= 1 || throw(InvalidConfiguration(
        "LinearAPDDetector excess_noise_factor must be >= 1"))
    dark_current >= 0 || throw(InvalidConfiguration(
        "LinearAPDDetector dark_current must be >= 0"))
    conversion_gain > 0 || throw(InvalidConfiguration(
        "LinearAPDDetector conversion_gain must be > 0"))

    converted_noise = validate_noise(convert_noise(noise, T))
    selector = _resolve_backend_selector(backend)
    storage = _resolve_array_backend(backend)
    n_channels = linear_apd_channel_count(topology)
    channels = storage{T}(undef, n_channels)
    noise_buffer = storage{T}(undef, n_channels)
    fill!(channels, zero(T))
    fill!(noise_buffer, zero(T))
    params = LinearAPDDetectorParams{T,typeof(topology)}(
        T(integration_time), T(qe), T(avalanche_gain), T(excess_noise_factor),
        T(dark_current), T(conversion_gain), topology)
    state = LinearAPDDetectorState{T,typeof(channels)}(channels, noise_buffer)
    return LinearAPDDetector{typeof(converted_noise),typeof(params),typeof(state),
        typeof(selector)}(converted_noise, params, state)
end

channel_output(det::LinearAPDDetector) = det.state.channels
readout_ready(::LinearAPDDetector) = true
reset_integration!(det::LinearAPDDetector) = det
supports_detector_thermal_model(::LinearAPDDetector) = false
supports_avalanche_gain(::LinearAPDDetector) = true

@inline linear_apd_photon_noise_enabled(::LinearAPDDetector{NoiseNone}) = false
@inline linear_apd_photon_noise_enabled(::LinearAPDDetector{NoisePhoton}) = true
@inline linear_apd_photon_noise_enabled(
    ::LinearAPDDetector{<:NoiseReadout}) = false
@inline linear_apd_photon_noise_enabled(
    ::LinearAPDDetector{<:NoisePhotonReadout}) = true

linear_apd_readout_sigma(::NoiseNone, ::Type{T}) where {T<:AbstractFloat} =
    zero(T)
linear_apd_readout_sigma(::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} =
    zero(T)
linear_apd_readout_sigma(noise::NoiseReadout,
    ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)
linear_apd_readout_sigma(noise::NoisePhotonReadout,
    ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)

function apply_linear_apd_avalanche!(det::LinearAPDDetector,
    rng::AbstractRNG)
    channels = det.state.channels
    factor = det.params.excess_noise_factor
    if factor > one(factor)
        randn_backend!(rng, det.state.noise_buffer)
        scale2 = factor * factor - one(factor)
        zero_t = zero(eltype(channels))
        @. channels = max(channels + sqrt(max(scale2 * channels, zero_t)) *
            det.state.noise_buffer, zero_t)
    end
    channels .*= det.params.avalanche_gain
    return channels
end

function finalize_linear_apd_capture!(det::LinearAPDDetector,
    rng::AbstractRNG)
    channels = det.state.channels
    T = eltype(channels)
    @. channels = channels * (det.params.qe * det.params.integration_time) +
        det.params.dark_current * det.params.integration_time
    linear_apd_photon_noise_enabled(det) && poisson_noise!(rng, channels)
    apply_linear_apd_avalanche!(det, rng)
    sigma = linear_apd_readout_sigma(det.noise, T)
    if sigma > zero(T)
        randn_backend!(rng, det.state.noise_buffer)
        @. channels += sigma * det.state.noise_buffer
    end
    channels .*= det.params.conversion_gain
    return channels
end

function capture!(det::LinearAPDDetector, photon_flux::AbstractVector;
    rng::AbstractRNG=Random.default_rng())
    length(photon_flux) == length(det.state.channels) ||
        throw(DimensionMismatchError(
            "linear APD input length must match its fixed channel topology"))
    copyto!(det.state.channels, photon_flux)
    return finalize_linear_apd_capture!(det, rng)
end

capture!(det::LinearAPDDetector, photon_flux::AbstractVector,
    rng::AbstractRNG) = capture!(det, photon_flux; rng=rng)

function capture!(det::LinearAPDDetector, photon_flux::Real;
    rng::AbstractRNG=Random.default_rng())
    length(det.state.channels) == 1 || throw(DimensionMismatchError(
        "scalar linear APD capture requires SingleElementAPD topology"))
    fill!(det.state.channels, photon_flux)
    return finalize_linear_apd_capture!(det, rng)
end

capture!(det::LinearAPDDetector, photon_flux::Real,
    rng::AbstractRNG) = capture!(det, photon_flux; rng=rng)

struct LinearAPDExportMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    avalanche_gain::T
    excess_noise_factor::T
    dark_current::T
    conversion_gain::T
    topology::Symbol
    n_channels::Int
    noise::Symbol
end

function detector_export_metadata(det::LinearAPDDetector;
    T::Type{<:AbstractFloat}=eltype(det.state.channels))
    return LinearAPDExportMetadata{T}(
        T(det.params.integration_time), T(det.params.qe),
        T(det.params.avalanche_gain), T(det.params.excess_noise_factor),
        T(det.params.dark_current), T(det.params.conversion_gain),
        linear_apd_topology_symbol(det.params.topology),
        length(det.state.channels), detector_noise_symbol(det.noise))
end
