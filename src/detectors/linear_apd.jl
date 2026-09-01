abstract type AbstractLinearAPDTopology end

"""
    SingleElementLinearAPD()

Select one linear-mode avalanche-photodiode element. The detector accepts a
scalar or one-element-vector photon flux and retains one-dimensional storage.
"""
struct SingleElementLinearAPD <: AbstractLinearAPDTopology end

"""
    LinearAPDChannelBank(n_channels)

Select a fixed bank of `n_channels` independent linear-mode avalanche-
photodiode channels. Use [`SingleElementLinearAPD`](@ref) for one channel.
"""
struct LinearAPDChannelBank <: AbstractLinearAPDTopology
    n_channels::Int
    function LinearAPDChannelBank(n_channels::Integer)
        n_channels >= 2 || throw(InvalidConfiguration(
            "LinearAPDChannelBank n_channels must be >= 2; use SingleElementLinearAPD for one channel"))
        return new(Int(n_channels))
    end
end

linear_apd_channel_count(::SingleElementLinearAPD) = 1
linear_apd_channel_count(topology::LinearAPDChannelBank) = topology.n_channels
linear_apd_topology_symbol(::SingleElementLinearAPD) = :single_element
linear_apd_topology_symbol(::LinearAPDChannelBank) = :channel_bank

struct LinearAPDDetectorParams{T<:AbstractFloat,TP<:AbstractLinearAPDTopology}
    exposure_duration::T
    qe::T
    avalanche_gain::T
    excess_noise_factor::T
    dark_current::T
    conversion_gain::T
    topology::TP
end

"""Replaceable stochastic scratch for linear-mode APD acquisition."""
struct LinearAPDDetectorWorkspace{T<:AbstractFloat,A<:AbstractVector{T}}
    noise_buffer::A
end

"""Caller-visible analog channel values from linear-mode APD acquisition."""
struct LinearAPDDetectorProducts{T<:AbstractFloat,A<:AbstractVector{T}}
    channels::A
end

"""
    LinearAPDDetector(;
        topology=SingleElementLinearAPD(),
        exposure_duration=1.0,
        qe=1.0,
        avalanche_gain=1.0,
        excess_noise_factor=1.0,
        dark_current=0.0,
        conversion_gain=1.0,
        noise=NoisePhoton(),
        T=Float64,
        backend=CPUBackend(),
    )

Construct a linear-mode single-element APD or fixed channel bank. Inputs are
photon fluxes in photons per channel per second and outputs are
one-dimensional analog channel values. The signal path applies quantum
efficiency and integration, dark current, optional photon noise, avalanche
multiplication, additive read noise, and conversion gain in that order.

`excess_noise_factor` uses the detector convention
`F = E[M^2] / E[M]^2`. This API is intentionally separate from both
area-detector frames and Geiger-mode SPAD counting arrays.
"""
struct LinearAPDDetector{N<:NoiseModel,P<:LinearAPDDetectorParams,
    W<:LinearAPDDetectorWorkspace,R<:LinearAPDDetectorProducts,
    B<:AbstractArrayBackend} <: AbstractDetector
    noise::N
    params::P
    workspace::W
    products::R
end

@inline backend(
    ::LinearAPDDetector{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

function LinearAPDDetector(;
    topology::AbstractLinearAPDTopology=SingleElementLinearAPD(),
    exposure_duration::Real=1.0,
    qe::Real=1.0,
    avalanche_gain::Real=1.0,
    excess_noise_factor::Real=1.0,
    dark_current::Real=0.0,
    conversion_gain::Real=1.0,
    noise::NoiseModel=NoisePhoton(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
)
    exposure_duration_t = T(exposure_duration)
    qe_t = T(qe)
    avalanche_gain_t = T(avalanche_gain)
    excess_noise_factor_t = T(excess_noise_factor)
    dark_current_t = T(dark_current)
    conversion_gain_t = T(conversion_gain)
    isfinite(exposure_duration_t) && exposure_duration_t > zero(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector exposure_duration must be finite and > 0"))
    isfinite(qe_t) && zero(T) <= qe_t <= one(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector qe must be finite and lie in [0, 1]"))
    isfinite(avalanche_gain_t) && avalanche_gain_t >= one(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector avalanche_gain must be finite and >= 1"))
    isfinite(excess_noise_factor_t) && excess_noise_factor_t >= one(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector excess_noise_factor must be finite and >= 1"))
    isfinite(dark_current_t) && dark_current_t >= zero(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector dark_current must be finite and >= 0"))
    isfinite(conversion_gain_t) && conversion_gain_t > zero(T) ||
        throw(InvalidConfiguration(
            "LinearAPDDetector conversion_gain must be finite and > 0"))

    converted_noise = validate_noise(convert_noise(noise, T))
    selector = _resolve_backend_selector(backend)
    storage = _resolve_array_backend(backend)
    n_channels = linear_apd_channel_count(topology)
    channels = storage{T}(undef, n_channels)
    noise_buffer = storage{T}(undef, n_channels)
    fill!(channels, zero(T))
    fill!(noise_buffer, zero(T))
    params = LinearAPDDetectorParams{T,typeof(topology)}(
        exposure_duration_t, qe_t, avalanche_gain_t, excess_noise_factor_t,
        dark_current_t, conversion_gain_t, topology)
    workspace = LinearAPDDetectorWorkspace{T,typeof(noise_buffer)}(
        noise_buffer)
    products = LinearAPDDetectorProducts{T,typeof(channels)}(channels)
    return LinearAPDDetector{typeof(converted_noise),typeof(params),
        typeof(workspace),typeof(products),typeof(selector)}(
        converted_noise, params, workspace, products)
end

channel_output(det::LinearAPDDetector) = det.products.channels
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
    channels = det.products.channels
    factor = det.params.excess_noise_factor
    if factor > one(factor)
        randn_backend!(rng, det.workspace.noise_buffer)
        scale2 = factor - one(factor)
        zero_t = zero(eltype(channels))
        @. channels = max(channels + sqrt(max(scale2 * channels, zero_t)) *
            det.workspace.noise_buffer, zero_t)
    end
    channels .*= det.params.avalanche_gain
    return channels
end

function finalize_linear_apd_capture!(det::LinearAPDDetector,
    rng::AbstractRNG)
    channels = det.products.channels
    T = eltype(channels)
    @. channels = channels * (det.params.qe * det.params.exposure_duration) +
        det.params.dark_current * det.params.exposure_duration
    linear_apd_photon_noise_enabled(det) && poisson_noise!(rng, channels)
    apply_linear_apd_avalanche!(det, rng)
    sigma = linear_apd_readout_sigma(det.noise, T)
    if sigma > zero(T)
        randn_backend!(rng, det.workspace.noise_buffer)
        @. channels += sigma * det.workspace.noise_buffer
    end
    channels .*= det.params.conversion_gain
    return channels
end

function capture!(det::LinearAPDDetector, photon_flux::AbstractVector;
    rng::AbstractRNG=runtime_rng())
    channels = det.products.channels
    length(photon_flux) == length(channels) ||
        throw(DimensionMismatchError(
            "linear APD input length must match its fixed channel topology"))
    photon_flux === channels || copyto!(channels, photon_flux)
    return finalize_linear_apd_capture!(det, rng)
end

capture!(det::LinearAPDDetector, photon_flux::AbstractVector,
    rng::AbstractRNG) = capture!(det, photon_flux; rng=rng)

function capture!(det::LinearAPDDetector, photon_flux::Real;
    rng::AbstractRNG=runtime_rng())
    length(det.products.channels) == 1 || throw(DimensionMismatchError(
        "scalar linear APD capture requires SingleElementLinearAPD topology"))
    fill!(det.products.channels, photon_flux)
    return finalize_linear_apd_capture!(det, rng)
end

capture!(det::LinearAPDDetector, photon_flux::Real,
    rng::AbstractRNG) = capture!(det, photon_flux; rng=rng)

struct LinearAPDExportMetadata{T<:AbstractFloat}
    exposure_duration::T
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
    T::Type{<:AbstractFloat}=eltype(det.products.channels))
    return LinearAPDExportMetadata{T}(
        T(det.params.exposure_duration), T(det.params.qe),
        T(det.params.avalanche_gain), T(det.params.excess_noise_factor),
        T(det.params.dark_current), T(det.params.conversion_gain),
        linear_apd_topology_symbol(det.params.topology),
        length(det.products.channels), detector_noise_symbol(det.noise))
end
