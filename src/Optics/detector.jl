using Random

abstract type NoiseModel end
struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end
struct NoiseReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end
struct NoisePhotonReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end

NoiseReadout(sigma::Real) = NoiseReadout{Float64}(float(sigma))
NoisePhotonReadout(sigma::Real) = NoisePhotonReadout{Float64}(float(sigma))

struct DetectorParams{T<:AbstractFloat}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    frame::A
    bin_buffer::A
    noise_buffer::A
end

struct Detector{N<:NoiseModel,P<:DetectorParams,S<:DetectorState} <: AbstractDetector
    noise::N
    params::P
    state::S
end

normalize_noise(noise::NoiseModel) = noise

function normalize_noise(noises::Tuple{Vararg{NoiseModel}})
    if isempty(noises)
        return NoiseNone()
    end
    has_photon = false
    has_readout = false
    sigma = 0.0
    sigma_set = false

    for noise in noises
        if noise isa NoiseNone
            continue
        elseif noise isa NoisePhoton
            has_photon = true
        elseif noise isa NoiseReadout
            has_readout = true
            sigma, sigma_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
        elseif noise isa NoisePhotonReadout
            has_photon = true
            has_readout = true
            sigma, sigma_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
        end
    end

    if has_photon && has_readout
        return NoisePhotonReadout(sigma)
    elseif has_photon
        return NoisePhoton()
    elseif has_readout
        return NoiseReadout(sigma)
    end
    return NoiseNone()
end

function merge_readout_sigma(current::Float64, set::Bool, new_sigma::Real)
    σ = float(new_sigma)
    if set && σ != current
        throw(InvalidConfiguration("conflicting readout noise values in noise tuple"))
    end
    return σ, true
end

convert_noise(noise::NoiseNone, ::Type{T}) where {T<:AbstractFloat} = NoiseNone()
convert_noise(noise::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} = NoisePhoton()
convert_noise(noise::NoiseReadout, ::Type{T}) where {T<:AbstractFloat} = NoiseReadout{T}(T(noise.sigma))
convert_noise(noise::NoisePhotonReadout, ::Type{T}) where {T<:AbstractFloat} = NoisePhotonReadout{T}(T(noise.sigma))

validate_noise(noise::NoiseNone) = noise
validate_noise(noise::NoisePhoton) = noise
function validate_noise(noise::NoiseReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoiseReadout"))
    end
    return noise
end
function validate_noise(noise::NoisePhotonReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoisePhotonReadout"))
    end
    return noise
end

function _build_detector(noise::NoiseModel; integration_time::Real, qe::Real,
    psf_sampling::Int, binning::Int, T::Type{<:AbstractFloat}, backend)
    params = DetectorParams{T}(T(integration_time), T(qe), psf_sampling, binning)
    frame = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    fill!(frame, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    state = DetectorState{T, typeof(frame)}(frame, bin_buffer, noise_buffer)
    return Detector{typeof(noise), typeof(params), typeof(state)}(noise, params, state)
end

function Detector(; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T}
    n_in, m_in = size(psf)
    sampling = det.params.psf_sampling
    binning = det.params.binning
    if sampling < 1 || binning < 1
        throw(InvalidConfiguration("psf_sampling and binning must be >= 1"))
    end
    if n_in % sampling != 0 || m_in % sampling != 0
        throw(DimensionMismatchError("psf_sampling must evenly divide input dimensions"))
    end
    n_mid = div(n_in, sampling)
    m_mid = div(m_in, sampling)
    if n_mid % binning != 0 || m_mid % binning != 0
        throw(DimensionMismatchError("binning must evenly divide sampled dimensions"))
    end
    n_out = div(n_mid, binning)
    m_out = div(m_mid, binning)
    ensure_buffers!(det, n_mid, m_mid, n_out, m_out)

    if sampling > 1
        bin2d!(det.state.bin_buffer, psf, sampling)
        if binning > 1
            bin2d!(det.state.frame, det.state.bin_buffer, binning)
        else
            det.state.frame .= det.state.bin_buffer
        end
    else
        if binning > 1
            bin2d!(det.state.frame, psf, binning)
        else
            det.state.frame .= psf
        end
    end
    @. det.state.frame *= det.params.qe * det.params.integration_time
    return det.state.frame
end

function capture!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    fill_frame!(det, psf)
    return det.state.frame
end

function capture!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    fill_frame!(det, psf)
    poisson_noise!(rng, det.state.frame)
    return det.state.frame
end

function capture!(det::Detector{NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    fill_frame!(det, psf)
    randn!(rng, det.state.noise_buffer)
    det.state.frame .+= det.noise.sigma .* det.state.noise_buffer
    return det.state.frame
end

function capture!(det::Detector{NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    fill_frame!(det, psf)
    poisson_noise!(rng, det.state.frame)
    randn!(rng, det.state.noise_buffer)
    det.state.frame .+= det.noise.sigma .* det.state.noise_buffer
    return det.state.frame
end

function capture!(det::Detector, psf::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T}
    return capture!(det, psf, rng)
end

function ensure_buffers!(det::Detector, n_mid::Int, m_mid::Int, n_out::Int, m_out::Int)
    if size(det.state.frame) != (n_out, m_out)
        det.state.frame = similar(det.state.frame, n_out, m_out)
    end
    if size(det.state.bin_buffer) != (n_mid, m_mid)
        det.state.bin_buffer = similar(det.state.bin_buffer, n_mid, m_mid)
    end
    if size(det.state.noise_buffer) != (n_out, m_out)
        det.state.noise_buffer = similar(det.state.noise_buffer, n_out, m_out)
    end
    return det
end
