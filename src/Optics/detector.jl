using Random

struct DetectorParams{T<:AbstractFloat}
    integration_time::T
    readout_noise::T
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
    for noise in noises
        has_photon |= noise isa NoisePhoton || noise isa NoisePhotonReadout
        has_readout |= noise isa NoiseReadout || noise isa NoisePhotonReadout
    end
    if has_photon && has_readout
        return NoisePhotonReadout()
    elseif has_photon
        return NoisePhoton()
    elseif has_readout
        return NoiseReadout()
    end
    return NoiseNone()
end

function _build_detector(noise::NoiseModel; integration_time::Real, readout_noise::Real, qe::Real,
    psf_sampling::Int, binning::Int, T::Type{<:AbstractFloat}, backend)
    params = DetectorParams{T}(T(integration_time), T(readout_noise), T(qe), psf_sampling, binning)
    frame = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    fill!(frame, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    state = DetectorState{T, typeof(frame)}(frame, bin_buffer, noise_buffer)
    return Detector{typeof(noise), typeof(params), typeof(state)}(params, state)
end

function Detector(; integration_time::Real=1.0, readout_noise::Real=0.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, readout_noise=readout_noise, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function Detector(::NoiseNone; integration_time::Real=1.0, readout_noise::Real=0.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    if readout_noise != 0
        throw(InvalidConfiguration("readout_noise must be 0 for NoiseNone"))
    end
    return _build_detector(NoiseNone(); integration_time=integration_time, readout_noise=readout_noise, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function Detector(::NoisePhoton; integration_time::Real=1.0, readout_noise::Real=0.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    if readout_noise != 0
        throw(InvalidConfiguration("readout_noise must be 0 for NoisePhoton"))
    end
    return _build_detector(NoisePhoton(); integration_time=integration_time, readout_noise=readout_noise, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function Detector(::NoiseReadout; integration_time::Real=1.0, readout_noise::Real=0.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    if readout_noise <= 0
        throw(InvalidConfiguration("readout_noise must be > 0 for NoiseReadout"))
    end
    return _build_detector(NoiseReadout(); integration_time=integration_time, readout_noise=readout_noise, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function Detector(::NoisePhotonReadout; integration_time::Real=1.0, readout_noise::Real=0.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    if readout_noise <= 0
        throw(InvalidConfiguration("readout_noise must be > 0 for NoisePhotonReadout"))
    end
    return _build_detector(NoisePhotonReadout(); integration_time=integration_time, readout_noise=readout_noise, qe=qe,
        psf_sampling=psf_sampling, binning=binning, T=T, backend=backend)
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T}
    n_in, m_in = size(psf)
    if det.params.binning > 1
        n_out = div(n_in, det.params.binning)
        m_out = div(m_in, det.params.binning)
        ensure_buffers!(det, n_in, m_in, n_out, m_out)
        @. det.state.bin_buffer = psf * det.params.qe * det.params.integration_time
        bin2d!(det.state.frame, det.state.bin_buffer, det.params.binning)
    else
        ensure_buffers!(det, n_in, m_in, n_in, m_in)
        @. det.state.frame = psf * det.params.qe * det.params.integration_time
    end
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
    det.state.frame .+= det.params.readout_noise .* det.state.noise_buffer
    return det.state.frame
end

function capture!(det::Detector{NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    fill_frame!(det, psf)
    poisson_noise!(rng, det.state.frame)
    randn!(rng, det.state.noise_buffer)
    det.state.frame .+= det.params.readout_noise .* det.state.noise_buffer
    return det.state.frame
end

function capture!(det::Detector, psf::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T}
    return capture!(det, psf, rng)
end

function ensure_buffers!(det::Detector, n_in::Int, m_in::Int, n_out::Int, m_out::Int)
    if size(det.state.frame) != (n_out, m_out)
        det.state.frame = similar(det.state.frame, n_out, m_out)
    end
    if size(det.state.bin_buffer) != (n_in, m_in)
        det.state.bin_buffer = similar(det.state.bin_buffer, n_in, m_in)
    end
    if size(det.state.noise_buffer) != (n_out, m_out)
        det.state.noise_buffer = similar(det.state.noise_buffer, n_out, m_out)
    end
    return det
end
