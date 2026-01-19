using Random

struct DetectorParams{T<:AbstractFloat}
    integration_time::T
    photon_noise::Bool
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

struct Detector{P<:DetectorParams,S<:DetectorState} <: AbstractDetector
    params::P
    state::S
end

function Detector(; integration_time::Real=1.0, photon_noise::Bool=true, readout_noise::Real=0.0,
    qe::Real=1.0, psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)

    params = DetectorParams{T}(T(integration_time), photon_noise, T(readout_noise), T(qe), psf_sampling, binning)
    frame = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    fill!(frame, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    state = DetectorState{T, typeof(frame)}(frame, bin_buffer, noise_buffer)
    return Detector(params, state)
end

function capture!(det::Detector, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
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

    if det.params.photon_noise
        poisson_noise!(rng, det.state.frame)
    end

    if det.params.readout_noise > 0
        randn!(rng, det.state.noise_buffer)
        det.state.frame .+= det.params.readout_noise .* det.state.noise_buffer
    end

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
