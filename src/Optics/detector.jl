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
end

struct Detector{P<:DetectorParams,S<:DetectorState} <: AbstractDetector
    params::P
    state::S
end

function Detector(; integration_time::Real=1.0, photon_noise::Bool=true, readout_noise::Real=0.0,
    qe::Real=1.0, psf_sampling::Int=1, binning::Int=1, T::Type{<:AbstractFloat}=Float64, backend=Array)

    params = DetectorParams{T}(T(integration_time), photon_noise, T(readout_noise), T(qe), psf_sampling, binning)
    frame = backend{T}(undef, 1, 1)
    fill!(frame, zero(T))
    state = DetectorState{T, typeof(frame)}(frame)
    return Detector(params, state)
end

function ensure_frame!(det::Detector, n::Int, m::Int=n)
    if size(det.state.frame) != (n, m)
        det.state.frame = similar(det.state.frame, n, m)
    end
    return det
end

function capture!(det::Detector, psf::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T}
    img = psf .* det.params.qe .* det.params.integration_time
    if det.params.binning > 1
        img = bin2d(img, det.params.binning)
    end

    ensure_frame!(det, size(img, 1), size(img, 2))
    det.state.frame .= img

    if det.params.photon_noise
        poisson_noise!(rng, det.state.frame)
    end

    if det.params.readout_noise > 0
        det.state.frame .+= det.params.readout_noise .* randn(rng, eltype(det.state.frame), size(det.state.frame))
    end

    return det.state.frame
end
