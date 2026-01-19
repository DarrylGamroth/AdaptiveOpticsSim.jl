using FFTW
using Random

struct KolmogorovParams{T<:AbstractFloat}
    r0::T
    L0::T
end

mutable struct KolmogorovState{T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractMatrix{Complex{T}},V<:AbstractVector{T}}
    opd::A
    psd::A
    spectrum::B
    noise_re::A
    noise_im::A
    freqs::V
end

struct KolmogorovAtmosphere{P<:KolmogorovParams,S<:KolmogorovState} <: AbstractAtmosphere
    params::P
    state::S
end

function KolmogorovAtmosphere(tel::Telescope; r0::Real, L0::Real=25.0, T::Type{<:AbstractFloat}=Float64, backend=Array)
    params = KolmogorovParams{T}(T(r0), T(L0))
    n = tel.params.resolution
    opd = backend{T}(undef, n, n)
    psd = backend{T}(undef, n, n)
    spectrum = backend{Complex{T}}(undef, n, n)
    noise_re = backend{T}(undef, n, n)
    noise_im = backend{T}(undef, n, n)
    freqs = backend{T}(undef, n)
    fill!(opd, zero(T))
    fill!(psd, zero(T))
    fill!(spectrum, zero(eltype(spectrum)))
    state = KolmogorovState{T, typeof(opd), typeof(spectrum), typeof(freqs)}(opd, psd, spectrum, noise_re, noise_im, freqs)
    return KolmogorovAtmosphere(params, state)
end

function update_psd!(atm::KolmogorovAtmosphere, delta::Real)
    n = size(atm.state.opd, 1)
    T = eltype(atm.state.opd)
    freqs = FFTW.fftfreq(n, delta)
    atm.state.freqs .= freqs

    r0 = atm.params.r0
    L0 = atm.params.L0
    @inbounds for i in 1:n, j in 1:n
        f = sqrt(freqs[i]^2 + freqs[j]^2)
        k = 2 * pi * f
        atm.state.psd[i, j] = T(0.023) * r0^(-T(5) / T(3)) * (k^2 + (T(1) / L0)^2)^(-T(11) / T(6))
    end
    return atm
end

function advance!(atm::KolmogorovAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng(), reuse_psd::Bool=true)
    n = tel.params.resolution
    delta = tel.params.diameter / n
    if reuse_psd
        update_psd!(atm, delta)
    end
    phase_screen_von_karman!(atm.state.opd, atm, delta, rng)
    return atm
end

function phase_screen_von_karman!(out::AbstractMatrix, atm::KolmogorovAtmosphere, delta::Real, rng::AbstractRNG)
    n = size(out, 1)
    if size(out) != (n, n)
        throw(DimensionMismatchError("output must be square"))
    end
    randn!(rng, atm.state.noise_re)
    randn!(rng, atm.state.noise_im)
    @. atm.state.spectrum = complex(atm.state.noise_re, atm.state.noise_im) * sqrt(atm.state.psd)
    FFTW.ifft!(atm.state.spectrum)
    @inbounds for i in 1:n, j in 1:n
        out[i, j] = real(atm.state.spectrum[i, j]) * (n * delta)
    end
    return out
end


function propagate!(atm::KolmogorovAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
