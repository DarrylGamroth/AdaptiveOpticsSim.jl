using Random

struct KolmogorovParams{T<:AbstractFloat}
    r0::T
    L0::T
end

mutable struct KolmogorovState{T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractMatrix{Complex{T}},V<:AbstractVector{T},P}
    opd::A
    psd::A
    spectrum::B
    noise_re::A
    noise_im::A
    freqs::V
    ifft_plan::P
    last_delta::T
    last_r0::T
    last_L0::T
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
    ifft_plan = plan_ifft_backend!(spectrum)
    fill!(opd, zero(T))
    fill!(psd, zero(T))
    fill!(spectrum, zero(eltype(spectrum)))
    state = KolmogorovState{T, typeof(opd), typeof(spectrum), typeof(freqs), typeof(ifft_plan)}(
        opd,
        psd,
        spectrum,
        noise_re,
        noise_im,
        freqs,
        ifft_plan,
        T(-1),
        T(-1),
        T(-1),
    )
    return KolmogorovAtmosphere(params, state)
end

function update_psd!(atm::KolmogorovAtmosphere, delta::Real)
    n = size(atm.state.opd, 1)
    T = eltype(atm.state.opd)
    fftfreq!(atm.state.freqs, n; d=delta)
    freqs = atm.state.freqs

    r0 = atm.params.r0
    L0 = atm.params.L0
    freq_x = reshape(freqs, n, 1)
    freq_y = reshape(freqs, 1, n)
    coeff = T(0.023) * r0^(-T(5) / T(3))
    inv_L0 = T(1) / L0
    two_pi_sq = T(2 * pi)^2
    @. atm.state.psd = coeff * (two_pi_sq * (freq_x^2 + freq_y^2) + inv_L0^2)^(-T(11) / T(6))
    return atm
end

function ensure_psd!(atm::KolmogorovAtmosphere, delta::Real)
    if atm.state.last_delta != delta || atm.state.last_r0 != atm.params.r0 || atm.state.last_L0 != atm.params.L0
        update_psd!(atm, delta)
        atm.state.last_delta = delta
        atm.state.last_r0 = atm.params.r0
        atm.state.last_L0 = atm.params.L0
    end
    return atm
end

function advance!(atm::KolmogorovAtmosphere, tel::Telescope, rng::AbstractRNG)
    n = tel.params.resolution
    delta = tel.params.diameter / n
    ensure_psd!(atm, delta)
    phase_screen_von_karman!(atm.state.opd, atm, delta, rng)
    return atm
end

function advance!(atm::KolmogorovAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    return advance!(atm, tel, rng)
end

function phase_screen_von_karman!(out::AbstractMatrix, atm::KolmogorovAtmosphere, delta::Real, rng::AbstractRNG)
    n = size(out, 1)
    if size(out) != (n, n)
        throw(DimensionMismatchError("output must be square"))
    end
    randn!(rng, atm.state.noise_re)
    randn!(rng, atm.state.noise_im)
    @. atm.state.spectrum = complex(atm.state.noise_re, atm.state.noise_im) * sqrt(atm.state.psd)
    mul!(atm.state.spectrum, atm.state.ifft_plan, atm.state.spectrum)
    @. out = real(atm.state.spectrum) * (n * delta)
    return out
end


function propagate!(atm::KolmogorovAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
