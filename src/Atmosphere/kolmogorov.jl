using FFTW
using Random

struct KolmogorovParams{T<:AbstractFloat}
    r0::T
    L0::T
end

mutable struct KolmogorovState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
end

struct KolmogorovAtmosphere{P<:KolmogorovParams,S<:KolmogorovState} <: AbstractAtmosphere
    params::P
    state::S
end

function KolmogorovAtmosphere(tel::Telescope; r0::Real, L0::Real=25.0, T=Float64, backend=Array)
    params = KolmogorovParams{T}(T(r0), T(L0))
    opd = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    fill!(opd, zero(T))
    state = KolmogorovState{T, typeof(opd)}(opd)
    return KolmogorovAtmosphere(params, state)
end

function phase_screen_von_karman(r0::Real, L0::Real, n::Int, delta::Real, rng::AbstractRNG, T)
    fx = FFTW.fftfreq(n, delta)
    fy = fx
    psd = Array{T}(undef, n, n)

    @inbounds for i in 1:n, j in 1:n
        f = sqrt(fx[i]^2 + fy[j]^2)
        k = 2 * pi * f
        psd[i, j] = 0.023 * r0^(-5 / 3) * (k^2 + (1 / L0)^2)^(-11 / 6)
    end

    noise_re = randn(rng, T, n, n)
    noise_im = randn(rng, T, n, n)
    noise = complex.(noise_re, noise_im)

    spectrum = noise .* sqrt.(psd)
    screen = real.(ifft(spectrum))
    return screen .* (n * delta)
end

function advance!(atm::KolmogorovAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    n = tel.params.resolution
    delta = tel.params.diameter / n
    screen = phase_screen_von_karman(atm.params.r0, atm.params.L0, n, delta, rng, eltype(atm.state.opd))
    atm.state.opd .= screen
    return atm
end

function propagate!(atm::KolmogorovAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
