using Random

@kernel function kolmogorov_psd_kernel!(psd, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds begin
            fx = freqs[i]
            fy = freqs[j]
            psd[i, j] = ifelse(
                i == 1 && j == 1,
                zero(eltype(psd)),
                coeff * (two_pi_sq * (fx * fx + fy * fy) + inv_L0_sq)^exponent,
            )
        end
    end
end

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
    noise_re_host::Matrix{T}
    noise_im_host::Matrix{T}
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
    noise_re_host = Matrix{T}(undef, n, n)
    noise_im_host = Matrix{T}(undef, n, n)
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
        noise_re_host,
        noise_im_host,
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
    coeff = T(0.023) * atm.params.r0^(-T(5) / T(3))
    inv_L0_sq = (T(1) / atm.params.L0)^2
    two_pi_sq = T(2 * pi)^2
    exponent = -T(11) / T(6)
    update_psd!(execution_style(atm.state.psd), atm.state.psd, atm.state.freqs, coeff, two_pi_sq, inv_L0_sq, exponent, n)
    return atm
end

function update_psd!(::ScalarCPUStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T, n::Int) where {T<:AbstractFloat}
    @inbounds for j in 1:n, i in 1:n
        fx = freqs[i]
        fy = freqs[j]
        psd[i, j] = i == 1 && j == 1 ? zero(T) :
            coeff * (two_pi_sq * (fx * fx + fy * fy) + inv_L0_sq)^exponent
    end
    return psd
end

function update_psd!(style::AcceleratorStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T, n::Int) where {T<:AbstractFloat}
    launch_kernel!(style, kolmogorov_psd_kernel!, psd, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, n;
        ndrange=size(psd))
    return psd
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
    atm.state.noise_re_host = randn_phase_noise!(rng, atm.state.noise_re, atm.state.noise_re_host)
    atm.state.noise_im_host = randn_phase_noise!(rng, atm.state.noise_im, atm.state.noise_im_host)
    @. atm.state.spectrum = complex(atm.state.noise_re, atm.state.noise_im) * sqrt(atm.state.psd)
    execute_fft_plan!(atm.state.spectrum, atm.state.ifft_plan)
    @. out = real(atm.state.spectrum) * (n * delta)
    return out
end

function randn_phase_noise!(rng::AbstractRNG, out::AbstractMatrix{T}, host::Matrix{T}) where {T<:AbstractFloat}
    randn_backend!(rng, out)
    return host
end


function propagate!(atm::KolmogorovAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
