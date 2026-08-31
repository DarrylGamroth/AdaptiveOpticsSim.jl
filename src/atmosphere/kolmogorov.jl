@kernel function kolmogorov_psd_kernel!(psd, freqs, coeff, inv_L0_sq, exponent, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds begin
            fx = freqs[i]
            fy = freqs[j]
            psd[i, j] = ifelse(
                i == 1 && j == 1,
                zero(eltype(psd)),
                coeff * (fx * fx + fy * fy + inv_L0_sq)^exponent,
            )
        end
    end
end

struct KolmogorovParams{T<:AbstractFloat}
    r0::T
    L0::T
    reference_wavelength_m::T
    opd_per_radian::T
    sampling_m::T
end

mutable struct KolmogorovState{T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractMatrix{Complex{T}},V<:AbstractVector{T},P}
    phase_rad::A
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
    timeline::AtmosphereTimelineState{T}
end

struct KolmogorovAtmosphere{P<:KolmogorovParams,S<:KolmogorovState,I<:AtmosphereIdentity,B<:AbstractArrayBackend} <: AbstractTimedAtmosphere
    params::P
    state::S
    identity::I
end

@inline backend(::KolmogorovAtmosphere{<:Any,<:Any,<:Any,B}) where {B} = B()
@inline atmosphere_layers(::KolmogorovAtmosphere) = ()
@inline atmosphere_numeric_type(atm::KolmogorovAtmosphere) =
    typeof(atm.params.r0)

"""
    KolmogorovAtmosphere(telescope; r0, reference_wavelength_m, L0=25, ...)

Prepare one finite von Kármán phase screen. `r0` is the Fried parameter in
metres at `reference_wavelength_m`. The persistent screen is phase in radians
at that reference wavelength. Atmosphere rendering converts it to OPD in
metres with `reference_wavelength_m / (2π)` before it enters an optical path.
"""
function KolmogorovAtmosphere(tel::Telescope;
    r0::Real,
    reference_wavelength_m::Real,
    L0::Real=25.0,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel),
)
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    r0_t = _converted_positive_finite(r0, T,
        "atmosphere Fried parameter r0")
    L0_t = _converted_positive_finite(L0, T,
        "atmosphere outer scale L0")
    reference_wavelength_t = _converted_positive_finite(
        reference_wavelength_m,
        T,
        "atmosphere reference wavelength in metres",
    )
    params = KolmogorovParams{T}(
        r0_t,
        L0_t,
        reference_wavelength_t,
        reference_wavelength_t / T(2 * pi),
        T(tel.aperture.sampling_m[1]),
    )
    n = tel.params.resolution
    phase_rad = backend{T}(undef, n, n)
    psd = backend{T}(undef, n, n)
    spectrum = backend{Complex{T}}(undef, n, n)
    noise_re = backend{T}(undef, n, n)
    noise_im = backend{T}(undef, n, n)
    noise_re_host = Matrix{T}(undef, n, n)
    noise_im_host = Matrix{T}(undef, n, n)
    freqs = backend{T}(undef, n)
    ifft_plan = plan_ifft_backend!(spectrum)
    fill!(phase_rad, zero(T))
    fill!(psd, zero(T))
    fill!(spectrum, zero(eltype(spectrum)))
    state = KolmogorovState{T, typeof(phase_rad), typeof(spectrum),
        typeof(freqs), typeof(ifft_plan)}(
        phase_rad,
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
        new_atmosphere_timeline(T),
    )
    identity = AtmosphereIdentity()
    return KolmogorovAtmosphere{typeof(params), typeof(state), typeof(identity),
        typeof(selector)}(params, state, identity)
end

function update_psd!(atm::KolmogorovAtmosphere, delta::Real)
    n = size(atm.state.phase_rad, 1)
    T = eltype(atm.state.phase_rad)
    fftfreq!(atm.state.freqs, n; d=delta)
    coeff = T(0.023) * atm.params.r0^(-T(5) / T(3))
    inv_L0_sq = (T(1) / atm.params.L0)^2
    exponent = -T(11) / T(6)
    update_psd!(execution_style(atm.state.psd), atm.state.psd,
        atm.state.freqs, coeff, inv_L0_sq, exponent, n)
    return atm
end

function update_psd!(::ScalarCPUStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, inv_L0_sq::T, exponent::T, n::Int) where {T<:AbstractFloat}
    @inbounds for j in 1:n, i in 1:n
        fx = freqs[i]
        fy = freqs[j]
        psd[i, j] = i == 1 && j == 1 ? zero(T) :
            coeff * (fx * fx + fy * fy + inv_L0_sq)^exponent
    end
    return psd
end

function update_psd!(style::AcceleratorStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, inv_L0_sq::T, exponent::T, n::Int) where {T<:AbstractFloat}
    launch_kernel!(style, kolmogorov_psd_kernel!, psd, freqs, coeff,
        inv_L0_sq, exponent, n;
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

function regenerate_phase_screen!(atm::KolmogorovAtmosphere,
    rng::AbstractRNG)
    delta = atm.params.sampling_m
    ensure_psd!(atm, delta)
    phase_screen_von_karman!(atm.state.phase_rad, atm, delta, rng)
    return atm
end

@inline initialize_atmosphere!(atm::KolmogorovAtmosphere, rng::AbstractRNG) =
    regenerate_phase_screen!(atm, rng)

# Initial publication already generated the independent screen. Subsequent
# positive-time publications refresh it once, preserving legacy semantics.
@inline evolve_initial_atmosphere!(atm::KolmogorovAtmosphere, ::Real,
    ::AbstractRNG) = atm

@inline evolve_atmosphere!(atm::KolmogorovAtmosphere, ::Real,
    rng::AbstractRNG) = regenerate_phase_screen!(atm, rng)

function phase_screen_von_karman!(out::AbstractMatrix, atm::KolmogorovAtmosphere, delta::Real, rng::AbstractRNG)
    n = size(out, 1)
    if size(out) != (n, n)
        throw(DimensionMismatchError("output must be square"))
    end
    atm.state.noise_re_host = randn_phase_noise!(rng, atm.state.noise_re, atm.state.noise_re_host)
    atm.state.noise_im_host = randn_phase_noise!(rng, atm.state.noise_im, atm.state.noise_im_host)
    @. atm.state.spectrum = complex(atm.state.noise_re, atm.state.noise_im) * sqrt(atm.state.psd)
    execute_fft_plan!(atm.state.spectrum, atm.state.ifft_plan)
    phase_scale = eltype(out)(n) / eltype(out)(delta)
    @. out = real(atm.state.spectrum) * phase_scale
    return out
end

function randn_phase_noise!(rng::AbstractRNG, out::AbstractMatrix{T}, host::Matrix{T}) where {T<:AbstractFloat}
    randn_backend!(rng, out)
    return host
end
function render_atmosphere_opd_impl!(dest::AbstractMatrix,
    renderer::AtmosphereDirectionRenderer, atm::KolmogorovAtmosphere)
    size(dest) == size(atm.state.phase_rad) || throw(DimensionMismatchError(
        "Kolmogorov screen dimensions do not match prepared output"))
    @. dest = atm.state.phase_rad * renderer.pupil *
        atm.params.opd_per_radian
    return dest
end
