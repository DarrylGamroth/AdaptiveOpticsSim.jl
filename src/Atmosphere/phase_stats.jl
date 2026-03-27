using Random
using SpecialFunctions

@kernel function phase_screen_psd_kernel!(psd, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, inv_fm_sq, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(psd)
        fx = freqs[i]
        fy = freqs[j]
        f_sq = fx * fx + fy * fy
        base = coeff * (two_pi_sq * f_sq + inv_L0_sq)^exponent
        with_inner = inv_fm_sq == zero(T) ? base : base * exp(-f_sq * inv_fm_sq)
        @inbounds psd[i, j] = ifelse(i == 1 && j == 1, zero(T), with_inner)
    end
end

@kernel function phase_spectrum_kernel!(out, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds out[i] = coeff * (two_pi_sq * (freqs[i] * freqs[i]) + inv_L0_sq)^exponent
    end
end

@kernel function add_subharmonics_kernel!(phs, coeff_re, coeff_im, freq_x, freq_y, two_pi, delta_t, offset::Int, n_terms::Int, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(phs)
        xi = (i - 1 - offset) * delta_t
        yj = (j - 1 - offset) * delta_t
        acc = zero(T)
        @inbounds for k in 1:n_terms
            phase = two_pi * (freq_x[k] * xi + freq_y[k] * yj)
            acc += coeff_re[k] * cos(phase) - coeff_im[k] * sin(phase)
        end
        @inbounds phs[i, j] += acc
    end
end

abstract type SubharmonicMode end
struct FastSubharmonics <: SubharmonicMode end
struct FidelitySubharmonics <: SubharmonicMode end

mutable struct PhaseStatsWorkspace{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    V<:AbstractVector{T},
    P}
    spectrum::C
    buffer::C
    psd::R
    noise_re::R
    noise_im::R
    freqs::V
    ifft_plan::P
end

function PhaseStatsWorkspace(n::Int; T::Type{<:AbstractFloat}=Float64, backend=Array)
    spectrum = backend{Complex{T}}(undef, n, n)
    buffer = similar(spectrum)
    psd = backend{T}(undef, n, n)
    noise_re = backend{T}(undef, n, n)
    noise_im = backend{T}(undef, n, n)
    freqs = backend{T}(undef, n)
    ifft_plan = plan_ifft_backend!(spectrum)
    return PhaseStatsWorkspace{T, typeof(spectrum), typeof(psd), typeof(freqs), typeof(ifft_plan)}(
        spectrum,
        buffer,
        psd,
        noise_re,
        noise_im,
        freqs,
        ifft_plan,
    )
end

function PhaseStatsWorkspace(ref::AbstractMatrix, n::Int; T::Type{<:AbstractFloat}=eltype(ref))
    spectrum = similar(ref, Complex{T}, n, n)
    buffer = similar(spectrum)
    psd = similar(ref, T, n, n)
    noise_re = similar(ref, T, n, n)
    noise_im = similar(ref, T, n, n)
    freqs = similar(ref, T, n)
    ifft_plan = plan_ifft_backend!(spectrum)
    return PhaseStatsWorkspace{T, typeof(spectrum), typeof(psd), typeof(freqs), typeof(ifft_plan)}(
        spectrum,
        buffer,
        psd,
        noise_re,
        noise_im,
        freqs,
        ifft_plan,
    )
end

function phase_variance(r0::Real, L0::Real; cn2::Real=1.0)
    L0r0 = (L0 / r0)^(5 / 3)
    cst = (24 * gamma(6 / 5))^(5 / 6) * (gamma(11 / 6) * gamma(5 / 6)) / (2 * pi^(8 / 3))
    return cn2 * cst * L0r0
end

phase_variance(atm::KolmogorovAtmosphere; cn2::Real=1.0) = phase_variance(atm.params.r0, atm.params.L0; cn2=cn2)
phase_variance(atm::MultiLayerAtmosphere; cn2::Real=sum(atm.params.cn2_fractions)) =
    phase_variance(atm.params.r0, atm.params.L0; cn2=cn2)

function _phase_covariance_cpu(rho::AbstractArray, r0::Real, L0::Real)
    L0r0 = (L0 / r0)^(5 / 3)
    cst = (24 * gamma(6 / 5) / 5)^(5 / 6) * (gamma(11 / 6) / ((2^(5 / 6)) * pi^(8 / 3))) * L0r0
    base = (24 * gamma(6 / 5) / 5)^(5 / 6) * (gamma(11 / 6) * gamma(5 / 6)) / (2 * pi^(8 / 3)) * L0r0
    out = fill(base, size(rho))
    u = 2 * pi .* rho ./ L0
    mask = rho .!= 0
    out[mask] .= cst .* (u[mask].^(5 / 6)) .* besselk.(5 / 6, u[mask])
    return out
end

function phase_covariance(rho::AbstractArray, r0::Real, L0::Real)
    return _phase_covariance(execution_style(rho), rho, r0, L0)
end

_phase_covariance(::ScalarCPUStyle, rho::AbstractArray, r0::Real, L0::Real) = _phase_covariance_cpu(rho, r0, L0)

function _phase_covariance(style::AcceleratorStyle, rho::AbstractArray, r0::Real, L0::Real)
    # `SpecialFunctions.besselk` is not GPU-compilable on current AMDGPU/CUDA paths,
    # so compute on the host and return the result on the original backend.
    out_host = _phase_covariance_cpu(Array(rho), r0, L0)
    out = similar(rho, eltype(out_host), size(out_host)...)
    copyto!(out, out_host)
    synchronize_backend!(style)
    return out
end

phase_covariance(rho::AbstractArray, atm::KolmogorovAtmosphere) = phase_covariance(rho, atm.params.r0, atm.params.L0)

function phase_spectrum(f::AbstractArray, r0::Real, L0::Real)
    T = promote_type(typeof(float(r0)), typeof(float(L0)), eltype(float.(f)))
    coeff = T(0.023) * T(r0)^(-T(5) / T(3))
    inv_L0_sq = inv(T(L0))^2
    two_pi_sq = T(2 * pi)^2
    exponent = -T(11) / T(6)
    out = similar(f, T, size(f)...)
    _phase_spectrum!(execution_style(out), out, f, coeff, two_pi_sq, inv_L0_sq, exponent)
    return out
end

phase_spectrum(f::AbstractArray, atm::KolmogorovAtmosphere) = phase_spectrum(f, atm.params.r0, atm.params.L0)

function _phase_spectrum!(::ScalarCPUStyle, out::AbstractArray{T}, freqs::AbstractArray{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T) where {T<:AbstractFloat}
    @inbounds for i in eachindex(out, freqs)
        f = freqs[i]
        out[i] = coeff * (two_pi_sq * (f * f) + inv_L0_sq)^exponent
    end
    return out
end

function _phase_spectrum!(style::AcceleratorStyle, out::AbstractArray{T}, freqs::AbstractArray{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T) where {T<:AbstractFloat}
    launch_kernel!(style, phase_spectrum_kernel!, out, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, length(out); ndrange=length(out))
    return out
end

function covariance_matrix(rho1::AbstractVector, rho2::AbstractVector, atm::KolmogorovAtmosphere)
    rho = abs.(rho1 .- transpose(rho2))
    return phase_covariance(rho, atm)
end

function ft_phase_screen(atm::KolmogorovAtmosphere, n::Int, delta::Real;
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false,
    ws::Union{Nothing,PhaseStatsWorkspace}=nothing)

    if ws === nothing
        ws = PhaseStatsWorkspace(atm.state.opd, n; T=eltype(atm.state.opd))
    end
    fftfreq!(ws.freqs, n; d=delta)
    T = eltype(ws.psd)
    coeff = T(0.023) * atm.params.r0^(-T(5) / T(3))
    inv_L0_sq = inv(T(atm.params.L0))^2
    two_pi_sq = T(2 * pi)^2
    exponent = -T(11) / T(6)
    fm = T(5.92) / T(l0) / T(2 * pi)
    inv_fm_sq = isfinite(fm) && fm > zero(T) ? inv(fm)^2 : zero(T)

    fill_phase_psd!(ws.psd, ws.freqs, coeff, two_pi_sq, inv_L0_sq, exponent, inv_fm_sq, n)

    randn_backend!(rng, ws.noise_re)
    randn_backend!(rng, ws.noise_im)
    @. ws.spectrum = complex(ws.noise_re, ws.noise_im) * sqrt(ws.psd)
    execute_fft_plan!(ws.spectrum, ws.ifft_plan)
    phs = real.(ws.spectrum) .* (n * T(delta))
    if return_psd
        return phs, ws.psd
    end
    return phs
end

function fill_phase_psd!(psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T, inv_fm_sq::T, n::Int) where {T<:AbstractFloat}
    _fill_phase_psd!(execution_style(psd), psd, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, inv_fm_sq, n)
    return psd
end

function _fill_phase_psd!(::ScalarCPUStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T, inv_fm_sq::T, n::Int) where {T<:AbstractFloat}
    @inbounds for i in 1:n, j in 1:n
        f_sq = freqs[i]^2 + freqs[j]^2
        base = coeff * (two_pi_sq * f_sq + inv_L0_sq)^exponent
        with_inner = inv_fm_sq == zero(T) ? base : base * exp(-f_sq * inv_fm_sq)
        psd[i, j] = i == 1 && j == 1 ? zero(T) : with_inner
    end
    return psd
end

function _fill_phase_psd!(style::AcceleratorStyle, psd::AbstractMatrix{T}, freqs::AbstractVector{T},
    coeff::T, two_pi_sq::T, inv_L0_sq::T, exponent::T, inv_fm_sq::T, n::Int) where {T<:AbstractFloat}
    launch_kernel!(style, phase_screen_psd_kernel!, psd, freqs, coeff, two_pi_sq, inv_L0_sq, exponent, inv_fm_sq, n; ndrange=size(psd))
    return psd
end

function ft_sh_phase_screen(atm::KolmogorovAtmosphere, n::Int, delta::Real;
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false,
    ws::Union{Nothing,PhaseStatsWorkspace}=nothing, subharmonics::Bool=true,
    profile::FidelityProfile=default_fidelity_profile(),
    mode::SubharmonicMode=default_subharmonic_mode(profile),
    n_levels::Union{Int,Nothing}=nothing, subharmonic_radius::Union{Int,Nothing}=nothing)
    phs, psd = ft_phase_screen(atm, n, delta; l0=l0, rng=rng, return_psd=true, ws=ws)
    if subharmonics
        levels = something(n_levels, default_subharmonic_levels(mode, atm.params.L0, n * delta))
        radius = something(subharmonic_radius, default_subharmonic_radius(mode))
        add_subharmonics!(phs, atm.params.r0, atm.params.L0, delta, l0;
            rng=rng, n_levels=levels, radius=radius)
    end
    if return_psd
        return phs, psd
    end
    return phs
end

default_subharmonic_levels(::FastSubharmonics, L0::Real, D::Real) = 3
default_subharmonic_levels(::FidelitySubharmonics, L0::Real, D::Real) =
    resolve_subharmonic_levels(L0, D)

default_subharmonic_radius(::FastSubharmonics) = 1
default_subharmonic_radius(::FidelitySubharmonics) = 2

default_subharmonic_mode(profile::FidelityProfile) = default_subharmonic_mode(atmosphere_profile(profile))
default_subharmonic_mode(::ScientificProfile) = FidelitySubharmonics()
default_subharmonic_mode(::FastProfile) = FastSubharmonics()

function add_subharmonics!(phs::AbstractMatrix{T}, r0::Real, L0::Real, delta::Real, l0::Real;
    rng::AbstractRNG=Random.default_rng(), n_levels::Int=3, radius::Int=2) where {T<:AbstractFloat}
    n_levels >= 1 || throw(ArgumentError("n_levels must be >= 1"))
    radius >= 1 || throw(ArgumentError("radius must be >= 1"))
    return _add_subharmonics!(execution_style(phs), phs, r0, L0, delta, l0; rng=rng, n_levels=n_levels, radius=radius)
end

function _subharmonic_terms(::Type{T}, n::Int, r0::Real, L0::Real, delta::Real, l0::Real, rng::AbstractRNG,
    n_levels::Int, radius::Int) where {T<:AbstractFloat}
    base_coeff = T(0.023) * T(r0)^(-T(5) / T(3))
    inv_L0_sq = inv(T(L0))^2
    two_pi_sq = T(2 * pi)^2
    exponent = -T(11) / T(6)
    fm = T(5.92) / T(l0) / T(2 * pi)
    inv_fm_sq = isfinite(fm) && fm > zero(T) ? inv(fm)^2 : zero(T)
    D = T(n) * T(delta)
    freq_x = T[]
    freq_y = T[]
    coeff_re = T[]
    coeff_im = T[]
    for p in 1:n_levels
        del_f = inv(D * T(3^p))
        for fx in -radius:radius, fy in -radius:radius
            if fx == 0 && fy == 0
                continue
            end
            fx_t = T(fx) * del_f
            fy_t = T(fy) * del_f
            f = sqrt(fx_t^2 + fy_t^2)
            psd = base_coeff * (two_pi_sq * (f^2) + inv_L0_sq)^exponent
            if inv_fm_sq != zero(T)
                psd *= exp(-(f^2) * inv_fm_sq)
            end
            amp = sqrt(psd) * del_f
            sample_coeff = (randn(rng) + im * randn(rng)) * amp
            push!(freq_x, fx_t)
            push!(freq_y, fy_t)
            push!(coeff_re, real(sample_coeff))
            push!(coeff_im, imag(sample_coeff))
        end
    end
    return (; freq_x, freq_y, coeff_re, coeff_im)
end

function _add_subharmonics!(::ScalarCPUStyle, phs::AbstractMatrix{T}, r0::Real, L0::Real, delta::Real, l0::Real;
    rng::AbstractRNG=Random.default_rng(), n_levels::Int=3, radius::Int=2) where {T<:AbstractFloat}
    n = size(phs, 1)
    offset = n ÷ 2
    delta_t = T(delta)
    terms = _subharmonic_terms(T, n, r0, L0, delta, l0, rng, n_levels, radius)
    @inbounds for k in eachindex(terms.freq_x)
        fx_t = terms.freq_x[k]
        fy_t = terms.freq_y[k]
        re = terms.coeff_re[k]
        im = terms.coeff_im[k]
        for i in 1:n, j in 1:n
            xi = (i - 1 - offset) * delta_t
            yj = (j - 1 - offset) * delta_t
            phase = T(2 * pi) * (fx_t * xi + fy_t * yj)
            phs[i, j] += re * cos(phase) - im * sin(phase)
        end
    end
    return phs
end

function _add_subharmonics!(style::AcceleratorStyle, phs::AbstractMatrix{T}, r0::Real, L0::Real, delta::Real, l0::Real;
    rng::AbstractRNG=Random.default_rng(), n_levels::Int=3, radius::Int=2) where {T<:AbstractFloat}
    n = size(phs, 1)
    terms = _subharmonic_terms(T, n, r0, L0, delta, l0, rng, n_levels, radius)
    n_terms = length(terms.freq_x)
    n_terms == 0 && return phs
    coeff_re = similar(phs, T, n_terms)
    coeff_im = similar(phs, T, n_terms)
    freq_x = similar(phs, T, n_terms)
    freq_y = similar(phs, T, n_terms)
    copyto!(coeff_re, terms.coeff_re)
    copyto!(coeff_im, terms.coeff_im)
    copyto!(freq_x, terms.freq_x)
    copyto!(freq_y, terms.freq_y)
    launch_kernel!(style, add_subharmonics_kernel!, phs, coeff_re, coeff_im, freq_x, freq_y,
        T(2 * pi), T(delta), n ÷ 2, n_terms, n; ndrange=size(phs))
    return phs
end

function resolve_subharmonic_levels(L0::Real, D::Real; min_levels::Int=3, max_levels::Int=6)
    min_levels >= 1 || throw(ArgumentError("min_levels must be >= 1"))
    max_levels >= min_levels || throw(ArgumentError("max_levels must be >= min_levels"))
    ratio = max(float(L0) / float(D), 1.0)
    levels = ceil(Int, log(ratio) / log(3)) + 2
    return clamp(levels, min_levels, max_levels)
end
