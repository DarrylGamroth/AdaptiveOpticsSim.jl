using Random
using SpecialFunctions

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
    ifft_plan = plan_ifft_backend!(buffer)
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
phase_variance(atm::MultiLayerAtmosphere; cn2::Real=sum(atm.params.r0_fractions)) =
    phase_variance(atm.params.r0, atm.params.L0; cn2=cn2)

function phase_covariance(rho::AbstractArray, r0::Real, L0::Real)
    L0r0 = (L0 / r0)^(5 / 3)
    cst = (24 * gamma(6 / 5) / 5)^(5 / 6) * (gamma(11 / 6) / ((2^(5 / 6)) * pi^(8 / 3))) * L0r0
    base = (24 * gamma(6 / 5) / 5)^(5 / 6) * (gamma(11 / 6) * gamma(5 / 6)) / (2 * pi^(8 / 3)) * L0r0
    out = fill(base, size(rho))
    u = 2 * pi .* rho ./ L0
    mask = rho .!= 0
    out[mask] .= cst .* (u[mask].^(5 / 6)) .* besselk.(5 / 6, u[mask])
    return out
end

phase_covariance(rho::AbstractArray, atm::KolmogorovAtmosphere) = phase_covariance(rho, atm.params.r0, atm.params.L0)

function phase_spectrum(f::AbstractArray, r0::Real, L0::Real)
    cst = (24 * gamma(6 / 5) / 5)^(5 / 6) * (gamma(11 / 6)^2) / (2 * pi^(11 / 3)) * r0^(-5 / 3)
    return cst .* (f .^ 2 .+ (1 / L0)^2) .^ (-11 / 6)
end

phase_spectrum(f::AbstractArray, atm::KolmogorovAtmosphere) = phase_spectrum(f, atm.params.r0, atm.params.L0)

function covariance_matrix(rho1::AbstractVector, rho2::AbstractVector, atm::KolmogorovAtmosphere)
    rho = abs.(rho1 .- transpose(rho2))
    return phase_covariance(rho, atm)
end

function ft_phase_screen(atm::KolmogorovAtmosphere, n::Int, delta::Real;
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false,
    ws::Union{Nothing,PhaseStatsWorkspace}=nothing)

    if ws === nothing
        ws = PhaseStatsWorkspace(n; T=eltype(atm.state.opd))
    end
    fftfreq!(ws.freqs, n; d=delta)
    del_f = 1 / (n * delta)

    r0 = atm.params.r0
    L0 = atm.params.L0
    fm = 5.92 / l0 / (2 * pi)
    f0 = 1 / L0

    @inbounds for i in 1:n, j in 1:n
        f = sqrt(ws.freqs[i]^2 + ws.freqs[j]^2)
        ws.psd[i, j] = 0.023 * r0^(-5 / 3) * exp(-((f / fm)^2)) / ((f^2 + f0^2)^(11 / 6))
    end
    ws.psd[div(n, 2) + 1, div(n, 2) + 1] = 0

    randn_backend!(rng, ws.noise_re)
    randn_backend!(rng, ws.noise_im)
    @. ws.spectrum = complex(ws.noise_re, ws.noise_im) * sqrt(ws.psd) * del_f
    fftshift2d!(ws.buffer, ws.spectrum)
    execute_fft_plan!(ws.buffer, ws.ifft_plan)
    fftshift2d!(ws.spectrum, ws.buffer)
    phs = real.(ws.spectrum) .* n
    if return_psd
        return phs, ws.psd
    end
    return phs
end

function ft_sh_phase_screen(atm::KolmogorovAtmosphere, n::Int, delta::Real;
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false,
    ws::Union{Nothing,PhaseStatsWorkspace}=nothing, subharmonics::Bool=true, n_levels::Int=3)
    phs, psd = ft_phase_screen(atm, n, delta; l0=l0, rng=rng, return_psd=true, ws=ws)
    if subharmonics
        add_subharmonics!(phs, atm.params.r0, atm.params.L0, delta, l0; rng=rng, n_levels=n_levels)
    end
    if return_psd
        return phs, psd
    end
    return phs
end

function add_subharmonics!(phs::AbstractMatrix{T}, r0::Real, L0::Real, delta::Real, l0::Real;
    rng::AbstractRNG=Random.default_rng(), n_levels::Int=3) where {T<:AbstractFloat}
    n = size(phs, 1)
    fm = 5.92 / l0 / (2 * pi)
    f0 = 1 / L0
    D = n * delta
    offset = n ÷ 2
    delta_t = T(delta)
    for p in 1:n_levels
        del_f = 1 / (D * 3^p)
        for fx in -1:1, fy in -1:1
            if fx == 0 && fy == 0
                continue
            end
            f = sqrt((fx * del_f)^2 + (fy * del_f)^2)
            psd = 0.023 * r0^(-5 / 3) * exp(-((f / fm)^2)) / ((f^2 + f0^2)^(11 / 6))
            amp = sqrt(psd) * del_f
            coeff = (randn(rng) + im * randn(rng)) * amp
            @inbounds for i in 1:n, j in 1:n
                xi = (i - 1 - offset) * delta_t
                yj = (j - 1 - offset) * delta_t
                phase = 2 * pi * ((fx * del_f) * xi + (fy * del_f) * yj)
                phs[i, j] += real(coeff * cis(phase))
            end
        end
    end
    return phs
end
