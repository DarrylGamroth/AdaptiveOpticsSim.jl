using FFTW
using Random
using SpecialFunctions

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
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false)

    freqs = zeros(Float64, n)
    fftfreq!(freqs, n; d=delta)
    del_f = 1 / (n * delta)

    r0 = atm.params.r0
    L0 = atm.params.L0
    fm = 5.92 / l0 / (2 * pi)
    f0 = 1 / L0
    psd = zeros(Float64, n, n)

    @inbounds for i in 1:n, j in 1:n
        f = sqrt(freqs[i]^2 + freqs[j]^2)
        psd[i, j] = 0.023 * r0^(-5 / 3) * exp(-((f / fm)^2)) / ((f^2 + f0^2)^(11 / 6))
    end
    psd[div(n, 2) + 1, div(n, 2) + 1] = 0

    noise_re = randn(rng, n, n)
    noise_im = randn(rng, n, n)
    cn = complex.(noise_re, noise_im) .* sqrt.(psd) .* del_f
    phs = real(ifft(FFTW.fftshift(cn))) * n
    if return_psd
        return phs, psd
    end
    return phs
end

function ft_sh_phase_screen(atm::KolmogorovAtmosphere, n::Int, delta::Real;
    l0::Real=1e-10, rng::AbstractRNG=Random.default_rng(), return_psd::Bool=false)
    # Sub-harmonics are omitted for now; this is a frequency-domain screen.
    return ft_phase_screen(atm, n, delta; l0=l0, rng=rng, return_psd=return_psd)
end
