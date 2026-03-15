const BAND_WAVELENGTHS = Dict(
    :U => 0.365e-6,
    :B => 0.445e-6,
    :V => 0.551e-6,
    :R => 0.658e-6,
    :I => 0.806e-6,
    :J => 1.22e-6,
    :H => 1.63e-6,
    :K => 2.19e-6,
)

const BAND_ZEROPOINTS = Dict(
    :U => 1.96e12 / 368,
    :B => 5.38e12 / 368,
    :V => 3.64e12 / 368,
    :R => 4.01e12 / 368,
    :I => 2.69e12 / 368,
    :J => 1.90e12 / 368,
    :H => 1.10e12 / 368,
    :K => 0.70e12 / 368,
)

abstract type LGSProfile end
struct LGSProfileNone <: LGSProfile end
struct LGSProfileNaProfile <: LGSProfile end

struct SourceParams{T<:AbstractFloat}
    band::Symbol
    magnitude::T
    coordinates::NTuple{2,T}
    wavelength::T
    n_photon::T
end

struct Source{P<:SourceParams} <: AbstractSource
    params::P
end

function Source(; band::Symbol=:I, magnitude::Real=0.0, coordinates=(0.0, 0.0), wavelength=nothing, T::Type{<:AbstractFloat}=Float64)
    if wavelength === nothing
        if haskey(BAND_WAVELENGTHS, band)
            wavelength = BAND_WAVELENGTHS[band]
        else
            throw(InvalidConfiguration("Unknown band $(band); provide wavelength."))
        end
    end
    coords = (T(coordinates[1]), T(coordinates[2]))
    n_photon = if haskey(BAND_ZEROPOINTS, band)
        T(BAND_ZEROPOINTS[band] * 10.0^(-0.4 * magnitude))
    else
        one(T)
    end
    params = SourceParams{T}(band, T(magnitude), coords, T(wavelength), n_photon)
    return Source(params)
end

wavelength(src::Source) = src.params.wavelength
photon_flux(src::Source) = src.params.n_photon

struct LGSSourceParams{T<:AbstractFloat,A}
    magnitude::T
    coordinates::NTuple{2,T}
    wavelength::T
    altitude::T
    elongation_factor::T
    laser_coordinates::NTuple{2,T}
    na_profile::A
    fwhm_spot_up::T
    n_photon::T
end

struct LGSSource{P<:LGSSourceParams} <: AbstractSource
    params::P
end

function LGSSource(; magnitude::Real=0.0, coordinates=(0.0, 0.0), wavelength::Real=589e-9,
    altitude::Real=90000.0, elongation_factor::Real=1.2, laser_coordinates=(0.0, 0.0),
    na_profile=nothing, fwhm_spot_up::Real=0.0, T::Type{<:AbstractFloat}=Float64)

    alt_val = T(altitude)
    if na_profile !== nothing
        if size(na_profile, 1) != 2
            throw(InvalidConfiguration("na_profile must be a 2xN matrix of (altitude, weight)"))
        end
        weights = na_profile[2, :]
        if sum(weights) > 0
            alt_val = T(sum(na_profile[1, :] .* weights) / sum(weights))
        end
    end

    params = LGSSourceParams{T, typeof(na_profile)}(
        T(magnitude),
        (T(coordinates[1]), T(coordinates[2])),
        T(wavelength),
        alt_val,
        T(elongation_factor),
        (T(laser_coordinates[1]), T(laser_coordinates[2])),
        na_profile,
        T(fwhm_spot_up),
        one(T),
    )
    return LGSSource(params)
end

wavelength(src::LGSSource) = src.params.wavelength
photon_flux(src::LGSSource) = src.params.n_photon

lgs_elongation_factor(src::LGSSource) = src.params.elongation_factor

function lgs_has_profile(src::LGSSource)
    return src.params.na_profile !== nothing
end

lgs_profile(::LGSSource{<:LGSSourceParams{<:AbstractFloat,Nothing}}) = LGSProfileNone()
lgs_profile(::LGSSource{<:LGSSourceParams{<:AbstractFloat,<:AbstractMatrix}}) = LGSProfileNaProfile()
