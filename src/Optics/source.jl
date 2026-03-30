const BAND_WAVELENGTHS = Dict(
    :U => 0.360e-6,
    :B => 0.440e-6,
    :V => 0.550e-6,
    :R => 0.640e-6,
    :I => 0.790e-6,
    :J => 1.215e-6,
    :H => 1.654e-6,
    :K => 2.179e-6,
)

const BAND_ZEROPOINTS = Dict(
    :U => 1.96e12 / 368,
    :B => 5.38e12 / 368,
    :V => 3.31e12 / 368,
    :R => 4.01e12 / 368,
    :I => 2.69e12 / 368,
    :J => 1.90e12 / 368,
    :H => 1.05e12 / 368,
    :K => 0.70e12 / 368,
)

abstract type LGSProfile end
struct LGSProfileNone <: LGSProfile end
struct LGSProfileNaProfile <: LGSProfile end

struct SourceParams{T<:AbstractFloat}
    band::Symbol
    magnitude::T
    coordinates_xy_arcsec::NTuple{2,T}
    wavelength::T
    n_photon::T
end

struct Source{P<:SourceParams} <: AbstractSource
    params::P
end

@inline function polar_arcsec_deg_to_xy_arcsec(radius_arcsec::Real, theta_deg::Real, ::Type{T}) where {T<:AbstractFloat}
    r = T(radius_arcsec)
    θ = T(deg2rad(theta_deg))
    sθ, cθ = sincos(θ)
    return (r * cθ, r * sθ)
end

function Source(; band::Symbol=:I, magnitude::Real=0.0, coordinates=(0.0, 0.0), wavelength=nothing, T::Type{<:AbstractFloat}=Float64)
    if wavelength === nothing
        if haskey(BAND_WAVELENGTHS, band)
            wavelength = BAND_WAVELENGTHS[band]
        else
            throw(InvalidConfiguration("Unknown band $(band); provide wavelength."))
        end
    end
    coords_xy_arcsec = polar_arcsec_deg_to_xy_arcsec(coordinates[1], coordinates[2], T)
    n_photon = if haskey(BAND_ZEROPOINTS, band)
        T(BAND_ZEROPOINTS[band] * 10.0^(-0.4 * magnitude))
    else
        one(T)
    end
    params = SourceParams{T}(band, T(magnitude), coords_xy_arcsec, T(wavelength), n_photon)
    return Source(params)
end

wavelength(src::Source) = src.params.wavelength
photon_flux(src::Source) = src.params.n_photon
optical_tag(src::Source) = "source($(src.params.band))"
source_height_m(::Source) = Inf

coordinates_xy_arcsec(src::Source) = src.params.coordinates_xy_arcsec

struct LGSSourceParams{T<:AbstractFloat,A}
    magnitude::T
    coordinates_xy_arcsec::NTuple{2,T}
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
        polar_arcsec_deg_to_xy_arcsec(coordinates[1], coordinates[2], T),
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
optical_tag(::LGSSource) = "lgs"
source_height_m(src::LGSSource) = src.params.altitude

coordinates_xy_arcsec(src::LGSSource) = src.params.coordinates_xy_arcsec

optical_tag(x) = lowercase(string(nameof(typeof(x))))

optical_path(parts...) = join(optical_tag.(parts), " -> ")
print_optical_path(io::IO, parts...) = print(io, optical_path(parts...))
print_optical_path(parts...) = print_optical_path(stdout, parts...)

lgs_elongation_factor(src::LGSSource) = src.params.elongation_factor

function lgs_has_profile(src::LGSSource)
    return src.params.na_profile !== nothing
end

lgs_profile(::LGSSource{<:LGSSourceParams{<:AbstractFloat,Nothing}}) = LGSProfileNone()
lgs_profile(::LGSSource{<:LGSSourceParams{<:AbstractFloat,<:AbstractMatrix}}) = LGSProfileNaProfile()
