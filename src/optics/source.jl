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

abstract type AbstractSourceRadiometry end

# Source expansions are deliberately flat alternatives. Concrete leaf source
# extensions inherit the leaf default; wrappers that expand wavelength,
# direction, or source extent opt into the expanded trait in their defining
# files. This lets optical models reject unsupported Cartesian products at a
# cold API boundary without an `isa` ladder.
abstract type AbstractSourceCompositionStyle end
struct LeafSourceComposition <: AbstractSourceCompositionStyle end
struct ExpandedSourceComposition <: AbstractSourceCompositionStyle end

@inline source_composition_style(::AbstractSource) = LeafSourceComposition()

@inline require_leaf_source(src::AbstractSource,
    label::AbstractString) = _require_leaf_source(
        source_composition_style(src), label)

@inline _require_leaf_source(::LeafSourceComposition,
    ::AbstractString) = nothing

@inline is_leaf_source(src::AbstractSource) =
    _is_leaf_source(source_composition_style(src))
@inline _is_leaf_source(::LeafSourceComposition) = true
@inline _is_leaf_source(::ExpandedSourceComposition) = false

function _require_leaf_source(::ExpandedSourceComposition,
    label::AbstractString)
    throw(UnsupportedAlgorithm(
        "$(label) does not support a composite source expansion; prepare " *
        "and combine its components explicitly"))
end

"""The source value is physical photon irradiance in photons·s⁻¹·m⁻²."""
struct PhysicalPhotonIrradianceSource <: AbstractSourceRadiometry end

"""The source value is dimensionless relative power for calibration or tests."""
struct NormalizedTestSource <: AbstractSourceRadiometry end

config_value(::PhysicalPhotonIrradianceSource) =
    "physical_photon_irradiance"
config_value(::NormalizedTestSource) = "normalized_test_power"

@inline function _converted_finite(value::Real, ::Type{T},
    label::AbstractString) where {T<:AbstractFloat}
    isfinite(value) || throw(InvalidConfiguration("$(label) must be finite"))
    converted = T(value)
    isfinite(converted) || throw(InvalidConfiguration(
        "$(label) must remain finite after conversion to $(T)"))
    return converted
end


function _converted_finite(::Any, ::Type{T},
    label::AbstractString) where {T<:AbstractFloat}
    throw(InvalidConfiguration("$(label) must be a real finite value"))
end

@inline function _converted_positive_finite(value::Real,
    ::Type{T}, label::AbstractString) where {T<:AbstractFloat}
    converted = T(value)
    isfinite(converted) && converted > zero(T) || throw(
        InvalidConfiguration("$(label) must be finite and positive after conversion to $(T)"))
    return converted
end

function _converted_positive_finite(::Any, ::Type{T},
    label::AbstractString) where {T<:AbstractFloat}
    throw(InvalidConfiguration("$(label) must be a real finite positive value"))
end

@inline function _converted_nonnegative_finite(value::Real,
    ::Type{T}, label::AbstractString) where {T<:AbstractFloat}
    isfinite(value) && value >= zero(value) || throw(
        InvalidConfiguration("$(label) must be finite and non-negative"))
    converted = T(value)
    isfinite(converted) && converted >= zero(T) || throw(
        InvalidConfiguration("$(label) must be finite and non-negative after conversion to $(T)"))
    return converted
end


function _converted_nonnegative_finite(::Any, ::Type{T},
    label::AbstractString) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "$(label) must be a real finite non-negative value"))
end

struct SourceParams{T<:AbstractFloat,R<:AbstractSourceRadiometry}
    band::Symbol
    magnitude::T
    coordinates_xy_arcsec::NTuple{2,T}
    wavelength::T
    radiometric_value::T
    radiometry::R
end

struct Source{P<:SourceParams} <: AbstractSource
    params::P
end

@inline function polar_arcsec_deg_to_xy_arcsec(radius_arcsec, theta_deg,
    ::Type{T}) where {T<:AbstractFloat}
    r = _converted_finite(radius_arcsec, T,
        "source angular-separation radius")
    theta = _converted_finite(theta_deg, T,
        "source position angle")
    θ = _converted_finite(deg2rad(theta), T,
        "source position angle in radians")
    sθ, cθ = sincos(θ)
    x = _converted_finite(r * cθ, T,
        "source x coordinate")
    y = _converted_finite(r * sθ, T,
        "source y coordinate")
    return (x, y)
end

@inline default_source_radiometry(band::Symbol, ::Nothing) =
    haskey(BAND_ZEROPOINTS, band) ? PhysicalPhotonIrradianceSource() :
    NormalizedTestSource()
@inline default_source_radiometry(::Symbol, ::Real) =
    PhysicalPhotonIrradianceSource()
@inline default_source_radiometry(::Symbol,
    radiometry::AbstractSourceRadiometry) = radiometry

function source_radiometric_value(::PhysicalPhotonIrradianceSource,
    band::Symbol, magnitude, ::Nothing, normalized_power, ::Type{T}) where {T}
    haskey(BAND_ZEROPOINTS, band) || throw(InvalidConfiguration(
        "physical custom-band sources require photon_irradiance"))
    value = BAND_ZEROPOINTS[band] * 10.0^(-0.4 * magnitude)
    return _converted_nonnegative_finite(value, T,
        "source photon irradiance")
end

function source_radiometric_value(::PhysicalPhotonIrradianceSource,
    ::Symbol, ::Any, irradiance::Real, normalized_power, ::Type{T}) where {T}
    return _converted_nonnegative_finite(irradiance, T,
        "source photon_irradiance")
end

function source_radiometric_value(::NormalizedTestSource, ::Symbol, ::Any,
    ::Nothing, normalized_power, ::Type{T}) where {T}
    return _converted_nonnegative_finite(normalized_power, T,
        "normalized source power")
end

function source_radiometric_value(::NormalizedTestSource, ::Symbol, ::Any,
    ::Real, normalized_power, ::Type)
    throw(InvalidConfiguration(
        "normalized sources cannot declare physical photon_irradiance"))
end

function Source(; band::Symbol=:I, magnitude::Real=0.0,
    coordinates=(0.0, 0.0), wavelength=nothing,
    photon_irradiance::Union{Nothing,Real}=nothing,
    normalized_power::Real=1.0,
    radiometry::Union{Nothing,AbstractSourceRadiometry}=nothing,
    T::Type{<:AbstractFloat}=Float64)
    if wavelength === nothing
        if haskey(BAND_WAVELENGTHS, band)
            wavelength = BAND_WAVELENGTHS[band]
        else
            throw(InvalidConfiguration("Unknown band $(band); provide wavelength."))
        end
    end
    wavelength_value = _converted_positive_finite(wavelength, T,
        "source wavelength")
    magnitude_value = _converted_finite(magnitude, T, "source magnitude")
    coords_xy_arcsec = polar_arcsec_deg_to_xy_arcsec(coordinates[1], coordinates[2], T)
    resolved_radiometry = default_source_radiometry(band,
        isnothing(radiometry) ? photon_irradiance : radiometry)
    radiometric_value = source_radiometric_value(resolved_radiometry, band,
        magnitude, photon_irradiance, normalized_power, T)
    params = SourceParams{T,typeof(resolved_radiometry)}(band,
        magnitude_value,
        coords_xy_arcsec, wavelength_value, radiometric_value,
        resolved_radiometry)
    return Source(params)
end

wavelength(src::Source) = src.params.wavelength

"""
    photon_irradiance(src)

Return source photon irradiance at the telescope entrance pupil in
photons·s⁻¹·m⁻².
"""
@inline source_radiometry(src::Source) = src.params.radiometry
@inline source_radiometric_value(src::Source) = src.params.radiometric_value
@inline photon_irradiance(src::Source) = _photon_irradiance(
    src.params.radiometry, src.params.radiometric_value)
@inline _photon_irradiance(::PhysicalPhotonIrradianceSource, value) = value
function _photon_irradiance(::NormalizedTestSource, value)
    throw(InvalidConfiguration(
        "normalized test sources do not declare physical photon irradiance"))
end

@inline _require_physical_photon_irradiance(src::AbstractSource,
    label::AbstractString) = _require_physical_photon_irradiance(
    source_radiometry(src), source_radiometric_value(src), label)
@inline _require_physical_photon_irradiance(
    ::PhysicalPhotonIrradianceSource, value, ::AbstractString) = value
function _require_physical_photon_irradiance(::NormalizedTestSource,
    ::Any, label::AbstractString)
    throw(InvalidConfiguration(
        "$(label) requires a source with physical photon irradiance"))
end
optical_tag(src::Source) = "source($(src.params.band))"
source_height_m(::Source) = Inf
is_lgs_source(::AbstractSource) = false

coordinates_xy_arcsec(src::Source) = src.params.coordinates_xy_arcsec

struct LGSSourceParams{T<:AbstractFloat,A,R<:AbstractSourceRadiometry}
    magnitude::T
    coordinates_xy_arcsec::NTuple{2,T}
    wavelength::T
    altitude::T
    elongation_factor::T
    laser_coordinates::NTuple{2,T}
    na_profile::A
    fwhm_spot_up::T
    radiometric_value::T
    radiometry::R
end

struct LGSSource{P<:LGSSourceParams} <: AbstractSource
    params::P
end

function _freeze_lgs_profile(profile::AbstractMatrix,
    ::Type{T}) where {T<:AbstractFloat}
    size(profile, 1) == 2 || throw(InvalidConfiguration(
        "na_profile must be a 2xN matrix of (altitude, weight)"))
    n_samples = size(profile, 2)
    n_samples > 0 || throw(InvalidConfiguration(
        "na_profile must contain at least one altitude sample"))
    raw_profile = Matrix(profile)
    frozen = Matrix{T}(undef, 2, n_samples)
    total_weight = zero(T)
    @inbounds for i in 1:n_samples
        frozen[1, i] = _converted_positive_finite(raw_profile[1, i], T,
            "na_profile altitude")
        weight = _converted_nonnegative_finite(raw_profile[2, i], T,
            "na_profile weight")
        frozen[2, i] = weight
        total_weight += weight
    end
    isfinite(total_weight) && total_weight > zero(T) || throw(
        InvalidConfiguration(
            "na_profile weights must have a finite positive sum"))
    mean_altitude = zero(T)
    @inbounds for i in 1:n_samples
        mean_altitude += frozen[1, i] * (frozen[2, i] / total_weight)
    end
    mean_altitude = _converted_positive_finite(mean_altitude, T,
        "na_profile weighted altitude")
    return frozen, mean_altitude
end

function _freeze_lgs_profile(::Any, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "na_profile must be a 2xN matrix of (altitude, weight)"))
end

function LGSSource(; magnitude::Real=0.0, coordinates=(0.0, 0.0), wavelength::Real=589e-9,
    altitude::Real=90000.0, elongation_factor::Real=1.2, laser_coordinates=(0.0, 0.0),
    na_profile=nothing, fwhm_spot_up::Real=0.0,
    photon_irradiance::Union{Nothing,Real}=nothing,
    normalized_power::Real=1.0,
    radiometry::Union{Nothing,AbstractSourceRadiometry}=nothing,
    T::Type{<:AbstractFloat}=Float64)

    altitude_value = _converted_positive_finite(altitude, T,
        "LGS altitude")
    frozen_profile, alt_val = isnothing(na_profile) ?
        (nothing, altitude_value) : _freeze_lgs_profile(na_profile, T)
    magnitude_value = _converted_finite(magnitude, T, "LGS magnitude")
    coordinates_value = polar_arcsec_deg_to_xy_arcsec(
        coordinates[1], coordinates[2], T)
    elongation_value = _converted_nonnegative_finite(elongation_factor, T,
        "LGS elongation factor")
    laser_coordinates_value = (
        _converted_finite(laser_coordinates[1], T,
            "LGS laser x coordinate"),
        _converted_finite(laser_coordinates[2], T,
            "LGS laser y coordinate"),
    )
    fwhm_value = _converted_nonnegative_finite(fwhm_spot_up, T,
        "LGS uplink spot FWHM")

    resolved_radiometry = isnothing(radiometry) ?
        (isnothing(photon_irradiance) ? NormalizedTestSource() :
         PhysicalPhotonIrradianceSource()) : radiometry
    radiometric_value = source_radiometric_value(resolved_radiometry,
        :LGS, magnitude, photon_irradiance, normalized_power, T)
    wavelength_value = _converted_positive_finite(wavelength, T,
        "LGS wavelength")
    params = LGSSourceParams{T, typeof(frozen_profile),
        typeof(resolved_radiometry)}(
        magnitude_value,
        coordinates_value,
        wavelength_value,
        alt_val,
        elongation_value,
        laser_coordinates_value,
        frozen_profile,
        fwhm_value,
        radiometric_value,
        resolved_radiometry,
    )
    return LGSSource(params)
end

wavelength(src::LGSSource) = src.params.wavelength
@inline source_radiometry(src::LGSSource) = src.params.radiometry
@inline source_radiometric_value(src::LGSSource) = src.params.radiometric_value
@inline photon_irradiance(src::LGSSource) = _photon_irradiance(
    src.params.radiometry, src.params.radiometric_value)
optical_tag(::LGSSource) = "lgs"
source_height_m(src::LGSSource) = src.params.altitude
is_lgs_source(::LGSSource) = true

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

"""
    freeze_source(source)

Return a run-owned source description. Maintained source types containing
mutable user inputs specialize this function and copy those inputs.
"""
@inline freeze_source(src::AbstractSource) = src
@inline freeze_source(src::Source) = src

function freeze_source(src::LGSSource)
    params = src.params
    profile = isnothing(params.na_profile) ? nothing : copy(params.na_profile)
    return LGSSource(LGSSourceParams{
        typeof(params.altitude),typeof(profile),typeof(params.radiometry),
    }(
        params.magnitude,
        params.coordinates_xy_arcsec,
        params.wavelength,
        params.altitude,
        params.elongation_factor,
        params.laser_coordinates,
        profile,
        params.fwhm_spot_up,
        params.radiometric_value,
        params.radiometry,
    ))
end
