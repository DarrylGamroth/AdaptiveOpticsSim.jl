struct SpectralSample{T<:AbstractFloat}
    wavelength::T
    weight::T
end

struct SpectralBundle{T<:AbstractFloat,V<:AbstractVector{SpectralSample{T}}}
    samples::V
end

function SpectralBundle(samples::AbstractVector{<:SpectralSample{T}}) where {T<:AbstractFloat}
    isempty(samples) && throw(InvalidConfiguration("SpectralBundle must contain at least one sample"))
    total_weight = sum(sample.weight for sample in samples)
    total_weight > zero(T) || throw(InvalidConfiguration("SpectralBundle weights must sum to a positive value"))
    normalized = Vector{SpectralSample{T}}(undef, length(samples))
    @inbounds for i in eachindex(samples)
        λ = samples[i].wavelength
        w = samples[i].weight
        λ > zero(T) || throw(InvalidConfiguration("SpectralSample wavelength must be > 0"))
        w >= zero(T) || throw(InvalidConfiguration("SpectralSample weight must be >= 0"))
        normalized[i] = SpectralSample{T}(λ, w / total_weight)
    end
    return SpectralBundle{T, typeof(normalized)}(normalized)
end

function SpectralBundle(wavelengths::AbstractVector, weights::AbstractVector; T::Type{<:AbstractFloat}=Float64)
    length(wavelengths) == length(weights) ||
        throw(DimensionMismatchError("SpectralBundle wavelengths and weights must have the same length"))
    samples = Vector{SpectralSample{T}}(undef, length(wavelengths))
    @inbounds for i in eachindex(wavelengths)
        samples[i] = SpectralSample{T}(T(wavelengths[i]), T(weights[i]))
    end
    return SpectralBundle(samples)
end

Base.length(bundle::SpectralBundle) = length(bundle.samples)
Base.getindex(bundle::SpectralBundle, i::Int) = bundle.samples[i]
Base.iterate(bundle::SpectralBundle, state...) = iterate(bundle.samples, state...)

weighted_wavelength(bundle::SpectralBundle{T}) where {T<:AbstractFloat} =
    sum(sample.wavelength * sample.weight for sample in bundle.samples)

function spectral_bundle_signature(bundle::SpectralBundle{T}) where {T<:AbstractFloat}
    sig = hash(length(bundle.samples), UInt(0))
    @inbounds for sample in bundle.samples
        sig = hash(sample.wavelength, sig)
        sig = hash(sample.weight, sig)
    end
    return sig
end

struct SpectralSource{S<:AbstractSource,B<:SpectralBundle} <: AbstractSource
    source::S
    bundle::B
end

function with_spectrum(src::AbstractSource, bundle::SpectralBundle)
    return SpectralSource(src, bundle)
end

wavelength(src::SpectralSource) = wavelength(src.source)
photon_flux(src::SpectralSource) = photon_flux(src.source)
coordinates_xy_arcsec(src::SpectralSource) = coordinates_xy_arcsec(src.source)
source_height_m(src::SpectralSource) = source_height_m(src.source)
optical_tag(src::SpectralSource) = optical_tag(src.source)
is_lgs_source(src::SpectralSource) = is_lgs_source(src.source)

has_spectral_bundle(::AbstractSource) = false
has_spectral_bundle(::SpectralSource) = true

spectral_bundle(::AbstractSource) = nothing
spectral_bundle(src::SpectralSource) = src.bundle

spectral_reference_source(src::AbstractSource) = src
spectral_reference_source(src::SpectralSource) = src.source

is_polychromatic(src::AbstractSource) = false
is_polychromatic(src::SpectralSource) = length(src.bundle) > 1

function source_measurement_signature(src::AbstractSource)
    return hash((typeof(src), wavelength(src)))
end

function source_measurement_signature(src::SpectralSource)
    return hash((typeof(src.source), wavelength(src.source), spectral_bundle_signature(src.bundle)))
end

function source_with_wavelength_and_flux(src::Source, λ::T, flux::T) where {T<:AbstractFloat}
    params = src.params
    return Source(SourceParams{T}(params.band, T(params.magnitude), (T(params.coordinates_xy_arcsec[1]), T(params.coordinates_xy_arcsec[2])), λ, flux))
end

function source_with_wavelength_and_flux(src::LGSSource, λ::T, flux::T) where {T<:AbstractFloat}
    params = src.params
    lcoords = (T(params.laser_coordinates[1]), T(params.laser_coordinates[2]))
    coords = (T(params.coordinates_xy_arcsec[1]), T(params.coordinates_xy_arcsec[2]))
    return LGSSource(LGSSourceParams{T, typeof(params.na_profile)}(
        T(params.magnitude),
        coords,
        λ,
        T(params.altitude),
        T(params.elongation_factor),
        lcoords,
        params.na_profile,
        T(params.fwhm_spot_up),
        flux,
    ))
end

source_with_wavelength_and_flux(src::SpectralSource, λ::T, flux::T) where {T<:AbstractFloat} =
    source_with_wavelength_and_flux(src.source, λ, flux)
