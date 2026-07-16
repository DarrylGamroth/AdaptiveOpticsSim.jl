struct SpectralSample{T<:AbstractFloat}
    wavelength::T
    weight::T
end

function _validated_spectral_weight_sum(
    samples::AbstractVector{<:SpectralSample{T}},
) where {T<:AbstractFloat}
    isempty(samples) && throw(InvalidConfiguration(
        "SpectralBundle must contain at least one sample"))
    total_weight = zero(T)
    @inbounds for sample in samples
        isfinite(sample.wavelength) && sample.wavelength > zero(T) ||
            throw(InvalidConfiguration(
                "SpectralSample wavelength must be finite and > 0"))
        isfinite(sample.weight) && sample.weight >= zero(T) ||
            throw(InvalidConfiguration(
                "SpectralSample weight must be finite and >= 0"))
        total_weight += sample.weight
    end
    isfinite(total_weight) && total_weight > zero(T) ||
        throw(InvalidConfiguration(
            "SpectralBundle weights must have a finite positive sum"))
    return total_weight
end

struct SpectralBundle{T<:AbstractFloat,V<:AbstractVector{SpectralSample{T}}}
    samples::V

    # The fully parameterized form is a validated, already-normalized storage
    # constructor. The ordinary outer constructor below accepts raw weights and
    # normalizes them before entering here.
    function SpectralBundle{T,V}(samples::V) where {
        T<:AbstractFloat,V<:AbstractVector{SpectralSample{T}},
    }
        total_weight = _validated_spectral_weight_sum(samples)
        isapprox(total_weight, one(T); atol=sqrt(eps(T)), rtol=sqrt(eps(T))) ||
            throw(InvalidConfiguration(
                "fully parameterized SpectralBundle weights must already sum to one"))
        return new{T,V}(samples)
    end
end

function SpectralBundle(samples::AbstractVector{<:SpectralSample{T}}) where {T<:AbstractFloat}
    total_weight = _validated_spectral_weight_sum(samples)
    normalized = Vector{SpectralSample{T}}(undef, length(samples))
    @inbounds for i in eachindex(samples)
        λ = samples[i].wavelength
        w = samples[i].weight
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

struct _SpectralSourceToken end
const _SPECTRAL_SOURCE_TOKEN = _SpectralSourceToken()

struct SpectralSource{
    S<:Union{Source,LGSSource},B<:SpectralBundle,
} <: AbstractSource
    source::S
    bundle::B

    function SpectralSource(::_SpectralSourceToken, source::S,
        bundle::B) where {S<:Union{Source,LGSSource},B<:SpectralBundle}
        return new{S,B}(source, bundle)
    end
end

function with_spectrum(src::AbstractSource, bundle::SpectralBundle)
    return _with_spectrum(source_composition_style(src), src, bundle)
end

function _with_spectrum(::LeafSourceComposition,
    src::Union{Source,LGSSource}, bundle::SpectralBundle)
    return SpectralSource(_SPECTRAL_SOURCE_TOKEN, freeze_source(src),
        freeze_spectral_bundle(bundle))
end

function _with_spectrum(::LeafSourceComposition,
    ::AbstractSource, ::SpectralBundle)
    throw(UnsupportedAlgorithm(
        "spectral expansion currently supports Source and LGSSource leaves only"))
end

function _with_spectrum(::ExpandedSourceComposition,
    ::AbstractSource, ::SpectralBundle)
    throw(UnsupportedAlgorithm(
        "nested source expansions are not implemented; prepare and combine " *
        "the spectral, spatial, or directional components explicitly"))
end

function SpectralSource(::AbstractSource, ::SpectralBundle)
    throw(UnsupportedAlgorithm(
        "construct spectral sources with with_spectrum from a Source or " *
        "LGSSource leaf"))
end

function freeze_spectral_bundle(bundle::SpectralBundle{T}) where {T<:AbstractFloat}
    samples = copy(bundle.samples)
    return SpectralBundle{T,typeof(samples)}(samples)
end

freeze_source(src::SpectralSource) = SpectralSource(_SPECTRAL_SOURCE_TOKEN,
    freeze_source(src.source), freeze_spectral_bundle(src.bundle))

@inline source_composition_style(::SpectralSource) =
    ExpandedSourceComposition()

wavelength(src::SpectralSource) = wavelength(src.source)
photon_irradiance(src::SpectralSource) = photon_irradiance(src.source)
source_radiometry(src::SpectralSource) = source_radiometry(src.source)
source_radiometric_value(src::SpectralSource) =
    source_radiometric_value(src.source)
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

function source_with_wavelength_and_radiometric_value(src::Source, λ::T,
    radiometric_value::T) where {T<:AbstractFloat}
    params = src.params
    coordinates = (T(params.coordinates_xy_arcsec[1]),
        T(params.coordinates_xy_arcsec[2]))
    return Source(SourceParams{T,typeof(params.radiometry)}(params.band,
        T(params.magnitude), coordinates, λ, radiometric_value,
        params.radiometry))
end

function source_with_wavelength_and_radiometric_value(src::LGSSource, λ::T,
    radiometric_value::T) where {T<:AbstractFloat}
    params = src.params
    lcoords = (T(params.laser_coordinates[1]), T(params.laser_coordinates[2]))
    coords = (T(params.coordinates_xy_arcsec[1]), T(params.coordinates_xy_arcsec[2]))
    profile = isnothing(params.na_profile) ? nothing : copy(params.na_profile)
    return LGSSource(LGSSourceParams{T, typeof(profile),
        typeof(params.radiometry)}(
        T(params.magnitude),
        coords,
        λ,
        T(params.altitude),
        T(params.elongation_factor),
        lcoords,
        profile,
        T(params.fwhm_spot_up),
        radiometric_value,
        params.radiometry,
    ))
end

source_with_wavelength_and_radiometric_value(src::SpectralSource, λ::T,
    radiometric_value::T) where {T<:AbstractFloat} =
    source_with_wavelength_and_radiometric_value(src.source, λ,
        radiometric_value)
