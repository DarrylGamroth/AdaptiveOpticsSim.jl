abstract type AbstractExtendedSourceModel end

@inline function _source_model_storage_matches(a::T, b::T) where {T<:AbstractFloat}
    return isapprox(a, b; atol=zero(T), rtol=T(8) * eps(T))
end

function _validate_normalized_source_model_storage(
    offsets_xy_arcsec::AbstractVector{<:NTuple{2,T}},
    weights::AbstractVector{T}, label::AbstractString,
) where {T<:AbstractFloat}
    length(offsets_xy_arcsec) == length(weights) || throw(
        DimensionMismatchError(
            "$(label) offsets and weights must have the same length"))
    isempty(offsets_xy_arcsec) && throw(InvalidConfiguration(
        "$(label) must contain at least one sample"))
    total = zero(T)
    @inbounds for i in eachindex(offsets_xy_arcsec, weights)
        x, y = offsets_xy_arcsec[i]
        isfinite(x) && isfinite(y) || throw(InvalidConfiguration(
            "$(label) offsets must be finite"))
        weight = weights[i]
        isfinite(weight) && weight >= zero(T) || throw(
            InvalidConfiguration(
                "$(label) weights must be finite and non-negative"))
        total += weight
    end
    normalization_tolerance = min(sqrt(eps(T)),
        T(8) * eps(T) * T(length(weights)))
    isfinite(total) && abs(total - one(T)) <= normalization_tolerance || throw(
        InvalidConfiguration(
            "fully parameterized $(label) weights must already sum to one"))
    return nothing
end

struct PointCloudSourceModel{T<:AbstractFloat,V<:AbstractVector{NTuple{2,T}},W<:AbstractVector{T}} <: AbstractExtendedSourceModel
    offsets_xy_arcsec::V
    weights::W

    function PointCloudSourceModel{T,V,W}(offsets_xy_arcsec::V,
        weights::W) where {
        T<:AbstractFloat,V<:AbstractVector{NTuple{2,T}},
        W<:AbstractVector{T},
    }
        _validate_normalized_source_model_storage(offsets_xy_arcsec, weights,
            "PointCloudSourceModel")
        return new{T,V,W}(offsets_xy_arcsec, weights)
    end
end

function _normalize_source_weights(offsets_xy_arcsec::AbstractVector{<:NTuple{2,T}}, weights::AbstractVector{<:Real}) where {T<:AbstractFloat}
    length(offsets_xy_arcsec) == length(weights) ||
        throw(DimensionMismatchError("source-model offsets and weights must have the same length"))
    isempty(offsets_xy_arcsec) && throw(InvalidConfiguration("source model must contain at least one sample"))
    total = zero(T)
    normalized_offsets = Vector{NTuple{2,T}}(undef, length(offsets_xy_arcsec))
    normalized_weights = Vector{T}(undef, length(weights))
    @inbounds for i in eachindex(offsets_xy_arcsec, weights)
        offset = offsets_xy_arcsec[i]
        w = T(weights[i])
        isfinite(w) && w >= zero(T) || throw(InvalidConfiguration(
            "source-model sample weights must be finite and >= 0"))
        x = T(offset[1])
        y = T(offset[2])
        isfinite(x) && isfinite(y) || throw(InvalidConfiguration(
            "source-model offsets must remain finite after conversion"))
        normalized_offsets[i] = (x, y)
        normalized_weights[i] = w
        total += w
    end
    isfinite(total) && total > zero(T) || throw(InvalidConfiguration(
        "source-model sample weights must have a finite positive sum"))
    @inbounds for i in eachindex(normalized_weights)
        normalized_weights[i] /= total
    end
    return normalized_offsets, normalized_weights
end

function PointCloudSourceModel(offsets_xy_arcsec::AbstractVector, weights::AbstractVector;
    T::Type{<:AbstractFloat}=Float64)
    offsets = Vector{NTuple{2,T}}(undef, length(offsets_xy_arcsec))
    @inbounds for i in eachindex(offsets_xy_arcsec)
        offset = offsets_xy_arcsec[i]
        length(offset) == 2 || throw(DimensionMismatchError("point-cloud source offsets must contain x/y pairs"))
        offsets[i] = (T(offset[1]), T(offset[2]))
    end
    normalized_offsets, normalized_weights = _normalize_source_weights(offsets, weights)
    return PointCloudSourceModel{T, typeof(normalized_offsets), typeof(normalized_weights)}(normalized_offsets, normalized_weights)
end

function PointCloudSourceModel(offsets_xy_arcsec::AbstractMatrix, weights::AbstractVector; T::Type{<:AbstractFloat}=Float64)
    size(offsets_xy_arcsec, 2) == 2 ||
        throw(DimensionMismatchError("point-cloud source offsets must be an N×2 matrix"))
    offsets = NTuple{2,T}[(T(offsets_xy_arcsec[i, 1]), T(offsets_xy_arcsec[i, 2])) for i in axes(offsets_xy_arcsec, 1)]
    return PointCloudSourceModel(offsets, weights; T=T)
end

struct GaussianDiskSourceModel{T<:AbstractFloat,V<:AbstractVector{NTuple{2,T}},W<:AbstractVector{T}} <: AbstractExtendedSourceModel
    sigma_arcsec::T
    radius_sigma::T
    n_side::Int
    offsets_xy_arcsec::V
    weights::W

    function GaussianDiskSourceModel{T,V,W}(sigma_arcsec::T,
        radius_sigma::T, n_side::Int, offsets_xy_arcsec::V,
        weights::W) where {
        T<:AbstractFloat,V<:AbstractVector{NTuple{2,T}},
        W<:AbstractVector{T},
    }
        isfinite(sigma_arcsec) && sigma_arcsec > zero(T) || throw(
            InvalidConfiguration(
                "GaussianDiskSourceModel sigma_arcsec must be finite and > 0"))
        isfinite(radius_sigma) && radius_sigma > zero(T) || throw(
            InvalidConfiguration(
                "GaussianDiskSourceModel radius_sigma must be finite and > 0"))
        n_side >= 3 && isodd(n_side) || throw(InvalidConfiguration(
            "GaussianDiskSourceModel n_side must be an odd integer >= 3"))
        n_side <= isqrt(typemax(Int)) || throw(InvalidConfiguration(
            "GaussianDiskSourceModel n_side is too large"))
        length(offsets_xy_arcsec) == n_side * n_side || throw(
            DimensionMismatchError(
                "GaussianDiskSourceModel must contain n_side^2 samples"))
        _validate_normalized_source_model_storage(offsets_xy_arcsec, weights,
            "GaussianDiskSourceModel")

        half = div(n_side - 1, 2)
        pitch = radius_sigma * sigma_arcsec / T(half)
        isfinite(pitch) || throw(InvalidConfiguration(
            "GaussianDiskSourceModel sample pitch must be finite"))
        raw_total = zero(T)
        @inbounds for iy in -half:half, ix in -half:half
            x = T(ix) * pitch
            y = T(iy) * pitch
            raw_total += exp(-T(0.5) * (x * x + y * y) /
                             (sigma_arcsec * sigma_arcsec))
        end
        isfinite(raw_total) && raw_total > zero(T) || throw(
            InvalidConfiguration(
                "GaussianDiskSourceModel canonical weights must have a finite positive sum"))

        idx = 1
        @inbounds for iy in -half:half, ix in -half:half
            expected_x = T(ix) * pitch
            expected_y = T(iy) * pitch
            expected_weight = exp(-T(0.5) *
                (expected_x * expected_x + expected_y * expected_y) /
                (sigma_arcsec * sigma_arcsec)) / raw_total
            x, y = offsets_xy_arcsec[idx]
            _source_model_storage_matches(x, expected_x) &&
                _source_model_storage_matches(y, expected_y) &&
                _source_model_storage_matches(weights[idx], expected_weight) ||
                throw(InvalidConfiguration(
                    "fully parameterized GaussianDiskSourceModel storage must match its declared grid"))
            idx += 1
        end
        return new{T,V,W}(sigma_arcsec, radius_sigma, n_side,
            offsets_xy_arcsec, weights)
    end
end

function GaussianDiskSourceModel(; sigma_arcsec::Real, radius_sigma::Real=2.0, n_side::Int=5,
    T::Type{<:AbstractFloat}=Float64)
    sigma = T(sigma_arcsec)
    sigma > zero(T) || throw(InvalidConfiguration("GaussianDiskSourceModel sigma_arcsec must be > 0"))
    radius = T(radius_sigma)
    radius > zero(T) || throw(InvalidConfiguration("GaussianDiskSourceModel radius_sigma must be > 0"))
    n_side >= 3 && isodd(n_side) ||
        throw(InvalidConfiguration("GaussianDiskSourceModel n_side must be an odd integer >= 3"))
    half = div(n_side - 1, 2)
    pitch = radius * sigma / T(half)
    offsets = Vector{NTuple{2,T}}(undef, n_side * n_side)
    weights = Vector{T}(undef, n_side * n_side)
    idx = 1
    @inbounds for iy in -half:half, ix in -half:half
        x = T(ix) * pitch
        y = T(iy) * pitch
        offsets[idx] = (x, y)
        weights[idx] = exp(-T(0.5) * (x * x + y * y) / (sigma * sigma))
        idx += 1
    end
    normalized_offsets, normalized_weights = _normalize_source_weights(offsets, weights)
    return GaussianDiskSourceModel{T, typeof(normalized_offsets), typeof(normalized_weights)}(
        sigma, radius, n_side, normalized_offsets, normalized_weights)
end

struct SampledImageSourceModel{T<:AbstractFloat,M<:AbstractMatrix{T},V<:AbstractVector{NTuple{2,T}},W<:AbstractVector{T}} <: AbstractExtendedSourceModel
    image::M
    pixel_scale_arcsec::T
    offsets_xy_arcsec::V
    weights::W

    function SampledImageSourceModel{T,M,V,W}(image::M,
        pixel_scale_arcsec::T, offsets_xy_arcsec::V, weights::W) where {
        T<:AbstractFloat,M<:AbstractMatrix{T},
        V<:AbstractVector{NTuple{2,T}},W<:AbstractVector{T},
    }
        size(image, 1) > 0 && size(image, 2) > 0 || throw(
            InvalidConfiguration(
                "SampledImageSourceModel image must be non-empty"))
        isfinite(pixel_scale_arcsec) && pixel_scale_arcsec > zero(T) || throw(
            InvalidConfiguration(
                "SampledImageSourceModel pixel_scale_arcsec must be finite and > 0"))
        _validate_normalized_source_model_storage(offsets_xy_arcsec, weights,
            "SampledImageSourceModel")

        total = zero(T)
        positive_count = 0
        @inbounds for value in image
            isfinite(value) && value >= zero(T) || throw(
                InvalidConfiguration(
                    "SampledImageSourceModel image weights must be finite and non-negative"))
            if value > zero(T)
                total += value
                positive_count += 1
            end
        end
        isfinite(total) && total > zero(T) || throw(InvalidConfiguration(
            "SampledImageSourceModel image weights must have a finite positive sum"))
        positive_count == length(offsets_xy_arcsec) || throw(
            DimensionMismatchError(
                "SampledImageSourceModel must contain one sample per positive image pixel"))

        row_center = (T(size(image, 1)) + one(T)) / T(2)
        col_center = (T(size(image, 2)) + one(T)) / T(2)
        idx = 1
        @inbounds for i in axes(image, 1), j in axes(image, 2)
            value = image[i, j]
            value > zero(T) || continue
            expected_x = (T(j) - col_center) * pixel_scale_arcsec
            expected_y = (T(i) - row_center) * pixel_scale_arcsec
            x, y = offsets_xy_arcsec[idx]
            _source_model_storage_matches(x, expected_x) &&
                _source_model_storage_matches(y, expected_y) &&
                _source_model_storage_matches(weights[idx], value / total) ||
                throw(InvalidConfiguration(
                    "fully parameterized SampledImageSourceModel storage must match its image grid"))
            idx += 1
        end
        return new{T,M,V,W}(image, pixel_scale_arcsec,
            offsets_xy_arcsec, weights)
    end
end

function SampledImageSourceModel(image::AbstractMatrix; pixel_scale_arcsec::Real, T::Type{<:AbstractFloat}=Float64)
    size(image, 1) > 0 && size(image, 2) > 0 ||
        throw(InvalidConfiguration("SampledImageSourceModel image must be non-empty"))
    scale = T(pixel_scale_arcsec)
    scale > zero(T) || throw(InvalidConfiguration("SampledImageSourceModel pixel_scale_arcsec must be > 0"))
    img = T.(image)
    any(x -> x < zero(T), img) && throw(InvalidConfiguration("SampledImageSourceModel image weights must be >= 0"))
    row_center = (T(size(img, 1)) + one(T)) / T(2)
    col_center = (T(size(img, 2)) + one(T)) / T(2)
    offsets = NTuple{2,T}[]
    weights = T[]
    @inbounds for i in axes(img, 1), j in axes(img, 2)
        w = img[i, j]
        if w > zero(T)
            push!(offsets, ((T(j) - col_center) * scale, (T(i) - row_center) * scale))
            push!(weights, w)
        end
    end
    normalized_offsets, normalized_weights = _normalize_source_weights(offsets, weights)
    return SampledImageSourceModel{T, typeof(img), typeof(normalized_offsets), typeof(normalized_weights)}(
        img, scale, normalized_offsets, normalized_weights)
end

@inline source_model_offsets(model::PointCloudSourceModel) = model.offsets_xy_arcsec
@inline source_model_offsets(model::GaussianDiskSourceModel) = model.offsets_xy_arcsec
@inline source_model_offsets(model::SampledImageSourceModel) = model.offsets_xy_arcsec

@inline source_model_weights(model::PointCloudSourceModel) = model.weights
@inline source_model_weights(model::GaussianDiskSourceModel) = model.weights
@inline source_model_weights(model::SampledImageSourceModel) = model.weights

function source_model_signature(model::AbstractExtendedSourceModel)
    sig = hash(typeof(model), UInt(0))
    @inbounds for offset in source_model_offsets(model)
        sig = hash(offset[1], sig)
        sig = hash(offset[2], sig)
    end
    @inbounds for weight in source_model_weights(model)
        sig = hash(weight, sig)
    end
    return sig
end

struct ExtendedSource{
    S<:Union{Source,LGSSource},
    M<:AbstractExtendedSourceModel,
    A<:Asterism,
} <: AbstractSource
    source::S
    model::M
    quadrature::A
end

function _require_asterism_leaf_source(::ExtendedSource)
    throw(UnsupportedAlgorithm(
        "Asterism does not support ExtendedSource children; explicit " *
        "spatial-by-directional Cartesian quadrature is not implemented"))
end

function with_spectrum(::ExtendedSource, ::SpectralBundle)
    throw(UnsupportedAlgorithm(
        "polychromatic extended-source composition is not implemented; " *
        "prepare and combine wavelength-by-spatial components explicitly"))
end

function with_extended_source(src::Union{Source,LGSSource}, model::AbstractExtendedSourceModel)
    frozen_source = freeze_source(src)
    frozen_model = freeze_extended_source_model(model)
    quadrature = _build_extended_source_asterism(frozen_source, frozen_model)
    return ExtendedSource{typeof(frozen_source),typeof(frozen_model),
        typeof(quadrature)}(frozen_source, frozen_model, quadrature)
end

function with_extended_source(src::AbstractSource, ::AbstractExtendedSourceModel)
    throw(UnsupportedAlgorithm("extended-source models currently support Source and LGSSource only"))
end

wavelength(src::ExtendedSource) = wavelength(src.source)
photon_irradiance(src::ExtendedSource) = photon_irradiance(src.source)
source_radiometry(src::ExtendedSource) = source_radiometry(src.source)
source_radiometric_value(src::ExtendedSource) =
    source_radiometric_value(src.source)
coordinates_xy_arcsec(src::ExtendedSource) = coordinates_xy_arcsec(src.source)
source_height_m(src::ExtendedSource) = source_height_m(src.source)
optical_tag(src::ExtendedSource) = "extended(" * optical_tag(src.source) * ")"
is_lgs_source(src::ExtendedSource) = is_lgs_source(src.source)
@inline source_composition_style(::ExtendedSource) =
    ExpandedSourceComposition()

has_extended_source_model(::AbstractSource) = false
has_extended_source_model(::ExtendedSource) = true

extended_source_model(::AbstractSource) = nothing
extended_source_model(src::ExtendedSource) = src.model

extended_source_reference_source(src::AbstractSource) = src
extended_source_reference_source(src::ExtendedSource) = src.source

function source_measurement_signature(src::ExtendedSource)
    return hash((source_measurement_signature(src.source), source_model_signature(src.model)))
end

function source_with_coordinates_and_radiometric_value(src::Source,
    coords_xy_arcsec::NTuple{2,T}, radiometric_value::T) where {T<:AbstractFloat}
    params = src.params
    λ = _converted_positive_finite(params.wavelength, T,
        "extended-source wavelength")
    value = _converted_nonnegative_finite(radiometric_value, T,
        "extended-source radiometric value")
    return Source(SourceParams{T,typeof(params.radiometry)}(params.band,
        T(params.magnitude), coords_xy_arcsec, λ, value, params.radiometry))
end

function source_with_coordinates_and_radiometric_value(src::LGSSource,
    coords_xy_arcsec::NTuple{2,T}, radiometric_value::T) where {T<:AbstractFloat}
    params = src.params
    profile = isnothing(params.na_profile) ? nothing : copy(params.na_profile)
    λ = _converted_positive_finite(params.wavelength, T,
        "extended LGS wavelength")
    value = _converted_nonnegative_finite(radiometric_value, T,
        "extended LGS radiometric value")
    return LGSSource(LGSSourceParams{T, typeof(profile),
        typeof(params.radiometry)}(
        T(params.magnitude),
        coords_xy_arcsec,
        λ,
        T(params.altitude),
        T(params.elongation_factor),
        (T(params.laser_coordinates[1]), T(params.laser_coordinates[2])),
        profile,
        T(params.fwhm_spot_up),
        value,
        params.radiometry,
    ))
end

function _build_extended_source_asterism(src::Union{Source,LGSSource}, model::AbstractExtendedSourceModel)
    offsets = source_model_offsets(model)
    weights = source_model_weights(model)
    raw_coords = coordinates_xy_arcsec(src)
    T = promote_type(typeof(raw_coords[1]), eltype(weights))
    base_coords = (T(raw_coords[1]), T(raw_coords[2]))
    total_value = T(source_radiometric_value(src))
    first_value = total_value * weights[1]
    samples = Vector{typeof(source_with_coordinates_and_radiometric_value(
        src, base_coords, first_value))}(undef, length(weights))
    @inbounds for i in eachindex(offsets, weights)
        offset = offsets[i]
        coords = (base_coords[1] + T(offset[1]),
            base_coords[2] + T(offset[2]))
        samples[i] = source_with_coordinates_and_radiometric_value(src,
            coords, total_value * T(weights[i]))
    end
    return Asterism(samples)
end

function freeze_extended_source_model(model::PointCloudSourceModel{T}) where {T<:AbstractFloat}
    offsets = copy(model.offsets_xy_arcsec)
    weights = copy(model.weights)
    return PointCloudSourceModel{T,typeof(offsets),typeof(weights)}(offsets, weights)
end

function freeze_extended_source_model(model::GaussianDiskSourceModel{T}) where {T<:AbstractFloat}
    offsets = copy(model.offsets_xy_arcsec)
    weights = copy(model.weights)
    return GaussianDiskSourceModel{T,typeof(offsets),typeof(weights)}(
        model.sigma_arcsec, model.radius_sigma, model.n_side, offsets, weights)
end

function freeze_extended_source_model(model::SampledImageSourceModel{T}) where {T<:AbstractFloat}
    image = copy(model.image)
    offsets = copy(model.offsets_xy_arcsec)
    weights = copy(model.weights)
    return SampledImageSourceModel{T,typeof(image),typeof(offsets),typeof(weights)}(
        image, model.pixel_scale_arcsec, offsets, weights)
end

freeze_source(src::ExtendedSource) = with_extended_source(src.source, src.model)

function extended_source_asterism(src::ExtendedSource)
    return src.quadrature
end
