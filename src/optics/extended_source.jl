abstract type AbstractExtendedSourceModel end

struct PointCloudSourceModel{T<:AbstractFloat,V<:AbstractVector{NTuple{2,T}},W<:AbstractVector{T}} <: AbstractExtendedSourceModel
    offsets_xy_arcsec::V
    weights::W
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
        w >= zero(T) || throw(InvalidConfiguration("source-model sample weights must be >= 0"))
        normalized_offsets[i] = (T(offset[1]), T(offset[2]))
        normalized_weights[i] = w
        total += w
    end
    total > zero(T) || throw(InvalidConfiguration("source-model sample weights must sum to a positive value"))
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

struct ExtendedSource{S<:Union{Source,LGSSource},M<:AbstractExtendedSourceModel} <: AbstractSource
    source::S
    model::M
end

function with_extended_source(src::Union{Source,LGSSource}, model::AbstractExtendedSourceModel)
    return ExtendedSource(src, model)
end

function with_extended_source(src::AbstractSource, ::AbstractExtendedSourceModel)
    throw(UnsupportedAlgorithm("extended-source models currently support Source and LGSSource only"))
end

wavelength(src::ExtendedSource) = wavelength(src.source)
photon_flux(src::ExtendedSource) = photon_flux(src.source)
coordinates_xy_arcsec(src::ExtendedSource) = coordinates_xy_arcsec(src.source)
source_height_m(src::ExtendedSource) = source_height_m(src.source)
optical_tag(src::ExtendedSource) = "extended(" * optical_tag(src.source) * ")"
is_lgs_source(src::ExtendedSource) = is_lgs_source(src.source)

has_extended_source_model(::AbstractSource) = false
has_extended_source_model(::ExtendedSource) = true

extended_source_model(::AbstractSource) = nothing
extended_source_model(src::ExtendedSource) = src.model

extended_source_reference_source(src::AbstractSource) = src
extended_source_reference_source(src::ExtendedSource) = src.source

function source_measurement_signature(src::ExtendedSource)
    return hash((source_measurement_signature(src.source), source_model_signature(src.model)))
end

function source_with_coordinates_and_flux(src::Source, coords_xy_arcsec::NTuple{2,T}, flux::T) where {T<:AbstractFloat}
    params = src.params
    return Source(SourceParams{T}(params.band, T(params.magnitude), coords_xy_arcsec, T(params.wavelength), flux))
end

function source_with_coordinates_and_flux(src::LGSSource, coords_xy_arcsec::NTuple{2,T}, flux::T) where {T<:AbstractFloat}
    params = src.params
    return LGSSource(LGSSourceParams{T, typeof(params.na_profile)}(
        T(params.magnitude),
        coords_xy_arcsec,
        T(params.wavelength),
        T(params.altitude),
        T(params.elongation_factor),
        (T(params.laser_coordinates[1]), T(params.laser_coordinates[2])),
        params.na_profile,
        T(params.fwhm_spot_up),
        flux,
    ))
end

function extended_source_asterism(src::ExtendedSource)
    base_coords = coordinates_xy_arcsec(src.source)
    total_flux = photon_flux(src.source)
    offsets = source_model_offsets(src.model)
    weights = source_model_weights(src.model)
    samples = Vector{typeof(source_with_coordinates_and_flux(src.source, base_coords, total_flux * eltype(weights)(weights[1])))}(undef, length(weights))
    @inbounds for i in eachindex(offsets, weights)
        offset = offsets[i]
        coords = (base_coords[1] + offset[1], base_coords[2] + offset[2])
        samples[i] = source_with_coordinates_and_flux(src.source, coords, total_flux * weights[i])
    end
    return Asterism(samples)
end
