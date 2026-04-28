"""
Internal marker for atmosphere layers that participate in the shared finite /
infinite rendering path.

Maintained layer implementations are expected to support:

- `sample_layer!(out, layer, tel, rng)` for evolving the layer and writing a
  pupil-sized sample into `out`
- `sample_layer_accumulate!(out, layer, tel, rng)` for evolving the layer and
  adding a pupil-sized sample into `out`
- `render_layer!(out, layer, shift_x, shift_y, footprint_scale)` for writing a
  source-aware pupil sample into `out` without evolving the layer state
- `render_layer_accumulate!(out, layer, shift_x, shift_y, footprint_scale)` for
  adding a source-aware pupil sample into `out` without evolving the layer
  state
- `layer_altitude(layer)` returning the layer altitude in meters

The shared atmosphere-container helpers in this file assume their `atm`
argument provides:

- `atm.layers`
- `atm.params.altitude`
- `atm.state.opd`
- `atm.state.source_geometry`
"""
abstract type AbstractAtmosphereLayer end

@inline layer_altitude(layer::AbstractAtmosphereLayer) = layer.params.altitude

mutable struct AtmosphereSourceGeometryCache{T<:AbstractFloat,V<:AbstractVector{T}}
    shift_x::V
    shift_y::V
    footprint_scale::V
    cached_x_arcsec::T
    cached_y_arcsec::T
    cached_height_m::T
    valid::Bool
end

function AtmosphereSourceGeometryCache(n_layers::Int, ::Type{T}=Float64) where {T<:AbstractFloat}
    return AtmosphereSourceGeometryCache(
        zeros(T, n_layers),
        zeros(T, n_layers),
        ones(T, n_layers),
        T(NaN),
        T(NaN),
        T(NaN),
        false,
    )
end

@inline function source_geometry_signature(::Nothing, ::Type{T}) where {T<:AbstractFloat}
    return zero(T), zero(T), T(Inf)
end

function source_geometry_signature(src::AbstractSource, ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec = coordinates_xy_arcsec(src)
    height_m = source_height_m(src)
    return T(x_arcsec), T(y_arcsec), isfinite(height_m) ? T(height_m) : T(Inf)
end

@inline function is_onaxis_infinite_source(src::AbstractSource, ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec, height_m = source_geometry_signature(src, T)
    return x_arcsec == zero(T) && y_arcsec == zero(T) && !isfinite(height_m)
end

@inline function is_onaxis_infinite_source(::Nothing, ::Type{T}) where {T<:AbstractFloat}
    return true
end

@inline function layer_source_geometry(::Nothing, layer_altitude::Real, tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    return zero(T), zero(T), one(T)
end

function layer_source_geometry(src::AbstractSource, layer_altitude::Real, tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec, height_m = source_geometry_signature(src, T)
    delta = T(tel.params.diameter / tel.params.resolution)
    altitude_t = T(layer_altitude)
    arcsec_to_rad = T(pi / (180 * 3600))
    shift_x = x_arcsec * arcsec_to_rad * altitude_t / delta
    shift_y = y_arcsec * arcsec_to_rad * altitude_t / delta

    footprint_scale = one(T)
    if isfinite(height_m)
        height_m > altitude_t || throw(InvalidConfiguration("source height must exceed layer altitude for finite-height propagation"))
        footprint_scale = (height_m - altitude_t) / height_m
    end
    return shift_x, shift_y, footprint_scale
end

function ensure_source_geometry_cache!(cache::AtmosphereSourceGeometryCache{T},
    src::AbstractSource, altitudes::AbstractVector, tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec, height_m = source_geometry_signature(src, T)
    if !cache.valid || cache.cached_x_arcsec != x_arcsec || cache.cached_y_arcsec != y_arcsec || cache.cached_height_m != height_m
        delta = T(tel.params.diameter / tel.params.resolution)
        arcsec_to_rad = T(pi / (180 * 3600))
        @inbounds for i in eachindex(altitudes, cache.shift_x, cache.shift_y, cache.footprint_scale)
            altitude_t = T(altitudes[i])
            cache.shift_x[i] = x_arcsec * arcsec_to_rad * altitude_t / delta
            cache.shift_y[i] = y_arcsec * arcsec_to_rad * altitude_t / delta
            if isfinite(height_m)
                height_m > altitude_t || throw(InvalidConfiguration("source height must exceed layer altitude for finite-height propagation"))
                cache.footprint_scale[i] = (height_m - altitude_t) / height_m
            else
                cache.footprint_scale[i] = one(T)
            end
        end
        cache.cached_x_arcsec = x_arcsec
        cache.cached_y_arcsec = y_arcsec
        cache.cached_height_m = height_m
        cache.valid = true
    end
    return cache
end

function accumulate_sampled_layers!(opd::AbstractMatrix, layers, tel::Telescope, rng::AbstractRNG)
    isempty(layers) && return opd
    sample_layer!(opd, layers[1], tel, rng)
    @inbounds for i in 2:length(layers)
        sample_layer_accumulate!(opd, layers[i], tel, rng)
    end
    return opd
end

function accumulate_rendered_layers!(opd::AbstractMatrix,
    layers, shift_x::AbstractVector, shift_y::AbstractVector, footprint_scale::AbstractVector)
    isempty(layers) && return opd
    render_layer!(opd, layers[1], shift_x[1], shift_y[1], footprint_scale[1])
    @inbounds for i in 2:length(layers)
        render_layer_accumulate!(opd, layers[i], shift_x[i], shift_y[i], footprint_scale[i])
    end
    return opd
end

function propagate_source_aware!(atm, tel::Telescope, src::AbstractSource)
    T = eltype(atm.state.opd)
    if is_onaxis_infinite_source(src, T)
        return propagate!(atm, tel)
    end
    cache = ensure_source_geometry_cache!(atm.state.source_geometry, src, atm.params.altitude, tel, T)
    accumulate_rendered_layers!(atm.state.opd, atm.layers,
        cache.shift_x, cache.shift_y, cache.footprint_scale)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
