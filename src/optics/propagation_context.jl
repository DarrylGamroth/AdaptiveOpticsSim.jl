struct LayerRenderContext{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    footprint_scale::T
end

@inline function layer_render_context(src::Union{AbstractSource,Nothing}, layer_altitude::Real, tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    shift_x, shift_y, footprint_scale = layer_source_geometry(src, layer_altitude, tel, T)
    return LayerRenderContext{T}(shift_x, shift_y, footprint_scale)
end

@inline layer_render_context(src::Union{AbstractSource,Nothing}, layer::AbstractAtmosphereLayer, tel::Telescope, ::Type{T}) where {T<:AbstractFloat} =
    layer_render_context(src, layer_altitude(layer), tel, T)

@inline render_layer!(out::AbstractMatrix{T}, layer::AbstractAtmosphereLayer, ctx::LayerRenderContext{T}) where {T<:AbstractFloat} =
    render_layer!(out, layer, ctx.shift_x, ctx.shift_y, ctx.footprint_scale)

@inline render_layer_accumulate!(out::AbstractMatrix{T}, layer::AbstractAtmosphereLayer, ctx::LayerRenderContext{T}) where {T<:AbstractFloat} =
    render_layer_accumulate!(out, layer, ctx.shift_x, ctx.shift_y, ctx.footprint_scale)
