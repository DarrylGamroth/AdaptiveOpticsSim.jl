abstract type AbstractAtmosphericFieldModel end

struct GeometricAtmosphericPropagation{T<:AbstractFloat} <: AbstractAtmosphericFieldModel
    chromatic_reference_wavelength::Union{Nothing,T}
end

struct LayeredFresnelAtmosphericPropagation{T<:AbstractFloat} <: AbstractAtmosphericFieldModel
    band_limit_factor::T
    chromatic_reference_wavelength::Union{Nothing,T}
end

struct AtmosphericFieldPropagationParams{T<:AbstractFloat,M<:AbstractAtmosphericFieldModel,V<:AbstractVector{Int},D<:AbstractVector{T}}
    zero_padding::Int
    model::M
    layer_order::V
    layer_distances_m::D
end

mutable struct AtmosphericFieldSlice{
    T<:AbstractFloat,
    S<:AbstractSource,
    F<:ElectricField,
    P<:AbstractMatrix{T},
    V,
}
    source::S
    field::F
    phase_buffer::P
    propagators::V
end

mutable struct AtmosphericFieldPropagationState{
    T<:AbstractFloat,
    S<:AbstractVector,
    R<:AbstractMatrix{T},
}
    slices::S
    intensity::R
end

struct AtmosphericFieldPropagation{P<:AtmosphericFieldPropagationParams,S<:AtmosphericFieldPropagationState}
    params::P
    state::S
end

function GeometricAtmosphericPropagation(; chromatic_reference_wavelength=nothing, T::Type{<:AbstractFloat}=Float64)
    ref = isnothing(chromatic_reference_wavelength) ? nothing : T(chromatic_reference_wavelength)
    return GeometricAtmosphericPropagation{T}(ref)
end

function LayeredFresnelAtmosphericPropagation(; band_limit_factor::Real=1.0,
    chromatic_reference_wavelength=nothing,
    T::Type{<:AbstractFloat}=Float64)
    zero(T) <= band_limit_factor <= one(T) ||
        throw(InvalidConfiguration("band_limit_factor must be between 0 and 1"))
    ref = isnothing(chromatic_reference_wavelength) ? nothing : T(chromatic_reference_wavelength)
    return LayeredFresnelAtmosphericPropagation{T}(T(band_limit_factor), ref)
end

@inline atmospheric_model_wavelength(model::AbstractAtmosphericFieldModel) = model.chromatic_reference_wavelength

function _layer_order_and_distances(atm, ::Type{T}) where {T<:AbstractFloat}
    order = sortperm(collect(1:length(atm.layers)); by=i -> atm.params.altitude[i], rev=true)
    distances = Vector{T}(undef, length(order))
    @inbounds for k in eachindex(order)
        idx = order[k]
        next_alt = k == length(order) ? zero(T) : T(atm.params.altitude[order[k + 1]])
        distances[k] = max(zero(T), T(atm.params.altitude[idx]) - next_alt)
    end
    return order, distances
end

function _build_field_slice(atm, tel::Telescope, src::AbstractSource, zero_padding::Int,
    ::Type{T}, model::GeometricAtmosphericPropagation{T}) where {T<:AbstractFloat}
    field = ElectricField(tel, src; zero_padding=zero_padding, T=T)
    phase_buffer = similar(tel.state.opd, T, tel.params.resolution, tel.params.resolution)
    fill!(phase_buffer, zero(T))
    return AtmosphericFieldSlice{T, typeof(src), typeof(field), typeof(phase_buffer), Nothing}(src, field, phase_buffer, nothing)
end

function _build_field_slice(atm, tel::Telescope, src::AbstractSource, zero_padding::Int,
    ::Type{T}, model::LayeredFresnelAtmosphericPropagation{T}) where {T<:AbstractFloat}
    field = ElectricField(tel, src; zero_padding=zero_padding, T=T)
    phase_buffer = similar(tel.state.opd, T, tel.params.resolution, tel.params.resolution)
    fill!(phase_buffer, zero(T))
    _, distances = _layer_order_and_distances(atm, T)
    propagators = map(distances) do dist
        dist == zero(T) ? nothing : FresnelPropagation(field; distance_m=dist)
    end
    return AtmosphericFieldSlice{T, typeof(src), typeof(field), typeof(phase_buffer), typeof(propagators)}(src, field, phase_buffer, propagators)
end

function _build_slices(atm, tel::Telescope, src::AbstractSource, zero_padding::Int,
    ::Type{T}, model::AbstractAtmosphericFieldModel) where {T<:AbstractFloat}
    slices = Vector{Any}(undef, 1)
    slices[1] = _build_field_slice(atm, tel, src, zero_padding, T, model)
    return slices
end

function _build_slices(atm, tel::Telescope, src::SpectralSource, zero_padding::Int,
    ::Type{T}, model::AbstractAtmosphericFieldModel) where {T<:AbstractFloat}
    bundle = spectral_bundle(src)
    slices = Vector{Any}(undef, length(bundle))
    @inbounds for i in eachindex(bundle.samples)
        sample = bundle.samples[i]
        sample_src = source_with_wavelength_and_flux(src, sample.wavelength, T(photon_flux(src) * sample.weight))
        slices[i] = _build_field_slice(atm, tel, sample_src, zero_padding, T, model)
    end
    return slices
end

function AtmosphericFieldPropagation(atm, tel::Telescope, src;
    model::AbstractAtmosphericFieldModel=GeometricAtmosphericPropagation(T=eltype(tel.state.opd)),
    zero_padding::Int=1,
    T::Type{<:AbstractFloat}=eltype(tel.state.opd))
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    order, distances = _layer_order_and_distances(atm, T)
    slices = _build_slices(atm, tel, src, zero_padding, T, model)
    intensity = similar(first(slices).field.state.intensity)
    fill!(intensity, zero(T))
    params = AtmosphericFieldPropagationParams{T, typeof(model), typeof(order), typeof(distances)}(
        zero_padding,
        model,
        order,
        distances,
    )
    state = AtmosphericFieldPropagationState{T, typeof(slices), typeof(intensity)}(slices, intensity)
    return AtmosphericFieldPropagation(params, state)
end

@inline function _chromatic_shift_scale(model::AbstractAtmosphericFieldModel, λ::T) where {T<:AbstractFloat}
    ref = atmospheric_model_wavelength(model)
    return isnothing(ref) ? one(T) : T(ref / λ)
end

function _render_layer_phase!(out::AbstractMatrix{T}, layer::AbstractAtmosphereLayer, tel::Telescope,
    src::AbstractSource, model::AbstractAtmosphericFieldModel, λ::T) where {T<:AbstractFloat}
    shift_x, shift_y, footprint_scale = layer_source_geometry(src, layer_altitude(layer), tel, T)
    shift_scale = _chromatic_shift_scale(model, λ)
    return render_layer!(out, layer, shift_x * shift_scale, shift_y * shift_scale, footprint_scale)
end

function _propagate_slice_geometric!(slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm, tel::Telescope) where {T<:AbstractFloat}
    fill_from_telescope!(slice.field, tel, slice.source)
    @inbounds for idx in prop.params.layer_order
        layer = atm.layers[idx]
        _render_layer_phase!(slice.phase_buffer, layer, tel, slice.source, prop.params.model, slice.field.params.wavelength)
        apply_phase!(slice.field, slice.phase_buffer; units=:opd)
    end
    return slice.field
end

function _propagate_slice_fresnel!(slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm, tel::Telescope) where {T<:AbstractFloat}
    fill_from_telescope!(slice.field, tel, slice.source)
    @inbounds for k in eachindex(prop.params.layer_order)
        idx = prop.params.layer_order[k]
        layer = atm.layers[idx]
        _render_layer_phase!(slice.phase_buffer, layer, tel, slice.source, prop.params.model, slice.field.params.wavelength)
        apply_phase!(slice.field, slice.phase_buffer; units=:opd)
        local propagator = slice.propagators[k]
        isnothing(propagator) || propagate_field!(slice.field, propagator)
    end
    return slice.field
end

function propagate_atmosphere_field!(prop::AtmosphericFieldPropagation, atm, tel::Telescope, src::AbstractSource)
    length(prop.state.slices) == 1 || throw(InvalidConfiguration("monochromatic atmosphere propagation requires a monochromatic workspace"))
    slice = prop.state.slices[1]
    if prop.params.model isa GeometricAtmosphericPropagation
        return _propagate_slice_geometric!(slice, prop, atm, tel)
    else
        return _propagate_slice_fresnel!(slice, prop, atm, tel)
    end
end

function propagate_atmosphere_field!(field::ElectricField, prop::AtmosphericFieldPropagation, atm, tel::Telescope, src::AbstractSource)
    propagated = propagate_atmosphere_field!(prop, atm, tel, src)
    size(field.state.field) == size(propagated.state.field) ||
        throw(DimensionMismatchError("ElectricField size must match propagated atmosphere field size"))
    copyto!(field.state.field, propagated.state.field)
    return field
end

function atmospheric_intensity!(out::AbstractMatrix{T}, prop::AtmosphericFieldPropagation, atm, tel::Telescope,
    src::AbstractSource) where {T<:AbstractFloat}
    propagated = propagate_atmosphere_field!(prop, atm, tel, src)
    size(out) == size(propagated.state.field) ||
        throw(DimensionMismatchError("atmospheric intensity output must match propagated field size"))
    @. out = abs2(propagated.state.field)
    return out
end

function atmospheric_intensity!(out::AbstractMatrix{T}, prop::AtmosphericFieldPropagation, atm, tel::Telescope,
    src::SpectralSource) where {T<:AbstractFloat}
    size(out) == size(prop.state.intensity) ||
        throw(DimensionMismatchError("spectral atmospheric intensity output must match propagation workspace size"))
    fill!(out, zero(T))
    @inbounds for slice_any in prop.state.slices
        slice = slice_any
        if prop.params.model isa GeometricAtmosphericPropagation
            _propagate_slice_geometric!(slice, prop, atm, tel)
        else
            _propagate_slice_fresnel!(slice, prop, atm, tel)
        end
        out .+= abs2.(slice.field.state.field)
    end
    return out
end

function atmospheric_intensity!(prop::AtmosphericFieldPropagation, atm, tel::Telescope, src)
    atmospheric_intensity!(prop.state.intensity, prop, atm, tel, src)
    return prop.state.intensity
end
