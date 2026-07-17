abstract type AbstractAtmosphericFieldModel end
abstract type AbstractAtmosphericFieldExecutionPlan end

struct GeometricFieldSynchronousPlan <: AbstractAtmosphericFieldExecutionPlan end
struct GeometricFieldAsyncPlan <: AbstractAtmosphericFieldExecutionPlan end
struct LayeredFresnelFieldSynchronousPlan <: AbstractAtmosphericFieldExecutionPlan end
struct LayeredFresnelFieldAsyncPlan <: AbstractAtmosphericFieldExecutionPlan end

struct GeometricAtmosphericPropagation{T<:AbstractFloat} <: AbstractAtmosphericFieldModel
    chromatic_reference_wavelength::Union{Nothing,T}
end

struct LayeredFresnelAtmosphericPropagation{T<:AbstractFloat} <: AbstractAtmosphericFieldModel
    band_limit_factor::T
    chromatic_reference_wavelength::Union{Nothing,T}
end

struct AtmosphericFieldPropagationParams{T<:AbstractFloat,M<:AbstractAtmosphericFieldModel,V<:AbstractVector{Int},D<:AbstractVector{T},I<:AtmosphereIdentity}
    zero_padding::Int
    model::M
    layer_order::V
    layer_distances_m::D
    atmosphere_identity::I
end

struct AtmosphericFieldSlice{
    T<:AbstractFloat,
    S<:AbstractSource,
    W<:PupilFunction,
    F<:ElectricField,
    FP<:PupilFieldFormationPlan,
    P<:AbstractMatrix{T},
    V,
    C<:AbstractVector{LayerRenderContext{T}},
}
    source::S
    wavefront::W
    field::F
    formation_plan::FP
    phase_buffer::P
    propagators::V
    contexts::C
end

struct AtmosphericFieldPropagationState{
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

@inline atmospheric_field_execution_plan(style::ExecutionStyle, model::AbstractAtmosphericFieldModel) =
    atmospheric_field_execution_plan(typeof(style), typeof(model))

@inline atmospheric_field_execution_plan(::Type{<:ScalarCPUStyle}, ::Type{<:GeometricAtmosphericPropagation}) =
    GeometricFieldSynchronousPlan()
@inline atmospheric_field_execution_plan(::Type{<:AcceleratorStyle}, ::Type{<:GeometricAtmosphericPropagation}) =
    GeometricFieldAsyncPlan()
@inline atmospheric_field_execution_plan(::Type{<:ScalarCPUStyle}, ::Type{<:LayeredFresnelAtmosphericPropagation}) =
    LayeredFresnelFieldSynchronousPlan()
@inline atmospheric_field_execution_plan(::Type{<:AcceleratorStyle}, ::Type{<:LayeredFresnelAtmosphericPropagation}) =
    LayeredFresnelFieldAsyncPlan()

function GeometricAtmosphericPropagation(; chromatic_reference_wavelength=nothing, T::Type{<:AbstractFloat}=Float64)
    ref = isnothing(chromatic_reference_wavelength) ? nothing :
        _converted_positive_finite(chromatic_reference_wavelength, T,
            "chromatic reference wavelength")
    return GeometricAtmosphericPropagation{T}(ref)
end

function LayeredFresnelAtmosphericPropagation(; band_limit_factor::Real=1.0,
    chromatic_reference_wavelength=nothing,
    T::Type{<:AbstractFloat}=Float64)
    factor = T(band_limit_factor)
    isfinite(factor) && zero(T) <= factor <= one(T) ||
        throw(InvalidConfiguration("band_limit_factor must be between 0 and 1"))
    ref = isnothing(chromatic_reference_wavelength) ? nothing :
        _converted_positive_finite(chromatic_reference_wavelength, T,
            "chromatic reference wavelength")
    return LayeredFresnelAtmosphericPropagation{T}(factor, ref)
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
    λ = _converted_positive_finite(wavelength(src), T,
        "atmospheric source wavelength")
    wavefront = PupilFunction(tel; T=T)
    copyto!(wavefront.opd, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=zero_padding, T=T)
    formation_plan = prepare_pupil_field(tel, wavefront, src, field)
    phase_buffer = similar(wavefront.opd)
    fill!(phase_buffer, zero(T))
    contexts = _build_layer_contexts(atm, tel, src, model, λ, T)
    return AtmosphericFieldSlice{
        T,typeof(src),typeof(wavefront),typeof(field),
        typeof(formation_plan),typeof(phase_buffer),Nothing,typeof(contexts),
    }(src, wavefront, field, formation_plan, phase_buffer, nothing, contexts)
end

function _build_field_slice(atm, tel::Telescope, src::AbstractSource, zero_padding::Int,
    ::Type{T}, model::LayeredFresnelAtmosphericPropagation{T}) where {T<:AbstractFloat}
    λ = _converted_positive_finite(wavelength(src), T,
        "atmospheric source wavelength")
    wavefront = PupilFunction(tel; T=T)
    copyto!(wavefront.opd, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=zero_padding, T=T)
    formation_plan = prepare_pupil_field(tel, wavefront, src, field)
    phase_buffer = similar(wavefront.opd)
    fill!(phase_buffer, zero(T))
    _, distances = _layer_order_and_distances(atm, T)
    propagators = map(distances) do dist
        dist == zero(T) ? nothing : FresnelPropagation(field;
            distance_m=dist, output_kind=field.metadata.kind)
    end
    contexts = _build_layer_contexts(atm, tel, src, model, λ, T)
    return AtmosphericFieldSlice{
        T,typeof(src),typeof(wavefront),typeof(field),
        typeof(formation_plan),typeof(phase_buffer),typeof(propagators),
        typeof(contexts),
    }(src, wavefront, field, formation_plan, phase_buffer, propagators, contexts)
end

function _build_slices(atm, tel::Telescope, src::AbstractSource, zero_padding::Int,
    ::Type{T}, model::AbstractAtmosphericFieldModel) where {T<:AbstractFloat}
    slice = _build_field_slice(atm, tel, src, zero_padding, T, model)
    slices = Vector{typeof(slice)}(undef, 1)
    slices[1] = slice
    return slices
end

function _build_slices(atm, tel::Telescope, src::SpectralSource, zero_padding::Int,
    ::Type{T}, model::AbstractAtmosphericFieldModel) where {T<:AbstractFloat}
    bundle = spectral_bundle(src)
    total_value = _converted_nonnegative_finite(
        source_radiometric_value(src), T,
        "spectral source radiometric value")
    first_sample = bundle.samples[1]
    first_wavelength = _converted_positive_finite(first_sample.wavelength,
        T, "spectral source wavelength")
    first_value = _converted_nonnegative_finite(
        total_value * first_sample.weight, T,
        "spectral source radiometric value")
    first_src = source_with_wavelength_and_radiometric_value(src,
        first_wavelength, first_value)
    first_slice = _build_field_slice(atm, tel, first_src, zero_padding, T, model)
    slices = Vector{typeof(first_slice)}(undef, length(bundle))
    slices[1] = first_slice
    @inbounds for i in 2:length(bundle)
        sample = bundle.samples[i]
        sample_wavelength = _converted_positive_finite(sample.wavelength,
            T, "spectral source wavelength")
        sample_value = _converted_nonnegative_finite(
            total_value * sample.weight, T,
            "spectral source radiometric value")
        sample_src = source_with_wavelength_and_radiometric_value(src,
            sample_wavelength, sample_value)
        slices[i] = _build_field_slice(atm, tel, sample_src, zero_padding, T, model)
    end
    return slices
end

function AtmosphericFieldPropagation(atm::AbstractTimedAtmosphere,
    tel::Telescope, src;
    model::AbstractAtmosphericFieldModel=GeometricAtmosphericPropagation(T=eltype(tel.state.opd)),
    zero_padding::Int=1,
    T::Type{<:AbstractFloat}=eltype(tel.state.opd))
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    require_same_backend(atm, tel)
    atmosphere_numeric_type(atm) === T || throw(InvalidConfiguration(
        "atmospheric field numeric type must match atmosphere layer storage"))
    frozen_source = freeze_source(src)
    order, distances = _layer_order_and_distances(atm, T)
    slices = _build_slices(atm, tel, frozen_source, zero_padding, T, model)
    intensity = similar(first(slices).field.values, T,
        size(first(slices).field.values)...)
    fill!(intensity, zero(T))
    identity = atmosphere_identity(atm)
    params = AtmosphericFieldPropagationParams{T, typeof(model), typeof(order),
        typeof(distances),typeof(identity)}(
        zero_padding,
        model,
        order,
        distances,
        identity,
    )
    state = AtmosphericFieldPropagationState{T, typeof(slices), typeof(intensity)}(slices, intensity)
    return AtmosphericFieldPropagation(params, state)
end

@inline function _chromatic_shift_scale(model::AbstractAtmosphericFieldModel, λ::T) where {T<:AbstractFloat}
    ref = atmospheric_model_wavelength(model)
    return isnothing(ref) ? one(T) : T(ref / λ)
end

function _build_layer_contexts(atm, tel::Telescope, src::AbstractSource,
    model::AbstractAtmosphericFieldModel, λ::Real,
    ::Type{T}) where {T<:AbstractFloat}
    shift_scale = _chromatic_shift_scale(model, T(λ))
    contexts = Vector{LayerRenderContext{T}}(undef, length(atm.layers))
    @inbounds for i in eachindex(atm.layers)
        ctx = layer_render_context(src, atm.layers[i], tel, T)
        contexts[i] = LayerRenderContext{T}(
            ctx.shift_x * shift_scale,
            ctx.shift_y * shift_scale,
            ctx.footprint_scale,
        )
    end
    return contexts
end

function _propagate_slice!(::GeometricFieldSynchronousPlan,
    slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    style = execution_style(slice.field.values)
    fill_electric_field_async!(slice.field, slice.wavefront,
        slice.formation_plan)
    @inbounds for idx in prop.params.layer_order
        layer = atm.layers[idx]
        render_layer!(slice.phase_buffer, layer, slice.contexts[idx])
        apply_phase_async!(slice.field, slice.phase_buffer,
            slice.formation_plan; units=:opd)
    end
    synchronize_backend!(style)
    return slice.field
end

function _propagate_slice!(::GeometricFieldAsyncPlan,
    slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    fill_electric_field_async!(slice.field, slice.wavefront,
        slice.formation_plan)
    @inbounds for idx in prop.params.layer_order
        layer = atm.layers[idx]
        render_layer!(slice.phase_buffer, layer, slice.contexts[idx])
        apply_phase_async!(slice.field, slice.phase_buffer,
            slice.formation_plan; units=:opd)
    end
    return slice.field
end

function _propagate_slice!(::LayeredFresnelFieldSynchronousPlan,
    slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    style = execution_style(slice.field.values)
    fill_electric_field_async!(slice.field, slice.wavefront,
        slice.formation_plan)
    @inbounds for k in eachindex(prop.params.layer_order)
        idx = prop.params.layer_order[k]
        layer = atm.layers[idx]
        render_layer!(slice.phase_buffer, layer, slice.contexts[idx])
        apply_phase_async!(slice.field, slice.phase_buffer,
            slice.formation_plan; units=:opd)
        synchronize_backend!(style)
        local propagator = slice.propagators[k]
        isnothing(propagator) || propagate_field!(slice.field.values,
            slice.field, propagator)
    end
    return slice.field
end

function _propagate_slice!(::LayeredFresnelFieldAsyncPlan,
    slice::AtmosphericFieldSlice{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    style = execution_style(slice.field.values)
    fill_electric_field_async!(slice.field, slice.wavefront,
        slice.formation_plan)
    @inbounds for k in eachindex(prop.params.layer_order)
        idx = prop.params.layer_order[k]
        layer = atm.layers[idx]
        render_layer!(slice.phase_buffer, layer, slice.contexts[idx])
        apply_phase_async!(slice.field, slice.phase_buffer,
            slice.formation_plan; units=:opd)
        synchronize_backend!(style)
        local propagator = slice.propagators[k]
        isnothing(propagator) || propagate_field!(slice.field.values,
            slice.field, propagator)
    end
    return slice.field
end

function _validate_field_epoch(prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    _validate_epoch_identity(prop.params.atmosphere_identity, atm, epoch)
    n_layers = length(atm.layers)
    length(prop.params.layer_order) == n_layers ||
        throw(AtmosphereEpochError(
            "atmosphere layer shape changed after field renderer preparation"))
    @inbounds for slice in prop.state.slices
        length(slice.contexts) == n_layers || throw(AtmosphereEpochError(
            "atmosphere layer shape changed after field renderer preparation"))
    end
    return epoch
end

function propagate_atmosphere_field!(prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    _validate_field_epoch(prop, atm, epoch)
    length(prop.state.slices) == 1 || throw(InvalidConfiguration(
        "monochromatic atmosphere propagation requires a monochromatic workspace"))
    slice = prop.state.slices[1]
    plan = atmospheric_field_execution_plan(execution_style(slice.field.values),
        prop.params.model)
    return _propagate_slice!(plan, slice, prop, atm)
end

function propagate_atmosphere_field!(field::ElectricField,
    prop::AtmosphericFieldPropagation, atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch)
    _validate_field_epoch(prop, atm, epoch)
    length(prop.state.slices) == 1 || throw(InvalidConfiguration(
        "monochromatic atmosphere propagation requires a monochromatic workspace"))
    prepared = prop.state.slices[1].field
    size(field.values) == size(prepared.values) ||
        throw(DimensionMismatchError("ElectricField size must match propagated atmosphere field size"))
    field.metadata == prepared.metadata || throw(InvalidConfiguration(
        "ElectricField metadata must match propagated atmosphere field metadata"))
    validate_plane_storage(field.metadata, field.values;
        label="atmospheric ElectricField destination")
    propagated = propagate_atmosphere_field!(prop, atm, epoch)
    copyto!(field.values, propagated.values)
    return field
end

function _validate_atmospheric_intensity_destination(out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation) where {T<:AbstractFloat}
    prepared = prop.state.intensity
    size(out) == size(prepared) ||
        throw(DimensionMismatchError("atmospheric intensity output must match propagated field size"))
    eltype(out) === eltype(prepared) || throw(InvalidConfiguration(
        "atmospheric intensity numeric type does not match prepared output"))
    typeof(backend(out)) === typeof(backend(prepared)) ||
        throw(InvalidConfiguration(
            "atmospheric intensity backend does not match prepared output"))
    plane_device(out) == plane_device(prepared) || throw(InvalidConfiguration(
        "atmospheric intensity destination is on a different physical device"))
    return out
end

function atmospheric_intensity!(out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation, atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch) where {T<:AbstractFloat}
    _validate_field_epoch(prop, atm, epoch)
    _validate_atmospheric_intensity_destination(out, prop)
    if length(prop.state.slices) == 1
        slice = prop.state.slices[1]
        plan = atmospheric_field_execution_plan(execution_style(slice.field.values),
            prop.params.model)
        _propagate_slice!(plan, slice, prop, atm)
        intensity!(out, slice.field)
        return out
    end
    return _spectral_atmospheric_intensity!(execution_style(out), out, prop,
        atm)
end

function _spectral_atmospheric_intensity!(::ScalarCPUStyle,
    out::AbstractMatrix{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    fill!(out, zero(T))
    plan = atmospheric_field_execution_plan(ScalarCPUStyle(), prop.params.model)
    _spectral_atmospheric_intensity_model!(plan, out, prop, atm)
    return out
end

function _spectral_atmospheric_intensity!(style::AcceleratorStyle,
    out::AbstractMatrix{T}, prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    fill!(out, zero(T))
    plan = atmospheric_field_execution_plan(style, prop.params.model)
    _spectral_atmospheric_intensity_model!(plan, out, prop, atm)
    synchronize_backend!(style)
    return out
end

function _spectral_atmospheric_intensity_model!(::GeometricFieldSynchronousPlan, out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation{<:AtmosphericFieldPropagationParams{T,<:GeometricAtmosphericPropagation}},
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    plan = GeometricFieldSynchronousPlan()
    @inbounds for slice in prop.state.slices
        _propagate_slice!(plan, slice, prop, atm)
        _accumulate_field_intensity!(out, slice.field)
    end
    return out
end

function _spectral_atmospheric_intensity_model!(::GeometricFieldAsyncPlan, out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation{<:AtmosphericFieldPropagationParams{T,<:GeometricAtmosphericPropagation}},
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    plan = GeometricFieldAsyncPlan()
    @inbounds for slice in prop.state.slices
        _propagate_slice!(plan, slice, prop, atm)
        _accumulate_field_intensity_async!(out, slice.field)
    end
    return out
end

function _spectral_atmospheric_intensity_model!(::LayeredFresnelFieldSynchronousPlan, out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation{<:AtmosphericFieldPropagationParams{T,<:LayeredFresnelAtmosphericPropagation}},
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    plan = LayeredFresnelFieldSynchronousPlan()
    @inbounds for slice in prop.state.slices
        _propagate_slice!(plan, slice, prop, atm)
        _accumulate_field_intensity!(out, slice.field)
    end
    return out
end

function _spectral_atmospheric_intensity_model!(::LayeredFresnelFieldAsyncPlan, out::AbstractMatrix{T},
    prop::AtmosphericFieldPropagation{<:AtmosphericFieldPropagationParams{T,<:LayeredFresnelAtmosphericPropagation}},
    atm::AbstractTimedAtmosphere) where {T<:AbstractFloat}
    plan = LayeredFresnelFieldAsyncPlan()
    @inbounds for slice in prop.state.slices
        _propagate_slice!(plan, slice, prop, atm)
        _accumulate_field_intensity_async!(out, slice.field)
    end
    return out
end

function atmospheric_intensity!(prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    atmospheric_intensity!(prop.state.intensity, prop, atm, epoch)
    return prop.state.intensity
end

# Transitional call shapes for optical consumers that have not yet migrated to
# the explicit epoch executor. Preparation has already frozen telescope and
# source state, so these arguments are intentionally not read on execution.
@inline propagate_atmosphere_field!(prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere, ::Telescope, ::AbstractSource) =
    propagate_atmosphere_field!(prop, atm, current_epoch(atm))

@inline propagate_atmosphere_field!(field::ElectricField,
    prop::AtmosphericFieldPropagation, atm::AbstractTimedAtmosphere,
    ::Telescope, ::AbstractSource) =
    propagate_atmosphere_field!(field, prop, atm, current_epoch(atm))

@inline atmospheric_intensity!(out::AbstractMatrix,
    prop::AtmosphericFieldPropagation, atm::AbstractTimedAtmosphere,
    ::Telescope, ::AbstractSource) =
    atmospheric_intensity!(out, prop, atm, current_epoch(atm))

@inline atmospheric_intensity!(prop::AtmosphericFieldPropagation,
    atm::AbstractTimedAtmosphere, ::Telescope, ::AbstractSource) =
    atmospheric_intensity!(prop, atm, current_epoch(atm))
