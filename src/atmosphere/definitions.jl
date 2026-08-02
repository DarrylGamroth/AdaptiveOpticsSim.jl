"""Cold scientific parameters and stable identity for one atmosphere layer."""
struct AtmosphereLayerDefinition{T<:AbstractFloat}
    id::AtmosphereLayerID
    cn2_fraction::T
    wind_speed::T
    wind_direction::T
    altitude::T
end

struct _FixedAtmosphereLayerDefinitions{T<:AbstractFloat} <:
    AbstractVector{AtmosphereLayerDefinition{T}}
    _storage::Tuple{Vararg{AtmosphereLayerDefinition{T}}}
end

Base.size(layers::_FixedAtmosphereLayerDefinitions) =
    (length(getfield(layers, :_storage)),)
Base.axes(layers::_FixedAtmosphereLayerDefinitions) =
    axes(getfield(layers, :_storage))
Base.length(layers::_FixedAtmosphereLayerDefinitions) =
    length(getfield(layers, :_storage))
Base.getindex(layers::_FixedAtmosphereLayerDefinitions, index::Int) =
    getfield(layers, :_storage)[index]
Base.IndexStyle(::Type{<:_FixedAtmosphereLayerDefinitions}) = IndexLinear()
Base.iterate(layers::_FixedAtmosphereLayerDefinitions, state...) =
    iterate(getfield(layers, :_storage), state...)
Base.copy(layers::_FixedAtmosphereLayerDefinitions) =
    collect(getfield(layers, :_storage))

function Base.getproperty(
    layers::_FixedAtmosphereLayerDefinitions,
    name::Symbol,
)
    name === :_storage && return collect(getfield(layers, :_storage))
    return getfield(layers, name)
end

"""Cold declaration of one finite Kolmogorov phase-screen atmosphere."""
struct KolmogorovAtmosphereDefinition{T<:AbstractFloat} <:
    AbstractTimedAtmosphereDefinition
    r0::T
    L0::T
end

function KolmogorovAtmosphereDefinition(;
    r0::Real,
    L0::Real=25.0,
    T::Type{<:AbstractFloat}=Float64,
)
    return KolmogorovAtmosphereDefinition{T}(
        _converted_positive_finite(r0, T, "atmosphere Fried parameter r0"),
        _converted_positive_finite(L0, T, "atmosphere outer scale L0"),
    )
end

"""
Cold declaration of one finite moving-screen multilayer atmosphere.

Layer parameters are copied into fixed-size homogeneous host configuration
storage. Explicit `layer_ids` are required because tuple position and altitude
are not stable stochastic-owner identities.
"""
struct MultiLayerAtmosphereDefinition{T<:AbstractFloat} <:
    AbstractTimedAtmosphereDefinition
    r0::T
    L0::T
    layers::_FixedAtmosphereLayerDefinitions{T}
end

"""Cold declaration of one infinite boundary-injection multilayer atmosphere."""
struct InfiniteMultiLayerAtmosphereDefinition{T<:AbstractFloat} <:
    AbstractTimedAtmosphereDefinition
    r0::T
    L0::T
    layers::_FixedAtmosphereLayerDefinitions{T}
    screen_resolution::Union{Nothing,Int}
    stencil_size::Union{Nothing,Int}
end

function _require_unique_cold_atmosphere_layer_ids(layers)
    @inbounds for index in eachindex(layers)
        id = layers[index].id
        for previous in firstindex(layers):(index - 1)
            layers[previous].id == id && throw(InvalidConfiguration(
                "atmosphere layer identity $(id) is duplicated"))
        end
    end
    return layers
end

function _cold_atmosphere_layer_definitions(
    ::Type{T};
    fractional_cn2,
    wind_speed,
    wind_direction,
    altitude,
    layer_ids,
) where {T<:AbstractFloat}
    isnothing(layer_ids) && throw(InvalidConfiguration(
        "cold multilayer atmosphere declarations require explicit layer_ids"))
    n_layers = length(fractional_cn2)
    n_layers > 0 || throw(InvalidConfiguration(
        "fractional_cn2 cannot be empty"))
    length(wind_speed) == n_layers || throw(InvalidConfiguration(
        "wind_speed length must match fractional_cn2"))
    length(wind_direction) == n_layers || throw(InvalidConfiguration(
        "wind_direction length must match fractional_cn2"))
    length(altitude) == n_layers || throw(InvalidConfiguration(
        "altitude length must match fractional_cn2"))
    length(layer_ids) == n_layers || throw(InvalidConfiguration(
        "layer_ids length must match fractional_cn2"))

    layers = Memory{AtmosphereLayerDefinition{T}}(undef, n_layers)
    @inbounds for index in eachindex(layers)
        layers[index] = AtmosphereLayerDefinition(
            _as_atmosphere_layer_id(layer_ids[index]),
            _converted_nonnegative_finite(fractional_cn2[index], T,
                "atmosphere Cn2 fraction"),
            _converted_nonnegative_finite(wind_speed[index], T,
                "atmosphere wind speed"),
            _converted_finite(wind_direction[index], T,
                "atmosphere wind direction"),
            _converted_nonnegative_finite(altitude[index], T,
                "atmosphere layer altitude"),
        )
    end
    isapprox(sum(layer -> layer.cn2_fraction, layers), one(T);
        atol=T(1e-6), rtol=T(1e-6)) || throw(InvalidConfiguration(
        "fractional_cn2 must sum to 1"))
    _require_unique_cold_atmosphere_layer_ids(layers)
    return _FixedAtmosphereLayerDefinitions{T}(Tuple(layers))
end

function MultiLayerAtmosphereDefinition(;
    r0::Real,
    L0::Real=25.0,
    fractional_cn2,
    wind_speed,
    wind_direction,
    altitude,
    layer_ids,
    T::Type{<:AbstractFloat}=Float64,
)
    layers = _cold_atmosphere_layer_definitions(T;
        fractional_cn2,
        wind_speed,
        wind_direction,
        altitude,
        layer_ids,
    )
    return MultiLayerAtmosphereDefinition{T}(
        _converted_positive_finite(r0, T, "atmosphere Fried parameter r0"),
        _converted_positive_finite(L0, T, "atmosphere outer scale L0"),
        layers,
    )
end

@inline function _optional_positive_resolution(
    ::Nothing,
    ::AbstractString,
)
    return nothing
end

function _optional_positive_resolution(::Bool, label::AbstractString)
    throw(InvalidConfiguration(
        "$label must be an integer count, not Bool"))
end

function _optional_positive_resolution(value::Integer, label::AbstractString)
    value > 0 || throw(InvalidConfiguration("$label must be positive"))
    value <= typemax(Int) || throw(InvalidConfiguration(
        "$label exceeds the supported Int range"))
    return Int(value)
end

function _optional_positive_resolution(value, label::AbstractString)
    throw(InvalidConfiguration(
        "$label must be an integer count or nothing; got $(typeof(value))"))
end

function InfiniteMultiLayerAtmosphereDefinition(;
    r0::Real,
    L0::Real=25.0,
    fractional_cn2,
    wind_speed,
    wind_direction,
    altitude,
    layer_ids,
    screen_resolution=nothing,
    stencil_size=nothing,
    T::Type{<:AbstractFloat}=Float64,
)
    layers = _cold_atmosphere_layer_definitions(T;
        fractional_cn2,
        wind_speed,
        wind_direction,
        altitude,
        layer_ids,
    )
    return InfiniteMultiLayerAtmosphereDefinition{T}(
        _converted_positive_finite(r0, T, "atmosphere Fried parameter r0"),
        _converted_positive_finite(L0, T, "atmosphere outer scale L0"),
        layers,
        _optional_positive_resolution(screen_resolution,
            "infinite atmosphere screen_resolution"),
        _optional_positive_resolution(stencil_size,
            "infinite atmosphere stencil_size"),
    )
end

@inline _definition_cn2_fractions(definition) =
    map(layer -> layer.cn2_fraction, definition.layers)

@inline _definition_wind_speeds(definition) =
    map(layer -> layer.wind_speed, definition.layers)

@inline _definition_wind_directions(definition) =
    map(layer -> layer.wind_direction, definition.layers)

@inline _definition_altitudes(definition) =
    map(layer -> layer.altitude, definition.layers)

@inline _definition_layer_ids(definition) =
    map(layer -> layer.id, definition.layers)

function _prepare_timed_atmosphere(
    definition::KolmogorovAtmosphereDefinition{T},
    telescope::Telescope,
    selector::AbstractArrayBackend,
) where {T}
    return KolmogorovAtmosphere(telescope;
        r0=definition.r0,
        L0=definition.L0,
        T,
        backend=selector,
    )
end

function _prepare_timed_atmosphere(
    definition::MultiLayerAtmosphereDefinition{T},
    telescope::Telescope,
    selector::AbstractArrayBackend,
) where {T}
    return MultiLayerAtmosphere(telescope;
        r0=definition.r0,
        L0=definition.L0,
        fractional_cn2=_definition_cn2_fractions(definition),
        wind_speed=_definition_wind_speeds(definition),
        wind_direction=_definition_wind_directions(definition),
        altitude=_definition_altitudes(definition),
        layer_ids=_definition_layer_ids(definition),
        T,
        backend=selector,
    )
end

function _prepare_timed_atmosphere(
    definition::InfiniteMultiLayerAtmosphereDefinition{T},
    telescope::Telescope,
    selector::AbstractArrayBackend,
) where {T}
    screen_resolution = something(definition.screen_resolution,
        default_infinite_screen_resolution(telescope.params.resolution))
    stencil_size = something(definition.stencil_size,
        default_infinite_stencil_size(telescope.params.resolution))
    return InfiniteMultiLayerAtmosphere(telescope;
        r0=definition.r0,
        L0=definition.L0,
        fractional_cn2=_definition_cn2_fractions(definition),
        wind_speed=_definition_wind_speeds(definition),
        wind_direction=_definition_wind_directions(definition),
        altitude=_definition_altitudes(definition),
        layer_ids=_definition_layer_ids(definition),
        screen_resolution,
        stencil_size,
        T,
        backend=selector,
    )
end

function _prepare_timed_atmosphere(
    definition::AbstractTimedAtmosphereDefinition,
    ::Telescope,
    ::AbstractArrayBackend,
)
    throw(InvalidConfiguration(
        "timed atmosphere definition $(typeof(definition)) does not " *
        "implement target-local preparation"))
end

function _require_atmosphere_target_array(array, target, label)
    actual = compute_device(array)
    actual == target || _throw_compute_device_error(
        :prepare_timed_atmosphere,
        :wrong_device,
        target,
        "prepared $label occupies $(actual)",
    )
    return array
end

function _require_kolmogorov_target(
    atmosphere::KolmogorovAtmosphere,
    target::AbstractComputeDevice,
)
    state = atmosphere.state
    _require_atmosphere_target_array(state.opd, target, "atmosphere OPD")
    _require_atmosphere_target_array(state.psd, target, "atmosphere PSD")
    _require_atmosphere_target_array(state.spectrum, target,
        "atmosphere spectrum")
    _require_atmosphere_target_array(state.noise_re, target,
        "atmosphere real noise workspace")
    _require_atmosphere_target_array(state.noise_im, target,
        "atmosphere imaginary noise workspace")
    _require_atmosphere_target_array(state.freqs, target,
        "atmosphere frequency workspace")
    return atmosphere
end

function _require_telescope_target(telescope::Telescope, target, label)
    _require_atmosphere_target_array(pupil_mask(telescope), target,
        "$label pupil mask")
    _require_atmosphere_target_array(pupil_reflectivity(telescope), target,
        "$label pupil reflectivity")
    return telescope
end

function validate_timed_atmosphere_target(
    atmosphere::KolmogorovAtmosphere,
    target::AbstractComputeDevice,
)
    return _require_kolmogorov_target(atmosphere, target)
end

function validate_timed_atmosphere_target(
    atmosphere::MultiLayerAtmosphere,
    target::AbstractComputeDevice,
)
    for layer in atmosphere.layers
        _require_kolmogorov_target(layer.generator, target)
        _require_telescope_target(layer.generator_telescope, target,
            "moving-layer generator telescope")
    end
    return atmosphere
end

function _require_infinite_boundary_target(model, target, label)
    stencil = model.stencil
    operator = model.operator
    for (name, array) in (
        (:stencil_coords, stencil.stencil_coords),
        (:boundary_coords, stencil.boundary_coords),
        (:stencil_positions, stencil.stencil_positions),
        (:boundary_positions, stencil.boundary_positions),
        (:predictor, operator.predictor),
        (:residual_factor, operator.residual_factor),
        (:cov_zz, operator.cov_zz),
        (:cov_xx, operator.cov_xx),
        (:cov_xz, operator.cov_xz),
        (:cov_zx, operator.cov_zx),
        (:residual_covariance, operator.residual_covariance),
        (:singular_values, operator.singular_values),
    )
        _require_atmosphere_target_array(array, target, "$label $name")
    end
    return model
end

function _require_infinite_screen_target(screen, target)
    state = screen.state
    for (name, array) in (
        (:screen, state.screen),
        (:screen_scratch, state.screen_scratch),
        (:extract_buffer, state.extract_buffer),
        (:stencil_buffer, state.stencil_buffer),
        (:boundary_buffer, state.boundary_buffer),
        (:noise_buffer, state.noise_buffer),
        (:column_positive_coords, state.column_positive_coords),
        (:column_negative_coords, state.column_negative_coords),
        (:row_positive_coords, state.row_positive_coords),
        (:row_negative_coords, state.row_negative_coords),
    )
        _require_atmosphere_target_array(array, target,
            "infinite atmosphere $name")
    end
    _require_infinite_boundary_target(state.column_positive, target,
        "positive-column boundary")
    _require_infinite_boundary_target(state.column_negative, target,
        "negative-column boundary")
    _require_infinite_boundary_target(state.row_positive, target,
        "positive-row boundary")
    _require_infinite_boundary_target(state.row_negative, target,
        "negative-row boundary")
    _require_kolmogorov_target(screen.generator, target)
    _require_telescope_target(screen.generator_telescope, target,
        "infinite-screen generator telescope")
    return screen
end

function validate_timed_atmosphere_target(
    atmosphere::InfiniteMultiLayerAtmosphere,
    target::AbstractComputeDevice,
)
    for layer in atmosphere.layers
        _require_infinite_screen_target(layer.screen, target)
    end
    return atmosphere
end

"""Qualified fail-closed exact-target validation seam for timed atmospheres."""
function validate_timed_atmosphere_target(
    prepared::AbstractTimedAtmosphere,
    ::AbstractComputeDevice,
)
    throw(InvalidConfiguration(
        "prepared timed atmosphere $(typeof(prepared)) does not implement " *
        "validate_timed_atmosphere_target"))
end

function validate_timed_atmosphere_target(prepared, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "timed atmosphere preparation must return AbstractTimedAtmosphere; " *
        "got $(typeof(prepared))"))
end

"""
    prepare_timed_atmosphere(definition, telescope, target)

Prepare one independent numerical timed atmosphere on an exact compute device.
The caller's previous accelerator context is restored on success and failure.
Host configuration and deliberate RNG staging buffers are not accelerator
data-plane arrays. This is a cold function-barrier boundary: backend FFT-plan
construction may be inference-opaque, while the returned runtime value is
concrete and subsequent numerical execution specializes on that value's type.
"""
function prepare_timed_atmosphere(
    definition::D,
    telescope::Telescope,
    target::AbstractComputeDevice,
) where {D<:AbstractTimedAtmosphereDefinition}
    ismutabletype(D) && throw(InvalidConfiguration(
        "timed atmosphere declarations must be immutable"))
    telescope_target = compute_device(pupil_mask(telescope))
    telescope_target == target || _throw_compute_device_error(
        :prepare_timed_atmosphere,
        :wrong_device,
        target,
        "prepared telescope occupies $(telescope_target)",
    )
    selector = compute_device_backend(target)
    return _with_compute_device(target) do
        prepared = _prepare_timed_atmosphere(
            definition, telescope, selector)
        validate_timed_atmosphere_target(prepared, target)
        return prepared
    end
end
