"""Internal marker for atmosphere layers rendered by prepared paths."""
abstract type AbstractAtmosphereLayer end

@inline layer_altitude(layer::AbstractAtmosphereLayer) = layer.params.altitude

# Static or extension-defined atmospheres without directional geometry retain
# their ordinary propagation semantics when a runtime supplies a source.
# Timed atmosphere models override this with the prepared-renderer path below.
@inline propagate!(atm::AbstractAtmosphere, pupil::PupilFunction,
    ::AbstractSource) = propagate!(atm, pupil)

"""Per-atmosphere identity shared by its epoch values and prepared renderers."""
mutable struct AtmosphereIdentity end

"""
    AtmosphereEpoch

Immutable identity token for one published atmosphere state. `model_time` is
expressed in seconds in the atmosphere's numeric type; `sequence` advances
once for every newly published state. The token does not retain layer storage:
the current renderer contract accepts it only while that state remains the
atmosphere's current publication.
"""
struct AtmosphereEpoch{T<:AbstractFloat,I<:AtmosphereIdentity}
    identity::I
    model_time::T
    sequence::UInt64
end

"""Return an atmosphere epoch's explicit model time in seconds."""
@inline epoch_time(epoch::AtmosphereEpoch) = epoch.model_time

"""Return an atmosphere epoch's monotonic publication sequence."""
@inline epoch_sequence(epoch::AtmosphereEpoch) = epoch.sequence

Base.:(==)(a::AtmosphereEpoch, b::AtmosphereEpoch) =
    a.identity === b.identity &&
    isequal(a.model_time, b.model_time) &&
    a.sequence == b.sequence

mutable struct AtmosphereTimelineState{T<:AbstractFloat}
    model_time::T
    sequence::UInt64
    initialized::Bool
end

@inline new_atmosphere_timeline(::Type{T}) where {T<:AbstractFloat} =
    AtmosphereTimelineState{T}(zero(T), UInt64(0), false)

@inline atmosphere_timeline(atm::AbstractTimedAtmosphere) = atm.state.timeline
@inline atmosphere_identity(atm::AbstractTimedAtmosphere) = atm.identity

"""Return the currently published epoch, or throw if none has been published."""
function current_epoch(atm::AbstractTimedAtmosphere)
    timeline = atmosphere_timeline(atm)
    timeline.initialized || throw(AtmosphereEpochError(
        "atmosphere has no published epoch; call advance_by! or advance_to! first"))
    return AtmosphereEpoch(atmosphere_identity(atm), timeline.model_time,
        timeline.sequence)
end

function _explicit_atmosphere_duration(duration::Real, ::Type{T}) where {T<:AbstractFloat}
    value = T(duration)
    isfinite(value) || throw(AtmosphereTimeError(
        "atmosphere elapsed duration must be finite"))
    value >= zero(T) || throw(AtmosphereTimeError(
        "atmosphere elapsed duration must be non-negative"))
    return value
end

function _explicit_atmosphere_target(target_time::Real, ::Type{T}) where {T<:AbstractFloat}
    value = T(target_time)
    isfinite(value) || throw(AtmosphereTimeError(
        "atmosphere target model time must be finite"))
    value >= zero(T) || throw(AtmosphereTimeError(
        "atmosphere target model time must be non-negative"))
    return value
end

# Concrete models own initialization, wind/refresh evolution, boundary
# extension, and all RNG consumption through these dispatch points.
function initialize_atmosphere!(atm::AbstractTimedAtmosphere, ::AbstractRNG)
    throw(UnsupportedAlgorithm(
        "$(typeof(atm)) does not implement explicit atmosphere initialization"))
end

function evolve_atmosphere!(atm::AbstractTimedAtmosphere, ::Real, ::AbstractRNG)
    throw(UnsupportedAlgorithm(
        "$(typeof(atm)) does not implement explicit atmosphere evolution"))
end


@inline evolve_initial_atmosphere!(atm::AbstractTimedAtmosphere, duration::Real,
    rng::AbstractRNG) = evolve_atmosphere!(atm, duration, rng)

"""
    advance_by!(atmosphere, elapsed_seconds; rng)

Advance one shared atmosphere writer by an explicit elapsed duration and return
the newly published epoch. The first call initializes the stochastic state;
thereafter a zero duration returns the existing epoch without consuming RNG or
changing its sequence.
"""
function advance_by!(atm::AbstractTimedAtmosphere, duration::Real,
    rng::AbstractRNG)
    timeline = atmosphere_timeline(atm)
    T = typeof(timeline.model_time)
    dt = _explicit_atmosphere_duration(duration, T)

    if timeline.initialized && iszero(dt)
        return current_epoch(atm)
    end

    next_time = timeline.model_time + dt
    isfinite(next_time) || throw(AtmosphereTimeError(
        "atmosphere model time overflowed its numeric representation"))
    timeline.sequence == typemax(UInt64) && throw(AtmosphereTimeError(
        "atmosphere epoch sequence is exhausted"))

    if timeline.initialized
        evolve_atmosphere!(atm, dt, rng)
    else
        initialize_atmosphere!(atm, rng)
        evolve_initial_atmosphere!(atm, dt, rng)
    end

    timeline.model_time = next_time
    timeline.sequence += UInt64(1)
    timeline.initialized = true
    return current_epoch(atm)
end

function advance_by!(atm::AbstractTimedAtmosphere, duration::Real;
    rng::AbstractRNG=Random.default_rng())
    return advance_by!(atm, duration, rng)
end

"""
    advance_to!(atmosphere, target_model_time; rng)

Advance to an explicit absolute model time. Equal-time calls are no-ops after
initialization. Backward targets fail before atmosphere state is mutated.
"""
function advance_to!(atm::AbstractTimedAtmosphere, target_time::Real,
    rng::AbstractRNG)
    timeline = atmosphere_timeline(atm)
    T = typeof(timeline.model_time)
    target = _explicit_atmosphere_target(target_time, T)
    target < timeline.model_time && throw(AtmosphereTimeError(
        "atmosphere target model time $target precedes current time $(timeline.model_time)"))
    return advance_by!(atm, target - timeline.model_time, rng)
end

function advance_to!(atm::AbstractTimedAtmosphere, target_time::Real;
    rng::AbstractRNG=Random.default_rng())
    return advance_to!(atm, target_time, rng)
end

@inline function source_geometry_signature(::Nothing,
    ::Type{T}) where {T<:AbstractFloat}
    return zero(T), zero(T), T(Inf)
end

function source_geometry_signature(src::AbstractSource,
    ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec = coordinates_xy_arcsec(src)
    height_m = source_height_m(src)
    return T(x_arcsec), T(y_arcsec),
        isfinite(height_m) ? T(height_m) : T(Inf)
end

@inline function is_onaxis_infinite_source(src::Union{AbstractSource,Nothing},
    ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec, height_m = source_geometry_signature(src, T)
    return iszero(x_arcsec) && iszero(y_arcsec) && !isfinite(height_m)
end

@inline function layer_source_geometry(::Nothing, layer_altitude::Real,
    sampling_m::NTuple{2,<:Real}, ::Type{T}) where {T<:AbstractFloat}
    return zero(T), zero(T), one(T)
end

function layer_source_geometry(src::AbstractSource, layer_altitude::Real,
    sampling_m::NTuple{2,<:Real}, ::Type{T}) where {T<:AbstractFloat}
    x_arcsec, y_arcsec, height_m = source_geometry_signature(src, T)
    altitude_t = T(layer_altitude)
    arcsec_to_rad = T(pi / (180 * 3600))
    shift_x = x_arcsec * arcsec_to_rad * altitude_t / T(sampling_m[2])
    shift_y = y_arcsec * arcsec_to_rad * altitude_t / T(sampling_m[1])

    footprint_scale = one(T)
    if isfinite(height_m)
        height_m > altitude_t || throw(InvalidConfiguration(
            "source height must exceed layer altitude for finite-height propagation"))
        footprint_scale = (height_m - altitude_t) / height_m
    end
    return shift_x, shift_y, footprint_scale
end

@inline layer_source_geometry(src::Union{AbstractSource,Nothing},
    layer_altitude::Real, tel::Telescope, ::Type{T}) where {T<:AbstractFloat} =
    layer_source_geometry(src, layer_altitude, tel.aperture.sampling_m, T)

struct LayerRenderContext{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    footprint_scale::T
end

@inline function layer_render_context(src::Union{AbstractSource,Nothing},
    layer_altitude::Real, sampling_m::NTuple{2,<:Real},
    ::Type{T}) where {T<:AbstractFloat}
    shift_x, shift_y, footprint_scale =
        layer_source_geometry(src, layer_altitude, sampling_m, T)
    return LayerRenderContext{T}(shift_x, shift_y, footprint_scale)
end

@inline layer_render_context(src::Union{AbstractSource,Nothing},
    layer_altitude::Real, tel::Telescope, ::Type{T}) where {T<:AbstractFloat} =
    layer_render_context(src, layer_altitude, tel.aperture.sampling_m, T)

@inline layer_render_context(src::Union{AbstractSource,Nothing},
    layer::AbstractAtmosphereLayer, tel::Telescope,
    ::Type{T}) where {T<:AbstractFloat} =
    layer_render_context(src, layer_altitude(layer), tel, T)

@inline render_layer!(out::AbstractMatrix{T}, layer::AbstractAtmosphereLayer,
    ctx::LayerRenderContext{T}) where {T<:AbstractFloat} =
    render_layer!(out, layer, ctx.shift_x, ctx.shift_y,
        ctx.footprint_scale)

@inline render_layer_accumulate!(out::AbstractMatrix{T},
    layer::AbstractAtmosphereLayer,
    ctx::LayerRenderContext{T}) where {T<:AbstractFloat} =
    render_layer_accumulate!(out, layer, ctx.shift_x, ctx.shift_y,
        ctx.footprint_scale)

@inline atmosphere_layers(atm::AbstractTimedAtmosphere) = atm.layers

function atmosphere_numeric_type(atm::AbstractTimedAtmosphere)
    throw(UnsupportedAlgorithm(
        "$(typeof(atm)) does not declare its atmosphere numeric type"))
end

"""
Prepared, path-local geometry for one source direction. Its source description,
geometry arrays, pupil support, output metadata, and atmosphere identity are
fixed for the lifetime of the prepared run.
"""
struct AtmosphereDirectionRenderer{I,S,V,M,A}
    identity::I
    source::S
    shift_x::V
    shift_y::V
    footprint_scale::V
    output_metadata::M
    pupil::A
end

struct PreparedAtmosphereDirections{S,R}
    source::S
    renderers::R
end

"""Return the frozen per-direction renderers in a plural preparation."""
@inline direction_renderers(prepared::PreparedAtmosphereDirections) =
    prepared.renderers

@inline freeze_source(::Nothing) = nothing

"""
    prepare_atmosphere_renderer(atmosphere, telescope, source=nothing)

Prepare immutable geometry for one frozen source direction. Use
`prepare_atmosphere_renderers` for an asterism or extended source.
"""
function prepare_atmosphere_renderer(atm::AbstractTimedAtmosphere,
    tel::Telescope, src::Union{AbstractSource,Nothing}=nothing;
    T::Type{<:AbstractFloat}=eltype(pupil_reflectivity(tel)))
    require_same_backend(atm, tel)
    atmosphere_numeric_type(atm) === T || throw(InvalidConfiguration(
        "atmosphere renderer numeric type must match atmosphere layer storage"))
    frozen_source = freeze_source(src)
    layers = atmosphere_layers(atm)
    shift_x = Vector{T}(undef, length(layers))
    shift_y = Vector{T}(undef, length(layers))
    footprint_scale = Vector{T}(undef, length(layers))
    sampling_m = (T(tel.aperture.sampling_m[1]),
        T(tel.aperture.sampling_m[2]))
    @inbounds for i in eachindex(layers)
        ctx = layer_render_context(frozen_source, layer_altitude(layers[i]),
            sampling_m, T)
        shift_x[i] = ctx.shift_x
        shift_y[i] = ctx.shift_y
        footprint_scale[i] = ctx.footprint_scale
    end

    template = PupilFunction(tel; T=T)
    frozen_pupil = copy(pupil_mask(tel))
    return AtmosphereDirectionRenderer(
        atmosphere_identity(atm),
        frozen_source,
        shift_x,
        shift_y,
        footprint_scale,
        template.metadata,
        frozen_pupil,
    )
end

function prepare_atmosphere_renderer(::AbstractTimedAtmosphere,
    ::Telescope, ::Asterism; kwargs...)
    throw(InvalidConfiguration(
        "an asterism contains multiple directions; use prepare_atmosphere_renderers"))
end

function prepare_atmosphere_renderer(::AbstractTimedAtmosphere,
    ::Telescope, ::ExtendedSource; kwargs...)
    throw(InvalidConfiguration(
        "an extended source contains multiple directions; use prepare_atmosphere_renderers"))
end

function _prepare_renderer_vector(atm::AbstractTimedAtmosphere,
    tel::Telescope, sources::AbstractVector{S};
    T::Type{<:AbstractFloat}) where {S<:AbstractSource}
    isempty(sources) && throw(InvalidConfiguration(
        "prepared atmosphere directions require at least one source"))
    first_renderer = prepare_atmosphere_renderer(atm, tel, first(sources);
        T=T)
    renderers = Vector{typeof(first_renderer)}(undef, length(sources))
    renderers[1] = first_renderer
    @inbounds for i in 2:length(sources)
        renderers[i] = prepare_atmosphere_renderer(atm, tel, sources[i]; T=T)
    end
    return renderers
end

"""
    prepare_atmosphere_renderers(atmosphere, telescope, source)

Freeze and prepare all directions represented by a point, asterism, or extended
source. Access the result with `direction_renderers`.
"""
function prepare_atmosphere_renderers(atm::AbstractTimedAtmosphere,
    tel::Telescope, src::AbstractSource;
    T::Type{<:AbstractFloat}=eltype(pupil_reflectivity(tel)))
    frozen = freeze_source(src)
    renderer = prepare_atmosphere_renderer(atm, tel, frozen; T=T)
    return PreparedAtmosphereDirections(frozen, (renderer,))
end

function prepare_atmosphere_renderers(atm::AbstractTimedAtmosphere,
    tel::Telescope, ast::Asterism;
    T::Type{<:AbstractFloat}=eltype(pupil_reflectivity(tel)))
    frozen = freeze_source(ast)
    renderers = map(src -> prepare_atmosphere_renderer(atm, tel, src; T=T),
        frozen.sources)
    return PreparedAtmosphereDirections(frozen, renderers)
end

function prepare_atmosphere_renderers(atm::AbstractTimedAtmosphere,
    tel::Telescope, src::ExtendedSource;
    T::Type{<:AbstractFloat}=eltype(pupil_reflectivity(tel)))
    frozen = freeze_source(src)
    asterism = extended_source_asterism(frozen)
    renderers = _prepare_renderer_vector(atm, tel, asterism.sources; T=T)
    return PreparedAtmosphereDirections(asterism, renderers)
end

function _validate_epoch_identity(identity::AtmosphereIdentity,
    atm::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    identity === atmosphere_identity(atm) || throw(AtmosphereEpochError(
        "prepared renderer belongs to a different atmosphere"))
    epoch.identity === identity || throw(AtmosphereEpochError(
        "epoch belongs to a different atmosphere"))
    published = current_epoch(atm)
    epoch == published || throw(AtmosphereEpochError(
        "epoch $(epoch.sequence) at model time $(epoch.model_time) is stale or unpublished"))
    return published
end

function _validate_atmosphere_destination(dest::AbstractMatrix,
    renderer::AtmosphereDirectionRenderer)
    metadata = renderer.output_metadata
    size(dest) == metadata.dimensions || throw(DimensionMismatchError(
        "atmosphere destination dimensions do not match prepared output geometry"))
    eltype(dest) === metadata.numeric_type || throw(InvalidConfiguration(
        "atmosphere destination numeric type does not match prepared output"))
    typeof(backend(dest)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "atmosphere destination backend does not match prepared output"))
    plane_device(dest) == metadata.device || throw(InvalidConfiguration(
        "atmosphere destination is on a different physical device"))
    return dest
end

function _validate_atmosphere_destination(dest::PupilFunction,
    renderer::AtmosphereDirectionRenderer)
    require_same_plane_grid(dest.metadata, renderer.output_metadata;
        label="atmosphere renderer and destination",
        require_spectral=false)
    validate_plane_storage(dest.metadata, dest.opd;
        label="atmosphere destination OPD")
    return dest.opd
end

function _validate_atmosphere_renderer_binding(
    renderer::AtmosphereDirectionRenderer,
    atm::AbstractTimedAtmosphere,
)
    renderer.identity === atmosphere_identity(atm) || throw(
        AtmosphereEpochError(
            "prepared renderer belongs to a different atmosphere"))
    n_layers = length(atmosphere_layers(atm))
    (length(renderer.shift_x) == n_layers &&
        length(renderer.shift_y) == n_layers &&
        length(renderer.footprint_scale) == n_layers) || throw(
        AtmosphereEpochError(
            "atmosphere layer shape changed after renderer preparation"))
    metadata = renderer.output_metadata
    size(renderer.pupil) == metadata.dimensions || throw(
        AtmosphereEpochError(
            "prepared renderer pupil shape changed after preparation"))
    typeof(backend(renderer.pupil)) === typeof(metadata.backend) || throw(
        AtmosphereEpochError(
            "prepared renderer pupil backend changed after preparation"))
    plane_device(renderer.pupil) == metadata.device || throw(
        AtmosphereEpochError(
            "prepared renderer pupil device changed after preparation"))
    return renderer
end

"""
    validate_atmosphere_rendering(destination, renderer, atmosphere, epoch)

Validate the complete current-epoch rendering binding without mutating the
destination. Serial plant execution uses this preflight operation for every
selected path before materializing any path input.
"""
function validate_atmosphere_rendering(dest,
    renderer::AtmosphereDirectionRenderer, atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch)
    _validate_epoch_identity(renderer.identity, atm, epoch)
    _validate_atmosphere_renderer_binding(renderer, atm)
    _validate_atmosphere_destination(dest, renderer)
    return dest
end

function accumulate_rendered_layers!(opd::AbstractMatrix,
    layers, shift_x::AbstractVector, shift_y::AbstractVector,
    footprint_scale::AbstractVector)
    isempty(layers) && return fill!(opd, zero(eltype(opd)))
    render_layer!(opd, layers[1], shift_x[1], shift_y[1],
        footprint_scale[1])
    @inbounds for i in 2:length(layers)
        render_layer_accumulate!(opd, layers[i], shift_x[i], shift_y[i],
            footprint_scale[i])
    end
    return opd
end

function render_atmosphere_opd_impl!(dest::AbstractMatrix,
    renderer::AtmosphereDirectionRenderer, atm::AbstractTimedAtmosphere)
    layers = atmosphere_layers(atm)
    accumulate_rendered_layers!(dest, layers, renderer.shift_x,
        renderer.shift_y, renderer.footprint_scale)
    dest .*= renderer.pupil
    return dest
end

"""
    render_atmosphere_opd!(destination, renderer, atmosphere, epoch)

Render one frozen epoch into a caller-owned OPD matrix. Epoch, shape, backend,
and physical-device compatibility are checked before output mutation.
"""
function render_atmosphere_opd!(dest::AbstractMatrix,
    renderer::AtmosphereDirectionRenderer, atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch)
    validate_atmosphere_rendering(dest, renderer, atm, epoch)
    return render_atmosphere_opd_impl!(dest, renderer, atm)
end

"""
    render_atmosphere!(destination, renderer, atmosphere, epoch)

Render the current compatible epoch through one prepared direction into a
caller-owned `PupilFunction`. Validation completes before output mutation.
"""
function render_atmosphere!(dest::PupilFunction,
    renderer::AtmosphereDirectionRenderer, atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch)
    validate_atmosphere_rendering(dest, renderer, atm, epoch)
    opd = dest.opd
    render_atmosphere_opd_impl!(opd, renderer, atm)
    return dest
end
