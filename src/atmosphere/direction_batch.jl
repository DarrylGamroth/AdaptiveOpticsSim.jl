"""
Capability trait for preparing a homogeneous atmosphere-direction batch.

Concrete atmosphere models opt in only when every direction can be rendered
from the same ordered set of sampled layer screens.
"""
abstract type AbstractAtmosphereDirectionBatchCapability end

"""The atmosphere exposes ordered sampled screens suitable for batching."""
struct ExtractedScreenDirectionBatchCapability <:
    AbstractAtmosphereDirectionBatchCapability end

"""The atmosphere retains its existing single-direction rendering path."""
struct UnsupportedAtmosphereDirectionBatchCapability <:
    AbstractAtmosphereDirectionBatchCapability end

struct AtmosphereSourceGeometryBinding{
    S,
    T<:AbstractFloat,
}
    source::S
    signature::NTuple{3,T}
end

Base.@noinline function _source_geometry_binding_matches(
    binding::AtmosphereSourceGeometryBinding{S,T},
) where {S,T<:AbstractFloat}
    return source_geometry_signature(binding.source, T) == binding.signature
end

@inline atmosphere_direction_batch_capability(
    ::Type{<:AbstractTimedAtmosphere},
) = UnsupportedAtmosphereDirectionBatchCapability()

@inline atmosphere_direction_batch_capability(
    ::Type{<:MultiLayerAtmosphere},
) = ExtractedScreenDirectionBatchCapability()

@inline atmosphere_direction_batch_capability(
    ::Type{<:InfiniteMultiLayerAtmosphere},
) = ExtractedScreenDirectionBatchCapability()

@inline _direction_batch_layer_screen(layer::MovingAtmosphereLayer) =
    layer.generator.state.opd
@inline _direction_batch_layer_screen(layer::InfiniteAtmosphereLayer) =
    layer.screen.state.screen

@inline _direction_batch_layer_sampling(layer::MovingAtmosphereLayer) =
    layer.params.sampling_m
@inline _direction_batch_layer_sampling(layer::InfiniteAtmosphereLayer) =
    layer.screen.params.pixel_scale

@inline _direction_batch_layer_offset_x(layer::AbstractAtmosphereLayer) =
    layer.state.offset_x
@inline _direction_batch_layer_offset_y(layer::AbstractAtmosphereLayer) =
    layer.state.offset_y
@inline _direction_batch_layer_amplitude(layer::AbstractAtmosphereLayer) =
    layer.params.amplitude_scale

"""
Immutable scientific and storage contract for a prepared atmosphere-direction
batch. Direction cardinality remains a runtime value rather than a type
parameter so large asterisms do not produce cardinality-specific method
instances.
"""
struct AtmosphereDirectionBatchParams{
    T<:AbstractFloat,
    I<:AtmosphereIdentity,
    C<:AbstractAtmosphereDirectionBatchCapability,
    S<:FixedSizeVector,
    G<:FixedSizeVector,
    R<:FixedSizeVector,
    L<:FixedSizeVector,
    M<:OpticalPlaneMetadata,
}
    identity::I
    capability::C
    sources::S
    source_geometry_bindings::G
    source_bindings::R
    layer_bindings::L
    output_metadata::M
    direction_count::Int
    direction_capacity::Int
    layer_count::Int
end

function AtmosphereDirectionBatchParams{T}(
    identity::I,
    capability::C,
    sources::S,
    source_geometry_bindings::G,
    source_bindings::R,
    layer_bindings::L,
    output_metadata::M,
    direction_count::Int,
    direction_capacity::Int,
    layer_count::Int,
) where {
    T<:AbstractFloat,
    I<:AtmosphereIdentity,
    C<:AbstractAtmosphereDirectionBatchCapability,
    S<:FixedSizeVector,
    G<:FixedSizeVector,
    R<:FixedSizeVector,
    L<:FixedSizeVector,
    M<:OpticalPlaneMetadata,
}
    return AtmosphereDirectionBatchParams{
        T,I,C,S,G,R,L,M,
    }(
        identity,
        capability,
        sources,
        source_geometry_bindings,
        source_bindings,
        layer_bindings,
        output_metadata,
        direction_count,
        direction_capacity,
        layer_count,
    )
end

"""Caller-owned output and prepared device-resident batch geometry."""
mutable struct AtmosphereDirectionBatchWorkspace{
    X<:AbstractMatrix,
    Y<:AbstractMatrix,
    F<:AbstractMatrix,
    P<:AbstractMatrix{Bool},
    O<:AbstractArray{<:AbstractFloat,3},
}
    shift_x::X
    shift_y::Y
    footprint_scale::F
    pupil::P
    output::O
end

"""Prepared immutable contract plus its path-local mutable workspace."""
struct PreparedAtmosphereDirectionBatch{
    P<:AtmosphereDirectionBatchParams,
    W<:AtmosphereDirectionBatchWorkspace,
}
    params::P
    workspace::W
end

"""Return the caller-owned atmospheric OPD stack."""
@inline atmosphere_direction_output(
    prepared::PreparedAtmosphereDirectionBatch,
) = prepared.workspace.output

"""Return the number of ordered source directions rendered by the batch."""
@inline atmosphere_direction_count(
    prepared::PreparedAtmosphereDirectionBatch,
) = prepared.params.direction_count

"""Return the fixed direction capacity of the caller-owned output."""
@inline atmosphere_direction_capacity(
    prepared::PreparedAtmosphereDirectionBatch,
) = prepared.params.direction_capacity

"""Return the pupil-plane metadata shared by every OPD-stack slice."""
@inline atmosphere_direction_metadata(
    prepared::PreparedAtmosphereDirectionBatch,
) = prepared.params.output_metadata

@inline plane_metadata(prepared::PreparedAtmosphereDirectionBatch) =
    atmosphere_direction_metadata(prepared)

@inline function _freeze_batch_sources(::Nothing)
    return _fixed_size_union_vector((nothing,))
end

function _freeze_batch_sources(src::AbstractSource)
    require_leaf_source(src, "atmosphere direction batch source")
    return _fixed_size_union_vector((freeze_source(src),))
end

function _freeze_batch_sources(
    sources::AbstractVector{<:AbstractSource},
)
    isempty(sources) && throw(InvalidConfiguration(
        "an atmosphere direction batch requires at least one source"))
    require_leaf_source(first(sources),
        "atmosphere direction batch source")
    first_source = freeze_source(first(sources))
    S = typeof(first_source)
    frozen = Vector{S}(undef, length(sources))
    frozen[1] = first_source
    @inbounds for direction in 2:length(sources)
        require_leaf_source(sources[direction],
            "atmosphere direction batch source")
        source = freeze_source(sources[direction])
        typeof(source) === S || throw(InvalidConfiguration(
            "atmosphere direction batch sources must have one concrete " *
            "numeric and radiometric storage contract"))
        frozen[direction] = source
    end
    return FixedSizeVectorDefault{S}(frozen)
end

function _freeze_batch_sources(ast::Asterism)
    frozen = freeze_source(ast)
    isempty(frozen.sources) && throw(InvalidConfiguration(
        "an atmosphere direction batch requires at least one source"))
    return frozen.sources
end

function _freeze_batch_sources(src::ExtendedSource)
    frozen = freeze_source(src)
    sources = extended_source_asterism(frozen).sources
    isempty(sources) && throw(InvalidConfiguration(
        "an atmosphere direction batch requires at least one source"))
    return sources
end

function _copy_source_geometry_bindings(
    sources,
    ::Type{T},
) where {T<:AbstractFloat}
    bindings = (
        AtmosphereSourceGeometryBinding(
            source,
            source_geometry_signature(source, T),
        )
        for source in sources
    )
    return _fixed_size_union_vector(bindings)
end

function _copy_source_bindings(sources)
    return _fixed_size_union_vector(sources)
end

function _copy_layer_bindings(layers::AbstractVector)
    return _fixed_size_union_vector(layers)
end

function _prepare_batch_geometry(
    layers,
    sources,
    tel::Telescope,
    ::Type{T},
) where {T<:AbstractFloat}
    n_layers = length(layers)
    n_directions = length(sources)
    shift_x = Matrix{T}(undef, n_layers, n_directions)
    shift_y = Matrix{T}(undef, n_layers, n_directions)
    footprint_scale = Matrix{T}(undef, n_layers, n_directions)
    sampling_m = (
        T(tel.aperture.sampling_m[1]),
        T(tel.aperture.sampling_m[2]),
    )
    @inbounds for direction in eachindex(sources)
        source = sources[direction]
        for layer_index in eachindex(layers)
            context = layer_render_context(
                source,
                layer_altitude(layers[layer_index]),
                sampling_m,
                T,
            )
            shift_x[layer_index, direction] = context.shift_x
            shift_y[layer_index, direction] = context.shift_y
            footprint_scale[layer_index, direction] =
                context.footprint_scale
        end
    end
    return shift_x, shift_y, footprint_scale
end

function _validate_batch_output(
    output::AbstractArray,
    metadata::OpticalPlaneMetadata,
    direction_count::Int,
)
    ndims(output) == 3 || throw(DimensionMismatchError(
        "an atmospheric OPD stack must have two pupil axes and one direction axis"))
    size(output, 1) == metadata.dimensions[1] &&
        size(output, 2) == metadata.dimensions[2] ||
        throw(DimensionMismatchError(
            "atmospheric OPD-stack pupil dimensions do not match the prepared grid"))
    size(output, 3) == direction_count || throw(DimensionMismatchError(
        "atmospheric OPD-stack direction capacity must equal the prepared direction count"))
    eltype(output) === metadata.numeric_type || throw(InvalidConfiguration(
        "atmospheric OPD-stack numeric type does not match the atmosphere"))
    typeof(backend(output)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "atmospheric OPD-stack backend does not match the prepared grid"))
    compute_device(output) == metadata.device || throw(InvalidConfiguration(
        "atmospheric OPD stack is on a different compute device"))
    return output
end

function _prepare_batch_output_metadata(
    tel::Telescope,
    output::AbstractArray{T,3},
    direction_count::Int,
) where {T<:AbstractFloat}
    n = tel.params.resolution
    size(output, 1) == n && size(output, 2) == n ||
        throw(DimensionMismatchError(
            "atmospheric OPD-stack pupil dimensions must match the telescope resolution"))
    size(output, 3) == direction_count || throw(DimensionMismatchError(
        "atmospheric OPD-stack direction capacity must equal the prepared direction count"))
    sampling = (
        T(tel.aperture.sampling_m[1]),
        T(tel.aperture.sampling_m[2]),
    )
    origin = (
        T(tel.aperture.origin_m[1]),
        T(tel.aperture.origin_m[2]),
    )
    storage = @view output[:, :, 1]
    return OpticalPlaneMetadata(
        PupilPlane(),
        storage;
        coordinate_domain=MetricCoordinates(),
        sampling=sampling,
        origin=origin,
        spectral=AchromaticSpectralCoordinate(),
        normalization=UnspecifiedNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
    )
end

function _validate_batch_layer_storage(
    layer::AbstractAtmosphereLayer,
    metadata::OpticalPlaneMetadata,
)
    screen = _direction_batch_layer_screen(layer)
    size(screen, 1) == size(screen, 2) || throw(DimensionMismatchError(
        "an atmosphere batch layer screen must be square"))
    size(screen, 1) >= metadata.dimensions[1] ||
        throw(DimensionMismatchError(
            "an atmosphere batch layer screen is smaller than the pupil grid"))
    eltype(screen) === metadata.numeric_type || throw(InvalidConfiguration(
        "an atmosphere batch layer screen has an incompatible numeric type"))
    typeof(backend(screen)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "an atmosphere batch layer screen has an incompatible backend"))
    compute_device(screen) == metadata.device || throw(InvalidConfiguration(
        "an atmosphere batch layer screen is on a different compute device"))
    sampling = _direction_batch_layer_sampling(layer)
    metadata.sampling[1] == sampling &&
        metadata.sampling[2] == sampling ||
        throw(InvalidConfiguration(
            "atmosphere layer and output pupil grid have incompatible physical sampling"))
    return screen
end

function _validate_batch_preparation(
    atm::AbstractTimedAtmosphere,
    tel::Telescope,
    output::AbstractArray,
    sources,
    ::Type{T},
) where {T<:AbstractFloat}
    isempty(sources) && throw(InvalidConfiguration(
        "an atmosphere direction batch requires at least one source"))
    layers = atmosphere_layers(atm)
    isempty(layers) && throw(UnsupportedAlgorithm(
        "$(typeof(atm)) does not expose sampled atmosphere layers for direction batching"))
    require_same_backend(atm, tel, output)
    atmosphere_numeric_type(atm) === T || throw(InvalidConfiguration(
        "atmospheric OPD-stack numeric type must match atmosphere layer storage"))
    compute_device(pupil_mask(tel)) == compute_device(output) ||
        throw(InvalidConfiguration(
            "telescope pupil and atmospheric OPD stack occupy different compute devices"))
    return layers
end

function _materialize_batch_geometry(
    output::AbstractArray{T,3},
    shift_x_host::AbstractMatrix{T},
    shift_y_host::AbstractMatrix{T},
    footprint_scale_host::AbstractMatrix{T},
) where {T<:AbstractFloat}
    geometry_dimensions = size(shift_x_host)
    shift_x = similar(output, T, geometry_dimensions...)
    shift_y = similar(output, T, geometry_dimensions...)
    footprint_scale = similar(output, T, geometry_dimensions...)
    copyto!(shift_x, shift_x_host)
    copyto!(shift_y, shift_y_host)
    copyto!(footprint_scale, footprint_scale_host)
    return shift_x, shift_y, footprint_scale
end

function _prepare_atmosphere_direction_batch(
    ::UnsupportedAtmosphereDirectionBatchCapability,
    atm::AbstractTimedAtmosphere,
    ::Telescope,
    ::Any,
    ::AbstractArray,
)
    throw(UnsupportedAlgorithm(
        "$(typeof(atm)) does not support homogeneous atmosphere-direction batching; use its single-direction renderer"))
end

function _prepare_atmosphere_direction_batch(
    capability::ExtractedScreenDirectionBatchCapability,
    atm::AbstractTimedAtmosphere,
    tel::Telescope,
    src,
    output::AbstractArray,
)
    T = atmosphere_numeric_type(atm)
    eltype(output) === T || throw(InvalidConfiguration(
        "atmospheric OPD-stack numeric type must match atmosphere layer storage"))
    ndims(output) == 3 || throw(DimensionMismatchError(
        "an atmospheric OPD stack must have two pupil axes and one direction axis"))

    sources = _freeze_batch_sources(src)
    layers = _validate_batch_preparation(atm, tel, output, sources, T)
    metadata = _prepare_batch_output_metadata(
        tel,
        output,
        length(sources),
    )
    _validate_batch_output(output, metadata, length(sources))
    @inbounds for layer in layers
        _validate_batch_layer_storage(layer, metadata)
    end

    shift_x_host, shift_y_host, footprint_scale_host =
        _prepare_batch_geometry(layers, sources, tel, T)
    shift_x, shift_y, footprint_scale = _materialize_batch_geometry(
        output,
        shift_x_host,
        shift_y_host,
        footprint_scale_host,
    )
    source_geometry_bindings =
        _copy_source_geometry_bindings(sources, T)
    source_bindings = _copy_source_bindings(sources)
    layer_bindings = _copy_layer_bindings(layers)
    pupil = similar(output, Bool, metadata.dimensions...)
    copyto!(pupil, pupil_mask(tel))
    direction_count = length(sources)
    params = AtmosphereDirectionBatchParams{T}(
        atmosphere_identity(atm),
        capability,
        sources,
        source_geometry_bindings,
        source_bindings,
        layer_bindings,
        metadata,
        direction_count,
        size(output, 3),
        length(layers),
    )
    workspace = AtmosphereDirectionBatchWorkspace(
        shift_x,
        shift_y,
        footprint_scale,
        pupil,
        output,
    )
    prepared = PreparedAtmosphereDirectionBatch(params, workspace)
    _validate_atmosphere_direction_batch_binding(prepared, atm)
    return prepared
end

"""
    prepare_atmosphere_direction_batch(atmosphere, telescope, source, output)
    prepare_atmosphere_direction_batch(atmosphere, telescope, sources, output)
    prepare_atmosphere_direction_batch(atmosphere, telescope, output)

Freeze one source, a homogeneous ordered vector of leaf sources, an ordered
`Asterism`, or an `ExtendedSource` quadrature and prepare its homogeneous
geometry for repeated current-epoch rendering into the caller-owned
three-dimensional `output`. The first two axes are the pupil-plane grid and the
third axis retains source order. The three-argument form prepares one on-axis
source direction.
"""
function prepare_atmosphere_direction_batch(
    atm::AbstractTimedAtmosphere,
    tel::Telescope,
    src::Union{AbstractSource,AbstractVector{<:AbstractSource},Nothing},
    output::AbstractArray,
)
    capability = atmosphere_direction_batch_capability(typeof(atm))
    return _prepare_atmosphere_direction_batch(
        capability,
        atm,
        tel,
        src,
        output,
    )
end

@inline function prepare_atmosphere_direction_batch(
    atm::AbstractTimedAtmosphere,
    tel::Telescope,
    output::AbstractArray,
)
    return prepare_atmosphere_direction_batch(atm, tel, nothing, output)
end

function _validate_batch_array(
    array::AbstractArray,
    dimensions::Tuple,
    ::Type{T},
    expected_backend::AbstractArrayBackend,
    expected_device::AbstractComputeDevice,
    label::AbstractString,
) where {T}
    size(array) == dimensions || throw(DimensionMismatchError(
        "$label dimensions do not match the prepared batch"))
    eltype(array) === T || throw(InvalidConfiguration(
        "$label numeric type does not match the prepared batch"))
    typeof(backend(array)) === typeof(expected_backend) ||
        throw(InvalidConfiguration(
            "$label backend does not match the prepared batch"))
    compute_device(array) == expected_device || throw(InvalidConfiguration(
        "$label is on a different compute device"))
    return array
end

function _validate_batch_sources(
    params::AtmosphereDirectionBatchParams,
)
    sources = params.sources
    length(sources) == params.direction_count || throw(
        AtmosphereEpochError(
            "prepared atmosphere direction order changed after preparation"))
    length(params.source_geometry_bindings) == params.direction_count &&
        length(params.source_bindings) == params.direction_count ||
        throw(AtmosphereEpochError(
            "prepared atmosphere source geometry changed after preparation"))
    @inbounds for direction in eachindex(params.source_bindings)
        sources[direction] === params.source_bindings[direction] || throw(
            AtmosphereEpochError(
                "prepared atmosphere source geometry or order changed after preparation"))
        _source_geometry_binding_matches(
            params.source_geometry_bindings[direction],
        ) || throw(
            AtmosphereEpochError(
                "prepared atmosphere source geometry changed after preparation"))
    end
    return sources
end

function _validate_batch_layers(
    params::AtmosphereDirectionBatchParams,
    atm::AbstractTimedAtmosphere,
)
    layers = atmosphere_layers(atm)
    length(layers) == params.layer_count || throw(AtmosphereEpochError(
        "atmosphere layer count changed after batch preparation"))
    length(params.layer_bindings) == params.layer_count || throw(
        AtmosphereEpochError(
            "prepared atmosphere layer bindings changed after preparation"))
    @inbounds for layer_index in eachindex(params.layer_bindings)
        layers[layer_index] === params.layer_bindings[layer_index] ||
            throw(AtmosphereEpochError(
                "atmosphere layer identity or order changed after batch preparation"))
        _validate_batch_layer_storage(
            layers[layer_index],
            params.output_metadata,
        )
    end
    return layers
end

function _validate_atmosphere_direction_batch_binding(
    prepared::PreparedAtmosphereDirectionBatch,
    atm::AbstractTimedAtmosphere,
)
    params = prepared.params
    workspace = prepared.workspace
    params.identity === atmosphere_identity(atm) || throw(
        AtmosphereEpochError(
            "prepared atmosphere direction batch belongs to a different atmosphere"))
    typeof(atmosphere_direction_batch_capability(typeof(atm))) ===
        typeof(params.capability) || throw(AtmosphereEpochError(
        "atmosphere direction-batch capability changed after preparation"))
    _validate_batch_sources(params)
    _validate_batch_layers(params, atm)
    metadata = params.output_metadata
    geometry_dimensions = (params.layer_count, params.direction_count)
    _validate_batch_array(
        workspace.shift_x,
        geometry_dimensions,
        metadata.numeric_type,
        metadata.backend,
        metadata.device,
        "atmosphere direction shift-x geometry",
    )
    _validate_batch_array(
        workspace.shift_y,
        geometry_dimensions,
        metadata.numeric_type,
        metadata.backend,
        metadata.device,
        "atmosphere direction shift-y geometry",
    )
    _validate_batch_array(
        workspace.footprint_scale,
        geometry_dimensions,
        metadata.numeric_type,
        metadata.backend,
        metadata.device,
        "atmosphere direction cone-scale geometry",
    )
    _validate_batch_array(
        workspace.pupil,
        metadata.dimensions,
        Bool,
        metadata.backend,
        metadata.device,
        "atmosphere direction pupil support",
    )
    _validate_batch_array(
        workspace.output,
        (
            metadata.dimensions[1],
            metadata.dimensions[2],
            params.direction_capacity,
        ),
        metadata.numeric_type,
        metadata.backend,
        metadata.device,
        "atmospheric OPD stack",
    )
    params.direction_capacity == params.direction_count || throw(
        AtmosphereEpochError(
            "atmosphere direction capacity changed after preparation"))
    return prepared
end

"""
    validate_atmosphere_direction_batch(prepared, atmosphere, epoch)

Validate epoch identity, frozen direction and layer order, shape, numeric type,
backend, and concrete compute-device residency for the whole batch. Validation
does not mutate the caller-owned output.
"""
function validate_atmosphere_direction_batch(
    prepared::PreparedAtmosphereDirectionBatch,
    atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch,
)
    _validate_epoch_identity(prepared.params.identity, atm, epoch)
    return _validate_atmosphere_direction_batch_binding(prepared, atm)
end

function _render_cpu_atmosphere_direction!(
    output::AbstractMatrix{T},
    layers,
    shift_x::AbstractMatrix{T},
    shift_y::AbstractMatrix{T},
    footprint_scale::AbstractMatrix{T},
    pupil::AbstractMatrix{Bool},
    direction::Int,
) where {T<:AbstractFloat}
    render_layer!(
        output,
        layers[1],
        shift_x[1, direction],
        shift_y[1, direction],
        footprint_scale[1, direction],
    )
    @inbounds for layer_index in 2:length(layers)
        render_layer_accumulate!(
            output,
            layers[layer_index],
            shift_x[layer_index, direction],
            shift_y[layer_index, direction],
            footprint_scale[layer_index, direction],
        )
    end
    output .*= pupil
    return output
end

function _render_atmosphere_direction_batch!(
    ::ScalarCPUStyle,
    prepared::PreparedAtmosphereDirectionBatch,
    atm::AbstractTimedAtmosphere,
)
    workspace = prepared.workspace
    output = workspace.output
    layers = atmosphere_layers(atm)
    @inbounds for direction in 1:prepared.params.direction_count
        _render_cpu_atmosphere_direction!(
            @view(output[:, :, direction]),
            layers,
            workspace.shift_x,
            workspace.shift_y,
            workspace.footprint_scale,
            workspace.pupil,
            direction,
        )
    end
    return output
end

@kernel function atmosphere_direction_layer_batch_kernel!(
    output,
    screen,
    shift_x,
    shift_y,
    footprint_scale,
    pupil,
    layer_offset_x,
    layer_offset_y,
    amplitude_scale,
    layer_index::Int,
    n::Int,
    m::Int,
    direction_count::Int,
    overwrite::Bool,
)
    i, j, direction = @index(Global, NTuple)
    if i <= n && j <= n && direction <= direction_count
        T = eltype(output)
        @inbounds begin
            scale = footprint_scale[layer_index, direction]
            start_x = T(m + 1) / T(2) -
                scale * T(n - 1) / T(2) -
                layer_offset_x +
                shift_x[layer_index, direction]
            start_y = T(m + 1) / T(2) -
                scale * T(n - 1) / T(2) -
                layer_offset_y +
                shift_y[layer_index, direction]
            start_x = normalize_start_coordinate(start_x, m)
            start_y = normalize_start_coordinate(start_y, m)

            y = start_y + scale * T(i - 1)
            y0 = unsafe_trunc(Int, floor(y))
            fy = y - T(y0)
            wy0 = one(T) - fy
            iy0 = wrap_upper_index(y0, m)
            iy1 = wrap_upper_index(y0 + 1, m)

            x = start_x + scale * T(j - 1)
            x0 = unsafe_trunc(Int, floor(x))
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_upper_index(x0, m)
            ix1 = wrap_upper_index(x0 + 1, m)

            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
            sample = amplitude_scale * (
                wy0 * (wx0 * v00 + fx * v01) +
                fy * (wx0 * v10 + fx * v11)
            )
            value = ifelse(pupil[i, j], sample, zero(T))
            if overwrite
                output[i, j, direction] = value
            else
                output[i, j, direction] += value
            end
        end
    end
end

function _queue_atmosphere_direction_layer!(
    phase::KernelLaunchPhase{<:AcceleratorStyle},
    prepared::PreparedAtmosphereDirectionBatch,
    layer::AbstractAtmosphereLayer,
    layer_index::Int,
)
    workspace = prepared.workspace
    output = workspace.output
    n = size(output, 1)
    screen = _direction_batch_layer_screen(layer)
    m = size(screen, 1)
    T = eltype(output)
    queue_kernel!(
        phase,
        atmosphere_direction_layer_batch_kernel!,
        output,
        screen,
        workspace.shift_x,
        workspace.shift_y,
        workspace.footprint_scale,
        workspace.pupil,
        T(_direction_batch_layer_offset_x(layer)),
        T(_direction_batch_layer_offset_y(layer)),
        T(_direction_batch_layer_amplitude(layer)),
        layer_index,
        n,
        m,
        prepared.params.direction_count,
        layer_index == 1;
        ndrange=(
            n,
            n,
            prepared.params.direction_count,
        ),
    )
    return nothing
end

function _render_atmosphere_direction_batch!(
    style::AcceleratorStyle,
    prepared::PreparedAtmosphereDirectionBatch,
    atm::AbstractTimedAtmosphere,
)
    phase = begin_kernel_phase(style)
    layers = atmosphere_layers(atm)
    @inbounds for layer_index in eachindex(layers)
        _queue_atmosphere_direction_layer!(
            phase,
            prepared,
            layers[layer_index],
            layer_index,
        )
    end
    finish_kernel_phase!(phase)
    return prepared.workspace.output
end

"""
    render_atmosphere_directions!(prepared, atmosphere, epoch)

Render every frozen direction from one current atmosphere epoch into the
prepared caller-owned atmospheric OPD stack. The complete batch contract is
validated before any output slice is mutated. CPU execution preserves the
single-direction serial algorithm; accelerator execution submits one
three-dimensional kernel per atmosphere layer and explicitly waits for
device-ready completion.
"""
function render_atmosphere_directions!(
    prepared::PreparedAtmosphereDirectionBatch,
    atm::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch,
)
    validate_atmosphere_direction_batch(prepared, atm, epoch)
    return _render_atmosphere_direction_batch!(
        execution_style(prepared.workspace.output),
        prepared,
        atm,
    )
end
