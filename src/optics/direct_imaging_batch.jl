"""
Capability trait for prepared direct-imaging batches.

Propagation models opt in only when compatible source samples can share one
stacked field workspace and one FFT plan without changing their public optical
products.
"""
abstract type AbstractDirectImagingBatchCapability end

"""Native Fraunhofer direct imaging supports stacked source samples."""
struct StackedFraunhoferDirectImagingBatchCapability <:
    AbstractDirectImagingBatchCapability end

"""The propagation model must use ordinary independent prepared paths."""
struct UnsupportedDirectImagingBatchCapability <:
    AbstractDirectImagingBatchCapability end

@inline direct_imaging_batch_capability(
    ::Type{<:AbstractPropagationModel},
) = UnsupportedDirectImagingBatchCapability()

@inline direct_imaging_batch_capability(
    ::Type{<:FraunhoferPropagation},
) = StackedFraunhoferDirectImagingBatchCapability()

"""
One ordered direct-imaging source sample in a compatibility signature.

The focal-plane sampling is wavelength dependent and therefore remains a
per-sample value even though every sample shares one FFT shape.
"""
struct DirectImagingBatchSampleParams{T<:AbstractFloat}
    coordinates_xy_arcsec::NTuple{2,T}
    wavelength_m::T
    photon_irradiance::T
    output_sampling_rad::T
    shift_samples::NTuple{2,Int}
end

"""The common detector-facing contract of every output view in a batch."""
struct DirectImagingBatchProductContract{
    K<:AbstractOpticalPlaneKind,
    C<:AbstractPlaneCoordinateDomain,
    N<:AbstractOpticalNormalization,
    M<:AbstractSpatialMeasure,
    H<:AbstractCombinationPolicy,
    Z,
    O<:PlaneAxisOrientation,
}
    kind::K
    coordinate_domain::C
    normalization::N
    spatial_measure::M
    coherence::H
    centering::Z
    orientation::O
    dimensions::NTuple{2,Int}
end

struct _FixedDirectImagingBatchVector{
    T,
    V<:Memory{T},
} <: AbstractVector{T}
    _storage::V
end

Base.size(values::_FixedDirectImagingBatchVector) =
    size(getfield(values, :_storage))
Base.axes(values::_FixedDirectImagingBatchVector) =
    axes(getfield(values, :_storage))
Base.length(values::_FixedDirectImagingBatchVector) =
    length(getfield(values, :_storage))
Base.getindex(values::_FixedDirectImagingBatchVector, index::Int) =
    @inbounds getfield(values, :_storage)[index]
Base.IndexStyle(::Type{<:_FixedDirectImagingBatchVector}) = IndexLinear()
Base.iterate(values::_FixedDirectImagingBatchVector, state...) =
    iterate(getfield(values, :_storage), state...)
Base.copy(values::_FixedDirectImagingBatchVector) =
    collect(getfield(values, :_storage))

function Base.getproperty(
    values::_FixedDirectImagingBatchVector,
    name::Symbol,
)
    name === :_storage && return collect(getfield(values, :_storage))
    return getfield(values, name)
end

"""
Immutable compatibility signature for one prepared direct-imaging batch.

Sample cardinality is deliberately a runtime value. The fixed vector wrapper
does not encode cardinality in its type, so large source registries do not
produce one method family per batch size.
"""
struct DirectImagingBatchCompatibilitySignature{
    T<:AbstractFloat,
    M<:AbstractPropagationModel,
    C<:AbstractDirectImagingBatchCapability,
    I<:OpticalPlaneMetadata,
    R<:AbstractSourceRadiometry,
    B<:AbstractArrayBackend,
    D<:AbstractComputeDevice,
    P<:DirectImagingBatchProductContract,
    S<:_FixedDirectImagingBatchVector,
}
    model_type::Type{M}
    capability::C
    input_metadata::I
    radiometry::R
    backend::B
    device::D
    product_contract::P
    samples::S
    padded_dimensions::NTuple{2,Int}
    fft_dimensions::NTuple{2,Int}
    numeric_type::Type{T}
    sample_count::Int
end

"""Single-writer stacked numerical storage and its exact FFT plan."""
mutable struct DirectImagingBatchWorkspace{
    F<:AbstractArray{<:Complex,3},
    O<:AbstractArray{<:AbstractFloat,3},
    X<:AbstractVector{Int},
    Y<:AbstractVector{Int},
    P,
}
    field_stack::F
    output_stack::O
    shift_axis1::X
    shift_axis2::Y
    fft_plan::P
end

"""
Exact run-immutable bindings for replaceable workspace fields and their
prepared field/product views.
"""
struct DirectImagingBatchWorkspaceBindings{F,O,X,Y,P,V,R}
    field_stack::F
    output_stack::O
    shift_axis1::X
    shift_axis2::Y
    fft_plan::P
    fields::V
    output::R
end

"""
Prepared native direct-imaging batch.

The fixed input, source, field-view, formation-plan, and output memberships are
separate from the mutable single-writer workspace. `output` preserves ordinary
`IntensityMap` products through `OpticalProductBundle`.
"""
struct PreparedDirectImagingBatch{
    S<:DirectImagingBatchCompatibilitySignature,
    R<:_FixedDirectImagingBatchVector,
    I<:_FixedDirectImagingBatchVector,
    F<:_FixedDirectImagingBatchVector,
    P<:_FixedDirectImagingBatchVector,
    O<:OpticalProductBundle,
    W<:DirectImagingBatchWorkspace,
    B<:DirectImagingBatchWorkspaceBindings,
}
    signature::S
    sources::R
    inputs::I
    fields::F
    formation_plans::P
    output::O
    workspace::W
    workspace_bindings::B
end

@inline direct_imaging_output(prepared::PreparedDirectImagingBatch) =
    prepared.output
@inline direct_imaging_batch_signature(
    prepared::PreparedDirectImagingBatch,
) = prepared.signature
@inline direct_imaging_batch_sources(
    prepared::PreparedDirectImagingBatch,
) = prepared.sources
@inline direct_imaging_batch_inputs(
    prepared::PreparedDirectImagingBatch,
) = prepared.inputs
@inline direct_imaging_batch_products(
    prepared::PreparedDirectImagingBatch,
) = prepared.output
@inline direct_imaging_batch_count(
    prepared::PreparedDirectImagingBatch,
) = prepared.signature.sample_count

@inline function _direct_imaging_batch_source(
    source::Source,
    ::Type{T},
) where {T<:AbstractFloat}
    _require_physical_photon_irradiance(source, "direct-imaging batch")
    return freeze_source(source)
end

function _direct_imaging_batch_source(
    source::AbstractSource,
    ::Type{<:AbstractFloat},
)
    throw(UnsupportedAlgorithm(
        "direct-imaging batching currently supports physical Source leaves; " *
        "$(typeof(source)) must use an ordinary independent prepared path",
    ))
end

function _direct_imaging_batch_sources(
    sources::AbstractVector{<:AbstractSource},
    ::Type{T},
) where {T<:AbstractFloat}
    isempty(sources) && throw(InvalidConfiguration(
        "a direct-imaging batch requires at least one source sample",
    ))
    first_source = _direct_imaging_batch_source(first(sources), T)
    P = typeof(first_source)
    frozen = Memory{P}(undef, length(sources))
    frozen[1] = first_source
    @inbounds for index in 2:length(sources)
        source = _direct_imaging_batch_source(sources[index], T)
        typeof(source) === P || throw(InvalidConfiguration(
            "direct-imaging batch sources must have one concrete numeric " *
            "and radiometric storage contract",
        ))
        frozen[index] = source
    end
    return _FixedDirectImagingBatchVector{P,typeof(frozen)}(frozen)
end

function _direct_imaging_batch_sources(
    source::Source,
    ::Type{T},
) where {T<:AbstractFloat}
    frozen_source = _direct_imaging_batch_source(source, T)
    P = typeof(frozen_source)
    frozen = Memory{P}(undef, 1)
    frozen[1] = frozen_source
    return _FixedDirectImagingBatchVector{P,typeof(frozen)}(frozen)
end

function _direct_imaging_batch_sources(
    asterism::Asterism,
    ::Type{T},
) where {T<:AbstractFloat}
    wavelength(asterism)
    frozen = freeze_source(asterism)
    return _direct_imaging_batch_sources(frozen.sources, T)
end

function _direct_imaging_batch_sources(
    source::SpectralSource{<:Source},
    ::Type{T},
) where {T<:AbstractFloat}
    frozen = freeze_source(source)
    leaves = _spectral_leaf_sources(frozen, T)
    return _direct_imaging_batch_sources(leaves, T)
end

function _direct_imaging_batch_sources(
    source::SpectralSource,
    ::Type{<:AbstractFloat},
)
    throw(UnsupportedAlgorithm(
        "direct-imaging spectral batching currently supports Source leaves; " *
        "$(typeof(source.source)) must use an ordinary independent prepared path",
    ))
end

function _direct_imaging_batch_sources(
    source::AbstractSource,
    ::Type{<:AbstractFloat},
)
    throw(UnsupportedAlgorithm(
        "direct-imaging batching does not support source composition " *
        "$(typeof(source)); prepare its physical Source leaves independently",
    ))
end

function _repeat_direct_imaging_batch_input(
    input::PupilFunction,
    count::Int,
)
    storage = Memory{typeof(input)}(undef, count)
    fill!(storage, input)
    return _FixedDirectImagingBatchVector{
        typeof(input),typeof(storage),
    }(storage)
end

function _direct_imaging_batch_inputs(
    inputs::AbstractVector{<:PupilFunction},
)
    isempty(inputs) && throw(InvalidConfiguration(
        "a direct-imaging batch requires at least one pupil input",
    ))
    P = typeof(first(inputs))
    storage = Memory{P}(undef, length(inputs))
    @inbounds for index in eachindex(inputs)
        input = inputs[index]
        typeof(input) === P || throw(InvalidConfiguration(
            "direct-imaging batch pupil inputs must have one concrete storage contract",
        ))
        storage[index] = input
    end
    return _FixedDirectImagingBatchVector{P,typeof(storage)}(storage)
end

@inline function _direct_imaging_batch_product_contract(
    metadata::OpticalPlaneMetadata,
)
    return DirectImagingBatchProductContract(
        metadata.kind,
        metadata.coordinate_domain,
        metadata.normalization,
        metadata.spatial_measure,
        metadata.coherence,
        metadata.centering,
        metadata.orientation,
        metadata.dimensions,
    )
end

function _require_direct_imaging_batch_product_contract(
    reference::DirectImagingBatchProductContract,
    metadata::OpticalPlaneMetadata,
)
    candidate = _direct_imaging_batch_product_contract(metadata)
    candidate == reference || throw(InvalidConfiguration(
        "direct-imaging batch outputs have incompatible product contracts",
    ))
    return metadata
end

function _require_direct_imaging_batch_input_grid(
    reference::PupilFunction,
    input::PupilFunction,
)
    require_same_plane_grid(
        reference.metadata,
        input.metadata;
        label="direct-imaging batch pupil inputs",
    )
    typeof(backend(input)) === typeof(backend(reference)) || throw(
        InvalidConfiguration(
            "direct-imaging batch pupil inputs use different backends",
        ),
    )
    input.metadata.device == reference.metadata.device || throw(
        InvalidConfiguration(
            "direct-imaging batch pupil inputs occupy different compute devices",
        ),
    )
    eltype(input.opd) === eltype(reference.opd) || throw(
        InvalidConfiguration(
            "direct-imaging batch pupil inputs use different numeric types",
        ),
    )
    return input
end

function _prepare_direct_imaging_batch_component(
    input::PupilFunction,
    source::Source,
    field_values::AbstractMatrix{Complex{T}},
    output_values::AbstractMatrix{T},
) where {T<:AbstractFloat}
    field_metadata = _pupil_field_metadata(
        input,
        source,
        field_values,
        source_field_normalization(source_radiometry(source)),
        source_field_measure(source_radiometry(source)),
        CoherentFieldCombination(),
    )
    field = ElectricField(field_metadata, field_values)
    formation_plan = prepare_pupil_field(input, source, field)
    output_sampling = _fraunhofer_output_sampling(field)
    output_metadata = _fraunhofer_output_metadata(
        field,
        output_values,
        output_sampling,
        IncoherentIntensityAddition(),
    )
    output = IntensityMap(output_metadata, output_values)
    validate_direct_imaging_output(output)
    shift_samples = _direct_imaging_shift(output, source)
    sample = DirectImagingBatchSampleParams{T}(
        (
            T(coordinates_xy_arcsec(source)[1]),
            T(coordinates_xy_arcsec(source)[2]),
        ),
        T(wavelength(source)),
        T(photon_irradiance(source)),
        T(output_sampling),
        shift_samples,
    )
    return (; field, formation_plan, output, sample)
end

function _prepare_direct_imaging_batch(
    ::UnsupportedDirectImagingBatchCapability,
    model_type::Type{<:AbstractPropagationModel},
    ::_FixedDirectImagingBatchVector,
    ::_FixedDirectImagingBatchVector;
    zero_padding::Int,
)
    throw(UnsupportedAlgorithm(
        "$(model_type) does not support direct-imaging batching; " *
        "prepare ordinary independent optical paths",
    ))
end

function _prepare_direct_imaging_batch(
    capability::StackedFraunhoferDirectImagingBatchCapability,
    model_type::Type{<:FraunhoferPropagation},
    inputs::_FixedDirectImagingBatchVector,
    sources::_FixedDirectImagingBatchVector;
    zero_padding::Int,
)
    zero_padding >= 1 || throw(InvalidConfiguration(
        "zero_padding must be >= 1",
    ))
    count = length(sources)
    length(inputs) == count || throw(DimensionMismatchError(
        "direct-imaging batch input and source counts must match",
    ))

    reference_input = first(inputs)
    T = eltype(reference_input.opd)
    T <: AbstractFloat || throw(InvalidConfiguration(
        "direct-imaging batch pupil OPD storage must use a floating-point type",
    ))
    n = reference_input.metadata.dimensions[1]
    reference_input.metadata.dimensions[2] == n || throw(
        DimensionMismatchError(
            "direct-imaging batch pupil inputs must be square",
        ),
    )
    padded_resolution = n * zero_padding
    field_stack = similar(
        reference_input.opd,
        Complex{T},
        padded_resolution,
        padded_resolution,
        count,
    )
    output_stack = similar(
        reference_input.opd,
        T,
        padded_resolution,
        padded_resolution,
        count,
    )
    fill!(field_stack, zero(eltype(field_stack)))
    fill!(output_stack, zero(T))

    first_component = _prepare_direct_imaging_batch_component(
        reference_input,
        first(sources),
        @view(field_stack[:, :, 1]),
        @view(output_stack[:, :, 1]),
    )
    Field = typeof(first_component.field)
    Formation = typeof(first_component.formation_plan)
    Product = typeof(first_component.output)
    fields = Memory{Field}(undef, count)
    formation_plans = Memory{Formation}(undef, count)
    products = Vector{Product}(undef, count)
    samples = Memory{DirectImagingBatchSampleParams{T}}(undef, count)
    shift_axis1_host = Vector{Int}(undef, count)
    shift_axis2_host = Vector{Int}(undef, count)
    fields[1] = first_component.field
    formation_plans[1] = first_component.formation_plan
    products[1] = first_component.output
    samples[1] = first_component.sample
    shift_axis1_host[1] = mod(
        first_component.sample.shift_samples[1],
        padded_resolution,
    )
    shift_axis2_host[1] = mod(
        first_component.sample.shift_samples[2],
        padded_resolution,
    )
    product_contract = _direct_imaging_batch_product_contract(
        first_component.output.metadata,
    )

    @inbounds for index in 2:count
        input = inputs[index]
        _require_direct_imaging_batch_input_grid(reference_input, input)
        component = _prepare_direct_imaging_batch_component(
            input,
            sources[index],
            @view(field_stack[:, :, index]),
            @view(output_stack[:, :, index]),
        )
        typeof(component.field) === Field || throw(InvalidConfiguration(
            "direct-imaging batch fields have incompatible storage contracts",
        ))
        typeof(component.formation_plan) === Formation || throw(
            InvalidConfiguration(
                "direct-imaging batch formation plans are incompatible",
            ),
        )
        typeof(component.output) === Product || throw(InvalidConfiguration(
            "direct-imaging batch outputs have incompatible storage contracts",
        ))
        _require_direct_imaging_batch_product_contract(
            product_contract,
            component.output.metadata,
        )
        fields[index] = component.field
        formation_plans[index] = component.formation_plan
        products[index] = component.output
        samples[index] = component.sample
        shift_axis1_host[index] = mod(
            component.sample.shift_samples[1],
            padded_resolution,
        )
        shift_axis2_host[index] = mod(
            component.sample.shift_samples[2],
            padded_resolution,
        )
    end

    shift_axis1 = similar(output_stack, Int, count)
    shift_axis2 = similar(output_stack, Int, count)
    copyto!(shift_axis1, shift_axis1_host)
    copyto!(shift_axis2, shift_axis2_host)
    fft_plan = plan_fft_backend!(field_stack, (1, 2))
    workspace = DirectImagingBatchWorkspace(
        field_stack,
        output_stack,
        shift_axis1,
        shift_axis2,
        fft_plan,
    )
    fixed_samples = _FixedDirectImagingBatchVector{
        DirectImagingBatchSampleParams{T},typeof(samples),
    }(samples)
    signature = DirectImagingBatchCompatibilitySignature(
        model_type,
        capability,
        reference_input.metadata,
        PhysicalPhotonIrradianceSource(),
        backend(reference_input),
        reference_input.metadata.device,
        product_contract,
        fixed_samples,
        (padded_resolution, padded_resolution),
        (1, 2),
        T,
        count,
    )
    fixed_fields =
        _FixedDirectImagingBatchVector{Field,typeof(fields)}(fields)
    fixed_formation_plans = _FixedDirectImagingBatchVector{
        Formation,typeof(formation_plans),
    }(formation_plans)
    output = OpticalProductBundle(products)
    workspace_bindings = DirectImagingBatchWorkspaceBindings(
        field_stack,
        output_stack,
        shift_axis1,
        shift_axis2,
        fft_plan,
        fixed_fields,
        output,
    )
    prepared = PreparedDirectImagingBatch(
        signature,
        sources,
        inputs,
        fixed_fields,
        fixed_formation_plans,
        output,
        workspace,
        workspace_bindings,
    )
    validate_direct_imaging_batch(prepared)
    return prepared
end

function _prepare_direct_imaging_batch(
    model_type::Type{<:AbstractPropagationModel},
    inputs::_FixedDirectImagingBatchVector,
    sources::_FixedDirectImagingBatchVector;
    zero_padding::Int,
)
    capability = direct_imaging_batch_capability(model_type)
    return _prepare_direct_imaging_batch(
        capability,
        model_type,
        inputs,
        sources;
        zero_padding=zero_padding,
    )
end

"""
    prepare_direct_imaging_batch([model_type], pupil, source; zero_padding=1)
    prepare_direct_imaging_batch([model_type], pupils, sources; zero_padding=1)

Prepare one fixed native direct-science batch. A `Source`, same-wavelength
ordered `Asterism`, or `SpectralSource{Source}` may share one pupil input. The
plural form binds one exact path-local pupil to each physical `Source` leaf.
The optional model type is an advanced capability seam and defaults to
`FraunhoferPropagation`.
"""
function prepare_direct_imaging_batch(
    model_type::Type{<:AbstractPropagationModel},
    input::PupilFunction,
    source::AbstractSource;
    zero_padding::Int=1,
)
    T = eltype(input.opd)
    sources = _direct_imaging_batch_sources(source, T)
    inputs = _repeat_direct_imaging_batch_input(input, length(sources))
    return _prepare_direct_imaging_batch(
        model_type,
        inputs,
        sources;
        zero_padding=zero_padding,
    )
end

function prepare_direct_imaging_batch(
    model_type::Type{<:AbstractPropagationModel},
    inputs::AbstractVector{<:PupilFunction},
    sources::AbstractVector{<:AbstractSource};
    zero_padding::Int=1,
)
    fixed_inputs = _direct_imaging_batch_inputs(inputs)
    T = eltype(first(fixed_inputs).opd)
    fixed_sources = _direct_imaging_batch_sources(sources, T)
    return _prepare_direct_imaging_batch(
        model_type,
        fixed_inputs,
        fixed_sources;
        zero_padding=zero_padding,
    )
end

@inline function prepare_direct_imaging_batch(
    input::PupilFunction,
    source::AbstractSource;
    zero_padding::Int=1,
)
    return prepare_direct_imaging_batch(
        FraunhoferPropagation,
        input,
        source;
        zero_padding=zero_padding,
    )
end

@inline function prepare_direct_imaging_batch(
    inputs::AbstractVector{<:PupilFunction},
    sources::AbstractVector{<:AbstractSource};
    zero_padding::Int=1,
)
    return prepare_direct_imaging_batch(
        FraunhoferPropagation,
        inputs,
        sources;
        zero_padding=zero_padding,
    )
end

function _validate_direct_imaging_batch_array(
    array::AbstractArray,
    dimensions::Tuple,
    ::Type{T},
    expected_backend::AbstractArrayBackend,
    expected_device::AbstractComputeDevice,
    label::AbstractString,
) where {T}
    size(array) == dimensions || throw(DimensionMismatchError(
        "$label dimensions do not match the prepared direct-imaging batch",
    ))
    eltype(array) === T || throw(InvalidConfiguration(
        "$label numeric type does not match the prepared direct-imaging batch",
    ))
    typeof(backend(array)) === typeof(expected_backend) || throw(
        InvalidConfiguration(
            "$label backend does not match the prepared direct-imaging batch",
        ),
    )
    compute_device(array) == expected_device || throw(InvalidConfiguration(
        "$label occupies a different compute device",
    ))
    return array
end

function _validate_direct_imaging_batch_workspace(
    prepared::PreparedDirectImagingBatch,
)
    signature = prepared.signature
    workspace = prepared.workspace
    bindings = prepared.workspace_bindings
    workspace.field_stack === bindings.field_stack || throw(
        InvalidConfiguration(
            "direct-imaging batch field stack does not match its prepared binding",
        ),
    )
    workspace.output_stack === bindings.output_stack || throw(
        InvalidConfiguration(
            "direct-imaging batch output stack does not match its prepared binding",
        ),
    )
    workspace.shift_axis1 === bindings.shift_axis1 || throw(
        InvalidConfiguration(
            "direct-imaging batch axis-1 shifts do not match their prepared binding",
        ),
    )
    workspace.shift_axis2 === bindings.shift_axis2 || throw(
        InvalidConfiguration(
            "direct-imaging batch axis-2 shifts do not match their prepared binding",
        ),
    )
    workspace.fft_plan === bindings.fft_plan || throw(InvalidConfiguration(
        "direct-imaging batch FFT plan does not match its prepared binding",
    ))
    prepared.fields === bindings.fields || throw(InvalidConfiguration(
        "direct-imaging batch field views do not match their prepared binding",
    ))
    prepared.output === bindings.output || throw(InvalidConfiguration(
        "direct-imaging batch product views do not match their prepared binding",
    ))
    n1, n2 = signature.padded_dimensions
    count = signature.sample_count
    T = signature.numeric_type
    _validate_direct_imaging_batch_array(
        workspace.field_stack,
        (n1, n2, count),
        Complex{T},
        signature.backend,
        signature.device,
        "direct-imaging batch field stack",
    )
    _validate_direct_imaging_batch_array(
        workspace.output_stack,
        (n1, n2, count),
        T,
        signature.backend,
        signature.device,
        "direct-imaging batch output stack",
    )
    _validate_direct_imaging_batch_array(
        workspace.shift_axis1,
        (count,),
        Int,
        signature.backend,
        signature.device,
        "direct-imaging batch axis-1 shifts",
    )
    _validate_direct_imaging_batch_array(
        workspace.shift_axis2,
        (count,),
        Int,
        signature.backend,
        signature.device,
        "direct-imaging batch axis-2 shifts",
    )
    return workspace
end

function _validate_direct_imaging_batch_components(
    prepared::PreparedDirectImagingBatch,
)
    signature = prepared.signature
    count = signature.sample_count
    length(prepared.sources) == count &&
        length(prepared.inputs) == count &&
        length(prepared.fields) == count &&
        length(prepared.formation_plans) == count &&
        length(prepared.output) == count ||
        throw(DimensionMismatchError(
            "direct-imaging batch membership changed after preparation",
        ))
    reference_input = prepared.inputs[1]
    reference_input.metadata === signature.input_metadata || throw(
        InvalidConfiguration(
            "direct-imaging batch reference input metadata changed after preparation",
        ),
    )
    eltype(reference_input.opd) === signature.numeric_type || throw(
        InvalidConfiguration(
            "direct-imaging batch input numeric type does not match its signature",
        ),
    )
    typeof(backend(reference_input)) === typeof(signature.backend) || throw(
        InvalidConfiguration(
            "direct-imaging batch input backend does not match its signature",
        ),
    )
    reference_input.metadata.device == signature.device || throw(
        InvalidConfiguration(
            "direct-imaging batch input device does not match its signature",
        ),
    )
    @inbounds for index in 1:count
        source = prepared.sources[index]
        input = prepared.inputs[index]
        field = prepared.fields[index]
        formation = prepared.formation_plans[index]
        product = prepared.output[index]
        sample = signature.samples[index]
        typeof(source_radiometry(source)) ===
            typeof(signature.radiometry) || throw(InvalidConfiguration(
            "direct-imaging batch source radiometry changed after preparation",
        ))
        coordinates = coordinates_xy_arcsec(source)
        (
            signature.numeric_type(coordinates[1]),
            signature.numeric_type(coordinates[2]),
        ) == sample.coordinates_xy_arcsec || throw(InvalidConfiguration(
            "direct-imaging batch source coordinates changed after preparation",
        ))
        signature.numeric_type(wavelength(source)) ==
            sample.wavelength_m || throw(InvalidConfiguration(
            "direct-imaging batch source wavelength changed after preparation",
        ))
        signature.numeric_type(photon_irradiance(source)) ==
            sample.photon_irradiance || throw(InvalidConfiguration(
            "direct-imaging batch source photon irradiance changed after preparation",
        ))
        input.metadata === formation.input_metadata || throw(
            InvalidConfiguration(
                "direct-imaging batch input does not match its formation plan",
            ),
        )
        field.metadata === formation.output_metadata || throw(
            InvalidConfiguration(
                "direct-imaging batch field does not match its formation plan",
            ),
        )
        _require_direct_imaging_batch_input_grid(
            reference_input,
            input,
        )
        validate_plane_storage(
            field.metadata,
            field.values;
            label="direct-imaging batch field",
        )
        validate_direct_imaging_output(product)
        _require_direct_imaging_batch_product_contract(
            signature.product_contract,
            product.metadata,
        )
        product.metadata.spectral ==
            MonochromaticChannel(sample.wavelength_m) || throw(
            InvalidConfiguration(
                "direct-imaging batch product wavelength changed after preparation",
            ),
        )
        product.metadata.sampling ==
            (sample.output_sampling_rad, sample.output_sampling_rad) ||
            throw(InvalidConfiguration(
                "direct-imaging batch product sampling changed after preparation",
            ))
    end
    return prepared
end

"""
    validate_direct_imaging_batch(prepared)

Validate fixed membership, input and output contracts, exact workspace
bindings, backend family, and concrete compute-device residency for the whole
batch before numerical storage is mutated.
"""
function validate_direct_imaging_batch(
    prepared::PreparedDirectImagingBatch,
)
    signature = prepared.signature
    signature.sample_count > 0 || throw(InvalidConfiguration(
        "a direct-imaging batch requires at least one source sample",
    ))
    length(signature.samples) == signature.sample_count || throw(
        DimensionMismatchError(
            "direct-imaging batch sample signature count changed after preparation",
        ),
    )
    typeof(signature.radiometry) === PhysicalPhotonIrradianceSource ||
        throw(InvalidConfiguration(
            "direct-imaging batch signature must use physical photon irradiance",
        ))
    signature.fft_dimensions == (1, 2) || throw(InvalidConfiguration(
        "direct-imaging batch FFT dimensions changed after preparation",
    ))
    capability = direct_imaging_batch_capability(signature.model_type)
    typeof(capability) === typeof(signature.capability) || throw(
        InvalidConfiguration(
        "direct-imaging batch model capability changed after preparation",
        ),
    )
    _validate_direct_imaging_batch_workspace(prepared)
    _validate_direct_imaging_batch_components(prepared)
    return prepared
end

function _require_exact_direct_imaging_target(
    prepared::PreparedDirectImagingBatch,
    target::AbstractComputeDevice,
)
    signature_device = prepared.signature.device
    signature_device == target || _throw_wrong_direct_imaging_target(
        target,
        "direct-imaging batch signature",
        signature_device,
    )
    validate_direct_imaging_batch(prepared)

    @inbounds for input in prepared.inputs
        _require_exact_direct_imaging_pupil_target(input, target)
    end
    @inbounds for field in prepared.fields
        _require_exact_direct_imaging_metadata_target(
            field.metadata, target, "direct-imaging batch field")
        _require_exact_direct_imaging_array_target(
            field.values, target, "direct-imaging batch field storage")
    end
    @inbounds for product in prepared.output
        _require_exact_direct_imaging_metadata_target(
            product.metadata, target, "direct-imaging batch output")
        _require_exact_direct_imaging_array_target(
            product.values, target,
            "direct-imaging batch output storage")
    end

    workspace = prepared.workspace
    _require_exact_direct_imaging_array_target(
        workspace.field_stack, target,
        "direct-imaging batch field stack")
    _require_exact_direct_imaging_array_target(
        workspace.output_stack, target,
        "direct-imaging batch output stack")
    _require_exact_direct_imaging_array_target(
        workspace.shift_axis1, target,
        "direct-imaging batch axis-1 shifts")
    _require_exact_direct_imaging_array_target(
        workspace.shift_axis2, target,
        "direct-imaging batch axis-2 shifts")
    return prepared
end

@kernel function direct_imaging_batch_intensity_kernel!(
    output,
    field_stack,
    shift_axis1,
    shift_axis2,
    intensity_scale,
    n1::Int,
    n2::Int,
    count::Int,
)
    i, j, sample = @index(Global, NTuple)
    if i <= n1 && j <= n2 && sample <= count
        @inbounds begin
            source_i = mod(i - shift_axis1[sample] - 1, n1) + 1
            source_j = mod(j - shift_axis2[sample] - 1, n2) + 1
            output[i, j, sample] =
                abs2(field_stack[source_i, source_j, sample]) *
                intensity_scale
        end
    end
end

function _direct_imaging_batch_intensity!(
    ::ScalarCPUStyle,
    prepared::PreparedDirectImagingBatch,
)
    workspace = prepared.workspace
    output = workspace.output_stack
    field_stack = workspace.field_stack
    n1, n2, count = size(output)
    T = eltype(output)
    intensity_scale = inv(T(n1) * T(n2))
    @inbounds for sample in 1:count
        shift1 = workspace.shift_axis1[sample]
        shift2 = workspace.shift_axis2[sample]
        for j in 1:n2
            source_j = mod(j - shift2 - 1, n2) + 1
            for i in 1:n1
                source_i = mod(i - shift1 - 1, n1) + 1
                output[i, j, sample] =
                    abs2(field_stack[source_i, source_j, sample]) *
                    intensity_scale
            end
        end
    end
    return output
end

function _direct_imaging_batch_intensity!(
    style::AcceleratorStyle,
    prepared::PreparedDirectImagingBatch,
)
    workspace = prepared.workspace
    output = workspace.output_stack
    n1, n2, count = size(output)
    T = eltype(output)
    phase = begin_kernel_phase(style)
    queue_kernel!(
        phase,
        direct_imaging_batch_intensity_kernel!,
        output,
        workspace.field_stack,
        workspace.shift_axis1,
        workspace.shift_axis2,
        inv(T(n1) * T(n2)),
        n1,
        n2,
        count;
        ndrange=size(output),
    )
    finish_kernel_phase!(phase)
    return output
end

function _form_direct_imaging_batch!(
    ::ScalarCPUStyle,
    prepared::PreparedDirectImagingBatch,
)
    @inbounds for index in 1:prepared.signature.sample_count
        fill_electric_field!(
            prepared.fields[index],
            prepared.inputs[index],
            prepared.formation_plans[index],
        )
    end
    execute_fft_plan!(
        prepared.workspace.field_stack,
        prepared.workspace.fft_plan,
    )
    _direct_imaging_batch_intensity!(ScalarCPUStyle(), prepared)
    return prepared.output
end

function _form_direct_imaging_batch!(
    style::AcceleratorStyle,
    prepared::PreparedDirectImagingBatch,
)
    @inbounds for index in 1:prepared.signature.sample_count
        fill_electric_field_async!(
            prepared.fields[index],
            prepared.inputs[index],
            prepared.formation_plans[index],
        )
    end
    synchronize_backend!(style)
    execute_fft_plan!(
        prepared.workspace.field_stack,
        prepared.workspace.fft_plan,
    )
    _direct_imaging_batch_intensity!(style, prepared)
    return prepared.output
end

"""
    form_direct_image!(prepared::PreparedDirectImagingBatch)

Validate and form every ordered direct-science product. Compatible samples
share one prepared FFT submission. Accelerator return establishes
device-ready completion; it performs no host observation or transfer.
"""
function form_direct_image!(prepared::PreparedDirectImagingBatch)
    validate_direct_imaging_batch(prepared)
    return _form_direct_imaging_batch!(
        execution_style(prepared.workspace.field_stack),
        prepared,
    )
end
