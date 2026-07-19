#
# Prepared native direct-imaging formation
#
# Direct imaging is an optical front end, not telescope state. Preparation
# binds explicit pupil/field and output storage; repeated execution writes one
# caller-owned focal-plane photon-arrival-rate product or a concrete bundle.
#

@kernel function apply_centering_phase_kernel!(field, phase_shift, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds field[i, j] *= cis(phase_shift * (i + j - 2))
    end
end

function apply_centering_phase!(::ScalarCPUStyle,
    field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    n, m = size(field)
    @inbounds for j in 1:m, i in 1:n
        field[i, j] *= cis(phase_shift * (i + j - 2))
    end
    return field
end

function apply_centering_phase!(style::AcceleratorStyle,
    field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    launch_kernel!(style, apply_centering_phase_kernel!, field, phase_shift,
        size(field, 1); ndrange=size(field))
    return field
end

abstract type AbstractDirectImagingInputPlan end

struct PreparedPupilImagingInput{
    F<:PupilFieldFormationPlan,
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    O<:AbstractMatrix,
} <: AbstractDirectImagingInputPlan
    formation::F
    metadata::M
    amplitude::A
    opd::O
end

struct PreparedFieldImagingInput{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
} <: AbstractDirectImagingInputPlan
    metadata::M
    values::A
end

"""
    DirectImagingPlan

Run-immutable direct-imaging contract bound to exact input, field-work, and
output storage. The stored integer shift is resolved against the declared axis
orientation during preparation and preserves the current nearest-sample,
periodic off-axis model; subpixel interpolation and finite-field loss require a
different prepared mapping.
"""
struct DirectImagingPlan{
    I<:AbstractDirectImagingInputPlan,
    FM<:OpticalPlaneMetadata,
    FA<:AbstractMatrix,
    OM<:OpticalPlaneMetadata,
    OA<:AbstractMatrix,
    P<:FraunhoferPropagation,
    R<:AbstractMatrix,
}
    input::I
    field_metadata::FM
    field_values::FA
    output_metadata::OM
    output_values::OA
    propagation::P
    unshifted_intensity::R
    shift_samples::NTuple{2,Int}
end

"""Single-writer propagation and off-axis scratch for direct imaging."""
struct DirectImagingWorkspace{
    P<:FraunhoferPropagation,
    R<:AbstractMatrix,
}
    propagation::P
    unshifted_intensity::R
end

"""
    PreparedDirectImaging

Prepared leaf direct-imaging stage. `input`, `field`, and `output` are explicit
products owned by the caller holding this value; `workspace` is single-writer.
For a preformed-field stage, `input` and `field` are the same `ElectricField`.
"""
struct PreparedDirectImaging{
    I<:Union{PupilFunction,ElectricField},
    F<:ElectricField,
    O<:IntensityMap,
    P<:DirectImagingPlan,
    W<:DirectImagingWorkspace,
}
    input::I
    field::F
    output::O
    plan::P
    workspace::W
end

"""Prepared same-grid incoherent direct-imaging composition."""
struct PreparedIncoherentDirectImaging{
    C<:AbstractVector,
    P<:AbstractVector,
    O<:IntensityMap,
    S<:PreparedIncoherentSum,
}
    components::C
    products::P
    output::O
    accumulation::S
end

"""Prepared direct-imaging composition whose physical grids remain separate."""
struct PreparedBundledDirectImaging{
    C<:AbstractVector,
    B<:OpticalProductBundle,
}
    components::C
    output::B
end

@inline direct_imaging_output(prepared::PreparedDirectImaging) = prepared.output
@inline direct_imaging_output(prepared::PreparedIncoherentDirectImaging) =
    prepared.output
@inline direct_imaging_output(prepared::PreparedBundledDirectImaging) =
    prepared.output

@inline direct_imaging_components(prepared::PreparedDirectImaging) = (prepared,)
@inline direct_imaging_components(prepared::PreparedIncoherentDirectImaging) =
    prepared.components
@inline direct_imaging_components(prepared::PreparedBundledDirectImaging) =
    prepared.components

"""
    focal_plane_pixel_scale_arcsec(map_or_metadata)

Return the scalar angular sample spacing in arcseconds for a focal-plane
product with equal sampling on both axes.
"""
@inline function focal_plane_pixel_scale_arcsec(metadata::OpticalPlaneMetadata)
    _require_direct_image_plane(metadata.kind)
    _require_direct_image_coordinates(metadata.coordinate_domain)
    metadata.sampling[1] == metadata.sampling[2] || throw(
        InvalidConfiguration(
            "scalar focal-plane pixel scale requires equal axis sampling"))
    return metadata.sampling[1] * (180 / pi) * 3600
end

@inline _handle_direct_shift_rounding_error(
    err::Exception, ::Symbol) = throw(err)

function _handle_direct_shift_rounding_error(::InexactError, axis::Symbol)
    throw(InvalidConfiguration(
        "direct-imaging $(axis) displacement is outside the supported integer sample range"))
end

@inline focal_plane_pixel_scale_arcsec(map::IntensityMap) =
    focal_plane_pixel_scale_arcsec(map.metadata)

function shift_direct_image!(out::AbstractMatrix, input::AbstractMatrix,
    shifts::NTuple{2,Int})
    circshift2d!(out, input, shifts)
    return out
end

@inline _require_direct_image_plane(::FocalPlane) = nothing
function _require_direct_image_plane(::AbstractOpticalPlaneKind)
    throw(InvalidConfiguration(
        "direct-imaging output must be declared on a focal plane"))
end

@inline _require_direct_image_coordinates(::AngularCoordinates) = nothing
function _require_direct_image_coordinates(::AbstractPlaneCoordinateDomain)
    throw(InvalidConfiguration(
        "native direct-imaging output must use angular coordinates"))
end

@inline _require_direct_image_rate(::PhotonRateNormalization) = nothing
function _require_direct_image_rate(::AbstractOpticalNormalization)
    throw(InvalidConfiguration(
        "native direct-imaging output must use photon-rate normalization"))
end

@inline _require_direct_image_measure(::CellIntegratedMeasure) = nothing
function _require_direct_image_measure(::AbstractSpatialMeasure)
    throw(InvalidConfiguration(
        "native direct-imaging output must use cell-integrated photon rate"))
end

@inline _require_direct_image_coherence(::IncoherentIntensityAddition) = nothing
function _require_direct_image_coherence(::AbstractCombinationPolicy)
    throw(InvalidConfiguration(
        "direct-imaging output must declare incoherent-intensity semantics"))
end

@inline _require_direct_image_spectral(::MonochromaticChannel) = nothing
function _require_direct_image_spectral(::AbstractSpectralCoordinate)
    throw(InvalidConfiguration(
        "native direct-imaging leaf output must declare one wavelength"))
end

function validate_direct_imaging_output(output::IntensityMap)
    metadata = validate_plane_storage(output.metadata, output.values;
        label="direct-imaging output")
    _require_direct_image_plane(metadata.kind)
    _require_direct_image_coordinates(metadata.coordinate_domain)
    _require_direct_image_rate(metadata.normalization)
    _require_direct_image_measure(metadata.spatial_measure)
    _require_direct_image_coherence(metadata.coherence)
    _require_direct_image_spectral(metadata.spectral)
    return output
end

@inline _direct_axis_coordinate(::Val{:x}, x, y) = x
@inline _direct_axis_coordinate(::Val{:y}, x, y) = y

function _direct_axis_coordinate(::Val{A}, x, y) where {A}
    throw(InvalidConfiguration(
        "direct-imaging axis orientation must use :x and :y; got $(A)"))
end

function _rounded_direct_shift(value::Real, axis::Symbol)
    isfinite(value) || throw(InvalidConfiguration(
        "direct-imaging $(axis) displacement is not finite in output samples"))
    return try
        round(Int, value)
    catch err
        _handle_direct_shift_rounding_error(err, axis)
    end
end

function _direct_imaging_shift(output::IntensityMap, src::AbstractSource)
    coordinates = coordinates_xy_arcsec(src)
    all(isfinite, coordinates) || throw(InvalidConfiguration(
        "direct-imaging source coordinates must be finite"))
    metadata = output.metadata
    orientation = metadata.orientation
    arcsec_per_radian = 180 * 3600 / pi
    return ntuple(2) do dimension
        axis = orientation.axes[dimension]
        coordinate_arcsec = _direct_axis_coordinate(Val(axis),
            coordinates[1], coordinates[2])
        scale_arcsec = metadata.sampling[dimension] * arcsec_per_radian
        displacement = orientation.signs[dimension] * coordinate_arcsec /
            scale_arcsec
        _rounded_direct_shift(displacement, axis)
    end
end

function _prepare_direct_imaging(input_plan::AbstractDirectImagingInputPlan,
    field::ElectricField, output::IntensityMap,
    propagation::FraunhoferPropagation, src::AbstractSource)
    require_same_plane_grid(output.metadata, propagation.output_metadata;
        label="direct-imaging output", require_numeric_type=false)
    require_compatible_radiometry(output.metadata,
        propagation.output_metadata; label="direct-imaging output")
    validate_direct_imaging_output(output)
    eltype(output.values) === real(eltype(field.values)) || throw(
        InvalidConfiguration(
            "direct-imaging field and output real numeric types must match"))
    plane_device(field.values) == plane_device(output.values) || throw(
        InvalidConfiguration(
            "direct-imaging field and output occupy different physical devices"))
    Base.mightalias(field.values, output.values) && throw(
        InvalidConfiguration(
            "direct-imaging field and output storage must not alias"))
    shift_samples = _direct_imaging_shift(output, src)
    unshifted = similar(output.values)
    fill!(unshifted, zero(eltype(unshifted)))
    plan = DirectImagingPlan{
        typeof(input_plan),typeof(field.metadata),typeof(field.values),
        typeof(output.metadata),typeof(output.values),
        typeof(propagation),typeof(unshifted),
    }(
        input_plan,
        field.metadata,
        field.values,
        output.metadata,
        output.values,
        propagation,
        unshifted,
        shift_samples,
    )
    workspace = DirectImagingWorkspace(propagation, unshifted)
    return (; plan, workspace)
end

"""
    prepare_direct_imaging(pupil, source, field, output)

Prepare direct imaging from an explicit `PupilFunction`, using caller-owned
`field` work storage and writing caller-owned focal-plane `output`. Physical
sources only are accepted; optical formation applies no elapsed time.
"""
function _prepare_direct_imaging(pupil::PupilFunction,
    src::Union{Source,LGSSource}, field::ElectricField,
    output::IntensityMap, propagation::FraunhoferPropagation)
    require_leaf_source(src, "direct-imaging source")
    _require_physical_photon_irradiance(src, "direct imaging")
    formation = prepare_pupil_field(pupil, src, field)
    input_plan = PreparedPupilImagingInput(formation, pupil.metadata,
        pupil.amplitude, pupil.opd)
    return _prepare_direct_imaging(input_plan, field, output, propagation,
        src)
end

function prepare_direct_imaging(pupil::PupilFunction,
    src::Union{Source,LGSSource}, field::ElectricField,
    output::IntensityMap)
    propagation = FraunhoferPropagation(field)
    return _prepare_direct_imaging(pupil, src, field, output, propagation)
end

"""
    prepare_direct_imaging(source, field, output)

Prepare direct imaging from an already formed physical pupil-plane electric
field. The field and output remain caller-owned and are bound exactly to the
returned plan.
"""
function _prepare_direct_imaging(src::Union{Source,LGSSource},
    field::ElectricField, output::IntensityMap,
    propagation::FraunhoferPropagation)
    require_leaf_source(src, "direct-imaging source")
    _require_physical_photon_irradiance(src, "direct imaging")
    validate_plane_storage(field.metadata, field.values;
        label="direct-imaging input field")
    typeof(field.metadata.kind) === PupilPlane || throw(InvalidConfiguration(
        "direct-imaging input ElectricField must be on a pupil plane"))
    _require_physical_field_normalization(field.metadata.normalization)
    _require_cell_integrated_field(field.metadata.spatial_measure)
    _require_coherent_field(field.metadata.coherence)
    field.metadata.spectral == MonochromaticChannel(
        real(eltype(field.values))(wavelength(src))) || throw(
        InvalidConfiguration(
            "direct-imaging source wavelength must match its input field"))
    input_plan = PreparedFieldImagingInput(field.metadata, field.values)
    return _prepare_direct_imaging(input_plan, field, output, propagation,
        src)
end


function prepare_direct_imaging(src::Union{Source,LGSSource},
    field::ElectricField, output::IntensityMap)
    propagation = FraunhoferPropagation(field)
    return _prepare_direct_imaging(src, field, output, propagation)
end

function _require_prepared_direct_output(output::IntensityMap,
    plan::DirectImagingPlan)
    output.metadata === plan.output_metadata || throw(InvalidConfiguration(
        "direct-imaging output does not match its prepared plan"))
    output.values === plan.output_values || throw(InvalidConfiguration(
        "direct-imaging output storage does not match its prepared plan"))
    validate_plane_storage(output.metadata, output.values;
        label="direct-imaging output")
    return nothing
end

function _require_prepared_direct_field(field::ElectricField,
    plan::DirectImagingPlan)
    field.metadata === plan.field_metadata || throw(InvalidConfiguration(
        "direct-imaging field does not match its prepared plan"))
    field.values === plan.field_values || throw(InvalidConfiguration(
        "direct-imaging field storage does not match its prepared plan"))
    validate_plane_storage(field.metadata, field.values;
        label="direct-imaging field")
    return nothing
end

function _require_prepared_direct_workspace(
    workspace::DirectImagingWorkspace, plan::DirectImagingPlan)
    workspace.propagation === plan.propagation || throw(
        InvalidConfiguration(
            "direct-imaging propagation does not match its prepared plan"))
    workspace.unshifted_intensity === plan.unshifted_intensity || throw(
        InvalidConfiguration(
            "direct-imaging scratch storage does not match its prepared plan"))
    size(workspace.unshifted_intensity) ==
        plan.output_metadata.dimensions || throw(DimensionMismatchError(
        "direct-imaging scratch dimensions do not match its prepared plan"))
    eltype(workspace.unshifted_intensity) ===
        plan.output_metadata.numeric_type || throw(InvalidConfiguration(
        "direct-imaging scratch numeric type does not match its prepared plan"))
    plane_device(workspace.unshifted_intensity) ==
        plan.output_metadata.device || throw(InvalidConfiguration(
        "direct-imaging scratch device does not match its prepared plan"))
    return nothing
end

function _prepare_direct_field!(field::ElectricField,
    pupil::PupilFunction, input::PreparedPupilImagingInput)
    pupil.metadata === input.metadata || throw(InvalidConfiguration(
        "direct-imaging pupil does not match its prepared plan"))
    pupil.amplitude === input.amplitude || throw(InvalidConfiguration(
        "direct-imaging pupil amplitude storage does not match its prepared plan"))
    pupil.opd === input.opd || throw(InvalidConfiguration(
        "direct-imaging pupil OPD storage does not match its prepared plan"))
    fill_electric_field!(field, pupil, input.formation)
    return field
end

function _prepare_direct_field!(field::ElectricField,
    input_field::ElectricField, input::PreparedFieldImagingInput)
    input_field === field || throw(InvalidConfiguration(
        "preformed direct-imaging input must be its prepared field"))
    input_field.metadata === input.metadata || throw(InvalidConfiguration(
        "preformed direct-imaging field metadata does not match its prepared plan"))
    input_field.values === input.values || throw(InvalidConfiguration(
        "preformed direct-imaging field storage does not match its prepared plan"))
    return field
end

function _form_direct_image!(output::IntensityMap, field::ElectricField,
    plan::DirectImagingPlan, workspace::DirectImagingWorkspace)
    shifts = plan.shift_samples
    if iszero(shifts[1]) && iszero(shifts[2])
        fraunhofer_intensity_from_field!(output.values, field,
            workspace.propagation)
    else
        fraunhofer_intensity_from_field!(workspace.unshifted_intensity,
            field, workspace.propagation)
        shift_direct_image!(output.values, workspace.unshifted_intensity,
            shifts)
    end
    return output
end

"""Form one prepared photon-arrival-rate direct image from a pupil."""
function form_direct_image!(output::IntensityMap, pupil::PupilFunction,
    field::ElectricField, plan::DirectImagingPlan{<:PreparedPupilImagingInput},
    workspace::DirectImagingWorkspace)
    _require_prepared_direct_output(output, plan)
    _require_prepared_direct_field(field, plan)
    _require_prepared_direct_workspace(workspace, plan)
    _prepare_direct_field!(field, pupil, plan.input)
    return _form_direct_image!(output, field, plan, workspace)
end

"""Form one prepared photon-arrival-rate direct image from a preformed field."""
function form_direct_image!(output::IntensityMap, field::ElectricField,
    plan::DirectImagingPlan{<:PreparedFieldImagingInput},
    workspace::DirectImagingWorkspace)
    _require_prepared_direct_output(output, plan)
    _require_prepared_direct_field(field, plan)
    _require_prepared_direct_workspace(workspace, plan)
    _prepare_direct_field!(field, field, plan.input)
    return _form_direct_image!(output, field, plan, workspace)
end

function prepare_direct_imaging(pupil::PupilFunction,
    src::Union{Source,LGSSource}; zero_padding::Int=1)
    field = ElectricField(pupil, src; zero_padding=zero_padding)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    prepared = _prepare_direct_imaging(pupil, src, field, output, propagation)
    return PreparedDirectImaging(pupil, field, output, prepared.plan,
        prepared.workspace)
end

function prepare_direct_imaging(src::Union{Source,LGSSource},
    field::ElectricField)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    prepared = _prepare_direct_imaging(src, field, output, propagation)
    return PreparedDirectImaging(field, field, output, prepared.plan,
        prepared.workspace)
end

@inline function _store_direct_component!(components::Vector{C}, index::Int,
    component::C) where {C}
    @inbounds components[index] = component
    return true
end

@inline _store_direct_component!(::Vector, ::Int, component) = false

function _widen_direct_components(components::Vector,
    sources::AbstractVector{<:AbstractSource}, mismatch_index::Int,
    mismatch::PreparedDirectImaging, pupil::PupilFunction,
    zero_padding::Int)
    widened = Vector{PreparedDirectImaging}(undef, length(sources))
    @inbounds for index in 1:(mismatch_index - 1)
        widened[index] = components[index]
    end
    widened[mismatch_index] = mismatch
    @inbounds for index in (mismatch_index + 1):length(sources)
        widened[index] = prepare_direct_imaging(pupil, sources[index];
            zero_padding=zero_padding)
    end
    return widened
end

function _prepare_direct_components(pupil::PupilFunction,
    sources::AbstractVector{<:AbstractSource}, zero_padding::Int)
    isempty(sources) && throw(InvalidConfiguration(
        "direct-imaging composition must contain at least one source"))
    # Preserve exact per-leaf field formation and component products here.
    # Reusing one same-pupil/common-wavelength FFT plus prepared weighted shifts
    # is a separate optimization and requires its own fidelity/GPU contract.
    first_component = prepare_direct_imaging(pupil, first(sources);
        zero_padding=zero_padding)
    components = Vector{typeof(first_component)}(undef, length(sources))
    components[1] = first_component
    @inbounds for index in 2:length(sources)
        component = prepare_direct_imaging(pupil, sources[index];
            zero_padding=zero_padding)
        if !_store_direct_component!(components, index, component)
            return _widen_direct_components(components, sources, index, component,
                pupil, zero_padding)
        end
    end
    return components
end

@inline function _store_direct_product!(products::Vector{P}, index::Int,
    product::P) where {P}
    @inbounds products[index] = product
    return true
end

@inline _store_direct_product!(::Vector, ::Int, product) = false

function _widen_direct_products(products::Vector,
    components::AbstractVector{<:PreparedDirectImaging}, mismatch_index::Int,
    mismatch::IntensityMap)
    widened = Vector{IntensityMap}(undef, length(components))
    @inbounds for index in 1:(mismatch_index - 1)
        widened[index] = products[index]
    end
    widened[mismatch_index] = mismatch
    @inbounds for index in (mismatch_index + 1):length(components)
        widened[index] = direct_imaging_output(components[index])
    end
    return widened
end

function _direct_imaging_products(
    components::AbstractVector{<:PreparedDirectImaging})
    isempty(components) && throw(InvalidConfiguration(
        "direct-imaging composition must contain at least one component"))
    first_product = direct_imaging_output(first(components))
    products = Vector{typeof(first_product)}(undef, length(components))
    products[1] = first_product
    @inbounds for index in 2:length(components)
        product = direct_imaging_output(components[index])
        if !_store_direct_product!(products, index, product)
            return _widen_direct_products(products, components, index, product)
        end
    end
    return products
end

function _similar_direct_output(reference::IntensityMap)
    values = similar(reference.values)
    fill!(values, zero(eltype(values)))
    return IntensityMap(reference.metadata, values)
end

function prepare_direct_imaging(pupil::PupilFunction, ast::Asterism;
    zero_padding::Int=1)
    isempty(ast.sources) && throw(InvalidConfiguration(
        "direct-imaging asterism must contain at least one source"))
    wavelength(ast)
    frozen = freeze_source(ast)
    components = _prepare_direct_components(pupil, frozen.sources, zero_padding)
    products = _direct_imaging_products(components)
    output = _similar_direct_output(first(products))
    accumulation = prepare_incoherent_sum(output, products)
    return PreparedIncoherentDirectImaging(components, products, output,
        accumulation)
end

function _spectral_leaf_source(reference::Union{Source,LGSSource}, sample,
    total, ::Type{T}) where {T<:AbstractFloat}
    sample_wavelength = _converted_positive_finite(sample.wavelength,
        T, "direct-imaging spectral wavelength")
    sample_value = _converted_nonnegative_finite(total * sample.weight, T,
        "direct-imaging spectral radiometric value")
    return source_with_wavelength_and_radiometric_value(reference,
        sample_wavelength, sample_value)
end

function _spectral_leaf_sources(src::SpectralSource,
    ::Type{T}) where {T<:AbstractFloat}
    reference = src.source
    total = source_radiometric_value(reference)
    samples = src.bundle.samples
    first_source = _spectral_leaf_source(reference, first(samples), total, T)
    sources = Vector{typeof(first_source)}(undef, length(samples))
    sources[1] = first_source
    @inbounds for index in 2:length(samples)
        sources[index] = _spectral_leaf_source(reference, samples[index],
            total, T)
    end
    return sources
end

function prepare_direct_imaging(pupil::PupilFunction, src::SpectralSource;
    zero_padding::Int=1)
    frozen = freeze_source(src)
    sources = _spectral_leaf_sources(frozen, eltype(pupil.opd))
    components = _prepare_direct_components(pupil, sources, zero_padding)
    products = _direct_imaging_products(components)
    output = OpticalProductBundle(products)
    return PreparedBundledDirectImaging(components, output)
end

function prepare_direct_imaging(::PupilFunction, src::AbstractSource;
    zero_padding::Int=1)
    throw(UnsupportedAlgorithm(
        "direct imaging does not support source composition $(typeof(src)); " *
        "prepare its physical components explicitly"))
end

function _require_direct_component_products(
    components::AbstractVector{<:PreparedDirectImaging},
    products::AbstractVector{<:IntensityMap})
    length(components) == length(products) || throw(DimensionMismatchError(
        "direct-imaging component and product counts must match"))
    @inbounds for index in eachindex(components, products)
        direct_imaging_output(components[index]) === products[index] || throw(
            InvalidConfiguration(
                "direct-imaging product does not match its prepared component"))
    end
    return nothing
end

function _form_direct_components!(
    components::AbstractVector{<:PreparedDirectImaging},
    products::AbstractVector{<:IntensityMap})
    _require_direct_component_products(components, products)
    @inbounds for component in components
        form_direct_image!(component)
    end
    return nothing
end

function form_direct_image!(prepared::PreparedDirectImaging)
    return _form_prepared_direct_image!(prepared, prepared.input)
end

@inline function _form_prepared_direct_image!(prepared::PreparedDirectImaging,
    input::PupilFunction)
    return form_direct_image!(prepared.output, input, prepared.field,
        prepared.plan, prepared.workspace)
end

@inline function _form_prepared_direct_image!(prepared::PreparedDirectImaging,
    input::ElectricField)
    return form_direct_image!(prepared.output, input, prepared.plan,
        prepared.workspace)
end

function form_direct_image!(prepared::PreparedIncoherentDirectImaging)
    _form_direct_components!(prepared.components, prepared.products)
    accumulate_intensity!(prepared.output, prepared.products,
        prepared.accumulation)
    return prepared.output
end

function form_direct_image!(prepared::PreparedBundledDirectImaging)
    _form_direct_components!(prepared.components, prepared.output.products)
    return prepared.output
end
