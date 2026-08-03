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

struct PupilDirectImagingInputPlan{F<:PupilFieldFormationPlan}
    formation::F
end

struct PreformedFieldDirectImagingStrategy end

"""
    DirectImagingPlan

Reusable run-immutable direct-imaging contract. Exact input, field-work,
output, and workspace identities belong to `PreparedDirectImaging`. The stored
integer shift is resolved against the declared axis orientation during
preparation and preserves the current nearest-sample, periodic off-axis model;
subpixel interpolation and finite-field loss require a different plan.
"""
struct DirectImagingPlan{
    I<:Union{PupilDirectImagingInputPlan,PreformedFieldDirectImagingStrategy},
    FM<:OpticalPlaneMetadata,
    OM<:OpticalPlaneMetadata,
    P<:FraunhoferPropagationPlan,
}
    input_plan::I
    field_metadata::FM
    output_metadata::OM
    propagation::P
    shift_samples::NTuple{2,Int}
end

"""Single-writer propagation and off-axis scratch for direct imaging."""
struct DirectImagingWorkspace{
    P<:FraunhoferPropagationWorkspace,
    R<:AbstractMatrix,
}
    propagation::P
    unshifted_intensity::R
end

struct _PreparedDirectImagingToken end
const _PREPARED_DIRECT_IMAGING_TOKEN = _PreparedDirectImagingToken()

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

    function PreparedDirectImaging(
        ::_PreparedDirectImagingToken,
        input::I,
        field::F,
        output::O,
        plan::P,
        workspace::W,
    ) where {
        I<:Union{PupilFunction,ElectricField},
        F<:ElectricField,
        O<:IntensityMap,
        P<:DirectImagingPlan,
        W<:DirectImagingWorkspace,
    }
        return new{I,F,O,P,W}(input, field, output, plan, workspace)
    end
end

struct PreparedDirectImagingCompositionBindings{
    C<:FixedSizeVector,
    P<:FixedSizeVector,
}
    components::C
    products::P
end

struct _PreparedDirectImagingCompositionToken end
const _PREPARED_DIRECT_IMAGING_COMPOSITION_TOKEN =
    _PreparedDirectImagingCompositionToken()

"""Prepared same-grid incoherent direct-imaging composition."""
struct PreparedIncoherentDirectImaging{
    C<:FixedSizeVector,
    P<:FixedSizeVector,
    O<:IntensityMap,
    S<:PreparedIncoherentSum,
    B<:PreparedDirectImagingCompositionBindings,
}
    components::C
    products::P
    output::O
    accumulation::S
    bindings::B

    function PreparedIncoherentDirectImaging(
        ::_PreparedDirectImagingCompositionToken,
        components::C,
        products::P,
        output::O,
        accumulation::S,
        bindings::B,
    ) where {
        C<:FixedSizeVector,
        P<:FixedSizeVector,
        O<:IntensityMap,
        S<:PreparedIncoherentSum,
        B<:PreparedDirectImagingCompositionBindings,
    }
        return new{C,P,O,S,B}(
            components, products, output, accumulation, bindings)
    end
end

"""Prepared direct-imaging composition whose physical grids remain separate."""
struct PreparedBundledDirectImaging{
    C<:FixedSizeVector,
    O<:OpticalProductBundle,
    B<:PreparedDirectImagingCompositionBindings,
}
    components::C
    output::O
    bindings::B

    function PreparedBundledDirectImaging(
        ::_PreparedDirectImagingCompositionToken,
        components::C,
        output::O,
        bindings::B,
    ) where {
        C<:FixedSizeVector,
        O<:OpticalProductBundle,
        B<:PreparedDirectImagingCompositionBindings,
    }
        return new{C,O,B}(components, output, bindings)
    end
end


"""Return the reusable plan from a prepared direct-imaging leaf."""
@inline direct_imaging_plan(prepared::PreparedDirectImaging) = prepared.plan

"""Return the single-writer workspace from a prepared direct-imaging leaf."""
@inline direct_imaging_workspace(prepared::PreparedDirectImaging) =
    prepared.workspace

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

function _prepare_direct_imaging(
    input_plan::Union{PupilDirectImagingInputPlan,PreformedFieldDirectImagingStrategy},
    field::ElectricField, output::IntensityMap,
    propagation::FraunhoferPropagation, src::AbstractSource)
    propagation_output = propagation_output_metadata(propagation)
    require_same_plane_grid(output.metadata, propagation_output;
        label="direct-imaging output", require_numeric_type=false)
    require_compatible_radiometry(output.metadata,
        propagation_output; label="direct-imaging output")
    validate_direct_imaging_output(output)
    eltype(output.values) === real(eltype(field.values)) || throw(
        InvalidConfiguration(
            "direct-imaging field and output real numeric types must match"))
    compute_device(field.values) == compute_device(output.values) || throw(
        InvalidConfiguration(
            "direct-imaging field and output occupy different compute devices"))
    Base.mightalias(field.values, output.values) && throw(
        InvalidConfiguration(
            "direct-imaging field and output storage must not alias"))
    shift_samples = _direct_imaging_shift(output, src)
    unshifted = similar(output.values)
    fill!(unshifted, zero(eltype(unshifted)))
    plan = DirectImagingPlan{
        typeof(input_plan),typeof(field.metadata),typeof(output.metadata),
        typeof(propagation.plan),
    }(
        input_plan,
        field.metadata,
        output.metadata,
        propagation.plan,
        shift_samples,
    )
    workspace = DirectImagingWorkspace(propagation.workspace, unshifted)
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
    input_plan = PupilDirectImagingInputPlan(formation)
    prepared = _prepare_direct_imaging(input_plan, field, output,
        propagation, src)
    return PreparedDirectImaging(
        _PREPARED_DIRECT_IMAGING_TOKEN,
        pupil,
        field,
        output,
        prepared.plan,
        prepared.workspace,
    )
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
returned `PreparedDirectImaging` owner.
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
    input_plan = PreformedFieldDirectImagingStrategy()
    prepared = _prepare_direct_imaging(input_plan, field, output,
        propagation, src)
    return PreparedDirectImaging(
        _PREPARED_DIRECT_IMAGING_TOKEN,
        field,
        field,
        output,
        prepared.plan,
        prepared.workspace,
    )
end


function prepare_direct_imaging(src::Union{Source,LGSSource},
    field::ElectricField, output::IntensityMap)
    propagation = FraunhoferPropagation(field)
    return _prepare_direct_imaging(src, field, output, propagation)
end

function _require_prepared_direct_output(output::IntensityMap,
    plan::DirectImagingPlan)
    output.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "direct-imaging output does not match its prepared plan"))
    validate_plane_storage(output.metadata, output.values;
        label="direct-imaging output")
    return nothing
end

function _require_prepared_direct_field(field::ElectricField,
    plan::DirectImagingPlan)
    field.metadata == plan.field_metadata || throw(InvalidConfiguration(
        "direct-imaging field does not match its prepared plan"))
    validate_plane_storage(field.metadata, field.values;
        label="direct-imaging field")
    return nothing
end

function _require_prepared_direct_workspace(
    workspace::DirectImagingWorkspace, plan::DirectImagingPlan)
    _require_propagation_plan_workspace(plan.propagation,
        workspace.propagation)
    size(workspace.unshifted_intensity) ==
        plan.output_metadata.dimensions || throw(DimensionMismatchError(
        "direct-imaging scratch dimensions do not match its prepared plan"))
    eltype(workspace.unshifted_intensity) ===
        plan.output_metadata.numeric_type || throw(InvalidConfiguration(
        "direct-imaging scratch numeric type does not match its prepared plan"))
    compute_device(workspace.unshifted_intensity) ==
        plan.output_metadata.device || throw(InvalidConfiguration(
        "direct-imaging scratch device does not match its prepared plan"))
    typeof(backend(workspace.unshifted_intensity)) ===
        typeof(plan.output_metadata.backend) || throw(InvalidConfiguration(
        "direct-imaging scratch backend does not match its prepared plan"))
    return nothing
end

@inline function _require_prepared_direct_input(
    pupil::PupilFunction,
    field::ElectricField,
    input_plan::PupilDirectImagingInputPlan,
)
    pupil.metadata == input_plan.formation.input_metadata || throw(
        InvalidConfiguration(
        "direct-imaging pupil does not match its prepared plan"))
    field.metadata == input_plan.formation.output_metadata || throw(
        InvalidConfiguration(
            "direct-imaging field does not match its formation plan"))
    return nothing
end

@inline function _require_prepared_direct_input(
    input_field::ElectricField,
    field::ElectricField,
    ::PreformedFieldDirectImagingStrategy,
)
    input_field === field || throw(InvalidConfiguration(
        "preformed direct-imaging input must be its prepared field"))
    return nothing
end

function _require_prepared_direct_bindings(
    prepared::PreparedDirectImaging,
    input::Union{PupilFunction,ElectricField},
    field::ElectricField,
    output::IntensityMap,
    workspace::DirectImagingWorkspace,
)
    input === prepared.input || throw(InvalidConfiguration(
        "direct-imaging input does not match its prepared owner"))
    field === prepared.field || throw(InvalidConfiguration(
        "direct-imaging field does not match its prepared owner"))
    output === prepared.output || throw(InvalidConfiguration(
        "direct-imaging output does not match its prepared owner"))
    workspace === prepared.workspace || throw(InvalidConfiguration(
        "direct-imaging workspace does not match its prepared owner"))
    return nothing
end

@noinline function _throw_wrong_direct_imaging_target(
    target::AbstractComputeDevice,
    label::AbstractString,
    actual::AbstractComputeDevice,
)
    _throw_compute_device_error(
        :validate_direct_imaging_target,
        :wrong_device,
        target,
        "$label occupies $(actual)",
    )
end

@inline function _require_exact_direct_imaging_array_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    actual = compute_device(storage)
    actual == target || _throw_wrong_direct_imaging_target(
        target, label, actual)
    return storage
end

@inline function _require_exact_direct_imaging_metadata_target(
    metadata::OpticalPlaneMetadata,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    actual = metadata.device
    actual == target || _throw_wrong_direct_imaging_target(
        target, "$label metadata", actual)
    return metadata
end

function _require_exact_direct_imaging_pupil_target(
    pupil::PupilFunction,
    target::AbstractComputeDevice,
)
    _require_exact_direct_imaging_metadata_target(
        pupil.metadata, target, "direct-imaging pupil")
    _require_exact_direct_imaging_array_target(
        pupil.support, target, "direct-imaging pupil support")
    _require_exact_direct_imaging_array_target(
        pupil.amplitude, target, "direct-imaging pupil amplitude")
    _require_exact_direct_imaging_array_target(
        pupil.opd, target, "direct-imaging pupil OPD")
    size(pupil.support) == pupil.metadata.dimensions || throw(
        DimensionMismatchError(
            "direct-imaging pupil support dimensions do not match its metadata"))
    typeof(backend(pupil.support)) === typeof(pupil.metadata.backend) || throw(
        InvalidConfiguration(
            "direct-imaging pupil support backend does not match its metadata"))
    validate_plane_storage(pupil.metadata, pupil.amplitude;
        label="direct-imaging pupil amplitude")
    validate_plane_storage(pupil.metadata, pupil.opd;
        label="direct-imaging pupil OPD")
    return pupil
end

function _require_exact_direct_imaging_input_target(
    input_plan::PupilDirectImagingInputPlan,
    target::AbstractComputeDevice,
)
    _require_exact_direct_imaging_metadata_target(
        input_plan.formation.input_metadata,
        target,
        "direct-imaging planned pupil",
    )
    _require_exact_direct_imaging_metadata_target(
        input_plan.formation.output_metadata,
        target,
        "direct-imaging planned field",
    )
    return input_plan
end

@inline _require_exact_direct_imaging_input_target(
    input_plan::PreformedFieldDirectImagingStrategy,
    ::AbstractComputeDevice,
) = input_plan

function _require_exact_direct_imaging_target(
    plan::DirectImagingPlan,
    target::AbstractComputeDevice,
)
    _require_exact_direct_imaging_input_target(plan.input_plan, target)
    _require_exact_direct_imaging_metadata_target(
        plan.field_metadata, target, "direct-imaging field")
    _require_exact_direct_imaging_metadata_target(
        plan.output_metadata, target, "direct-imaging output")
    propagation = plan.propagation
    propagation.input_metadata == plan.field_metadata || throw(
        InvalidConfiguration(
            "direct-imaging propagation does not match its prepared field"))
    require_same_plane_grid(
        propagation.output_metadata,
        plan.output_metadata;
        label="direct-imaging propagation and output",
        require_numeric_type=false,
    )
    require_compatible_radiometry(
        propagation.output_metadata,
        plan.output_metadata;
        label="direct-imaging propagation and output",
    )
    _require_exact_direct_imaging_metadata_target(
        propagation.output_metadata, target,
        "direct-imaging propagation output")
    return plan
end

function _require_exact_direct_imaging_target(
    prepared::PreparedDirectImaging,
    target::AbstractComputeDevice,
)
    _require_prepared_direct_bindings(prepared, prepared.input,
        prepared.field, prepared.output, prepared.workspace)
    _require_prepared_direct_input(
        prepared.input, prepared.field, prepared.plan.input_plan)
    _require_prepared_direct_field(prepared.field, prepared.plan)
    _require_prepared_direct_output(prepared.output, prepared.plan)
    _require_prepared_direct_workspace(prepared.workspace, prepared.plan)
    _require_exact_direct_imaging_target(prepared.plan, target)
    _require_exact_direct_imaging_owner_input_target(prepared.input, target)
    _require_exact_direct_imaging_owner_input_target(prepared.field, target)
    _require_exact_direct_imaging_metadata_target(
        prepared.output.metadata, target, "direct-imaging output")
    _require_exact_direct_imaging_array_target(
        prepared.output.values, target, "direct-imaging output storage")
    _require_exact_direct_imaging_array_target(
        prepared.workspace.propagation.scratch,
        target,
        "direct-imaging propagation scratch",
    )
    _require_exact_direct_imaging_array_target(
        prepared.workspace.unshifted_intensity,
        target,
        "direct-imaging unshifted intensity scratch",
    )
    return prepared
end

@inline _require_exact_direct_imaging_owner_input_target(
    pupil::PupilFunction,
    target::AbstractComputeDevice,
) = _require_exact_direct_imaging_pupil_target(pupil, target)

function _require_exact_direct_imaging_owner_input_target(
    field::ElectricField,
    target::AbstractComputeDevice,
)
    _require_exact_direct_imaging_metadata_target(
        field.metadata, target, "direct-imaging input field")
    _require_exact_direct_imaging_array_target(
        field.values, target, "direct-imaging input field storage")
    validate_plane_storage(field.metadata, field.values;
        label="direct-imaging input field")
    return field
end

function _require_exact_direct_imaging_target(
    prepared::PreparedIncoherentDirectImaging,
    target::AbstractComputeDevice,
)
    _require_direct_composition_bindings(prepared)
    _require_direct_component_products(prepared.components, prepared.products)
    prepared.output.metadata === prepared.accumulation.output_metadata || throw(
        InvalidConfiguration(
            "direct-imaging incoherent output does not match its accumulation plan"))
    length(prepared.products) == length(prepared.accumulation.inputs.metadata) ||
        throw(DimensionMismatchError(
            "direct-imaging incoherent product count does not match its accumulation plan"))
    _require_prepared_intensity_inputs(
        prepared.output.values,
        prepared.products,
        prepared.accumulation.inputs.metadata,
        prepared.accumulation.inputs.values,
    )
    @inbounds for component in prepared.components
        _require_exact_direct_imaging_target(component, target)
    end
    _require_exact_direct_imaging_metadata_target(
        prepared.output.metadata, target, "direct-imaging incoherent output")
    _require_exact_direct_imaging_array_target(
        prepared.output.values, target,
        "direct-imaging incoherent output storage")
    validate_plane_storage(prepared.output.metadata, prepared.output.values;
        label="direct-imaging incoherent output")
    return prepared
end

function _require_exact_direct_imaging_target(
    prepared::PreparedBundledDirectImaging,
    target::AbstractComputeDevice,
)
    _require_direct_composition_bindings(prepared)
    _require_direct_component_products(
        prepared.components, prepared.output.products)
    @inbounds for component in prepared.components
        _require_exact_direct_imaging_target(component, target)
    end
    return prepared
end

function _require_exact_direct_imaging_target(
    prepared,
    ::AbstractComputeDevice,
)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for prepared direct-imaging owner $(typeof(prepared))"))
end

function _prepare_direct_field!(field::ElectricField,
    pupil::PupilFunction, input_plan::PupilDirectImagingInputPlan)
    _require_prepared_direct_input(pupil, field, input_plan)
    fill_electric_field!(field, pupil, input_plan.formation)
    return field
end

function _prepare_direct_field!(field::ElectricField,
    input_field::ElectricField,
    input_plan::PreformedFieldDirectImagingStrategy,
)
    _require_prepared_direct_input(input_field, field, input_plan)
    return field
end

function _form_direct_image!(output::IntensityMap, field::ElectricField,
    plan::DirectImagingPlan, workspace::DirectImagingWorkspace)
    shifts = plan.shift_samples
    if iszero(shifts[1]) && iszero(shifts[2])
        fraunhofer_intensity_from_field!(output.values, field,
            plan.propagation, workspace.propagation)
    else
        fraunhofer_intensity_from_field!(workspace.unshifted_intensity,
            field, plan.propagation, workspace.propagation)
        shift_direct_image!(output.values, workspace.unshifted_intensity,
            shifts)
    end
    return output
end

function _form_prepared_direct_image!(
    prepared::PreparedDirectImaging,
    input::Union{PupilFunction,ElectricField},
    field::ElectricField,
    output::IntensityMap,
    workspace::DirectImagingWorkspace,
)
    _require_prepared_direct_execution(
        prepared, input, field, output, workspace)
    return _execute_prepared_direct_image!(prepared)
end

function _require_prepared_direct_execution(
    prepared::PreparedDirectImaging,
    input::Union{PupilFunction,ElectricField},
    field::ElectricField,
    output::IntensityMap,
    workspace::DirectImagingWorkspace,
)
    _require_prepared_direct_bindings(
        prepared, input, field, output, workspace)
    plan = prepared.plan
    _require_prepared_direct_input(input, field, plan.input_plan)
    _require_prepared_direct_output(output, plan)
    _require_prepared_direct_field(field, plan)
    _require_prepared_direct_workspace(workspace, plan)
    return prepared
end

function _execute_prepared_direct_image!(prepared::PreparedDirectImaging)
    plan = prepared.plan
    _prepare_direct_field!(prepared.field, prepared.input, plan.input_plan)
    return _form_direct_image!(
        prepared.output, prepared.field, plan, prepared.workspace)
end

function prepare_direct_imaging(pupil::PupilFunction,
    src::Union{Source,LGSSource}; zero_padding::Int=1)
    field = ElectricField(pupil, src; zero_padding=zero_padding)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    return _prepare_direct_imaging(pupil, src, field, output, propagation)
end

function prepare_direct_imaging(src::Union{Source,LGSSource},
    field::ElectricField)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    return _prepare_direct_imaging(src, field, output, propagation)
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
    C = typeof(first_component)
    components = Vector{C}(undef, length(sources))
    @inbounds components[1] = first_component
    slot = 2
    for source in Iterators.drop(sources, 1)
        component = prepare_direct_imaging(pupil, source;
            zero_padding=zero_padding)
        typeof(component) === C || throw(InvalidConfiguration(
            "direct-imaging composition must prepare one concrete leaf type; got $(C) and $(typeof(component))"))
        @inbounds components[slot] = component
        slot += 1
    end
    return FixedSizeVectorDefault{C}(components)
end

function _direct_imaging_products(
    components::AbstractVector{<:PreparedDirectImaging})
    isempty(components) && throw(InvalidConfiguration(
        "direct-imaging composition must contain at least one component"))
    first_product = direct_imaging_output(first(components))
    P = typeof(first_product)
    products = Vector{P}(undef, length(components))
    @inbounds products[1] = first_product
    slot = 2
    for component in Iterators.drop(components, 1)
        product = direct_imaging_output(component)
        typeof(product) === P || throw(InvalidConfiguration(
            "direct-imaging composition must produce one concrete product type; got $(P) and $(typeof(product))"))
        @inbounds products[slot] = product
        slot += 1
    end
    return FixedSizeVectorDefault{P}(products)
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
    bindings = PreparedDirectImagingCompositionBindings(
        _fixed_concrete_vector(
            components, "direct-imaging composition components"),
        _fixed_concrete_vector(
            products, "direct-imaging composition products"),
    )
    return PreparedIncoherentDirectImaging(
        _PREPARED_DIRECT_IMAGING_COMPOSITION_TOKEN,
        components,
        products,
        output,
        accumulation,
        bindings,
    )
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
    bindings = PreparedDirectImagingCompositionBindings(
        _fixed_concrete_vector(
            components, "direct-imaging composition components"),
        _fixed_concrete_vector(
            products, "direct-imaging composition products"),
    )
    return PreparedBundledDirectImaging(
        _PREPARED_DIRECT_IMAGING_COMPOSITION_TOKEN,
        components,
        output,
        bindings,
    )
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
        _require_prepared_direct_execution(
            component,
            component.input,
            component.field,
            component.output,
            component.workspace,
        )
    end
    @inbounds for component in components
        _execute_prepared_direct_image!(component)
    end
    return nothing
end

function _require_direct_composition_bindings(
    prepared::PreparedIncoherentDirectImaging,
)
    _require_exact_fixed_membership(
        prepared.components,
        prepared.bindings.components,
        "direct-imaging composition component",
    )
    _require_exact_fixed_membership(
        prepared.products,
        prepared.bindings.products,
        "direct-imaging composition product",
    )
    return prepared
end

function _require_direct_composition_bindings(
    prepared::PreparedBundledDirectImaging,
)
    _require_exact_fixed_membership(
        prepared.components,
        prepared.bindings.components,
        "direct-imaging composition component",
    )
    _require_exact_fixed_membership(
        prepared.output.products,
        prepared.bindings.products,
        "direct-imaging composition product",
    )
    return prepared
end

function form_direct_image!(prepared::PreparedDirectImaging)
    return _form_prepared_direct_image!(
        prepared,
        prepared.input,
        prepared.field,
        prepared.output,
        prepared.workspace,
    )
end

function form_direct_image!(prepared::PreparedIncoherentDirectImaging)
    _require_direct_composition_bindings(prepared)
    _form_direct_components!(prepared.components, prepared.products)
    accumulate_intensity!(prepared.output, prepared.products,
        prepared.accumulation)
    return prepared.output
end

function form_direct_image!(prepared::PreparedBundledDirectImaging)
    _require_direct_composition_bindings(prepared)
    _form_direct_components!(prepared.components, prepared.output.products)
    return prepared.output
end
