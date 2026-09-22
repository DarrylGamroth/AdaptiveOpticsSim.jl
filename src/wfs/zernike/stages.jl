#
# Prepared Zernike WFS stages
#

"""
Run-immutable physical and numerical contract for one Zernike detector-plane
photon-arrival-rate map.
"""
struct ZernikeOpticsPlan{P,S} <: AbstractWFSOpticsPlan
    propagation::P
    source::S
end

"""Exact live owner for one prepared Zernike optics execution."""
struct PreparedZernikeOptics{P,F,W,I,O,R,B,D}
    plan::P
    front_end::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    backend::B
    device::D
end

@inline wfs_optical_products(prepared::PreparedZernikeOptics) =
    prepared.output

@inline function _zernike_propagation_workspace_binding(workspace)
    return (workspace.field, workspace.focal_field, workspace.pupil_field,
        workspace.phasor, workspace.phase_mask, workspace.pupil_intensity,
        workspace.nominal_frame, workspace.fft_plan, workspace.ifft_plan)
end

@inline _zernike_input_storages(input::PupilFunction) =
    (input.amplitude, input.opd, input.support)
@inline _zernike_input_storages(input::ElectricField) = (input.values,)

@inline _zernike_mightalias_any(value, ::Tuple{}) = false
@inline function _zernike_mightalias_any(value, values::Tuple)
    return _wfs_storage_mightalias(value, first(values)) ||
        _zernike_mightalias_any(value, Base.tail(values))
end

@inline _zernike_any_alias(::Tuple{}) = false
@inline function _zernike_any_alias(values::Tuple)
    remaining = Base.tail(values)
    return _zernike_mightalias_any(first(values), remaining) ||
        _zernike_any_alias(remaining)
end

function _require_zernike_optics_aliases(input, output, workspace)
    storages = (_zernike_input_storages(input)..., output.values,
        _zernike_propagation_workspace_binding(workspace)...)
    _zernike_any_alias(storages) && throw(WFSPreparationError(
        :wfs_optics, :aliasing,
        "Zernike input, rate product, and propagation workspace must not alias"))
    return nothing
end

@inline modulated_wfs_propagation_storage(
    front_end::ZernikeOpticalFrontEnd) =
    zernike_propagation_workspace(front_end).field

function ZernikeOpticalFrontEnd(sensor::ZernikeWFS, source=nothing)
    front_end = sensor.front_end
    return ZernikeOpticalFrontEnd(front_end.phase_spot,
        front_end.propagation, front_end.binning, source)
end

@inline function zernike_rate_dimensions(front_end::ZernikeOpticalFrontEnd)
    pupil_samples = zernike_propagation_plan(
        front_end.propagation).pupil_samples
    return (div(pupil_samples, front_end.binning),
        div(pupil_samples, front_end.binning))
end

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "Zernike WFS optics require a source for PupilFunction input"))
    require_leaf_source(source, "prepared Zernike optics")
    return source
end

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a Zernike source"))
    return nothing
end

@inline _zernike_front_end_wavelength(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction) = modulated_input_wavelength(input,
        front_end.source)
@inline _zernike_front_end_wavelength(::ZernikeOpticalFrontEnd,
    input::ElectricField) = modulated_input_wavelength(input)

@inline _require_zernike_rate_coordinates(
    ::NormalizedPupilCoordinates) = nothing

function _require_zernike_rate_coordinates(::AbstractPlaneCoordinateDomain)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Zernike detector output must use normalized pupil coordinates"))
end

@inline _require_zernike_rate_measure(::CellIntegratedMeasure) = nothing

function _require_zernike_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "Zernike detector output must carry cell-integrated rate"))
end

function _require_zernike_rate_wavelength(channel::MonochromaticChannel,
    wavelength_m)
    channel.wavelength_m == wavelength_m || throw(
        WFSPreparationError(:wfs_optics, :plane_metadata,
            "Zernike detector output wavelength differs from its input"))
    return nothing
end

function _require_zernike_rate_wavelength(
    ::AbstractSpectralCoordinate, ::Any)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Zernike detector output wavelength differs from its input"))
end

function _require_zernike_rate_map(output::IntensityMap,
    expected_dimensions, wavelength_m)
    validate_wfs_optical_products(output)
    _require_zernike_rate_coordinates(output.metadata.coordinate_domain)
    _require_zernike_rate_measure(output.metadata.spatial_measure)
    size(output.values) == expected_dimensions || throw(
        WFSPreparationError(:wfs_optics, :shape,
            "Zernike detector output has the wrong prepared dimensions"))
    _require_zernike_rate_wavelength(output.metadata.spectral, wavelength_m)
    return output
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    pupil_resolution = zernike_propagation_plan(
        front_end.propagation).pupil_resolution
    input.metadata.dimensions == (pupil_resolution,
        pupil_resolution) || throw(WFSPreparationError(
        :wfs_optics, :shape,
        "Zernike pupil input dimensions differ from the prepared relay"))
    return nothing
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::ElectricField)
    workspace = zernike_propagation_workspace(front_end)
    input.metadata.dimensions == size(workspace.field) || throw(
        WFSPreparationError(:wfs_optics, :shape,
            "Zernike ElectricField dimensions differ from the prepared diffraction grid"))
    return nothing
end

function prepare_wfs_optics(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_zernike_front_end_source(front_end, input)
    _require_zernike_input_geometry(front_end, input)
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    _require_zernike_rate_map(output, zernike_rate_dimensions(front_end),
        wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    workspace = zernike_propagation_workspace(front_end)
    eltype(workspace.pupil_intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :wfs_optics, :numeric_type,
            "Zernike output precision differs from prepared propagation"))
    _require_zernike_optics_aliases(input, output, workspace)
    propagation = front_end.propagation
    propagation_plan = zernike_propagation_plan(propagation)
    plan = ZernikeOpticsPlan(propagation_plan, front_end.source)
    return PreparedZernikeOptics(plan, front_end, workspace, input, output,
        _zernike_propagation_workspace_binding(workspace),
        input.metadata.backend, input.metadata.device)
end

function zernike_rate_map(sensor::ZernikeWFS,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return zernike_rate_map(ZernikeOpticalFrontEnd(sensor, source), input)
end

function zernike_rate_map(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    dimensions = zernike_rate_dimensions(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    T = eltype(zernike_propagation_workspace(front_end).pupil_intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    normalized_sampling = T(front_end.binning / propagation_plan.pupil_samples)
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(normalized_sampling, normalized_sampling),
        spectral=MonochromaticChannel(T(wavelength_m)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

function _require_zernike_optical_binding(
    prepared::PreparedZernikeOptics, input, output)
    input === prepared.input && output === prepared.output || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike optical products do not match their prepared plan"))
    workspace = prepared.workspace
    prepared.backend === input.metadata.backend &&
        prepared.device == input.metadata.device || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike optical input target changed after preparation"))
    prepared.plan.propagation ===
        zernike_propagation_plan(prepared.front_end.propagation) &&
        prepared.front_end.phase_spot ===
            prepared.plan.propagation.phase_spot &&
        prepared.front_end.source === prepared.plan.source &&
        workspace === zernike_propagation_workspace(prepared.front_end) &&
        _zernike_propagation_workspace_binding(workspace) ===
            prepared.workspace_binding || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike propagation storage changed after preparation"))
    return nothing
end

@inline validate_wfs_optics_binding(output::IntensityMap, input,
    plan::PreparedZernikeOptics) =
    _require_zernike_optical_binding(plan, input, output)

function _form_zernike_input_field!(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    propagation = zernike_propagation_workspace(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    T = eltype(propagation.pupil_intensity)
    n = propagation_plan.pupil_resolution
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    cell_area = T(input.metadata.sampling[1] * input.metadata.sampling[2])
    amplitude_scale = sqrt(T(photon_irradiance(front_end.source)) *
        cell_area)
    opd_to_cycles = T(2) / T(wavelength(front_end.source))
    fill!(propagation.field, zero(eltype(propagation.field)))
    @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] =
        amplitude_scale * input.amplitude * cispi(opd_to_cycles * input.opd)
    return propagation.field
end

function _form_zernike_input_field!(front_end::ZernikeOpticalFrontEnd,
    input::ElectricField)
    workspace = zernike_propagation_workspace(front_end)
    copyto!(workspace.field, input.values)
    return workspace.field
end

function _form_zernike_rate!(output::AbstractMatrix,
    front_end::ZernikeOpticalFrontEnd, input)
    propagation = zernike_propagation_workspace(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    _form_zernike_input_field!(front_end, input)
    copyto!(propagation.focal_field, propagation.field)
    @. propagation.focal_field *= propagation.phasor
    execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
    @. propagation.focal_field *= propagation.phase_mask
    copyto!(propagation.pupil_field, propagation.focal_field)
    execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    n = propagation_plan.pupil_resolution
    pad = size(propagation.pupil_field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    @views @. propagation.pupil_intensity =
        abs2(propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    sampling = div(n, propagation_plan.pupil_samples)
    bin2d!(propagation.nominal_frame, propagation.pupil_intensity, sampling)
    if front_end.binning == 1
        copyto!(output, propagation.nominal_frame)
    else
        bin2d!(output, propagation.nominal_frame, front_end.binning)
    end
    return output
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedZernikeOptics)
    validate_wfs_optics_binding(output, input, plan)
    _form_zernike_rate!(output.values, plan.front_end, input)
    return output
end
