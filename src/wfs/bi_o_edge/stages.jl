#
# Prepared Bi-O-edge WFS stages
#

"""
Run-immutable physical and numerical contract for one Bi-O-edge detector-plane
photon-arrival-rate map.
"""
struct BiOEdgeOpticsPlan{P,O,S,L<:AbstractPreparedFourPupilLGS} <:
        AbstractWFSOpticsPlan
    propagation::P
    operating_modulation::O
    source::S
    lgs_model::L
    propagation_revision::UInt
end
"""Exact live owner for one prepared Bi-O-edge optics execution."""
struct PreparedBiOEdgeOptics{P,F,W,I,O,R,B,D}
    plan::P
    front_end::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    backend::B
    device::D
end

"""Run-immutable contract for one fixed Bi-O-edge optics-product bundle."""
struct BiOEdgeOpticsBundlePlan{P<:Tuple} <: AbstractWFSOpticsPlan
    plans::P
end

"""Exact live owner for one prepared Bi-O-edge optics-product bundle."""
struct PreparedBiOEdgeOpticsBundle{P,C<:Tuple,I,O}
    plan::P
    components::C
    input::I
    output::O
end

@inline wfs_optical_products(prepared::PreparedBiOEdgeOptics) =
    prepared.output
@inline wfs_optical_products(prepared::PreparedBiOEdgeOpticsBundle) =
    prepared.output

@inline function _bi_o_edge_propagation_workspace_binding(workspace)
    return (workspace.field, workspace.focal_field, workspace.pupil_field,
        workspace.bi_o_edge_masks, workspace.phasor, workspace.intensity,
        workspace.temp, workspace.scratch, workspace.asterism_stack,
        workspace.fft_buffer, workspace.fft_plan, workspace.ifft_plan,
        workspace.elongation_kernel, workspace.lgs_kernel_fft)
end

@inline modulated_wfs_propagation_storage(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end).field

@inline function bi_o_edge_output_sampling_factor(
    front_end::BiOEdgeOpticalFrontEnd, pupil_resolution::Int)
    pupil_resolution % front_end.pupil_samples == 0 ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge pupil resolution must be divisible by pupil_samples"))
    pupil_sample = div(pupil_resolution, front_end.pupil_samples)
    return front_end.binning == 1 ? pupil_sample :
        2 * pupil_sample * front_end.binning
end

function bi_o_edge_output_dimensions(front_end::BiOEdgeOpticalFrontEnd,
    pupil_resolution::Int)
    factor = bi_o_edge_output_sampling_factor(front_end, pupil_resolution)
    side = size(bi_o_edge_propagation_workspace(front_end).intensity, 1)
    side % factor == 0 || throw(WFSPreparationError(:wfs_optics,
        :shape, "Bi-O-edge sampling does not evenly divide the detector plane"))
    output_side = div(side, factor)
    return (output_side, output_side)
end

@inline function _bi_o_edge_front_end_wavelength(
    front_end::BiOEdgeOpticalFrontEnd, input::PupilFunction)
    return modulated_input_wavelength(input, front_end.source)
end

@inline function _bi_o_edge_front_end_wavelength(
    ::BiOEdgeOpticalFrontEnd, input::ElectricField)
    return modulated_input_wavelength(input)
end

function _require_bi_o_edge_source(front_end::BiOEdgeOpticalFrontEnd,
    ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "Bi-O-edge WFS optics require a source for PupilFunction input"))
    return _require_single_bi_o_edge_source(source)
end

@inline _require_single_bi_o_edge_source(source) = source

function _require_single_bi_o_edge_source(source::SpectralSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "spectral Bi-O-edge optics require an OpticalProductBundle"))
end

function _require_single_bi_o_edge_source(source::Asterism)
    throw(WFSPreparationError(:wfs_optics,
        :plane_count,
        "asterism Bi-O-edge optics require path-local pupil inputs"))
end

function _require_single_bi_o_edge_source(source::ExtendedSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "extended Bi-O-edge optics require path-local pupil inputs"))
end

function _require_bi_o_edge_source(front_end::BiOEdgeOpticalFrontEnd,
    ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    return nothing
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_bi_o_edge_source(front_end, input)
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge pupil input must be square"))
    size(front_end.modulation.phases, 1) == resolution ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge modulation was prepared for another pupil resolution"))
    expected = bi_o_edge_output_dimensions(front_end, resolution)
    wavelength_m = _bi_o_edge_front_end_wavelength(front_end, input)
    require_four_pupil_rate_map(output, expected, wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    eltype(bi_o_edge_propagation_workspace(front_end).intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :wfs_optics, :numeric_type,
            "Bi-O-edge output precision differs from prepared propagation"))
    lgs_model = prepare_four_pupil_lgs(front_end.source, input, front_end)
    propagation = front_end.propagation
    propagation_plan = bi_o_edge_propagation_plan(propagation)
    workspace = bi_o_edge_propagation_workspace(propagation)
    plan = BiOEdgeOpticsPlan(propagation_plan, front_end.modulation,
        front_end.source, lgs_model, workspace.revision)
    return PreparedBiOEdgeOptics(plan, front_end, workspace, input, output,
        _bi_o_edge_propagation_workspace_binding(workspace),
        input.metadata.backend, input.metadata.device)
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    output::OpticalProductBundle)
    return prepare_bi_o_edge_optical_bundle(front_end, input, output,
        front_end.source)
end

function prepare_bi_o_edge_optical_bundle(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::OpticalProductBundle,
    source::SpectralSource)
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge spectral output count does not match the source"))
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    plans = ntuple(length(samples)) do index
        sample = samples[index]
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        prepare_wfs_optics(
            bi_o_edge_front_end_with_source(front_end, component), input,
            output[index])
    end
    return PreparedBiOEdgeOpticsBundle(
        BiOEdgeOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, input, output)
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle)
    return prepare_bi_o_edge_optical_bundle(front_end, inputs, output,
        front_end.source)
end

function prepare_bi_o_edge_optical_bundle(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle,
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local pupil count does not match the source count"))
    length(output) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local output count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "Bi-O-edge path-local source collection is empty"))
    plans = ntuple(length(sources)) do index
        prepare_wfs_optics(
            bi_o_edge_front_end_with_source(front_end, sources[index]),
            inputs[index], output[index])
    end
    return PreparedBiOEdgeOpticsBundle(
        BiOEdgeOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, inputs, output)
end

function prepare_bi_o_edge_optical_bundle(front_end, input, output, source)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge product bundles require a spectral or path-expanded source"))
end

function _bi_o_edge_native_rate!(front_end::BiOEdgeOpticalFrontEnd,
    input::PupilFunction, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = bi_o_edge_propagation_workspace(front_end)
    source = front_end.source
    resolution = size(input.opd, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    amplitude_scale = sqrt(T(photon_irradiance(source)) *
        T(input.metadata.sampling[1] * input.metadata.sampling[2]))
    opd_to_cycles = T(2) / T(wavelength(source))
    fill!(propagation.intensity, zero(T))
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = amplitude_scale * weight *
            input.amplitude * front_end.modulation.phases[:, :, point] *
            cispi(opd_to_cycles * input.opd)
        copyto!(propagation.focal_field, propagation.field)
        accumulate_bi_o_edge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function _bi_o_edge_native_rate!(front_end::BiOEdgeOpticalFrontEnd,
    input::ElectricField, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = bi_o_edge_propagation_workspace(front_end)
    resolution = size(input.values, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    fill!(propagation.intensity, zero(T))
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = weight * input.values *
            front_end.modulation.phases[:, :, point]
        copyto!(propagation.focal_field, propagation.field)
        accumulate_bi_o_edge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedBiOEdgeOptics)
    validate_wfs_optics_binding(output, input, plan)
    native = _bi_o_edge_native_rate!(plan.front_end, input,
        plan.plan.lgs_model)
    factor = bi_o_edge_output_sampling_factor(plan.front_end,
        input.metadata.dimensions[1])
    bin2d!(output.values, native, factor)
    return output
end

function validate_wfs_optics_binding(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedBiOEdgeOptics)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge optical products do not match prepared storage"))
    workspace = bi_o_edge_propagation_workspace(plan.front_end)
    workspace === plan.workspace &&
        _bi_o_edge_propagation_workspace_binding(workspace) ===
            plan.workspace_binding &&
        workspace.revision == plan.plan.propagation_revision ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge propagation workspace changed after preparation"))
    plan.front_end.propagation.plan === plan.plan.propagation &&
        plan.front_end.amplitude_mask ===
            plan.plan.propagation.amplitude_mask &&
        plan.front_end.modulation === plan.plan.operating_modulation &&
        plan.front_end.source === plan.plan.source ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge optics definition changed after preparation"))
    return nothing
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedBiOEdgeOpticsBundle)
    validate_wfs_optics_binding(output, input, plan)
    return form_four_pupil_bundle!(output, input, plan.components)
end

function validate_wfs_optics_binding(
    output::OpticalProductBundle, input,
    plan::PreparedBiOEdgeOpticsBundle)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge spectral products do not match prepared storage"))
    validate_four_pupil_bundle_binding(output, input, plan.components)
    return nothing
end

function bi_o_edge_rate_map(sensor::BiOEdgeWFS,
    inputs::Union{Tuple,AbstractVector}, source)
    return bi_o_edge_rate_map(BiOEdgeOpticalFrontEnd(sensor, source), inputs)
end

function bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector})
    return bi_o_edge_path_rate_bundle(front_end, inputs, front_end.source)
end

function bi_o_edge_path_rate_bundle(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector},
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local pupil count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "Bi-O-edge path-local source collection is empty"))
    first_map = bi_o_edge_rate_map(
        bi_o_edge_front_end_with_source(front_end, sources[1]), inputs[1])
    maps = Vector{typeof(first_map)}(undef, length(sources))
    maps[1] = first_map
    @inbounds for index in 2:length(sources)
        maps[index] = bi_o_edge_rate_map(
            bi_o_edge_front_end_with_source(front_end, sources[index]),
            inputs[index])
    end
    return OpticalProductBundle(maps)
end

function bi_o_edge_path_rate_bundle(front_end, inputs, source)
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-local Bi-O-edge inputs require an Asterism or ExtendedSource"))
end

function bi_o_edge_rate_map(sensor::BiOEdgeWFS,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return bi_o_edge_rate_map(BiOEdgeOpticalFrontEnd(sensor, source), input)
end

function bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    return _bi_o_edge_rate_map(front_end, input, front_end.source)
end

@inline _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input,
    source::SpectralSource) =
    _bi_o_edge_spectral_rate_bundle(front_end, input, source)

function _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input,
    source::Union{Asterism,ExtendedSource})
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-expanded Bi-O-edge sources require path-local pupil inputs"))
end

function _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input, source)
    wavelength_m = _bi_o_edge_front_end_wavelength(front_end, input)
    dimensions = bi_o_edge_output_dimensions(front_end,
        input.metadata.dimensions[1])
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    factor = bi_o_edge_output_sampling_factor(front_end,
        input.metadata.dimensions[1])
    normalized_sampling = T(factor / input.metadata.dimensions[1])
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(normalized_sampling, normalized_sampling),
        spectral=MonochromaticChannel(T(wavelength_m)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

function _bi_o_edge_spectral_rate_bundle(front_end::BiOEdgeOpticalFrontEnd,
    input, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    function component_map(sample)
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = bi_o_edge_front_end_with_source(front_end,
            component)
        return bi_o_edge_rate_map(component_front_end, input)
    end
    first_map = component_map(first(samples))
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        maps[index] = component_map(samples[index])
    end
    return OpticalProductBundle(maps)
end
