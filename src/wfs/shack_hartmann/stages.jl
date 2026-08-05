#
# Prepared Shack-Hartmann stage composition
#

"""Run-immutable physical and numerical contract for one SH rate plane."""
struct ShackHartmannOpticsPlan{M,S,T<:AbstractFloat,SS} <:
        AbstractWFSOpticsPlan
    microlens_array::M
    source::S
    threshold_convolution::T
    wavelength_m::T
    sampling_signature::SS
end

"""Exact live owner for one prepared Shack-Hartmann optics execution."""
struct PreparedShackHartmannOptics{P,F,W,I,O,R,L,B,D}
    plan::P
    optics::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    layout_binding::L
    backend::B
    device::D
end

"""Fixed-membership spectral SH optics plan."""
struct ShackHartmannOpticsBundlePlan{P<:FixedSizeVector} <:
        AbstractWFSOpticsPlan
    plans::P
end

"""Exact live owner for fixed-membership spectral SH optics."""
struct PreparedShackHartmannOpticsBundle{
    P,C<:FixedSizeVector,I,O}
    plan::P
    components::C
    input::I
    output::O
end

@inline wfs_optical_products(prepared::PreparedShackHartmannOptics) =
    prepared.output
@inline wfs_optical_products(prepared::PreparedShackHartmannOpticsBundle) =
    prepared.output

struct ShackHartmannSpectralComponent{S,T<:AbstractFloat} <: AbstractSource
    source::S
    wavelength_m::T
    photon_rate_m2_s::T
end

@inline wavelength(source::ShackHartmannSpectralComponent) =
    source.wavelength_m
@inline photon_irradiance(source::ShackHartmannSpectralComponent) =
    source.photon_rate_m2_s

"""Run-immutable Shack-Hartmann estimation and calibration contract."""
struct ShackHartmannEstimationPlan{
    E,P<:AbstractWFSMeasurementPath,C} <: AbstractWFSEstimationPlan
    extraction::E
    path::P
    calibration_binding::C
    n_lenslets::Int
    n_pixels_per_lenslet::Int
end

"""Exact live owner for one prepared Shack-Hartmann estimator."""
struct PreparedShackHartmannEstimator{P,W,WS,PR,I,M,WB,PB,B,D}
    plan::P
    sensor::W
    workspace::WS
    products::PR
    input::I
    measurement::M
    workspace_binding::WB
    products_binding::PB
    backend::B
    device::D
end

struct ShackHartmannCalibrationBinding{
    T<:AbstractFloat,R<:AbstractMatrix{T},U}
    layout_revision::UInt
    revision::UInt
    wavelength_m::T
    signature::UInt
    centroid_response::T
    output_units::U
    reference_signal::R
end

struct ShackHartmannLayoutBinding
    revision::UInt
end

@inline wfs_measurement_path(prepared::PreparedShackHartmannEstimator) =
    prepared.plan.path

@inline function _sh_microlens_workspace_binding(workspace)
    return (
        workspace.field,
        workspace.phasor,
        workspace.fft_buffer,
        workspace.fft_stack,
        workspace.intensity,
        workspace.intensity_stack,
        workspace.intensity_tmp_stack,
        workspace.temp,
        workspace.bin_buffer,
        workspace.spot,
        workspace.sampled_spot_cube,
        workspace.spot_cube_accum,
        workspace.fft_plan,
        workspace.fft_stack_plan,
        workspace.ifft_plan,
        workspace.ifft_stack_plan,
        workspace.elongation_kernel,
        workspace.lgs_kernel_fft,
        workspace.fft_asterism_stack,
        workspace.fft_asterism_plan,
        workspace.amp_scales,
        workspace.amp_scales_host,
        workspace.opd_to_cycles,
        workspace.opd_to_cycles_host,
    )
end

@inline function _sh_layout_storage_binding(layout::SubapertureLayout)
    return (layout.valid_mask, layout.valid_mask_host,
        layout.valid_indices_host)
end

@inline function _sh_estimator_workspace_binding(workspace)
    return (workspace.spot_cube, workspace.detector_noise_cube,
        workspace.spot_stats, workspace.spot_stats_accum,
        workspace.slopes_host, workspace.centroid_host)
end

@inline _sh_estimator_workspace_binding(::Nothing) = ()

@inline function _sh_estimator_owner_binding(sensor::ShackHartmannWFS)
    return (_sh_estimator_workspace_binding(sensor.workspace),
        _sh_layout_storage_binding(sensor.front_end.layout))
end

@inline function _sh_products_binding(products)
    return (products.slopes, products.legacy_spot_cube)
end

@inline _sh_input_storages(input::PupilFunction) =
    (input.amplitude, input.opd)
@inline _sh_input_storages(input::ElectricField) = (input.values,)
@inline _sh_input_storages(input::WFSObservation) = (input.storage,)

@inline _sh_mightalias_any(::AbstractArray, ::Tuple{}) = false
@inline function _sh_mightalias_any(value::AbstractArray, values::Tuple)
    return _wfs_storage_mightalias(value, first(values)) ||
        _sh_mightalias_any(value, Base.tail(values))
end

@inline _sh_any_alias(::Tuple{}, ::Tuple) = false
@inline function _sh_any_alias(values::Tuple, resources::Tuple)
    return _sh_mightalias_any(first(values), resources) ||
        _sh_any_alias(Base.tail(values), resources)
end

@inline function _require_sh_optics_aliases(input, output::IntensityMap,
    workspace, layout::SubapertureLayout)
    input_storages = _sh_input_storages(input)
    resources = (_sh_microlens_workspace_binding(workspace)...,
        layout.valid_mask, layout.valid_mask_host)
    (_sh_mightalias_any(output.values, input_storages) ||
        _sh_mightalias_any(output.values, resources) ||
        _sh_any_alias(input_storages, resources)) &&
        throw(WFSPreparationError(:wfs_optics, :aliasing,
            "Shack-Hartmann input, rate product, layout, and workspace storage must not alias"))
    return nothing
end

@inline function _require_sh_estimation_aliases(sensor::ShackHartmannWFS,
    input, measurement::WFSMeasurement)
    input_storages = _sh_input_storages(input)
    resources = (_sh_estimator_workspace_binding(sensor.workspace)...,
        sensor.products.slopes,
        sensor.products.legacy_spot_cube,
        sensor.front_end.layout.valid_mask,
        sensor.front_end.layout.valid_mask_host,
        sensor.calibration.reference_signal_2d,
        sensor.calibration.reference_signal_host)
    (_sh_mightalias_any(measurement.storage, input_storages) ||
        _sh_mightalias_any(measurement.storage, resources) ||
        _sh_any_alias(input_storages, resources)) &&
        throw(WFSPreparationError(:estimation, :aliasing,
            "Shack-Hartmann estimator input, measurement, calibration, workspace, and products must not alias"))
    return nothing
end

@inline function _sh_pupil_diameter(metadata::OpticalPlaneMetadata)
    return metadata.sampling[1] * metadata.dimensions[1]
end

@inline function _sh_sampling_signature(
    optics::ShackHartmannOptics)
    workspace = microlens_propagation_workspace(optics.propagation)
    return (size(workspace.fft_stack), workspace.effective_padding,
        workspace.binning_pixel_scale, workspace.sampled_n_pix_subap,
        workspace.phasor_ratio,
        subaperture_layout_revision(optics.front_end.layout))
end

@inline function _sh_output_wavelength(source::AbstractSource)
    return wavelength(source)
end

@inline function _sh_output_wavelength(source::SpectralSource)
    return require_sh_common_spectral_grid(source)
end

function require_sh_common_spectral_grid(source::SpectralSource)
    samples = spectral_bundle(source).samples
    isempty(samples) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "spectral source must contain at least one sample"))
    wavelength_ref = first(samples).wavelength
    @inbounds for index in 2:length(samples)
        samples[index].wavelength == wavelength_ref ||
            throw(WFSPreparationError(:wfs_optics, :plane_count,
                "distinct Shack-Hartmann spectral grids require an OpticalProductBundle"))
    end
    return wavelength_ref
end

@inline function _sh_front_end_wavelength(optics::ShackHartmannOptics,
    input::PupilFunction)
    source = optics.front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry,
        "dimensionless PupilFunction input requires an illumination source"))
    return _sh_output_wavelength(source)
end

@inline function _sh_front_end_wavelength(::ShackHartmannOptics,
    input::ElectricField)
    return _sh_monochromatic_wavelength(input.metadata.spectral)
end


@inline _sh_monochromatic_wavelength(channel::MonochromaticChannel) =
    channel.wavelength_m

function _sh_monochromatic_wavelength(::AbstractSpectralCoordinate)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Shack-Hartmann electric-field input must declare a monochromatic channel"))
end

function _require_sh_front_end_input(optics::ShackHartmannOptics,
    input::PupilFunction)
    optics.front_end.source === nothing && throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "Shack-Hartmann WFS optics require a source for PupilFunction input"))
    _require_sh_pupil_semantics(input, :wfs_optics)
    return nothing
end

function _require_sh_pupil_semantics(input::PupilFunction, stage::Symbol)
    _require_sh_metric_coordinates(input.metadata.coordinate_domain, stage)
    _require_sh_default_orientation(input.metadata.orientation, stage)
    _require_sh_centered_geometry(input.metadata, "pupil input", stage)
    _require_sh_pupil_normalization(input.metadata.normalization, stage)
    _require_sh_pupil_measure(input.metadata.spatial_measure, stage)
    _require_sh_coherent_input(input.metadata.coherence, stage)
    _require_sh_achromatic_pupil(input.metadata.spectral, stage)
    return nothing
end

function _require_sh_front_end_input(optics::ShackHartmannOptics,
    input::ElectricField)
    optics.front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    _require_sh_metric_coordinates(input.metadata.coordinate_domain,
        :wfs_optics)
    _require_sh_default_orientation(input.metadata.orientation,
        :wfs_optics)
    _require_sh_centered_geometry(input.metadata, "electric-field input",
        :wfs_optics)
    _require_sh_field_normalization(input.metadata.normalization)
    _require_sh_field_measure(input.metadata.spatial_measure)
    _require_sh_coherent_input(input.metadata.coherence, :wfs_optics)
    return nothing
end


@inline _require_sh_metric_coordinates(::MetricCoordinates, ::Symbol) =
    nothing

function _require_sh_metric_coordinates(::AbstractPlaneCoordinateDomain,
    stage::Symbol)
    throw(WFSPreparationError(stage, :plane_metadata,
        "Shack-Hartmann pupil input must use metric coordinates"))
end

function _require_sh_default_orientation(orientation::PlaneAxisOrientation,
    stage::Symbol)
    orientation.axes == (:x, :y) && orientation.signs == (1, 1) ||
        throw(WFSPreparationError(stage, :plane_metadata,
            "Shack-Hartmann currently requires (:x, :y) axes with positive signs"))
    return nothing
end

function _require_sh_centered_geometry(metadata::OpticalPlaneMetadata,
    label::AbstractString, stage::Symbol)
    expected_centering = (
        axis_centering(metadata.dimensions[1]),
        axis_centering(metadata.dimensions[2]),
    )
    metadata.centering == expected_centering ||
        throw(WFSPreparationError(stage, :plane_metadata,
            "Shack-Hartmann $label centering does not match its dimensions"))
    metadata.origin == centered_grid_origin(metadata.dimensions,
        metadata.sampling) ||
        throw(WFSPreparationError(stage, :plane_metadata,
            "Shack-Hartmann $label must use a centered-grid origin"))
    return nothing
end

@inline _require_sh_pupil_normalization(::DimensionlessNormalization,
    ::Symbol) = nothing

function _require_sh_pupil_normalization(::AbstractOpticalNormalization,
    stage::Symbol)
    throw(WFSPreparationError(stage, :radiometry,
        "Shack-Hartmann PupilFunction amplitude must be dimensionless"))
end

@inline _require_sh_pupil_measure(::PointSampledMeasure, ::Symbol) = nothing

function _require_sh_pupil_measure(::AbstractSpatialMeasure, stage::Symbol)
    throw(WFSPreparationError(stage, :radiometry,
        "Shack-Hartmann PupilFunction amplitude must be point sampled"))
end

@inline _require_sh_coherent_input(::CoherentFieldCombination, ::Symbol) =
    nothing

function _require_sh_coherent_input(::AbstractCombinationPolicy,
    stage::Symbol)
    throw(WFSPreparationError(stage, :radiometry,
        "Shack-Hartmann pupil input must declare coherent field combination"))
end

@inline _require_sh_achromatic_pupil(::AchromaticSpectralCoordinate,
    ::Symbol) = nothing

function _require_sh_achromatic_pupil(::AbstractSpectralCoordinate,
    stage::Symbol)
    throw(WFSPreparationError(stage, :plane_metadata,
        "Shack-Hartmann PupilFunction input must be achromatic"))
end


@inline _require_sh_field_normalization(::PhotonRateNormalization) = nothing

function _require_sh_field_normalization(::AbstractOpticalNormalization)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "electric-field input must carry photon-rate normalization"))
end

@inline _require_sh_field_measure(::CellIntegratedMeasure) = nothing

function _require_sh_field_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "electric-field input must carry cell-integrated photon rate"))
end

function _require_sh_stage_domains(optics::ShackHartmannOptics,
    input, output::IntensityMap)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:wfs_optics, :backend,
            "Shack-Hartmann input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:wfs_optics, :device,
            "Shack-Hartmann input and output occupy different devices"))
    workspace = microlens_propagation_workspace(optics.propagation)
    propagation_storage = workspace.fft_stack
    typeof(input.metadata.backend) === typeof(backend(propagation_storage)) ||
        throw(WFSPreparationError(:wfs_optics, :backend,
            "Shack-Hartmann input and microlens propagation backends differ"))
    input.metadata.device == compute_device(propagation_storage) ||
        throw(WFSPreparationError(:wfs_optics, :device,
            "Shack-Hartmann input and microlens propagation occupy different devices"))
    typeof(input.metadata.backend) ===
        typeof(backend(optics.front_end.layout.valid_mask)) ||
        throw(WFSPreparationError(:wfs_optics, :backend,
            "Shack-Hartmann input and subaperture-layout backends differ"))
    input.metadata.device == compute_device(optics.front_end.layout.valid_mask) ||
        throw(WFSPreparationError(:wfs_optics, :device,
            "Shack-Hartmann input and subaperture layout occupy different devices"))
    propagation_type = eltype(workspace.intensity)
    _sh_input_numeric_type(input) === propagation_type ||
        throw(WFSPreparationError(:wfs_optics, :numeric_type,
            "Shack-Hartmann input type must match prepared propagation precision"))
    output.metadata.numeric_type === propagation_type ||
        throw(WFSPreparationError(:wfs_optics, :numeric_type,
            "Shack-Hartmann rate output type must match prepared propagation precision"))
    return nothing
end

function _require_sh_layout_geometry(layout::SubapertureLayout,
    n_sub::Int, input, stage::Symbol)
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions[2] == resolution ||
        throw(WFSPreparationError(stage, :shape,
            "Shack-Hartmann pupil input must be square"))
    input.metadata.sampling[1] == input.metadata.sampling[2] ||
        throw(WFSPreparationError(stage, :plane_metadata,
            "Shack-Hartmann pupil sampling must be square"))
    resolution % n_sub == 0 ||
        throw(WFSPreparationError(stage, :shape,
            "Shack-Hartmann pupil resolution must be divisible by the lenslet count"))
    layout.subap_pixels == div(resolution, n_sub) ||
        throw(WFSPreparationError(stage, :shape,
            "SubapertureLayout sampling does not match the pupil grid"))
    layout_diameter = layout.pitch_m * n_sub
    pupil_diameter = _sh_pupil_diameter(input.metadata)
    T_geometry = promote_type(typeof(layout_diameter),
        typeof(pupil_diameter))
    isapprox(T_geometry(layout_diameter), T_geometry(pupil_diameter);
        rtol=T_geometry(8) * eps(T_geometry), atol=zero(T_geometry)) ||
        throw(WFSPreparationError(stage, :plane_metadata,
            "SubapertureLayout diameter does not match the pupil grid"))
    return resolution, pupil_diameter
end

function _require_sh_rate_mosaic(optics::ShackHartmannOptics,
    output::IntensityMap, wavelength_m::Real, pupil_diameter_m::Real)
    workspace = microlens_propagation_workspace(optics.propagation)
    n_sub = n_lenslets(optics)
    n_pix = workspace.sampled_n_pix_subap
    size(output.values) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Shack-Hartmann rate output must be an n_lenslets-by-n_lenslets tiled spot mosaic"))
    _require_sh_rate_coordinates(output.metadata.coordinate_domain)
    _require_sh_rate_measure(output.metadata.spatial_measure)
    _require_sh_default_orientation(output.metadata.orientation,
        :wfs_optics)
    _require_sh_centered_geometry(output.metadata, "rate output",
        :wfs_optics)
    T = eltype(workspace.intensity)
    expected_scale_arcsec = sh_pixel_scale_init(
        pupil_diameter_m / n_sub,
        workspace.effective_padding, wavelength_m) *
        workspace.binning_pixel_scale
    expected_scale_rad = T(expected_scale_arcsec / ARCSEC_PER_RAD)
    output.metadata.sampling == (expected_scale_rad, expected_scale_rad) ||
        throw(WFSPreparationError(:wfs_optics, :plane_metadata,
            "Shack-Hartmann rate output sampling does not match prepared microlens propagation"))
    _require_sh_output_wavelength(output.metadata.spectral, wavelength_m)
    return nothing
end


@inline _require_sh_rate_coordinates(::AngularCoordinates) = nothing

function _require_sh_rate_coordinates(::AbstractPlaneCoordinateDomain)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Shack-Hartmann rate output must use angular coordinates"))
end

@inline _require_sh_rate_measure(::CellIntegratedMeasure) = nothing

function _require_sh_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "Shack-Hartmann rate output must contain cell-integrated photon rate"))
end


function _require_sh_output_wavelength(channel::MonochromaticChannel,
    wavelength_m::Real)
    channel.wavelength_m == typeof(channel.wavelength_m)(wavelength_m) ||
        throw(WFSPreparationError(
        :wfs_optics, :plane_metadata,
        "Shack-Hartmann rate output wavelength does not match its optical input"))
    return nothing
end


function _require_sh_output_wavelength(::AbstractSpectralCoordinate,
    ::Real)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Shack-Hartmann rate output must declare a monochromatic channel"))
end

function _prepare_sh_optics_candidate(
    optics::ShackHartmannOptics, input,
    output::IntensityMap)
    validate_wfs_optical_input(input)
    validate_wfs_optical_products(output)
    _require_sh_front_end_input(optics, input)
    _require_sh_stage_domains(optics, input, output)
    resolution, pupil_diameter = _require_sh_layout_geometry(
        optics.front_end.layout, n_lenslets(optics), input,
        :wfs_optics)
    wavelength_m = _sh_front_end_wavelength(optics, input)
    _prepare_microlens_sampling_wavelength!(optics, resolution,
        _sh_pupil_diameter(input.metadata), wavelength_m)
    _prepare_sh_source_workspace!(execution_style(_sh_input_storage(input)),
        optics, input, optics.front_end.source, wavelength_m)
    _require_sh_rate_mosaic(optics, output, wavelength_m,
        pupil_diameter)
    sampling_signature = _sh_sampling_signature(optics)
    T = eltype(microlens_propagation_workspace(
        optics.propagation).intensity)
    plan = ShackHartmannOpticsPlan(
        optics.front_end.microlens_array, optics.front_end.source,
        T(optics.front_end.threshold_convolution), T(wavelength_m),
        sampling_signature)
    workspace = microlens_propagation_workspace(optics.propagation)
    _require_sh_optics_aliases(input, output, workspace,
        optics.front_end.layout)
    return PreparedShackHartmannOptics(plan, optics,
        workspace, input, output,
        _sh_microlens_workspace_binding(workspace),
        _sh_layout_storage_binding(optics.front_end.layout),
        input.metadata.backend, input.metadata.device)
end

function _sh_independent_optics(
    optics::ShackHartmannOptics{F,PR}, source) where {F,PR}
    sourced = _sh_optics_with_source(optics, source)
    propagation = _prepare_microlens_propagation_like(sourced.propagation)
    return ShackHartmannOptics(sourced.front_end,
        propagation)
end

function _sh_commit_optics_candidate!(
    optics::ShackHartmannOptics,
    candidate::PreparedShackHartmannOptics{P}) where {P}
    workspace = microlens_propagation_workspace(optics.propagation)
    _replace_microlens_propagation_workspace!(workspace,
        candidate.workspace)
    return PreparedShackHartmannOptics(candidate.plan, optics,
        workspace, candidate.input, candidate.output,
        _sh_microlens_workspace_binding(workspace),
        candidate.layout_binding, candidate.backend, candidate.device)
end

function prepare_wfs_optics(
    optics::ShackHartmannOptics, input,
    output::IntensityMap)
    candidate_optics = _sh_independent_optics(optics,
        optics.front_end.source)
    candidate = _prepare_sh_optics_candidate(
        candidate_optics, input, output)
    return _sh_commit_optics_candidate!(optics, candidate)
end

@inline _prepare_sh_source_workspace!(::ExecutionStyle,
    ::ShackHartmannOptics, input, source, wavelength_m) =
    nothing

@inline _prepare_sh_asterism_workspace!(::ExecutionStyle,
    ::ShackHartmannOptics, ::Asterism) = nothing

function _prepare_sh_asterism_workspace!(::AcceleratorStyle,
    optics::ShackHartmannOptics, source::Asterism)
    sh_stacked_asterism_compatible(source) &&
        ensure_sh_asterism_buffers!(optics, length(source.sources))
    return nothing
end

function _prepare_sh_source_workspace!(style::ExecutionStyle,
    optics::ShackHartmannOptics, input::PupilFunction,
    source::Asterism, wavelength_m::Real)
    _prepare_sh_asterism_workspace!(style, optics, source)
    @inbounds for component in source.sources
        _prepare_sh_source_workspace!(style, optics, input, component,
            wavelength(component))
    end
    return nothing
end

function _prepare_sh_source_workspace!(::ExecutionStyle,
    optics::ShackHartmannOptics, input::PupilFunction,
    source::LGSSource, wavelength_m::Real)
    return _prepare_sh_lgs_workspace!(lgs_profile(source), optics, input,
        source, wavelength_m)
end

function _prepare_sh_source_workspace!(::ExecutionStyle,
    optics::ShackHartmannOptics, input::PupilFunction,
    source::ShackHartmannSpectralComponent{<:LGSSource}, wavelength_m::Real)
    return _prepare_sh_lgs_workspace!(lgs_profile(source.source), optics,
        input, source.source, wavelength_m)
end

function _prepare_sh_lgs_workspace!(::LGSProfileNone,
    optics::ShackHartmannOptics, ::PupilFunction,
    source::LGSSource, ::Real)
    workspace = microlens_propagation_workspace(optics.propagation)
    T = eltype(workspace.elongation_kernel)
    factor = T(lgs_elongation_factor(source))
    workspace.elongation_kernel =
        prepare_elongation_kernel!(workspace.elongation_kernel, factor)
    return nothing
end

function _prepare_sh_lgs_workspace!(::LGSProfileNaProfile,
    optics::ShackHartmannOptics, input::PupilFunction,
    source::LGSSource, wavelength_m::Real)
    metadata = input.metadata
    ensure_lgs_kernels!(optics, source, metadata.dimensions,
        _sh_pupil_diameter(metadata), metadata.sampling, metadata.origin,
        wavelength_m)
    return nothing
end

function prepare_wfs_optics(
    optics::ShackHartmannOptics{
        <:ShackHartmannOpticalFrontEnd{
            M,L,FT,<:SpectralSource}}, input,
    output::OpticalProductBundle) where {M,L,FT}
    validate_wfs_optical_input(input)
    validate_wfs_optical_products(output)
    source = optics.front_end.source
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "spectral Shack-Hartmann output count must match its spectral samples"))
    input_type = _sh_input_numeric_type(input)
    T = promote_type(input_type, typeof(photon_irradiance(source)))
    first_sample = first(samples)
    first_component = ShackHartmannSpectralComponent(source.source,
        T(first_sample.wavelength),
        T(photon_irradiance(source)) * T(first_sample.weight))
    first_optics = _sh_independent_optics(optics, first_component)
    first_plan = _prepare_sh_optics_candidate(first_optics,
        input, output[1])
    components = Vector{typeof(first_plan)}(undef, length(samples))
    components[1] = first_plan
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_optics = _sh_independent_optics(optics, component)
        components[index] = _prepare_sh_optics_candidate(
            component_optics, input, output[index])
    end
    fixed_components = FixedSizeVectorDefault{typeof(first_plan)}(components)
    plan_values = FixedSizeVectorDefault{
        typeof(first_plan.plan)}(map(component -> component.plan,
            components))
    plan = ShackHartmannOpticsBundlePlan(plan_values)
    return PreparedShackHartmannOpticsBundle(plan,
        fixed_components, input, output)
end

@inline function _sh_optics_with_source(
    optics::ShackHartmannOptics, source)
    return ShackHartmannOptics(optics, source)
end

@inline function _require_sh_optics_binding(output::IntensityMap, input,
    prepared::PreparedShackHartmannOptics)
    output.metadata === prepared.output.metadata &&
        output.values === prepared.output.values ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann rate output does not match prepared storage"))
    input === prepared.input || throw(WFSPreparationError(:wfs_optics,
        :prepared_binding,
        "Shack-Hartmann pupil input does not match prepared storage"))
    optics = prepared.optics
    plan = prepared.plan
    workspace = microlens_propagation_workspace(optics.propagation)
    workspace === prepared.workspace &&
        _sh_microlens_workspace_binding(workspace) ===
            prepared.workspace_binding ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann microlens workspace was replaced after preparation"))
    _sh_layout_storage_binding(optics.front_end.layout) ===
        prepared.layout_binding ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann layout storage was replaced after preparation"))
    optics.front_end.microlens_array === plan.microlens_array &&
        optics.front_end.source === plan.source &&
        isequal(optics.front_end.threshold_convolution,
            plan.threshold_convolution) ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann optical definition changed after preparation"))
    _sh_sampling_signature(optics) == plan.sampling_signature ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann microlens sampling no longer matches its prepared plan"))
    input.metadata.backend === prepared.backend &&
        output.metadata.backend === prepared.backend &&
        input.metadata.device == prepared.device &&
        output.metadata.device == prepared.device ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Shack-Hartmann backend or device binding changed after preparation"))
    _require_sh_optics_aliases(input, output, workspace,
        optics.front_end.layout)
    return nothing
end

@inline validate_wfs_optics_binding(output::IntensityMap, input,
    plan::PreparedShackHartmannOptics) =
    _require_sh_optics_binding(output, input, plan)

@kernel function sh_explicit_pupil_stack_kernel!(fft_stack, valid_mask,
    amplitude, opd, phasor, amp_scale, opd_to_cycles, n_sub::Int, sub::Int,
    ox::Int, oy::Int, n::Int, pad::Int)
    x, y, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && i <= n_sub && j <= n_sub
        index = sh_lenslet_index(i, j, n_sub)
        xi = x - ox
        yi = y - oy
        value = zero(eltype(fft_stack))
        if 1 <= xi <= sub && 1 <= yi <= sub && @inbounds(valid_mask[i, j])
            px = (i - 1) * sub + xi
            py = (j - 1) * sub + yi
            if px <= n && py <= n
                @inbounds value = amp_scale * amplitude[px, py] *
                    cispi(opd_to_cycles * opd[px, py]) * phasor[x, y]
            end
        end
        @inbounds fft_stack[x, y, index] = value
    end
end

@kernel function sh_explicit_field_stack_kernel!(fft_stack, valid_mask,
    pupil_field, phasor, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int,
    pad::Int)
    x, y, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && i <= n_sub && j <= n_sub
        index = sh_lenslet_index(i, j, n_sub)
        xi = x - ox
        yi = y - oy
        value = zero(eltype(fft_stack))
        if 1 <= xi <= sub && 1 <= yi <= sub && @inbounds(valid_mask[i, j])
            px = (i - 1) * sub + xi
            py = (j - 1) * sub + yi
            if px <= n && py <= n
                @inbounds value = pupil_field[px, py] * phasor[x, y]
            end
        end
        @inbounds fft_stack[x, y, index] = value
    end
end

@kernel function sh_explicit_pupil_asterism_stack_kernel!(fft_stack,
    valid_mask, amplitude, opd, phasor, amp_scales, opd_to_cycles,
    n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int,
    n_spots::Int, n_src::Int)
    x, y, source_index, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && source_index <= n_src &&
            i <= n_sub && j <= n_sub
        spot_index = sh_lenslet_index(i, j, n_sub)
        stack_index = (source_index - 1) * n_spots + spot_index
        xi = x - ox
        yi = y - oy
        value = zero(eltype(fft_stack))
        if 1 <= xi <= sub && 1 <= yi <= sub &&
                @inbounds(valid_mask[i, j])
            px = (i - 1) * sub + xi
            py = (j - 1) * sub + yi
            if px <= n && py <= n
                @inbounds value = amp_scales[source_index] *
                    amplitude[px, py] *
                    cispi(opd_to_cycles[source_index] * opd[px, py]) *
                    phasor[x, y]
            end
        end
        @inbounds fft_stack[x, y, stack_index] = value
    end
end

function _form_sh_explicit_stack!(::ScalarCPUStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = microlens_propagation_workspace(sensor.propagation)
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(workspace.intensity_stack)
    pupil_cell_area = T(input.metadata.sampling[1] *
        input.metadata.sampling[2])
    amp_scale = sqrt(T(photon_irradiance(source)) * pupil_cell_area)
    opd_to_cycles = T(2) / T(wavelength_m)
    fill!(workspace.fft_stack, zero(eltype(workspace.fft_stack)))
    index = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if sensor.front_end.layout.valid_mask[i, j]
            for y in 1:sub, x in 1:sub
                px = (i - 1) * sub + x
                py = (j - 1) * sub + y
                workspace.fft_stack[ox + x, oy + y, index] =
                    amp_scale * input.amplitude[px, py] *
                    cispi(opd_to_cycles * input.opd[px, py]) *
                    workspace.phasor[ox + x, oy + y]
            end
        end
        index += 1
    end
    return _finish_sh_explicit_stack!(ScalarCPUStyle(), sensor, input,
        source)
end

function _form_sh_explicit_stack!(style::AcceleratorStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = microlens_propagation_workspace(sensor.propagation)
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(workspace.intensity_stack)
    pupil_cell_area = T(input.metadata.sampling[1] *
        input.metadata.sampling[2])
    amp_scale = sqrt(T(photon_irradiance(source)) * pupil_cell_area)
    opd_to_cycles = T(2) / T(wavelength_m)
    launch_kernel_async!(style, sh_explicit_pupil_stack_kernel!,
        workspace.fft_stack, sensor.front_end.layout.valid_mask, input.amplitude,
        input.opd, workspace.phasor, amp_scale, opd_to_cycles, n_sub, sub,
        ox, oy, n, pad; ndrange=(pad, pad, n_sub, n_sub))
    synchronize_backend!(style)
    return _finish_sh_explicit_stack!(style, sensor, input, source)
end

function _form_sh_explicit_stack!(::ScalarCPUStyle,
    sensor::ShackHartmannOptics, input::ElectricField, ::Nothing,
    wavelength_m)
    workspace = microlens_propagation_workspace(sensor.propagation)
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    fill!(workspace.fft_stack, zero(eltype(workspace.fft_stack)))
    index = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if sensor.front_end.layout.valid_mask[i, j]
            for y in 1:sub, x in 1:sub
                px = (i - 1) * sub + x
                py = (j - 1) * sub + y
                workspace.fft_stack[ox + x, oy + y, index] =
                    input.values[px, py] * workspace.phasor[ox + x, oy + y]
            end
        end
        index += 1
    end
    return _finish_sh_explicit_stack!(ScalarCPUStyle(), sensor, input,
        nothing)
end

function _form_sh_explicit_stack!(style::AcceleratorStyle,
    sensor::ShackHartmannOptics, input::ElectricField, ::Nothing,
    wavelength_m)
    workspace = microlens_propagation_workspace(sensor.propagation)
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    launch_kernel_async!(style, sh_explicit_field_stack_kernel!,
        workspace.fft_stack, sensor.front_end.layout.valid_mask, input.values,
        workspace.phasor, n_sub, sub, ox, oy, n, pad;
        ndrange=(pad, pad, n_sub, n_sub))
    synchronize_backend!(style)
    return _finish_sh_explicit_stack!(style, sensor, input, nothing)
end

function _form_sh_explicit_asterism_serial!(style::ExecutionStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::Asterism)
    wavelength(source)
    workspace = microlens_propagation_workspace(sensor.propagation)
    fill!(workspace.spot_cube_accum,
        zero(eltype(workspace.spot_cube_accum)))
    @inbounds for component in source.sources
        _form_sh_explicit_stack!(style, sensor, input, component,
            wavelength(component))
        @. workspace.spot_cube_accum = workspace.spot_cube_accum +
            workspace.sampled_spot_cube
    end
    copyto!(workspace.sampled_spot_cube, workspace.spot_cube_accum)
    return workspace.sampled_spot_cube
end

@inline _form_sh_explicit_asterism!(style::ScalarCPUStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::Asterism) =
    _form_sh_explicit_asterism_serial!(style, sensor, input, source)

function _form_sh_explicit_asterism!(style::AcceleratorStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::Asterism)
    sh_stacked_asterism_compatible(source) ||
        return _form_sh_explicit_asterism_serial!(style, sensor, input,
            source)
    workspace = microlens_propagation_workspace(sensor.propagation)
    n_src = length(source.sources)
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    n_spots = n_sub * n_sub
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(workspace.intensity_stack)
    pupil_cell_area = T(input.metadata.sampling[1] *
        input.metadata.sampling[2])
    @inbounds for index in eachindex(source.sources)
        component = source.sources[index]
        workspace.amp_scales_host[index] =
            sqrt(T(photon_irradiance(component)) * pupil_cell_area)
        workspace.opd_to_cycles_host[index] = T(2) / T(wavelength(component))
    end
    copyto!(workspace.amp_scales, workspace.amp_scales_host)
    copyto!(workspace.opd_to_cycles, workspace.opd_to_cycles_host)
    total = n_spots * n_src
    fft_view = @view workspace.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view workspace.intensity_tmp_stack[:, :, 1:total]
    launch_kernel_async!(style, sh_explicit_pupil_asterism_stack_kernel!,
        fft_view, sensor.front_end.layout.valid_mask, input.amplitude, input.opd,
        workspace.phasor, workspace.amp_scales, workspace.opd_to_cycles,
        n_sub, sub, ox, oy, n, pad, n_spots, n_src;
        ndrange=(pad, pad, n_src, n_sub, n_sub))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, workspace.fft_asterism_plan)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, complex_abs2_stack_kernel!, intensity_view,
        fft_view, intensity_scale, pad, total;
        ndrange=size(intensity_view))
    queue_kernel!(phase, reduce_grouped_blocks_kernel!,
        workspace.intensity_stack, intensity_view, n_spots, n_src,
        size(workspace.intensity_stack, 1),
        size(workspace.intensity_stack, 2);
        ndrange=size(workspace.intensity_stack))
    finish_kernel_phase!(phase)
    sample_spot_stack!(style, sensor)
    return workspace.sampled_spot_cube
end

@inline _form_sh_explicit_stack!(style::ScalarCPUStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

@inline _form_sh_explicit_stack!(style::AcceleratorStyle,
    sensor::ShackHartmannOptics, input::PupilFunction,
    source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

function _finish_sh_explicit_stack!(style::ExecutionStyle,
    sensor::ShackHartmannOptics, input, source)
    workspace = microlens_propagation_workspace(sensor.propagation)
    execute_fft_plan!(workspace.fft_stack, workspace.fft_stack_plan)
    T = eltype(workspace.intensity_stack)
    pad = size(workspace.fft_stack, 1)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    _sh_stack_intensity!(style, sensor, intensity_scale, pad)
    _apply_sh_source_spot_model!(style, sensor, input, source)
    sample_spot_stack!(style, sensor)
    return workspace.sampled_spot_cube
end

@inline _apply_sh_source_spot_model!(::ExecutionStyle,
    ::ShackHartmannOptics, input, source) = nothing

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannOptics, input, source::LGSSource)
    return _apply_sh_lgs_spot_model!(lgs_profile(source), style, sensor,
        input, source, wavelength(source))
end

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannOptics, input,
    source::ShackHartmannSpectralComponent{<:LGSSource})
    return _apply_sh_lgs_spot_model!(lgs_profile(source.source), style,
        sensor, input, source.source, wavelength(source))
end

function _apply_sh_lgs_spot_model!(::LGSProfileNone,
    ::ExecutionStyle, sensor::ShackHartmannOptics, input,
    source::LGSSource,
    wavelength_m::Real)
    workspace = microlens_propagation_workspace(sensor.propagation)
    workspace.elongation_kernel = apply_prepared_elongation_stack!(
        workspace.intensity_stack, lgs_elongation_factor(source),
        workspace.intensity_tmp_stack, workspace.elongation_kernel)
    return nothing
end

function _apply_sh_lgs_spot_model!(::LGSProfileNaProfile,
    ::ExecutionStyle, sensor::ShackHartmannOptics, input,
    source::LGSSource,
    wavelength_m::Real)
    metadata = input.metadata
    ensure_lgs_kernels!(sensor, source, metadata.dimensions,
        _sh_pupil_diameter(metadata), metadata.sampling, metadata.origin,
        wavelength_m)
    workspace = microlens_propagation_workspace(sensor.propagation)
    apply_lgs_convolution_stack!(workspace.intensity_stack,
        workspace.lgs_kernel_fft, workspace.fft_stack,
        workspace.fft_stack_plan, workspace.ifft_stack_plan)
    return nothing
end

function _sh_stack_intensity!(::ScalarCPUStyle,
    sensor::ShackHartmannOptics,
    intensity_scale, ::Int)
    workspace = microlens_propagation_workspace(sensor.propagation)
    @inbounds for index in axes(workspace.fft_stack, 3),
            y in axes(workspace.fft_stack, 2),
            x in axes(workspace.fft_stack, 1)
        workspace.intensity_stack[x, y, index] =
            abs2(workspace.fft_stack[x, y, index]) * intensity_scale
    end
    return workspace.intensity_stack
end

function _sh_stack_intensity!(style::AcceleratorStyle,
    sensor::ShackHartmannOptics, intensity_scale, pad::Int)
    workspace = microlens_propagation_workspace(sensor.propagation)
    launch_kernel!(style, complex_abs2_stack_kernel!,
        workspace.intensity_stack, workspace.fft_stack,
        intensity_scale, pad, n_lenslets(sensor)^2;
        ndrange=size(workspace.intensity_stack))
    return workspace.intensity_stack
end

function form_wfs_optical_products!(output::IntensityMap, input,
    prepared::PreparedShackHartmannOptics)
    _require_sh_optics_binding(output, input, prepared)
    optics = prepared.optics
    wavelength_m = _sh_front_end_wavelength(optics, input)
    _form_sh_explicit_stack!(execution_style(_sh_input_storage(input)),
        optics, input, optics.front_end.source,
        wavelength_m)
    _tile_shack_hartmann_spot_cube!(output.values,
        prepared.workspace.sampled_spot_cube, n_lenslets(optics))
    return output
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    prepared::PreparedShackHartmannOpticsBundle)
    validate_wfs_optics_binding(output, input, prepared)
    @inbounds for index in eachindex(prepared.components)
        form_wfs_optical_products!(output[index], input,
            prepared.components[index])
    end
    return output
end

function validate_wfs_optics_binding(
    output::OpticalProductBundle, input,
    prepared::PreparedShackHartmannOpticsBundle)
    output === prepared.output && input === prepared.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "spectral Shack-Hartmann products do not match prepared storage"))
    length(prepared.plan.plans) == length(prepared.components) ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "spectral Shack-Hartmann plan membership changed"))
    @inbounds for index in eachindex(prepared.components)
        prepared.components[index].plan === prepared.plan.plans[index] ||
            throw(WFSPreparationError(:wfs_optics,
                :prepared_binding,
                "spectral Shack-Hartmann component plan was replaced"))
        validate_wfs_optics_binding(output[index], input,
            prepared.components[index])
    end
    return nothing
end

function _require_sh_observation_semantics(observation::WFSObservation)
    isequal(observation.metadata.layout, :lenslet_mosaic) ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "Shack-Hartmann detector observations require :lenslet_mosaic layout"))
    return nothing
end

function _require_sh_real_observation(observation::WFSObservation,
    stage::Symbol)
    observation.metadata.numeric_type <: Real ||
        throw(WFSPreparationError(stage, :numeric_type,
            "Shack-Hartmann observations require real detector samples"))
    return nothing
end

function _require_sh_floating_measurement(measurement::WFSMeasurement)
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "Shack-Hartmann measurements require floating-point storage"))
    return nothing
end

function _require_sh_measurement_semantics(
    sensor::ShackHartmannWFS{<:Diffractive},
    measurement::WFSMeasurement)
    isequal(measurement.units, sensor.calibration.output_units) ||
        throw(WFSPreparationError(:estimation, :units,
            "diffractive Shack-Hartmann measurement units must match the calibration output units"))
    isequal(measurement.metadata.kind, :centroid_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "diffractive Shack-Hartmann measurements require :centroid_slopes kind"))
    return nothing
end

function _require_sh_measurement_semantics(
    ::ShackHartmannWFS{<:Geometric}, measurement::WFSMeasurement)
    isequal(measurement.units, :radian) ||
        throw(WFSPreparationError(:estimation, :units,
            "geometric Shack-Hartmann measurements require :radian units"))
    isequal(measurement.metadata.kind, :geometric_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "geometric Shack-Hartmann measurements require :geometric_slopes kind"))
    return nothing
end

function _require_sh_storage_domain(stage::Symbol, metadata, storage,
    label::AbstractString)
    return _require_wfs_storage_domain(stage, metadata, storage, label)
end

function _prepare_sh_calibration_binding(sensor::ShackHartmannWFS)
    calibration = sensor.calibration
    calibration.calibrated || throw(WFSPreparationError(:estimation,
        :estimator,
        "diffractive Shack-Hartmann estimation requires explicit calibration"))
    isfinite(calibration.centroid_response) &&
        calibration.centroid_response != zero(calibration.centroid_response) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "Shack-Hartmann slope calibration must be finite and nonzero"))
    isfinite(calibration.wavelength) &&
        calibration.wavelength > zero(calibration.wavelength) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "Shack-Hartmann calibration wavelength must be finite and positive"))
    return ShackHartmannCalibrationBinding(
        subaperture_layout_revision(sensor.front_end.layout),
        calibration.revision,
        calibration.wavelength, calibration.signature,
        calibration.centroid_response, calibration.output_units,
        calibration.reference_signal_2d)
end

function _require_sh_calibration_binding(sensor::ShackHartmannWFS,
    binding::ShackHartmannCalibrationBinding)
    calibration = sensor.calibration
    subaperture_layout_revision(sensor.front_end.layout) ==
        binding.layout_revision &&
        calibration.calibrated &&
        calibration.revision == binding.revision &&
        isequal(calibration.wavelength, binding.wavelength_m) &&
        calibration.signature == binding.signature &&
        isequal(calibration.centroid_response, binding.centroid_response) &&
        isequal(calibration.output_units, binding.output_units) &&
        calibration.reference_signal_2d === binding.reference_signal ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann layout or calibration changed after estimator preparation"))
    return nothing
end


@inline function _require_sh_layout_binding(sensor::ShackHartmannWFS,
    binding::ShackHartmannLayoutBinding)
    subaperture_layout_revision(sensor.front_end.layout) ==
        binding.revision ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann subaperture layout changed after estimator preparation"))
    return nothing
end

function prepare_wfs_estimation(sensor::ShackHartmannWFS{<:Diffractive},
    observation::WFSObservation, measurement::WFSMeasurement)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    _require_sh_observation_semantics(observation)
    _require_sh_real_observation(observation, :estimation)
    _require_sh_floating_measurement(measurement)
    _require_sh_measurement_semantics(sensor, measurement)
    _require_sh_storage_domain(:estimation, observation.metadata,
        sensor.workspace.spot_cube, "observation")
    _require_sh_storage_domain(:estimation, measurement.metadata,
        sensor.products.slopes, "measurement")
    _require_sh_storage_domain(:estimation, observation.metadata,
        sensor.front_end.layout.valid_mask, "observation/layout")
    n_sub = n_lenslets(sensor)
    n_pix = sensor.optics.propagation.workspace.sampled_n_pix_subap
    size(observation.storage) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann estimator requires a tiled lenslet mosaic"))
    size(measurement.storage) == size(sensor.products.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann measurement storage has the wrong slope shape"))
    calibration_binding = _prepare_sh_calibration_binding(sensor)
    ensure_sh_acquisition_buffers!(sensor, n_pix)
    _require_sh_estimation_aliases(sensor, observation, measurement)
    plan = ShackHartmannEstimationPlan(slope_extraction_model(sensor),
        AcquiredObservationPath(), calibration_binding, n_sub, n_pix)
    return PreparedShackHartmannEstimator(plan, sensor, sensor.workspace,
        sensor.products, observation, measurement,
        _sh_estimator_owner_binding(sensor),
        _sh_products_binding(sensor.products), measurement.metadata.backend,
        measurement.metadata.device)
end

@kernel function sh_unpack_mosaic_kernel!(spot_cube, mosaic, n_sub::Int,
    n_pix::Int)
    i, j, x, y = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub && x <= n_pix && y <= n_pix
        index = sh_lenslet_index(i, j, n_sub)
        @inbounds spot_cube[index, x, y] =
            mosaic[(i - 1) * n_pix + x, (j - 1) * n_pix + y]
    end
end

function _unpack_sh_mosaic!(::ScalarCPUStyle, spot_cube, mosaic,
    n_sub::Int, n_pix::Int)
    @inbounds for y in 1:n_pix, x in 1:n_pix, j in 1:n_sub, i in 1:n_sub
        index = sh_lenslet_index(i, j, n_sub)
        spot_cube[index, x, y] =
            mosaic[(i - 1) * n_pix + x, (j - 1) * n_pix + y]
    end
    return spot_cube
end

function _unpack_sh_mosaic!(style::AcceleratorStyle, spot_cube, mosaic,
    n_sub::Int, n_pix::Int)
    launch_kernel!(style, sh_unpack_mosaic_kernel!, spot_cube, mosaic,
        n_sub, n_pix; ndrange=(n_sub, n_sub, n_pix, n_pix))
    return spot_cube
end

function _require_sh_estimation_binding(measurement::WFSMeasurement,
    input, prepared::PreparedShackHartmannEstimator)
    measurement === prepared.measurement && input === prepared.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator storage does not match its prepared owner"))
    sensor = prepared.sensor
    sensor.workspace === prepared.workspace &&
        sensor.products === prepared.products &&
        _sh_estimator_owner_binding(sensor) ===
            prepared.workspace_binding &&
        _sh_products_binding(prepared.products) ===
            prepared.products_binding ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator workspace or products were replaced after preparation"))
    measurement.metadata.backend === prepared.backend &&
        input.metadata.backend === prepared.backend &&
        measurement.metadata.device == prepared.device &&
        input.metadata.device == prepared.device ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator backend or device binding changed after preparation"))
    _require_sh_estimation_aliases(sensor, input, measurement)
    return sensor
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation,
    prepared::PreparedShackHartmannEstimator{
        <:ShackHartmannEstimationPlan{
            E,<:AcquiredObservationPath}}) where {E}
    sensor = _require_sh_estimation_binding(measurement, observation,
        prepared)
    plan = prepared.plan
    _require_sh_calibration_binding(sensor, plan.calibration_binding)
    n_sub = plan.n_lenslets
    n_pix = plan.n_pixels_per_lenslet
    style = execution_style(observation.storage)
    _unpack_sh_mosaic!(style, prepared.workspace.spot_cube,
        observation.storage, n_sub, n_pix)
    peak = sh_safe_peak_value(prepared.workspace.spot_cube)
    sh_signal_from_spots_calibrated!(sensor, peak,
        plan.extraction)
    copyto!(measurement.storage, prepared.products.slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    prepared::PreparedShackHartmannEstimator)
    sensor = _require_sh_estimation_binding(measurement, input, prepared)
    binding = prepared.plan.calibration_binding
    _require_sh_estimation_state_binding(sensor, binding)
    return nothing
end


@inline _require_sh_estimation_state_binding(sensor::ShackHartmannWFS,
    binding::ShackHartmannCalibrationBinding) =
    _require_sh_calibration_binding(sensor, binding)

@inline _require_sh_estimation_state_binding(sensor::ShackHartmannWFS,
    binding::ShackHartmannLayoutBinding) =
    _require_sh_layout_binding(sensor, binding)

function prepare_wfs_estimation(sensor::ShackHartmannWFS{<:Geometric},
    input::PupilFunction, measurement::WFSMeasurement)
    validate_wfs_optical_input(input)
    validate_wfs_measurement(measurement)
    _require_sh_pupil_semantics(input, :estimation)
    _require_sh_layout_geometry(sensor.front_end.layout,
        n_lenslets(sensor), input,
        :estimation)
    _require_sh_measurement_semantics(sensor, measurement)
    _require_sh_floating_measurement(measurement)
    _require_sh_storage_domain(:estimation, input.metadata,
        sensor.products.slopes, "geometric input")
    _require_sh_storage_domain(:estimation, input.metadata,
        sensor.front_end.layout.valid_mask, "geometric input/layout")
    _require_sh_storage_domain(:estimation, measurement.metadata,
        sensor.products.slopes, "geometric measurement")
    size(measurement.storage) == size(sensor.products.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric Shack-Hartmann measurement storage has the wrong slope shape"))
    binding = ShackHartmannLayoutBinding(
        subaperture_layout_revision(sensor.front_end.layout))
    _require_sh_estimation_aliases(sensor, input, measurement)
    plan = ShackHartmannEstimationPlan(slope_extraction_model(sensor),
        DirectMeasurementPath(), binding, n_lenslets(sensor), 0)
    return PreparedShackHartmannEstimator(plan, sensor, sensor.workspace,
        sensor.products, input, measurement,
        _sh_estimator_owner_binding(sensor),
        _sh_products_binding(sensor.products), measurement.metadata.backend,
        measurement.metadata.device)
end


function prepare_wfs_estimation(::ShackHartmannWFS{<:Geometric},
    ::ElectricField, ::WFSMeasurement)
    throw(WFSPreparationError(:estimation, :unsupported,
        "geometric Shack-Hartmann estimation requires OPD-bearing PupilFunction input"))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    prepared::PreparedShackHartmannEstimator{
        <:ShackHartmannEstimationPlan{
            E,<:DirectMeasurementPath}}) where {E}
    sensor = _require_sh_estimation_binding(measurement, input, prepared)
    _require_sh_layout_binding(sensor,
        prepared.plan.calibration_binding)
    geometric_wavefront_slopes!(prepared.products.slopes, input.opd,
        sensor.front_end.layout.valid_mask, input.metadata.sampling)
    copyto!(measurement.storage, prepared.products.slopes)
    return measurement
end

"""
    shack_hartmann_rate_map(sensor, input, source=nothing)

Allocate a detector-plane photon-rate mosaic with metadata compatible with a
prepared Shack-Hartmann optical front end. Execution remains allocation-free;
applications may instead construct and retain an equivalent `IntensityMap`.
"""
function shack_hartmann_rate_map(
    sensor::ShackHartmannWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source::SpectralSource)
    return shack_hartmann_rate_map(
        ShackHartmannOptics(sensor.optics, source), input,
        source)
end

function shack_hartmann_rate_map(
    model::ShackHartmannOptics,
    input::Union{PupilFunction,ElectricField}, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(model.propagation.workspace.intensity)
    first_sample = first(samples)
    first_component = ShackHartmannSpectralComponent(source.source,
        T(first_sample.wavelength),
        T(photon_irradiance(source)) * T(first_sample.weight))
    first_map = shack_hartmann_rate_map(
        _sh_optics_with_source(model, first_component), input)
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        maps[index] = shack_hartmann_rate_map(
            _sh_optics_with_source(model, component), input)
    end
    return OpticalProductBundle(maps)
end


function shack_hartmann_rate_map(sensor::ShackHartmannWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    optics = source === nothing ? sensor.optics :
        ShackHartmannOptics(sensor.optics, source)
    return shack_hartmann_rate_map(optics, input)
end

function shack_hartmann_rate_map(model::ShackHartmannOptics,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    resolved_model = source === nothing ? model :
        _sh_optics_with_source(model, source)
    wavelength_m = _sh_front_end_wavelength(resolved_model, input)
    sampling = _sh_microlens_sampling_configuration(resolved_model,
        input.metadata.dimensions[1],
        _sh_pupil_diameter(input.metadata), wavelength_m)
    propagation = microlens_propagation_workspace(resolved_model.propagation)
    n = n_lenslets(resolved_model) * sampling.spot_samples_per_axis
    T = eltype(propagation.intensity)
    values = similar(_sh_input_storage(input), T, n, n)
    fill!(values, zero(T))
    pixel_scale_rad = T(sampling.pixel_scale_arcsec / ARCSEC_PER_RAD)
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(pixel_scale_rad, pixel_scale_rad),
        spectral=MonochromaticChannel(T(wavelength_m)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end


@inline _sh_input_storage(input::PupilFunction) = input.opd
@inline _sh_input_storage(input::ElectricField) = input.values
@inline _sh_input_numeric_type(input::PupilFunction) = eltype(input.opd)
@inline _sh_input_numeric_type(input::ElectricField) =
    typeof(real(zero(eltype(input.values))))
