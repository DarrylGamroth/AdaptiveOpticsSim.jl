#
# Prepared Shack-Hartmann stage composition
#

struct PreparedShackHartmannOpticalFormation{F,I,O,S}
    front_end::F
    input::I
    output::O
    sampling_signature::S
end

struct PreparedShackHartmannOpticalBundleFormation{P<:AbstractVector,I,O}
    plans::P
    input::I
    output::O
end

struct ShackHartmannSpectralComponent{S,T<:AbstractFloat} <: AbstractSource
    source::S
    wavelength_m::T
    photon_rate_m2_s::T
end

@inline wavelength(source::ShackHartmannSpectralComponent) =
    source.wavelength_m
@inline photon_irradiance(source::ShackHartmannSpectralComponent) =
    source.photon_rate_m2_s

struct PreparedShackHartmannEstimator{W,I,M,P<:AbstractWFSMeasurementPath,C}
    sensor::W
    input::I
    measurement::M
    path::P
    calibration_binding::C
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

@inline wfs_measurement_path(plan::PreparedShackHartmannEstimator) = plan.path

@inline function _sh_pupil_diameter(metadata::OpticalPlaneMetadata)
    return metadata.sampling[1] * metadata.dimensions[1]
end

@inline function _sh_sampling_signature(
    front_end::ShackHartmannOpticalFrontEnd)
    workspace = front_end.propagation
    return (size(workspace.fft_stack), workspace.effective_padding,
        workspace.binning_pixel_scale, workspace.sampled_n_pix_subap,
        workspace.phasor_ratio,
        subaperture_layout_revision(front_end.layout))
end

@inline function _sh_output_wavelength(source::AbstractSource)
    return wavelength(source)
end

@inline function _sh_output_wavelength(source::SpectralSource)
    return require_sh_common_spectral_grid(source)
end

function require_sh_common_spectral_grid(source::SpectralSource)
    samples = spectral_bundle(source).samples
    isempty(samples) && throw(WFSPreparationError(:optical_formation,
        :plane_count, "spectral source must contain at least one sample"))
    wavelength_ref = first(samples).wavelength
    @inbounds for index in 2:length(samples)
        samples[index].wavelength == wavelength_ref ||
            throw(WFSPreparationError(:optical_formation, :plane_count,
                "distinct Shack-Hartmann spectral grids require an OpticalProductBundle"))
    end
    return wavelength_ref
end

@inline function _sh_front_end_wavelength(front_end::ShackHartmannOpticalFrontEnd,
    input::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:optical_formation,
        :radiometry,
        "dimensionless PupilFunction input requires an illumination source"))
    return _sh_output_wavelength(source)
end

@inline function _sh_front_end_wavelength(::ShackHartmannOpticalFrontEnd,
    input::ElectricField)
    return _sh_monochromatic_wavelength(input.metadata.spectral)
end


@inline _sh_monochromatic_wavelength(channel::MonochromaticChannel) =
    channel.wavelength_m

function _sh_monochromatic_wavelength(::AbstractSpectralCoordinate)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "Shack-Hartmann electric-field input must declare a monochromatic channel"))
end

function _require_sh_front_end_input(front_end::ShackHartmannOpticalFrontEnd,
    input::PupilFunction)
    front_end.source === nothing && throw(WFSPreparationError(
        :optical_formation, :radiometry,
        "PupilFunction formation requires a source"))
    _require_sh_pupil_semantics(input, :optical_formation)
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

function _require_sh_front_end_input(front_end::ShackHartmannOpticalFrontEnd,
    input::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :optical_formation, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    _require_sh_metric_coordinates(input.metadata.coordinate_domain,
        :optical_formation)
    _require_sh_default_orientation(input.metadata.orientation,
        :optical_formation)
    _require_sh_centered_geometry(input.metadata, "electric-field input",
        :optical_formation)
    _require_sh_field_normalization(input.metadata.normalization)
    _require_sh_field_measure(input.metadata.spatial_measure)
    _require_sh_coherent_input(input.metadata.coherence, :optical_formation)
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
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "electric-field input must carry photon-rate normalization"))
end

@inline _require_sh_field_measure(::CellIntegratedMeasure) = nothing

function _require_sh_field_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "electric-field input must carry cell-integrated photon rate"))
end

function _require_sh_stage_domains(front_end::ShackHartmannOpticalFrontEnd,
    input, output::IntensityMap)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "Shack-Hartmann input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:optical_formation, :device,
            "Shack-Hartmann input and output occupy different devices"))
    propagation_storage = front_end.propagation.fft_stack
    typeof(input.metadata.backend) === typeof(backend(propagation_storage)) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "Shack-Hartmann input and microlens propagation backends differ"))
    input.metadata.device == plane_device(propagation_storage) ||
        throw(WFSPreparationError(:optical_formation, :device,
            "Shack-Hartmann input and microlens propagation occupy different devices"))
    typeof(input.metadata.backend) === typeof(backend(front_end.layout.valid_mask)) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "Shack-Hartmann input and subaperture-layout backends differ"))
    input.metadata.device == plane_device(front_end.layout.valid_mask) ||
        throw(WFSPreparationError(:optical_formation, :device,
            "Shack-Hartmann input and subaperture layout occupy different devices"))
    propagation_type = eltype(front_end.propagation.intensity)
    _sh_input_numeric_type(input) === propagation_type ||
        throw(WFSPreparationError(:optical_formation, :numeric_type,
            "Shack-Hartmann input type must match prepared propagation precision"))
    output.metadata.numeric_type === propagation_type ||
        throw(WFSPreparationError(:optical_formation, :numeric_type,
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

function _require_sh_rate_mosaic(front_end::ShackHartmannOpticalFrontEnd,
    output::IntensityMap, wavelength_m::Real, pupil_diameter_m::Real)
    n_sub = n_lenslets(front_end)
    n_pix = front_end.propagation.sampled_n_pix_subap
    size(output.values) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "Shack-Hartmann rate output must be an n_lenslets-by-n_lenslets tiled spot mosaic"))
    _require_sh_rate_coordinates(output.metadata.coordinate_domain)
    _require_sh_rate_measure(output.metadata.spatial_measure)
    _require_sh_default_orientation(output.metadata.orientation,
        :optical_formation)
    _require_sh_centered_geometry(output.metadata, "rate output",
        :optical_formation)
    T = eltype(front_end.propagation.intensity)
    expected_scale_arcsec = sh_pixel_scale_init(
        pupil_diameter_m / n_sub,
        front_end.propagation.effective_padding, wavelength_m) *
        front_end.propagation.binning_pixel_scale
    expected_scale_rad = T(expected_scale_arcsec / ARCSEC_PER_RAD)
    output.metadata.sampling == (expected_scale_rad, expected_scale_rad) ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "Shack-Hartmann rate output sampling does not match prepared microlens propagation"))
    _require_sh_output_wavelength(output.metadata.spectral, wavelength_m)
    return nothing
end


@inline _require_sh_rate_coordinates(::AngularCoordinates) = nothing

function _require_sh_rate_coordinates(::AbstractPlaneCoordinateDomain)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "Shack-Hartmann rate output must use angular coordinates"))
end

@inline _require_sh_rate_measure(::CellIntegratedMeasure) = nothing

function _require_sh_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "Shack-Hartmann rate output must contain cell-integrated photon rate"))
end


function _require_sh_output_wavelength(channel::MonochromaticChannel,
    wavelength_m::Real)
    channel.wavelength_m == typeof(channel.wavelength_m)(wavelength_m) ||
        throw(WFSPreparationError(
        :optical_formation, :plane_metadata,
        "Shack-Hartmann rate output wavelength does not match its optical input"))
    return nothing
end


function _require_sh_output_wavelength(::AbstractSpectralCoordinate,
    ::Real)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "Shack-Hartmann rate output must declare a monochromatic channel"))
end

function prepare_wfs_optical_formation(
    front_end::ShackHartmannOpticalFrontEnd, input,
    output::IntensityMap)
    validate_wfs_optical_input(input)
    validate_wfs_optical_products(output)
    _require_sh_front_end_input(front_end, input)
    _require_sh_stage_domains(front_end, input, output)
    resolution, pupil_diameter = _require_sh_layout_geometry(
        front_end.layout, n_lenslets(front_end), input,
        :optical_formation)
    wavelength_m = _sh_front_end_wavelength(front_end, input)
    _prepare_microlens_sampling_wavelength!(front_end, resolution,
        _sh_pupil_diameter(input.metadata), wavelength_m)
    _prepare_sh_source_workspace!(execution_style(_sh_input_storage(input)),
        front_end, front_end.source)
    _require_sh_rate_mosaic(front_end, output, wavelength_m,
        pupil_diameter)
    return PreparedShackHartmannOpticalFormation(front_end, input, output,
        _sh_sampling_signature(front_end))
end

@inline _prepare_sh_source_workspace!(::ExecutionStyle,
    ::ShackHartmannOpticalFrontEnd, source) = nothing

function _prepare_sh_source_workspace!(::AcceleratorStyle,
    front_end::ShackHartmannOpticalFrontEnd, source::Asterism)
    sh_stacked_asterism_compatible(source) &&
        ensure_sh_asterism_buffers!(front_end, length(source.sources))
    return nothing
end

function prepare_wfs_optical_formation(
    front_end::ShackHartmannOpticalFrontEnd{
        <:Any,<:Any,<:Any,<:Any,<:SpectralSource}, input,
    output::OpticalProductBundle)
    validate_wfs_optical_input(input)
    validate_wfs_optical_products(output)
    source = front_end.source
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "spectral Shack-Hartmann output count must match its spectral samples"))
    input_type = _sh_input_numeric_type(input)
    T = promote_type(input_type, typeof(photon_irradiance(source)))
    first_sample = first(samples)
    first_component = ShackHartmannSpectralComponent(source.source,
        T(first_sample.wavelength),
        T(photon_irradiance(source)) * T(first_sample.weight))
    first_plan = prepare_wfs_optical_formation(
        _sh_front_end_with_source(front_end, first_component),
        input, output[1])
    plans = Vector{typeof(first_plan)}(undef, length(samples))
    plans[1] = first_plan
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = _sh_front_end_with_source(front_end, component)
        plans[index] = prepare_wfs_optical_formation(component_front_end,
            input, output[index])
        plans[index].sampling_signature == first_plan.sampling_signature ||
            throw(WFSPreparationError(:optical_formation, :plane_count,
                "spectral Shack-Hartmann components require separate microlens workspaces when their prepared FFT sampling differs"))
    end
    return PreparedShackHartmannOpticalBundleFormation(plans, input, output)
end

@inline function _sh_front_end_with_source(
    front_end::ShackHartmannOpticalFrontEnd, source)
    return ShackHartmannOpticalFrontEnd(front_end.microlens_array,
        front_end.propagation, front_end.layout, source;
        threshold_convolution=front_end.threshold_convolution)
end

@inline function _require_sh_optical_binding(output::IntensityMap, input,
    plan::PreparedShackHartmannOpticalFormation)
    output.metadata === plan.output.metadata &&
        output.values === plan.output.values ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "Shack-Hartmann rate output does not match prepared storage"))
    input === plan.input || throw(WFSPreparationError(:optical_formation,
        :prepared_binding,
        "Shack-Hartmann pupil input does not match prepared storage"))
    _sh_sampling_signature(plan.front_end) == plan.sampling_signature ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "Shack-Hartmann microlens sampling no longer matches its prepared plan"))
    return nothing
end

@inline validate_wfs_optical_formation_binding(output::IntensityMap, input,
    plan::PreparedShackHartmannOpticalFormation) =
    _require_sh_optical_binding(output, input, plan)

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
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = sensor.propagation
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
        if sensor.layout.valid_mask[i, j]
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
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = sensor.propagation
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
        workspace.fft_stack, sensor.layout.valid_mask, input.amplitude,
        input.opd, workspace.phasor, amp_scale, opd_to_cycles, n_sub, sub,
        ox, oy, n, pad; ndrange=(pad, pad, n_sub, n_sub))
    synchronize_backend!(style)
    return _finish_sh_explicit_stack!(style, sensor, input, source)
end

function _form_sh_explicit_stack!(::ScalarCPUStyle,
    sensor::ShackHartmannOpticalFrontEnd, input::ElectricField, ::Nothing,
    wavelength_m)
    workspace = sensor.propagation
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    fill!(workspace.fft_stack, zero(eltype(workspace.fft_stack)))
    index = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if sensor.layout.valid_mask[i, j]
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
    sensor::ShackHartmannOpticalFrontEnd, input::ElectricField, ::Nothing,
    wavelength_m)
    workspace = sensor.propagation
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    launch_kernel_async!(style, sh_explicit_field_stack_kernel!,
        workspace.fft_stack, sensor.layout.valid_mask, input.values,
        workspace.phasor, n_sub, sub, ox, oy, n, pad;
        ndrange=(pad, pad, n_sub, n_sub))
    synchronize_backend!(style)
    return _finish_sh_explicit_stack!(style, sensor, input, nothing)
end

function _form_sh_explicit_asterism_serial!(style::ExecutionStyle,
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::Asterism)
    wavelength(source)
    workspace = sensor.propagation
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
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::Asterism) =
    _form_sh_explicit_asterism_serial!(style, sensor, input, source)

function _form_sh_explicit_asterism!(style::AcceleratorStyle,
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::Asterism)
    sh_stacked_asterism_compatible(source) ||
        return _form_sh_explicit_asterism_serial!(style, sensor, input,
            source)
    workspace = sensor.propagation
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
        fft_view, sensor.layout.valid_mask, input.amplitude, input.opd,
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
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

@inline _form_sh_explicit_stack!(style::AcceleratorStyle,
    sensor::ShackHartmannOpticalFrontEnd, input::PupilFunction,
    source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

function _finish_sh_explicit_stack!(style::ExecutionStyle,
    sensor::ShackHartmannOpticalFrontEnd, input, source)
    workspace = sensor.propagation
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
    ::ShackHartmannOpticalFrontEnd, input, source) = nothing

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannOpticalFrontEnd, input, source::LGSSource)
    return _apply_sh_lgs_spot_model!(lgs_profile(source), style, sensor,
        input, source, wavelength(source))
end

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannOpticalFrontEnd, input,
    source::ShackHartmannSpectralComponent{<:LGSSource})
    return _apply_sh_lgs_spot_model!(lgs_profile(source.source), style,
        sensor, input, source.source, wavelength(source))
end

function _apply_sh_lgs_spot_model!(::LGSProfileNone,
    ::ExecutionStyle, sensor::ShackHartmannOpticalFrontEnd, input,
    source::LGSSource,
    wavelength_m::Real)
    workspace = sensor.propagation
    workspace.elongation_kernel = apply_elongation_stack!(
        workspace.intensity_stack, lgs_elongation_factor(source),
        workspace.intensity_tmp_stack, workspace.elongation_kernel)
    return nothing
end

function _apply_sh_lgs_spot_model!(::LGSProfileNaProfile,
    ::ExecutionStyle, sensor::ShackHartmannOpticalFrontEnd, input,
    source::LGSSource,
    wavelength_m::Real)
    metadata = input.metadata
    ensure_lgs_kernels!(sensor, source, metadata.dimensions,
        _sh_pupil_diameter(metadata), metadata.sampling, metadata.origin,
        wavelength_m)
    workspace = sensor.propagation
    apply_lgs_convolution_stack!(workspace.intensity_stack,
        workspace.lgs_kernel_fft, workspace.fft_stack,
        workspace.fft_stack_plan, workspace.ifft_stack_plan)
    return nothing
end

function _sh_stack_intensity!(::ScalarCPUStyle,
    sensor::ShackHartmannOpticalFrontEnd,
    intensity_scale, ::Int)
    workspace = sensor.propagation
    @inbounds for index in axes(workspace.fft_stack, 3),
            y in axes(workspace.fft_stack, 2),
            x in axes(workspace.fft_stack, 1)
        workspace.intensity_stack[x, y, index] =
            abs2(workspace.fft_stack[x, y, index]) * intensity_scale
    end
    return workspace.intensity_stack
end

function _sh_stack_intensity!(style::AcceleratorStyle,
    sensor::ShackHartmannOpticalFrontEnd, intensity_scale, pad::Int)
    workspace = sensor.propagation
    launch_kernel!(style, complex_abs2_stack_kernel!,
        workspace.intensity_stack, workspace.fft_stack,
        intensity_scale, pad, n_lenslets(sensor)^2;
        ndrange=size(workspace.intensity_stack))
    return workspace.intensity_stack
end

function form_wfs_optical_products!(output::IntensityMap, input,
    plan::PreparedShackHartmannOpticalFormation)
    _require_sh_optical_binding(output, input, plan)
    front_end = plan.front_end
    wavelength_m = _sh_front_end_wavelength(front_end, input)
    _form_sh_explicit_stack!(execution_style(_sh_input_storage(input)),
        front_end, input, front_end.source,
        wavelength_m)
    shack_hartmann_detector_image!(output.values,
        front_end.propagation.sampled_spot_cube, n_lenslets(front_end))
    return output
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedShackHartmannOpticalBundleFormation)
    validate_wfs_optical_formation_binding(output, input, plan)
    @inbounds for index in eachindex(plan.plans)
        form_wfs_optical_products!(output[index], input, plan.plans[index])
    end
    return output
end

function validate_wfs_optical_formation_binding(
    output::OpticalProductBundle, input,
    plan::PreparedShackHartmannOpticalBundleFormation)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "spectral Shack-Hartmann products do not match prepared storage"))
    @inbounds for index in eachindex(plan.plans)
        validate_wfs_optical_formation_binding(output[index], input,
            plan.plans[index])
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
        sensor.acquisition.spot_cube, "observation")
    _require_sh_storage_domain(:estimation, measurement.metadata,
        sensor.estimator.slopes, "measurement")
    _require_sh_storage_domain(:estimation, observation.metadata,
        sensor.front_end.layout.valid_mask, "observation/layout")
    n_sub = n_lenslets(sensor)
    n_pix = sensor.front_end.propagation.sampled_n_pix_subap
    size(observation.storage) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann estimator requires a tiled lenslet mosaic"))
    size(measurement.storage) == size(sensor.estimator.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann measurement storage has the wrong slope shape"))
    calibration_binding = _prepare_sh_calibration_binding(sensor)
    ensure_sh_acquisition_buffers!(sensor, n_pix)
    return PreparedShackHartmannEstimator(sensor, observation, measurement,
        AcquiredObservationPath(), calibration_binding)
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

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation,
    plan::PreparedShackHartmannEstimator{
        <:Any,<:Any,<:Any,<:AcquiredObservationPath,<:Any})
    measurement === plan.measurement && observation === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator storage does not match its plan"))
    sensor = plan.sensor
    _require_sh_calibration_binding(sensor, plan.calibration_binding)
    n_sub = n_lenslets(sensor)
    n_pix = sensor.front_end.propagation.sampled_n_pix_subap
    style = execution_style(observation.storage)
    _unpack_sh_mosaic!(style, sensor.acquisition.spot_cube,
        observation.storage, n_sub, n_pix)
    peak = sh_safe_peak_value(sensor.acquisition.spot_cube)
    sh_signal_from_spots_calibrated!(sensor, peak,
        slope_extraction_model(sensor))
    copyto!(measurement.storage, sensor.estimator.slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedShackHartmannEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator storage does not match its plan"))
    return nothing
end

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
        sensor.estimator.slopes, "geometric input")
    _require_sh_storage_domain(:estimation, input.metadata,
        sensor.front_end.layout.valid_mask, "geometric input/layout")
    _require_sh_storage_domain(:estimation, measurement.metadata,
        sensor.estimator.slopes, "geometric measurement")
    size(measurement.storage) == size(sensor.estimator.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric Shack-Hartmann measurement storage has the wrong slope shape"))
    return PreparedShackHartmannEstimator(sensor, input, measurement,
        DirectMeasurementPath(), ShackHartmannLayoutBinding(
            subaperture_layout_revision(sensor.front_end.layout)))
end


function prepare_wfs_estimation(::ShackHartmannWFS{<:Geometric},
    ::ElectricField, ::WFSMeasurement)
    throw(WFSPreparationError(:estimation, :unsupported,
        "geometric Shack-Hartmann estimation requires OPD-bearing PupilFunction input"))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    plan::PreparedShackHartmannEstimator{
        <:Any,<:Any,<:Any,<:DirectMeasurementPath,<:Any})
    measurement === plan.measurement && input === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "geometric Shack-Hartmann estimator storage does not match its plan"))
    sensor = plan.sensor
    _require_sh_layout_binding(sensor, plan.calibration_binding)
    geometric_wavefront_slopes!(sensor.estimator.slopes, input.opd,
        sensor.front_end.layout.valid_mask, input.metadata.sampling)
    copyto!(measurement.storage, sensor.estimator.slopes)
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
        ShackHartmannOpticalFrontEnd(sensor.front_end, source), input,
        source)
end

function shack_hartmann_rate_map(
    front_end::ShackHartmannOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(front_end.propagation.intensity)
    first_sample = first(samples)
    first_component = ShackHartmannSpectralComponent(source.source,
        T(first_sample.wavelength),
        T(photon_irradiance(source)) * T(first_sample.weight))
    first_map = shack_hartmann_rate_map(
        _sh_front_end_with_source(front_end, first_component), input)
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        maps[index] = shack_hartmann_rate_map(
            _sh_front_end_with_source(front_end, component), input)
    end
    return OpticalProductBundle(maps)
end


function shack_hartmann_rate_map(sensor::ShackHartmannWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    front_end = source === nothing ? sensor.front_end :
        ShackHartmannOpticalFrontEnd(sensor.front_end, source)
    return shack_hartmann_rate_map(front_end, input)
end

function shack_hartmann_rate_map(front_end::ShackHartmannOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    resolved_front_end = source === nothing ? front_end :
        _sh_front_end_with_source(front_end, source)
    wavelength_m = _sh_front_end_wavelength(resolved_front_end, input)
    _prepare_microlens_sampling_wavelength!(resolved_front_end,
        input.metadata.dimensions[1],
        _sh_pupil_diameter(input.metadata), wavelength_m)
    propagation = resolved_front_end.propagation
    n = n_lenslets(resolved_front_end) * propagation.sampled_n_pix_subap
    T = eltype(propagation.intensity)
    values = similar(_sh_input_storage(input), T, n, n)
    fill!(values, zero(T))
    pixel_scale_arcsec = sh_pixel_scale_init(
        _sh_pupil_diameter(input.metadata) / n_lenslets(resolved_front_end),
        propagation.effective_padding, wavelength_m) *
        propagation.binning_pixel_scale
    pixel_scale_rad = T(pixel_scale_arcsec / ARCSEC_PER_RAD)
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
