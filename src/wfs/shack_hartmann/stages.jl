#
# Prepared Shack-Hartmann stage composition
#

"""Bind a diffractive Shack-Hartmann optical front end to its illumination."""
struct ShackHartmannOpticalFrontEnd{W<:ShackHartmannWFS,S}
    sensor::W
    source::S
end

function ShackHartmannOpticalFrontEnd(
    sensor::ShackHartmannWFS{<:Diffractive}, source::AbstractSource)
    return ShackHartmannOpticalFrontEnd{typeof(sensor),typeof(source)}(
        sensor, source)
end

function ShackHartmannOpticalFrontEnd(
    sensor::ShackHartmannWFS{<:Diffractive})
    return ShackHartmannOpticalFrontEnd{typeof(sensor),Nothing}(sensor,
        nothing)
end

function ShackHartmannOpticalFrontEnd(
    ::ShackHartmannWFS{<:Geometric}, source=nothing)
    throw(WFSPreparationError(:optical_formation, :unsupported,
        "geometric Shack-Hartmann sensing uses DirectMeasurementPath and has no optical front end"))
end

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

struct PreparedShackHartmannAcquisition{D,P,O}
    detector::D
    detector_plan::P
    observation::O
end

struct PreparedShackHartmannEstimator{W,I,M,P<:AbstractWFSMeasurementPath}
    sensor::W
    input::I
    measurement::M
    path::P
end

@inline wfs_measurement_path(plan::PreparedShackHartmannEstimator) = plan.path

@inline function _sh_pupil_diameter(metadata::OpticalPlaneMetadata)
    return metadata.sampling[1] * metadata.dimensions[1]
end

@inline function _sh_sampling_signature(sensor::ShackHartmannWFS)
    workspace = sensor.optical_workspace
    return (size(workspace.fft_stack), workspace.effective_padding,
        workspace.binning_pixel_scale, workspace.sampled_n_pix_subap,
        workspace.phasor_ratio)
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
    return nothing
end

function _require_sh_front_end_input(front_end::ShackHartmannOpticalFrontEnd,
    input::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :optical_formation, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    _require_sh_field_normalization(input.metadata.normalization)
    _require_sh_field_measure(input.metadata.spatial_measure)
    return nothing
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

function _require_sh_stage_domains(input, output::IntensityMap)
    typeof(input.metadata.backend) === typeof(output.metadata.backend) ||
        throw(WFSPreparationError(:optical_formation, :backend,
            "Shack-Hartmann input and output backends differ"))
    input.metadata.device == output.metadata.device ||
        throw(WFSPreparationError(:optical_formation, :device,
            "Shack-Hartmann input and output occupy different devices"))
    return nothing
end

function _require_sh_rate_mosaic(sensor::ShackHartmannWFS,
    output::IntensityMap, wavelength_m::Real)
    n_sub = n_lenslets(sensor)
    n_pix = sensor.optical_workspace.sampled_n_pix_subap
    size(output.values) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "Shack-Hartmann rate output must be an n_lenslets-by-n_lenslets tiled spot mosaic"))
    _require_sh_output_wavelength(output.metadata.spectral, wavelength_m)
    return nothing
end


function _require_sh_output_wavelength(channel::MonochromaticChannel,
    wavelength_m::Real)
    channel.wavelength_m == wavelength_m || throw(WFSPreparationError(
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
    _require_sh_stage_domains(input, output)
    sensor = front_end.sensor
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions[2] == resolution ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "Shack-Hartmann pupil input must be square"))
    input.metadata.sampling[1] == input.metadata.sampling[2] ||
        throw(WFSPreparationError(:optical_formation, :plane_metadata,
            "Shack-Hartmann pupil sampling must be square"))
    wavelength_m = _sh_front_end_wavelength(front_end, input)
    _prepare_microlens_sampling_wavelength!(sensor, resolution,
        _sh_pupil_diameter(input.metadata), wavelength_m)
    _prepare_sh_source_workspace!(execution_style(_sh_input_storage(input)),
        sensor, front_end.source)
    _require_sh_rate_mosaic(sensor, output, wavelength_m)
    return PreparedShackHartmannOpticalFormation(front_end, input, output,
        _sh_sampling_signature(sensor))
end

@inline _prepare_sh_source_workspace!(::ExecutionStyle,
    ::ShackHartmannWFS, source) = nothing

function _prepare_sh_source_workspace!(::AcceleratorStyle,
    sensor::ShackHartmannWFS, source::Asterism)
    sh_stacked_asterism_compatible(source) &&
        ensure_sh_asterism_buffers!(sensor, length(source.sources))
    return nothing
end

function prepare_wfs_optical_formation(
    front_end::ShackHartmannOpticalFrontEnd{<:Any,<:SpectralSource}, input,
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
        ShackHartmannOpticalFrontEnd(front_end.sensor, first_component),
        input, output[1])
    plans = Vector{typeof(first_plan)}(undef, length(samples))
    plans[1] = first_plan
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = ShackHartmannOpticalFrontEnd(
            front_end.sensor, component)
        plans[index] = prepare_wfs_optical_formation(component_front_end,
            input, output[index])
        plans[index].sampling_signature == first_plan.sampling_signature ||
            throw(WFSPreparationError(:optical_formation, :plane_count,
                "spectral Shack-Hartmann components require separate microlens workspaces when their prepared FFT sampling differs"))
    end
    return PreparedShackHartmannOpticalBundleFormation(plans, input, output)
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
    _sh_sampling_signature(plan.front_end.sensor) == plan.sampling_signature ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "Shack-Hartmann microlens sampling no longer matches its prepared plan"))
    return nothing
end

@kernel function sh_explicit_pupil_stack_kernel!(fft_stack, valid_mask,
    amplitude, opd, phasor, amp_scale, opd_to_cycles, n_sub::Int, sub::Int,
    ox::Int, oy::Int, n::Int, pad::Int)
    x, y, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && i <= n_sub && j <= n_sub
        index = (i - 1) * n_sub + j
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
        index = (i - 1) * n_sub + j
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
        spot_index = (i - 1) * n_sub + j
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
    sensor::ShackHartmannWFS, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = sensor.optical_workspace
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
    @inbounds for i in 1:n_sub, j in 1:n_sub
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
    sensor::ShackHartmannWFS, input::PupilFunction,
    source::AbstractSource, wavelength_m)
    workspace = sensor.optical_workspace
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
    sensor::ShackHartmannWFS, input::ElectricField, ::Nothing, wavelength_m)
    workspace = sensor.optical_workspace
    n = input.metadata.dimensions[1]
    n_sub = n_lenslets(sensor)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    fill!(workspace.fft_stack, zero(eltype(workspace.fft_stack)))
    index = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
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
    sensor::ShackHartmannWFS, input::ElectricField, ::Nothing, wavelength_m)
    workspace = sensor.optical_workspace
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
    sensor::ShackHartmannWFS, input::PupilFunction, source::Asterism)
    wavelength(source)
    workspace = sensor.optical_workspace
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
    sensor::ShackHartmannWFS, input::PupilFunction, source::Asterism) =
    _form_sh_explicit_asterism_serial!(style, sensor, input, source)

function _form_sh_explicit_asterism!(style::AcceleratorStyle,
    sensor::ShackHartmannWFS, input::PupilFunction, source::Asterism)
    sh_stacked_asterism_compatible(source) ||
        return _form_sh_explicit_asterism_serial!(style, sensor, input,
            source)
    workspace = sensor.optical_workspace
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
    sensor::ShackHartmannWFS, input::PupilFunction, source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

@inline _form_sh_explicit_stack!(style::AcceleratorStyle,
    sensor::ShackHartmannWFS, input::PupilFunction, source::Asterism,
    wavelength_m) = _form_sh_explicit_asterism!(style, sensor, input, source)

function _finish_sh_explicit_stack!(style::ExecutionStyle,
    sensor::ShackHartmannWFS, input, source)
    workspace = sensor.optical_workspace
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
    ::ShackHartmannWFS, input, source) = nothing

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannWFS, input, source::LGSSource)
    return _apply_sh_lgs_spot_model!(lgs_profile(source), style, sensor,
        input, source, wavelength(source))
end

@inline function _apply_sh_source_spot_model!(style::ExecutionStyle,
    sensor::ShackHartmannWFS, input,
    source::ShackHartmannSpectralComponent{<:LGSSource})
    return _apply_sh_lgs_spot_model!(lgs_profile(source.source), style,
        sensor, input, source.source, wavelength(source))
end

function _apply_sh_lgs_spot_model!(::LGSProfileNone,
    ::ExecutionStyle, sensor::ShackHartmannWFS, input, source::LGSSource,
    wavelength_m::Real)
    workspace = sensor.optical_workspace
    workspace.elongation_kernel = apply_elongation_stack!(
        workspace.intensity_stack, lgs_elongation_factor(source),
        workspace.intensity_tmp_stack, workspace.elongation_kernel)
    return nothing
end

function _apply_sh_lgs_spot_model!(::LGSProfileNaProfile,
    ::ExecutionStyle, sensor::ShackHartmannWFS, input, source::LGSSource,
    wavelength_m::Real)
    metadata = input.metadata
    ensure_lgs_kernels!(sensor, source, metadata.dimensions,
        _sh_pupil_diameter(metadata), metadata.sampling, metadata.origin,
        wavelength_m)
    workspace = sensor.optical_workspace
    apply_lgs_convolution_stack!(workspace.intensity_stack,
        workspace.lgs_kernel_fft, workspace.fft_stack,
        workspace.fft_stack_plan, workspace.ifft_stack_plan)
    return nothing
end

function _sh_stack_intensity!(::ScalarCPUStyle, sensor::ShackHartmannWFS,
    intensity_scale, ::Int)
    workspace = sensor.optical_workspace
    @inbounds for index in axes(workspace.fft_stack, 3),
            y in axes(workspace.fft_stack, 2),
            x in axes(workspace.fft_stack, 1)
        workspace.intensity_stack[x, y, index] =
            abs2(workspace.fft_stack[x, y, index]) * intensity_scale
    end
    return workspace.intensity_stack
end

function _sh_stack_intensity!(style::AcceleratorStyle,
    sensor::ShackHartmannWFS, intensity_scale, pad::Int)
    workspace = sensor.optical_workspace
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
    sensor = front_end.sensor
    wavelength_m = _sh_front_end_wavelength(front_end, input)
    _form_sh_explicit_stack!(execution_style(_sh_input_storage(input)),
        sensor, input, front_end.source,
        wavelength_m)
    shack_hartmann_detector_image!(output.values,
        sensor.optical_workspace.sampled_spot_cube, n_lenslets(sensor))
    return output
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedShackHartmannOpticalBundleFormation)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "spectral Shack-Hartmann products do not match prepared storage"))
    @inbounds for index in eachindex(plan.plans)
        form_wfs_optical_products!(output[index], input, plan.plans[index])
    end
    return output
end

function prepare_wfs_acquisition(detector::Detector,
    optical_product::IntensityMap, observation::WFSObservation)
    validate_wfs_optical_products(optical_product)
    validate_wfs_observation(observation)
    plan = prepare_detector_acquisition(detector, optical_product)
    size(observation.storage) == size(output_frame(detector)) ||
        throw(WFSPreparationError(:acquisition, :shape,
            "WFS observation storage must match the prepared detector output"))
    typeof(backend(observation.storage)) === typeof(backend(detector)) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "WFS observation and detector backends differ"))
    plane_device(observation.storage) == plane_device(output_frame(detector)) ||
        throw(WFSPreparationError(:acquisition, :device,
            "WFS observation and detector output occupy different devices"))
    return PreparedShackHartmannAcquisition(detector, plan, observation)
end

function acquire_wfs_observation!(observation::WFSObservation,
    optical_product::IntensityMap,
    plan::PreparedShackHartmannAcquisition, rng::AbstractRNG)
    observation === plan.observation || throw(WFSPreparationError(
        :acquisition, :prepared_binding,
        "WFS observation does not match prepared storage"))
    frame = capture!(plan.detector, optical_product, plan.detector_plan, rng)
    copyto!(observation.storage, frame)
    return observation
end

function prepare_wfs_estimation(sensor::ShackHartmannWFS{<:Diffractive},
    observation::WFSObservation, measurement::WFSMeasurement)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    n_sub = n_lenslets(sensor)
    n_pix = sensor.optical_workspace.sampled_n_pix_subap
    ensure_sh_acquisition_buffers!(sensor, n_pix)
    size(observation.storage) == (n_sub * n_pix, n_sub * n_pix) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann estimator requires a tiled lenslet mosaic"))
    size(measurement.storage) == size(sensor.estimator.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Shack-Hartmann measurement storage has the wrong slope shape"))
    return PreparedShackHartmannEstimator(sensor, observation,
        measurement, AcquiredObservationPath())
end

@kernel function sh_unpack_mosaic_kernel!(spot_cube, mosaic, n_sub::Int,
    n_pix::Int)
    i, j, x, y = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub && x <= n_pix && y <= n_pix
        index = (i - 1) * n_sub + j
        @inbounds spot_cube[index, x, y] =
            mosaic[(i - 1) * n_pix + x, (j - 1) * n_pix + y]
    end
end

function _unpack_sh_mosaic!(::ScalarCPUStyle, spot_cube, mosaic,
    n_sub::Int, n_pix::Int)
    @inbounds for y in 1:n_pix, x in 1:n_pix, i in 1:n_sub, j in 1:n_sub
        index = (i - 1) * n_sub + j
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
    plan::PreparedShackHartmannEstimator{<:Any,<:Any,<:Any,<:AcquiredObservationPath})
    measurement === plan.measurement && observation === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Shack-Hartmann estimator storage does not match its plan"))
    sensor = plan.sensor
    n_sub = n_lenslets(sensor)
    n_pix = sensor.optical_workspace.sampled_n_pix_subap
    style = execution_style(observation.storage)
    _unpack_sh_mosaic!(style, sensor.acquisition.spot_cube,
        observation.storage, n_sub, n_pix)
    peak = sh_safe_peak_value(sensor.acquisition.spot_cube)
    sh_signal_from_spots_calibrated!(sensor, peak,
        slope_extraction_model(sensor))
    copyto!(measurement.storage, sensor.estimator.slopes)
    return measurement
end

function prepare_wfs_estimation(sensor::ShackHartmannWFS{<:Geometric},
    input::PupilFunction, measurement::WFSMeasurement)
    validate_wfs_optical_input(input)
    validate_wfs_measurement(measurement)
    size(measurement.storage) == size(sensor.estimator.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric Shack-Hartmann measurement storage has the wrong slope shape"))
    return PreparedShackHartmannEstimator(sensor, input, measurement,
        DirectMeasurementPath())
end


function prepare_wfs_estimation(::ShackHartmannWFS{<:Geometric},
    ::ElectricField, ::WFSMeasurement)
    throw(WFSPreparationError(:estimation, :unsupported,
        "geometric Shack-Hartmann estimation requires OPD-bearing PupilFunction input"))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    plan::PreparedShackHartmannEstimator{<:Any,<:Any,<:Any,<:DirectMeasurementPath})
    measurement === plan.measurement && input === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "geometric Shack-Hartmann estimator storage does not match its plan"))
    sensor = plan.sensor
    geometric_slopes!(sensor.estimator.slopes, input.opd,
        sensor.layout.valid_mask)
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
    samples = spectral_bundle(source).samples
    T = promote_type(_sh_input_numeric_type(input),
        typeof(photon_irradiance(source)))
    first_sample = first(samples)
    first_component = ShackHartmannSpectralComponent(source.source,
        T(first_sample.wavelength),
        T(photon_irradiance(source)) * T(first_sample.weight))
    first_map = shack_hartmann_rate_map(sensor, input, first_component)
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        sample = samples[index]
        component = ShackHartmannSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        maps[index] = shack_hartmann_rate_map(sensor, input, component)
    end
    return OpticalProductBundle(maps)
end


function shack_hartmann_rate_map(sensor::ShackHartmannWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    front_end = source === nothing ? ShackHartmannOpticalFrontEnd(sensor) :
        ShackHartmannOpticalFrontEnd(sensor, source)
    wavelength_m = _sh_front_end_wavelength(front_end, input)
    _prepare_microlens_sampling_wavelength!(sensor,
        input.metadata.dimensions[1],
        _sh_pupil_diameter(input.metadata), wavelength_m)
    n = n_lenslets(sensor) * sensor.optical_workspace.sampled_n_pix_subap
    T = _sh_input_numeric_type(input)
    values = similar(_sh_input_storage(input), T, n, n)
    fill!(values, zero(T))
    pixel_scale_arcsec = sh_pixel_scale_init(
        _sh_pupil_diameter(input.metadata) / n_lenslets(sensor),
        sensor.optical_workspace.effective_padding, wavelength_m) *
        sensor.optical_workspace.binning_pixel_scale
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
