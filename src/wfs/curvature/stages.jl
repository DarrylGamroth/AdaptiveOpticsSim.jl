#
# Prepared Curvature WFS stages
#

@kernel function curvature_branch_field_from_pupil_kernel!(field_stack,
    amplitude, opd, defocus_stack, phasor, amplitude_scale, opd_to_cycles,
    ox::Int, oy::Int, n::Int, pad::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= size(field_stack, 3)
        xi = x - ox
        yi = y - oy
        value = zero(eltype(field_stack))
        if 1 <= xi <= n && 1 <= yi <= n
            value = amplitude_scale * @inbounds(amplitude[xi, yi]) *
                cispi(opd_to_cycles * @inbounds(opd[xi, yi]))
        end
        @inbounds field_stack[x, y, branch] = value *
            defocus_stack[x, y, branch] * phasor[x, y]
    end
end

@kernel function curvature_reduce_observation_pair_kernel!(plus_out,
    minus_out, plus_in, minus_in, plus_factor::Int, minus_factor::Int,
    plus_scale, minus_scale, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        T = eltype(plus_out)
        plus = zero(T)
        minus = zero(T)
        @inbounds for jj in 1:plus_factor, ii in 1:plus_factor
            plus += T(plus_in[(i - 1) * plus_factor + ii,
                (j - 1) * plus_factor + jj])
        end
        @inbounds for jj in 1:minus_factor, ii in 1:minus_factor
            minus += T(minus_in[(i - 1) * minus_factor + ii,
                (j - 1) * minus_factor + jj])
        end
        @inbounds begin
            plus_out[i, j] = plus * T(plus_scale)
            minus_out[i, j] = minus * T(minus_scale)
        end
    end
end

@kernel function curvature_unpack_channel_pair_kernel!(plus_out,
    minus_out, channels, plus_scale, minus_scale, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        T = eltype(plus_out)
        @inbounds begin
            plus_out[i, j] = T(channels[1, idx]) * T(plus_scale)
            minus_out[i, j] = T(channels[2, idx]) * T(minus_scale)
        end
    end
end

struct CurvaturePropagationBinding{H,F,D,I,C,P,M,Q}
    phasor::H
    field_stack::F
    defocus_stack::D
    intensity_stack::I
    cropped_plus::C
    cropped_minus::C
    fft_plan::P
    frame_plus::M
    frame_minus::Q
    revision::UInt
end

struct PreparedCurvatureOpticalFormation{F,I,O,B}
    front_end::F
    input::I
    output::O
    binding::B
end

@inline modulated_wfs_propagation_storage(
    front_end::CurvatureOpticalFrontEnd) = front_end.propagation.field_stack

"""
    CurvaturePackedAcquisition(detector; readout_model, source,
        branch_durations)

Declare a single-detector mapping for the two Curvature branch rate planes.
`branch_durations` is a preparation-time assertion only; the detector remains
the sole owner that applies elapsed time.
"""
struct CurvaturePackedAcquisition{D,R,S,T<:AbstractFloat}
    detector::D
    readout_model::R
    source::S
    branch_durations::NTuple{2,T}
end

@inline _curvature_detector_duration(detector::Detector) =
    detector.params.integration_time
@inline _curvature_detector_duration(detector::AbstractCountingDetector) =
    counting_integration_time(detector)

@inline _curvature_branch_durations(::Nothing, duration::T) where {
    T<:AbstractFloat,
} = (duration, duration)

function _curvature_branch_durations(values::Tuple{A,B}, duration::T) where {
    A<:Real,B<:Real,T<:AbstractFloat,
}
    return (T(values[1]), T(values[2]))
end

function _curvature_branch_durations(::Any, ::AbstractFloat)
    throw(InvalidConfiguration(
        "branch_durations must be a two-element tuple of real durations"))
end

function CurvaturePackedAcquisition(detector::AbstractDetector;
    readout_model::CurvatureReadoutModel=CurvatureFrameReadout(),
    source=nothing, branch_durations=nothing)
    duration = _curvature_detector_duration(detector)
    T = typeof(duration)
    durations = _curvature_branch_durations(branch_durations, duration)
    return CurvaturePackedAcquisition{typeof(detector),
        typeof(readout_model),typeof(source),T}(detector, readout_model,
        source, durations)
end

struct PreparedCurvaturePackedFrameAcquisition{M,I,O,R,P,T}
    model::M
    optical_products::I
    observation::O
    packed_rate::R
    detector_plan::P
    detector_duration::T
end

struct PreparedCurvaturePackedCountingAcquisition{M,I,O,C,A,F,T}
    model::M
    optical_products::I
    observation::O
    channels::C
    detector_input::A
    detector_output::F
    detector_duration::T
end

struct CurvatureCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

abstract type AbstractPreparedCurvatureObservationMapping end

struct CurvatureImagePairMapping{P,M} <:
    AbstractPreparedCurvatureObservationMapping
    plus::P
    minus::M
    plus_reduction::Int
    minus_reduction::Int
end

struct CurvatureChannelPairMapping{C} <:
    AbstractPreparedCurvatureObservationMapping
    channels::C
end

struct PreparedCurvatureEstimator{W,I,M,P<:AbstractWFSMeasurementPath,C,G,T}
    sensor::W
    input::I
    measurement::M
    path::P
    calibration_binding::C
    mapping::G
    branch_rate_scales::NTuple{2,T}
end

@inline wfs_measurement_path(plan::PreparedCurvatureEstimator) = plan.path

function CurvatureOpticalFrontEnd(sensor::CurvatureWFS, source=nothing)
    front_end = sensor.front_end
    return CurvatureOpticalFrontEnd(front_end.defocus_pair,
        front_end.propagation, front_end.pupil_resolution,
        front_end.pupil_samples, front_end.readout_crop_resolution,
        front_end.readout_pixels_per_sample, front_end.branch_response,
        source)
end

@inline function curvature_branch_dimensions(
    front_end::CurvatureOpticalFrontEnd)
    side = front_end.pupil_samples * front_end.readout_pixels_per_sample
    return (side, side)
end

function _require_curvature_front_end_source(
    front_end::CurvatureOpticalFrontEnd, ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:optical_formation,
        :radiometry, "Curvature PupilFunction formation requires a source"))
    require_leaf_source(source, "prepared Curvature optical formation")
    return source
end

function _require_curvature_front_end_source(
    front_end::CurvatureOpticalFrontEnd, ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :optical_formation, :radiometry,
        "photon-rate ElectricField input must not also supply a Curvature source"))
    return nothing
end

@inline _curvature_front_end_wavelength(front_end::CurvatureOpticalFrontEnd,
    input::PupilFunction) = modulated_input_wavelength(input,
        front_end.source)
@inline _curvature_front_end_wavelength(::CurvatureOpticalFrontEnd,
    input::ElectricField) = modulated_input_wavelength(input)

function _require_curvature_input_geometry(
    front_end::CurvatureOpticalFrontEnd, input::PupilFunction)
    input.metadata.dimensions == (front_end.pupil_resolution,
        front_end.pupil_resolution) || throw(WFSPreparationError(
        :optical_formation, :shape,
        "Curvature pupil input dimensions differ from the prepared relay"))
    return nothing
end

function _require_curvature_input_geometry(
    front_end::CurvatureOpticalFrontEnd, input::ElectricField)
    input.metadata.dimensions == size(front_end.propagation.phasor) || throw(
        WFSPreparationError(:optical_formation, :shape,
            "Curvature ElectricField dimensions differ from the prepared diffraction grid"))
    return nothing
end

@inline _require_curvature_rate_coordinates(
    ::NormalizedPupilCoordinates, ::AbstractString) = nothing

function _require_curvature_rate_coordinates(
    ::AbstractPlaneCoordinateDomain, label::AbstractString)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "$label must use normalized pupil coordinates"))
end

@inline _require_curvature_rate_measure(
    ::CellIntegratedMeasure, ::AbstractString) = nothing

function _require_curvature_rate_measure(
    ::AbstractSpatialMeasure, label::AbstractString)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "$label must carry cell-integrated photon rate"))
end

function _require_curvature_rate_wavelength(channel::MonochromaticChannel,
    wavelength_m, label::AbstractString)
    channel.wavelength_m == wavelength_m || throw(
        WFSPreparationError(:optical_formation, :plane_metadata,
            "$label wavelength differs from its input"))
    return nothing
end

function _require_curvature_rate_wavelength(
    ::AbstractSpectralCoordinate, ::Any, label::AbstractString)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "$label wavelength differs from its input"))
end

function _require_curvature_rate_product(product::IntensityMap,
    expected_dimensions, wavelength_m, label::AbstractString)
    validate_wfs_optical_products(product)
    _require_curvature_rate_coordinates(product.metadata.coordinate_domain,
        label)
    _require_curvature_rate_measure(product.metadata.spatial_measure, label)
    size(product.values) == expected_dimensions || throw(
        WFSPreparationError(:optical_formation, :shape,
            "$label has the wrong prepared dimensions"))
    _require_curvature_rate_wavelength(product.metadata.spectral,
        wavelength_m, label)
    return product
end

function prepare_wfs_optical_formation(front_end::CurvatureOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    output::Tuple{<:IntensityMap,<:IntensityMap})
    require_modulated_wfs_input(input)
    _require_curvature_front_end_source(front_end, input)
    _require_curvature_input_geometry(front_end, input)
    wavelength_m = _curvature_front_end_wavelength(front_end, input)
    dimensions = curvature_branch_dimensions(front_end)
    _require_curvature_rate_product(output[1], dimensions, wavelength_m,
        "positive-defocus Curvature product")
    _require_curvature_rate_product(output[2], dimensions, wavelength_m,
        "negative-defocus Curvature product")
    output[1].values === output[2].values && throw(WFSPreparationError(
        :optical_formation, :prepared_binding,
        "Curvature branch products require distinct storage"))
    require_modulated_wfs_domains(front_end, input, output[1])
    require_modulated_wfs_domains(front_end, input, output[2])
    T = eltype(front_end.propagation.intensity_stack)
    output[1].metadata.numeric_type === T &&
        output[2].metadata.numeric_type === T || throw(
        WFSPreparationError(:optical_formation, :numeric_type,
            "Curvature output precision differs from prepared propagation"))
    propagation = front_end.propagation
    binding = CurvaturePropagationBinding(propagation.phasor,
        propagation.field_stack, propagation.defocus_stack,
        propagation.intensity_stack, propagation.cropped_plus,
        propagation.cropped_minus, propagation.fft_stack_plan,
        propagation.frame_plus, propagation.frame_minus,
        propagation.revision)
    return PreparedCurvatureOpticalFormation(front_end, input, output,
        binding)
end

function curvature_rate_maps(sensor::CurvatureWFS,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return curvature_rate_maps(CurvatureOpticalFrontEnd(sensor, source),
        input)
end

function curvature_rate_maps(front_end::CurvatureOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    wavelength_m = _curvature_front_end_wavelength(front_end, input)
    dimensions = curvature_branch_dimensions(front_end)
    T = eltype(front_end.propagation.intensity_stack)
    sampling_value = T(front_end.readout_crop_resolution /
        (front_end.pupil_resolution * dimensions[1]))
    function branch_map()
        values = similar(_modulated_input_storage(input), T, dimensions...)
        fill!(values, zero(T))
        metadata = OpticalPlaneMetadata(DetectorPlane(), values;
            coordinate_domain=NormalizedPupilCoordinates(),
            sampling=(sampling_value, sampling_value),
            spectral=MonochromaticChannel(T(wavelength_m)),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition())
        return IntensityMap(metadata, values)
    end
    return (branch_map(), branch_map())
end

function _require_curvature_optical_binding(
    plan::PreparedCurvatureOpticalFormation, input, output)
    input === plan.input && output === plan.output || throw(
        WFSPreparationError(:optical_formation, :prepared_binding,
            "Curvature optical products do not match their prepared plan"))
    propagation = plan.front_end.propagation
    binding = plan.binding
    propagation.phasor === binding.phasor &&
        propagation.field_stack === binding.field_stack &&
        propagation.defocus_stack === binding.defocus_stack &&
        propagation.intensity_stack === binding.intensity_stack &&
        propagation.cropped_plus === binding.cropped_plus &&
        propagation.cropped_minus === binding.cropped_minus &&
        propagation.fft_stack_plan === binding.fft_plan &&
        propagation.frame_plus === binding.frame_plus &&
        propagation.frame_minus === binding.frame_minus &&
        propagation.revision == binding.revision || throw(
        WFSPreparationError(:optical_formation, :prepared_binding,
            "Curvature propagation storage changed after preparation"))
    return nothing
end

@inline validate_wfs_optical_formation_binding(
    output::Tuple{<:IntensityMap,<:IntensityMap}, input,
    plan::PreparedCurvatureOpticalFormation) =
    _require_curvature_optical_binding(plan, input, output)

function _form_curvature_branch_fields!(::ScalarCPUStyle,
    front_end::CurvatureOpticalFrontEnd, input::PupilFunction)
    propagation = front_end.propagation
    T = eltype(propagation.intensity_stack)
    n = front_end.pupil_resolution
    pad = size(propagation.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    cell_area = T(input.metadata.sampling[1] * input.metadata.sampling[2])
    amplitude_scale = sqrt(T(photon_irradiance(front_end.source)) *
        cell_area)
    opd_to_cycles = T(2) / T(wavelength(front_end.source))
    fill!(propagation.field_stack, zero(eltype(propagation.field_stack)))
    @inbounds for y in 1:n, x in 1:n
        value = amplitude_scale * input.amplitude[x, y] *
            cispi(opd_to_cycles * input.opd[x, y])
        xx = ox + x
        yy = oy + y
        common = value * propagation.phasor[xx, yy]
        propagation.field_stack[xx, yy, 1] = common *
            propagation.defocus_stack[xx, yy, 1]
        propagation.field_stack[xx, yy, 2] = common *
            propagation.defocus_stack[xx, yy, 2]
    end
    return propagation.field_stack
end

function _form_curvature_branch_fields!(style::AcceleratorStyle,
    front_end::CurvatureOpticalFrontEnd, input::PupilFunction)
    propagation = front_end.propagation
    T = eltype(propagation.intensity_stack)
    n = front_end.pupil_resolution
    pad = size(propagation.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    cell_area = T(input.metadata.sampling[1] * input.metadata.sampling[2])
    amplitude_scale = sqrt(T(photon_irradiance(front_end.source)) *
        cell_area)
    opd_to_cycles = T(2) / T(wavelength(front_end.source))
    launch_kernel!(style, curvature_branch_field_from_pupil_kernel!,
        propagation.field_stack, input.amplitude, input.opd,
        propagation.defocus_stack, propagation.phasor, amplitude_scale,
        opd_to_cycles, ox, oy, n, pad;
        ndrange=size(propagation.field_stack))
    return propagation.field_stack
end

function _form_curvature_branch_fields!(style::ExecutionStyle,
    front_end::CurvatureOpticalFrontEnd, input::ElectricField)
    propagation = front_end.propagation
    return _form_curvature_field_input!(style, propagation, input)
end

function _form_curvature_field_input!(::ScalarCPUStyle, propagation,
    input::ElectricField)
    pad = size(propagation.field_stack, 1)
    @inbounds for branch in axes(propagation.field_stack, 3), y in 1:pad,
            x in 1:pad
        propagation.field_stack[x, y, branch] = input.values[x, y] *
            propagation.defocus_stack[x, y, branch] *
            propagation.phasor[x, y]
    end
    return propagation.field_stack
end

function _form_curvature_field_input!(style::AcceleratorStyle,
    propagation, input::ElectricField)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_branch_field_from_input_kernel!,
        propagation.field_stack, input.values, propagation.defocus_stack,
        propagation.phasor, size(propagation.field_stack, 1),
        size(propagation.field_stack, 3);
        ndrange=size(propagation.field_stack))
    finish_kernel_phase!(phase)
    return propagation.field_stack
end

function _sample_curvature_rate_planes!(::ScalarCPUStyle, output,
    front_end::CurvatureOpticalFrontEnd)
    propagation = front_end.propagation
    crop_n = front_end.readout_crop_resolution
    pad = size(propagation.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    factor = div(crop_n, size(output[1].values, 1))
    copyto!(propagation.cropped_plus,
        @view(propagation.intensity_stack[
            ox+1:ox+crop_n, oy+1:oy+crop_n, 1]))
    copyto!(propagation.cropped_minus,
        @view(propagation.intensity_stack[
            ox+1:ox+crop_n, oy+1:oy+crop_n, 2]))
    bin2d!(output[1].values, propagation.cropped_plus, factor)
    bin2d!(output[2].values, propagation.cropped_minus, factor)
    response = front_end.branch_response
    @. output[1].values = response.plus_throughput * output[1].values +
        response.plus_background
    @. output[2].values = response.minus_throughput * output[2].values +
        response.minus_background
    return output
end

function _sample_curvature_rate_planes!(style::AcceleratorStyle, output,
    front_end::CurvatureOpticalFrontEnd)
    propagation = front_end.propagation
    crop_n = front_end.readout_crop_resolution
    pad = size(propagation.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    factor = div(crop_n, size(output[1].values, 1))
    response = front_end.branch_response
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_sample_branch_kernel!,
        output[1].values, propagation.intensity_stack, ox, oy, factor, 1,
        response.plus_throughput, response.plus_background,
        size(output[1].values, 1), size(output[1].values, 2);
        ndrange=size(output[1].values))
    queue_kernel!(phase, curvature_sample_branch_kernel!,
        output[2].values, propagation.intensity_stack, ox, oy, factor, 2,
        response.minus_throughput, response.minus_background,
        size(output[2].values, 1), size(output[2].values, 2);
        ndrange=size(output[2].values))
    finish_kernel_phase!(phase)
    return output
end

function form_wfs_optical_products!(
    output::Tuple{<:IntensityMap,<:IntensityMap},
    input::Union{PupilFunction,ElectricField},
    plan::PreparedCurvatureOpticalFormation)
    validate_wfs_optical_formation_binding(output, input, plan)
    style = execution_style(output[1].values)
    _form_curvature_branch_fields!(style, plan.front_end, input)
    propagation = plan.front_end.propagation
    fraunhofer_intensity_stack!(propagation.intensity_stack,
        propagation.field_stack, propagation.fft_stack_plan)
    _sample_curvature_rate_planes!(style, output, plan.front_end)
    return output
end

function _require_curvature_branch_compatibility(products)
    plus = products[1].metadata
    minus = products[2].metadata
    plus.dimensions == minus.dimensions || throw(WFSPreparationError(
        :acquisition, :shape,
        "packed Curvature branches require identical dimensions"))
    typeof(plus.kind) === typeof(minus.kind) &&
        typeof(plus.coordinate_domain) === typeof(minus.coordinate_domain) &&
        plus.sampling == minus.sampling && plus.origin == minus.origin &&
        plus.centering == minus.centering &&
        plus.orientation == minus.orientation &&
        plus.spectral == minus.spectral || throw(WFSPreparationError(
        :acquisition, :plane_metadata,
        "packed Curvature branches require compatible plane geometry"))
    typeof(plus.normalization) === typeof(minus.normalization) &&
        typeof(plus.spatial_measure) === typeof(minus.spatial_measure) &&
        typeof(plus.coherence) === typeof(minus.coherence) || throw(
        WFSPreparationError(:acquisition, :radiometry,
            "packed Curvature branches require compatible radiometry"))
    plus.numeric_type === minus.numeric_type || throw(WFSPreparationError(
        :acquisition, :numeric_type,
        "packed Curvature branches require identical numeric types"))
    typeof(plus.backend) === typeof(minus.backend) || throw(
        WFSPreparationError(:acquisition, :backend,
            "packed Curvature branches occupy different backends"))
    plus.device == minus.device || throw(WFSPreparationError(:acquisition,
        :device, "packed Curvature branches occupy different devices"))
    products[1].values === products[2].values && throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "packed Curvature branches require distinct storage"))
    return nothing
end

function _require_curvature_packed_duration(model::CurvaturePackedAcquisition)
    duration = _curvature_detector_duration(model.detector)
    d1, d2 = model.branch_durations
    isfinite(d1) && isfinite(d2) && d1 > zero(d1) && d2 > zero(d2) ||
        throw(WFSPreparationError(:acquisition, :duration,
            "packed Curvature branch durations must be finite and positive"))
    d1 == d2 == duration || throw(WFSPreparationError(:acquisition,
        :duration,
        "packed Curvature branches must share the detector integration duration"))
    return duration
end

function _curvature_packed_rate_map(products)
    plus = products[1]
    n, m = size(plus.values)
    T = eltype(plus.values)
    values = similar(plus.values, T, 2 * n, m)
    fill!(values, zero(T))
    metadata = plus.metadata
    packed_metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=metadata.coordinate_domain,
        sampling=metadata.sampling,
        orientation=metadata.orientation,
        spectral=metadata.spectral,
        normalization=metadata.normalization,
        spatial_measure=metadata.spatial_measure,
        coherence=metadata.coherence)
    return IntensityMap(packed_metadata, values)
end

@inline _require_curvature_counting_measure(::CellIntegratedMeasure) = nothing

function _require_curvature_counting_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:acquisition, :radiometry,
        "packed counting Curvature inputs must carry cell-integrated photon rate"))
end

function prepare_wfs_acquisition(
    model::CurvaturePackedAcquisition{<:Detector,<:CurvatureFrameReadout},
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    observation::WFSObservation)
    validate_wfs_optical_products(optical_products)
    validate_wfs_observation(observation)
    _require_curvature_branch_compatibility(optical_products)
    duration = _require_curvature_packed_duration(model)
    isequal(observation.metadata.layout, :curvature_branch_regions) ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "packed frame Curvature observation requires :curvature_branch_regions layout"))
    packed_rate = _curvature_packed_rate_map(optical_products)
    detector_plan = prepare_wfs_acquisition(model.detector, packed_rate,
        observation; source=model.source)
    return PreparedCurvaturePackedFrameAcquisition(model,
        optical_products, observation, packed_rate, detector_plan, duration)
end

function prepare_wfs_acquisition(
    model::CurvaturePackedAcquisition{
        <:AbstractCountingDetector,<:CurvatureCountingReadout},
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    observation::WFSObservation)
    validate_wfs_optical_products(optical_products)
    validate_wfs_observation(observation)
    _require_curvature_branch_compatibility(optical_products)
    duration = _require_curvature_packed_duration(model)
    _require_curvature_counting_measure(
        optical_products[1].metadata.spatial_measure)
    isequal(observation.metadata.layout, :curvature_branch_channels) ||
        throw(WFSPreparationError(:acquisition, :detector_mapping,
            "packed counting Curvature observation requires :curvature_branch_channels layout"))
    n, m = size(optical_products[1].values)
    n == m || throw(WFSPreparationError(:acquisition, :shape,
        "packed counting Curvature branches must be square"))
    T = eltype(optical_products[1].values)
    channels = similar(optical_products[1].values, T, 2, n * n)
    fill!(channels, zero(T))
    detector = model.detector
    _require_counting_wfs_source(detector, model.source)
    _require_counting_wfs_spectral_match(optical_products[1].metadata,
        model.source)
    eltype(channels) === eltype(counting_array(detector)) || throw(
        WFSPreparationError(:acquisition, :numeric_type,
            "counting detector and packed Curvature precision differ"))
    typeof(backend(detector)) === typeof(backend(channels)) || throw(
        WFSPreparationError(:acquisition, :backend,
            "counting detector and packed Curvature channels use different backends"))
    plane_device(counting_array(detector)) == plane_device(channels) ||
        throw(WFSPreparationError(:acquisition, :device,
            "counting detector and packed Curvature channels occupy different devices"))
    ensure_buffers!(detector, size(channels))
    output = output_frame(detector)
    size(observation.storage) == size(output) || throw(WFSPreparationError(
        :acquisition, :shape,
        "packed counting Curvature observation must match detector output"))
    observation.metadata.numeric_type === eltype(output) || throw(
        WFSPreparationError(:acquisition, :numeric_type,
            "packed counting Curvature observation type must match detector output"))
    typeof(backend(observation.storage)) === typeof(backend(detector)) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "packed counting Curvature observation and detector backends differ"))
    plane_device(observation.storage) == plane_device(output) || throw(
        WFSPreparationError(:acquisition, :device,
            "packed counting Curvature observation and detector output occupy different devices"))
    return PreparedCurvaturePackedCountingAcquisition(model,
        optical_products, observation, channels, counting_array(detector),
        output, duration)
end

function prepare_wfs_acquisition(model::CurvaturePackedAcquisition,
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    observation::WFSObservation)
    throw(WFSPreparationError(:acquisition, :detector_mapping,
        "CurvatureFrameReadout requires a frame Detector and CurvatureCountingReadout requires a counting detector"))
end

function _pack_curvature_regions!(::ScalarCPUStyle, packed, plus, minus)
    n = size(plus, 1)
    @views copyto!(packed[1:n, :], plus)
    @views copyto!(packed[n+1:2*n, :], minus)
    return packed
end

function _pack_curvature_regions!(style::AcceleratorStyle, packed, plus,
    minus)
    launch_kernel!(style, curvature_frame_pack_kernel!, packed, plus,
        minus, size(plus, 1); ndrange=size(plus))
    return packed
end

function _pack_curvature_channels!(::ScalarCPUStyle, channels, plus,
    minus)
    n = size(plus, 1)
    @inbounds for i in 1:n, j in 1:n
        index = (i - 1) * n + j
        channels[1, index] = plus[i, j]
        channels[2, index] = minus[i, j]
    end
    return channels
end

function _pack_curvature_channels!(style::AcceleratorStyle, channels,
    plus, minus)
    launch_kernel!(style, curvature_channel_pack_kernel!, channels, plus,
        minus, size(plus, 1); ndrange=size(plus))
    return channels
end

function acquire_wfs_observation!(observation::WFSObservation,
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    plan::PreparedCurvaturePackedFrameAcquisition, rng::AbstractRNG)
    validate_wfs_acquisition_binding(observation, optical_products, plan)
    _pack_curvature_regions!(execution_style(plan.packed_rate.values),
        plan.packed_rate.values, optical_products[1].values,
        optical_products[2].values)
    acquire_wfs_observation!(observation, plan.packed_rate,
        plan.detector_plan, rng)
    return observation
end

function validate_wfs_acquisition_binding(observation::WFSObservation,
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    plan::PreparedCurvaturePackedFrameAcquisition)
    observation === plan.observation &&
        optical_products === plan.optical_products || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "packed Curvature frame storage does not match its plan"))
    _curvature_detector_duration(plan.model.detector) ==
        plan.detector_duration || throw(WFSPreparationError(:acquisition,
        :prepared_binding,
        "packed Curvature detector duration changed after preparation"))
    validate_wfs_acquisition_binding(observation, plan.packed_rate,
        plan.detector_plan)
    return nothing
end

function acquire_wfs_observation!(observation::WFSObservation,
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    plan::PreparedCurvaturePackedCountingAcquisition, rng::AbstractRNG)
    validate_wfs_acquisition_binding(observation, optical_products, plan)
    detector = plan.model.detector
    _pack_curvature_channels!(execution_style(plan.channels), plan.channels,
        optical_products[1].values, optical_products[2].values)
    frame = _capture_counting_wfs!(detector, plan.channels,
        plan.model.source, rng)
    copyto!(observation.storage, frame)
    return observation
end

function validate_wfs_acquisition_binding(observation::WFSObservation,
    optical_products::Tuple{<:IntensityMap,<:IntensityMap},
    plan::PreparedCurvaturePackedCountingAcquisition)
    observation === plan.observation &&
        optical_products === plan.optical_products || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "packed Curvature channel storage does not match its plan"))
    detector = plan.model.detector
    counting_array(detector) === plan.detector_input &&
        output_frame(detector) === plan.detector_output || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "packed Curvature counting storage changed after preparation"))
    _curvature_detector_duration(detector) == plan.detector_duration ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "packed Curvature detector duration changed after preparation"))
    return nothing
end

function _curvature_calibration_binding(sensor::CurvatureWFS)
    state = sensor.estimator.state
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "Curvature estimation requires explicit calibration"))
    return CurvatureCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength, state.calibration_signature,
        state.reference_signal_2d, state.valid_mask)
end

function _require_curvature_calibration(sensor::CurvatureWFS,
    binding::CurvatureCalibrationBinding)
    state = sensor.estimator.state
    state.calibrated &&
        state.calibration_revision == binding.revision &&
        state.calibration_wavelength == binding.wavelength_m &&
        state.calibration_signature == binding.signature &&
        state.reference_signal_2d === binding.reference_signal &&
        state.valid_mask === binding.valid_support || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Curvature calibration changed after estimator preparation"))
    return nothing
end

function _require_curvature_measurement(sensor::CurvatureWFS,
    measurement::WFSMeasurement)
    validate_wfs_measurement(measurement)
    isequal(measurement.metadata.kind, :curvature_signal) || throw(
        WFSPreparationError(:estimation, :estimator,
            "Curvature measurement kind must be :curvature_signal"))
    isequal(measurement.units, :dimensionless) || throw(
        WFSPreparationError(:estimation, :units,
            "Curvature differential signals are dimensionless"))
    measurement.metadata.numeric_type <: AbstractFloat || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Curvature measurement storage must be floating point"))
    measurement.metadata.numeric_type === eltype(sensor.estimator.state.slopes) ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "Curvature measurement precision differs from its estimator"))
    size(measurement.storage) == size(sensor.estimator.state.slopes) || throw(
        WFSPreparationError(:estimation, :shape,
            "Curvature measurement has the wrong signal-vector dimensions"))
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        sensor.estimator.state.slopes, "Curvature measurement")
    return measurement
end

function _require_curvature_image_observation(sensor::CurvatureWFS,
    observation::WFSObservation, label::AbstractString)
    validate_wfs_observation(observation)
    observation.metadata.numeric_type <: Real || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "$label requires real detector samples"))
    dimensions = observation.metadata.dimensions
    length(dimensions) == 2 && dimensions[1] == dimensions[2] || throw(
        WFSPreparationError(:estimation, :shape,
            "$label must be a square detector image"))
    dimensions[1] % sensor.params.pupil_samples == 0 || throw(
        WFSPreparationError(:estimation, :shape,
            "$label sampling must evenly reduce to the Curvature pupil grid"))
    _require_wfs_storage_domain(:estimation, observation.metadata,
        sensor.estimator.state.reduced_plus, label)
    return div(dimensions[1], sensor.params.pupil_samples)
end

function _curvature_rate_scales(sensor::CurvatureWFS,
    values::Tuple{A,B}) where {A<:Real,B<:Real}
    T = eltype(sensor.estimator.state.slopes)
    scales = (T(values[1]), T(values[2]))
    all(value -> isfinite(value) && value > zero(T), scales) || throw(
        WFSPreparationError(:estimation, :radiometry,
            "Curvature branch rate scales must be finite and positive"))
    return scales
end

function _curvature_rate_scales(::CurvatureWFS, ::Any)
    throw(WFSPreparationError(:estimation, :radiometry,
        "Curvature branch rate scales must be a two-element tuple of real values"))
end

function prepare_wfs_estimation(sensor::CurvatureWFS,
    observations::Tuple{<:WFSObservation,<:WFSObservation},
    measurement::WFSMeasurement; branch_rate_scales=(1, 1))
    _require_curvature_measurement(sensor, measurement)
    isequal(observations[1].metadata.layout, :curvature_branch_image) &&
        isequal(observations[2].metadata.layout,
            :curvature_branch_image) || throw(WFSPreparationError(
        :estimation, :detector_mapping,
        "separate Curvature observations require :curvature_branch_image layout"))
    plus_reduction = _require_curvature_image_observation(sensor,
        observations[1], "positive-defocus Curvature observation")
    minus_reduction = _require_curvature_image_observation(sensor,
        observations[2], "negative-defocus Curvature observation")
    mapping = CurvatureImagePairMapping(observations[1].storage,
        observations[2].storage, plus_reduction, minus_reduction)
    scales = _curvature_rate_scales(sensor, branch_rate_scales)
    binding = _curvature_calibration_binding(sensor)
    return PreparedCurvatureEstimator(sensor, observations, measurement,
        AcquiredObservationPath(), binding, mapping, scales)
end

function prepare_wfs_estimation(sensor::CurvatureWFS,
    observation::WFSObservation, measurement::WFSMeasurement;
    branch_rate_scales=(1, 1))
    _require_curvature_measurement(sensor, measurement)
    validate_wfs_observation(observation)
    observation.metadata.numeric_type <: Real || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Curvature observations require real detector samples"))
    layout = observation.metadata.layout
    if isequal(layout, :curvature_branch_regions)
        dimensions = observation.metadata.dimensions
        length(dimensions) == 2 && dimensions[1] == 2 * dimensions[2] ||
            throw(WFSPreparationError(:estimation, :shape,
                "packed Curvature regions require a 2N-by-N frame"))
        side = dimensions[2]
        side % sensor.params.pupil_samples == 0 || throw(
            WFSPreparationError(:estimation, :shape,
                "packed Curvature regions must reduce to the pupil grid"))
        _require_wfs_storage_domain(:estimation, observation.metadata,
            sensor.estimator.state.reduced_plus,
            "packed Curvature observation")
        plus = @view observation.storage[1:side, :]
        minus = @view observation.storage[side+1:2*side, :]
        reduction = div(side, sensor.params.pupil_samples)
        mapping = CurvatureImagePairMapping(plus, minus, reduction,
            reduction)
    elseif isequal(layout, :curvature_branch_channels)
        dimensions = observation.metadata.dimensions
        dimensions == (2, sensor.params.pupil_samples^2) || throw(
            WFSPreparationError(:estimation, :shape,
                "packed Curvature channels require a 2-by-N² observation"))
        _require_wfs_storage_domain(:estimation, observation.metadata,
            sensor.estimator.state.reduced_plus,
            "packed Curvature observation")
        mapping = CurvatureChannelPairMapping(observation.storage)
    else
        throw(WFSPreparationError(:estimation, :detector_mapping,
            "Curvature estimator requires branch regions, channels, or a separate observation pair"))
    end
    scales = _curvature_rate_scales(sensor, branch_rate_scales)
    binding = _curvature_calibration_binding(sensor)
    return PreparedCurvatureEstimator(sensor, observation, measurement,
        AcquiredObservationPath(), binding, mapping, scales)
end

function _reduce_curvature_observations!(::ScalarCPUStyle,
    sensor::CurvatureWFS, mapping::CurvatureImagePairMapping, scales)
    state = sensor.estimator.state
    T = eltype(state.reduced_plus)
    n_sub = sensor.params.pupil_samples
    @inbounds for j in 1:n_sub, i in 1:n_sub
        plus = zero(T)
        minus = zero(T)
        for jj in 1:mapping.plus_reduction,
                ii in 1:mapping.plus_reduction
            plus += T(mapping.plus[(i - 1) * mapping.plus_reduction + ii,
                (j - 1) * mapping.plus_reduction + jj])
        end
        for jj in 1:mapping.minus_reduction,
                ii in 1:mapping.minus_reduction
            minus += T(mapping.minus[(i - 1) * mapping.minus_reduction + ii,
                (j - 1) * mapping.minus_reduction + jj])
        end
        state.reduced_plus[i, j] = plus * scales[1]
        state.reduced_minus[i, j] = minus * scales[2]
    end
    return nothing
end

function _reduce_curvature_observations!(style::AcceleratorStyle,
    sensor::CurvatureWFS, mapping::CurvatureImagePairMapping, scales)
    state = sensor.estimator.state
    launch_kernel!(style, curvature_reduce_observation_pair_kernel!,
        state.reduced_plus, state.reduced_minus, mapping.plus, mapping.minus,
        mapping.plus_reduction, mapping.minus_reduction, scales[1],
        scales[2], sensor.params.pupil_samples;
        ndrange=size(state.reduced_plus))
    return nothing
end

function _reduce_curvature_observations!(::ScalarCPUStyle,
    sensor::CurvatureWFS, mapping::CurvatureChannelPairMapping, scales)
    state = sensor.estimator.state
    T = eltype(state.reduced_plus)
    n_sub = sensor.params.pupil_samples
    @inbounds for i in 1:n_sub, j in 1:n_sub
        index = (i - 1) * n_sub + j
        state.reduced_plus[i, j] = T(mapping.channels[1, index]) *
            scales[1]
        state.reduced_minus[i, j] = T(mapping.channels[2, index]) *
            scales[2]
    end
    return nothing
end

function _reduce_curvature_observations!(style::AcceleratorStyle,
    sensor::CurvatureWFS, mapping::CurvatureChannelPairMapping, scales)
    state = sensor.estimator.state
    launch_kernel!(style, curvature_unpack_channel_pair_kernel!,
        state.reduced_plus, state.reduced_minus, mapping.channels,
        scales[1], scales[2], sensor.params.pupil_samples;
        ndrange=size(state.reduced_plus))
    return nothing
end

function estimate_wfs_measurement!(measurement::WFSMeasurement, input,
    plan::PreparedCurvatureEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Curvature estimator storage does not match its plan"))
    sensor = plan.sensor
    _require_curvature_calibration(sensor, plan.calibration_binding)
    style = execution_style(sensor.estimator.state.reduced_plus)
    _reduce_curvature_observations!(style, sensor, plan.mapping,
        plan.branch_rate_scales)
    _estimate_curvature_from_reduced!(style, sensor)
    copyto!(measurement.storage, sensor.estimator.state.slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedCurvatureEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Curvature estimator storage does not match its plan"))
    return nothing
end

@inline _estimate_curvature_from_reduced!(::ScalarCPUStyle,
    sensor::CurvatureWFS) = curvature_signal_from_planes!(sensor)

@inline _estimate_curvature_from_reduced!(style::AcceleratorStyle,
    sensor::CurvatureWFS) = curvature_signal_from_planes!(style, sensor)

function set_curvature_calibration!(sensor::CurvatureWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0))
    state = sensor.estimator.state
    size(reference) == size(state.reference_signal_2d) || throw(
        InvalidConfiguration(
            "Curvature calibration reference has the wrong dimensions"))
    eltype(reference) <: Real || throw(InvalidConfiguration(
        "Curvature calibration reference must contain real values"))
    typeof(backend(reference)) === typeof(backend(state.reference_signal_2d)) ||
        throw(InvalidConfiguration(
            "Curvature calibration reference backend differs from the estimator"))
    plane_device(reference) == plane_device(state.reference_signal_2d) ||
        throw(InvalidConfiguration(
            "Curvature calibration reference occupies another device"))
    all(isfinite, host_array(reference)) || throw(InvalidConfiguration(
        "Curvature calibration reference must contain finite values"))
    T = eltype(state.reference_signal_2d)
    wavelength_value = T(wavelength_m)
    isfinite(wavelength_value) && wavelength_value > zero(T) || throw(
        InvalidConfiguration(
            "Curvature calibration wavelength must be finite and positive"))
    copyto!(state.reference_signal_2d, reference)
    fill!(state.signal_2d, zero(T))
    fill!(state.slopes, zero(T))
    state.calibrated = true
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibration_revision += UInt(1)
    return sensor
end
