const _DETECTOR_STATE_TARGET_VALIDATION_FIELDS = (
    :accum_buffer,
    :latent_buffer,
    :thermal_state,
    :integrated_time,
    :readout_ready,
)

const _DETECTOR_WORKSPACE_TARGET_VALIDATION_FIELDS = (
    :presampling_buffer,
    :presampling_scratch,
    :response_buffer,
    :bin_buffer,
    :temporal_buffer,
    :noise_buffer,
    :noise_buffer_host,
    :batched_buffer_host,
    :output_buffer_host,
    :readout,
)

const _DETECTOR_PRODUCTS_TARGET_VALIDATION_FIELDS = (
    :frame,
    :output_buffer,
    :readout,
)

const _DETECTOR_ACQUISITION_PLAN_TARGET_VALIDATION_FIELDS = (
    :detector_params,
    :input_metadata,
    :input_shape,
    :frame_shape,
    :output_shape,
    :rate_scale,
    :quantum_efficiency,
)

const _PREPARED_DETECTOR_ACQUISITION_TARGET_VALIDATION_FIELDS = (
    :detector,
    :input,
    :plan,
    :state,
    :workspace,
    :products,
    :state_binding,
    :workspace_binding,
    :product_binding,
)

const _SKIPPER_READOUT_TARGET_VALIDATION_FIELDS = (
    :mean_frame,
    :sample_count,
)

const _SAMPLED_READOUT_TARGET_VALIDATION_FIELDS = (
    :reference_frame,
    :signal_frame,
    :read_cube,
)

const _MULTI_READ_TARGET_VALIDATION_FIELDS = (
    :reference_frame,
    :signal_frame,
    :combined_frame,
    :reference_cube,
    :signal_cube,
    :read_cube,
    :read_offsets_s,
)

const _RAMP_READOUT_TARGET_VALIDATION_FIELDS = (
    :slope_frame,
    :intercept_frame,
    :integrated_frame,
    :read_cube,
    :read_offsets_s,
    :acquisition_kind,
)

const _COUNTING_DETECTOR_STATE_TARGET_VALIDATION_FIELDS = (
    :thermal_state,
)

const _COUNTING_DETECTOR_WORKSPACE_TARGET_VALIDATION_FIELDS = (
    :noise_buffer,
    :host_buffer,
    :output_buffer_host,
)

const _COUNTING_DETECTOR_PRODUCTS_TARGET_VALIDATION_FIELDS = (
    :counts,
    :output_buffer,
)

@noinline function _throw_wrong_detector_target(
    target::AbstractComputeDevice,
    label::AbstractString,
    actual::AbstractComputeDevice,
)
    _throw_compute_device_error(
        :validate_detector_acquisition_target,
        :wrong_device,
        target,
        "$label occupies $(actual)",
    )
end

@inline function _require_exact_detector_array_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    actual = compute_device(storage)
    actual == target || _throw_wrong_detector_target(target, label, actual)
    return storage
end

@inline _require_exact_optional_detector_array_target(
    ::Nothing, ::AbstractComputeDevice, ::AbstractString) = nothing

@inline function _require_exact_optional_detector_array_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    return _require_exact_detector_array_target(storage, target, label)
end

function _require_exact_optional_detector_array_target(
    storage,
    ::AbstractComputeDevice,
    label::AbstractString,
)
    throw(InvalidConfiguration(
        "$label must be target-local array storage or nothing; got " *
        "$(typeof(storage))"))
end

@inline function _require_exact_detector_response_target(
    ::NullFrameResponse, ::AbstractComputeDevice)
    return nothing
end

@inline function _require_exact_detector_response_target(
    response::GaussianPixelResponse, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        response.kernel, target, "detector Gaussian response kernel")
    return response
end

@inline function _require_exact_detector_response_target(
    response::SampledFrameResponse, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        response.kernel, target, "detector sampled response kernel")
    return response
end

@inline function _require_exact_detector_response_target(
    response::RectangularPixelAperture, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        response.kernel_x, target, "detector pixel-aperture x kernel")
    _require_exact_detector_array_target(
        response.kernel_y, target, "detector pixel-aperture y kernel")
    return response
end

function _require_exact_detector_response_target(
    response::AbstractFrameResponse, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for detector response " *
        "$(typeof(response))"))
end

@inline function _require_exact_detector_coupling_target(
    ::NullChargeCoupling, ::AbstractComputeDevice)
    return nothing
end

@inline function _require_exact_detector_coupling_target(
    coupling::InterpixelCapacitance, target::AbstractComputeDevice)
    _require_exact_detector_response_target(coupling.response, target)
    return coupling
end

function _require_exact_detector_coupling_target(
    coupling::AbstractChargeCouplingModel, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for detector charge coupling " *
        "$(typeof(coupling))"))
end

@inline function _require_exact_detector_defect_target(
    ::NullDetectorDefectModel, ::AbstractComputeDevice)
    return nothing
end

@inline function _require_exact_detector_defect_target(
    model::PixelResponseNonuniformity, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        model.gain_map, target, "detector PRNU gain map")
    return model
end

@inline function _require_exact_detector_defect_target(
    model::DarkSignalNonuniformity, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        model.dark_map, target, "detector DSNU dark map")
    return model
end

@inline function _require_exact_detector_defect_target(
    model::BadPixelMask, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        model.mask, target, "detector bad-pixel mask")
    return model
end

@inline _require_exact_detector_defect_stages_target(
    ::Tuple{}, ::AbstractComputeDevice) = nothing

@inline function _require_exact_detector_defect_stages_target(
    stages::Tuple, target::AbstractComputeDevice)
    _require_exact_detector_defect_target(first(stages), target)
    _require_exact_detector_defect_stages_target(Base.tail(stages), target)
    return stages
end

@inline function _require_exact_detector_defect_target(
    model::CompositeDetectorDefectModel, target::AbstractComputeDevice)
    _require_exact_detector_defect_stages_target(model.stages, target)
    return model
end

function _require_exact_detector_defect_target(
    model::AbstractDetectorDefectModel, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for detector defect model " *
        "$(typeof(model))"))
end

@inline _require_exact_cmos_read_noise_target(
    ::NullCMOSReadNoise, ::AbstractComputeDevice) = nothing

@inline function _require_exact_cmos_read_noise_target(
    model::CMOSReadNoiseMap, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        model.sigma, target, "detector CMOS read-noise map")
    return model
end

function _require_exact_cmos_read_noise_target(
    model::AbstractCMOSReadNoiseModel, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for CMOS read-noise model " *
        "$(typeof(model))"))
end

@inline _require_exact_cmos_output_target(
    ::NullCMOSOutputModel, ::AbstractComputeDevice) = nothing

@inline function _require_exact_cmos_output_target(
    model::StaticCMOSOutputPattern, target::AbstractComputeDevice)
    _require_exact_detector_array_target(
        model.gains, target, "detector CMOS output gains")
    _require_exact_detector_array_target(
        model.offsets, target, "detector CMOS output offsets")
    return model
end

function _require_exact_cmos_output_target(
    model::AbstractCMOSOutputModel, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for CMOS output model " *
        "$(typeof(model))"))
end

@inline _require_exact_detector_sensor_target(
    ::CCDSensor, ::AbstractComputeDevice) = nothing
@inline _require_exact_detector_sensor_target(
    ::EMCCDSensor, ::AbstractComputeDevice) = nothing
@inline _require_exact_detector_sensor_target(
    ::InGaAsSensor, ::AbstractComputeDevice) = nothing
@inline _require_exact_detector_sensor_target(
    ::HgCdTeSensor, ::AbstractComputeDevice) = nothing
@inline _require_exact_detector_sensor_target(
    ::HgCdTeAvalancheArraySensor, ::AbstractComputeDevice) = nothing

@inline function _require_exact_detector_sensor_target(
    sensor::CMOSSensor, target::AbstractComputeDevice)
    _require_exact_cmos_read_noise_target(sensor.readout_noise_model, target)
    _require_exact_cmos_output_target(sensor.output_model, target)
    return sensor
end

function _require_exact_detector_sensor_target(
    sensor::AbstractFrameSensor, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for frame-detector sensor " *
        "$(typeof(sensor))"))
end

@inline _require_exact_detector_background_target(
    ::NoBackground, ::AbstractComputeDevice, ::AbstractString) = nothing
@inline _require_exact_detector_background_target(
    ::ScalarBackground, ::AbstractComputeDevice, ::AbstractString) = nothing

@inline function _require_exact_detector_background_target(
    background::BackgroundFrame,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_detector_array_target(background.map, target, label)
    return background
end

function _require_exact_detector_background_target(
    background::BackgroundModel,
    ::AbstractComputeDevice,
    label::AbstractString,
)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for $label model " *
        "$(typeof(background))"))
end

@inline function _require_detector_target_validation_layout(
    detector::Detector)
    fieldnames(typeof(detector.state)) ==
        _DETECTOR_STATE_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "DetectorState fields changed without updating exact-target validation"))
    fieldnames(typeof(detector.workspace)) ==
        _DETECTOR_WORKSPACE_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "DetectorWorkspace fields changed without updating exact-target validation"))
    fieldnames(typeof(detector.products)) ==
        _DETECTOR_PRODUCTS_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "DetectorProducts fields changed without updating exact-target validation"))
    return detector
end

@inline function _require_detector_readout_target_validation_layout(
    products::FrameReadoutProducts,
    expected_fields::Tuple,
)
    fieldnames(typeof(products)) == expected_fields || throw(
        InvalidConfiguration(
            "$(nameof(typeof(products))) fields changed without updating " *
            "exact-target validation"))
    return products
end

function _require_exact_detector_readout_products_target(
    products::NoFrameReadoutProducts, ::AbstractComputeDevice)
    _require_detector_readout_target_validation_layout(products, ())
    return products
end

function _require_exact_detector_readout_products_target(
    products::SkipperReadoutProducts, target::AbstractComputeDevice)
    _require_detector_readout_target_validation_layout(
        products, _SKIPPER_READOUT_TARGET_VALIDATION_FIELDS)
    _require_exact_detector_array_target(
        products.mean_frame, target, "detector Skipper mean frame")
    # `sample_count` is host scalar configuration.
    return products
end

function _require_exact_detector_readout_products_target(
    products::SampledFrameReadoutProducts, target::AbstractComputeDevice)
    _require_detector_readout_target_validation_layout(
        products, _SAMPLED_READOUT_TARGET_VALIDATION_FIELDS)
    _require_exact_optional_detector_array_target(products.reference_frame,
        target, "detector sampled reference frame")
    _require_exact_detector_array_target(
        products.signal_frame, target, "detector sampled signal frame")
    _require_exact_optional_detector_array_target(
        products.read_cube, target, "detector sampled read cube")
    return products
end

function _require_exact_detector_readout_products_target(
    products::MultiReadFrameReadoutProducts, target::AbstractComputeDevice)
    _require_detector_readout_target_validation_layout(
        products, _MULTI_READ_TARGET_VALIDATION_FIELDS)
    _require_exact_optional_detector_array_target(products.reference_frame,
        target, "detector multi-read reference frame")
    _require_exact_detector_array_target(
        products.signal_frame, target, "detector multi-read signal frame")
    _require_exact_detector_array_target(products.combined_frame,
        target, "detector multi-read combined frame")
    _require_exact_optional_detector_array_target(products.reference_cube,
        target, "detector multi-read reference cube")
    _require_exact_optional_detector_array_target(products.signal_cube,
        target, "detector multi-read signal cube")
    _require_exact_optional_detector_array_target(products.read_cube,
        target, "detector multi-read read cube")
    # `read_offsets_s` is a host read-schedule vector, not a data-plane array.
    return products
end

function _require_exact_detector_readout_products_target(
    products::UpTheRampReadoutProducts, target::AbstractComputeDevice)
    _require_detector_readout_target_validation_layout(
        products, _RAMP_READOUT_TARGET_VALIDATION_FIELDS)
    _require_exact_detector_array_target(
        products.slope_frame, target, "detector ramp slope frame")
    _require_exact_detector_array_target(
        products.intercept_frame, target, "detector ramp intercept frame")
    _require_exact_detector_array_target(products.integrated_frame,
        target, "detector ramp integrated frame")
    _require_exact_detector_array_target(
        products.read_cube, target, "detector ramp read cube")
    # `read_offsets_s` and `acquisition_kind` are host schedule configuration.
    return products
end

function _require_exact_detector_readout_products_target(
    products::FrameReadoutProducts, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for detector readout products " *
        "$(typeof(products))"))
end

@inline function _require_exact_detector_readout_workspace_target(
    workspace::NoFrameReadoutWorkspace, ::AbstractComputeDevice)
    return workspace
end

function _require_exact_detector_readout_workspace_target(
    workspace::SkipperReadoutWorkspace, target::AbstractComputeDevice)
    _require_exact_detector_array_target(workspace.baseline_frame, target,
        "detector Skipper baseline workspace")
    _require_exact_detector_array_target(workspace.sample_sum, target,
        "detector Skipper sample-sum workspace")
    return workspace
end

function _require_exact_detector_readout_workspace_target(
    workspace::MultiReadFrameReadoutWorkspace,
    target::AbstractComputeDevice)
    _require_exact_detector_array_target(workspace.reference_average, target,
        "detector multi-read reference-average workspace")
    _require_exact_detector_array_target(workspace.signal_average, target,
        "detector multi-read signal-average workspace")
    _require_exact_optional_detector_array_target(workspace.reference_cube,
        target, "detector multi-read reference-cube workspace")
    _require_exact_detector_array_target(workspace.signal_cube, target,
        "detector multi-read signal-cube workspace")
    return workspace
end

function _require_exact_detector_readout_workspace_target(
    workspace::UpTheRampReadoutWorkspace, target::AbstractComputeDevice)
    _require_exact_detector_array_target(workspace.slope, target,
        "detector ramp slope workspace")
    _require_exact_detector_array_target(workspace.intercept, target,
        "detector ramp intercept workspace")
    _require_exact_detector_array_target(workspace.integrated, target,
        "detector ramp integrated workspace")
    _require_exact_detector_array_target(workspace.cube, target,
        "detector ramp cube workspace")
    return workspace
end

function _require_exact_detector_readout_workspace_target(
    workspace::FrameReadoutWorkspace, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for detector readout workspace " *
        "$(typeof(workspace))"))
end

"""
    _require_exact_detector_target(detector, target)

Internal cold-path validation that every detector-owned numerical data-plane
array occupies `target`. The detector's sampled-QE vectors are host
configuration. `noise_buffer_host`, `batched_buffer_host`, and
`output_buffer_host` are deliberate host staging buffers and are excluded.
"""
function _require_exact_detector_target(
    detector::Detector, target::AbstractComputeDevice)
    _require_detector_target_validation_layout(detector)
    params = detector.params
    state = detector.state
    workspace = detector.workspace
    products = detector.products

    _require_exact_detector_sensor_target(params.sensor, target)
    _require_exact_detector_response_target(params.response_model, target)
    _require_exact_detector_coupling_target(
        params.charge_coupling_model, target)
    _require_exact_detector_defect_target(params.defect_model, target)
    _require_exact_detector_background_target(detector.background_flux,
        target, "detector background-flux frame")
    _require_exact_detector_background_target(detector.background_map,
        target, "detector background-subtraction frame")

    _require_exact_detector_array_target(
        products.frame, target, "detector frame")
    _require_exact_detector_array_target(workspace.presampling_buffer,
        target, "detector presampling buffer")
    _require_exact_detector_array_target(workspace.presampling_scratch,
        target, "detector presampling scratch")
    _require_exact_detector_array_target(
        workspace.response_buffer, target, "detector response buffer")
    _require_exact_detector_array_target(
        workspace.bin_buffer, target, "detector bin buffer")
    _require_exact_detector_array_target(
        workspace.temporal_buffer, target, "detector temporal buffer")
    _require_exact_detector_array_target(
        workspace.noise_buffer, target, "detector noise buffer")
    _require_exact_detector_array_target(
        state.accum_buffer, target, "detector accumulation buffer")
    _require_exact_detector_array_target(
        state.latent_buffer, target, "detector latent buffer")
    _require_exact_optional_detector_array_target(
        products.output_buffer, target, "detector output buffer")
    _require_exact_detector_readout_products_target(
        products.readout, target)
    _require_exact_detector_readout_workspace_target(
        workspace.readout, target)

    # Deliberate host staging: workspace.noise_buffer_host,
    # workspace.batched_buffer_host, and workspace.output_buffer_host. Thermal
    # state, exposure duration, readiness, and all remaining params are host
    # configuration or scalars.
    return detector
end

"""
    _require_exact_counting_detector_target(detector, target)

Internal cold-path validation for a maintained accumulated-count detector.
The count, noise, and optional converted-output arrays must occupy `target`.
`host_buffer` and `output_buffer_host` are deliberate RNG/output staging
mirrors; thermal state and detector parameters are host scalar configuration.
"""
function _require_exact_counting_detector_target(
    detector::Union{SPADArrayDetector,MKIDArrayDetector},
    target::AbstractComputeDevice,
)
    state = detector.state
    workspace = detector.workspace
    products = detector.products
    fieldnames(typeof(state)) ==
        _COUNTING_DETECTOR_STATE_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "counting-detector state fields changed without updating " *
            "exact-target validation"))
    fieldnames(typeof(workspace)) ==
        _COUNTING_DETECTOR_WORKSPACE_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "counting-detector workspace fields changed without updating " *
            "exact-target validation"))
    fieldnames(typeof(products)) ==
        _COUNTING_DETECTOR_PRODUCTS_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "counting-detector product fields changed without updating " *
            "exact-target validation"))
    _require_exact_detector_array_target(
        products.counts, target, "counting detector accumulation array")
    _require_exact_detector_array_target(
        workspace.noise_buffer, target, "counting detector noise buffer")
    _require_exact_optional_detector_array_target(
        products.output_buffer, target, "counting detector output buffer")
    return detector
end

function _require_exact_counting_detector_target(
    detector::AbstractCountingDetector,
    ::AbstractComputeDevice,
)
    throw(InvalidConfiguration(
        "no exact-target validator is defined for counting detector " *
        "$(typeof(detector))"))
end

"""
    _require_exact_detector_acquisition_target(prepared, target)

Internal fail-closed exact-target validation for a prepared frame-detector
acquisition owner.
"""
function _require_exact_detector_acquisition_target(
    prepared::PreparedDetectorAcquisition,
    target::AbstractComputeDevice,
)
    fieldnames(typeof(prepared)) ==
        _PREPARED_DETECTOR_ACQUISITION_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "PreparedDetectorAcquisition fields changed without updating " *
            "exact-target validation"))
    plan = prepared.plan
    fieldnames(typeof(plan)) ==
        _DETECTOR_ACQUISITION_PLAN_TARGET_VALIDATION_FIELDS || throw(
        InvalidConfiguration(
            "DetectorAcquisitionPlan fields changed without updating " *
            "exact-target validation"))
    _require_prepared_detector_binding(prepared)
    detector = prepared.detector
    _require_exact_detector_target(detector, target)

    metadata_device = plan.input_metadata.device
    metadata_device == target || _throw_wrong_detector_target(
        target, "detector acquisition input metadata", metadata_device)
    _require_exact_detector_array_target(prepared.input.values, target,
        "detector acquisition input storage")
    validate_plane_storage(plan.input_metadata, prepared.input.values;
        label="detector acquisition input")

    # Radiometric scale and quantum efficiency are host scalars. Exact state,
    # workspace, product, input, and detector bindings were checked above.
    return prepared
end

@inline function _require_exact_detector_acquisition_target(
    detector::Detector, prepared::PreparedDetectorAcquisition,
    target::AbstractComputeDevice)
    detector === prepared.detector || throw(InvalidConfiguration(
        "detector does not match its prepared acquisition owner"))
    return _require_exact_detector_acquisition_target(prepared, target)
end
