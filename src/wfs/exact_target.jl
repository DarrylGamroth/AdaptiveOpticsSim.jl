#
# Exact-target validation for prepared WFS stage plans.
#
# These cold-path validators enumerate WFS-owned numerical data-plane storage.
# Host configuration, scalar state, fixed registries, and documented host
# mirrors/staging buffers are deliberately excluded.
#

@noinline function _throw_wrong_wfs_target(
    stage::Symbol,
    target::AbstractComputeDevice,
    label::AbstractString,
    actual::AbstractComputeDevice,
)
    throw(WFSPreparationError(
        stage,
        :device,
        "$label occupies $actual instead of the requested exact target $target",
    ))
end

@inline function _require_exact_wfs_storage_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    stage::Symbol,
    label::AbstractString,
)
    actual = _wfs_storage_device(storage)
    actual == target ||
        _throw_wrong_wfs_target(stage, target, label, actual)
    return storage
end


@inline function _require_exact_wfs_storage_target(
    storage::Base.RefValue,
    ::AbstractComputeDevice,
    ::Symbol,
    ::AbstractString,
)
    # WFS stage contracts define `Ref{T}` as a host-resident scalar product,
    # not as an accelerator data-plane array.
    return storage
end

@inline function _require_exact_wfs_metadata_target(
    metadata,
    target::AbstractComputeDevice,
    stage::Symbol,
    label::AbstractString,
)
    actual = metadata.device
    actual == target ||
        _throw_wrong_wfs_target(stage, target, "$label metadata", actual)
    return metadata
end


@inline function _require_exact_wfs_product_metadata_target(
    metadata,
    ::AbstractArray,
    target::AbstractComputeDevice,
    stage::Symbol,
    label::AbstractString,
)
    return _require_exact_wfs_metadata_target(
        metadata, target, stage, label)
end


@inline function _require_exact_wfs_product_metadata_target(
    metadata,
    ::Base.RefValue,
    ::AbstractComputeDevice,
    ::Symbol,
    ::AbstractString,
)
    # The ordinary WFS validator has already bound this metadata to the
    # canonical host scalar domain.
    return metadata
end

@inline _require_exact_wfs_array_targets(
    ::Tuple{},
    ::Tuple{},
    ::AbstractComputeDevice,
    ::Symbol,
) = nothing

@inline function _require_exact_wfs_array_targets(
    storages::Tuple,
    labels::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_storage_target(
        first(storages), target, stage, first(labels))
    return _require_exact_wfs_array_targets(
        Base.tail(storages), Base.tail(labels), target, stage)
end

function _require_exact_wfs_input_target(
    input::PupilFunction,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    validate_wfs_optical_input(input)
    _require_exact_wfs_metadata_target(
        input.metadata, target, stage, "WFS pupil input")
    _require_exact_wfs_array_targets(
        (input.support, input.amplitude, input.opd),
        ("WFS pupil support", "WFS pupil amplitude", "WFS pupil OPD"),
        target,
        stage,
    )
    size(input.support) == input.metadata.dimensions ||
        throw(WFSPreparationError(stage, :shape,
            "WFS pupil support dimensions do not match its metadata"))
    typeof(backend(input.support)) === typeof(input.metadata.backend) ||
        throw(WFSPreparationError(stage, :backend,
            "WFS pupil support backend does not match its metadata"))
    return input
end

function _require_exact_wfs_input_target(
    input::ElectricField,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    validate_wfs_optical_input(input)
    _require_exact_wfs_metadata_target(
        input.metadata, target, stage, "WFS electric-field input")
    _require_exact_wfs_storage_target(
        input.values, target, stage, "WFS electric-field storage")
    return input
end

@inline _require_exact_wfs_input_targets(
    ::Tuple{},
    ::AbstractComputeDevice,
    ::Symbol,
) = nothing

@inline function _require_exact_wfs_input_targets(
    inputs::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_input_target(first(inputs), target, stage)
    return _require_exact_wfs_input_targets(
        Base.tail(inputs), target, stage)
end

function _require_exact_wfs_input_target(
    inputs::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_input_targets(inputs, target, stage)
    return inputs
end

function _require_exact_wfs_input_target(
    inputs::AbstractVector,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    @inbounds for input in inputs
        _require_exact_wfs_input_target(input, target, stage)
    end
    return inputs
end

function _require_exact_wfs_input_target(
    input,
    ::AbstractComputeDevice,
    stage::Symbol,
)
    throw(WFSPreparationError(stage, :unsupported_target_validation,
        "no exact-target validator is defined for WFS input $(typeof(input))"))
end

function _require_exact_wfs_product_target(
    product::IntensityMap,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    validate_wfs_optical_products(product)
    _require_exact_wfs_metadata_target(
        product.metadata, target, stage, "WFS optical product")
    _require_exact_wfs_storage_target(
        product.values, target, stage, "WFS optical-product storage")
    return product
end

@inline _require_exact_wfs_product_targets(
    ::Tuple{},
    ::AbstractComputeDevice,
    ::Symbol,
) = nothing

@inline function _require_exact_wfs_product_targets(
    products::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_product_target(first(products), target, stage)
    return _require_exact_wfs_product_targets(
        Base.tail(products), target, stage)
end

function _require_exact_wfs_product_target(
    products::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_product_targets(products, target, stage)
    return products
end

function _require_exact_wfs_product_target(
    products::AbstractVector,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    @inbounds for product in products
        _require_exact_wfs_product_target(product, target, stage)
    end
    return products
end

function _require_exact_wfs_product_target(
    bundle::OpticalProductBundle,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_product_target(bundle.products, target, stage)
    return bundle
end

function _require_exact_wfs_observation_target(
    observation::WFSObservation,
    target::AbstractComputeDevice,
    stage::Symbol=:acquisition,
)
    validate_wfs_observation(observation)
    _require_exact_wfs_product_metadata_target(
        observation.metadata, observation.storage, target, stage,
        "WFS observation")
    _require_exact_wfs_storage_target(
        observation.storage, target, stage,
        "WFS observation storage")
    return observation
end

@inline _require_exact_wfs_observation_targets(
    ::Tuple{},
    ::AbstractComputeDevice,
    ::Symbol,
) = nothing

@inline function _require_exact_wfs_observation_targets(
    observations::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_observation_target(
        first(observations), target, stage)
    return _require_exact_wfs_observation_targets(
        Base.tail(observations), target, stage)
end

function _require_exact_wfs_observation_target(
    observations::Tuple,
    target::AbstractComputeDevice,
    stage::Symbol=:acquisition,
)
    _require_exact_wfs_observation_targets(observations, target, stage)
    return observations
end

@inline _require_exact_wfs_estimator_input_target(
    input::Union{PupilFunction,ElectricField},
    target::AbstractComputeDevice,
) = _require_exact_wfs_input_target(input, target, :estimation)

@inline _require_exact_wfs_estimator_input_target(
    input::WFSObservation,
    target::AbstractComputeDevice,
) = _require_exact_wfs_observation_target(input, target, :estimation)

@inline _require_exact_wfs_estimator_input_target(
    input::Tuple,
    target::AbstractComputeDevice,
) = _require_exact_wfs_observation_target(input, target, :estimation)

function _require_exact_wfs_measurement_target(
    measurement::WFSMeasurement,
    target::AbstractComputeDevice,
)
    validate_wfs_measurement(measurement)
    _require_exact_wfs_product_metadata_target(
        measurement.metadata, measurement.storage, target, :estimation,
        "WFS measurement")
    _require_exact_wfs_storage_target(
        measurement.storage, target, :estimation,
        "WFS measurement storage")
    return measurement
end

function _require_exact_sh_propagation_target(
    propagation::PreparedMicrolensPropagation,
    target::AbstractComputeDevice,
)
    workspace = microlens_propagation_workspace(propagation)
    _require_exact_wfs_array_targets(
        (
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
            workspace.elongation_kernel,
            workspace.lgs_kernel_fft,
            workspace.fft_asterism_stack,
            workspace.amp_scales,
            workspace.opd_to_cycles,
        ),
        (
            "Shack-Hartmann field",
            "Shack-Hartmann phasor",
            "Shack-Hartmann FFT buffer",
            "Shack-Hartmann FFT stack",
            "Shack-Hartmann intensity",
            "Shack-Hartmann intensity stack",
            "Shack-Hartmann temporary intensity stack",
            "Shack-Hartmann temporary image",
            "Shack-Hartmann bin buffer",
            "Shack-Hartmann spot",
            "Shack-Hartmann sampled-spot cube",
            "Shack-Hartmann accumulated-spot cube",
            "Shack-Hartmann elongation kernel",
            "Shack-Hartmann LGS Fourier kernel",
            "Shack-Hartmann asterism FFT stack",
            "Shack-Hartmann amplitude scales",
            "Shack-Hartmann OPD-to-cycle scales",
        ),
        target,
        :wfs_optics,
    )
    return propagation
end

function _require_exact_sh_layout_target(
    layout::SubapertureLayout,
    target::AbstractComputeDevice,
    stage::Symbol,
)
    _require_exact_wfs_storage_target(
        layout.valid_mask, target, stage,
        "Shack-Hartmann valid-subaperture mask")
    return layout
end

function _require_exact_sh_front_end_target(
    optics::ShackHartmannOptics,
    target::AbstractComputeDevice,
)
    _require_exact_sh_propagation_target(optics.propagation, target)
    _require_exact_sh_layout_target(
        optics.front_end.layout, target, :wfs_optics)
    return optics
end

function _require_exact_wfs_target(
    plan::PreparedShackHartmannOptics,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :wfs_optics)
    _require_exact_wfs_product_target(
        plan.output, target, :wfs_optics)
    _require_exact_sh_front_end_target(plan.optics, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedShackHartmannOpticsBundle,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.components
        _require_exact_wfs_target(component, target)
    end
    return plan
end

@inline _require_exact_prepared_four_pupil_lgs_target(
    ::NoPreparedFourPupilLGS,
    ::AbstractComputeDevice,
) = nothing

function _require_exact_prepared_four_pupil_lgs_target(
    model::PreparedFourPupilElongation,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(
        model.kernel, target, :wfs_optics,
        "four-pupil LGS elongation kernel")
    return nothing
end

function _require_exact_prepared_four_pupil_lgs_target(
    model::PreparedFourPupilSodiumLayerProfile,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(
        model.kernel_fft, target, :wfs_optics,
        "four-pupil sodium-layer-profile Fourier kernel")
    return nothing
end

function _require_exact_focal_plane_modulation_target(
    modulation::PreparedFocalPlaneModulation,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_wfs_storage_target(
        modulation.phases, target, :wfs_optics, "$label phases")
    return modulation
end

@inline function _require_exact_pyramid_modulation_batch_target(
    ::NoPyramidModulationBatchWorkspace,
    ::AbstractComputeDevice,
)
    return nothing
end

function _require_exact_pyramid_modulation_batch_target(
    batch::PyramidModulationBatchWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            batch.field_stack,
            batch.operating_weights,
            batch.calibration_weights,
        ),
        (
            "Pyramid modulation field stack",
            "Pyramid operating modulation weights",
            "Pyramid calibration modulation weights",
        ),
        target,
        :wfs_optics,
    )
    return nothing
end


function _require_exact_pyramid_modulation_batch_target(
    batch::PyramidShiftedMaskModulationWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            batch.field_stack,
            batch.shifted_masks,
            batch.operating_weights,
            batch.axis_1_shifts_rad,
            batch.axis_2_shifts_rad,
        ),
        (
            "Pyramid shifted-mask field stack",
            "Pyramid shifted focal masks",
            "Pyramid shifted-mask operating weights",
            "Pyramid shifted-mask axis-1 coordinate shifts",
            "Pyramid shifted-mask axis-2 coordinate shifts",
        ),
        target,
        :wfs_optics,
    )
    return nothing
end

function _require_exact_pyramid_front_end_target(
    front_end::PyramidOpticalFrontEnd,
    target::AbstractComputeDevice,
)
    propagation = pyramid_propagation_workspace(front_end)
    _require_exact_wfs_array_targets(
        (
            propagation.field,
            propagation.focal_field,
            propagation.pupil_field,
            propagation.pyramid_mask,
            propagation.phasor,
            propagation.intensity,
            propagation.temp,
            propagation.scratch,
            propagation.asterism_stack,
            propagation.elongation_kernel,
            propagation.lgs_kernel_fft,
        ),
        (
            "Pyramid field",
            "Pyramid focal field",
            "Pyramid pupil field",
            "Pyramid phase mask",
            "Pyramid phasor",
            "Pyramid intensity",
            "Pyramid temporary image",
            "Pyramid scratch image",
            "Pyramid asterism stack",
            "Pyramid elongation kernel",
            "Pyramid LGS Fourier kernel",
        ),
        target,
        :wfs_optics,
    )
    _require_exact_focal_plane_modulation_target(
        front_end.modulation, target, "Pyramid operating modulation")
    _require_exact_focal_plane_modulation_target(
        front_end.calibration_modulation, target,
        "Pyramid calibration modulation")
    _require_exact_pyramid_modulation_batch_target(
        propagation.modulation_batch, target)
    return front_end
end

function _require_exact_wfs_target(
    plan::PreparedPyramidOptics,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :wfs_optics)
    _require_exact_wfs_product_target(
        plan.output, target, :wfs_optics)
    _require_exact_pyramid_front_end_target(plan.front_end, target)
    _require_exact_prepared_four_pupil_lgs_target(
        plan.plan.lgs_model, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedPyramidOpticsBundle,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.components
        _require_exact_wfs_target(component, target)
    end
    return plan
end

function _require_exact_bi_o_edge_front_end_target(
    front_end::BiOEdgeOpticalFrontEnd,
    target::AbstractComputeDevice,
)
    propagation = bi_o_edge_propagation_workspace(front_end)
    _require_exact_wfs_array_targets(
        (
            propagation.field,
            propagation.focal_field,
            propagation.pupil_field,
            propagation.bi_o_edge_masks,
            propagation.phasor,
            propagation.intensity,
            propagation.temp,
            propagation.scratch,
            propagation.asterism_stack,
            propagation.fft_buffer,
            propagation.elongation_kernel,
            propagation.lgs_kernel_fft,
        ),
        (
            "Bi-O-edge field",
            "Bi-O-edge focal field",
            "Bi-O-edge pupil field",
            "Bi-O-edge amplitude masks",
            "Bi-O-edge phasor",
            "Bi-O-edge intensity",
            "Bi-O-edge temporary image",
            "Bi-O-edge scratch image",
            "Bi-O-edge asterism stack",
            "Bi-O-edge FFT buffer",
            "Bi-O-edge elongation kernel",
            "Bi-O-edge LGS Fourier kernel",
        ),
        target,
        :wfs_optics,
    )
    _require_exact_focal_plane_modulation_target(
        front_end.modulation, target, "Bi-O-edge operating modulation")
    _require_exact_focal_plane_modulation_target(
        front_end.calibration_modulation, target,
        "Bi-O-edge calibration modulation")
    return front_end
end

function _require_exact_wfs_target(
    plan::PreparedBiOEdgeOptics,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :wfs_optics)
    _require_exact_wfs_product_target(
        plan.output, target, :wfs_optics)
    _require_exact_bi_o_edge_front_end_target(plan.front_end, target)
    _require_exact_prepared_four_pupil_lgs_target(
        plan.plan.lgs_model, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedBiOEdgeOpticsBundle,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.components
        _require_exact_wfs_target(component, target)
    end
    return plan
end

function _require_exact_zernike_propagation_target(
    propagation::ZernikePropagationWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            propagation.field,
            propagation.focal_field,
            propagation.pupil_field,
            propagation.phasor,
            propagation.phase_mask,
            propagation.pupil_intensity,
            propagation.nominal_frame,
        ),
        (
            "Zernike field",
            "Zernike focal field",
            "Zernike pupil field",
            "Zernike phasor",
            "Zernike phase mask",
            "Zernike pupil intensity",
            "Zernike nominal frame",
        ),
        target,
        :wfs_optics,
    )
    return propagation
end

function _require_exact_wfs_target(
    plan::PreparedZernikeOptics,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :wfs_optics)
    _require_exact_wfs_product_target(
        plan.output, target, :wfs_optics)
    _require_exact_zernike_propagation_target(
        plan.workspace, target)
    return plan
end

function _require_exact_curvature_propagation_target(
    workspace::CurvaturePropagationWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            workspace.phasor,
            workspace.field_stack,
            workspace.defocus_stack,
            workspace.intensity_stack,
            workspace.cropped_plus,
            workspace.cropped_minus,
            workspace.frame_plus,
            workspace.frame_minus,
        ),
        (
            "Curvature phasor",
            "Curvature field stack",
            "Curvature defocus stack",
            "Curvature intensity stack",
            "Curvature positive-defocus crop",
            "Curvature negative-defocus crop",
            "Curvature positive-defocus frame",
            "Curvature negative-defocus frame",
        ),
        target,
        :wfs_optics,
    )
    return workspace
end

function _require_exact_wfs_target(
    plan::PreparedCurvatureOptics,
    target::AbstractComputeDevice,
)
    validate_wfs_optics_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :wfs_optics)
    _require_exact_wfs_product_target(
        plan.output, target, :wfs_optics)
    _require_exact_curvature_propagation_target(
        plan.workspace, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedWFSDetectorAcquisition,
    target::AbstractComputeDevice,
)
    optical_product = detector_acquisition_input(plan.acquisition)
    validate_wfs_acquisition_binding(
        plan.observation, optical_product, plan)
    _require_exact_wfs_product_target(
        optical_product, target, :acquisition)
    _require_exact_wfs_observation_target(plan.observation, target)
    _require_exact_detector_acquisition_target(plan.acquisition, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedWFSCountingAcquisition,
    target::AbstractComputeDevice,
)
    validate_wfs_acquisition_binding(
        plan.observation, plan.optical_product, plan)
    _require_exact_wfs_product_target(
        plan.optical_product, target, :acquisition)
    _require_exact_wfs_observation_target(plan.observation, target)
    _require_exact_counting_detector_target(plan.detector, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedWFSMultipleDetectorAcquisition,
    target::AbstractComputeDevice,
)
    validate_wfs_acquisition_binding(
        plan.observations, plan.optical_products, plan)
    @inbounds for component in plan.acquisitions
        _require_exact_wfs_target(component, target)
    end
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedCurvaturePackedFrameAcquisition,
    target::AbstractComputeDevice,
)
    validate_wfs_acquisition_binding(
        plan.observation, plan.optical_products, plan)
    _require_exact_wfs_product_target(
        plan.optical_products, target, :acquisition)
    _require_exact_wfs_product_target(
        plan.packed_rate, target, :acquisition)
    _require_exact_wfs_observation_target(plan.observation, target)
    _require_exact_wfs_target(plan.detector_acquisition, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedCurvaturePackedChannelAcquisition,
    target::AbstractComputeDevice,
)
    validate_wfs_acquisition_binding(
        plan.observation, plan.optical_products, plan)
    _require_exact_wfs_product_target(
        plan.optical_products, target, :acquisition)
    _require_exact_wfs_observation_target(plan.observation, target)
    _require_exact_wfs_array_targets(
        (plan.channels, plan.detector_input, plan.detector_output),
        (
            "packed Curvature channels",
            "packed Curvature detector input",
            "packed Curvature detector output",
        ),
        target,
        :acquisition,
    )
    return plan
end

function _require_exact_sh_estimator_state_target(
    sensor::ShackHartmannWFS,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            sensor.products.slopes,
            sensor.workspace.spot_stats,
            sensor.workspace.spot_stats_accum,
            sensor.calibration.reference_signal_2d,
            sensor.front_end.layout.valid_mask,
        ),
        (
            "Shack-Hartmann slopes",
            "Shack-Hartmann spot statistics",
            "Shack-Hartmann accumulated spot statistics",
            "Shack-Hartmann calibration reference",
            "Shack-Hartmann valid-subaperture mask",
        ),
        target,
        :estimation,
    )
    return sensor
end

function _require_exact_sh_mode_target(
    sensor::ShackHartmannWFS{<:Diffractive},
    target::AbstractComputeDevice,
)
    _require_exact_sh_front_end_target(sensor.optics, target)
    _require_exact_wfs_array_targets(
        (
            sensor.workspace.spot_cube,
            sensor.products.legacy_spot_cube,
            sensor.workspace.detector_noise_cube,
        ),
        (
            "Shack-Hartmann spot cube",
            "Shack-Hartmann internal legacy spot diagnostic",
            "Shack-Hartmann detector-noise cube",
        ),
        target,
        :estimation,
    )
    return sensor
end

function _require_exact_sh_mode_target(
    sensor::ShackHartmannWFS{<:Geometric},
    target::AbstractComputeDevice,
)
    _require_exact_sh_layout_target(
        sensor.front_end.layout, target, :estimation)
    return sensor
end

function _require_exact_wfs_target(
    plan::PreparedShackHartmannEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_estimator_input_target(plan.input, target)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    _require_exact_sh_estimator_state_target(plan.sensor, target)
    _require_exact_sh_mode_target(plan.sensor, target)
    return plan
end

function _require_exact_pyramid_estimator_state_target(
    state::PyramidEstimatorState,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.optical_gain,
            state.valid_i4q,
            state.reference_signal_2d,
        ),
        (
            "Pyramid valid mask",
            "Pyramid optical gain",
            "Pyramid valid I4Q mask",
            "Pyramid calibration reference",
        ),
        target,
        :estimation,
    )
    return state
end

function _require_exact_pyramid_estimator_workspace_target(
    workspace::PyramidEstimatorWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            workspace.valid_signal,
            workspace.valid_signal_indices,
            workspace.valid_flux_sum_buffer,
            workspace.flux_i4q,
            workspace.signal_2d,
        ),
        (
            "Pyramid valid-signal mask",
            "Pyramid valid-signal indices",
            "Pyramid flux-sum buffer",
            "Pyramid I4Q flux",
            "Pyramid signal",
        ),
        target,
        :estimation,
    )
    return workspace
end

function _require_exact_pyramid_estimator_products_target(
    products::PyramidEstimatorProducts,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(products.slopes, target, :estimation,
        "Pyramid slopes")
    return products
end

function _require_exact_wfs_target(
    plan::PreparedPyramidEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_estimator_input_target(plan.input, target)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    _require_exact_pyramid_estimator_state_target(
        plan.state, target)
    _require_exact_pyramid_estimator_workspace_target(
        plan.workspace, target)
    _require_exact_pyramid_estimator_products_target(
        plan.products, target)
    if plan.sensor.front_end !== nothing
        _require_exact_pyramid_front_end_target(
            plan.sensor.front_end, target)
        acquisition = pyramid_acquisition_products(plan.sensor)
        _require_exact_wfs_storage_target(acquisition.frame, target,
            :estimation, "Pyramid acquisition frame")
    end
    return plan
end

function _require_exact_bi_o_edge_estimator_state_target(
    state::BiOEdgeEstimatorState,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.edge_mask,
            state.optical_gain,
            state.valid_i4q,
            state.reference_signal_2d,
        ),
        (
            "Bi-O-edge valid mask",
            "Bi-O-edge edge mask",
            "Bi-O-edge optical gain",
            "Bi-O-edge valid I4Q mask",
            "Bi-O-edge calibration reference",
        ),
        target,
        :estimation,
    )
    return state
end

function _require_exact_bi_o_edge_estimator_workspace_target(
    workspace::BiOEdgeEstimatorWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            workspace.valid_signal,
            workspace.valid_signal_indices,
            workspace.valid_flux_sum_buffer,
            workspace.flux_i4q,
            workspace.signal_2d,
            workspace.binned_phase,
            workspace.edge_mask_binned,
        ),
        (
            "Bi-O-edge valid-signal mask",
            "Bi-O-edge valid-signal indices",
            "Bi-O-edge flux-sum buffer",
            "Bi-O-edge I4Q flux",
            "Bi-O-edge signal",
            "Bi-O-edge binned phase",
            "Bi-O-edge binned edge mask",
        ),
        target,
        :estimation,
    )
    return workspace
end

function _require_exact_bi_o_edge_estimator_products_target(
    products::BiOEdgeEstimatorProducts,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(products.slopes, target, :estimation,
        "Bi-O-edge slopes")
    return products
end

function _require_exact_wfs_target(
    plan::PreparedBiOEdgeEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_estimator_input_target(plan.input, target)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    _require_exact_bi_o_edge_estimator_state_target(
        plan.state, target)
    _require_exact_bi_o_edge_estimator_workspace_target(
        plan.workspace, target)
    _require_exact_bi_o_edge_estimator_products_target(
        plan.products, target)
    if plan.sensor.front_end !== nothing
        _require_exact_bi_o_edge_front_end_target(
            plan.sensor.front_end, target)
        acquisition = bi_o_edge_acquisition_products(plan.sensor)
        _require_exact_wfs_storage_target(acquisition.frame, target,
            :estimation, "Bi-O-edge acquisition frame")
    end
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedZernikeEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_observation_target(plan.input, target, :estimation)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    state = plan.state
    workspace = plan.workspace
    products = plan.products
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.reference_signal_2d,
        ),
        (
            "Zernike valid mask",
            "Zernike calibration reference",
        ),
        target,
        :estimation,
    )
    _require_exact_wfs_array_targets(
        (
            workspace.valid_signal_indices,
            workspace.signal_2d,
            workspace.normalization_frame,
            workspace.normalization_partials,
            workspace.normalization_sum,
        ),
        (
            "Zernike valid-signal indices",
            "Zernike signal",
            "Zernike normalization frame",
            "Zernike normalization partials",
            "Zernike normalization sum",
        ),
        target,
        :estimation,
    )
    _require_exact_wfs_storage_target(products.signal, target,
        :estimation, "Zernike signal product")
    return plan
end

function _require_exact_curvature_mapping_target(
    mapping::CurvatureImagePairMapping,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (mapping.plus, mapping.minus),
        (
            "Curvature positive-defocus estimator input",
            "Curvature negative-defocus estimator input",
        ),
        target,
        :estimation,
    )
    return mapping
end

function _require_exact_curvature_mapping_target(
    mapping::CurvatureChannelPairMapping,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(
        mapping.channels, target, :estimation,
        "Curvature channel estimator input")
    return mapping
end

function _require_exact_wfs_target(
    plan::PreparedCurvatureEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_observation_target(plan.input, target, :estimation)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    _require_exact_curvature_mapping_target(plan.plan.mapping, target)
    state = plan.state
    workspace = plan.workspace
    products = plan.products
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.reference_signal_2d,
        ),
        (
            "Curvature valid mask",
            "Curvature calibration reference",
        ),
        target,
        :estimation,
    )
    _require_exact_wfs_array_targets(
        (
            workspace.reduced_plus,
            workspace.reduced_minus,
            workspace.signal_2d,
        ),
        (
            "Curvature reduced positive-defocus image",
            "Curvature reduced negative-defocus image",
            "Curvature signal",
        ),
        target,
        :estimation,
    )
    _require_exact_wfs_storage_target(products.signal, target,
        :estimation, "Curvature signal product")
    return plan
end

function _require_exact_wfs_target(
    plan,
    ::AbstractComputeDevice,
)
    throw(WFSPreparationError(
        :preparation,
        :unsupported_target_validation,
        "no exact-target validator is defined for prepared WFS plan $(typeof(plan))",
    ))
end

"""
    validate_wfs_target(plan, target)

Qualified cold extension seam for exact-target validation of a prepared WFS
WFS optics, acquisition, or estimation plan. Maintained plans delegate to the
fail-closed validators above; custom prepared plans must add a more-specific
method that enumerates their numerical data-plane storage.
"""
@inline validate_wfs_target(
    plan, target::AbstractComputeDevice) =
    _require_exact_wfs_target(plan, target)
