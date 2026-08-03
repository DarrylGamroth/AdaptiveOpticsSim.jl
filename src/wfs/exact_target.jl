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
    _require_exact_wfs_array_targets(
        (
            propagation.field,
            propagation.phasor,
            propagation.fft_buffer,
            propagation.fft_stack,
            propagation.intensity,
            propagation.intensity_stack,
            propagation.intensity_tmp_stack,
            propagation.temp,
            propagation.bin_buffer,
            propagation.spot,
            propagation.sampled_spot_cube,
            propagation.spot_cube_accum,
            propagation.elongation_kernel,
            propagation.lgs_kernel_fft,
            propagation.fft_asterism_stack,
            propagation.amp_scales,
            propagation.opd_to_cycles,
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
        :optical_formation,
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
    front_end::ShackHartmannOpticalFrontEnd,
    target::AbstractComputeDevice,
)
    _require_exact_sh_propagation_target(front_end.propagation, target)
    _require_exact_sh_layout_target(
        front_end.layout, target, :optical_formation)
    return front_end
end

function _require_exact_wfs_target(
    plan::PreparedShackHartmannOpticalFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :optical_formation)
    _require_exact_wfs_product_target(
        plan.output, target, :optical_formation)
    _require_exact_sh_front_end_target(plan.front_end, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedShackHartmannOpticalBundleFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.plans
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
        model.kernel, target, :optical_formation,
        "four-pupil LGS elongation kernel")
    return nothing
end

function _require_exact_prepared_four_pupil_lgs_target(
    model::PreparedFourPupilSodiumProfile,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_storage_target(
        model.kernel_fft, target, :optical_formation,
        "four-pupil sodium-profile Fourier kernel")
    return nothing
end

function _require_exact_focal_plane_modulation_target(
    modulation::PreparedFocalPlaneModulation,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_wfs_storage_target(
        modulation.phases, target, :optical_formation, "$label phases")
    return modulation
end

function _require_exact_pyramid_front_end_target(
    front_end::PyramidOpticalFrontEnd,
    target::AbstractComputeDevice,
)
    propagation = front_end.propagation
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
        :optical_formation,
    )
    _require_exact_focal_plane_modulation_target(
        front_end.modulation, target, "Pyramid operating modulation")
    _require_exact_focal_plane_modulation_target(
        front_end.calibration_modulation, target,
        "Pyramid calibration modulation")
    return front_end
end

function _require_exact_wfs_target(
    plan::PreparedPyramidOpticalFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :optical_formation)
    _require_exact_wfs_product_target(
        plan.output, target, :optical_formation)
    _require_exact_pyramid_front_end_target(plan.front_end, target)
    _require_exact_prepared_four_pupil_lgs_target(
        plan.lgs_model, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedPyramidOpticalBundleFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.plans
        _require_exact_wfs_target(component, target)
    end
    return plan
end

function _require_exact_bioedge_front_end_target(
    front_end::BioEdgeOpticalFrontEnd,
    target::AbstractComputeDevice,
)
    propagation = front_end.propagation
    _require_exact_wfs_array_targets(
        (
            propagation.field,
            propagation.focal_field,
            propagation.pupil_field,
            propagation.bioedge_masks,
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
            "BioEdge field",
            "BioEdge focal field",
            "BioEdge pupil field",
            "BioEdge amplitude masks",
            "BioEdge phasor",
            "BioEdge intensity",
            "BioEdge temporary image",
            "BioEdge scratch image",
            "BioEdge asterism stack",
            "BioEdge FFT buffer",
            "BioEdge elongation kernel",
            "BioEdge LGS Fourier kernel",
        ),
        target,
        :optical_formation,
    )
    _require_exact_focal_plane_modulation_target(
        front_end.modulation, target, "BioEdge operating modulation")
    _require_exact_focal_plane_modulation_target(
        front_end.calibration_modulation, target,
        "BioEdge calibration modulation")
    return front_end
end

function _require_exact_wfs_target(
    plan::PreparedBioEdgeOpticalFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :optical_formation)
    _require_exact_wfs_product_target(
        plan.output, target, :optical_formation)
    _require_exact_bioedge_front_end_target(plan.front_end, target)
    _require_exact_prepared_four_pupil_lgs_target(
        plan.lgs_model, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedBioEdgeOpticalBundleFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    @inbounds for component in plan.plans
        _require_exact_wfs_target(component, target)
    end
    return plan
end

function _require_exact_zernike_propagation_target(
    propagation::PreparedZernikePropagation,
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
        :optical_formation,
    )
    return propagation
end

function _require_exact_wfs_target(
    plan::PreparedZernikeOpticalFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :optical_formation)
    _require_exact_wfs_product_target(
        plan.output, target, :optical_formation)
    _require_exact_zernike_propagation_target(
        plan.front_end.propagation, target)
    return plan
end

function _require_exact_curvature_propagation_target(
    propagation::PreparedCurvaturePropagation,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            propagation.phasor,
            propagation.field_stack,
            propagation.defocus_stack,
            propagation.intensity_stack,
            propagation.cropped_plus,
            propagation.cropped_minus,
            propagation.frame_plus,
            propagation.frame_minus,
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
        :optical_formation,
    )
    _require_exact_curvature_atmospheric_target(
        propagation.atmospheric_propagation, target)
    return propagation
end

@inline _require_exact_curvature_atmospheric_target(
    ::Nothing,
    ::AbstractComputeDevice,
) = nothing

function _require_exact_curvature_atmospheric_target(
    propagation::AtmosphericFieldPropagation,
    ::AbstractComputeDevice,
)
    throw(WFSPreparationError(
        :optical_formation,
        :unsupported_target_validation,
        "cached Curvature atmospheric-field propagation requires its Atmospheres-owned exact-target validator",
    ))
end

function _require_exact_wfs_target(
    plan::PreparedCurvatureOpticalFormation,
    target::AbstractComputeDevice,
)
    validate_wfs_optical_formation_binding(
        plan.output, plan.input, plan)
    _require_exact_wfs_input_target(
        plan.input, target, :optical_formation)
    _require_exact_wfs_product_target(
        plan.output, target, :optical_formation)
    _require_exact_curvature_propagation_target(
        plan.front_end.propagation, target)
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
    _require_exact_wfs_array_targets(
        (plan.detector_input, plan.detector_output),
        ("counting WFS detector input", "counting WFS detector output"),
        target,
        :acquisition,
    )
    _require_exact_counting_detector_target(plan.detector, target)
    return plan
end

function _require_exact_wfs_target(
    plan::PreparedWFSMultipleDetectorAcquisition,
    target::AbstractComputeDevice,
)
    validate_wfs_acquisition_binding(
        plan.observations, plan.optical_products, plan)
    @inbounds for component in plan.plans
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
    estimator = sensor.estimator
    _require_exact_wfs_array_targets(
        (
            estimator.slopes,
            estimator.spot_stats,
            estimator.spot_stats_accum,
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
    _require_exact_sh_front_end_target(sensor.front_end, target)
    acquisition = sensor.acquisition
    _require_exact_wfs_array_targets(
        (
            acquisition.spot_cube,
            acquisition.exported_spot_cube,
            acquisition.detector_noise_cube,
        ),
        (
            "Shack-Hartmann spot cube",
            "Shack-Hartmann exported spot cube",
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
            state.slopes,
            state.optical_gain,
            state.valid_i4q,
            state.valid_signal,
            state.valid_signal_indices,
            state.valid_flux_sum_buffer,
            state.flux_i4q,
            state.signal_2d,
            state.reference_signal_2d,
        ),
        (
            "Pyramid valid mask",
            "Pyramid slopes",
            "Pyramid optical gain",
            "Pyramid valid I4Q mask",
            "Pyramid valid-signal mask",
            "Pyramid valid-signal indices",
            "Pyramid flux-sum buffer",
            "Pyramid I4Q flux",
            "Pyramid signal",
            "Pyramid calibration reference",
        ),
        target,
        :estimation,
    )
    return state
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
        plan.sensor.estimator.state, target)
    if plan.sensor.front_end !== nothing
        _require_exact_pyramid_front_end_target(
            plan.sensor.front_end, target)
        acquisition = plan.sensor.acquisition.state
        _require_exact_wfs_array_targets(
            (acquisition.binned_intensity, acquisition.camera_frame),
            ("Pyramid binned intensity", "Pyramid camera frame"),
            target,
            :estimation,
        )
    end
    return plan
end

function _require_exact_bioedge_estimator_state_target(
    state::BioEdgeEstimatorState,
    target::AbstractComputeDevice,
)
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.edge_mask,
            state.slopes,
            state.optical_gain,
            state.valid_i4q,
            state.valid_signal,
            state.valid_signal_indices,
            state.valid_flux_sum_buffer,
            state.flux_i4q,
            state.signal_2d,
            state.reference_signal_2d,
            state.binned_phase,
            state.edge_mask_binned,
        ),
        (
            "BioEdge valid mask",
            "BioEdge edge mask",
            "BioEdge slopes",
            "BioEdge optical gain",
            "BioEdge valid I4Q mask",
            "BioEdge valid-signal mask",
            "BioEdge valid-signal indices",
            "BioEdge flux-sum buffer",
            "BioEdge I4Q flux",
            "BioEdge signal",
            "BioEdge calibration reference",
            "BioEdge binned phase",
            "BioEdge binned edge mask",
        ),
        target,
        :estimation,
    )
    return state
end

function _require_exact_wfs_target(
    plan::PreparedBioEdgeEstimator,
    target::AbstractComputeDevice,
)
    validate_wfs_estimation_binding(
        plan.measurement, plan.input, plan)
    _require_exact_wfs_estimator_input_target(plan.input, target)
    _require_exact_wfs_measurement_target(plan.measurement, target)
    _require_exact_bioedge_estimator_state_target(
        plan.sensor.estimator.state, target)
    if plan.sensor.front_end !== nothing
        _require_exact_bioedge_front_end_target(
            plan.sensor.front_end, target)
        acquisition = plan.sensor.acquisition.state
        _require_exact_wfs_array_targets(
            (acquisition.binned_intensity, acquisition.camera_frame),
            ("BioEdge binned intensity", "BioEdge camera frame"),
            target,
            :estimation,
        )
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
    state = plan.estimator.estimator.state
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.valid_signal_indices,
            state.slopes,
            state.signal_2d,
            state.reference_signal_2d,
            state.reference_frame,
            state.normalization_frame,
            state.normalization_partials,
            state.normalization_sum,
        ),
        (
            "Zernike valid mask",
            "Zernike valid-signal indices",
            "Zernike slopes",
            "Zernike signal",
            "Zernike calibration reference",
            "Zernike reference frame",
            "Zernike normalization frame",
            "Zernike normalization partials",
            "Zernike normalization sum",
        ),
        target,
        :estimation,
    )
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
    _require_exact_curvature_mapping_target(plan.mapping, target)
    _require_exact_curvature_propagation_target(
        plan.sensor.front_end.propagation, target)
    _require_exact_wfs_storage_target(
        plan.sensor.acquisition.state.camera_frame,
        target,
        :estimation,
        "Curvature acquisition camera frame",
    )
    state = plan.sensor.estimator.state
    _require_exact_wfs_array_targets(
        (
            state.valid_mask,
            state.slopes,
            state.reduced_plus,
            state.reduced_minus,
            state.signal_2d,
            state.reference_signal_2d,
        ),
        (
            "Curvature valid mask",
            "Curvature slopes",
            "Curvature reduced positive-defocus image",
            "Curvature reduced negative-defocus image",
            "Curvature signal",
            "Curvature calibration reference",
        ),
        target,
        :estimation,
    )
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
formation, acquisition, or estimation plan. Maintained plans delegate to the
fail-closed validators above; custom prepared plans must add a more-specific
method that enumerates their numerical data-plane storage.
"""
@inline validate_wfs_target(
    plan, target::AbstractComputeDevice) =
    _require_exact_wfs_target(plan, target)
