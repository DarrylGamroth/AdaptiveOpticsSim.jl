@inline sh_resource_bytes(array) = UInt64(sizeof(eltype(array)) * length(array))

function sh_resource_bytes(arrays::Tuple)
    total = UInt64(0)
    for array in arrays
        total += sh_resource_bytes(array)
    end
    return total
end

function sh_resource_fixture()
    telescope = Telescope(resolution=16, diameter=8.0, T=Float32)
    pupil = PupilFunction(telescope; T=Float32)
    source = Source(band=:I, magnitude=7.0, T=Float32)
    sensor = ShackHartmannWFS(telescope;
        n_lenslets=4,
        n_pix_subap=4,
        mode=Diffractive(),
        T=Float32,
    )
    front_end = ShackHartmannOpticalFrontEnd(sensor.front_end, source)
    output = shack_hartmann_rate_map(front_end, pupil)
    optical_formation = prepare_wfs_optical_formation(front_end, pupil,
        output)
    set_subaperture_calibration!(sensor.calibration,
        zeros(Float32, size(sensor.calibration.reference_signal_2d));
        centroid_response=1f0, wavelength=500f-9)
    detector = Detector(integration_time=1f0, noise=NoiseNone(), qe=1f0,
        response_model=NullFrameResponse(), T=Float32)
    observation = WFSObservation(zeros(Float32, size(output.values));
        units=:detected_electrons, layout=:lenslet_mosaic)
    acquisition = prepare_wfs_acquisition(detector, output, observation)
    measurement = WFSMeasurement(zeros(Float32,
        size(sensor.estimator.slopes)); units=:pixel, kind=:centroid_slopes)
    estimator = prepare_wfs_estimation(sensor, observation, measurement)
    products = AcquisitionProducts(observation, measurement;
        metadata=(kind=:shack_hartmann,))
    return (; sensor, output, optical_formation, detector, observation,
        acquisition, measurement, estimator, products)
end

@testset "Gate 9A Shack-Hartmann structural byte reporters" begin
    target = HostComputeDevice()
    fixture = sh_resource_fixture()
    sensor = fixture.sensor
    propagation = sensor.front_end.propagation
    detector_state = fixture.detector.state

    wrappers = (
        (fixture.optical_formation,
            StructuralResourceOwnerID(:wfs_optical_formation, :sh)),
        (fixture.acquisition,
            StructuralResourceOwnerID(:detector_acquisition, :sh)),
        (Detectors.detector_acquisition_plan(
                fixture.acquisition.acquisition),
            StructuralResourceOwnerID(:detector_acquisition, :plan)),
        (fixture.estimator,
            StructuralResourceOwnerID(:wfs_estimator, :prepared)),
    )
    for (owner, id) in wrappers
        fact = structural_resource_fact(owner, id, target)
        @test structural_resource_known(fact)
        @test structural_resource_owner_id(fact) == id
        @test structural_resident_bytes(fact) == UInt64(0)
        @test structural_workspace_bytes(fact) == UInt64(0)
    end

    output_fact = structural_resource_fact(fixture.output,
        StructuralResourceOwnerID(:acquisition_product, :optical_rate), target)
    observation_fact = structural_resource_fact(fixture.observation,
        StructuralResourceOwnerID(:acquisition_product, :observation), target)
    measurement_fact = structural_resource_fact(fixture.measurement,
        StructuralResourceOwnerID(:acquisition_product, :measurement), target)
    @test structural_resident_bytes(output_fact) ==
        sh_resource_bytes(fixture.output.values)
    @test structural_resident_bytes(observation_fact) ==
        sh_resource_bytes(fixture.observation.storage)
    @test structural_resident_bytes(measurement_fact) ==
        sh_resource_bytes(fixture.measurement.storage)
    products_fact = structural_resource_fact(fixture.products,
        StructuralResourceOwnerID(:acquisition_product, :container), target)
    @test structural_resident_bytes(products_fact) == sh_resource_bytes((
        fixture.observation.storage,
        fixture.measurement.storage,
    ))
    @test structural_workspace_bytes(products_fact) == UInt64(0)

    acquisition_fact = structural_resource_fact(sensor.acquisition,
        StructuralResourceOwnerID(:wfs_estimator, :acquisition_state), target)
    @test structural_resident_bytes(acquisition_fact) == sh_resource_bytes((
        sensor.acquisition.exported_spot_cube,
    ))
    @test structural_workspace_bytes(acquisition_fact) == sh_resource_bytes((
        sensor.acquisition.spot_cube,
        sensor.acquisition.detector_noise_cube,
    ))

    estimator_fact = structural_resource_fact(sensor.estimator,
        StructuralResourceOwnerID(:wfs_estimator, :state), target)
    @test structural_resident_bytes(estimator_fact) == sh_resource_bytes((
        sensor.estimator.slopes,
    ))
    @test structural_workspace_bytes(estimator_fact) == sh_resource_bytes((
        sensor.estimator.spot_stats,
        sensor.estimator.spot_stats_accum,
        sensor.estimator.slopes_host,
        sensor.estimator.centroid_host,
    ))

    calibration_fact = structural_resource_fact(sensor.calibration,
        StructuralResourceOwnerID(:wfs_estimator, :calibration), target)
    @test structural_resident_bytes(calibration_fact) == sh_resource_bytes((
        sensor.calibration.reference_signal_2d,
        sensor.calibration.reference_signal_host,
    ))

    layout_fact = structural_resource_fact(sensor.front_end.layout,
        StructuralResourceOwnerID(:wfs_estimator, :layout), target)
    @test structural_resident_bytes(layout_fact) == sh_resource_bytes((
        sensor.front_end.layout.valid_mask,
        sensor.front_end.layout.valid_mask_host,
        sensor.front_end.layout.valid_indices_host,
    ))

    detector_fact = structural_resource_fact(detector_state,
        StructuralResourceOwnerID(:detector_state, :sh_camera), target)
    @test structural_resident_bytes(detector_fact) == sh_resource_bytes((
        detector_state.accum_buffer,
        detector_state.latent_buffer,
    ))
    @test iszero(structural_workspace_bytes(detector_fact))

    workspace_fact = structural_resource_fact(propagation,
        StructuralResourceOwnerID(:workspace, :microlens), target)
    @test structural_resident_bytes(workspace_fact) == UInt64(0)
    @test structural_workspace_bytes(workspace_fact) == sh_resource_bytes((
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
        propagation.amp_scales_host,
        propagation.opd_to_cycles,
        propagation.opd_to_cycles_host,
    ))

    unsupported = Detector(integration_time=1f0, noise=NoiseNone(),
        qe=1f0, response_model=NullFrameResponse(), output_type=UInt16,
        T=Float32)
    unsupported_fact = structural_resource_fact(unsupported.state,
        StructuralResourceOwnerID(:detector_state, :unsupported), target)
    @test structural_resource_known(unsupported_fact)
    @test structural_resident_bytes(unsupported_fact) == sh_resource_bytes((
        unsupported.state.accum_buffer,
        unsupported.state.latent_buffer,
    ))

    wrong_target = AcceleratorComputeDevice(CUDABackend(), 0)
    for (owner, id) in (
        (fixture.output,
            StructuralResourceOwnerID(:acquisition_product, :wrong_output)),
        (sensor.acquisition,
            StructuralResourceOwnerID(:wfs_estimator, :wrong_acquisition)),
        (detector_state,
            StructuralResourceOwnerID(:detector_state, :wrong_camera)),
        (propagation,
            StructuralResourceOwnerID(:workspace, :wrong_microlens)),
    )
        fact = structural_resource_fact(owner, id, wrong_target)
        @test !structural_resource_known(fact)
        @test structural_resource_unknown_reason(fact) ==
            :owner_not_on_device
    end
end
