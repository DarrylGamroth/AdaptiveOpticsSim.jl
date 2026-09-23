struct CommonContractWFS <: WavefrontSensors.AbstractWFS end

@testset "Common WavefrontSensors ownership" begin
    for name in (
        :WFSPreparationError,
        :Diffractive,
        :Geometric,
        :WFSObservation,
        :WFSMeasurement,
        :prepare_wfs_optics,
        :form_wfs_optical_products!,
        :prepare_wfs_acquisition,
        :acquire_wfs_observation!,
        :prepare_wfs_estimation,
        :estimate_wfs_measurement!,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end

    @test !isdefined(WavefrontSensors, :measure!)
    @test !isdefined(WavefrontSensors, :slopes)
    @test !isdefined(WavefrontSensors, :WFSNormalization)
    @test !isdefined(WavefrontSensors, :MeanValidFluxNormalization)
    @test !isdefined(WavefrontSensors, :IncidenceFluxNormalization)
    @test !isdefined(WavefrontSensors, :FluxThresholdValidSubapertures)


    for name in (
        :AbstractPyramidModulationPropagationStrategy,
        :PyramidPupilTiltStrategy,
        :PyramidShiftedMaskStrategy,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test Base.ispublic(WavefrontSensors, name)
        @test !Base.isexported(WavefrontSensors, name)
    end

    for name in (
        :PreparedLiFTForward,
        :LiFTObservation,
        :LiFTIdentityMapping,
        :LiFTFrameMapping,
        :LiFTPhotonRate,
        :LiFTExpectedCounts,
        :LiFTNormalizedIntensity,
        :prepare_lift_forward_model,
        :evaluate_lift_forward!,
        :predict_lift_observation!,
        :lift_forward_output,
        :lift_observation_contract,
        :LiFTForwardModel,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test Base.isexported(WavefrontSensors, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end
    for name in (:LiFT, :PreparedLiFTEstimator,
        :prepare_lift_estimator, :LiFTAnalyticJacobian,
        :LiFTNumericalJacobian, :LiFTSolveNormalEquations,
        :LiFTLevenbergMarquardt, :LiFTVarianceMapWeighting,
        :LiFTEstimationPlan, :reconstruct!, :reconstruct)
        @test !isdefined(WavefrontSensors, name)
    end
    for name in (
        :ShackHartmannWFS,
        :ShackHartmannOpticalFrontEnd,
        :SubapertureLayout,
        :PyramidWFS,
        :BiOEdgeWFS,
        :PyramidOpticalFrontEnd,
        :BiOEdgeOpticalFrontEnd,
        :pyramid_rate_map,
        :bi_o_edge_rate_map,
        :pyramid_modulation_frame,
        :pyramid_modulation_frame!,
        :ZernikeWFS,
        :ZernikeOpticalFrontEnd,
        :zernike_rate_map,
        :CurvatureWFS,
        :CurvatureOpticalFrontEnd,
        :curvature_rate_maps,
        :CurvaturePackedAcquisition,
        :CurvatureReadoutModel,
        :CurvatureFrameReadout,
        :CurvatureChannelReadout,
        :CurvatureBranchResponse,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end
    @test parentmodule(MicrolensArray) === Optics
    @test parentmodule(PyramidPhaseMask) === Optics
    @test parentmodule(BiOEdgeAmplitudeMask) === Optics
    @test MicrolensArray(; n_lenslets=2, n_pix_subap=2) isa
        MicrolensArray

    sensor = CommonContractWFS()
    @test @inferred(WavefrontSensors.sensing_mode(sensor)) isa Diffractive
    @test !(@inferred supports_prepared_runtime(sensor, nothing))
    @test !(@inferred supports_stacked_sources(sensor, nothing))
    @test !(@inferred supports_grouped_execution(sensor, nothing))
    @test @inferred(valid_subaperture_mask(sensor)) === nothing
    @test Base.ispublic(WavefrontSensors, :wfs_optical_products)
    @test !Base.isexported(WavefrontSensors, :wfs_optical_products)
    @test !isdefined(WavefrontSensors, :camera_frame)
    @test !isdefined(WavefrontSensors, :shack_hartmann_detector_image)
    @test !isdefined(WavefrontSensors, :shack_hartmann_detector_image!)
    @test !isdefined(WavefrontSensors, :shack_hartmann_spot_cube)

    observation = @inferred WFSObservation(zeros(Float32, 2, 3);
        units=:electron_count, layout=:detector_frame)
    measurement = @inferred WFSMeasurement(zeros(Float32, 4);
        units=:radian, kind=:slope)
    @test observation_metadata(observation).dimensions == (2, 3)
    @test measurement_metadata(measurement).dimensions == (4,)
    observation_storage(observation)
    measurement_storage(measurement)
    @test @allocated(observation_storage(observation)) == 0
    @test @allocated(measurement_storage(measurement)) == 0

    common_entry = read(joinpath(dirname(pathof(AdaptiveOpticsSim)), "wfs",
        "wavefront_sensors.jl"), String)
    @test occursin("include(\"shack_hartmann.jl\")", common_entry)
    @test occursin("include(\"pyramid.jl\")", common_entry)
    @test occursin("include(\"bi_o_edge.jl\")", common_entry)
    @test occursin("include(\"zernike.jl\")", common_entry)
    @test occursin("include(\"curvature.jl\")", common_entry)
    @test occursin("include(\"lift.jl\")", common_entry)
end

@testset "WFS convention controls" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i + j / 10
    end

    pyr_auto = PyramidWFS(tel; pupil_samples=4, modulation=1.0)
    @test size(pyr_auto.front_end.modulation.phases, 3) == 8

    pyr_path = PyramidWFS(tel; pupil_samples=4, modulation=0.0,
        user_modulation_path=((1.0, 0.0), (0.0, 1.0)))
    @test size(pyr_path.front_end.modulation.phases, 3) == 2

    pyr_default = PyramidWFS(tel; pupil_samples=4, modulation=1.0)
    pyr_rooftop = PyramidWFS(tel; pupil_samples=4, modulation=1.0,
        rooftop=0.5, phase_mask_rotation_rad=0.2)
    pyr_old = PyramidWFS(tel; pupil_samples=4, modulation=1.0, old_mask=true)
    @test pyramid_propagation_workspace(pyr_default).pyramid_mask !=
        pyramid_propagation_workspace(pyr_rooftop).pyramid_mask
    @test pyramid_propagation_workspace(pyr_default).pyramid_mask !=
        pyramid_propagation_workspace(pyr_old).pyramid_mask
    @test pyr_rooftop.front_end.phase_mask.rotation_rad == 0.2
    @test WavefrontSensors.pyramid_propagation_plan(
        pyr_default).modulation_propagation_strategy isa
        WavefrontSensors.PyramidPupilTiltStrategy
    @test WavefrontSensors._pyramid_modulation_batch_size(
        32, 8 * 1024 * 1024) == 32
    @test WavefrontSensors._pyramid_modulation_batch_size(
        32, 16 * 1024 * 1024) == 16
    @test WavefrontSensors._pyramid_modulation_batch_size(
        32, 64 * 1024 * 1024) == 4
    @test WavefrontSensors._pyramid_modulation_batch_size(
        12, 40 * 1024 * 1024) == 6

    shifted_path = ((1.3, 0.7),)
    shifted_common = (
        pupil_samples=4,
        modulation=0.0,
        user_modulation_path=shifted_path,
        diffraction_padding=4,
    )
    pupil_tilt = PyramidWFS(
        tel;
        shifted_common...,
        modulation_propagation_strategy=
            WavefrontSensors.PyramidPupilTiltStrategy(),
    )
    shifted_mask = PyramidWFS(
        tel;
        shifted_common...,
        modulation_propagation_strategy=
            WavefrontSensors.PyramidShiftedMaskStrategy(),
    )
    separable_batch = pyramid_propagation_workspace(
        shifted_mask).modulation_batch
    @test separable_batch isa
        WavefrontSensors.PyramidSeparableShiftedMaskModulationWorkspace
    point_count = size(separable_batch.axis_1_factors, 2)
    pad = size(separable_batch.axis_1_factors, 1)
    full_masks = similar(
        pyramid_propagation_workspace(shifted_mask).pyramid_mask,
        pad,
        pad,
        point_count,
    )
    full_batch = WavefrontSensors.PyramidShiftedMaskModulationWorkspace(
        similar(separable_batch.field_stack),
        full_masks,
        separable_batch.operating_weights,
        separable_batch.axis_1_shifts_rad,
        separable_batch.axis_2_shifts_rad,
        separable_batch.bfft_plan,
        separable_batch.batch_size,
    )
    WavefrontSensors._build_pyramid_shifted_masks!(
        AdaptiveOpticsSim.Backends.ScalarCPUStyle(),
        full_batch,
        shifted_mask,
        PupilFunction(tel),
    )
    for point in 1:point_count
        separable_mask = separable_batch.axis_1_factors[:, point] *
            transpose(separable_batch.axis_2_factors[:, point])
        @test separable_mask ≈ full_masks[:, :, point] rtol = 2e-14
    end
    shifted_general_mask = PyramidWFS(
        tel;
        shifted_common...,
        rooftop=0.25,
        phase_mask_rotation_rad=0.1,
        modulation_propagation_strategy=
            WavefrontSensors.PyramidShiftedMaskStrategy(),
    )
    @test pyramid_propagation_workspace(
        shifted_general_mask).modulation_batch isa
        WavefrontSensors.PyramidShiftedMaskModulationWorkspace
    strategy_pupil = PupilFunction(tel)
    strategy_resolution = tel.params.resolution
    @inbounds for i in 1:strategy_resolution, j in 1:strategy_resolution
        x = 2 * (i - 1) / (strategy_resolution - 1) - 1
        y = 2 * (j - 1) / (strategy_resolution - 1) - 1
        strategy_pupil.opd[i, j] = 40e-9 *
            (sinpi(x) + 0.4 * cospi(y) + 0.2 * x * y)
    end
    shifted_general_front_end = PyramidOpticalFrontEnd(
        shifted_general_mask, src)
    shifted_general_rate = pyramid_rate_map(
        shifted_general_front_end, strategy_pupil)
    shifted_general_prepared = prepare_wfs_optics(
        shifted_general_front_end,
        strategy_pupil,
        shifted_general_rate,
    )
    form_wfs_optical_products!(
        shifted_general_rate,
        strategy_pupil,
        shifted_general_prepared,
    )
    @test all(isfinite, shifted_general_rate.values)
    pupil_tilt_front_end = PyramidOpticalFrontEnd(pupil_tilt, src)
    shifted_mask_front_end = PyramidOpticalFrontEnd(shifted_mask, src)
    pupil_tilt_rate = pyramid_rate_map(pupil_tilt_front_end, strategy_pupil)
    shifted_mask_rate = pyramid_rate_map(shifted_mask_front_end, strategy_pupil)
    pupil_tilt_prepared = prepare_wfs_optics(
        pupil_tilt_front_end,
        strategy_pupil,
        pupil_tilt_rate,
    )
    shifted_mask_prepared = prepare_wfs_optics(
        shifted_mask_front_end,
        strategy_pupil,
        shifted_mask_rate,
    )
    form_wfs_optical_products!(
        pupil_tilt_rate,
        strategy_pupil,
        pupil_tilt_prepared,
    )
    form_wfs_optical_products!(
        shifted_mask_rate,
        strategy_pupil,
        shifted_mask_prepared,
    )
    relative_shifted_mask_error = norm(
        shifted_mask_rate.values .- pupil_tilt_rate.values,
    ) / norm(pupil_tilt_rate.values)
    @test sum(shifted_mask_rate.values) ≈ sum(pupil_tilt_rate.values) rtol = 1e-12
    @test 0 < relative_shifted_mask_error < 0.02
    form_wfs_optical_products!(
        shifted_mask_rate,
        strategy_pupil,
        shifted_mask_prepared,
    )
    @test @allocated(form_wfs_optical_products!(
        shifted_mask_rate,
        strategy_pupil,
        shifted_mask_prepared,
    )) == 0

    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=4, phase_mask_rotation_rad=NaN)
    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=4, modulation_phase_offset_rad=NaN)
    @test_throws InvalidConfiguration PyramidWFS(
        tel;
        pupil_samples=4,
        old_mask=true,
        modulation_propagation_strategy=
            WavefrontSensors.PyramidShiftedMaskStrategy(),
    )
    @test_throws InvalidConfiguration PyramidWFS(
        tel;
        pupil_samples=4,
        psf_centering=false,
        modulation_propagation_strategy=
            WavefrontSensors.PyramidShiftedMaskStrategy(),
    )
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=4, modulation_phase_offset_rad=NaN)

    bio_plain = BiOEdgeWFS(tel; pupil_samples=4, modulation=1.0)
    bio_gray = BiOEdgeWFS(tel; pupil_samples=4, modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bi_o_edge_propagation_workspace(
        bio_gray).bi_o_edge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bi_o_edge_propagation_workspace(bio_plain).bi_o_edge_masks !=
        bi_o_edge_propagation_workspace(bio_gray).bi_o_edge_masks

    bio_gray_front_end = BiOEdgeOpticalFrontEnd(bio_gray, src)
    bio_gray_rate = bi_o_edge_rate_map(bio_gray_front_end, pupil)
    bio_gray_plan = prepare_wfs_optics(
        bio_gray_front_end, pupil, bio_gray_rate)
    form_wfs_optical_products!(bio_gray_rate, pupil, bio_gray_plan)
    @test all(isfinite, bio_gray_rate.values)
end

@testset "WFS asterism and pupil-image geometry" begin
    tel = Telescope(resolution=20, diameter=8.0,
        central_obstruction=0.0)
    pupil = PupilFunction(tel)

    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=5, binning=2)
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=5, binning=2)
    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=0)
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=0)

    pyramid = PyramidWFS(tel; pupil_samples=4, diffraction_padding=3)
    WavefrontSensors.prepare_pyramid_sampling!(pyramid, pupil)
    @test size(pyramid_acquisition_products(pyramid).frame) == (12, 12)

    bi_o_edge = BiOEdgeWFS(tel; pupil_samples=4, diffraction_padding=3)
    WavefrontSensors.prepare_bi_o_edge_sampling!(bi_o_edge, pupil)
    @test size(WavefrontSensors.bi_o_edge_acquisition_products(
        bi_o_edge).frame) == (24, 24)

    ngs = Source(wavelength=589e-9, photon_irradiance=1.0)
    lgs = LGSSource(wavelength=589e-9, elongation_factor=1.4,
        photon_irradiance=1.0)
    heterogeneous = Asterism(AdaptiveOpticsSim.Optics.AbstractSource[ngs, lgs])
    bio_sensor = BiOEdgeWFS(tel; pupil_samples=4)
    heterogeneous_front_end = BiOEdgeOpticalFrontEnd(
        bio_sensor, heterogeneous)
    @test_throws WFSPreparationError bi_o_edge_rate_map(
        heterogeneous_front_end, pupil)

    common_lgs = Asterism([
        LGSSource(wavelength=589e-9, elongation_factor=1.4,
            separation_arcsec=0.0, position_angle_deg=0.0, photon_irradiance=1.0),
        LGSSource(wavelength=589e-9, elongation_factor=1.4,
            separation_arcsec=3.0, position_angle_deg=90.0, photon_irradiance=2.0),
    ])
    @test AdaptiveOpticsSim.WavefrontSensors.common_wfs_calibration_source(
        common_lgs, "test WFS") === first(common_lgs.sources)
end
