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
        :measure!,
        :slopes,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end

    for name in (
        :LiFT,
        :PreparedLiFTForwardModel,
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
        :diagnostics,
        :LiFTSolveAuto,
        :LiFTSolveQR,
        :LiFTSolveNormalEquations,
        :LiFTLevenbergMarquardt,
        :LiFTAdaptiveLevenbergMarquardt,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test Base.isexported(WavefrontSensors, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end
    @test isdefined(WavefrontSensors, :reconstruct!)
    @test isdefined(WavefrontSensors, :reconstruct)
    @test !Base.isexported(WavefrontSensors, :reconstruct!)
    @test !Base.isexported(WavefrontSensors, :reconstruct)
    @test WavefrontSensors.reconstruct! !== Control.reconstruct!
    @test WavefrontSensors.reconstruct !== Control.reconstruct
    for name in (
        :ShackHartmannWFS,
        :ShackHartmannDirectFrontEnd,
        :ShackHartmannOpticalFrontEnd,
        :SubapertureLayout,
        :SubapertureCalibration,
        :PyramidWFS,
        :BiOEdgeWFS,
        :PyramidOpticalFrontEnd,
        :BiOEdgeOpticalFrontEnd,
        :pyramid_rate_map,
        :bi_o_edge_rate_map,
        :set_pyramid_calibration!,
        :set_bi_o_edge_calibration!,
        :pyramid_modulation_frame,
        :pyramid_modulation_frame!,
        :ZernikeWFS,
        :ZernikeOpticalFrontEnd,
        :zernike_rate_map,
        :set_zernike_calibration!,
        :CurvatureWFS,
        :CurvatureOpticalFrontEnd,
        :curvature_rate_maps,
        :CurvaturePackedAcquisition,
        :set_curvature_calibration!,
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

@testset "OOPAO parity knobs" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i + j / 10
    end

    sh_plain = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8)
    sh_shift = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8,
        half_pixel_shift=true)
    sh_thresh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8,
        threshold_cog=0.2)
    @test measure!(sh_plain, pupil, src) != measure!(sh_shift, pupil, src)
    @test measure!(sh_plain, pupil, src) != measure!(sh_thresh, pupil, src)

    pyr_auto = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0)
    @test size(pyr_auto.front_end.modulation.phases, 3) == 8

    pyr_path = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=0.0,
        user_modulation_path=((1.0, 0.0), (0.0, 1.0)))
    @test size(pyr_path.front_end.modulation.phases, 3) == 2

    pyr_default = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0)
    pyr_rooftop = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0,
        rooftop=0.5, phase_mask_rotation_rad=0.2)
    pyr_old = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0, old_mask=true)
    @test pyramid_propagation_workspace(pyr_default).pyramid_mask !=
        pyramid_propagation_workspace(pyr_rooftop).pyramid_mask
    @test pyramid_propagation_workspace(pyr_default).pyramid_mask !=
        pyramid_propagation_workspace(pyr_old).pyramid_mask
    @test pyr_rooftop.front_end.phase_mask.rotation_rad == 0.2

    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=4, phase_mask_rotation_rad=NaN)
    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=4, modulation_phase_offset_rad=NaN)
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=4, modulation_phase_offset_rad=NaN)

    bio_plain = BiOEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0)
    bio_gray = BiOEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bi_o_edge_propagation_workspace(
        bio_gray).bi_o_edge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bi_o_edge_propagation_workspace(bio_plain).bi_o_edge_masks !=
        bi_o_edge_propagation_workspace(bio_gray).bi_o_edge_masks

    bio_gray_slopes = measure!(bio_gray, pupil, src)
    @test length(bio_gray_slopes) == 2 * 4 * 4
    @test all(isfinite, bio_gray_slopes)
end

@testset "WFS asterism calibration and pupil-image geometry" begin
    tel = Telescope(resolution=20, diameter=8.0,
        central_obstruction=0.0)
    pupil = PupilFunction(tel)

    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=5, binning=2, mode=Diffractive())
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=5, binning=2, mode=Diffractive())
    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=0, mode=Diffractive())
    @test_throws InvalidConfiguration BiOEdgeWFS(tel;
        pupil_samples=0, mode=Diffractive())

    pyramid = PyramidWFS(tel; pupil_samples=4,
        diffraction_padding=3, mode=Diffractive())
    WavefrontSensors.prepare_pyramid_sampling!(pyramid, pupil)
    @test_throws InvalidConfiguration begin
        WavefrontSensors.resize_pyramid_signal_buffers!(pyramid, 3)
    end
    @test_throws DimensionMismatchError begin
        WavefrontSensors.pyramid_signal!(pyramid, pupil, zeros(8, 6))
    end

    bi_o_edge = BiOEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration begin
        WavefrontSensors.resize_bi_o_edge_signal_buffers!(bi_o_edge, 7, 1)
    end
    @test_throws DimensionMismatchError begin
        WavefrontSensors.bi_o_edge_signal!(bi_o_edge, pupil, zeros(8, 6))
    end

    compact_bi_o_edge = BiOEdgeWFS(tel; pupil_samples=2,
        mode=Diffractive())
    WavefrontSensors.bi_o_edge_acquisition_workspace(
        compact_bi_o_edge).nominal_detector_resolution = 4
    WavefrontSensors.resize_bi_o_edge_signal_buffers!(compact_bi_o_edge, 4)
    fill!(compact_bi_o_edge.estimator.state.valid_i4q, true)
    WavefrontSensors.update_bi_o_edge_valid_signal!(compact_bi_o_edge)
    WavefrontSensors.update_bi_o_edge_valid_signal_indices!(compact_bi_o_edge)
    WavefrontSensors.resize_bi_o_edge_slope_buffers!(compact_bi_o_edge)
    fill!(compact_bi_o_edge.estimator.state.reference_signal_2d, 0.0)
    compact_frame = [4.0 4.0 1.0 1.0;
                     4.0 4.0 1.0 1.0;
                     3.0 3.0 2.0 2.0;
                     3.0 3.0 2.0 2.0]
    compact_slopes = copy(WavefrontSensors.bi_o_edge_signal!(
        compact_bi_o_edge, pupil, compact_frame))

    padded_bi_o_edge = BiOEdgeWFS(tel; pupil_samples=2,
        mode=Diffractive())
    WavefrontSensors.bi_o_edge_acquisition_workspace(
        padded_bi_o_edge).nominal_detector_resolution = 4
    WavefrontSensors.resize_bi_o_edge_signal_buffers!(padded_bi_o_edge, 8)
    fill!(padded_bi_o_edge.estimator.state.valid_i4q, true)
    WavefrontSensors.update_bi_o_edge_valid_signal!(padded_bi_o_edge)
    WavefrontSensors.update_bi_o_edge_valid_signal_indices!(padded_bi_o_edge)
    WavefrontSensors.resize_bi_o_edge_slope_buffers!(padded_bi_o_edge)
    fill!(padded_bi_o_edge.estimator.state.reference_signal_2d, 0.0)
    padded_frame = zeros(8, 8)
    @views padded_frame[3:6, 3:6] .= compact_frame
    @test WavefrontSensors.bi_o_edge_signal!(padded_bi_o_edge, pupil,
        padded_frame) ≈ compact_slopes

    ngs = Source(wavelength=589e-9, photon_irradiance=1.0)
    lgs = LGSSource(wavelength=589e-9, elongation_factor=1.4,
        photon_irradiance=1.0)
    heterogeneous = Asterism(AdaptiveOpticsSim.Optics.AbstractSource[ngs, lgs])
    detector = Detector(noise=NoiseNone(), exposure_duration=1.0,
        qe=1.0, binning=1)
    sensors = (
        ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive()),
        PyramidWFS(tel; pupil_samples=4, mode=Diffractive()),
        BiOEdgeWFS(tel; pupil_samples=4, mode=Diffractive()),
    )
    for wfs in sensors
        @test_throws InvalidConfiguration measure!(wfs, pupil,
            heterogeneous)
        @test_throws InvalidConfiguration measure!(wfs, pupil,
            heterogeneous, detector)
    end
    @test !sensors[1].calibration.calibrated
    @test !sensors[2].estimator.state.calibrated
    @test !sensors[3].estimator.state.calibrated

    common_lgs = Asterism([
        LGSSource(wavelength=589e-9, elongation_factor=1.4,
            coordinates=(0.0, 0.0), photon_irradiance=1.0),
        LGSSource(wavelength=589e-9, elongation_factor=1.4,
            coordinates=(3.0, 90.0), photon_irradiance=2.0),
    ])
    @test AdaptiveOpticsSim.WavefrontSensors.common_wfs_calibration_source(
        common_lgs, "test WFS") === first(common_lgs.sources)
end
