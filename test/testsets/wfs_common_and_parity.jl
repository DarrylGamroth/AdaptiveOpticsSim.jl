struct CommonContractWFS <: WavefrontSensors.AbstractWFS end

@testset "Common WavefrontSensors ownership" begin
    for name in (
        :WFSPreparationError,
        :Diffractive,
        :Geometric,
        :WFSObservation,
        :WFSMeasurement,
        :prepare_wfs_optical_formation,
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
        @test getfield(AdaptiveOpticsSim, name) ===
            getfield(WavefrontSensors, name)
    end

    for family in (
        :ZernikeWFS,
        :CurvatureWFS,
        :LiFT,
    )
        @test !isdefined(WavefrontSensors, family)
    end
    for name in (
        :ShackHartmannWFS,
        :ShackHartmannDirectFrontEnd,
        :ShackHartmannOpticalFrontEnd,
        :SubapertureLayout,
        :SubapertureCalibration,
        :shack_hartmann_detector_image,
        :PyramidWFS,
        :BioEdgeWFS,
        :PyramidOpticalFrontEnd,
        :BioEdgeOpticalFrontEnd,
        :pyramid_rate_map,
        :bioedge_rate_map,
        :set_pyramid_calibration!,
        :set_bioedge_calibration!,
        :pyramid_modulation_frame,
        :pyramid_modulation_frame!,
    )
        @test parentmodule(getfield(WavefrontSensors, name)) ===
            WavefrontSensors
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
    end
    @test parentmodule(MicrolensArray) === Optics
    @test parentmodule(PyramidPhaseMask) === Optics
    @test parentmodule(BioEdgeAmplitudeMask) === Optics
    @test MicrolensArray(; n_lenslets=2, n_pix_subap=2) isa
        MicrolensArray

    sensor = CommonContractWFS()
    @test @inferred(WavefrontSensors.sensing_mode(sensor)) isa Diffractive
    @test !(@inferred supports_prepared_runtime(sensor, nothing))
    @test !(@inferred supports_stacked_sources(sensor, nothing))
    @test !(@inferred supports_grouped_execution(sensor, nothing))
    @test @inferred(valid_subaperture_mask(sensor)) === nothing
    @test @inferred(camera_frame(sensor)) === nothing

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
    @test occursin("include(\"bioedge.jl\")", common_entry)
    for family_source in (
        "zernike.jl",
        "curvature.jl",
        "lift.jl",
    )
        @test !occursin("include(\"$family_source\")", common_entry)
    end
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
        rooftop=0.5, theta_rotation=0.2)
    pyr_old = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0, old_mask=true)
    @test pyr_default.front_end.propagation.pyramid_mask != pyr_rooftop.front_end.propagation.pyramid_mask
    @test pyr_default.front_end.propagation.pyramid_mask != pyr_old.front_end.propagation.pyramid_mask

    bio_plain = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0)
    bio_gray = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bio_gray.front_end.propagation.bioedge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bio_plain.front_end.propagation.bioedge_masks != bio_gray.front_end.propagation.bioedge_masks

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
    @test_throws InvalidConfiguration BioEdgeWFS(tel;
        pupil_samples=5, binning=2, mode=Diffractive())
    @test_throws InvalidConfiguration PyramidWFS(tel;
        pupil_samples=0, mode=Diffractive())
    @test_throws InvalidConfiguration BioEdgeWFS(tel;
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

    bioedge = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration begin
        WavefrontSensors.resize_bioedge_signal_buffers!(bioedge, 7, 1)
    end
    @test_throws DimensionMismatchError begin
        WavefrontSensors.bioedge_signal!(bioedge, pupil, zeros(8, 6))
    end

    compact_bioedge = BioEdgeWFS(tel; pupil_samples=2,
        mode=Diffractive())
    compact_bioedge.acquisition.state.nominal_detector_resolution = 4
    WavefrontSensors.resize_bioedge_signal_buffers!(compact_bioedge, 4)
    fill!(compact_bioedge.estimator.state.valid_i4q, true)
    WavefrontSensors.update_bioedge_valid_signal!(compact_bioedge)
    WavefrontSensors.update_bioedge_valid_signal_indices!(compact_bioedge)
    WavefrontSensors.resize_bioedge_slope_buffers!(compact_bioedge)
    fill!(compact_bioedge.estimator.state.reference_signal_2d, 0.0)
    compact_frame = [4.0 4.0 1.0 1.0;
                     4.0 4.0 1.0 1.0;
                     3.0 3.0 2.0 2.0;
                     3.0 3.0 2.0 2.0]
    compact_slopes = copy(WavefrontSensors.bioedge_signal!(
        compact_bioedge, pupil, compact_frame))

    padded_bioedge = BioEdgeWFS(tel; pupil_samples=2,
        mode=Diffractive())
    padded_bioedge.acquisition.state.nominal_detector_resolution = 4
    WavefrontSensors.resize_bioedge_signal_buffers!(padded_bioedge, 8)
    fill!(padded_bioedge.estimator.state.valid_i4q, true)
    WavefrontSensors.update_bioedge_valid_signal!(padded_bioedge)
    WavefrontSensors.update_bioedge_valid_signal_indices!(padded_bioedge)
    WavefrontSensors.resize_bioedge_slope_buffers!(padded_bioedge)
    fill!(padded_bioedge.estimator.state.reference_signal_2d, 0.0)
    padded_frame = zeros(8, 8)
    @views padded_frame[3:6, 3:6] .= compact_frame
    @test WavefrontSensors.bioedge_signal!(padded_bioedge, pupil,
        padded_frame) ≈ compact_slopes

    ngs = Source(wavelength=589e-9, photon_irradiance=1.0)
    lgs = LGSSource(wavelength=589e-9, elongation_factor=1.4,
        photon_irradiance=1.0)
    heterogeneous = Asterism(AdaptiveOpticsSim.Optics.AbstractSource[ngs, lgs])
    detector = Detector(noise=NoiseNone(), integration_time=1.0,
        qe=1.0, binning=1)
    sensors = (
        ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive()),
        PyramidWFS(tel; pupil_samples=4, mode=Diffractive()),
        BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive()),
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
    @test AdaptiveOpticsSim.common_wfs_calibration_source(
        common_lgs, "test WFS") === first(common_lgs.sources)
end
