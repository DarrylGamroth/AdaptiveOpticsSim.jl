@testset "Interface conformance" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    lgs = LGSSource()
    atm = KolmogorovAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    det = Detector(noise=NoiseNone())
    spad = SPADArrayDetector((8, 8); noise=NoisePhoton())
    mkid = MKIDArrayDetector(noise=NoisePhoton())
    psf = fill(1.0, 8, 8)
    linear_apd = LinearAPDDetector(
        topology=LinearAPDChannelBank(length(psf)), noise=NoisePhoton())
    opd_map = OPDMap(fill(0.1, size(pupil.opd)))
    ncpa = NCPA(fill(0.01, size(pupil.opd)))
    poly = with_spectrum(src, SpectralBundle([wavelength(src), 1.1 * wavelength(src)], [0.7, 0.3]))
    poly_common = with_spectrum(src, SpectralBundle(
        fill(wavelength(src), 2), [0.7, 0.3]))
    pyr = PyramidWFS(tel; pupil_samples=2)
    bio = BiOEdgeWFS(tel; pupil_samples=2)
    zwfs = ZernikeWFS(tel; pupil_samples=2)
    curv = CurvatureWFS(tel; pupil_samples=2)
    curv_count = CurvatureWFS(tel; pupil_samples=2,
        readout_model=CurvatureChannelReadout())
    ast = Asterism([src, Source(band=:I, magnitude=1.0, coordinates=(1.0, -45.0))])
    moving_atm = MultiLayerAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction_deg=[0.0], altitude=[0.0])
    infinite_atm = InfiniteMultiLayerAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction_deg=[0.0], altitude=[0.0], screen_resolution=33, stencil_size=35)
    @test CCDSensor <: AbstractFrameSensor
    @test CMOSSensor <: AbstractFrameSensor
    @test AbstractAvalancheFrameSensor <: AbstractFrameSensor
    @test AbstractHgCdTeSensor <: AbstractFrameSensor
    @test AbstractHgCdTeAvalancheArraySensor <: AbstractHgCdTeSensor
    @test EMCCDSensor <: AbstractAvalancheFrameSensor
    @test InGaAsSensor <: AbstractFrameSensor
    @test HgCdTeSensor <: AbstractHgCdTeSensor
    @test HgCdTeAvalancheArraySensor <: AbstractHgCdTeAvalancheArraySensor
    for removed_sensor_root in (
        :SensorType,
        :FrameSensorType,
        :CountingSensorType,
        :AvalancheFrameSensorType,
        :HgCdTeSensorType,
        :HgCdTeAvalancheArraySensorType,
        :SPADArraySensorType,
        :MKIDArraySensorType,
    )
        @test !isdefined(AdaptiveOpticsSim.Detectors, removed_sensor_root)
    end
    @test !supports_avalanche_gain(CCDSensor())
    @test !supports_sensor_glow(CMOSSensor())
    @test supports_detector_defect_maps(CMOSSensor())
    @test supports_detector_defect_maps(InGaAsSensor())
    @test supports_shutter_timing(CMOSSensor())
    @test !supports_shutter_timing(CCDSensor())
    @test !supports_detector_persistence(CMOSSensor())
    @test supports_detector_persistence(InGaAsSensor(persistence_model=ExponentialPersistence(0.1, 0.9)))
    @test !supports_detector_nonlinearity(CMOSSensor())
    @test supports_detector_nonlinearity(InGaAsSensor())
    @test !supports_nondestructive_reads(CCDSensor())
    @test supports_nondestructive_reads(
        CCDSensor(sampling_mode=SkipperSampling(4)))
    @test supports_nondestructive_reads(HgCdTeSensor())
    @test !supports_reference_read_subtraction(EMCCDSensor())
    @test supports_reference_read_subtraction(HgCdTeSensor())
    @test !supports_readout_correction(EMCCDSensor())
    @test supports_readout_correction(HgCdTeSensor())
    @test supports_read_cube(HgCdTeSensor())
    @test AdaptiveOpticsSim.Detectors.readout_correction_symbol(ReferenceRowCommonModeCorrection()) == :reference_row_common_mode
    @test AdaptiveOpticsSim.Detectors.readout_correction_symbol(ReferenceColumnCommonModeCorrection()) == :reference_column_common_mode
    @test AdaptiveOpticsSim.Detectors.readout_correction_symbol(ReferenceOutputCommonModeCorrection(4)) == :reference_output_common_mode
    @test AbstractSPADArraySensor <: AbstractCountingSensor
    @test AbstractMKIDArraySensor <: AbstractCountingSensor
    @test supports_photon_counting(mkid.params.sensor)
    @test !supports_energy_resolving(mkid.params.sensor)
    @test !supports_photon_number_resolving(mkid.params.sensor)
    @test curv_count.acquisition.plan.readout_model isa
        CurvatureChannelReadout

    # IF-SRC
    assert_source_interface(src)
    assert_source_interface(lgs)
    # IF-ATM
    assert_atmosphere_interface(atm, tel)
    assert_atmosphere_interface(moving_atm, tel)
    assert_atmosphere_interface(infinite_atm, tel)
    @test prepare_atmosphere_renderer(moving_atm, tel, src) isa
        AtmosphereDirectionRenderer
    @test prepare_atmosphere_renderer(infinite_atm, tel, src) isa
        AtmosphereDirectionRenderer
    assert_atmosphere_layer_interface(moving_atm.layers[1], tel, MersenneTwister(11), src)
    assert_atmosphere_layer_interface(infinite_atm.layers[1], tel, MersenneTwister(12), src)
    # IF-WFS
    @test applicable(update_valid_mask!, wfs, pupil)
    @test supports_valid_subaperture_mask(wfs)
    @test !supports_reference_signal(wfs)
    @test !applicable(update_valid_mask!, pyr, pupil)
    @test !applicable(measure!, pyr, pupil)
    @test !applicable(slopes, pyr)
    @test !supports_valid_subaperture_mask(pyr)
    @test !supports_reference_signal(pyr)
    pyramid_front_end = PyramidOpticalFrontEnd(pyr, src)
    pyramid_rate = pyramid_rate_map(pyramid_front_end, pupil)
    @test applicable(prepare_wfs_optics, pyramid_front_end, pupil,
        pyramid_rate)
    pyramid_optics = prepare_wfs_optics(pyramid_front_end, pupil,
        pyramid_rate)
    @test applicable(form_wfs_optical_products!, pyramid_rate, pupil,
        pyramid_optics)
    @test !applicable(update_valid_mask!, bio, pupil)
    @test !applicable(measure!, bio, pupil)
    @test !applicable(slopes, bio)
    @test !supports_valid_subaperture_mask(bio)
    @test !supports_reference_signal(bio)
    bio_front_end = BiOEdgeOpticalFrontEnd(bio, src)
    bio_rate = bi_o_edge_rate_map(bio_front_end, pupil)
    bio_optics = prepare_wfs_optics(bio_front_end, pupil, bio_rate)
    @test applicable(form_wfs_optical_products!, bio_rate, pupil,
        bio_optics)
    @test !applicable(update_valid_mask!, zwfs, pupil)
    @test !applicable(measure!, zwfs, pupil)
    @test !applicable(slopes, zwfs)
    @test !supports_valid_subaperture_mask(zwfs)
    @test !supports_reference_signal(zwfs)
    @test !applicable(update_valid_mask!, curv, pupil)
    @test !applicable(measure!, curv, pupil)
    @test !applicable(slopes, curv)
    @test !supports_valid_subaperture_mask(curv)
    @test !supports_reference_signal(curv)
    @test !applicable(update_valid_mask!, curv_count, pupil)
    @test !applicable(measure!, curv_count, pupil)
    @test !applicable(slopes, curv_count)
    @test !supports_valid_subaperture_mask(curv_count)
    @test !supports_reference_signal(curv_count)
    curvature_front_end = CurvatureOpticalFrontEnd(curv, src)
    curvature_rates = curvature_rate_maps(curvature_front_end, pupil)
    curvature_optics = prepare_wfs_optics(curvature_front_end, pupil,
        curvature_rates)
    @test applicable(form_wfs_optical_products!, curvature_rates, pupil,
        curvature_optics)
    @test supports_valid_subaperture_mask(wfs)
    @test valid_subaperture_mask(wfs) === wfs.front_end.layout.valid_mask
    @test !isdefined(WavefrontSensors, :camera_frame)
    @test !isdefined(WavefrontSensors, :shack_hartmann_detector_image)
    @test !isdefined(WavefrontSensors, :shack_hartmann_detector_image!)
    @test !applicable(wfs_detector_image, pyr)
    @test !applicable(wfs_detector_image, bio)
    @test !applicable(wfs_detector_image, zwfs)
    @test !applicable(wfs_detector_image, curv)
    @test !isdefined(WavefrontSensors, :shack_hartmann_spot_cube)
    @test !applicable(wfs_detector_image, wfs)
    # IF-DM
    assert_dm_interface(dm, tel)
    # IF-DET
    assert_detector_interface(det, psf)
    assert_detector_interface(linear_apd, vec(psf))
    assert_detector_interface(spad, psf)
    assert_detector_interface(mkid, psf)
    # IF-OPT
    assert_optical_element_interface(opd_map, tel)
    assert_optical_element_interface(ncpa, tel)
    # WFS execution capabilities
    @test !supports_prepared_runtime(wfs, src)
    @test !supports_prepared_runtime(wfs, poly)
    @test !supports_prepared_runtime(wfs, poly_common)
    @test !supports_prepared_runtime(wfs, ast)
    @test supports_prepared_runtime(zwfs, src)
    @test supports_prepared_runtime(curv, src)
    @test supports_prepared_runtime(PyramidWFS(tel; pupil_samples=2), src)
    @test supports_prepared_runtime(bio, src)
    @test !supports_detector_output(wfs, det)
    @test supports_detector_output(pyr, det)
    @test supports_detector_output(bio, det)
    @test supports_detector_output(zwfs, det)
    @test supports_detector_output(curv, det)
    @test !supports_detector_output(curv_count, det)
    @test supports_detector_output(curv_count, linear_apd)
    @test supports_detector_output(curv_count, spad)
    @test !supports_stacked_sources(wfs, src)
    @test supports_stacked_sources(wfs, ast)
    @test !supports_stacked_sources(wfs, poly)
    @test supports_stacked_sources(wfs, poly_common)
    @test !supports_grouped_execution(wfs, src)
    @test supports_grouped_execution(wfs, ast)
    @test !supports_grouped_execution(wfs, poly)
    @test supports_grouped_execution(wfs, poly_common)
    @test supports_grouped_execution(pyr, ast)
    @test supports_grouped_execution(pyr, poly)
    @test supports_grouped_execution(bio, ast)
    prepare_runtime_wfs!(zwfs, pupil, src)
    zernike_front_end = ZernikeOpticalFrontEnd(zwfs, src)
    zernike_rate = zernike_rate_map(zernike_front_end, pupil)
    zernike_optics = prepare_wfs_optics(zernike_front_end, pupil, zernike_rate)
    @test applicable(form_wfs_optical_products!, zernike_rate, pupil, zernike_optics)
    @test prepare_runtime_wfs!(curv, pupil, src) === curv
end
