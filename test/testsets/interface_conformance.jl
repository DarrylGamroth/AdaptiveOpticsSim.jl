@testset "Interface conformance" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    lgs = LGSSource()
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    det = Detector(noise=NoiseNone())
    spad = SPADArrayDetector((8, 8); noise=NoisePhoton())
    mkid = MKIDArrayDetector(noise=NoisePhoton())
    psf = fill(1.0, 8, 8)
    linear_apd = LinearAPDDetector(
        topology=LinearAPDChannelBank(length(psf)), noise=NoisePhoton())
    opd_map = OPDMap(fill(0.1, size(pupil.opd)))
    ncpa = NCPA(tel, dm, atm; coefficients=[0.01, -0.02])
    imat = interaction_matrix(dm, wfs, pupil; amplitude=0.1)
    modal = ModalReconstructor(imat; gain=1.0)
    factorized = FactorizedReconstructor(imat; gain=1.0)
    mapped = MappedReconstructor(Matrix{Float64}(I, length(dm.state.coefs), length(dm.state.coefs)), imat; gain=0.5)
    ctrl = DiscreteIntegratorController(length(slopes(wfs)); gain=0.1, tau=0.02)
    controlled = ControlledReconstructor(
        factorized,
        DiscreteIntegratorController(length(dm.state.coefs);
            gain=0.1, tau=0.02);
        dt=TEST_ATMOSPHERE_STEP,
    )
    wfs_diffractive = ShackHartmannWFS(tel; n_lenslets=2, mode=Diffractive())
    poly = with_spectrum(src, SpectralBundle([wavelength(src), 1.1 * wavelength(src)], [0.7, 0.3]))
    poly_common = with_spectrum(src, SpectralBundle(
        fill(wavelength(src), 2), [0.7, 0.3]))
    pyr = PyramidWFS(tel; pupil_samples=2, mode=Diffractive())
    bio = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    zwfs = ZernikeWFS(tel; pupil_samples=2)
    curv = CurvatureWFS(tel; pupil_samples=2)
    curv_count = CurvatureWFS(tel; pupil_samples=2,
        readout_model=CurvatureChannelReadout())
    ast = Asterism([src, Source(band=:I, magnitude=1.0, coordinates=(1.0, -45.0))])
    moving_atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0])
    infinite_atm = InfiniteMultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0], screen_resolution=33, stencil_size=35)
    @test CCDSensor <: FrameSensorType
    @test CMOSSensor <: FrameSensorType
    @test AvalancheFrameSensorType <: FrameSensorType
    @test HgCdTeSensorType <: FrameSensorType
    @test HgCdTeAvalancheArraySensorType <: HgCdTeSensorType
    @test EMCCDSensor <: AvalancheFrameSensorType
    @test InGaAsSensor <: FrameSensorType
    @test HgCdTeSensor <: HgCdTeSensorType
    @test HgCdTeAvalancheArraySensor <: HgCdTeAvalancheArraySensorType
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
    @test SPADArraySensorType <: CountingSensorType
    @test MKIDArraySensorType <: CountingSensorType
    @test supports_photon_counting(mkid.params.sensor)
    @test !supports_energy_resolving(mkid.params.sensor)
    @test !supports_photon_number_resolving(mkid.params.sensor)
    @test curv_count.params.readout_model isa CurvatureChannelReadout

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
    assert_wfs_interface(wfs, tel)
    assert_wfs_interface(wfs_diffractive, tel)
    assert_wfs_interface(pyr, tel)
    assert_wfs_interface(bio, tel)
    assert_wfs_interface(zwfs, tel)
    assert_wfs_interface(curv, tel)
    assert_wfs_interface(curv_count, tel)
    @test !supports_camera_frame(wfs)
    @test supports_valid_subaperture_mask(wfs)
    @test supports_reference_signal(wfs)
    @test supports_camera_frame(pyr)
    @test supports_camera_frame(bio)
    @test supports_camera_frame(zwfs)
    @test supports_camera_frame(curv)
    @test valid_subaperture_mask(wfs) === wfs.front_end.layout.valid_mask
    @test reference_signal(wfs) === wfs.calibration.reference_signal_2d
    @test camera_frame(pyr) === pyr.acquisition.state.camera_frame
    @test camera_frame(bio) === bio.acquisition.state.camera_frame
    @test camera_frame(zwfs) === zwfs.acquisition.state.camera_frame
    @test camera_frame(curv) === curv.acquisition.state.camera_frame
    @test wfs_detector_image(pyr) === pyr.acquisition.state.camera_frame
    @test wfs_detector_image(bio) === bio.acquisition.state.camera_frame
    @test wfs_detector_image(zwfs) === zwfs.acquisition.state.camera_frame
    @test wfs_detector_image(curv) === curv.acquisition.state.camera_frame
    @test_throws MethodError shack_hartmann_spot_cube(wfs)
    @test_throws InvalidConfiguration wfs_detector_image(wfs)
    sh_image = wfs_detector_image(wfs_diffractive; gap=1)
    sh_cube = shack_hartmann_spot_cube(wfs_diffractive)
    @test ndims(sh_image) == 2
    n_lenslets = microlens_array(wfs_diffractive.front_end).params.n_lenslets
    @test size(sh_image) == (n_lenslets * size(sh_cube, 2) + n_lenslets - 1,
                             n_lenslets * size(sh_cube, 3) + n_lenslets - 1)
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
    # IF-REC
    assert_reconstructor_interface(modal, slopes(wfs), length(dm.state.coefs))
    assert_reconstructor_interface(factorized, slopes(wfs), length(dm.state.coefs))
    assert_reconstructor_interface(mapped, slopes(wfs), length(dm.state.coefs))
    assert_reconstructor_interface(controlled, slopes(wfs), length(dm.state.coefs))
    @test reset_controller!(controlled) === controlled
    # IF-CTRL
    assert_controller_interface(ctrl, slopes(wfs), 0.01)
    @test supports_controller_reset(ctrl)
    # WFS execution capabilities
    @test !supports_prepared_runtime(wfs, src)
    @test supports_prepared_runtime(wfs_diffractive, src)
    @test !supports_prepared_runtime(wfs_diffractive, poly)
    @test supports_prepared_runtime(wfs_diffractive, poly_common)
    @test supports_prepared_runtime(wfs_diffractive, ast)
    @test supports_prepared_runtime(zwfs, src)
    @test supports_prepared_runtime(curv, src)
    @test supports_prepared_runtime(PyramidWFS(tel; pupil_samples=2, mode=Diffractive()), src)
    @test !supports_detector_output(wfs, det)
    @test supports_detector_output(wfs_diffractive, det)
    @test supports_detector_output(pyr, det)
    @test supports_detector_output(bio, det)
    @test supports_detector_output(zwfs, det)
    @test supports_detector_output(curv, det)
    @test !supports_detector_output(curv_count, det)
    @test supports_detector_output(curv_count, linear_apd)
    @test supports_detector_output(curv_count, spad)
    @test !supports_stacked_sources(wfs, src)
    @test supports_stacked_sources(wfs, ast)
    @test supports_stacked_sources(wfs_diffractive, ast)
    @test !supports_stacked_sources(wfs_diffractive, poly)
    @test supports_stacked_sources(wfs_diffractive, poly_common)
    @test !supports_grouped_execution(wfs_diffractive, src)
    @test supports_grouped_execution(wfs_diffractive, ast)
    @test !supports_grouped_execution(wfs_diffractive, poly)
    @test supports_grouped_execution(wfs_diffractive, poly_common)
    @test supports_grouped_execution(pyr, ast)
    @test supports_grouped_execution(pyr, poly)
    @test supports_grouped_execution(bio, ast)
    prepare_runtime_wfs!(wfs_diffractive, pupil, src)
    @test wfs_diffractive.calibration.calibrated
    prepare_runtime_wfs!(zwfs, pupil, src)
    @test zwfs.estimator.state.calibrated
    prepare_runtime_wfs!(curv, pupil, src)
    @test curv.estimator.state.calibrated
end
