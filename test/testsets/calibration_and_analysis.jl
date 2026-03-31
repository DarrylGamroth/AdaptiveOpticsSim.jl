@testset "OOPAO parity knobs" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i + j / 10
    end

    sh_plain = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_shift = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        half_pixel_shift=true)
    sh_thresh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        threshold_cog=0.2)
    @test measure!(sh_plain, tel, src) != measure!(sh_shift, tel, src)
    @test measure!(sh_plain, tel, src) != measure!(sh_thresh, tel, src)

    pyr_auto = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    @test size(pyr_auto.state.modulation_phases, 3) == 8

    pyr_path = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=0.0,
        user_modulation_path=((1.0, 0.0), (0.0, 1.0)))
    @test size(pyr_path.state.modulation_phases, 3) == 2

    pyr_default = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    pyr_rooftop = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        rooftop=0.5, theta_rotation=0.2)
    pyr_old = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0, old_mask=true)
    @test pyr_default.state.pyramid_mask != pyr_rooftop.state.pyramid_mask
    @test pyr_default.state.pyramid_mask != pyr_old.state.pyramid_mask

    bio_plain = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    bio_gray = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bio_gray.state.bioedge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bio_plain.state.bioedge_masks != bio_gray.state.bioedge_masks

    bio_gray_slopes = measure!(bio_gray, tel, src)
    @test length(bio_gray_slopes) == 2 * 4 * 4
    @test all(isfinite, bio_gray_slopes)
end
@testset "Calibration vault and modal basis" begin
    D = rand(4, 3)
    vault = CalibrationVault(D)
    @test size(vault.M) == (3, 4)
    @test vault.cond > 0
    @test vault.effective_rank == 3
    vault_trunc = with_truncation(vault, 1)
    @test vault_trunc.n_trunc == 1

    D_sing = [1.0 0.0; 0.0 1e-12]
    vault_exact = CalibrationVault(D_sing; policy=ExactPseudoInverse())
    vault_tsvd = CalibrationVault(D_sing; policy=TSVDInverse(rtol=1e-9))
    @test vault_exact.effective_rank == 2
    @test vault_tsvd.effective_rank == 1
    @test maximum(abs, vault_exact.M .- vault_tsvd.M) > 0
    @test CalibrationVault(D_sing).policy isa TSVDInverse

    imat = InteractionMatrix(D_sing, 0.1)
    recon_exact = ModalReconstructor(imat; policy=ExactPseudoInverse())
    recon_tsvd = ModalReconstructor(imat; policy=TSVDInverse(rtol=1e-9))
    @test recon_exact.effective_rank == 2
    @test recon_tsvd.effective_rank == 1
    @test ModalReconstructor(imat).policy isa TSVDInverse

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    basis = modal_basis(dm, tel; n_modes=2)
    @test size(basis.M2C, 2) == 2
    opd = rand(16, 16)
    fit, corr, turb = fitting_error(opd, basis.projector, basis.basis)
    @test size(fit) == size(opd)
    @test size(corr) == size(opd)
    @test size(turb) == size(opd)

    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    M2C, basis_hht = kl_modal_basis(KLHHtPSD(), dm, tel, atm; n_modes=2)
    @test size(M2C, 2) == 2
    @test size(basis_hht, 3) == 2
    basis2 = modal_basis(dm, tel; n_modes=2, method=KLHHtPSD(), atm=atm)
    @test size(basis2.M2C, 2) == 2
end

@testset "OPD maps and NCPA" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    map = OPDMap(fill(1.0, 8, 8))
    apply!(map, tel, DMReplace())
    @test sum(tel.state.opd) ≈ 64.0

    basis_default = AdaptiveOpticsSim.ncpa_basis(KLBasis(), tel, dm, atm; n_modes=2)
    basis_hht = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLHHtPSD()), tel, dm, atm; n_modes=2)
    basis_dm = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLDMModes()), tel, dm, atm; n_modes=2)
    @test basis_default ≈ basis_hht
    @test sum(abs.(basis_default .- basis_dm)) > 0

    coeffs = [1e-9, 2e-9]
    ncpa_default_kl = NCPA(tel, dm, atm; basis=KLBasis(), coefficients=coeffs)
    ncpa_hht = NCPA(tel, dm, atm; basis=KLBasis(KLHHtPSD()), coefficients=coeffs)
    ncpa_dm = NCPA(tel, dm, atm; basis=KLBasis(KLDMModes()), coefficients=coeffs)
    @test ncpa_default_kl.opd ≈ ncpa_hht.opd
    @test sum(abs.(ncpa_default_kl.opd .- ncpa_dm.opd)) > 0

    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 1e-9, 2e-9])
    @test size(ncpa.opd) == (8, 8)
end

@testset "Spatial filter" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    sf = SpatialFilter(tel; shape=CircularFilter(), diameter=4, zero_padding=2)
    phase, amp = filter!(sf, tel, src)
    @test size(phase) == (8, 8)
    @test size(amp) == (8, 8)
end

@testset "Gain sensing camera" begin
    mask = ones(8, 8)
    basis = rand(8, 8, 3)
    gsc = GainSensingCamera(mask, basis)
    frame = abs.(randn(8, 8))
    calibrate!(gsc, frame)
    og = compute_optical_gains!(gsc, frame)
    @test length(og) == 3
    @test length(AdaptiveOpticsSim.weak_mode_mask(gsc)) == 3
    @test all(isfinite, og)
    @test AdaptiveOpticsSim.detector_metadata(gsc) === nothing

    weak_gsc = GainSensingCamera(mask, zeros(8, 8, 2); sensitivity_floor=1e-6)
    calibrate!(weak_gsc, frame)
    weak_og = compute_optical_gains!(weak_gsc, frame)
    @test all(AdaptiveOpticsSim.weak_mode_mask(weak_gsc))
    @test weak_og == ones(2)

    det = Detector(noise=NoiseReadout(1e-3), integration_time=2.0, qe=0.8, psf_sampling=2, binning=4)
    gsc_with_det = GainSensingCamera(mask, basis; detector=det)
    metadata = AdaptiveOpticsSim.detector_metadata(gsc_with_det)
    @test metadata isa AdaptiveOpticsSim.GSCDetectorMetadata
    @test metadata.integration_time == 2.0
    @test metadata.qe == 0.8
    @test metadata.psf_sampling == 2
    @test metadata.binning == 4
    @test metadata.noise == :readout
    @test metadata.readout_sigma == 1e-3
    @test occursin("psf_sampling=2", sprint(show, MIME"text/plain"(), gsc_with_det))

    AdaptiveOpticsSim.detach_detector!(gsc_with_det)
    @test AdaptiveOpticsSim.detector_metadata(gsc_with_det) === nothing
    AdaptiveOpticsSim.attach_detector!(gsc_with_det, det)
    @test AdaptiveOpticsSim.detector_metadata(gsc_with_det) isa AdaptiveOpticsSim.GSCDetectorMetadata
end

@testset "LiFT" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(8, 8, 3)
    diversity = zeros(8, 8)
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=true)
    H = lift_interaction_matrix(lift, zeros(3), [1, 2])
    @test size(H) == (64, 2)
    kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    lift_analytic = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=kernel)
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3), [1, 2])
    @test size(H_analytic) == (64, 2)
    psf = compute_psf!(tel, src; zero_padding=1)
    coeffs = reconstruct(lift, psf, [1, 2])
    @test length(coeffs) == 2
    @test AdaptiveOpticsSim.effective_solve_mode(AdaptiveOpticsSim.ScalarCPUStyle(), LiFTSolveAuto()) isa LiFTSolveQR
    diag = diagnostics(lift)
    @test diag.used_qr isa Bool
    @test isfinite(diag.residual_norm)
    @test isfinite(diag.weighted_residual_norm)
    @test isfinite(diag.update_norm)
    @test isfinite(diag.condition_ratio) || isinf(diag.condition_ratio)
    lift_normal = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, solve_mode=LiFTSolveNormalEquations())
    coeffs_normal = reconstruct(lift_normal, psf, [1, 2])
    @test length(coeffs_normal) == 2
    @test all(isfinite, coeffs_normal)
    @test !diagnostics(lift_normal).used_qr
    lift_damped = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTLevenbergMarquardt())
    coeffs_damped = reconstruct(lift_damped, psf, [1, 2])
    @test length(coeffs_damped) == 2
    @test all(isfinite, coeffs_damped)
    @test diagnostics(lift_damped).regularization >= 0
    lift_adaptive = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTAdaptiveLevenbergMarquardt())
    coeffs_adaptive = reconstruct(lift_adaptive, psf, [1, 2])
    @test length(coeffs_adaptive) == 2
    @test all(isfinite, coeffs_adaptive)
    @test diagnostics(lift_adaptive).regularization >= 0
    det_binned = Detector(noise=NoiseNone(), psf_sampling=1, binning=2)
    frame_binned = capture!(det_binned, psf; rng=MersenneTwister(3))
    lift_binned = LiFT(tel, src, basis, det_binned; diversity_opd=diversity, iterations=2,
        img_resolution=size(frame_binned, 1), numerical=true)
    coeffs_binned = reconstruct(lift_binned, frame_binned, [1, 2])
    @test length(coeffs_binned) == 2
    @test all(isfinite, coeffs_binned)
    det_readout = Detector(noise=NoiseReadout(1e-3), psf_sampling=1)
    lift_readout = LiFT(tel, src, basis, det_readout; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true)
    coeffs_readout = reconstruct(lift_readout, psf, [1, 2])
    @test length(coeffs_readout) == 2
    @test all(isfinite, coeffs_readout)

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    lift_sep = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=sep_kernel)
    @test lift_sep.params.object_kernel isa AdaptiveOpticsSim.LiFTSeparableObjectKernel
    dense_conv = similar(psf)
    sep_conv = similar(psf)
    tmp_conv = similar(psf)
    AdaptiveOpticsSim.conv2d_same!(dense_conv, psf, sep_kernel)
    AdaptiveOpticsSim.conv2d_same_separable!(
        sep_conv,
        tmp_conv,
        psf,
        lift_sep.params.object_kernel.row,
        lift_sep.params.object_kernel.col,
    )
    @test isapprox(sep_conv, dense_conv; rtol=1e-6, atol=1e-6)
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    rho = [0.0, 1e-6, 0.1, 1.0]
    cov = phase_covariance(rho, atm)
    @test length(cov) == length(rho)
    @test all(isfinite, cov)
    @test cov[1] >= cov[2]
    @test cov[2] <= cov[1]
    @test cov[3] <= cov[2]
    @test cov[4] <= cov[3]
    @test abs(cov[1] - cov[2]) / cov[1] < 1e-3
    var = phase_variance(atm)
    @test var > 0
    psd = phase_spectrum([0.1], atm)
    @test length(psd) == 1
    @test psd[1] > 0
    delta = tel.params.diameter / tel.params.resolution
    screen = ft_phase_screen(atm, 8, 0.1; rng=MersenneTwister(1))
    @test size(screen) == (8, 8)
    ensure_psd!(atm, delta)
    runtime_screen_rng = MersenneTwister(7)
    helper_screen_rng = MersenneTwister(7)
    advance!(atm, tel; rng=runtime_screen_rng)
    helper_screen, helper_psd = ft_phase_screen(atm, tel.params.resolution, delta; rng=helper_screen_rng, return_psd=true)
    @test helper_screen ≈ atm.state.opd
    @test helper_psd ≈ atm.state.psd

    for z in (1e-6, 1e-3, 0.1, 1.0, 4.0, 10.0, 40.0, 140.0)
        ref = SpecialFunctions.besselk(5 / 6, z)
        approx = AdaptiveOpticsSim._kv56_scalar(z)
        scaled_ref = z^(5 / 6) * ref
        scaled_approx = AdaptiveOpticsSim._scaled_kv56_scalar(z)
        @test isapprox(real(approx), ref; rtol=2e-4, atol=1e-10)
        @test isapprox(scaled_approx, scaled_ref; rtol=1e-7, atol=1e-12)
    end
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    basis = modal_basis(dm, tel; n_modes=2)
    meta = AdaptiveOpticsSim.compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    est = AdaptiveOpticsSim.estimate_misregistration(meta, meta.calib0.D; misregistration_zero=Misregistration())
    @test est.shift_x ≈ 0.0
    @test est.shift_y ≈ 0.0
end

@testset "Interface conformance" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    lgs = LGSSource()
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    wfs = ShackHartmann(tel; n_subap=2)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    det = Detector(noise=NoiseNone())
    psf = fill(1.0, 8, 8)
    opd_map = OPDMap(fill(0.1, size(tel.state.opd)))
    ncpa = NCPA(tel, dm, atm; coefficients=[0.01, -0.02])
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    modal = ModalReconstructor(imat; gain=1.0)
    mapped = MappedReconstructor(Matrix{Float64}(I, length(dm.state.coefs), length(dm.state.coefs)), imat; gain=0.5)
    ctrl = DiscreteIntegratorController(length(wfs.state.slopes); gain=0.1, tau=0.02)
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    runtime = ClosedLoopRuntime(sim, modal; rng=MersenneTwister(9))
    wfs_diffractive = ShackHartmann(tel; n_subap=2, mode=Diffractive())
    zwfs = ZernikeWFS(tel; n_subap=2)
    curv_count = CurvatureWFS(tel; n_subap=2, readout_model=CurvatureCountingReadout())
    ast = Asterism([src, Source(band=:I, magnitude=1.0, coordinates=(1.0, -45.0))])
    moving_atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0])
    infinite_atm = InfiniteMultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0], screen_resolution=33, stencil_size=35)
    @test CCDSensor <: FrameSensorType
    @test CMOSSensor <: FrameSensorType
    @test AvalancheFrameSensorType <: FrameSensorType
    @test HgCdTeAvalancheArraySensorType <: AvalancheFrameSensorType
    @test EMCCDSensor <: AvalancheFrameSensorType
    @test InGaAsSensor <: FrameSensorType
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
    @test supports_nondestructive_reads(HgCdTeAvalancheArraySensor())
    @test !supports_reference_read_subtraction(EMCCDSensor())
    @test supports_reference_read_subtraction(HgCdTeAvalancheArraySensor())
    @test !supports_readout_correction(EMCCDSensor())
    @test supports_readout_correction(HgCdTeAvalancheArraySensor())
    @test supports_read_cube(HgCdTeAvalancheArraySensor())
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceRowCommonModeCorrection()) == :reference_row_common_mode
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceColumnCommonModeCorrection()) == :reference_column_common_mode
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceOutputCommonModeCorrection(4)) == :reference_output_common_mode
    @test APDSensor <: CountingSensorType
    @test curv_count.params.readout_model isa CurvatureCountingReadout

    # IF-SRC
    assert_source_interface(src)
    assert_source_interface(lgs)
    # IF-ATM
    assert_atmosphere_interface(atm, tel)
    @test applicable(propagate!, moving_atm, tel, src)
    @test applicable(propagate!, infinite_atm, tel, src)
    assert_atmosphere_layer_interface(moving_atm.layers[1], tel, MersenneTwister(11), src)
    assert_atmosphere_layer_interface(infinite_atm.layers[1], tel, MersenneTwister(12), src)
    # IF-WFS
    assert_wfs_interface(wfs, tel)
    # IF-DM
    assert_dm_interface(dm, tel)
    # IF-DET
    assert_detector_interface(det, psf)
    # IF-OPT
    assert_optical_element_interface(opd_map, tel)
    assert_optical_element_interface(ncpa, tel)
    # IF-REC
    assert_reconstructor_interface(modal, wfs.state.slopes, length(dm.state.coefs))
    assert_reconstructor_interface(mapped, wfs.state.slopes, length(dm.state.coefs))
    # IF-CTRL
    assert_controller_interface(ctrl, wfs.state.slopes, 0.01)
    # IF-SIM
    iface = assert_control_simulation_interface(runtime)
    @test !supports_prepared_runtime(runtime)
    @test !supports_detector_output(runtime)
    @test !supports_stacked_sources(runtime)
    @test !supports_grouped_execution(runtime)
    @test simulation_interface(iface) === iface
    @test !supports_prepared_runtime(wfs, src)
    @test supports_prepared_runtime(wfs_diffractive, src)
    @test supports_prepared_runtime(wfs_diffractive, ast)
    @test supports_prepared_runtime(zwfs, src)
    @test supports_prepared_runtime(CurvatureWFS(tel; n_subap=2), src)
    @test supports_prepared_runtime(PyramidWFS(tel; n_subap=2, mode=Diffractive()), src)
    @test !supports_stacked_sources(wfs, src)
    @test supports_stacked_sources(wfs, ast)
    @test supports_stacked_sources(wfs_diffractive, ast)
    prepare_runtime_wfs!(wfs_diffractive, tel, src)
    @test wfs_diffractive.state.calibrated
    prepare_runtime_wfs!(zwfs, tel, src)
    @test zwfs.state.calibrated
    curv = CurvatureWFS(tel; n_subap=2)
    prepare_runtime_wfs!(curv, tel, src)
    @test curv.state.calibrated
end

@testset "Calibration workflow contracts" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

    basis = modal_basis(dm, tel; n_modes=2)
    assert_modal_basis_contract(basis, length(dm.state.coefs), 2)

    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    assert_interaction_matrix_contract(imat, length(wfs.state.slopes), length(dm.state.coefs), 0.1)

    imat_basis = interaction_matrix(dm, wfs, tel, basis.M2C; amplitude=0.1)
    assert_interaction_matrix_contract(imat_basis, length(wfs.state.slopes), size(basis.M2C, 2), 0.1)

    vault = CalibrationVault(imat.matrix)
    assert_calibration_vault_contract(vault, imat.matrix)
    vault_noinv = CalibrationVault(imat.matrix; invert=false)
    assert_calibration_vault_contract(vault_noinv, imat.matrix; inverted=false)
    vault_trunc = with_truncation(vault, 0)
    assert_calibration_vault_contract(vault_trunc, imat.matrix)
    @test vault_trunc.n_trunc == 0

    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1, basis=basis)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)
    @test calib.calibration.D == imat_basis.matrix

    meta = AdaptiveOpticsSim.compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    assert_meta_sensitivity_contract(meta, 2)

    sprint = AdaptiveOpticsSim.SPRINT(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test sprint.meta isa AdaptiveOpticsSim.MetaSensitivity
    est = AdaptiveOpticsSim.estimate!(sprint, meta.calib0.D)
    @test est isa Misregistration

    diversity = fill(eltype(tel.state.opd)(1e-9), size(tel.state.opd))
    lift_basis = basis_from_m2c(dm, tel, basis.M2C)
    lift = LiFT(tel, src, lift_basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=true)
    psf_in = compute_psf!(tel, src; zero_padding=1)
    coeffs = zeros(eltype(psf_in), 2)
    reconstruct!(coeffs, lift, psf_in, [1, 2]; check_convergence=false)
    @test length(coeffs) == 2
    @test diagnostics(lift).residual_norm >= 0
end

