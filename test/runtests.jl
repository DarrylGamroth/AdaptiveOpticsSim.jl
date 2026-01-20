using Test
using AdaptiveOptics
using Random

function assert_source_interface(src)
    @test hasmethod(wavelength, Tuple{typeof(src)})
end

function assert_atmosphere_interface(atm, tel)
    @test applicable(advance!, atm, tel)
    @test applicable(propagate!, atm, tel)
end

function assert_wfs_interface(wfs, tel)
    @test applicable(update_valid_mask!, wfs, tel)
    @test applicable(measure!, wfs, tel)
end

function assert_detector_interface(det, psf)
    @test applicable(capture!, det, psf)
end

function assert_dm_interface(dm, tel)
    @test applicable(build_influence_functions!, dm, tel)
    @test applicable(apply!, dm, tel, DMAdditive())
end

@test AdaptiveOptics.PROJECT_STATUS == :in_development

@testset "Telescope and PSF" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    psf = compute_psf!(tel, src; zero_padding=2)
    @test size(psf) == (64, 64)
    @test maximum(psf) > 0
    @test isfinite(sum(psf))
end

@testset "Zernike basis" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end

@testset "Atmosphere propagation" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (32, 32)
    @test sum(abs.(tel.state.opd)) > 0

    delta = tel.params.diameter / tel.params.resolution
    ensure_psd!(atm, delta)
    psd_snapshot = copy(atm.state.psd)
    ensure_psd!(atm, delta)
    @test psd_snapshot == atm.state.psd

    atm = KolmogorovAtmosphere(tel; r0=0.1, L0=25.0)
    ensure_psd!(atm, delta)
    @test psd_snapshot != atm.state.psd
end

@testset "Sub-harmonic phase screens" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    delta = tel.params.diameter / tel.params.resolution
    rng = MersenneTwister(11)
    ws = PhaseStatsWorkspace(32; T=Float64)
    phs_base = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=false)
    rng = MersenneTwister(11)
    phs_sh = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=true)
    @test sum(abs.(phs_sh .- phs_base)) > 0
end

@testset "Multi-layer atmosphere" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_r0=[0.5, 0.5],
        wind_speed=[5.0, 10.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (16, 16)
    @test sum(abs.(tel.state.opd)) > 0
end

@testset "Deformable mirror and WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    dm.state.coefs .= 0.1
    apply!(dm, tel, DMReplace())
    @test sum(abs.(tel.state.opd)) > 0

    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    wfs = ShackHartmann(tel; n_subap=4)
    slopes = measure!(wfs, tel)
    @test length(slopes) == 2 * 4 * 4
    @test maximum(slopes) > 0
end

@testset "Calibration and control" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    @test size(imat.matrix) == (length(wfs.state.slopes), length(dm.state.coefs))

    recon = ModalReconstructor(imat; gain=1.0)
    cmd = reconstruct(recon, wfs.state.slopes)
    @test length(cmd) == length(dm.state.coefs)

    ctrl = DiscreteIntegratorController(length(wfs.state.slopes); gain=0.1, tau=0.02)
    dm_cmd = update!(ctrl, wfs.state.slopes, 0.01)
    @test length(dm_cmd) == length(wfs.state.slopes)
end

@testset "Detector" begin
    psf = fill(1.0, 8, 8)
    det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=2)
    frame = capture!(det, psf; rng=MersenneTwister(2))
    @test size(frame) == (4, 4)
    @test sum(frame) == sum(psf)

    det_tuple = Detector(integration_time=1.0, noise=(NoisePhoton(), NoiseReadout(0.5)),
        qe=1.0, binning=1)
    @test det_tuple.noise isa NoisePhotonReadout
end

@testset "Asterism PSF" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    ast = Asterism([src1, src2])
    psf = compute_psf!(tel, ast; zero_padding=2)
    @test length(tel.state.psf_list) == 2
    @test size(psf) == (32, 32)
    @test sum(psf) >= sum(tel.state.psf_list[1])
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    base_sum = sum(tel.state.pupil)
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(tel.state.pupil) < base_sum

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(tel.state.pupil) == sum(custom)

    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; n_subap=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    sh = ShackHartmann(tel; n_subap=4)
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = sh.params.n_subap * sh.params.n_subap
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5)

    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_profile = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    pyr = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4
end
@testset "Calibration vault and modal basis" begin
    D = rand(4, 3)
    vault = CalibrationVault(D)
    @test size(vault.M) == (3, 4)
    @test vault.cond > 0
    vault_trunc = with_truncation(vault, 1)
    @test vault_trunc.n_trunc == 1

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
    lift_analytic = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=false)
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3), [1, 2])
    @test size(H_analytic) == (64, 2)
    psf = compute_psf!(tel, src; zero_padding=1)
    coeffs = reconstruct(lift, psf, [1, 2])
    @test length(coeffs) == 2
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    rho = [0.0, 0.1]
    cov = phase_covariance(rho, atm)
    @test length(cov) == 2
    var = phase_variance(atm)
    @test var > 0
    psd = phase_spectrum([0.1], atm)
    @test length(psd) == 1
    screen = ft_phase_screen(atm, 8, 0.1; rng=MersenneTwister(1))
    @test size(screen) == (8, 8)
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    basis = modal_basis(dm, tel; n_modes=2)
    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    est = estimate_misregistration(meta, meta.calib0.D; misregistration_zero=Misregistration())
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

    assert_source_interface(src)
    assert_source_interface(lgs)
    assert_atmosphere_interface(atm, tel)
    assert_wfs_interface(wfs, tel)
    assert_dm_interface(dm, tel)
    assert_detector_interface(det, psf)
end
