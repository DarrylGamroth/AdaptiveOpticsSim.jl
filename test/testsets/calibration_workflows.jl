@testset "Modal bases and fitting" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
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

@testset "Interaction-matrix WFS output finalization" begin
    T = Float32
    tel = Telescope(
        resolution=480,
        diameter=8.0,
        central_obstruction=0.0,
        T=T,
    )
    src = Source(band=:R, magnitude=3.0, T=T)
    dm = DeformableMirror(tel; n_act=2, influence_width=T(0.2), T=T)
    wfs = PyramidWFS(
        tel;
        pupil_samples=20,
        threshold=T(0.1),
        modulation=zero(T),
        modulation_points=1,
        light_ratio=T(0.1),
        n_pix_separation=4,
        n_pix_edge=2,
        psf_centering=true,
        mode=Diffractive(),
        T=T,
    )
    initial_rows = length(slopes(wfs))
    imat = interaction_matrix(
        dm,
        wfs,
        PupilFunction(tel; T=T),
        src;
        amplitude=T(5e-9),
    )

    @test length(slopes(wfs)) < initial_rows
    @test size(imat.matrix) ==
        (length(slopes(wfs)), length(dm.state.coefs))
    @test all(isfinite, imat.matrix)

    stale_wfs = PyramidWFS(
        tel;
        pupil_samples=20,
        threshold=T(0.1),
        modulation=zero(T),
        modulation_points=1,
        light_ratio=T(0.1),
        n_pix_separation=4,
        n_pix_edge=2,
        psf_centering=true,
        mode=Diffractive(),
        T=T,
    )
    stale_out = zeros(T, initial_rows, length(dm.state.coefs))
    @test_throws DimensionMismatchError interaction_matrix!(
        stale_out,
        dm,
        stale_wfs,
        PupilFunction(tel; T=T),
        src;
        amplitude=T(5e-9),
    )
    @test length(slopes(stale_wfs)) == size(imat.matrix, 1)
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    basis = modal_basis(dm, tel; n_modes=2)
    fields = collect(Calibration.MISREG_FIELDS)
    meta, meta_fd, meta_ad = mktempdir() do root
        cd(root) do
            @test isempty(readdir())
            local_meta = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
            local_fd = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2];
                n_mis_reg=length(fields), field_order=fields,
                sensitivity=:finite_difference)
            local_ad = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2];
                n_mis_reg=length(fields), field_order=fields)
            @test isempty(readdir())
            return local_meta, local_fd, local_ad
        end
    end
    est = Calibration.estimate_misregistration(
        meta, meta.calib0.D; misregistration_zero=Misregistration())
    @test est.shift_x ≈ 0.0
    @test est.shift_y ≈ 0.0

    @test meta_ad.field_order == fields
    @test isapprox(meta_ad.calib0.D, meta_fd.calib0.D; rtol=1e-12, atol=1e-12)
    @test isapprox(meta_ad.meta.D, meta_fd.meta.D; rtol=2e-3, atol=1e-9)

    mktempdir() do root
        cache_path = joinpath(root, "meta-sensitivity.bin")
        @test_throws MethodError Calibration.compute_meta_sensitivity_matrix(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            cache_path=cache_path)
        @test_throws MethodError Calibration.SPRINT(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            save_sensitivity=false)
        @test_throws MethodError Calibration.SPRINT(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            recompute_sensitivity=true)
        @test !ispath(cache_path)
        @test isempty(readdir(root))
    end

    sampled_topology = SampledActuatorTopology(actuator_coordinates(dm)[:, 1:2])
    measured_dm = DeformableMirror(tel; topology=sampled_topology,
        influence_model=MeasuredInfluenceFunctions(Array(dm.state.modes[:, 1:2])))
    @test_throws UnsupportedAlgorithm Calibration.compute_meta_sensitivity_matrix(
        tel, measured_dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test_throws UnsupportedAlgorithm Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2, wfs_mis_registered=true)
    @test_throws InvalidConfiguration Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2, wfs_mis_registered=true,
        sensitivity=:finite_difference)
end

@testset "Calibration workflow contracts" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1)

    basis = modal_basis(dm, tel; n_modes=2)
    assert_modal_basis_contract(basis, length(dm.state.coefs), 2)

    imat = interaction_matrix(dm, wfs, pupil; amplitude=0.1)
    assert_interaction_matrix_contract(imat, length(slopes(wfs)), length(dm.state.coefs), 0.1)

    imat_basis = interaction_matrix(dm, wfs, pupil, basis.M2C;
        amplitude=0.1)
    assert_interaction_matrix_contract(imat_basis, length(slopes(wfs)), size(basis.M2C, 2), 0.1)

    aliased_commands = reshape(dm.state.coefs, :, 1)
    coefficients_before_alias_rejection = copy(dm.state.coefs)
    @test_throws InvalidConfiguration interaction_matrix(
        dm, wfs, pupil, aliased_commands; amplitude=0.1)
    @test dm.state.coefs == coefficients_before_alias_rejection

    control_matrix = ControlMatrix(imat.matrix)
    assert_control_matrix_contract(control_matrix, imat.matrix)
    noninverted_control_matrix = ControlMatrix(imat.matrix; invert=false)
    assert_control_matrix_contract(noninverted_control_matrix, imat.matrix; inverted=false)
    truncated_control_matrix = with_truncation(control_matrix, 0)
    assert_control_matrix_contract(truncated_control_matrix, imat.matrix)
    @test truncated_control_matrix.n_trunc == 0
    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)

    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1, basis=basis)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)
    @test calib.calibration.D == imat_basis.matrix

    meta = Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    assert_meta_sensitivity_contract(meta, 2)

    sprint = Calibration.SPRINT(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test sprint.meta isa Calibration.MetaSensitivity
    @test !hasfield(typeof(sprint), :cache_path)
    @test !hasfield(typeof(sprint), :save_sensitivity)
    @test !hasfield(typeof(sprint), :recompute_sensitivity)
    est = Calibration.estimate!(sprint, meta.calib0.D)
    @test est isa Misregistration
    mktempdir() do root
        cd(root) do
            refreshed = Calibration.estimate!(sprint, meta.calib0.D;
                n_update_zero_point=1, tel=tel, dm=dm, wfs=wfs,
                basis=basis.M2C[:, 1:2])
            @test refreshed isa Misregistration
            @test sprint.meta isa Calibration.MetaSensitivity
            @test isempty(readdir())
        end
    end

    diversity = fill(eltype(pupil.opd)(1e-9), size(pupil.opd))
    lift_basis = basis_from_m2c(dm, tel, basis.M2C)
    lift_forward = prepare_lift_forward_model(tel, src, lift_basis;
        diversity_opd=diversity, focal_resolution=8)
    lift = LiFT(lift_forward; iterations=2, numerical=true)
    psf_in = reference_direct_image(tel, src; zero_padding=1)
    lift_observation = LiFTObservation(lift_forward, copy(psf_in))
    coeffs = zeros(eltype(psf_in), 2)
    WavefrontSensors.reconstruct!(coeffs, lift, lift_observation;
        check_convergence=false)
    @test length(coeffs) == 2
    @test diagnostics(lift).residual_norm >= 0
end
