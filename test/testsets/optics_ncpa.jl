@testset "OPD maps and NCPA" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    atm = KolmogorovAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0)
    map = OPDMap(fill(1.0, 8, 8))
    pupil = PupilFunction(tel)
    apply_surface!(pupil, map, DMReplace())
    @test sum(pupil.opd) ≈ 64.0
    sampled_opd = fill(3e-9, 8, 8)
    physical_ncpa = @inferred NCPA(sampled_opd)
    @test surface_opd(physical_ncpa) === sampled_opd
    @test fieldnames(typeof(physical_ncpa)) == (:opd,)
    @test parentmodule(typeof(physical_ncpa)) === Optics
    @test parentmodule(typeof(KLBasis())) === Calibration

    basis_default = Calibration.ncpa_basis(KLBasis(), tel, dm, atm; n_modes=2)
    basis_hht = Calibration.ncpa_basis(KLBasis(KLHHtPSD()), tel, dm, atm; n_modes=2)
    basis_dm = Calibration.ncpa_basis(KLBasis(KLDMModes()), tel, dm, atm; n_modes=2)
    basis_dm_without_atmosphere =
        Calibration.ncpa_basis(KLBasis(KLDMModes()), tel, dm; n_modes=2)
    @test basis_default ≈ basis_hht
    @test sum(abs.(basis_default .- basis_dm)) > 0
    @test basis_dm_without_atmosphere ≈ basis_dm
    @test_throws InvalidConfiguration Calibration.ncpa_basis(
        KLBasis(KLHHtPSD()), tel, dm; n_modes=2)

    modal_to_command, _ =
        kl_modal_basis(KLDMModes(), dm, tel; n_modes=2)
    @test_throws InvalidConfiguration Calibration.ncpa_basis(
        M2CBasis(), tel, dm; n_modes=2)
    external_basis = Calibration.ncpa_basis(
        M2CBasis(), tel, dm, atm; n_modes=2, M2C=modal_to_command)
    @test external_basis ≈
        basis_from_m2c(dm, tel, modal_to_command)

    coeffs = [1e-9, 2e-9]
    ncpa_default_kl = @inferred NCPA(
        tel, dm, atm; basis=KLBasis(), coefficients=coeffs)
    ncpa_hht = NCPA(tel, dm, atm; basis=KLBasis(KLHHtPSD()), coefficients=coeffs)
    ncpa_dm = NCPA(tel, dm, atm; basis=KLBasis(KLDMModes()), coefficients=coeffs)
    ncpa_zero = NCPA(tel, dm, atm)
    @test eltype(ncpa_zero.opd) == eltype(pupil_reflectivity(tel))
    @test all(iszero, ncpa_zero.opd)
    @test ncpa_default_kl.opd ≈ ncpa_hht.opd
    @test sum(abs.(ncpa_default_kl.opd .- ncpa_dm.opd)) > 0

    expanded_opd = similar(pupil.opd)
    @test @inferred(Calibration.combine_basis!(
        expanded_opd, basis_dm, coeffs, pupil_mask(tel))) === expanded_opd
    expansion_plan = ModalOPDExpansionPlan(basis_dm, pupil_mask(tel))
    @test expansion_plan.basis !== basis_dm
    @test expansion_plan.pupil_support !== pupil_mask(tel)
    @test @inferred(combine_basis!(
        expanded_opd, expansion_plan, coeffs)) === expanded_opd
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(Calibration.combine_basis!(
            expanded_opd, basis_dm, coeffs, pupil_mask(tel))) == 0
        @test @allocated(combine_basis!(
            expanded_opd, expansion_plan, coeffs)) == 0
    end

    amplitude = 2e-9
    random_ncpa = NCPA(tel, dm, atm;
        basis=KLBasis(KLDMModes()),
        f2=(amplitude, 1, 2, 1.0),
        seed=17)
    repeated_random_ncpa = NCPA(tel, dm, atm;
        basis=KLBasis(KLDMModes()),
        f2=(amplitude, 1, 2, 1.0),
        seed=17)
    @test random_ncpa.opd == repeated_random_ncpa.opd
    @test std(random_ncpa.opd[pupil_mask(tel)]) ≈ amplitude
    @test_throws InvalidConfiguration NCPA(tel, dm, atm;
        f2=(amplitude, 1, 2))

    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 1e-9, 2e-9])
    @test size(ncpa.opd) == (8, 8)
end
