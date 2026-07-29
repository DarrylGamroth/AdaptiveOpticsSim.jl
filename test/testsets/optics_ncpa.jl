@testset "OPD maps and NCPA" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    map = OPDMap(fill(1.0, 8, 8))
    pupil = PupilFunction(tel)
    apply_surface!(pupil, map, DMReplace())
    @test sum(pupil.opd) ≈ 64.0
    sampled_opd = fill(3e-9, 8, 8)
    physical_ncpa = NCPA(sampled_opd)
    @test surface_opd(physical_ncpa) === sampled_opd
    @test fieldnames(typeof(physical_ncpa)) == (:opd,)

    basis_default = AdaptiveOpticsSim.ncpa_basis(KLBasis(), tel, dm, atm; n_modes=2)
    basis_hht = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLHHtPSD()), tel, dm, atm; n_modes=2)
    basis_dm = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLDMModes()), tel, dm, atm; n_modes=2)
    @test basis_default ≈ basis_hht
    @test sum(abs.(basis_default .- basis_dm)) > 0

    coeffs = [1e-9, 2e-9]
    ncpa_default_kl = NCPA(tel, dm, atm; basis=KLBasis(), coefficients=coeffs)
    ncpa_hht = NCPA(tel, dm, atm; basis=KLBasis(KLHHtPSD()), coefficients=coeffs)
    ncpa_dm = NCPA(tel, dm, atm; basis=KLBasis(KLDMModes()), coefficients=coeffs)
    ncpa_zero = NCPA(tel, dm, atm)
    @test eltype(ncpa_zero.opd) == eltype(pupil_reflectivity(tel))
    @test all(iszero, ncpa_zero.opd)
    @test ncpa_default_kl.opd ≈ ncpa_hht.opd
    @test sum(abs.(ncpa_default_kl.opd .- ncpa_dm.opd)) > 0

    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 1e-9, 2e-9])
    @test size(ncpa.opd) == (8, 8)
end
