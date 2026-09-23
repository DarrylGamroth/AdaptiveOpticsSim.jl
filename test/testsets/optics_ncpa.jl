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

    _, basis_kl = Calibration.modal_basis_components(
        KarhunenLoeveBasis(), dm, tel, atm; n_modes=2)
    _, basis_dm = Calibration.modal_basis_components(
        AOCModalBases.InfluenceFunctionEigenbasis(), dm, tel, nothing; n_modes=2)
    @test sum(abs.(basis_kl .- basis_dm)) > 0
    @test_throws InvalidConfiguration Calibration.modal_basis_components(
        KarhunenLoeveBasis(), dm, tel, nothing; n_modes=2)

    sampled_basis = modal_basis(dm, tel; n_modes=2, projector=false,
        method=AOCModalBases.InfluenceFunctionEigenbasis())
    external_basis = basis_from_m2c(dm, tel, sampled_basis.M2C)
    for mode in axes(external_basis, 3)
        raw_mode = @view external_basis[:, :, mode]
        centered_mode = raw_mode .- mean(raw_mode[pupil_mask(tel)])
        @test centered_mode ≈ @view(basis_dm[:, :, mode])
    end

    coefficients = [1e-9, 2e-9]
    specification = AOCModalBases.ModalOPDExpansionSpecification(
        size(basis_kl, 1), size(basis_kl, 2), size(basis_kl, 3), pupil_mask(tel), Float64)
    plan = AdaptiveOpticsCalibration.prepare(AOCModalBases.ModalOPDExpansion(), specification)
    product = AdaptiveOpticsCalibration.process(
        plan, AOCModalBases.ModalOPDExpansionInputs(basis_kl, coefficients))
    ncpa_kl = NCPA(product.opd)
    dm_plan = AdaptiveOpticsCalibration.prepare(
        AOCModalBases.ModalOPDExpansion(),
        AOCModalBases.ModalOPDExpansionSpecification(8, 8, 2, pupil_mask(tel), Float64))
    ncpa_dm = NCPA(AdaptiveOpticsCalibration.process(
        dm_plan, AOCModalBases.ModalOPDExpansionInputs(basis_dm, coefficients)).opd)
    @test sum(abs.(ncpa_kl.opd .- ncpa_dm.opd)) > 0
    @test all(iszero, NCPA(zeros(eltype(pupil_reflectivity(tel)), 8, 8)).opd)

    # The plant graph retains its separately qualified prepared executor until
    # graph ownership moves. It must agree with the cold AOC product.
    expanded_opd = similar(pupil.opd)
    expansion_plan = ModalOPDExpansionPlan(basis_dm, pupil_mask(tel))
    @test expansion_plan.basis !== basis_dm
    @test expansion_plan.pupil_support !== pupil_mask(tel)
    @test @inferred(combine_basis!(expanded_opd, expansion_plan, coefficients)) === expanded_opd
    @test expanded_opd ≈ ncpa_dm.opd
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(combine_basis!(expanded_opd, expansion_plan, coefficients)) == 0
    end

    apply_surface!(pupil, ncpa_kl, DMReplace())
    @test pupil.opd ≈ ncpa_kl.opd
end

@testset "Frozen AOS NCPA synthesis parity" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "..", "fixtures", "aos_ncpa_opd.toml"))
    @test fixture["schema"] == "test.adaptive-optics-sim/ncpa-opd/1"
    @test fixture["source_revision"] == "6fb8617e172a29affc82631df79cfc8dc303220e"
    @test fixture["matrix_order"] == "row-major"
    @test fixture["opd_unit"] == "m"

    rows = fixture["resolution"]
    row_major_matrix(values) = permutedims(reshape(values, rows, rows))
    tel = Telescope(resolution=rows, diameter=fixture["diameter_m"],
        central_obstruction=fixture["central_obstruction"])
    dm = DeformableMirror(tel; n_act=fixture["actuator_count_per_axis"],
        influence_width=fixture["influence_width"])
    zernike = ZernikeBasis(tel, fixture["mode_count"])
    compute_zernike!(zernike, tel)
    basis = zernike.modes
    support = pupil_mask(tel)
    @test support == row_major_matrix(fixture["pupil_support"])

    plan = AdaptiveOpticsCalibration.prepare(
        AOCModalBases.ModalOPDExpansion(),
        AOCModalBases.ModalOPDExpansionSpecification(
            rows, rows, fixture["mode_count"], support, Float64))
    explicit = fixture["explicit_coefficients"]
    product = AdaptiveOpticsCalibration.process(plan,
        AOCModalBases.ModalOPDExpansionInputs(basis, explicit["coefficients"]))
    @test product.opd ≈ row_major_matrix(explicit["expected_opd"]) atol=1e-22 rtol=1e-12
    @test surface_opd(NCPA(product.opd)) === product.opd

    # The historical `f2` tuple is now an explicit composing-layer policy:
    # seeded coefficients followed by a target supported-pupil sample std.
    seeded = fixture["seeded_power_law"]
    rng = runtime_rng(seeded["seed"])
    coefficients = zeros(Float64, fixture["mode_count"])
    for mode in seeded["start_mode"]:seeded["end_mode"]
        coefficients[mode] = randn(rng) / sqrt(mode + seeded["cutoff"])
    end
    @test coefficients ≈ seeded["expected_unscaled_coefficients"] atol=1e-15 rtol=1e-14
    seeded_product = AdaptiveOpticsCalibration.process(plan,
        AOCModalBases.ModalOPDExpansionInputs(basis, coefficients))
    sampled_std = std(seeded_product.opd[support])
    @test sampled_std > 0
    seeded_product.opd .*= seeded["target_supported_pupil_std_m"] / sampled_std
    @test seeded_product.opd ≈ row_major_matrix(seeded["expected_opd"]) atol=1e-21 rtol=1e-12
    @test std(seeded_product.opd[support]) ≈ seeded["target_supported_pupil_std_m"]
end
