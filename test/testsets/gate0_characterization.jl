@testset "Gate 0 pre-refactor characterization" begin
    root = default_gate0_reference_root()
    @test has_gate0_reference_bundle(root)
    bundle = load_gate0_reference_bundle(root)
    @test bundle.metadata["baseline"] == "julia_pre_hil_gate0"
    @test bundle.metadata["source_revision"] ==
        "7cf8564a63508ec94209dde41c7af0261abee592"

    expected_kinds = Set((
        :gate0_telescope_planes,
        :gate0_electric_field,
        :gate0_spatial_filter,
        :gate0_optic_surface,
        :gate0_radiometric_chain,
        :gate0_spectral_psf,
        :gate0_direct_science,
        :gate0_atmosphere_directions,
    ))
    @test Set(case.kind for case in bundle.cases) == expected_kinds
    for case in bundle.cases
        @test case.baseline === :julia_gate0
        result = validate_gate0_case(case)
        @test size(result.actual) == size(result.expected)
        @test result.ok
    end

    radiometric = only(case for case in bundle.cases
        if case.kind === :gate0_radiometric_chain)
    radiometric_actual = compute_gate0_actual(radiometric)
    flux = @view radiometric_actual[:, :, 1]
    field_intensity = @view radiometric_actual[:, :, 2]
    psf = @view radiometric_actual[:, :, 3]
    @test field_intensity ≈ flux rtol=1e-12 atol=1e-12
    @test sum(psf) ≈ sum(field_intensity) rtol=1e-12
    detector_cfgs = radiometric.config["detectors"]
    @test length(detector_cfgs) == 2
    for (index, detector_cfg) in enumerate(detector_cfgs)
        frame = @view radiometric_actual[:, :, index + 3]
        exposure = Float64(detector_cfg["integration_time"])
        qe = Float64(detector_cfg["qe"])
        @test frame ≈ psf .* (exposure * qe) rtol=1e-12 atol=1e-12
        @test radiometric.config["telescope"]["sampling_time"] != exposure
    end
    exposure_ratio = Float64(detector_cfgs[2]["integration_time"]) /
        Float64(detector_cfgs[1]["integration_time"])
    @views @test radiometric_actual[:, :, 5] ≈
        radiometric_actual[:, :, 4] .* exposure_ratio

    spectral_cases = [case for case in bundle.cases
        if case.kind === :gate0_spectral_psf]
    @test length(spectral_cases) == 2
    for case in spectral_cases
        actual = compute_gate0_actual(case)
        @views @test actual[:, :, end] ≈
            sum(actual[:, :, 1:end-1]; dims=3)[:, :, 1]
        wavelengths = Float64.(case.config["spectrum"]["wavelengths"])
        grid_class = String(case.config["compute"]["physical_grid_class"])
        if grid_class == "compatible"
            @test all(==(first(wavelengths)), wavelengths)
        elseif grid_class == "incompatible"
            @test !all(==(first(wavelengths)), wavelengths)
        else
            @test false
        end
    end

    science = only(case for case in bundle.cases
        if case.kind === :gate0_direct_science)
    science_actual = compute_gate0_actual(science)
    @views begin
        @test science_actual[:, :, 1] ≈ science_actual[:, :, 3]
        @test science_actual[:, :, 2] ≈ science_actual[:, :, 4]
        @test science_actual[:, :, 5] ≈
            science_actual[:, :, 3] .+ science_actual[:, :, 4]
    end

    atmosphere_cases = [case for case in bundle.cases
        if case.kind === :gate0_atmosphere_directions]
    @test length(atmosphere_cases) == 2
    for case in atmosphere_cases
        telescope_dt = Float64(case.config["telescope"]["sampling_time"])
        steps = Int(case.config["atmosphere"]["advance_steps"])
        declared_duration = Float64(
            case.config["atmosphere"]["legacy_elapsed_duration"])
        @test declared_duration == telescope_dt * steps
        actual = compute_gate0_actual(case)
        tel = gate0_telescope(case)
        @views begin
            @test actual[:, :, 2] == actual[:, :, 1] .* tel.state.pupil
            @test actual[:, :, 3] != actual[:, :, 2]
            @test actual[:, :, 4] != actual[:, :, 3]
        end
    end

    reference_ids = Set{String}()
    union!(reference_ids, (case.id for case in
        load_reference_bundle(default_reference_root()).cases))
    union!(reference_ids, (case.id for case in
        load_reference_bundle(default_specula_reference_root()).cases))
    required_ids = String.(bundle.metadata["required_reference_cases"])
    @test all(in(reference_ids), required_ids)
end
