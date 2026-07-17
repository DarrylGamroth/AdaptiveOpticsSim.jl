@testset "Telemetry and config" begin
    telemetry = AdaptiveOpticsSim.Telemetry()
    AdaptiveOpticsSim.record!(telemetry, 1, 0.0; wfe_rms=1.0, strehl=0.5, loop_gain=0.2)
    @test length(telemetry) == 1
    row = telemetry[1]
    @test row.iter == 1
    @test row.strehl ≈ 0.5

    @test Tables.istable(AdaptiveOpticsSim.Telemetry)
    @test Tables.rowaccess(AdaptiveOpticsSim.Telemetry)
    rows = collect(Tables.rows(telemetry))
    @test length(rows) == 1
    @test rows[1].loop_gain ≈ 0.2

    closed_loop = AdaptiveOpticsSim.ClosedLoopTrace(Float32[
        100 80 0.50 3 4
        120 90 0.45 5 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(closed_loop) == 2
    @test eltype(typeof(closed_loop)) == AdaptiveOpticsSim.ClosedLoopTraceRow{Float32}
    @test closed_loop[2].iter == 2
    @test closed_loop[2].t ≈ 0.102f0
    @test closed_loop[2].residual_rms_nm ≈ 90.0f0
    @test Tables.istable(typeof(closed_loop))
    @test Tables.rowaccess(typeof(closed_loop))
    closed_rows = collect(Tables.rows(closed_loop))
    @test length(closed_rows) == 2
    @test closed_rows[1].command_norm ≈ 4.0f0

    gsc_closed_loop = AdaptiveOpticsSim.GSCClosedLoopTrace(Float32[
        100 80 0.50 3 0.9 4
        120 90 0.45 5 0.8 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(gsc_closed_loop) == 2
    @test eltype(typeof(gsc_closed_loop)) == AdaptiveOpticsSim.GSCClosedLoopTraceRow{Float32}
    @test gsc_closed_loop[2].mean_optical_gain ≈ 0.8f0
    @test Tables.istable(typeof(gsc_closed_loop))
    @test Tables.rowaccess(typeof(gsc_closed_loop))
    gsc_closed_rows = collect(Tables.rows(gsc_closed_loop))
    @test length(gsc_closed_rows) == 2
    @test gsc_closed_rows[1].slope_norm ≈ 3.0f0

    replay = AdaptiveOpticsSim.GSCAtmosphereReplayTrace(Float32[
        140 100 110 0.50 0.45 3 0.9
        150 105 115 0.48 0.42 5 0.8
    ]; dt=0.002f0, t0=0.1f0)
    @test length(replay) == 2
    @test eltype(typeof(replay)) == AdaptiveOpticsSim.GSCAtmosphereReplayTraceRow{Float32}
    @test replay[2].sci_strehl ≈ 0.42f0
    @test replay[2].slope_norm ≈ 5.0f0
    @test Tables.istable(typeof(replay))
    @test Tables.rowaccess(typeof(replay))
    replay_rows = collect(Tables.rows(replay))
    @test length(replay_rows) == 2
    @test replay_rows[1].ngs_forcing_rms_nm ≈ 140.0f0

    @test_throws DimensionMismatchError AdaptiveOpticsSim.ClosedLoopTrace(zeros(2, 4))
    @test_throws DimensionMismatchError AdaptiveOpticsSim.GSCClosedLoopTrace(zeros(2, 5))
    @test_throws DimensionMismatchError AdaptiveOpticsSim.GSCAtmosphereReplayTrace(zeros(2, 6))

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    cfg = AdaptiveOpticsSim.snapshot_config(tel=tel, src=src)
    @test haskey(cfg, "tel")
    @test haskey(cfg, "src")
    path, io = mktemp()
    close(io)
    AdaptiveOpticsSim.write_config_toml(path, cfg)
    parsed = TOML.parsefile(path)
    @test parsed["tel"]["resolution"] == tel.params.resolution
    @test parsed["src"]["radiometry"] == "physical_photon_irradiance"
    @test parsed["src"]["radiometric_value"] ==
        source_radiometric_value(src)
end

@testset "Reference compare conventions" begin
    convention = parse_reference_compare_convention(Dict(
        "convention" => "oopao_geometric_sh_signal_2d",
    ))
    # AdaptiveOpticsSim stores [axis 1; axis 2], with each 2-by-2 block in
    # Julia column-major order. OOPAO stores [axis 2; axis 1], with each block
    # in NumPy row-major order. Keep the fixture asymmetric so that this test
    # detects both a block swap and an accidental transpose.
    raw = Float64[1, 2, 3, 4, 11, 12, 13, 14]
    adapted = adapt_compare_convention(convention, raw)
    @test adapted ≈ OOPAO_GEOMETRIC_SH_SLOPE_SCALE .*
        Float64[11, 13, 12, 14, 1, 3, 2, 4]

    legacy = parse_reference_compare_convention(Dict(
        "swap_halves" => true,
        "scale" => OOPAO_GEOMETRIC_SH_SLOPE_SCALE,
    ))
    @test adapt_compare_convention(legacy, raw) == adapted
    @test parse_reference_compare_convention(nothing) isa IdentityCompareConvention
    @test adapt_compare_convention(IdentityCompareConvention(), raw) === raw
    @test_throws InvalidConfiguration parse_reference_compare_convention(Dict(
        "swap_halves" => true,
    ))
    @test_throws InvalidConfiguration adapt_compare_convention(convention,
        Float64[1, 2, 3, 4, 5, 6])
end

@testset "Reference storage and detector conventions" begin
    col_data = Float64[1, 2, 3, 4]
    row_data = Float64[1, 2, 3, 4, 5, 6]
    @test reshape_reference_data(col_data, (2, 2), JuliaColumnMajorStorage()) == [1 3; 2 4]
    @test reshape_reference_data(row_data, (2, 3), NumPyRowMajorStorage()) == [1 2 3; 4 5 6]
    @test parse_reference_storage_convention("F") isa JuliaColumnMajorStorage
    @test parse_reference_storage_convention("numpy_row_major") isa NumPyRowMajorStorage
    @test_throws InvalidConfiguration parse_reference_storage_convention("weird")

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=2, binning=1)
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}()) == 16
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}("img_resolution" => 12)) == 12
end

@testset "Reference harness fixture" begin
    root = mktempdir()
    create_reference_fixture(root)
    bundle = load_reference_bundle(root)
    @test length(bundle.cases) == 6
    for case in bundle.cases
        result = validate_reference_case(case)
        @test size(result.actual) == size(result.expected)
        @test result.ok
    end
end

@testset "Reference bundle baselines" begin
    oopao_bundle = load_reference_bundle(default_reference_root())
    @test !isempty(reference_cases(oopao_bundle, :oopao))
    @test isempty(reference_cases(oopao_bundle, :specula))

    specula_root = default_specula_reference_root()
    if has_specula_reference_bundle(specula_root)
        specula_bundle = load_reference_bundle(specula_root)
        @test specula_bundle.metadata["baseline"] == "specula"
        @test !isempty(reference_cases(specula_bundle, :specula))
        @test isempty(reference_cases(specula_bundle, :oopao))
    else
        @info "Skipping SPECULA reference baseline inventory; no manifest found" root=specula_root
        @test true
    end
end

@testset "OOPAO reference regression" begin
    root = default_reference_root()
    if has_reference_bundle(root)
        bundle = load_reference_bundle(root)
        @test !isempty(bundle.cases)
        for case in bundle.cases
            result = validate_reference_case(case)
            @test size(result.actual) == size(result.expected)
            @test result.ok
        end
    else
        @info "Skipping OOPAO reference regression; no manifest found" root=root
        @test true
    end
end

@testset "SPECULA reference regression" begin
    root = default_specula_reference_root()
    if has_specula_reference_bundle(root)
        bundle = load_reference_bundle(root)
        cases = reference_cases(bundle, :specula)
        @test !isempty(cases)
        for case in cases
            result = validate_reference_case(case)
            @test size(result.actual) == size(result.expected)
            @test result.ok
        end
    else
        @info "Skipping SPECULA reference regression; no manifest found" root=root
        @test true
    end
end

@testset "Tutorial examples" begin
    image = run_tutorial_example("image_formation.jl")
    @test size(image.nominal_image) == size(image.aberrated_image)
    @test maximum(image.nominal_image) > maximum(image.aberrated_image)

    detector = run_tutorial_example("detector.jl")
    @test sum(detector.frame_native) ≈ sum(detector.photon_rate_image)
    @test size(detector.frame_sampled, 1) < size(detector.frame_native, 1)

    asterism = run_tutorial_example("asterism.jl")
    @test length(asterism.component_images) == 4
    @test sum(asterism.combined_image) >= sum(first(asterism.component_images))

    extended = run_tutorial_example("extended_source_sensing.jl")
    @test extended.n_samples == 25
    @test extended.sh_extended_rate ≈ extended.sh_point_rate rtol=1e-12
    @test extended.pyramid_extended_rate ≈
        extended.pyramid_point_rate rtol=1e-12
    @test extended.sh_extended_peak <= extended.sh_point_peak * (1 + 1e-12)
    # Direct WFS propagation does not yet apply each quadrature component's
    # angular offset to detector-plane morphology.
    @test_broken extended.sh_relative_morphology > 1e-6
    @test_broken extended.pyramid_relative_morphology > 1e-6

    sh_subaps = run_tutorial_example("shack_hartmann_subapertures.jl")
    @test sh_subaps.n_valid > 0
    @test sh_subaps.calibrated
    @test sh_subaps.metadata.n_valid_subap == sh_subaps.n_valid
    @test all(isfinite, sh_subaps.slopes)

    spatial = run_tutorial_example("spatial_filter.jl")
    @test size(spatial.filtered_phase) == size(spatial.filtered_amplitude)
    @test all(isfinite, spatial.filtered_phase)

    ncpa = run_tutorial_example("ncpa.jl")
    @test size(ncpa.image, 1) == 2 * size(ncpa.ncpa_opd, 1)
    @test size(ncpa.image, 2) == 2 * size(ncpa.ncpa_opd, 2)
    @test maximum(ncpa.image) > 0

    lift = run_tutorial_example("lift.jl")
    @test length(lift.coeffs_true) == length(lift.coeffs_fit)
    @test all(isfinite, lift.coeffs_fit)

    sprint = run_tutorial_example("sprint.jl")
    @test isfinite(sprint.estimate.shift_x)
    @test isfinite(sprint.estimate.shift_y)

    gsc = run_tutorial_example("gain_sensing_camera.jl")
    @test length(gsc.optical_gains) == 4
    @test all(isfinite, gsc.optical_gains)
    @test sum(gsc.calibration_frame) > 0
    @test sum(gsc.frame) > 0
    @test size(gsc.atmosphere_trace, 2) == 7
    @test all(isfinite, gsc.atmosphere_trace)

    transfer = run_tutorial_example("transfer_function.jl")
    @test size(transfer.rejection_db, 2) == length(transfer.loop_gains)
    @test all(isfinite, transfer.rejection_db)
    @test all(isfinite, transfer.closed_loop_db)

    tomography = run_tutorial_example("tomography.jl")
    @test all(isfinite, tomography.wavefront[.!isnan.(tomography.wavefront)])
    @test length(tomography.commands) == 4
    @test all(isfinite, tomography.commands)

    for name in ("closed_loop_shack_hartmann.jl", "closed_loop_pyramid.jl", "closed_loop_bioedge.jl", "closed_loop_zernike.jl")
        loop = run_tutorial_example(name)
        @test length(loop.residual_before) == length(loop.residual_after)
        @test all(isfinite, loop.residual_before)
        @test all(isfinite, loop.residual_after)
        if hasproperty(loop, :final_image)
            @test maximum(loop.final_image) > 0
        else
            @test maximum(loop.final_frame) > 0
            @test all(isfinite, loop.final_slopes)
            @test all(isfinite, loop.final_command)
        end
    end
end
