function aoc_covariance_reconstructor(
    projection::AbstractMatrix{T},
    phase_covariance::AbstractMatrix{T},
    fit_phase_covariance::AbstractMatrix{T},
    measurement_noise_covariance::AbstractMatrix{T},
) where {T<:AbstractFloat}
    tomography = AdaptiveOpticsCalibration.Tomography
    specification = tomography.CovarianceReconstructorSpecification(
        size(projection, 2), size(projection, 1), size(fit_phase_covariance, 1), T)
    plan = AdaptiveOpticsCalibration.prepare(
        tomography.CovarianceReconstructor(), specification)
    inputs = tomography.CovarianceReconstructorInputs(
        projection, phase_covariance, fit_phase_covariance,
        measurement_noise_covariance)
    return tomography.reconstructor(AdaptiveOpticsCalibration.process(plan, inputs))
end

@testset "Tomography Parameters and Geometry" begin
    atm = TomographyAtmosphereParams(
        zenith_angle_deg=30.0,
        layer_altitudes_m=[5_000.0, 10_000.0, 15_000.0],
        L0=30.0,
        r0_zenith=0.15,
        fractional_cn2=[0.5, 0.3, 0.2],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[90.0, 45.0, 180.0],
        wind_speed=[10.0, 20.0, 15.0],
    )
    @test isapprox(zenith_angle_rad(atm), deg2rad(30.0))
    @test zenith_angle_deg(atm) ≈ 30.0
    @test isapprox(airmass(atm), inv(cosd(30.0)))
    @test atm.layer_altitudes_m == [5_000.0, 10_000.0, 15_000.0]
    @test atm.reference_wavelength_m == 500e-9
    @test layer_slant_ranges_m(atm) ≈ [5000.0, 10000.0, 15000.0] .* airmass(atm)
    @test wind_direction_rad(atm) ≈ deg2rad.([90.0, 45.0, 180.0])
    @test wind_direction_deg(atm) ≈ [90.0, 45.0, 180.0]
    vx, vy = wind_velocity_components(atm)
    @test vx ≈ [0.0, 20cosd(45.0), -15.0]
    @test vy ≈ [10.0, 20sind(45.0), 0.0]

    lgs = LGSAsterismParams(radius_arcsec=30.0, wavelength_m=589e-9, base_height_m=90_000.0, n_lgs=4)
    dirs = lgs_directions(lgs)
    @test size(dirs) == (4, 2)
    @test all(isapprox.(dirs[:, 1], fill(30.0 * π / (180 * 3600), 4)))
    vectors = direction_vectors(view(dirs, :, 1), view(dirs, :, 2))
    @test size(vectors) == (3, 4)
    @test all(vectors[3, :] .== 1.0)
    @test isapprox(lgs_height_m(lgs, atm), 90_000.0 * airmass(atm))
    @test lgs.wavelength_m == 589e-9

    tomo = TomographyParams(n_fit_src=3, fov_optimization_arcsec=4.0)
    zenith, azimuth = optimization_geometry(tomo)
    @test length(zenith) == 9
    @test length(azimuth) == 9
    @test maximum(zenith) > 0

    wfs = LGSWFSParams(
        pupil_diameter_m=8.2,
        n_lenslets=40,
        n_px=16,
        field_stop_size_arcsec=2.5,
        valid_lenslet_map=Bool[
            1 0 1
            0 1 0
            1 0 1
        ],
        lenslet_grid_rotations_rad=zeros(4),
        lenslet_grid_offsets_fraction=zeros(2, 4),
    )
    @test n_valid_subapertures(wfs) == 5
    @test size(valid_lenslet_support(wfs)) == (7, 7)
    @test lenslet_grid_support_diameter_m(wfs) ≈ 8.2 * 7 / 40
    @test wfs.n_lenslets == 40
    @test wfs.pupil_diameter_m == 8.2
    gamma, grid_mask = sparse_gradient_matrix(valid_lenslet_support(wfs))
    @test size(gamma, 1) == 2 * n_valid_subapertures(wfs)
    @test count(grid_mask) > 0

    registered_wfs = LGSWFSParams(
        pupil_diameter_m=8.0,
        n_lenslets=2,
        n_px=4,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(2, 2),
        lenslet_grid_rotations_rad=[0.0],
        lenslet_grid_offsets_fraction=reshape([0.125, -0.25], 2, 1),
    )
    registration_diameter_m = lenslet_grid_support_diameter_m(registered_wfs)
    registered_x, registered_y = AdaptiveOpticsSim.Tomography._guide_star_grid(
        1,
        registration_diameter_m,
        only(registered_wfs.lenslet_grid_rotations_rad),
        registered_wfs.lenslet_grid_offsets_fraction[1, 1],
        registered_wfs.lenslet_grid_offsets_fraction[2, 1],
    )
    @test only(registered_x) == -0.125 * registration_diameter_m
    @test only(registered_y) == 0.25 * registration_diameter_m
    @test_throws DimensionMismatchError AdaptiveOpticsSim.Tomography._active_guide_grid_params(
        zeros(2), zeros(2), zeros(2), 1)

    @test :altitude_km ∉ propertynames(atm)
    @test :wavelength ∉ propertynames(lgs)
    @test :diameter ∉ propertynames(wfs)
    @test :n_lenslet ∉ propertynames(wfs)
    @test :lenslet_rotation_rad ∉ propertynames(wfs)
    @test :lenslet_offset ∉ propertynames(wfs)

    dm = TomographyDMParams(
        heights_m=[0.0, 1000.0],
        pitch_m=[0.5, 0.5],
        cross_coupling=0.15,
        n_actuators=[20, 20],
        valid_actuators=Bool[
            1 0 1 0
            0 1 0 1
            1 0 1 0
            0 1 0 1
        ],
    )
    @test size(dm_valid_support(dm)) == (8, 8)

    @test_throws InvalidConfiguration TomographyParams(n_fit_src=2, fov_optimization_arcsec=0.0)
    @test_throws InvalidConfiguration TomographyAtmosphereParams(
        zenith_angle_deg=0.0,
        layer_altitudes_m=[0.0, 1_000.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_cn2=[0.6, 0.3],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0, 90.0],
        wind_speed=[5.0, 10.0],
    )
end

@testset "Tomography Fitting and Reconstruction" begin
    influence = Matrix{Float64}(I, 3, 3)
    fitting = TomographyFitting(influence; regularization=0.0, resolution=3)
    opd = [1.0, 2.0, 3.0]
    @test fit_commands(fitting, opd) ≈ opd

    atm = TomographyAtmosphereParams(
        zenith_angle_deg=0.0,
        layer_altitudes_m=[0.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_cn2=[1.0],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0],
        wind_speed=[10.0],
    )
    lgs = LGSAsterismParams(radius_arcsec=7.6, wavelength_m=589e-9, base_height_m=90_000.0, n_lgs=1)
    wfs = LGSWFSParams(
        pupil_diameter_m=8.0,
        n_lenslets=1,
        n_px=8,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=zeros(1),
        lenslet_grid_offsets_fraction=zeros(2, 1),
    )
    tomo = TomographyParams(n_fit_src=1, fov_optimization_arcsec=0.0)
    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[1],
        valid_actuators=trues(1, 1),
    )
    imat = reshape([1.0, 0.5], 2, 1)
    grid_mask = trues(1, 1)
    recon = build_reconstructor(
        InteractionMatrixTomography(),
        imat,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        fitting=fitting,
    )
    slopes = [1.0, 2.0]
    expected = (recon.operators.cox * transpose(imat)) / Matrix(imat * recon.operators.cxx * transpose(imat) .+ recon.operators.cnz) * slopes
    @test reconstruct_wavefront(recon, slopes) ≈ expected
    out = zeros(1)
    @test @inferred(reconstruct_wavefront!(out, recon, slopes)) === out
    @test out ≈ expected
    if !coverage_instrumented()
        @test @allocated(reconstruct_wavefront!(out, recon, slopes)) == 0
    end
    mapped = reconstruct_wavefront_map(recon, slopes)
    @test size(mapped) == (1, 1)
    @test mapped[1, 1] ≈ expected[1]
    @test recon.fitting === fitting
    @test recon.operators.cxx isa AbstractMatrix
    @test recon.operators.cox isa AbstractMatrix
    @test size(recon.operators.cnz, 1) == 2

    recon_noise = build_reconstructor(
        InteractionMatrixTomography(),
        imat,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        fitting=fitting,
        noise_model=ScalarMeasurementNoise(1e-2),
    )
    @test diag(recon_noise.operators.cnz) == fill(1e-2, 2)

    recon_cpu = build_reconstructor(
        InteractionMatrixTomography(),
        imat,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        fitting=fitting,
        build_backend=Calibration.CPUBuildBackend(),
    )
    @test recon_cpu.reconstructor isa Matrix
    @test recon_cpu.grid_mask isa Matrix{Bool}
    @test recon_cpu.operators.recstat ≈ aoc_covariance_reconstructor(
        imat, recon_cpu.operators.cxx, recon_cpu.operators.cox,
        recon_cpu.operators.cnz)
    cpu_system = imat * recon_cpu.operators.cxx * transpose(imat) .+
        recon_cpu.operators.cnz
    @test (@inferred AdaptiveOpticsSim.Tomography._tomographic_covariance_reconstructor(
        ScalarCPUStyle(), Calibration.CPUBuildBackend(), imat,
        recon_cpu.operators.cxx, recon_cpu.operators.cox,
        recon_cpu.operators.cnz, cpu_system)) ≈ recon_cpu.operators.recstat

    det = Detector(noise=NoiseReadout(0.2), qe=0.8, binning=2)
    detector_noise = PhotonReadoutSlopeNoise(det; photons_per_subaperture=1000.0, excess_noise=1.2)
    recon_detector_noise = build_reconstructor(
        InteractionMatrixTomography(),
        imat,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        fitting=fitting,
        noise_model=detector_noise,
    )
    @test all(diag(recon_detector_noise.operators.cnz) .> 0)

    model = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm,
    )
    @test size(model.operators.gamma, 1) == 2
    @test size(model.reconstructor, 2) == 2
    model_map = reconstruct_wavefront_map(model, [0.1, -0.2])
    @test size(model_map) == size(model.grid_mask)
    @test count(isnan, model_map) > 0

    model_cpu = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        build_backend=Calibration.CPUBuildBackend(),
    )
    @test model_cpu.reconstructor isa Matrix
    @test model_cpu.grid_mask isa Matrix{Bool}
    @test model_cpu.operators.recstat ≈ aoc_covariance_reconstructor(
        model_cpu.operators.gamma, model_cpu.operators.cxx,
        model_cpu.operators.cox, model_cpu.operators.cnz)
    @test model_cpu.reconstructor ≈
        (lgs.wavelength_m / 2) .* model_cpu.operators.recstat

    model_noise = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        noise_model=DiagonalMeasurementNoise([1e-2, 2e-2]),
    )
    @test diag(model_noise.operators.cnz) == [1e-2, 2e-2]

    native_mask = Bool[
        1 0
        1 1
    ]
    native_recon = TomographicReconstructor(
        InteractionMatrixTomography(),
        Matrix{Float64}(I, 3, 3),
        native_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm,
        nothing,
        nothing,
    )
    native_map = reconstruct_wavefront_map(native_recon, [1.0, 2.0, 3.0])
    @test native_map ≈ [
        1.0 NaN
        2.0 3.0
    ] nans=true

    cmds = dm_commands(recon, slopes)
    @test length(cmds) == size(recon.reconstructor, 1)
end

@testset "Tomography Command Assembly" begin
    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[2],
        valid_actuators=trues(2, 2),
    )
    modes = influence_functions(dm; resolution=5)
    @test size(modes) == (25, 4)

    mat = reshape(1.0:16.0, 2, 8)
    swapped = swap_xy_blocks(mat, 2; n_channels=2)
    @test swapped == mat[:, [3, 4, 1, 2, 7, 8, 5, 6]]
    interleaved = interleave_xy_columns(swapped, 2; n_channels=2)
    @test interleaved == mat[:, [3, 1, 4, 2, 7, 5, 8, 6]]

    atm = TomographyAtmosphereParams(
        zenith_angle_deg=0.0,
        layer_altitudes_m=[0.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_cn2=[1.0],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0],
        wind_speed=[10.0],
    )
    lgs = LGSAsterismParams(radius_arcsec=7.6, wavelength_m=589e-9, base_height_m=90_000.0, n_lgs=1)
    wfs = LGSWFSParams(
        pupil_diameter_m=8.0,
        n_lenslets=1,
        n_px=8,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=zeros(1),
        lenslet_grid_offsets_fraction=zeros(2, 1),
    )
    tomo = TomographyParams(n_fit_src=1, fov_optimization_arcsec=0.0)
    model_recon = build_reconstructor(ModelBasedTomography(), atm, lgs, wfs, tomo, dm)
    cmd_recon = assemble_reconstructor_and_fitting(
        model_recon,
        dm;
        n_channels=1,
        slope_order=SimulationSlopes(),
        scaling_factor=2.0,
    )
    @test size(cmd_recon.matrix, 2) == 2
    cmd_recon_cpu = assemble_reconstructor_and_fitting(
        model_recon,
        dm;
        build_backend=Calibration.CPUBuildBackend(),
    )
    @test cmd_recon_cpu.matrix isa Matrix
    commands = dm_commands(cmd_recon, [0.1, -0.2])
    @test length(commands) == count(dm.valid_actuators)
    command_out = similar(commands)
    command_input = [0.1, -0.2]
    @test @inferred(dm_commands!(command_out, cmd_recon, command_input)) ===
        command_out
    if !coverage_instrumented()
        @test @allocated(dm_commands!(
            command_out, cmd_recon, command_input)) == 0
    end
    original = copy(cmd_recon.matrix)
    mask_actuators!(cmd_recon, 1)
    @test all(iszero, @view cmd_recon.matrix[1, :])
    @test cmd_recon.matrix[2:end, :] == original[2:end, :]

    imat = reshape([1.0, 0.5], 2, 1)
    im_recon = build_reconstructor(
        InteractionMatrixTomography(),
        imat,
        trues(1, 1),
        atm,
        lgs,
        wfs,
        tomo,
        dm,
    )
    @test dm_commands(im_recon, [0.1, -0.2]) ≈ reconstruct_wavefront(im_recon, [0.1, -0.2])
end

@testset "Frozen S6 tomography CPU source characterization" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "fixtures",
        "aos_s6_tomography_cpu.toml"))
    @test fixture["schema"] == "test.adaptive-optics-sim/tomography-cpu/1"
    @test fixture["source_repository"] == "AdaptiveOpticsSim.jl"
    @test fixture["source_revision"] ==
        "429bf423f2bd19a9e893dbc767a5e86908ed84ea"
    @test fixture["source_paths"] == [
        "src/tomography/fitting.jl",
        "src/tomography/parameters.jl",
        "src/tomography/reconstructors.jl",
    ]
    @test fixture["source_test_path"] == "test/tomography.jl"
    @test fixture["array_backend"] == "CPU"
    @test fixture["float_type"] == "Float64"
    @test fixture["storage_order"] == "Julia column-major vec order"

    fixture_array(section) = reshape(Float64.(section["values"]),
        Tuple(Int.(section["shape"])))
    fixture_bool_array(section) = reshape(Bool.(section["values"]),
        Tuple(Int.(section["shape"])))
    fixture_vector(section) = Float64.(section["values"])
    tolerance(section) = (
        rtol=Float64(fixture["tolerances"][section]["rtol"]),
        atol=Float64(fixture["tolerances"][section]["atol"]),
    )
    fixture_matches(actual, expected, section; nans=false) =
        isapprox(actual, expected; tolerance(section)..., nans=nans)

    grid_mask = fixture_bool_array(fixture["grid_mask"])
    gamma = fixture_array(fixture["gamma"])
    cxx = fixture_array(fixture["cxx"])
    cox = fixture_array(fixture["cox"])
    cnz = fixture_array(fixture["cnz"])
    recstat = fixture_array(fixture["recstat"])
    r = fixture_array(fixture["r"])
    h = fixture_array(fixture["h"])
    f = fixture_array(fixture["f"])
    ordered_r = fixture_array(fixture["ordered_r"])
    k = fixture_array(fixture["k"])
    s_native = fixture_vector(fixture["slopes"]["native"])
    s_sim = fixture_vector(fixture["slopes"]["simulation"])
    wavefront = fixture_vector(fixture["wavefront"])
    wavefront_map = fixture_array(fixture["wavefront_map"])
    command = fixture_vector(fixture["command"])

    @test size(grid_mask) == (11, 11)
    @test count(grid_mask) == 9
    @test fixture["grid_mask"]["active_linear_indices"] == findall(vec(grid_mask))
    @test fixture["grid_mask"]["unit"] ==
        "dimensionless Boolean padded phase-sample support"
    @test fixture["gamma"]["unit"] ==
        "dimensionless finite-difference slope coefficient per phase-radian sample"
    @test fixture["cxx"]["unit"] == "phase-radian² covariance"
    @test fixture["cox"]["unit"] == "phase-radian² covariance"
    @test fixture["cnz"]["unit"] == "phase-radian² covariance"
    @test fixture["recstat"]["unit"] ==
        "phase-radian sample per source slope coordinate"
    @test fixture["r"]["unit"] ==
        "metres OPD per source slope coordinate"
    @test fixture["h"]["unit"] ==
        "source sampled influence-function coefficient per actuator coordinate"
    @test fixture["f"]["unit"] ==
        "actuator coordinate per source sampled influence-function coefficient"
    @test fixture["ordered_r"]["unit"] ==
        "metres OPD per simulation slope coordinate"
    @test fixture["k"]["unit"] ==
        "scaled actuator coordinate per simulation slope coordinate"
    @test fixture["wavefront"]["unit"] == "metres OPD"
    @test fixture["wavefront_map"]["unit"] == "metres OPD; NaN outside grid mask"
    @test fixture["command"]["unit"] == "scaled actuator coordinate"
    @test fixture["fitting"]["pinv_rtol"] == 1e-15
    @test fixture["command_assembly"]["scaling_factor"] == 2.0
    @test fixture["scaling"]["reconstructor_scale_m"] == 589e-9 / 2
    @test fixture["constructor"]["lgs"]["n_lgs"] == 1
    @test fixture["constructor"]["wfs"]["n_lenslets"] == 1
    @test fixture["constructor"]["dm"]["n_actuators"] == [2]
    @test fixture["constructor"]["dm"]["valid_actuators"] ==
        [true, true, true, true]
    @test fixture["constructor"]["noise"] == Dict(
        "model" => "RelativeSignalNoise",
        "fraction" => 0.1,
    )
    @test s_native == [0.1, -0.2]
    @test s_sim == [-0.2, 0.1]

    atm = TomographyAtmosphereParams(
        zenith_angle_deg=0.0,
        layer_altitudes_m=[0.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_cn2=[1.0],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0],
        wind_speed=[10.0],
    )
    lgs = LGSAsterismParams(
        radius_arcsec=7.6,
        wavelength_m=589e-9,
        base_height_m=90_000.0,
        n_lgs=1,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=8.0,
        n_lenslets=1,
        n_px=8,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=zeros(1),
        lenslet_grid_offsets_fraction=zeros(2, 1),
    )
    tomo = TomographyParams(n_fit_src=1, fov_optimization_arcsec=0.0)
    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[2],
        valid_actuators=trues(2, 2),
    )
    model = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        build_backend=Calibration.CPUBuildBackend(),
    )
    command_reconstructor = assemble_reconstructor_and_fitting(
        model,
        dm;
        n_channels=1,
        slope_order=SimulationSlopes(),
        scaling_factor=fixture["command_assembly"]["scaling_factor"],
        build_backend=Calibration.CPUBuildBackend(),
    )

    @test model.grid_mask == grid_mask
    @test fixture_matches(Matrix(model.operators.gamma), gamma, "gamma")
    @test fixture_matches(model.operators.cxx, cxx, "cxx")
    @test fixture_matches(model.operators.cox, cox, "cox")
    @test fixture_matches(Matrix(model.operators.cnz), cnz, "cnz")
    @test fixture_matches(model.operators.recstat, recstat, "recstat")
    @test fixture_matches(model.reconstructor, r, "r")
    @test fixture_matches(model.reconstructor,
        fixture["scaling"]["reconstructor_scale_m"] .* model.operators.recstat,
        "r")

    sampled_h = influence_functions(dm; resolution=size(model.grid_mask, 1))
    @test fixture_matches(sampled_h, h, "h")
    @test fixture_matches(command_reconstructor.fitting.fitting_matrix, f, "f")
    @test fixture_matches(
        prepare_slope_order(SimulationSlopes(), model.reconstructor, 1),
        ordered_r,
        "ordered_r",
    )
    @test fixture_matches(ordered_r, r[:, [2, 1]], "ordered_r")
    @test fixture_matches(command_reconstructor.matrix, k, "k")
    @test fixture_matches(command_reconstructor.matrix,
        -(command_reconstructor.fitting.fitting_matrix * ordered_r) .* 2.0,
        "k")

    actual_wavefront = reconstruct_wavefront(model, s_native)
    actual_wavefront_map = reconstruct_wavefront_map(model, s_native)
    actual_command = dm_commands(command_reconstructor, s_sim)
    @test fixture_matches(actual_wavefront, wavefront, "wavefront")
    @test fixture_matches(actual_wavefront_map, wavefront_map, "wavefront_map";
        nans=true)
    @test fixture_matches(actual_command, command, "command")
    @test fixture_matches(actual_wavefront, model.reconstructor * s_native,
        "wavefront")
    @test fixture_matches(actual_command, command_reconstructor.matrix * s_sim,
        "command")
    @test !isapprox(actual_command, command_reconstructor.matrix * s_native)
end
