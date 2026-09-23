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
    @test fitting.fitting_matrix ≈ influence

    for T in (Float32, Float64)
        sampled = T[1 0; 0 1e-8; 0 0]
        for rtol in (zero(T), T(1e-6))
            prepared = @inferred TomographyFitting(
                sampled; regularization=rtol, resolution=3)
            @test prepared.influence_functions == sampled
            @test prepared.fitting_matrix ≈ pinv(sampled; rtol=rtol)
        end
        default_fitting = TomographyFitting(sampled; resolution=3)
        @test default_fitting.fitting_matrix ≈ pinv(sampled; rtol=T(1e-15))

        padded = zeros(T, 6, 4)
        sampled_view = @view padded[1:2:5, 1:2:3]
        copyto!(sampled_view, sampled)
        view_fitting = @inferred TomographyFitting(
            sampled_view; regularization=T(1e-6), resolution=3)
        @test view_fitting.influence_functions == sampled
        @test view_fitting.fitting_matrix ≈ pinv(sampled; rtol=T(1e-6))

        rank_deficient = T[1 0; 0 0; 0 0]
        rank_fitting = TomographyFitting(
            rank_deficient; regularization=zero(T), resolution=3)
        @test rank_fitting.fitting_matrix == pinv(rank_deficient; rtol=zero(T))

        zero_fitting = TomographyFitting(
            zeros(T, 3, 2); regularization=zero(T), resolution=3)
        @test all(iszero, zero_fitting.fitting_matrix)

        cutoff_matrix = T[1 0; 0 0.25; 0 0]
        cutoff_fitting = TomographyFitting(
            cutoff_matrix; regularization=T(0.25), resolution=3)
        @test cutoff_fitting.fitting_matrix ==
            pinv(cutoff_matrix; rtol=T(0.25))
    end
    @test_throws InvalidConfiguration TomographyFitting(
        influence; regularization=Inf, resolution=3)

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
    expected_matrix = (recon.operators.cox * transpose(imat)) /
        Matrix(imat * recon.operators.cxx * transpose(imat) .+ recon.operators.cnz)
    @test recon.reconstructor ≈ expected_matrix
    @test size(recon.reconstructor, 1) == count(grid_mask)
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

    # The measured interaction matrix may have a different precision from the
    # physical covariance model. The source builder promotes the cold solve.
    imat32 = Float32.(imat)
    recon_mixed = build_reconstructor(
        InteractionMatrixTomography(), imat32, grid_mask, atm, lgs, wfs,
        tomo, dm)
    @test eltype(recon_mixed.reconstructor) === Float64
    @test recon_mixed.operators.recstat ≈ aoc_covariance_reconstructor(
        Float64.(imat32), recon_mixed.operators.cxx,
        recon_mixed.operators.cox, recon_mixed.operators.cnz)
    recon_mixed_cpu = build_reconstructor(
        InteractionMatrixTomography(), imat32, grid_mask, atm, lgs, wfs,
        tomo, dm; build_backend=Calibration.CPUBuildBackend())
    @test eltype(recon_mixed_cpu.reconstructor) === Float64
    @test recon_mixed_cpu.operators.recstat ≈ recon_mixed.operators.recstat

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
    @test size(model.reconstructor, 1) == count(model.grid_mask)
    @test count(model.grid_mask) < length(model.grid_mask)

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
    masked_cross = cross_correlation(atm, lgs, wfs, tomo;
        grid_mask=model_recon.grid_mask)
    @test dropdims(masked_cross; dims=1) ≈ model_recon.operators.cxx
    @test model_recon.operators.cox ≈ model_recon.operators.cxx
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
    @test size(cmd_recon.matrix, 1) == count(dm.valid_actuators)
    original = copy(cmd_recon.matrix)
    mask_actuators!(cmd_recon, 1)
    @test all(iszero, @view cmd_recon.matrix[1, :])
    @test cmd_recon.matrix[2:end, :] == original[2:end, :]

end

@testset "Frozen S6 multi-source finite-height covariance adapter" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "fixtures",
        "aos_s6_von_karman_multisource.toml"))
    source = fixture["physical_source"]
    @test fixture["schema"] == "test.adaptive-optics-sim/von-karman-covariance/1"
    @test source["n_lgs"] == 2
    @test source["fit_src_height_m"] < Inf
    @test source["fractional_cn2"] == [0.4, 0.6]

    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=source["zenith_angle_deg"],
        layer_altitudes_m=source["layer_altitudes_m"],
        L0=source["L0_m"],
        r0_zenith=source["r0_zenith_m"],
        fractional_cn2=source["fractional_cn2"],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0, 90.0],
        wind_speed=[5.0, 10.0],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=source["lgs_radius_arcsec"],
        wavelength_m=589e-9,
        base_height_m=source["lgs_base_height_m"],
        n_lgs=source["n_lgs"],
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=source["pupil_diameter_m"],
        n_lenslets=source["n_lenslets"],
        n_px=4,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(2, 2),
        lenslet_grid_rotations_rad=source["lenslet_grid_rotations_rad"],
        lenslet_grid_offsets_fraction=reshape(
            source["lenslet_grid_offsets_fraction"], 2, source["n_lgs"]),
    )
    tomography = TomographyParams(
        n_fit_src=source["n_fit_src"],
        fov_optimization_arcsec=source["fov_optimization_arcsec"],
        fit_src_height_m=source["fit_src_height_m"],
    )
    grid_mask = reshape(Bool.(source["grid_mask"]),
        Tuple(Int.(source["grid_mask_shape"])))
    expected_cxx = reshape(fixture["cxx"]["values"],
        Tuple(Int.(fixture["cxx"]["shape"])))
    expected_cox = reshape(fixture["cox"]["values"],
        Tuple(Int.(fixture["cox"]["shape"])))

    cxx = @inferred auto_correlation(atmosphere, asterism, wfs, grid_mask)
    cox = @inferred cross_correlation(atmosphere, asterism, wfs, tomography;
        grid_mask=grid_mask)
    @test cxx ≈ expected_cxx rtol=2e-12 atol=2e-12
    @test cox ≈ expected_cox rtol=2e-12 atol=2e-12
    @test size(cox) == (4, 2, 4)

    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[2],
        valid_actuators=trues(2, 2),
    )
    model = build_reconstructor(ModelBasedTomography(), atmosphere, asterism,
        wfs, tomography, dm; build_backend=Calibration.CPUBuildBackend())
    model_cross = cross_correlation(atmosphere, asterism, wfs, tomography;
        grid_mask=model.grid_mask)
    @test model.operators.cxx ≈
        auto_correlation(atmosphere, asterism, wfs, model.grid_mask)
    @test model.operators.cox ≈
        dropdims(sum(model_cross; dims=1) ./ size(model_cross, 1); dims=1)
    @test size(model.operators.cox, 2) == size(model.operators.cxx, 1)

    # AOS accepts this nearly normalized profile and forwards its raw layer
    # strengths; AOC must not impose a narrower Float64 preparation boundary.
    near_atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=source["zenith_angle_deg"],
        layer_altitudes_m=source["layer_altitudes_m"],
        L0=source["L0_m"],
        r0_zenith=source["r0_zenith_m"],
        fractional_cn2=[0.4, 0.6000001],
        reference_wavelength_m=500e-9,
        wind_direction_deg=[0.0, 90.0],
        wind_speed=[5.0, 10.0],
    )
    @test auto_correlation(near_atmosphere, asterism, wfs, grid_mask)[1, 1] >
        cxx[1, 1]
    @test cross_correlation(near_atmosphere, asterism, wfs, tomography;
        grid_mask=grid_mask)[1, 1, 1] > cox[1, 1, 1]

    big_atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=big"0.0",
        layer_altitudes_m=BigFloat[0],
        L0=big"25.0",
        r0_zenith=big"0.2",
        fractional_cn2=BigFloat[1],
        reference_wavelength_m=big"5e-7",
        wind_direction_deg=BigFloat[0],
        wind_speed=BigFloat[10],
    )
    big_asterism = LGSAsterismParams(
        radius_arcsec=big"0.0",
        wavelength_m=big"5.89e-7",
        base_height_m=big"90000.0",
        n_lgs=1,
    )
    big_wfs = LGSWFSParams(
        pupil_diameter_m=big"8.0",
        n_lenslets=1,
        n_px=4,
        field_stop_size_arcsec=big"2.0",
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=BigFloat[0],
        lenslet_grid_offsets_fraction=zeros(BigFloat, 2, 1),
    )
    @test_throws UnsupportedAlgorithm auto_correlation(
        big_atmosphere, big_asterism, big_wfs, trues(1, 1))

    false_mask = falses(size(grid_mask))
    @test size(auto_correlation(atmosphere, asterism, wfs, false_mask)) == (0, 0)
    @test size(cross_correlation(atmosphere, asterism, wfs, tomography;
        grid_mask=false_mask)) == (4, 0, 0)

    zero_asterism = LGSAsterismParams(
        radius_arcsec=source["lgs_radius_arcsec"],
        wavelength_m=589e-9,
        base_height_m=source["lgs_base_height_m"],
        n_lgs=0,
    )
    zero_wfs = LGSWFSParams(
        pupil_diameter_m=source["pupil_diameter_m"],
        n_lenslets=source["n_lenslets"],
        n_px=4,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(2, 2),
    )
    @test size(auto_correlation(atmosphere, zero_asterism, zero_wfs, grid_mask)) ==
        (0, 0)
    @test size(cross_correlation(atmosphere, zero_asterism, zero_wfs, tomography;
        grid_mask=grid_mask)) == (4, 2, 0)

    atmosphere32 = TomographyAtmosphereParams(
        zenith_angle_deg=0f0,
        layer_altitudes_m=Float32[0],
        L0=25f0,
        r0_zenith=0.2f0,
        fractional_cn2=Float32[1],
        reference_wavelength_m=5f-7,
        wind_direction_deg=Float32[0],
        wind_speed=Float32[10],
    )
    asterism32 = LGSAsterismParams(
        radius_arcsec=0f0,
        wavelength_m=5.89f-7,
        base_height_m=90_000f0,
        n_lgs=1,
    )
    wfs32 = LGSWFSParams(
        pupil_diameter_m=8f0,
        n_lenslets=1,
        n_px=4,
        field_stop_size_arcsec=2f0,
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=Float32[0],
        lenslet_grid_offsets_fraction=zeros(Float32, 2, 1),
    )
    tomography32 = TomographyParams(
        n_fit_src=1,
        fov_optimization_arcsec=0f0,
        fit_src_height_m=Inf32,
    )
    @test eltype(@inferred auto_correlation(
        atmosphere32, asterism32, wfs32, trues(1, 1))) === Float32
    @test eltype(@inferred cross_correlation(
        atmosphere32, asterism32, wfs32, tomography32;
        grid_mask=trues(1, 1))) === Float32
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
    @test occursin("pre-adoption", fixture["provenance"]["characterization"])

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
    @test !isapprox(cox, cxx)
    @test fixture_matches(recstat,
        aoc_covariance_reconstructor(gamma, cxx, cox, cnz), "recstat")
    @test fixture_matches(r,
        fixture["scaling"]["reconstructor_scale_m"] .* recstat, "r")
    @test fixture_matches(ordered_r, r[:, [2, 1]], "ordered_r")
    @test fixture_matches(k, -(f * ordered_r) .* 2.0, "k")
    @test fixture_matches(wavefront, r * s_native, "wavefront")
    historical_map = fill(NaN, size(grid_mask))
    historical_map[grid_mask] .= wavefront
    @test fixture_matches(historical_map, wavefront_map, "wavefront_map";
        nans=true)
    @test fixture_matches(command, k * s_sim, "command")

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
    aligned = TOML.parsefile(joinpath(@__DIR__, "fixtures",
        "aos_s6_tomography_grid_aligned.toml"))
    @test aligned["schema"] == "test.adaptive-optics-sim/tomography-grid-aligned/1"
    @test aligned["source_fixture"] == "aos_s6_tomography_cpu.toml"
    @test aligned["array_backend"] == "CPU"
    @test aligned["float_type"] == "Float64"
    @test aligned["storage_order"] == "Julia column-major vec order"
    @test Tuple(aligned["grid_mask_shape"]) == size(model.grid_mask)
    @test aligned["active_linear_indices"] == findall(vec(model.grid_mask))
    @test aligned["guide_stars"] == lgs.n_lgs
    @test aligned["fit_sources"] == tomo.n_fit_src^2
    @test aligned["simulation_slopes"] == s_sim

    @test model.grid_mask == grid_mask
    @test fixture_matches(Matrix(model.operators.gamma), gamma, "gamma")
    @test fixture_matches(model.operators.cxx, cxx, "cxx")
    @test fixture_matches(Matrix(model.operators.cnz), cnz, "cnz")
    @test model.operators.cox ≈ model.operators.cxx
    @test model.operators.recstat ≈ aoc_covariance_reconstructor(
        model.operators.gamma, model.operators.cxx,
        model.operators.cox, model.operators.cnz)
    @test model.reconstructor ≈
        fixture["scaling"]["reconstructor_scale_m"] .* model.operators.recstat

    sampled_h = influence_functions(dm; resolution=size(model.grid_mask, 1))
    @test fixture_matches(sampled_h, h, "h")
    @test fixture_matches(command_reconstructor.fitting.fitting_matrix, f, "f")
    ordered_live = prepare_slope_order(SimulationSlopes(), model.reconstructor, 1)
    @test ordered_live ≈ model.reconstructor[:, [2, 1]]
    @test command_reconstructor.matrix ≈
        -(command_reconstructor.fitting.fitting_matrix * ordered_live) .* 2.0

    actual_wavefront = model.reconstructor * s_native
    actual_wavefront_map = fill(NaN, size(model.grid_mask))
    actual_wavefront_map[model.grid_mask] .= actual_wavefront
    actual_command = command_reconstructor.matrix * s_sim
    aligned_matrix = fixture_array(aligned["command_matrix"])
    aligned_command = fixture_vector(aligned["command"])
    @test aligned["command_matrix"]["unit"] ==
        "scaled actuator coordinate per simulation slope coordinate"
    @test aligned["command"]["unit"] == "scaled actuator coordinate"
    @test isapprox(command_reconstructor.matrix, aligned_matrix;
        rtol=aligned["rtol"], atol=aligned["atol"])
    @test isapprox(actual_command, aligned_command;
        rtol=aligned["rtol"], atol=aligned["atol"])
    @test actual_wavefront ≈ model.reconstructor * s_native
    @test actual_wavefront_map[grid_mask] ≈ actual_wavefront
    @test all(isnan, actual_wavefront_map[.!grid_mask])
    @test actual_command ≈ command_reconstructor.matrix * s_sim
    @test !isapprox(actual_command, command_reconstructor.matrix * s_native)
end
