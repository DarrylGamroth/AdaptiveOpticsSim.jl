@testset "Tomography Parameters and Geometry" begin
    atm = TomographyAtmosphereParams(
        zenith_angle_deg=30.0,
        altitude_km=[5.0, 10.0, 15.0],
        L0=30.0,
        r0_zenith=0.15,
        fractional_r0=[0.5, 0.3, 0.2],
        wavelength=500e-9,
        wind_direction_deg=[90.0, 45.0, 180.0],
        wind_speed=[10.0, 20.0, 15.0],
    )
    @test isapprox(airmass(atm), inv(cosd(30.0)))
    @test layer_altitude_m(atm) ≈ [5000.0, 10000.0, 15000.0] .* airmass(atm)
    vx, vy = wind_velocity_components(atm)
    @test vx ≈ [0.0, 20cosd(45.0), -15.0]
    @test vy ≈ [10.0, 20sind(45.0), 0.0]

    lgs = LGSAsterismParams(radius_arcsec=30.0, wavelength=589e-9, base_height_m=90_000.0, n_lgs=4)
    dirs = lgs_directions(lgs)
    @test size(dirs) == (4, 2)
    @test all(isapprox.(dirs[:, 1], fill(30.0 * π / (180 * 3600), 4)))
    vectors = direction_vectors(view(dirs, :, 1), view(dirs, :, 2))
    @test size(vectors) == (3, 4)
    @test all(vectors[3, :] .== 1.0)
    @test isapprox(lgs_height_m(lgs, atm), 90_000.0 * airmass(atm))

    tomo = TomographyParams(n_fit_src=3, fov_optimization_arcsec=4.0)
    zenith, azimuth = optimization_geometry(tomo)
    @test length(zenith) == 9
    @test length(azimuth) == 9
    @test maximum(zenith) > 0

    wfs = LGSWFSParams(
        diameter=8.2,
        n_lenslet=40,
        n_px=16,
        field_stop_size_arcsec=2.5,
        valid_lenslet_map=Bool[
            1 0 1
            0 1 0
            1 0 1
        ],
        lenslet_rotation_rad=zeros(4),
        lenslet_offset=zeros(2, 4),
    )
    @test n_valid_subapertures(wfs) == 5
    @test size(valid_lenslet_support(wfs)) == (7, 7)

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
        altitude_km=[0.0, 1.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_r0=[0.6, 0.3],
        wavelength=500e-9,
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
        altitude_km=[0.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_r0=[1.0],
        wavelength=500e-9,
        wind_direction_deg=[0.0],
        wind_speed=[10.0],
    )
    lgs = LGSAsterismParams(radius_arcsec=7.6, wavelength=589e-9, base_height_m=90_000.0, n_lgs=4)
    wfs = LGSWFSParams(
        diameter=8.0,
        n_lenslet=2,
        n_px=8,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=trues(2, 2),
        lenslet_rotation_rad=zeros(4),
        lenslet_offset=zeros(2, 4),
    )
    tomo = TomographyParams(n_fit_src=1, fov_optimization_arcsec=0.0)
    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[2],
        valid_actuators=trues(2, 2),
    )
    imat = [1.0 0.0 1.0; 0.0 1.0 1.0]
    grid_mask = trues(2, 2)
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
        rcond=0.0,
    )
    slopes = [1.0, 2.0]
    @test reconstruct_wavefront(recon, slopes) ≈ pinv(imat) * slopes
    out = zeros(3)
    reconstruct_wavefront!(out, recon, slopes)
    @test out ≈ pinv(imat) * slopes
    @test recon.fitting === fitting

    @test_throws UnsupportedAlgorithm build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm,
    )
end
