using AdaptiveOpticsSim
using Profile

function _sync_backend_array!(A)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(A))
    return nothing
end

function run_gpu_model_tomography_profile(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    policy = AdaptiveOpticsSim.default_gpu_precision_policy(B)
    high_accuracy = AdaptiveOpticsSim.high_accuracy_gpu_precision_policy(B)
    TB = AdaptiveOpticsSim.gpu_build_type(policy)
    TH = AdaptiveOpticsSim.gpu_build_type(high_accuracy)
    build_backend = AdaptiveOpticsSim.GPUArrayBuildBackend(B)

    n_lenslet = 3
    n_lgs = 2
    n_fit_src = 2
    n_dm = 2
    grid_side = 2 * n_lenslet

    atm_tomo = TomographyAtmosphereParams(
        zenith_angle_deg=TB(0.0),
        altitude_km=TB[0.0, 10.0],
        L0=TB(25.0),
        r0_zenith=TB(0.2),
        fractional_r0=TB[0.6, 0.4],
        wavelength=TB(500e-9),
        wind_direction_deg=TB[0.0, 45.0],
        wind_speed=TB[10.0, 20.0],
    )
    lgs = LGSAsterismParams(
        radius_arcsec=TB(7.6),
        wavelength=TB(589e-9),
        base_height_m=TB(90_000.0),
        n_lgs=n_lgs,
    )
    lgswfs = LGSWFSParams(
        diameter=TB(8.0),
        n_lenslet=n_lenslet,
        n_px=8,
        field_stop_size_arcsec=TB(2.0),
        valid_lenslet_map=trues(n_lenslet, n_lenslet),
        lenslet_rotation_rad=zeros(TB, n_lenslet^2),
        lenslet_offset=zeros(TB, 2, n_lenslet^2),
    )
    tomo = TomographyParams(
        n_fit_src=n_fit_src,
        fov_optimization_arcsec=TB(15.0),
        fit_src_height_m=TB(Inf),
    )
    tdm = TomographyDMParams(
        heights_m=collect(TB, range(TB(0.0), length=n_dm, step=TB(6000.0))),
        pitch_m=fill(TB(8.0 / n_lenslet), n_dm),
        cross_coupling=TB(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=trues(grid_side, grid_side),
    )
    noise_model = RelativeSignalNoise(TB(0.1))

    function build_model()
        recon = build_reconstructor(
            ModelBasedTomography(),
            atm_tomo,
            lgs,
            lgswfs,
            tomo,
            tdm;
            noise_model=noise_model,
            build_backend=build_backend,
        )
        _sync_backend_array!(recon.reconstructor)
        return recon
    end

    atm_tomo_hi = TomographyAtmosphereParams(
        zenith_angle_deg=TH(0.0),
        altitude_km=TH[0.0, 10.0],
        L0=TH(25.0),
        r0_zenith=TH(0.2),
        fractional_r0=TH[0.6, 0.4],
        wavelength=TH(500e-9),
        wind_direction_deg=TH[0.0, 45.0],
        wind_speed=TH[10.0, 20.0],
    )
    lgs_hi = LGSAsterismParams(
        radius_arcsec=TH(7.6),
        wavelength=TH(589e-9),
        base_height_m=TH(90_000.0),
        n_lgs=n_lgs,
    )
    lgswfs_hi = LGSWFSParams(
        diameter=TH(8.0),
        n_lenslet=n_lenslet,
        n_px=8,
        field_stop_size_arcsec=TH(2.0),
        valid_lenslet_map=trues(n_lenslet, n_lenslet),
        lenslet_rotation_rad=zeros(TH, n_lenslet^2),
        lenslet_offset=zeros(TH, 2, n_lenslet^2),
    )
    tomo_hi = TomographyParams(
        n_fit_src=n_fit_src,
        fov_optimization_arcsec=TH(15.0),
        fit_src_height_m=TH(Inf),
    )
    tdm_hi = TomographyDMParams(
        heights_m=collect(TH, range(TH(0.0), length=n_dm, step=TH(6000.0))),
        pitch_m=fill(TH(8.0 / n_lenslet), n_dm),
        cross_coupling=TH(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=trues(grid_side, grid_side),
    )

    function build_model_high_accuracy()
        recon = build_reconstructor(
            ModelBasedTomography(),
            atm_tomo_hi,
            lgs_hi,
            lgswfs_hi,
            tomo_hi,
            tdm_hi;
            noise_model=RelativeSignalNoise(TH(0.1)),
            build_backend=build_backend,
        )
        _sync_backend_array!(recon.reconstructor)
        return recon
    end

    build_model()
    build_model_high_accuracy()

    alloc_default = @allocated build_model()
    t0 = time_ns()
    recon = build_model()
    dt_default = time_ns() - t0
    default_type = typeof(recon.reconstructor)

    alloc_high = @allocated build_model_high_accuracy()
    t1 = time_ns()
    recon_hi = build_model_high_accuracy()
    dt_high = time_ns() - t1
    high_type = typeof(recon_hi.reconstructor)

    Profile.clear()
    Profile.init(n = 10^7, delay = 0.0005)
    Profile.@profile begin
        for _ in 1:3
            build_model()
        end
    end

    println("GPU model tomography profile")
    println("  backend: ", string(something(AdaptiveOpticsSim.gpu_backend_name(B), B)))
    println("  case: medium")
    println("  default_build_ns: ", dt_default)
    println("  default_build_alloc_bytes: ", alloc_default)
    println("  default_reconstructor_type: ", default_type)
    println("  high_accuracy_build_ns: ", dt_high)
    println("  high_accuracy_build_alloc_bytes: ", alloc_high)
    println("  high_accuracy_reconstructor_type: ", high_type)
    println("  profile_flat:")
    Profile.print(stdout; format=:flat, sortedby=:count, mincount=6, noisefloor=2.0)
    println("  profile_tree:")
    Profile.print(stdout; mincount=6, noisefloor=2.0)
    return nothing
end
