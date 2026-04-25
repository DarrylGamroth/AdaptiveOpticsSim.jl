using AdaptiveOpticsSim
using Random

_sync_cpu_array!(A) = nothing

function _time_block_ns(f::F) where {F<:Function}
    t0 = time_ns()
    result = f()
    return time_ns() - t0, result
end

function _warm_builder!(f::F) where {F<:Function}
    return f()
end

function run_cpu_sync_audit()
    policy = SplitGPUPrecision(Float32, Float32)
    high_accuracy = SplitGPUPrecision(Float32, Float64)
    T = gpu_runtime_type(policy)
    TB = gpu_build_type(policy)
    TH = gpu_build_type(high_accuracy)
    rng = MersenneTwister(1)
    build_backend = AdaptiveOpticsSim.CPUBuildBackend()

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=CPUBackend())
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3, T=T, backend=CPUBackend())
    wfs = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    runtime = ClosedLoopRuntime(sim, recon; rng=rng)
    step!(runtime)

    runtime_stats = runtime_timing(runtime; warmup=10, samples=100, gc_before=false)
    phase_stats = runtime_phase_timing(runtime; warmup=10, samples=100, gc_before=false)

    _warm_builder!() do
        ModalReconstructor(imat; build_backend=build_backend)
    end
    modal_build_ns, modal_recon = _time_block_ns() do
        ModalReconstructor(imat; build_backend=build_backend)
    end

    atm_tomo = TomographyAtmosphereParams(
        zenith_angle_deg=TB(0.0),
        altitude_km=TB[0.0],
        L0=TB(25.0),
        r0_zenith=TB(0.2),
        fractional_r0=TB[1.0],
        wavelength=TB(500e-9),
        wind_direction_deg=TB[0.0],
        wind_speed=TB[10.0],
    )
    lgs = LGSAsterismParams(radius_arcsec=TB(7.6), wavelength=TB(589e-9),
        base_height_m=TB(90_000.0), n_lgs=1)
    lgswfs = LGSWFSParams(diameter=TB(8.0), n_lenslet=1, n_px=8,
        field_stop_size_arcsec=TB(2.0), valid_lenslet_map=trues(1, 1),
        lenslet_rotation_rad=zeros(TB, 1), lenslet_offset=zeros(TB, 2, 1))
    tomo = TomographyParams(n_fit_src=1, fov_optimization_arcsec=TB(0.0),
        fit_src_height_m=TB(Inf))
    tdm = TomographyDMParams(heights_m=TB[0.0], pitch_m=TB[0.5],
        cross_coupling=TB(0.2), n_actuators=[1], valid_actuators=trues(1, 1))
    grid_mask = trues(1, 1)
    tomo_noise = RelativeSignalNoise(TB(0.1))
    imat_t = reshape(TB[1.0, 0.5], 2, 1)

    _warm_builder!() do
        build_reconstructor(
            InteractionMatrixTomography(),
            imat_t,
            grid_mask,
            atm_tomo,
            lgs,
            lgswfs,
            tomo,
            tdm;
            noise_model=tomo_noise,
            build_backend=build_backend,
        )
    end
    im_tomo_build_ns, im_tomo = _time_block_ns() do
        build_reconstructor(
            InteractionMatrixTomography(),
            imat_t,
            grid_mask,
            atm_tomo,
            lgs,
            lgswfs,
            tomo,
            tdm;
            noise_model=tomo_noise,
            build_backend=build_backend,
        )
    end

    _warm_builder!() do
        build_reconstructor(
            ModelBasedTomography(),
            atm_tomo,
            lgs,
            lgswfs,
            tomo,
            tdm;
            noise_model=tomo_noise,
            build_backend=build_backend,
        )
    end
    model_tomo_build_ns, model_tomo = _time_block_ns() do
        build_reconstructor(
            ModelBasedTomography(),
            atm_tomo,
            lgs,
            lgswfs,
            tomo,
            tdm;
            noise_model=tomo_noise,
            build_backend=build_backend,
        )
    end

    atm_tomo_hi = TomographyAtmosphereParams(
        zenith_angle_deg=TH(0.0),
        altitude_km=TH[0.0],
        L0=TH(25.0),
        r0_zenith=TH(0.2),
        fractional_r0=TH[1.0],
        wavelength=TH(500e-9),
        wind_direction_deg=TH[0.0],
        wind_speed=TH[10.0],
    )
    lgs_hi = LGSAsterismParams(radius_arcsec=TH(7.6), wavelength=TH(589e-9),
        base_height_m=TH(90_000.0), n_lgs=1)
    lgswfs_hi = LGSWFSParams(diameter=TH(8.0), n_lenslet=1, n_px=8,
        field_stop_size_arcsec=TH(2.0), valid_lenslet_map=trues(1, 1),
        lenslet_rotation_rad=zeros(TH, 1), lenslet_offset=zeros(TH, 2, 1))
    tomo_hi = TomographyParams(n_fit_src=1, fov_optimization_arcsec=TH(0.0),
        fit_src_height_m=TH(Inf))
    tdm_hi = TomographyDMParams(heights_m=TH[0.0], pitch_m=TH[0.5],
        cross_coupling=TH(0.2), n_actuators=[1], valid_actuators=trues(1, 1))

    _warm_builder!() do
        build_reconstructor(
            ModelBasedTomography(),
            atm_tomo_hi,
            lgs_hi,
            lgswfs_hi,
            tomo_hi,
            tdm_hi;
            noise_model=RelativeSignalNoise(TH(0.1)),
            build_backend=build_backend,
        )
    end
    model_tomo_high_ns, model_tomo_high = _time_block_ns() do
        build_reconstructor(
            ModelBasedTomography(),
            atm_tomo_hi,
            lgs_hi,
            lgswfs_hi,
            tomo_hi,
            tdm_hi;
            noise_model=RelativeSignalNoise(TH(0.1)),
            build_backend=build_backend,
        )
    end

    println("CPU sync audit")
    println("  runtime_precision_policy: ", typeof(policy))
    println("  runtime_eltype: ", T)
    println("  build_eltype: ", TB)
    println("  runtime_step_mean_ns: ", runtime_stats.mean_ns)
    println("  runtime_step_p95_ns: ", runtime_stats.p95_ns)
    println("  runtime_phase_sense_mean_ns: ", phase_stats.sense_mean_ns)
    println("  runtime_phase_reconstruct_mean_ns: ", phase_stats.reconstruct_mean_ns)
    println("  runtime_phase_apply_mean_ns: ", phase_stats.apply_mean_ns)
    println("  runtime_phase_total_p95_ns: ", phase_stats.total_p95_ns)
    println("  modal_build_ns: ", modal_build_ns)
    println("  interaction_tomography_build_ns: ", im_tomo_build_ns)
    println("  model_tomography_build_ns: ", model_tomo_build_ns)
    println("  model_tomography_high_accuracy_build_ns: ", model_tomo_high_ns)
    println("  runtime_command_type: ", typeof(runtime.command))
    println("  modal_reconstructor_type: ", typeof(modal_recon.reconstructor))
    println("  interaction_tomography_type: ", typeof(im_tomo.reconstructor))
    println("  model_tomography_type: ", typeof(model_tomo.reconstructor))
    println("  model_tomography_high_accuracy_type: ", typeof(model_tomo_high.reconstructor))
    return nothing
end

run_cpu_sync_audit()
