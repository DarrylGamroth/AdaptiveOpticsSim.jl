using AdaptiveOpticsSim
using Random

abstract type SweepExecutionTarget end
struct CPUSweepTarget <: SweepExecutionTarget end
struct GPUSweepTarget{B<:GPUBackendTag} <: SweepExecutionTarget end

_sync_target!(::CPUSweepTarget, _) = nothing

function _sync_target!(::GPUSweepTarget, A)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(A))
    return nothing
end

_gpu_target_type(::CPUSweepTarget) = nothing
_gpu_target_type(::GPUSweepTarget{B}) where {B<:GPUBackendTag} = B

_sweep_policy(::CPUSweepTarget) = SplitGPUPrecision(Float32, Float32)
_sweep_policy(::GPUSweepTarget{B}) where {B<:GPUBackendTag} =
    AdaptiveOpticsSim.default_gpu_precision_policy(B)

_high_accuracy_policy(::CPUSweepTarget) = SplitGPUPrecision(Float32, Float64)
_high_accuracy_policy(::GPUSweepTarget{B}) where {B<:GPUBackendTag} =
    AdaptiveOpticsSim.high_accuracy_gpu_precision_policy(B)

_sweep_backend_array(::CPUSweepTarget) = Array

function _sweep_backend_array(::GPUSweepTarget{B}) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    return BackendArray
end

_build_backend(::CPUSweepTarget) = AdaptiveOpticsSim.CPUBuildBackend()
_build_backend(::GPUSweepTarget{B}) where {B<:GPUBackendTag} = AdaptiveOpticsSim.GPUArrayBuildBackend(B)

_backend_name(::CPUSweepTarget) = "cpu"
_backend_name(::GPUSweepTarget{B}) where {B<:GPUBackendTag} =
    string(something(AdaptiveOpticsSim.gpu_backend_name(B), B))

function _time_block_ns(f::F) where {F<:Function}
    t0 = time_ns()
    result = f()
    return time_ns() - t0, result
end

function _runtime_case(target::SweepExecutionTarget; resolution::Int, n_subap::Int, n_act::Int)
    policy = _sweep_policy(target)
    T = AdaptiveOpticsSim.gpu_runtime_type(policy)
    BackendArray = _sweep_backend_array(target)
    rng = MersenneTwister(1)
    tel = Telescope(resolution=resolution, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.3, T=T, backend=BackendArray)
    wfs = ShackHartmann(tel; n_subap=n_subap, mode=Diffractive(), T=T, backend=BackendArray)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    runtime = ClosedLoopRuntime(sim, recon; rng=rng)
    step!(runtime)
    _sync_target!(target, runtime.command)
    return runtime
end

function _timed_runtime_case(target::SweepExecutionTarget; label::AbstractString, resolution::Int, n_subap::Int, n_act::Int)
    runtime = _runtime_case(target; resolution, n_subap, n_act)
    stats = runtime_timing(() -> begin
        step!(runtime)
        _sync_target!(target, runtime.command)
    end; warmup=10, samples=50, gc_before=false)
    return (
        kind="runtime",
        label=String(label),
        backend=_backend_name(target),
        resolution=resolution,
        n_subap=n_subap,
        n_act=n_act,
        mean_ns=stats.mean_ns,
        p95_ns=stats.p95_ns,
    )
end

function _tomography_case_params(target::SweepExecutionTarget; n_lenslet::Int, n_lgs::Int, n_fit_src::Int, n_dm::Int)
    policy = _sweep_policy(target)
    high_accuracy = _high_accuracy_policy(target)
    TB = AdaptiveOpticsSim.gpu_build_type(policy)
    TH = AdaptiveOpticsSim.gpu_build_type(high_accuracy)
    build_backend = _build_backend(target)
    lgswfs = LGSWFSParams(
        diameter=TB(8.0),
        n_lenslet=n_lenslet,
        n_px=8,
        field_stop_size_arcsec=TB(2.0),
        valid_lenslet_map=trues(n_lenslet, n_lenslet),
        lenslet_rotation_rad=zeros(TB, n_lenslet^2),
        lenslet_offset=zeros(TB, 2, n_lenslet^2),
    )
    lgs = LGSAsterismParams(
        radius_arcsec=TB(7.6),
        wavelength=TB(589e-9),
        base_height_m=TB(90_000.0),
        n_lgs=n_lgs,
    )
    tomo = TomographyParams(
        n_fit_src=n_fit_src,
        fov_optimization_arcsec=TB(15.0),
        fit_src_height_m=TB(Inf),
    )
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
    grid_side = max(1, 2 * n_lenslet)
    pitch = TB(8.0 / max(1, n_lenslet))
    valid = trues(grid_side, grid_side)
    tdm = TomographyDMParams(
        heights_m=collect(TB, range(TB(0.0), length=n_dm, step=TB(6000.0))),
        pitch_m=fill(pitch, n_dm),
        cross_coupling=TB(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=valid,
    )
    noise_model = RelativeSignalNoise(TB(0.1))
    imat_rows = 2 * n_lenslet^2 * n_lgs
    grid_mask = trues(grid_side, grid_side)
    imat_cols = n_lgs * count(grid_mask)
    imat_t = reshape(TB.(range(TB(0.1), length=imat_rows * imat_cols, step=TB(0.01))), imat_rows, imat_cols)

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
        pitch_m=fill(TH(8.0 / max(1, n_lenslet)), n_dm),
        cross_coupling=TH(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=valid,
    )
    return (; build_backend, lgswfs, lgs, tomo, atm_tomo, tdm, noise_model, imat_t, grid_mask,
        lgswfs_hi, lgs_hi, tomo_hi, atm_tomo_hi, tdm_hi, high_accuracy)
end

function _timed_builder_case(target::SweepExecutionTarget; label::AbstractString, n_lenslet::Int, n_lgs::Int, n_fit_src::Int, n_dm::Int)
    p = _tomography_case_params(target; n_lenslet, n_lgs, n_fit_src, n_dm)

    modal_resolution = 4 * n_lenslet
    modal_act = max(4, n_lenslet + 2)
    runtime = _runtime_case(target; resolution=modal_resolution, n_subap=n_lenslet, n_act=modal_act)
    imat = interaction_matrix(runtime.dm, runtime.wfs, runtime.tel, runtime.src;
        amplitude=eltype(runtime.command)(0.05))
    ModalReconstructor(imat; build_backend=p.build_backend)
    modal_build_ns, modal_recon = _time_block_ns() do
        recon_local = ModalReconstructor(imat; build_backend=p.build_backend)
        _sync_target!(target, recon_local.reconstructor)
        recon_local
    end

    build_reconstructor(
        InteractionMatrixTomography(),
        p.imat_t,
        p.grid_mask,
        p.atm_tomo,
        p.lgs,
        p.lgswfs,
        p.tomo,
        p.tdm;
        noise_model=p.noise_model,
        build_backend=p.build_backend,
    )
    interaction_ns, interaction_recon = _time_block_ns() do
        recon_local = build_reconstructor(
            InteractionMatrixTomography(),
            p.imat_t,
            p.grid_mask,
            p.atm_tomo,
            p.lgs,
            p.lgswfs,
            p.tomo,
            p.tdm;
            noise_model=p.noise_model,
            build_backend=p.build_backend,
        )
        _sync_target!(target, recon_local.reconstructor)
        recon_local
    end

    build_reconstructor(
        ModelBasedTomography(),
        p.atm_tomo,
        p.lgs,
        p.lgswfs,
        p.tomo,
        p.tdm;
        noise_model=p.noise_model,
        build_backend=p.build_backend,
    )
    model_ns, model_recon = _time_block_ns() do
        recon_local = build_reconstructor(
            ModelBasedTomography(),
            p.atm_tomo,
            p.lgs,
            p.lgswfs,
            p.tomo,
            p.tdm;
            noise_model=p.noise_model,
            build_backend=p.build_backend,
        )
        _sync_target!(target, recon_local.reconstructor)
        recon_local
    end

    build_reconstructor(
        ModelBasedTomography(),
        p.atm_tomo_hi,
        p.lgs_hi,
        p.lgswfs_hi,
        p.tomo_hi,
        p.tdm_hi;
        noise_model=RelativeSignalNoise(AdaptiveOpticsSim.gpu_build_type(p.high_accuracy)(0.1)),
        build_backend=p.build_backend,
    )
    model_hi_ns, model_hi_recon = _time_block_ns() do
        recon_local = build_reconstructor(
            ModelBasedTomography(),
            p.atm_tomo_hi,
            p.lgs_hi,
            p.lgswfs_hi,
            p.tomo_hi,
            p.tdm_hi;
            noise_model=RelativeSignalNoise(AdaptiveOpticsSim.gpu_build_type(p.high_accuracy)(0.1)),
            build_backend=p.build_backend,
        )
        _sync_target!(target, recon_local.reconstructor)
        recon_local
    end

    return (
        kind="builder",
        label=String(label),
        backend=_backend_name(target),
        n_lenslet=n_lenslet,
        n_lgs=n_lgs,
        n_fit_src=n_fit_src,
        n_dm=n_dm,
        modal_build_ns=modal_build_ns,
        interaction_build_ns=interaction_ns,
        model_build_ns=model_ns,
        model_high_accuracy_build_ns=model_hi_ns,
        modal_type=typeof(modal_recon.reconstructor),
        interaction_type=typeof(interaction_recon.reconstructor),
        model_type=typeof(model_recon.reconstructor),
        model_hi_type=typeof(model_hi_recon.reconstructor),
    )
end

function _print_runtime_result(r)
    println("runtime_case label=$(r.label) backend=$(r.backend) resolution=$(r.resolution) n_subap=$(r.n_subap) n_act=$(r.n_act) mean_ns=$(r.mean_ns) p95_ns=$(r.p95_ns)")
end

function _print_builder_result(r)
    println("builder_case label=$(r.label) backend=$(r.backend) n_lenslet=$(r.n_lenslet) n_lgs=$(r.n_lgs) n_fit_src=$(r.n_fit_src) n_dm=$(r.n_dm) modal_build_ns=$(r.modal_build_ns) interaction_build_ns=$(r.interaction_build_ns) model_build_ns=$(r.model_build_ns) model_high_accuracy_build_ns=$(r.model_high_accuracy_build_ns)")
end

function run_backend_crossover_sweep(target::SweepExecutionTarget)
    runtime_cases = (
        (label="compact", resolution=16, n_subap=4, n_act=4),
        (label="small", resolution=32, n_subap=8, n_act=8),
        (label="medium", resolution=64, n_subap=16, n_act=16),
    )
    builder_cases = (
        (label="compact", n_lenslet=1, n_lgs=1, n_fit_src=1, n_dm=1),
        (label="small", n_lenslet=2, n_lgs=2, n_fit_src=2, n_dm=1),
        (label="medium", n_lenslet=3, n_lgs=2, n_fit_src=2, n_dm=2),
    )

    println("backend_crossover_sweep backend=$(_backend_name(target))")
    for case in runtime_cases
        _print_runtime_result(_timed_runtime_case(target; case...))
    end
    for case in builder_cases
        _print_builder_result(_timed_builder_case(target; case...))
    end
    return nothing
end
