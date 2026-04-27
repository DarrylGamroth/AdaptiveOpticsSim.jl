using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using AdaptiveOpticsSim
using BenchmarkTools
using Random

abstract type BenchmarkExecutionTarget end
struct CPUBenchmarkTarget <: BenchmarkExecutionTarget end
struct GPUBenchmarkTarget{B<:GPUBackendTag} <: BenchmarkExecutionTarget end

_sync_target!(::CPUBenchmarkTarget, _) = nothing

function _sync_target!(::GPUBenchmarkTarget, A)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(A))
    return nothing
end

_gpu_target_type(::CPUBenchmarkTarget) = nothing
_gpu_target_type(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag} = B

_benchmark_policy(::CPUBenchmarkTarget) = SplitGPUPrecision(Float32, Float32)
_benchmark_policy(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag} =
    AdaptiveOpticsSim.default_gpu_precision_policy(B)

_high_accuracy_policy(::CPUBenchmarkTarget) = SplitGPUPrecision(Float32, Float64)
_high_accuracy_policy(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag} =
    AdaptiveOpticsSim.high_accuracy_gpu_precision_policy(B)

_benchmark_backend_array(::CPUBenchmarkTarget) = Array

function _benchmark_backend_array(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    return BackendArray
end

_build_backend(::CPUBenchmarkTarget) = CPUBuildBackend()
_build_backend(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag} = GPUArrayBuildBackend(B)

_backend_name(::CPUBenchmarkTarget) = "cpu"
_backend_name(::GPUBenchmarkTarget{B}) where {B<:GPUBackendTag} =
    string(something(AdaptiveOpticsSim.gpu_backend_name(B), B))

function _configure_benchmarks!()
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.0
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
    BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
    BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
    BenchmarkTools.DEFAULT_PARAMETERS.gcsample = false
    return nothing
end

function _runtime_case(target::BenchmarkExecutionTarget; resolution::Int, n_lenslets::Int, n_act::Int)
    policy = _benchmark_policy(target)
    T = AdaptiveOpticsSim.gpu_runtime_type(policy)
    BackendArray = _benchmark_backend_array(target)
    rng = runtime_rng(1)
    tel = Telescope(resolution=resolution, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.3, T=T, backend=BackendArray)
    wfs = ShackHartmann(tel; n_lenslets=n_lenslets, mode=Diffractive(), T=T, backend=BackendArray)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    runtime = ClosedLoopRuntime(sim, recon; rng=rng)
    step!(runtime)
    _sync_target!(target, runtime.command)
    return runtime
end

function _tomography_case_params(target::BenchmarkExecutionTarget; n_lenslet::Int, n_lgs::Int, n_fit_src::Int, n_dm::Int)
    policy = _benchmark_policy(target)
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

function _canonical_suite(target::BenchmarkExecutionTarget)
    _configure_benchmarks!()

    runtime = _runtime_case(target; resolution=16, n_lenslets=4, n_act=4)
    runtime_trial = run(@benchmarkable begin
        step!($runtime)
        _sync_target!($target, $runtime.command)
    end)

    p = _tomography_case_params(target; n_lenslet=3, n_lgs=2, n_fit_src=2, n_dm=2)
    builder_label, builder_trial, builder_high_accuracy_label, builder_high_accuracy_trial =
        _builder_benchmarks(target, p)

    return (
        backend=_backend_name(target),
        runtime_trial=runtime_trial,
        builder_label=builder_label,
        builder_trial=builder_trial,
        builder_high_accuracy_label=builder_high_accuracy_label,
        builder_high_accuracy_trial=builder_high_accuracy_trial,
    )
end

function _builder_benchmarks(target::BenchmarkExecutionTarget, p)
    model_builder_trial = run(@benchmarkable begin
        recon_local = build_reconstructor(
            ModelBasedTomography(),
            $(p.atm_tomo),
            $(p.lgs),
            $(p.lgswfs),
            $(p.tomo),
            $(p.tdm);
            noise_model=$(p.noise_model),
            build_backend=$(p.build_backend),
        )
        _sync_target!($target, recon_local.reconstructor)
    end)

    model_builder_high_accuracy_trial = run(@benchmarkable begin
        recon_local = build_reconstructor(
            ModelBasedTomography(),
            $(p.atm_tomo_hi),
            $(p.lgs_hi),
            $(p.lgswfs_hi),
            $(p.tomo_hi),
            $(p.tdm_hi);
            noise_model=RelativeSignalNoise($(AdaptiveOpticsSim.gpu_build_type(p.high_accuracy))(0.1)),
            build_backend=$(p.build_backend),
        )
        _sync_target!($target, recon_local.reconstructor)
    end)

    return (
        "model_tomography_build",
        model_builder_trial,
        "model_tomography_high_accuracy_build",
        model_builder_high_accuracy_trial,
    )
end

function _print_trial(label::AbstractString, trial::BenchmarkTools.Trial)
    println(label)
    println("  median_ns: ", BenchmarkTools.median(trial).time)
    println("  mean_ns: ", BenchmarkTools.mean(trial).time)
    println("  memory_bytes: ", BenchmarkTools.memory(trial))
    println("  allocs: ", BenchmarkTools.allocs(trial))
end

function _print_trial(label::Nothing, trial::Nothing)
    return nothing
end

function run_backend_benchmark_suite(target::BenchmarkExecutionTarget)
    suite = _canonical_suite(target)
    println("backend_benchmark_suite backend=", suite.backend)
    _print_trial("runtime_step", suite.runtime_trial)
    _print_trial(suite.builder_label, suite.builder_trial)
    _print_trial(suite.builder_high_accuracy_label, suite.builder_high_accuracy_trial)
    return suite
end

function _builder_benchmarks(::GPUBenchmarkTarget{AMDGPUBackendTag}, p)
    interaction_trial = run(@benchmarkable begin
        recon_local = build_reconstructor(
            InteractionMatrixTomography(),
            $(p.imat_t),
            $(p.grid_mask),
            $(p.atm_tomo),
            $(p.lgs),
            $(p.lgswfs),
            $(p.tomo),
            $(p.tdm);
            noise_model=$(p.noise_model),
            build_backend=$(p.build_backend),
        )
        _sync_target!(GPUBenchmarkTarget{AMDGPUBackendTag}(), recon_local.reconstructor)
    end)

    return (
        "interaction_tomography_build",
        interaction_trial,
        nothing,
        nothing,
    )
end
