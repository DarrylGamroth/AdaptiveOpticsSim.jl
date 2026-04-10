using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _scale_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "compact"

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_zernike_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_zernike_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_zernike_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_zernike_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_interface!(::Nothing, _) = nothing

function _sync_interface!(::Type{B}, interface::SimulationInterface) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(interface.command))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(interface.slopes))
    frame = simulation_wfs_frame(interface)
    isnothing(frame) || AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(frame))
    return nothing
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    return lowered
end

function _zernike_scale_config(scale_name::AbstractString, T::Type{<:AbstractFloat})
    scale = _resolve_scale(scale_name)
    if scale == "compact"
        return (
            scale=scale,
            resolution=16,
            n_subap=4,
            n_act=3,
            diffraction_padding=2,
            interaction_amplitude=T(1e-8),
            warmup=5,
            samples=20,
        )
    elseif scale == "medium"
        return (
            scale=scale,
            resolution=64,
            n_subap=16,
            n_act=9,
            diffraction_padding=2,
            interaction_amplitude=T(5e-9),
            warmup=3,
            samples=12,
        )
    end
    return (
        scale=scale,
        resolution=128,
        n_subap=32,
        n_act=17,
        diffraction_padding=2,
        interaction_amplitude=T(2.5e-9),
        warmup=2,
        samples=6,
    )
end

function run_profile(; backend_name::AbstractString="cpu", scale_name::AbstractString="compact",
    samples::Union{Int,Nothing}=nothing, warmup::Union{Int,Nothing}=nothing)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    T = Float32
    rng = MersenneTwister(5)
    cfg = _zernike_scale_config(scale_name, T)
    resolved_samples = something(samples, cfg.samples)
    resolved_warmup = something(warmup, cfg.warmup)

    tel = Telescope(
        resolution=cfg.resolution,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        T=T,
        backend=BackendArray,
    )
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=cfg.n_act, influence_width=0.35, T=T, backend=BackendArray)
    wfs = ZernikeWFS(tel; n_subap=cfg.n_subap, diffraction_padding=cfg.diffraction_padding, T=T, backend=BackendArray)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=BackendArray)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=cfg.interaction_amplitude)
    recon = ModalReconstructor(imat; gain=T(0.4))
    runtime = ClosedLoopRuntime(sim, recon; rng=rng, wfs_detector=det)
    interface = simulation_interface(runtime)

    t0 = time_ns()
    prepare!(interface)
    _sync_interface!(backend_tag, interface)
    build_time_ns = time_ns() - t0

    step!(interface)
    _sync_interface!(backend_tag, interface)
    timing = runtime_timing(() -> begin
        step!(interface)
        _sync_interface!(backend_tag, interface)
    end; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)
    phase = runtime_phase_timing(interface; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)

    println("zernike_runtime_profile")
    println("  backend: ", backend_label)
    println("  scale: ", cfg.scale)
    println("  build_time_ns: ", build_time_ns)
    println("  runtime_step_mean_ns: ", timing.mean_ns)
    println("  runtime_step_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  sense_mean_ns: ", phase.sense_mean_ns)
    println("  reconstruct_mean_ns: ", phase.reconstruct_mean_ns)
    println("  apply_mean_ns: ", phase.apply_mean_ns)
    println("  snapshot_mean_ns: ", phase.snapshot_mean_ns)
    println("  total_phase_mean_ns: ", phase.total_mean_ns)
    println("  total_phase_p95_ns: ", phase.total_p95_ns)
    println("  pupil_resolution: ", cfg.resolution)
    println("  n_subap: ", cfg.n_subap)
    println("  dm_n_act: ", cfg.n_act)
    println("  wfs_frame_shape: ", size(simulation_wfs_frame(interface)))
    println("  slope_length: ", length(simulation_slopes(interface)))
    println("  command_length: ", length(simulation_command(interface)))
    return nothing
end

run_profile(; backend_name=_backend_arg, scale_name=_scale_arg)
