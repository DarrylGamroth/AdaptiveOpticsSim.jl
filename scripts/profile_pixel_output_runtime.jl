using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _branch_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "sequential"
const _replay_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "direct"
const _scale_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "medium"

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
using .SubaruAO188Simulation

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_pixel_output_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_pixel_output_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_pixel_output_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_pixel_output_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function _resolve_branch_mode(name::AbstractString)
    lowered = lowercase(name)
    lowered == "sequential" && return SequentialExecution()
    lowered == "task" && return ThreadedExecution()
    lowered == "stream" && return BackendStreamExecution()
    error("unsupported branch mode '$name'; use sequential, task, or stream")
end

function _resolve_replay_mode(name::AbstractString)
    lowered = lowercase(name)
    lowered == "direct" && return DirectReplayMode()
    lowered == "prepared" && return PreparedReplayMode()
    error("unsupported replay mode '$name'; use direct or prepared")
end

_sync_runtime!(::Nothing, runtime) = nothing

function _sync_runtime!(::Type{B}, runtime) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(runtime.command))
    return nothing
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    return lowered
end

function _allocated_bytes(f!::F; warmup::Int=2, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
end

function _pixel_scale_params(scale_name::AbstractString)
    scale = _resolve_scale(scale_name)
    if scale == "compact"
        return (
            scale=scale,
            kwargs=(
                resolution=64,
                n_act=32,
                n_active_actuators=768,
                n_control_modes=64,
                control_grid_side=10,
                high_order_samples=8,
                low_order_lenslets=2,
                n_low_order_modes=4,
                r0=0.18,
                source_magnitude=7.0,
            ),
            warmup=3,
            samples=12,
        )
    elseif scale == "medium"
        return (
            scale=scale,
            kwargs=NamedTuple(),
            warmup=5,
            samples=20,
        )
    end
    return (
        scale=scale,
        kwargs=(
            resolution=160,
            n_act=64,
            n_active_actuators=3228,
            n_control_modes=512,
            control_grid_side=32,
            high_order_samples=20,
            low_order_lenslets=4,
            n_low_order_modes=16,
            r0=0.14,
            source_magnitude=9.0,
        ),
        warmup=2,
        samples=6,
    )
end

function run_profile(; backend_name::AbstractString="cpu", branch_name::AbstractString="sequential",
    replay_name::AbstractString="direct", scale_name::AbstractString="medium",
    samples::Union{Int,Nothing}=nothing, warmup::Union{Int,Nothing}=nothing)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    branch_mode = _resolve_branch_mode(branch_name)
    replay_mode = _resolve_replay_mode(replay_name)
    cfg = _pixel_scale_params(scale_name)
    params = AO188SimulationParams(; cfg.kwargs..., branch_execution=branch_mode, replay_mode=replay_mode)
    resolved_samples = something(samples, cfg.samples)
    resolved_warmup = something(warmup, cfg.warmup)

    t0 = time_ns()
    scenario = subaru_ao188_simulation(; params=params, backend=BackendArray, rng=runtime_rng(1))
    _sync_runtime!(backend_tag, scenario)
    build_time_ns = time_ns() - t0

    step!(scenario)
    _sync_runtime!(backend_tag, scenario)
    timing = runtime_timing(() -> begin
        step!(scenario)
        _sync_runtime!(backend_tag, scenario)
    end; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)
    runtime_alloc_bytes = _allocated_bytes(() -> begin
        step!(scenario)
        _sync_runtime!(backend_tag, scenario)
    end; warmup=resolved_warmup, gc_before=false)
    phase = subaru_ao188_phase_timing(scenario; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)

    println("pixel_output_runtime_profile")
    println("  backend: ", label)
    println("  branch_mode: ", branch_name)
    println("  replay_mode: ", replay_name)
    println("  scale: ", cfg.scale)
    println("  build_time_ns: ", build_time_ns)
    println("  runtime_step_mean_ns: ", timing.mean_ns)
    println("  runtime_step_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  runtime_alloc_bytes: ", runtime_alloc_bytes)
    println("  high_sense_mean_ns: ", phase.high_sense_mean_ns)
    println("  low_sense_mean_ns: ", phase.low_sense_mean_ns)
    println("  high_reconstruct_mean_ns: ", phase.high_reconstruct_mean_ns)
    println("  low_reconstruct_mean_ns: ", phase.low_reconstruct_mean_ns)
    println("  delay_mean_ns: ", phase.delay_mean_ns)
    println("  apply_mean_ns: ", phase.apply_mean_ns)
    println("  total_phase_mean_ns: ", phase.total_mean_ns)
    println("  total_phase_p95_ns: ", phase.total_p95_ns)
    println("  pupil_resolution: ", params.resolution)
    println("  high_order_samples: ", params.high_order_samples)
    println("  low_order_lenslets: ", params.low_order_lenslets)
    println("  dm_n_act: ", params.n_act)
    println("  n_control_modes: ", params.n_control_modes)
    println("  wfs_frame_shape: ", size(scenario.high_wfs.state.spot_cube))
    return nothing
end

run_profile(; backend_name=_backend_arg, branch_name=_branch_arg, replay_name=_replay_arg, scale_name=_scale_arg)
