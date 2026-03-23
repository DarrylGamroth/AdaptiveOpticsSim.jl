using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _branch_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "sequential"
const _replay_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "direct"

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
        isdefined(Main, :CUDA) || error("profile_pixel_output_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_pixel_output_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_pixel_output_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_pixel_output_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function _resolve_branch_mode(name::AbstractString)
    lowered = lowercase(name)
    lowered == "sequential" && return SequentialBranchExecution()
    lowered == "task" && return TaskParallelBranchExecution()
    lowered == "stream" && return BackendStreamBranchExecution()
    error("unsupported branch mode '$name'; use sequential, task, or stream")
end

function _resolve_replay_mode(name::AbstractString)
    lowered = lowercase(name)
    lowered == "direct" && return DirectReplayMode()
    lowered == "prepared" && return PreparedReplayMode()
    error("unsupported replay mode '$name'; use direct or prepared")
end

_sync_runtime!(::Nothing, runtime) = nothing

function _sync_runtime!(::Type{B}, runtime) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(runtime.command))
    return nothing
end

function run_profile(; backend_name::AbstractString="cpu", branch_name::AbstractString="sequential",
    replay_name::AbstractString="direct", samples::Int=20, warmup::Int=5)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    branch_mode = _resolve_branch_mode(branch_name)
    replay_mode = _resolve_replay_mode(replay_name)
    params = AO1883kSurrogateParams(branch_execution=branch_mode, replay_mode=replay_mode)

    t0 = time_ns()
    scenario = ao188_3k_surrogate(; params=params, backend=BackendArray, rng=MersenneTwister(1))
    _sync_runtime!(backend_tag, scenario)
    build_time_ns = time_ns() - t0

    step!(scenario)
    _sync_runtime!(backend_tag, scenario)
    timing = runtime_timing(() -> begin
        step!(scenario)
        _sync_runtime!(backend_tag, scenario)
    end; warmup=warmup, samples=samples, gc_before=false)
    phase = ao188_3k_phase_timing(scenario; warmup=warmup, samples=samples, gc_before=false)

    println("pixel_output_runtime_profile")
    println("  backend: ", label)
    println("  branch_mode: ", branch_name)
    println("  replay_mode: ", replay_name)
    println("  build_time_ns: ", build_time_ns)
    println("  runtime_step_mean_ns: ", timing.mean_ns)
    println("  runtime_step_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  high_sense_mean_ns: ", phase.high_sense_mean_ns)
    println("  low_sense_mean_ns: ", phase.low_sense_mean_ns)
    println("  high_reconstruct_mean_ns: ", phase.high_reconstruct_mean_ns)
    println("  low_reconstruct_mean_ns: ", phase.low_reconstruct_mean_ns)
    println("  delay_mean_ns: ", phase.delay_mean_ns)
    println("  apply_mean_ns: ", phase.apply_mean_ns)
    println("  total_phase_mean_ns: ", phase.total_mean_ns)
    println("  total_phase_p95_ns: ", phase.total_p95_ns)
    println("  wfs_frame_shape: ", size(scenario.high_wfs.state.spot_cube))
    return nothing
end

run_profile(; backend_name=_backend_arg, branch_name=_branch_arg, replay_name=_replay_arg)
