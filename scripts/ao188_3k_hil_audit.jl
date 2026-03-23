using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])

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
        isdefined(Main, :CUDA) || error("ao188_3k_hil_audit.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("ao188_3k_hil_audit.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("ao188_3k_hil_audit.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("ao188_3k_hil_audit.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_runtime!(::Nothing, runtime) = nothing

function _sync_runtime!(::Type{B}, runtime) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(runtime.command))
    return nothing
end

function _build_scenario(params, BackendArray, rng, backend_tag)
    t0 = time_ns()
    scenario = ao188_3k_surrogate(; params=params, backend=BackendArray, rng=rng)
    _sync_runtime!(backend_tag, scenario)
    build_time_ns = time_ns() - t0
    return scenario, build_time_ns
end

function run_ao188_3k_hil_audit(; backend_name::AbstractString="cpu", samples::Int=20, warmup::Int=5)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    params = AO1883kSurrogateParams()
    scenario, build_time_ns = _build_scenario(params, BackendArray, MersenneTwister(1), backend_tag)
    step!(scenario)
    _sync_runtime!(backend_tag, scenario)
    timing = runtime_timing(() -> begin
        step!(scenario)
        _sync_runtime!(backend_tag, scenario)
    end; warmup=warmup, samples=samples, gc_before=false)
    phase = ao188_3k_phase_timing(scenario; warmup=warmup, samples=samples, gc_before=false)
    frame_rate_hz = 1.0e9 / timing.mean_ns

    println("AO188/3k surrogate HIL audit")
    println("  backend: ", label)
    println("  wfs_proxy: diffractive_shack_hartmann")
    println("  diameter_m: ", params.diameter)
    println("  resolution: ", params.resolution)
    println("  dm_grid: ", params.n_act, "x", params.n_act)
    println("  active_actuators: ", count(scenario.active_mask))
    println("  control_modes: ", params.n_control_modes)
    println("  low_order_modes: ", params.n_low_order_modes)
    println("  n_subap_high: ", params.n_subap)
    println("  n_subap_low: ", params.n_low_order_subap)
    println("  slopes_high: ", length(scenario.high_wfs.state.slopes))
    println("  slopes_low: ", length(scenario.low_wfs.state.slopes))
    println("  latency_high_measurement_frames: ", params.latency.high_measurement_delay_frames)
    println("  latency_low_measurement_frames: ", params.latency.low_measurement_delay_frames)
    println("  latency_reconstruction_frames: ", params.latency.reconstruction_delay_frames)
    println("  latency_dm_frames: ", params.latency.dm_delay_frames)
    println("  high_detector_noise_type: ", typeof(scenario.high_detector.noise))
    println("  low_detector_noise_type: ", typeof(scenario.low_detector.noise))
    println("  build_time_ns: ", build_time_ns)
    println("  runtime_step_mean_ns: ", timing.mean_ns)
    println("  runtime_step_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", frame_rate_hz)
    println("  high_sense_mean_ns: ", phase.high_sense_mean_ns)
    println("  low_sense_mean_ns: ", phase.low_sense_mean_ns)
    println("  high_reconstruct_mean_ns: ", phase.high_reconstruct_mean_ns)
    println("  low_reconstruct_mean_ns: ", phase.low_reconstruct_mean_ns)
    println("  delay_mean_ns: ", phase.delay_mean_ns)
    println("  apply_mean_ns: ", phase.apply_mean_ns)
    println("  total_phase_mean_ns: ", phase.total_mean_ns)
    println("  total_phase_p95_ns: ", phase.total_p95_ns)
    println("  command_type: ", typeof(scenario.command))
    println("  high_modal_reconstructor_type: ", typeof(scenario.high_reconstructor.reconstructor))
    println("  high_command_basis_type: ", typeof(scenario.high_reconstructor.command_basis))
    println("  low_modal_reconstructor_type: ", typeof(scenario.low_reconstructor.reconstructor))
    println("  low_command_basis_type: ", typeof(scenario.low_reconstructor.command_basis))
    return nothing
end

run_ao188_3k_hil_audit(; backend_name=isempty(ARGS) ? "cpu" : ARGS[1])
