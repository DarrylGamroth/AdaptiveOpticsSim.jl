using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _scale_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "medium"
const _response_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "default"
const _sampling_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "default"
const _correction_arg = length(ARGS) >= 5 ? lowercase(ARGS[5]) : "default"
const _thermal_arg = length(ARGS) >= 6 ? lowercase(ARGS[6]) : "default"

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao3k_simulation.jl"))
using .SubaruAO3kSimulation
using .SubaruAO3kSimulation.SubaruAO188Simulation: AO188WFSDetectorConfig, subaru_ao188_phase_timing

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_ao3k_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_ao3k_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_ao3k_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_ao3k_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function _sync_simulation!(::Nothing, simulation)
    return nothing
end

function _sync_simulation!(::Type{B}, simulation) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(simulation.command))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(simulation.high_wfs.state.slopes))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(simulation.low_wfs.state.slopes))
    high_frame = simulation_wfs_frame(simulation)[1]
    low_frame = simulation_wfs_frame(simulation)[2]
    isnothing(high_frame) || AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(high_frame))
    isnothing(low_frame) || AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(low_frame))
    return nothing
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    return lowered
end

function _resolve_response(name::AbstractString)
    lowered = lowercase(name)
    lowered == "default" && return nothing, "default"
    lowered == "null" && return NullFrameResponse(), "null"
    error("unsupported response mode '$name'; use default or null")
end

function _resolve_sampling(name::AbstractString, T::Type{<:AbstractFloat})
    lowered = lowercase(name)
    lowered == "default" && return HgCdTeAvalancheArraySensor(read_time=T(2.5e-4), sampling_mode=CorrelatedDoubleSampling(), T=T), "cds"
    lowered == "single" && return HgCdTeAvalancheArraySensor(read_time=T(2.5e-4), sampling_mode=SingleRead(), T=T), "single"
    lowered == "ndr4" && return HgCdTeAvalancheArraySensor(read_time=T(2.5e-4), sampling_mode=AveragedNonDestructiveReads(4), T=T), "ndr4"
    lowered == "cds" && return HgCdTeAvalancheArraySensor(read_time=T(2.5e-4), sampling_mode=CorrelatedDoubleSampling(), T=T), "cds"
    lowered == "fowler8" && return HgCdTeAvalancheArraySensor(read_time=T(2.5e-4), sampling_mode=FowlerSampling(8), T=T), "fowler8"
    error("unsupported sampling mode '$name'; use default, single, ndr4, cds, or fowler8")
end

function _resolve_correction(name::AbstractString)
    lowered = lowercase(name)
    lowered == "default" && return ReferencePixelCommonModeCorrection(4, 4), "reference_pixel"
    lowered == "none" && return NullFrameReadoutCorrection(), "none"
    lowered == "row" && return ReferenceRowCommonModeCorrection(4), "row"
    lowered == "column" && return ReferenceColumnCommonModeCorrection(4), "column"
    lowered == "output" && return ReferenceOutputCommonModeCorrection(32; edge_rows=4, edge_cols=4), "output"
    error("unsupported correction mode '$name'; use default, none, row, column, or output")
end

function _resolve_thermal(name::AbstractString, T::Type{<:AbstractFloat})
    lowered = lowercase(name)
    lowered == "default" && return FixedTemperature(
        temperature_K=80.0,
        dark_current_law=ArrheniusRateLaw(293.0, 4500.0),
        glow_rate_law=ArrheniusRateLaw(293.0, 3000.0),
        T=T), "fixed80"
    lowered == "none" && return nothing, "none"
    lowered == "fixed80" && return FixedTemperature(
        temperature_K=80.0,
        dark_current_law=ArrheniusRateLaw(293.0, 4500.0),
        glow_rate_law=ArrheniusRateLaw(293.0, 3000.0),
        T=T), "fixed80"
    lowered == "dynamic80" && return FirstOrderThermalModel(
        ambient_temperature_K=293.0,
        setpoint_temperature_K=80.0,
        initial_temperature_K=293.0,
        time_constant_s=0.05,
        min_temperature_K=70.0,
        max_temperature_K=300.0,
        dark_current_law=ArrheniusRateLaw(293.0, 4500.0),
        glow_rate_law=ArrheniusRateLaw(293.0, 3000.0),
        T=T), "dynamic80"
    error("unsupported thermal mode '$name'; use default, none, fixed80, or dynamic80")
end

function _ao3k_scale_params(scale_name::AbstractString)
    scale = _resolve_scale(scale_name)
    if scale == "compact"
        return (
            scale=scale,
            kwargs=(
                resolution=96,
                n_subap=16,
                n_control_modes=256,
                control_grid_side=20,
                source_magnitude=8.5,
            ),
            warmup=3,
            samples=12,
        )
    elseif scale == "medium"
        return (
            scale=scale,
            kwargs=(;
                control_grid_side=40,
            ),
            warmup=3,
            samples=10,
        )
    end
    return (
        scale=scale,
        kwargs=(
            resolution=240,
            n_subap=40,
            n_control_modes=2048,
            control_grid_side=56,
            source_magnitude=11.0,
        ),
        warmup=2,
        samples=5,
    )
end

function _allocated_bytes(f!::F; warmup::Int=2, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
end

function run_profile(; backend_name::AbstractString="cpu", scale_name::AbstractString="medium",
    response_name::AbstractString="default", sampling_name::AbstractString="default",
    correction_name::AbstractString="default", thermal_name::AbstractString="default",
    samples::Union{Int,Nothing}=nothing,
    warmup::Union{Int,Nothing}=nothing)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    response_model, response_label = _resolve_response(response_name)
    sensor, sampling_label = _resolve_sampling(sampling_name, Float32)
    correction_model, correction_label = _resolve_correction(correction_name)
    thermal_model, thermal_label = _resolve_thermal(thermal_name, Float32)
    cfg = _ao3k_scale_params(scale_name)
    resolved_samples = something(samples, cfg.samples)
    resolved_warmup = something(warmup, cfg.warmup)
    T = Float32
    params = AO3kSimulationParams(; cfg.kwargs..., high_detector=AO188WFSDetectorConfig(
        T=T,
        integration_time=1e-3,
        qe=0.9,
        psf_sampling=1,
        binning=1,
        gain=1.0,
        dark_current=0.02,
        noise=NoiseReadout(0.1),
        sensor=HgCdTeAvalancheArraySensor(
            glow_rate=T(0.02),
            read_time=sensor.read_time,
            sampling_mode=sensor.sampling_mode,
            avalanche_gain=sensor.avalanche_gain,
            excess_noise_factor=sensor.excess_noise_factor,
            T=T),
        response_model=response_model,
        correction_model=correction_model,
        thermal_model=thermal_model,
    ))

    t0 = time_ns()
    simulation = subaru_ao3k_simulation(; params=params, backend=BackendArray, rng=MersenneTwister(1))
    _sync_simulation!(backend_tag, simulation)
    build_time_ns = time_ns() - t0

    step!(simulation)
    _sync_simulation!(backend_tag, simulation)
    timing = runtime_timing(() -> begin
        step!(simulation)
        _sync_simulation!(backend_tag, simulation)
    end; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)
    runtime_alloc_bytes = _allocated_bytes(() -> begin
        step!(simulation)
        _sync_simulation!(backend_tag, simulation)
    end; warmup=resolved_warmup, gc_before=false)
    phase = subaru_ao188_phase_timing(simulation; warmup=resolved_warmup, samples=resolved_samples, gc_before=false)
    high_metadata, low_metadata = simulation_wfs_metadata(simulation)
    high_frame, low_frame = simulation_wfs_frame(simulation)
    high_slopes, low_slopes = simulation_slopes(simulation)
    high_detector = simulation.high_detector
    high_reference_frame = isnothing(high_detector) ? nothing : detector_reference_frame(high_detector)
    high_signal_frame = isnothing(high_detector) ? nothing : detector_signal_frame(high_detector)
    high_combined_frame = isnothing(high_detector) ? nothing : detector_combined_frame(high_detector)
    high_reference_cube = isnothing(high_detector) ? nothing : detector_reference_cube(high_detector)
    high_signal_cube = isnothing(high_detector) ? nothing : detector_signal_cube(high_detector)
    high_read_cube = isnothing(high_detector) ? nothing : detector_read_cube(high_detector)
    high_read_times = isnothing(high_detector) ? nothing : detector_read_times(high_detector)

    println("ao3k_runtime_profile")
    println("  backend: ", backend_label)
    println("  scale: ", cfg.scale)
    println("  response_mode: ", response_label)
    println("  sampling_mode: ", sampling_label)
    println("  correction_mode: ", correction_label)
    println("  thermal_mode: ", thermal_label)
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
    println("  n_subap: ", params.n_subap)
    println("  n_control_modes: ", params.n_control_modes)
    println("  high_frame_shape: ", isnothing(high_frame) ? nothing : size(high_frame))
    println("  low_frame_shape: ", isnothing(low_frame) ? nothing : size(low_frame))
    println("  high_slopes_length: ", length(high_slopes))
    println("  low_slopes_length: ", length(low_slopes))
    println("  command_length: ", length(simulation_command(simulation)))
    println("  high_response_family: ", isnothing(high_metadata) ? nothing : high_metadata.frame_response)
    println("  low_response_family: ", isnothing(low_metadata) ? nothing : low_metadata.frame_response)
    println("  high_sensor: ", isnothing(high_metadata) ? nothing : high_metadata.sensor)
    println("  high_thermal_model: ", isnothing(high_metadata) ? nothing : high_metadata.thermal_model)
    println("  high_detector_temperature_K: ", isnothing(high_metadata) ? nothing : high_metadata.detector_temperature_K)
    println("  high_cooling_setpoint_K: ", isnothing(high_metadata) ? nothing : high_metadata.cooling_setpoint_K)
    println("  high_ambient_temperature_K: ", isnothing(high_metadata) ? nothing : high_metadata.ambient_temperature_K)
    println("  high_sampling_mode: ", isnothing(high_metadata) ? nothing : high_metadata.sampling_mode)
    println("  high_readout_correction: ", isnothing(high_metadata) ? nothing : high_metadata.readout_correction)
    println("  high_reference_frame_shape: ", isnothing(high_reference_frame) ? nothing : size(high_reference_frame))
    println("  high_signal_frame_shape: ", isnothing(high_signal_frame) ? nothing : size(high_signal_frame))
    println("  high_combined_frame_shape: ", isnothing(high_combined_frame) ? nothing : size(high_combined_frame))
    println("  high_reference_cube_shape: ", isnothing(high_reference_cube) ? nothing : size(high_reference_cube))
    println("  high_signal_cube_shape: ", isnothing(high_signal_cube) ? nothing : size(high_signal_cube))
    println("  high_read_cube_shape: ", isnothing(high_read_cube) ? nothing : size(high_read_cube))
    println("  high_read_times_length: ", isnothing(high_read_times) ? nothing : length(high_read_times))
    return nothing
end

run_profile(; backend_name=_backend_arg, scale_name=_scale_arg, response_name=_response_arg,
    sampling_name=_sampling_arg, correction_name=_correction_arg, thermal_name=_thermal_arg)
