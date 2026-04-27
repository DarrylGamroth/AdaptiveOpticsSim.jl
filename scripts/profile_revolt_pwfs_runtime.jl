using AdaptiveOpticsSim
using Random
using Statistics: mean
using KernelAbstractions: @index, @kernel

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _model_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "pwfs_unmod"
const _response_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "default"
const _samples_arg = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4
const _warmup_arg = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 1

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

@kernel function _fill_dm_command_kernel!(coefs, phase, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(coefs)
        @inbounds coefs[i] = T(3e-8) * (sin(T(0.017) * T(i) + phase) + T(0.25) * cos(T(0.031) * T(i) - phase))
    end
end

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_revolt_pwfs_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_revolt_pwfs_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_revolt_pwfs_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_revolt_pwfs_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_backend!(::Nothing, _) = nothing

function _sync_backend!(::Type{B}, array) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(array))
    return nothing
end

function _timed_stats!(f!::F; warmup::Int=1, samples::Int=4) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    GC.gc()
    timings = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        f!()
        timings[i] = time_ns() - t0
    end
    sorted = sort(timings)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return mean(timings), sorted[p95_idx]
end

function _allocated_bytes(f!::F; warmup::Int=1, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
end

function _resolve_model(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "pwfs"
        return (
            model=:pwfs,
            label="pwfs_modulated",
            resolution=40,
            source_magnitude=5.752575f0,
            calibration_magnitude=3.252575f0,
            pupil_samples=20,
            modulation=5.0f0,
            modulation_points=32,
            light_ratio=0.1f0,
            n_pix_separation=4,
            n_pix_edge=2,
            psf_centering=true,
        )
    elseif lowered in ("pwfs_unmod", "pwfs-unmod", "unmod")
        return (
            model=:pwfs_unmod,
            label="pwfs_unmodulated",
            resolution=480,
            source_magnitude=3.252575f0,
            calibration_magnitude=3.252575f0,
            pupil_samples=20,
            modulation=0.0f0,
            modulation_points=1,
            light_ratio=0.1f0,
            n_pix_separation=4,
            n_pix_edge=2,
            psf_centering=true,
        )
    end
    error("unsupported REVOLT PWFS model '$name'; use pwfs or pwfs_unmod")
end

function _resolve_response(name::AbstractString)
    lowered = lowercase(name)
    lowered == "default" && return nothing, :default
    lowered == "null" && return NullFrameResponse(), :null
    error("unsupported response mode '$name'; use default or null")
end

@inline function _fill_dm_command!(::AdaptiveOpticsSim.ScalarCPUStyle, coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat}
    @inbounds for i in eachindex(coefs)
        coefs[i] = T(3e-8) * (sin(T(0.017) * T(i) + phase) + T(0.25) * cos(T(0.031) * T(i) - phase))
    end
    return coefs
end

@inline function _fill_dm_command!(style::AdaptiveOpticsSim.AcceleratorStyle, coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat}
    AdaptiveOpticsSim.launch_kernel!(style, _fill_dm_command_kernel!, coefs, phase, length(coefs); ndrange=length(coefs))
    return coefs
end

@inline _fill_dm_command!(coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat} =
    _fill_dm_command!(AdaptiveOpticsSim.execution_style(coefs), coefs, phase)

function _phase_step!(runtime, backend_tag; phase_index::Int=1)
    advance!(runtime.atm, runtime.tel; rng=runtime.rng)
    propagate!(runtime.atm, runtime.tel)
    _fill_dm_command!(runtime.dm.state.coefs, eltype(runtime.dm.state.coefs)(phase_index))
    apply!(runtime.dm, runtime.tel, DMAdditive())
    sense!(runtime)
    _sync_backend!(backend_tag, simulation_wfs_frame(runtime))
    return nothing
end

function _camera_specs(::Type{T}) where {T<:AbstractFloat}
    return (
        ixon=(
            resolution=128,
            bits=14,
            full_well=T(16000.0),
            qe=T(0.95),
            dark_current=T(20.0),
            readout_noise=T(1.0),
            sensor=EMCCDSensor(T=T),
        ),
        cred2=(
            resolution=480,
            bits=14,
            full_well=T(33643.0),
            qe=T(0.80),
            dark_current=T(334.0),
            readout_noise=T(25.0),
            sensor=InGaAsSensor(T=T),
        ),
        cs165cu=(
            resolution=1080,
            bits=10,
            full_well=T(11000.0),
            qe=T(0.65),
            dark_current=T(4.0),
            readout_noise=T(4.0),
            sensor=CMOSSensor(T=T),
        ),
    )
end

function _detector_from_spec(spec; integration_time, response_model, T::Type{<:AbstractFloat}, backend)
    return Detector(
        NoisePhotonReadout(spec.readout_noise);
        integration_time=T(integration_time),
        qe=spec.qe,
        gain=T(1.0),
        dark_current=spec.dark_current,
        bits=spec.bits,
        full_well=spec.full_well,
        sensor=spec.sensor,
        response_model=response_model,
        T=T,
        backend=backend,
    )
end

function run_profile(; backend_name::AbstractString="cpu", model_name::AbstractString="pwfs_unmod",
    response_name::AbstractString="default", samples::Int=4, warmup::Int=1)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    cfg = _resolve_model(model_name)
    response_model, response_label = _resolve_response(response_name)
    T = Float32
    cameras = _camera_specs(T)

    tel = Telescope(
        resolution=cfg.resolution,
        diameter=8.0,
        sampling_time=0.002,
        central_obstruction=0.0,
        T=T,
        backend=BackendArray,
    )
    src = Source(band=:R, magnitude=cfg.source_magnitude, T=T)
    calibration_source = Source(band=:R, magnitude=cfg.calibration_magnitude, T=T)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.15),
        L0=T(25.0),
        fractional_cn2=T[0.5, 0.3, 0.2],
        wind_speed=T[10.0, 15.0, 20.0],
        wind_direction=T[0.0, 45.0, 90.0],
        altitude=T[0.0, 6000.0, 12000.0],
        T=T,
        backend=BackendArray,
    )
    dm = DeformableMirror(tel; n_act=16, influence_width=T(0.2), T=T, backend=BackendArray)
    wfs = PyramidWFS(tel;
        pupil_samples=cfg.pupil_samples,
        threshold=T(0.1),
        modulation=cfg.modulation,
        modulation_points=cfg.modulation_points,
        light_ratio=cfg.light_ratio,
        n_pix_separation=cfg.n_pix_separation,
        n_pix_edge=cfg.n_pix_edge,
        psf_centering=cfg.psf_centering,
        mode=Diffractive(),
        T=T,
        backend=BackendArray,
    )
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)

    wfs_detector = _detector_from_spec(cameras.ixon;
        integration_time=tel.params.sampling_time,
        response_model=response_model,
        T=T,
        backend=BackendArray,
    )
    gain_detector = _detector_from_spec(cameras.cs165cu;
        integration_time=tel.params.sampling_time,
        response_model=response_model,
        T=T,
        backend=BackendArray,
    )
    science_detector = _detector_from_spec(cameras.cred2;
        integration_time=tel.params.sampling_time,
        response_model=response_model,
        T=T,
        backend=BackendArray,
    )

    t0 = time_ns()
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    runtime = ClosedLoopRuntime(sim, recon; rng=runtime_rng(0), wfs_detector=wfs_detector)
    prepare!(runtime)
    _phase_step!(runtime, backend_tag; phase_index=1)
    build_time_ns = time_ns() - t0

    phase_index = Ref(1)
    atmosphere_mean_ns, atmosphere_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        advance!(runtime.atm, runtime.tel; rng=runtime.rng)
        propagate!(runtime.atm, runtime.tel)
        _sync_backend!(backend_tag, runtime.tel.state.opd)
    end; warmup=warmup, samples=samples)
    dm_apply_mean_ns, dm_apply_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        advance!(runtime.atm, runtime.tel; rng=runtime.rng)
        propagate!(runtime.atm, runtime.tel)
        _fill_dm_command!(runtime.dm.state.coefs, T(phase_index[]))
        apply!(runtime.dm, runtime.tel, DMAdditive())
        _sync_backend!(backend_tag, runtime.tel.state.opd)
    end; warmup=warmup, samples=samples)
    sense_mean_ns, sense_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup, samples=samples)

    atmosphere_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        advance!(runtime.atm, runtime.tel; rng=runtime.rng)
        propagate!(runtime.atm, runtime.tel)
        _sync_backend!(backend_tag, runtime.tel.state.opd)
    end; warmup=warmup)
    dm_apply_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        advance!(runtime.atm, runtime.tel; rng=runtime.rng)
        propagate!(runtime.atm, runtime.tel)
        _fill_dm_command!(runtime.dm.state.coefs, T(phase_index[]))
        apply!(runtime.dm, runtime.tel, DMAdditive())
        _sync_backend!(backend_tag, runtime.tel.state.opd)
    end; warmup=warmup)
    sense_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup)
    total_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup)

    wfs_metadata = detector_export_metadata(wfs_detector)
    gain_metadata = detector_export_metadata(gain_detector)
    science_metadata = detector_export_metadata(science_detector)
    wfs_frame = simulation_wfs_frame(runtime)
    wfs_reference_frame = detector_reference_frame(wfs_detector)
    wfs_signal_frame = detector_signal_frame(wfs_detector)
    wfs_combined_frame = detector_combined_frame(wfs_detector)
    wfs_reference_cube = detector_reference_cube(wfs_detector)
    wfs_signal_cube = detector_signal_cube(wfs_detector)
    wfs_read_cube = detector_read_cube(wfs_detector)

    println("revolt_pwfs_runtime_profile")
    println("  backend: ", backend_label)
    println("  model: ", cfg.label)
    println("  response_mode: ", response_label)
    println("  wfs_family: pyramid")
    println("  wfs_detector: ixon")
    println("  wfs_detector_sensor: ", wfs_metadata.sensor)
    println("  wfs_detector_nominal_resolution: ", cameras.ixon.resolution)
    println("  wfs_frame_response: ", wfs_metadata.frame_response)
    println("  science_detector: cred2")
    println("  science_detector_sensor: ", science_metadata.sensor)
    println("  science_detector_nominal_resolution: ", cameras.cred2.resolution)
    println("  gain_detector: cs165cu")
    println("  gain_detector_sensor: ", gain_metadata.sensor)
    println("  gain_detector_nominal_resolution: ", cameras.cs165cu.resolution)
    println("  pupil_resolution: ", cfg.resolution)
    println("  pupil_samples: ", cfg.pupil_samples)
    println("  dm_grid_shape: ", (dm.params.n_act, dm.params.n_act))
    println("  dm_command_length: ", length(dm.state.coefs))
    println("  slope_length: ", length(wfs.state.slopes))
    println("  source_band: ", src.params.band)
    println("  effective_sky_magnitude: ", src.params.magnitude)
    println("  effective_calibration_magnitude: ", calibration_source.params.magnitude)
    println("  wfs_frame_shape: ", size(wfs_frame))
    println("  detector_output_shape: ", size(output_frame(wfs_detector)))
    println("  wfs_reference_frame_shape: ", isnothing(wfs_reference_frame) ? nothing : size(wfs_reference_frame))
    println("  wfs_signal_frame_shape: ", isnothing(wfs_signal_frame) ? nothing : size(wfs_signal_frame))
    println("  wfs_combined_frame_shape: ", isnothing(wfs_combined_frame) ? nothing : size(wfs_combined_frame))
    println("  wfs_reference_cube_shape: ", isnothing(wfs_reference_cube) ? nothing : size(wfs_reference_cube))
    println("  wfs_signal_cube_shape: ", isnothing(wfs_signal_cube) ? nothing : size(wfs_signal_cube))
    println("  wfs_read_cube_shape: ", isnothing(wfs_read_cube) ? nothing : size(wfs_read_cube))
    println("  build_time_ns: ", build_time_ns)
    println("  atmosphere_mean_ns: ", atmosphere_mean_ns)
    println("  atmosphere_p95_ns: ", atmosphere_p95_ns)
    println("  dm_apply_mean_ns: ", dm_apply_mean_ns)
    println("  dm_apply_p95_ns: ", dm_apply_p95_ns)
    println("  sense_mean_ns: ", sense_mean_ns)
    println("  sense_p95_ns: ", sense_p95_ns)
    println("  total_mean_ns: ", total_mean_ns)
    println("  total_p95_ns: ", total_p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / total_mean_ns)
    println("  atmosphere_alloc_bytes: ", atmosphere_alloc_bytes)
    println("  dm_apply_alloc_bytes: ", dm_apply_alloc_bytes)
    println("  sense_alloc_bytes: ", sense_alloc_bytes)
    println("  total_alloc_bytes: ", total_alloc_bytes)
    return nothing
end

run_profile(; backend_name=_backend_arg, model_name=_model_arg, response_name=_response_arg,
    samples=_samples_arg, warmup=_warmup_arg)
