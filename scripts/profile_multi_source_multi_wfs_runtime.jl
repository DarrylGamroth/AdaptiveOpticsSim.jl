using AdaptiveOpticsSim
using Random
using Statistics

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
        isdefined(Main, :CUDA) || error("profile_multi_source_multi_wfs_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_multi_source_multi_wfs_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_backend!(::Nothing, _) = nothing

function _sync_backend!(::Type{B}, array) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(array))
    return nothing
end

function _timed_stats!(f!::F; warmup::Int=5, samples::Int=20) where {F<:Function}
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
    return mean(timings), sort(timings)[clamp(round(Int, 0.95 * samples), 1, samples)]
end

function _build_recon(dm, wfs, tel, src)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=eltype(dm.state.coefs)(0.05))
    return ModalReconstructor(imat; gain=eltype(dm.state.coefs)(0.5))
end

function _seed_opd!(tel::Telescope{<:Any,<:Any}, ::Type{T}) where {T<:AbstractFloat}
    host = Matrix{T}(undef, tel.params.resolution, tel.params.resolution)
    @inbounds for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        host[i, j] = T(i + j / 10)
    end
    copyto!(tel.state.opd, host)
    return tel.state.opd
end

function run_profile(; backend_name::AbstractString="cpu", T::Type{<:AbstractFloat}=Float32)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    rng = MersenneTwister(7)

    tel = Telescope(resolution=32, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T)
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T)
    ast = Asterism([src1, src2])
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)

    _seed_opd!(tel, T)

    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
    sh_mean_ns, sh_p95_ns = _timed_stats!(() -> begin
        measure!(sh, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, sh.state.slopes)
    end)

    pyr = PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    pyr_mean_ns, pyr_p95_ns = _timed_stats!(() -> begin
        measure!(pyr, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, pyr.state.slopes)
    end)

    bio = BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    bio_mean_ns, bio_p95_ns = _timed_stats!(() -> begin
        measure!(bio, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, bio.state.slopes)
    end)

    tel1 = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    tel2 = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm1 = KolmogorovAtmosphere(tel1; r0=T(0.2), L0=T(25.0), T=T, backend=BackendArray)
    atm2 = KolmogorovAtmosphere(tel2; r0=T(0.2), L0=T(25.0), T=T, backend=BackendArray)
    dm1 = DeformableMirror(tel1; n_act=4, influence_width=T(0.3), T=T, backend=BackendArray)
    dm2 = DeformableMirror(tel2; n_act=4, influence_width=T(0.3), T=T, backend=BackendArray)
    sh1 = ShackHartmann(tel1; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
    sh2 = ShackHartmann(tel2; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
    pyr2 = PyramidWFS(tel2; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    det1 = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)
    det2 = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)
    sim_sh1 = AOSimulation(tel1, atm1, src, dm1, sh1)
    sim_sh2 = AOSimulation(tel2, atm2, src, dm2, sh2)
    sim_mixed2 = AOSimulation(tel2, atm2, src, dm2, pyr2)
    rt_sh1 = ClosedLoopRuntime(sim_sh1, _build_recon(dm1, sh1, tel1, src); rng=MersenneTwister(11), wfs_detector=det1)
    rt_sh2 = ClosedLoopRuntime(sim_sh2, _build_recon(dm2, sh2, tel2, src); rng=MersenneTwister(12), wfs_detector=det2)
    rt_mixed2 = ClosedLoopRuntime(sim_mixed2, _build_recon(dm2, pyr2, tel2, src); rng=MersenneTwister(13), wfs_detector=det2)
    step!(rt_sh1)
    step!(rt_sh2)
    step!(rt_mixed2)
    composite_compatible = CompositeSimulationInterface(SimulationInterface(rt_sh1), SimulationInterface(rt_sh2))
    composite_mixed = CompositeSimulationInterface(SimulationInterface(rt_sh1), SimulationInterface(rt_mixed2))

    comp_mean_ns, comp_p95_ns = _timed_stats!(() -> begin
        step!(composite_compatible)
        _sync_backend!(backend_tag, simulation_command(composite_compatible))
    end)
    mixed_mean_ns, mixed_p95_ns = _timed_stats!(() -> begin
        step!(composite_mixed)
        _sync_backend!(backend_tag, simulation_command(composite_mixed))
    end)

    println("multi_source_multi_wfs_profile")
    println("  backend: ", label)
    println("  sh_asterism_mean_ns: ", sh_mean_ns)
    println("  sh_asterism_p95_ns: ", sh_p95_ns)
    println("  pyramid_asterism_mean_ns: ", pyr_mean_ns)
    println("  pyramid_asterism_p95_ns: ", pyr_p95_ns)
    println("  bioedge_asterism_mean_ns: ", bio_mean_ns)
    println("  bioedge_asterism_p95_ns: ", bio_p95_ns)
    println("  composite_compatible_mean_ns: ", comp_mean_ns)
    println("  composite_compatible_p95_ns: ", comp_p95_ns)
    println("  composite_mixed_mean_ns: ", mixed_mean_ns)
    println("  composite_mixed_p95_ns: ", mixed_p95_ns)
    return nothing
end

run_profile(; backend_name=_backend_arg)
