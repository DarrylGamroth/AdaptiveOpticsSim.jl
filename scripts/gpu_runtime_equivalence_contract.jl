using AdaptiveOpticsSim
using Random

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
using .SubaruAO188Simulation

const _RUNTIME_EQ_RTOL = 1f-4
const _RUNTIME_EQ_ATOL = 5f-5

struct StaticAtmosphere{A,B<:AdaptiveOpticsSim.AbstractArrayBackend} <: AdaptiveOpticsSim.AbstractAtmosphere
    screen::A
end

AdaptiveOpticsSim.backend(::StaticAtmosphere{<:Any,B}) where {B} = B()

function _deterministic_phase_screen(tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    n = tel.params.resolution
    host = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        x = T((i - 1) / max(n - 1, 1))
        y = T((j - 1) / max(n - 1, 1))
        host[i, j] = T(3e-8) * sin(T(2π) * x) +
                     T(2e-8) * cos(T(4π) * y) +
                     T(1e-8) * sin(T(2π) * (x + y))
    end
    host .*= Array(tel.state.pupil)
    return host
end

function StaticAtmosphere(tel::Telescope; T::Type{<:AbstractFloat}=Float32, backend::AdaptiveOpticsSim.AbstractArrayBackend=backend(tel))
    selector = AdaptiveOpticsSim.require_same_backend(tel, AdaptiveOpticsSim._resolve_backend_selector(backend))
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    host = _deterministic_phase_screen(tel, T)
    screen = array_backend{T}(undef, size(host)...)
    copyto!(screen, host)
    return StaticAtmosphere{typeof(screen),typeof(selector)}(screen)
end

AdaptiveOpticsSim.advance!(atm::StaticAtmosphere, tel::Telescope, rng::AbstractRNG) = atm
AdaptiveOpticsSim.advance!(atm::StaticAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng()) = atm

function AdaptiveOpticsSim.propagate!(atm::StaticAtmosphere, tel::Telescope)
    copyto!(tel.state.opd, atm.screen)
    return tel
end

function _gpu_backend_selector(::Type{AdaptiveOpticsSim.CUDABackendTag})
    return AdaptiveOpticsSim.CUDABackend()
end

function _gpu_backend_selector(::Type{AdaptiveOpticsSim.AMDGPUBackendTag})
    return AdaptiveOpticsSim.AMDGPUBackend()
end

function _ao188_noise_free_params(::Type{T}=Float32;
    branch_execution::AbstractExecutionPolicy=SequentialExecution()) where {T<:AbstractFloat}
    return AO188SimulationParams(
        T=T,
        source_magnitude=0.0,
        branch_execution=branch_execution,
        replay_mode=PreparedReplayMode(),
        latency=AO188LatencyModel(
            high_measurement_delay_frames=0,
            low_measurement_delay_frames=0,
            reconstruction_delay_frames=0,
            dm_delay_frames=0,
        ),
        high_detector=AO188WFSDetectorConfig(T=T, noise=NoiseNone()),
        low_detector=AO188WFSDetectorConfig(T=T, noise=NoiseNone()),
    )
end

function _set_deterministic_opd!(tel::Telescope)
    copyto!(tel.state.opd, _deterministic_phase_screen(tel, eltype(tel.state.opd)))
    return tel
end

function _evaluate_ao188!(surrogate::AO188Simulation)
    fill!(surrogate.dm.state.coefs, zero(eltype(surrogate.dm.state.coefs)))
    fill!(surrogate.high_command, zero(eltype(surrogate.high_command)))
    fill!(surrogate.low_command, zero(eltype(surrogate.low_command)))
    fill!(surrogate.combined_command, zero(eltype(surrogate.combined_command)))
    fill!(surrogate.command, zero(eltype(surrogate.command)))
    _set_deterministic_opd!(surrogate.tel)
    SubaruAO188Simulation._measure_branches!(surrogate.params.branch_execution, surrogate)
    reconstruct!(surrogate.high_command, surrogate.high_reconstructor, surrogate.high_wfs.state.slopes)
    reconstruct!(surrogate.low_command, surrogate.low_reconstructor, surrogate.low_wfs.state.slopes)
    surrogate.combined_command .= surrogate.high_command .+ surrogate.low_command
    copyto!(surrogate.command, surrogate.combined_command)
    return surrogate
end

function _max_relative_error(actual::AbstractArray, expected::AbstractArray)
    num = abs.(Array(actual) .- Array(expected))
    den = max.(abs.(Array(expected)), eps(real(eltype(Array(expected)))))
    return maximum(num ./ den)
end

function _assert_close(label::AbstractString, actual, expected; rtol=_RUNTIME_EQ_RTOL, atol=_RUNTIME_EQ_ATOL)
    actual_cpu = Array(actual)
    expected_cpu = Array(expected)
    ok = isapprox(actual_cpu, expected_cpu; rtol=rtol, atol=atol)
    max_abs = maximum(abs.(actual_cpu .- expected_cpu))
    max_rel = _max_relative_error(actual_cpu, expected_cpu)
    println("  ", label, ": max_abs=", max_abs, " max_rel=", max_rel)
    @assert ok "$label mismatch (max_abs=$max_abs, max_rel=$max_rel, rtol=$rtol, atol=$atol)"
end

function _assert_max_abs(label::AbstractString, actual, expected; atol::Real)
    actual_cpu = Array(actual)
    expected_cpu = Array(expected)
    max_abs = maximum(abs.(actual_cpu .- expected_cpu))
    max_rel = _max_relative_error(actual_cpu, expected_cpu)
    println("  ", label, ": max_abs=", max_abs, " max_rel=", max_rel)
    @assert max_abs <= atol "$label mismatch (max_abs=$max_abs, atol=$atol)"
end

function _run_ao188_equivalence(::Type{B}, branch_mode::AbstractExecutionPolicy) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    backend = _gpu_backend_selector(B)

    cpu = subaru_ao188_simulation(; params=_ao188_noise_free_params(), backend=CPUBackend(), rng=MersenneTwister(1))
    gpu = subaru_ao188_simulation(; params=_ao188_noise_free_params(Float32; branch_execution=branch_mode), backend=backend, rng=MersenneTwister(1))

    _evaluate_ao188!(cpu)
    _evaluate_ao188!(gpu)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu.command))

    println("ao188_runtime_equivalence")
    _assert_close("high_spot_cube", gpu.high_wfs.state.spot_cube, cpu.high_wfs.state.spot_cube)
    _assert_close("low_spot_cube", gpu.low_wfs.state.spot_cube, cpu.low_wfs.state.spot_cube)
    _assert_close("high_slopes", gpu.high_wfs.state.slopes, cpu.high_wfs.state.slopes)
    _assert_close("low_slopes", gpu.low_wfs.state.slopes, cpu.low_wfs.state.slopes)
    _assert_max_abs("command", gpu.command, cpu.command; atol=1f-3)
end

function _post_command_observation!(surrogate::AO188Simulation, host_command::AbstractVector)
    fill!(surrogate.dm.state.coefs, zero(eltype(surrogate.dm.state.coefs)))
    fill!(surrogate.high_command, zero(eltype(surrogate.high_command)))
    fill!(surrogate.low_command, zero(eltype(surrogate.low_command)))
    fill!(surrogate.combined_command, zero(eltype(surrogate.combined_command)))
    fill!(surrogate.command, zero(eltype(surrogate.command)))
    _set_deterministic_opd!(surrogate.tel)
    copyto!(surrogate.command, host_command)
    copyto!(surrogate.dm.state.coefs, host_command)
    apply!(surrogate.dm, surrogate.tel, DMAdditive())
    SubaruAO188Simulation._measure_branches!(SequentialExecution(), surrogate)
    return surrogate
end

function _run_ao188_post_command_equivalence(::Type{B}, branch_mode::AbstractExecutionPolicy;
    T::Type{<:AbstractFloat}=Float64) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    backend = _gpu_backend_selector(B)

    cmd_src = subaru_ao188_simulation(; params=_ao188_noise_free_params(T), backend=CPUBackend(), rng=MersenneTwister(1))
    _evaluate_ao188!(cmd_src)
    host_command = Array(cmd_src.command)

    cpu = subaru_ao188_simulation(; params=_ao188_noise_free_params(T), backend=CPUBackend(), rng=MersenneTwister(2))
    gpu = subaru_ao188_simulation(; params=_ao188_noise_free_params(T; branch_execution=branch_mode), backend=backend, rng=MersenneTwister(2))
    _post_command_observation!(cpu, host_command)
    _post_command_observation!(gpu, host_command)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu.command))

    println("ao188_post_command_equivalence T=", T)
    _assert_close("tel_opd", gpu.tel.state.opd, cpu.tel.state.opd; rtol=T(1e-8), atol=T(1e-12))
    _assert_close("post_high_spot_cube", gpu.high_wfs.state.spot_cube, cpu.high_wfs.state.spot_cube; rtol=T(1e-8), atol=T(1e-4))
    _assert_close("post_low_spot_cube", gpu.low_wfs.state.spot_cube, cpu.low_wfs.state.spot_cube; rtol=T(1e-8), atol=T(1e-4))
    _assert_close("post_high_slopes", gpu.high_wfs.state.slopes, cpu.high_wfs.state.slopes; rtol=T(1e-8), atol=T(1e-8))
    _assert_close("post_low_slopes", gpu.low_wfs.state.slopes, cpu.low_wfs.state.slopes; rtol=T(1e-8), atol=T(1e-8))
end

function _na_profile(T::Type{<:AbstractFloat})
    return T[
        89_500 90_000 90_500 91_000 91_500;
        0.10   0.25   0.30   0.22   0.13
    ]
end

function _build_lgs_case(backend, ::Type{T}, profile::Symbol) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=112,
        diameter=8.2,
        sampling_time=1e-3,
        central_obstruction=0.30,
        T=T,
        backend=backend,
    )
    src = if profile == :none
        LGSSource(
            magnitude=T(0),
            wavelength=T(589e-9),
            laser_coordinates=(T(5.0), T(0.0)),
            elongation_factor=T(1.6),
            T=T,
        )
    else
        LGSSource(
            magnitude=T(0),
            wavelength=T(589e-9),
            laser_coordinates=(T(5.0), T(0.0)),
            na_profile=_na_profile(T),
            fwhm_spot_up=T(1.0),
            T=T,
        )
    end
    wfs = ShackHartmann(tel; n_subap=14, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1e-3), qe=T(1), binning=1, T=T, backend=backend)
    _set_deterministic_opd!(tel)
    return tel, src, wfs, det
end

function _run_lgs_equivalence(::Type{B}, profile::Symbol) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    T = Float32

    tel_cpu, src_cpu, wfs_cpu, det_cpu = _build_lgs_case(CPUBackend(), T, profile)
    tel_gpu, src_gpu, wfs_gpu, det_gpu = _build_lgs_case(_gpu_backend_selector(B), T, profile)

    measure!(wfs_cpu, tel_cpu, src_cpu, det_cpu; rng=MersenneTwister(2))
    measure!(wfs_gpu, tel_gpu, src_gpu, det_gpu; rng=MersenneTwister(2))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_gpu.state.slopes))

    println("lgs_sh_equivalence profile=", profile)
    _assert_close("spot_cube", wfs_gpu.state.spot_cube, wfs_cpu.state.spot_cube)
    _assert_close("slopes", wfs_gpu.state.slopes, wfs_cpu.state.slopes)
end

function _build_mixed_sh_asterism_case(backend, ::Type{T}) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=112,
        diameter=8.2,
        sampling_time=1e-3,
        central_obstruction=0.30,
        T=T,
        backend=backend,
    )
    ngs = Source(wavelength=T(589e-9), magnitude=T(0), coordinates=(T(0.0), T(0.0)), T=T)
    lgs = LGSSource(
        magnitude=T(0),
        wavelength=T(589e-9),
        laser_coordinates=(T(5.0), T(0.0)),
        na_profile=_na_profile(T),
        fwhm_spot_up=T(1.0),
        T=T,
    )
    ast = Asterism([ngs, lgs])
    wfs = ShackHartmann(tel; n_subap=14, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1e-3), qe=T(1), binning=1, T=T, backend=backend)
    _set_deterministic_opd!(tel)
    return tel, ast, wfs, det
end

function _run_mixed_sh_asterism_equivalence(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    T = Float32

    tel_cpu, ast_cpu, wfs_cpu, det_cpu = _build_mixed_sh_asterism_case(CPUBackend(), T)
    tel_gpu, ast_gpu, wfs_gpu, det_gpu = _build_mixed_sh_asterism_case(_gpu_backend_selector(B), T)

    measure!(wfs_cpu, tel_cpu, ast_cpu, det_cpu; rng=MersenneTwister(3))
    measure!(wfs_gpu, tel_gpu, ast_gpu, det_gpu; rng=MersenneTwister(3))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_gpu.state.slopes))

    println("mixed_sh_asterism_equivalence")
    _assert_close("spot_cube", wfs_gpu.state.spot_cube, wfs_cpu.state.spot_cube)
    _assert_close("slopes", wfs_gpu.state.slopes, wfs_cpu.state.slopes)
end

function _build_zernike_case(backend, ::Type{T}) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=32,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        T=T,
        backend=backend,
    )
    src = Source(band=:I, magnitude=T(0), T=T)
    wfs = ZernikeWFS(tel; n_subap=8, diffraction_padding=2, T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1e-3), qe=T(1), binning=1, T=T, backend=backend)
    _set_deterministic_opd!(tel)
    return tel, src, wfs, det
end

function _run_zernike_equivalence(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    T = Float32

    tel_cpu, src_cpu, wfs_cpu, det_cpu = _build_zernike_case(CPUBackend(), T)
    tel_gpu, src_gpu, wfs_gpu, det_gpu = _build_zernike_case(_gpu_backend_selector(B), T)

    measure!(wfs_cpu, tel_cpu, src_cpu, det_cpu; rng=MersenneTwister(4))
    measure!(wfs_gpu, tel_gpu, src_gpu, det_gpu; rng=MersenneTwister(4))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_gpu.state.slopes))

    println("zernike_equivalence")
    _assert_close("camera_frame", wfs_gpu.state.camera_frame, wfs_cpu.state.camera_frame; rtol=1f-5, atol=8f0)
    _assert_close("reference_signal_2d", wfs_gpu.state.reference_signal_2d, wfs_cpu.state.reference_signal_2d; rtol=5f-5, atol=2f-6)
    _assert_close("slopes", wfs_gpu.state.slopes, wfs_cpu.state.slopes; rtol=5f-5, atol=2f-6)
    _assert_close("detector_frame", output_frame(det_gpu), output_frame(det_cpu); rtol=1f-5, atol=1f-2)
end

@inline _low_order_label(::Val{:tiptilt}) = :tiptilt
@inline _low_order_label(::Val{:steering}) = :steering
@inline _low_order_label(::Val{:focus}) = :focus

@inline _wfs_case_label(::Val{:sh}) = "sh"
@inline _wfs_case_label(::Val{:pyr}) = "pyr"
@inline _wfs_case_label(::Val{:bio}) = "bio"

function _build_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:sh}) where {T<:AbstractFloat}
    return ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
end

function _build_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:pyr}) where {T<:AbstractFloat}
    return PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:bio}) where {T<:AbstractFloat}
    return BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_multi_optic_hil_case(backend, ::Type{T}, ::Val{:tiptilt}, wfs_case::Val{W}=Val(:sh)) where {T<:AbstractFloat,W}
    tel = Telescope(
        resolution=16,
        diameter=T(8.0),
        sampling_time=T(1e-3),
        central_obstruction=T(0.0),
        T=T,
        backend=backend,
    )
    src = Source(band=:I, magnitude=T(0), T=T)
    atm = StaticAtmosphere(tel; T=T, backend=backend)
    tiptilt = TipTiltMirror(tel; scale=T(0.1), T=T, backend=backend, label=:tiptilt)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = _build_low_order_wfs(tel, backend, T, wfs_case)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=Symbol(:multi_optic_equivalence_, _wfs_case_label(wfs_case)), branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, AOSimulation(tel, src, atm, optic, wfs), NullReconstructor();
            wfs_detector=det,
            rng=MersenneTwister(91)),
    )
    prepare!(scenario)
    return scenario
end

function _build_multi_optic_hil_case(backend, ::Type{T}, ::Val{:steering}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=16,
        diameter=T(8.0),
        sampling_time=T(1e-3),
        central_obstruction=T(0.0),
        T=T,
        backend=backend,
    )
    src = Source(band=:I, magnitude=T(0), T=T)
    atm = StaticAtmosphere(tel; T=T, backend=backend)
    steering = LowOrderMirror(tel, (
        (x, y) -> T(0.1) * x,
        (x, y) -> T(0.1) * y,
    ); labels=:steering, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:steering => steering, :dm => dm)
    wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:steering_multi_optic_equivalence, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, AOSimulation(tel, src, atm, optic, wfs), NullReconstructor();
            wfs_detector=det,
            rng=MersenneTwister(91)),
    )
    prepare!(scenario)
    return scenario
end

function _build_multi_optic_hil_case(backend, ::Type{T}, ::Val{:focus}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=16,
        diameter=T(8.0),
        sampling_time=T(1e-3),
        central_obstruction=T(0.0),
        T=T,
        backend=backend,
    )
    src = Source(band=:I, magnitude=T(0), T=T)
    atm = StaticAtmosphere(tel; T=T, backend=backend)
    focus = FocusStage(tel; scale=T(0.1), T=T, backend=backend, label=:focus)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:focus => focus, :dm => dm)
    wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:focus_multi_optic_equivalence, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, AOSimulation(tel, src, atm, optic, wfs), NullReconstructor();
            wfs_detector=det,
            rng=MersenneTwister(91)),
    )
    prepare!(scenario)
    return scenario
end

function _multi_optic_case_commands(::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return fill(T(0.0125), 1), fill(T(0.025), 1), fill(T(0.02), 16)
end

function _multi_optic_case_commands(::Type{T}, ::Union{Val{:tiptilt},Val{:steering}}) where {T<:AbstractFloat}
    return fill(T(0.0125), 2), fill(T(0.025), 2), fill(T(0.02), 16)
end

function _multi_optic_case_tolerances(::Val{:focus}, ::Val{:sh}=Val(:sh))
    return (slopes_rtol=5f-3, slopes_atol=6f-3, frame_rtol=6f-3, frame_atol=1f6)
end

function _multi_optic_case_tolerances(::Union{Val{:tiptilt},Val{:steering}}, ::Union{Val{:sh},Val{:pyr}})
    return (slopes_rtol=3f-3, slopes_atol=4f-3, frame_rtol=4f-3, frame_atol=1f6)
end

function _multi_optic_case_tolerances(::Val{:tiptilt}, ::Val{:bio})
    return (slopes_rtol=1.5f-1, slopes_atol=4f-3, frame_rtol=6f-1, frame_atol=1f6)
end

function _run_multi_optic_hil_equivalence(::Type{B}, case::Val{K}, wfs_case::Val{W}=Val(:sh)) where {B<:AdaptiveOpticsSim.GPUBackendTag,K,W}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    T = Float32

    cpu = _build_multi_optic_hil_case(CPUBackend(), T, case, wfs_case)
    gpu = _build_multi_optic_hil_case(_gpu_backend_selector(B), T, case, wfs_case)
    gpu_repeat = _build_multi_optic_hil_case(_gpu_backend_selector(B), T, case, wfs_case)
    initial_low_order, updated_low_order, initial_dm = _multi_optic_case_commands(T, case)
    label = _low_order_label(case)
    tol = _multi_optic_case_tolerances(case, wfs_case)

    set_command!(cpu, NamedTuple{(label, :dm)}((initial_low_order, initial_dm)))
    set_command!(gpu, NamedTuple{(label, :dm)}((initial_low_order, initial_dm)))
    set_command!(gpu_repeat, NamedTuple{(label, :dm)}((initial_low_order, initial_dm)))
    sense!(cpu)
    sense!(gpu)
    sense!(gpu_repeat)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu_repeat)))

    println("$(K)_dm_", _wfs_case_label(wfs_case), "_hil_equivalence initial")
    _assert_close("command", command(gpu), command(cpu); rtol=1f-6, atol=1f-6)
    _assert_close("slopes", slopes(gpu), slopes(cpu); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    _assert_close("wfs_frame", wfs_frame(gpu), wfs_frame(cpu); rtol=tol.frame_rtol, atol=tol.frame_atol)
    _assert_close("gpu_repeat_command", command(gpu_repeat), command(gpu); rtol=1f-6, atol=1f-6)
    _assert_close("gpu_repeat_slopes", slopes(gpu_repeat), slopes(gpu); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    _assert_close("gpu_repeat_wfs_frame", wfs_frame(gpu_repeat), wfs_frame(gpu); rtol=tol.frame_rtol, atol=tol.frame_atol)

    update_command!(cpu, NamedTuple{(label,)}((updated_low_order,)))
    update_command!(gpu, NamedTuple{(label,)}((updated_low_order,)))
    update_command!(gpu_repeat, NamedTuple{(label,)}((updated_low_order,)))
    sense!(cpu)
    sense!(gpu)
    sense!(gpu_repeat)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu_repeat)))

    println("$(K)_dm_", _wfs_case_label(wfs_case), "_hil_equivalence updated")
    _assert_close("command", command(gpu), command(cpu); rtol=1f-6, atol=1f-6)
    _assert_close("slopes", slopes(gpu), slopes(cpu); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    _assert_close("wfs_frame", wfs_frame(gpu), wfs_frame(cpu); rtol=tol.frame_rtol, atol=tol.frame_atol)
    _assert_close("gpu_repeat_command", command(gpu_repeat), command(gpu); rtol=1f-6, atol=1f-6)
    _assert_close("gpu_repeat_slopes", slopes(gpu_repeat), slopes(gpu); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    _assert_close("gpu_repeat_wfs_frame", wfs_frame(gpu_repeat), wfs_frame(gpu); rtol=tol.frame_rtol, atol=tol.frame_atol)
end

function run_gpu_runtime_equivalence(::Type{B}; branch_mode::AbstractExecutionPolicy=SequentialExecution()) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    _run_ao188_equivalence(B, branch_mode)
    _run_lgs_equivalence(B, :none)
    _run_lgs_equivalence(B, :na)
    _run_mixed_sh_asterism_equivalence(B)
    _run_zernike_equivalence(B)
    _run_multi_optic_hil_equivalence(B, Val(:tiptilt), Val(:sh))
    _run_multi_optic_hil_equivalence(B, Val(:tiptilt), Val(:pyr))
    _run_multi_optic_hil_equivalence(B, Val(:tiptilt), Val(:bio))
    _run_multi_optic_hil_equivalence(B, Val(:steering), Val(:sh))
    _run_multi_optic_hil_equivalence(B, Val(:focus), Val(:sh))
    println("gpu_runtime_equivalence complete")
    return nothing
end

function run_gpu_runtime_equivalence_high_accuracy(::Type{B};
    branch_mode::AbstractExecutionPolicy=SequentialExecution()) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    _run_ao188_post_command_equivalence(B, branch_mode; T=Float64)
    println("gpu_runtime_equivalence_high_accuracy complete")
    return nothing
end
