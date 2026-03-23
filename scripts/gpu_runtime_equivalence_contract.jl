using AdaptiveOpticsSim
using Random

const _RUNTIME_EQ_RTOL = 1f-4
const _RUNTIME_EQ_ATOL = 5f-5

function _ao188_noise_free_params(::Type{T}=Float32;
    branch_execution::AO188BranchExecutionMode=SequentialBranchExecution()) where {T<:AbstractFloat}
    return AO1883kSurrogateParams(
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
    n = tel.params.resolution
    T = eltype(tel.state.opd)
    host = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        x = T((i - 1) / max(n - 1, 1))
        y = T((j - 1) / max(n - 1, 1))
        host[i, j] = T(3e-8) * sin(T(2π) * x) +
                     T(2e-8) * cos(T(4π) * y) +
                     T(1e-8) * sin(T(2π) * (x + y))
    end
    host .*= Array(tel.state.pupil)
    copyto!(tel.state.opd, host)
    return tel
end

function _evaluate_ao188!(surrogate::AO1883kSurrogate)
    fill!(surrogate.dm.state.coefs, zero(eltype(surrogate.dm.state.coefs)))
    fill!(surrogate.high_command, zero(eltype(surrogate.high_command)))
    fill!(surrogate.low_command, zero(eltype(surrogate.low_command)))
    fill!(surrogate.combined_command, zero(eltype(surrogate.combined_command)))
    fill!(surrogate.command, zero(eltype(surrogate.command)))
    _set_deterministic_opd!(surrogate.tel)
    AdaptiveOpticsSim._measure_branches!(surrogate.params.branch_execution, surrogate)
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

function _run_ao188_equivalence(::Type{B}, branch_mode::AO188BranchExecutionMode) where {B<:GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    cpu = ao188_3k_surrogate(; params=_ao188_noise_free_params(), backend=Array, rng=MersenneTwister(1))
    gpu = ao188_3k_surrogate(; params=_ao188_noise_free_params(Float32; branch_execution=branch_mode), backend=BackendArray, rng=MersenneTwister(1))

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

function _post_command_observation!(surrogate::AO1883kSurrogate, host_command::AbstractVector)
    fill!(surrogate.dm.state.coefs, zero(eltype(surrogate.dm.state.coefs)))
    fill!(surrogate.high_command, zero(eltype(surrogate.high_command)))
    fill!(surrogate.low_command, zero(eltype(surrogate.low_command)))
    fill!(surrogate.combined_command, zero(eltype(surrogate.combined_command)))
    fill!(surrogate.command, zero(eltype(surrogate.command)))
    _set_deterministic_opd!(surrogate.tel)
    copyto!(surrogate.command, host_command)
    copyto!(surrogate.dm.state.coefs, host_command)
    apply!(surrogate.dm, surrogate.tel, DMAdditive())
    AdaptiveOpticsSim._measure_branches!(SequentialBranchExecution(), surrogate)
    return surrogate
end

function _run_ao188_post_command_equivalence(::Type{B}, branch_mode::AO188BranchExecutionMode;
    T::Type{<:AbstractFloat}=Float64) where {B<:GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    cmd_src = ao188_3k_surrogate(; params=_ao188_noise_free_params(T), backend=Array, rng=MersenneTwister(1))
    _evaluate_ao188!(cmd_src)
    host_command = Array(cmd_src.command)

    cpu = ao188_3k_surrogate(; params=_ao188_noise_free_params(T), backend=Array, rng=MersenneTwister(2))
    gpu = ao188_3k_surrogate(; params=_ao188_noise_free_params(T; branch_execution=branch_mode), backend=BackendArray, rng=MersenneTwister(2))
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

function _run_lgs_equivalence(::Type{B}, profile::Symbol) where {B<:GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")
    T = Float32

    tel_cpu, src_cpu, wfs_cpu, det_cpu = _build_lgs_case(Array, T, profile)
    tel_gpu, src_gpu, wfs_gpu, det_gpu = _build_lgs_case(BackendArray, T, profile)

    measure!(wfs_cpu, tel_cpu, src_cpu, det_cpu; rng=MersenneTwister(2))
    measure!(wfs_gpu, tel_gpu, src_gpu, det_gpu; rng=MersenneTwister(2))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_gpu.state.slopes))

    println("lgs_sh_equivalence profile=", profile)
    _assert_close("spot_cube", wfs_gpu.state.spot_cube, wfs_cpu.state.spot_cube)
    _assert_close("slopes", wfs_gpu.state.slopes, wfs_cpu.state.slopes)
end

function run_gpu_runtime_equivalence(::Type{B}; branch_mode::AO188BranchExecutionMode=SequentialBranchExecution()) where {B<:GPUBackendTag}
    _run_ao188_equivalence(B, branch_mode)
    _run_lgs_equivalence(B, :none)
    _run_lgs_equivalence(B, :na)
    println("gpu_runtime_equivalence complete")
    return nothing
end

function run_gpu_runtime_equivalence_high_accuracy(::Type{B};
    branch_mode::AO188BranchExecutionMode=SequentialBranchExecution()) where {B<:GPUBackendTag}
    _run_ao188_post_command_equivalence(B, branch_mode; T=Float64)
    println("gpu_runtime_equivalence_high_accuracy complete")
    return nothing
end
