module AdaptiveOpticsSimAMDGPUExt

import AdaptiveOpticsSim
import AdaptiveOpticsSim: Backends, Calibration, Tomography,
    WavefrontSensors
using AMDGPU
using AbstractFFTs
using KernelAbstractions
using LinearAlgebra
using Random

#
# AMDGPU backend extension
#
# This extension supplies backend-native dense linear algebra and FFT plumbing
# for the maintained ROCArray execution paths. The main mathematical surfaces
# implemented here are:
#
# - pseudoinverse construction from SVD
# - stable Hermitian right division used by tomography/calibration
# - normal-equation solves and SVD fallback for LiFT
#
# The key rule is that the algorithms match the core implementation, while the
# execution is specialized to rocBLAS / rocSOLVER / rocFFT where that improves
# performance or avoids host fallback.
#
Backends.gpu_backend_loaded(::Type{Backends.AMDGPUBackendTag}) = true
Backends.gpu_backend_array_type(::Type{Backends.AMDGPUBackendTag}) = AMDGPU.ROCArray
Backends.gpu_backend_name(::Type{Backends.AMDGPUBackendTag}) = :amdgpu
Backends.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
Backends.array_backend_selector(::Type{<:AMDGPU.ROCArray}) = Backends.AMDGPUBackend()
Backends.disable_scalar_backend!(::Type{Backends.AMDGPUBackendTag}) = AMDGPU.allowscalar(false)
Backends.backend_rand(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
Backends.backend_randn(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
Backends.backend_zeros(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
Backends.backend_fill(::Type{Backends.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)
Backends.compute_device_identifier(array::AMDGPU.ROCArray) =
    AMDGPU.device_id(AMDGPU.device(array))

function Backends.compute_device_availability(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend,I},
) where {I<:Integer}
    identifier = try
        Int(Backends.compute_device_identifier(device))
    catch
        return Backends.ComputeDeviceUnavailable(:invalid_device_identifier)
    end
    identifier >= 1 || return Backends.ComputeDeviceUnavailable(
        :invalid_device_identifier)
    AMDGPU.functional() || return Backends.ComputeDeviceUnavailable(
        :backend_runtime_unavailable)
    try
        AMDGPU.HIPDevice(identifier)
    catch
        return Backends.ComputeDeviceUnavailable(:device_unavailable)
    end
    return Backends.ComputeDeviceAvailable()
end

Backends.compute_device_availability(
    ::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
) = Backends.ComputeDeviceUnavailable(:invalid_device_identifier)

@noinline function _require_amdgpu_device(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend,I},
) where {I<:Integer}
    availability = Backends.compute_device_availability(device)
    Backends.compute_device_is_available(availability) ||
        Backends._throw_compute_device_error(
            :select,
            Backends.compute_device_unavailable_reason(availability),
            device,
            "AMDGPU cannot address the requested device identifier",
        )
    return AMDGPU.HIPDevice(Int(Backends.compute_device_identifier(device)))
end

@noinline function _require_amdgpu_device(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
)
    Backends._throw_compute_device_error(
        :select,
        :invalid_device_identifier,
        device,
        "AMDGPU device identifiers must be positive integer identifiers",
    )
end

function Backends._with_compute_device(
    f::F,
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
) where {F}
    return AMDGPU.device!(f, _require_amdgpu_device(device))
end

struct AMDGPUPreparedDeviceExecutionContext{
    D<:Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
} <: Backends._AbstractPreparedDeviceExecutionContext
    device::AMDGPU.HIPDevice
    stream::AMDGPU.HIPStream
    compute_device::D
end

function Backends._prepare_device_execution_context(
    storage::AMDGPU.ROCArray,
)
    device = AMDGPU.device(storage)
    stream = AMDGPU.device!(device) do
        AMDGPU.HIPStream()
    end
    return AMDGPUPreparedDeviceExecutionContext(
        device,
        stream,
        Backends.compute_device(storage),
    )
end

function Backends._prepare_device_execution_context(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
)
    runtime_device = _require_amdgpu_device(device)
    stream = AMDGPU.device!(runtime_device) do
        AMDGPU.HIPStream()
    end
    return AMDGPUPreparedDeviceExecutionContext(
        runtime_device,
        stream,
        device,
    )
end

@inline Backends._prepared_device_execution_compute_device(
    context::AMDGPUPreparedDeviceExecutionContext,
) = context.compute_device

function Backends._with_prepared_device_execution_context(
    f::F,
    context::AMDGPUPreparedDeviceExecutionContext,
) where {F}
    return AMDGPU.device!(context.device) do
        # AMDGPU 2.7's scoped `stream!(f, stream)` restores the old stream but
        # does not select the requested stream before invoking `f`. Select and
        # restore explicitly so the prepared stream is authoritative even for
        # nested and multi-device execution contexts.
        old_stream = AMDGPU.stream()
        AMDGPU.stream!(context.stream)
        try
            f()
        finally
            AMDGPU.stream!(old_stream)
        end
    end
end

@inline function Backends._synchronize_prepared_device_execution_context!(
    context::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.synchronize(context.stream)
    return nothing
end

function Backends.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AMDGPU.rocFFT.ROCFFTPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
function Backends.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AbstractFFTs.ScaledPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
Calibration.default_build_backend(::AMDGPU.ROCArray) =
    Calibration.GPUArrayBuildBackend(Backends.AMDGPUBackendTag)
Calibration.prepare_build_matrix(
    ::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AbstractMatrix,
) = Matrix(A)
WavefrontSensors.grouped_accumulation_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:WavefrontSensors.PyramidWFS},
) = WavefrontSensors.GroupedStaged2DStrategy()
WavefrontSensors.grouped_accumulation_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:WavefrontSensors.BiOEdgeWFS},
) = WavefrontSensors.GroupedStaged2DStrategy()
function WavefrontSensors.sh_sensing_execution_strategy(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::WavefrontSensors.ShackHartmannWFS,
)
    return WavefrontSensors.ShackHartmannWFSROCmHostStatsStrategy()
end

AdaptiveOpticsSim.Detectors.detector_execution_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detectors.Detector},
) = AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy()
AdaptiveOpticsSim.Detectors._detector_value_strategy(
    strategy::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
) = strategy
AdaptiveOpticsSim.Detectors.can_apply_device_readout_correction(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AdaptiveOpticsSim.Detectors.FrameReadoutCorrectionModel,
) = false
AdaptiveOpticsSim.Detectors.counting_output_execution_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detectors.AbstractCountingDetector},
    ::Type{<:AMDGPU.ROCArray{T,2}},
) where {T<:Integer} = AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy()
Backends.reduction_execution_strategy(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AMDGPU.ROCArray,
) = Backends.HostMirrorReductionStrategy()
Backends.randn_backend_async!(::Backends.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
Backends._randn_backend!(::Backends.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    out::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    Backends.randn_backend!(rng, out)
    return out
end
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    Backends.randn_backend!(rng, cube)
    return cube
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    img::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.Detectors.detector_host_frame!(det, img)
    Backends._poisson_noise!(Backends.ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.Detectors.detector_host_cube!(det, cube)
    Backends._poisson_noise!(Backends.ScalarCPUStyle(), rng, host)
    copyto!(cube, host)
    return cube
end
function AdaptiveOpticsSim.Atmospheres.randn_phase_noise!(rng::AbstractRNG, out::AMDGPU.ROCArray{T,2}, host::Matrix{T}) where {T<:AbstractFloat}
    if size(host) != size(out)
        host = Matrix{T}(undef, size(out)...)
    end
    randn!(rng, host)
    copyto!(out, host)
    return host
end
function AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    psd::AMDGPU.ROCArray{T,2},
    freqs::AMDGPU.ROCArray{T,1},
    coeff::T,
    inv_L0_sq::T,
    exponent::T,
    inv_fm_sq::T,
    n::Int,
) where {T<:AbstractFloat}
    host_psd = Matrix{T}(undef, size(psd))
    AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(Backends.ScalarCPUStyle(), host_psd,
        Array(freqs), coeff, inv_L0_sq, exponent, inv_fm_sq, n)
    copyto!(psd, host_psd)
    return psd
end

# AMDGPU 2.7/GPUCompiler currently fails IR validation for the variable-trip
# KernelAbstractions slope kernels on gfx1030. Keep these backend-specific
# fallbacks explicit so other accelerator backends retain the device kernels.
function WavefrontSensors._geometric_slopes!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    WavefrontSensors._geometric_slopes!(Backends.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), sub, n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

function WavefrontSensors._edge_geometric_slopes!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    edge_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    WavefrontSensors._edge_geometric_slopes!(Backends.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), Array(edge_mask), sub,
        n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

Backends.backend_matmul(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'N', A, B)
Backends.backend_matmul_transpose_right(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'T', A, B)

function dense_copy_to_roc(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    out = AMDGPU.ROCArray{T}(undef, size(A)...)
    copyto!(out, A)
    return out
end

function dense_host_matrix(A::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat}
    host_parent = Array(parent(A))
    return Matrix(@view host_parent[parentindices(A)...])
end

dense_copy_to_roc(A::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat} =
    AMDGPU.ROCArray(dense_host_matrix(A))

copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::AMDGPU.ROCArray) = copyto!(dest, src)
copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat} =
    copyto!(dest, dense_host_matrix(src))
copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::AbstractMatrix) = copyto!(dest, src)

function roc_svd(A::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    F = copy(A)
    U, S, Vt = AMDGPU.rocSOLVER.gesvd!('S', 'S', F)
    return (; U, S, Vt, s_host=Calibration.singular_values_host(S))
end

function roc_lu_solve!(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    n = size(A, 1)
    ipiv = AMDGPU.ROCArray{Cint}(undef, n)
    AMDGPU.rocSOLVER.getrf!(A, ipiv)
    AMDGPU.rocSOLVER.getrs!('N', A, ipiv, B)
    return B
end

function roc_cholesky_solve!(
    A::AMDGPU.ROCArray{T,2},
    rhs::AMDGPU.ROCArray{T,N};
    check::Bool=false,
) where {T<:AbstractFloat,N}
    rhs_mat = if N == 1
        out = AMDGPU.ROCArray{T}(undef, length(rhs), 1)
        copyto!(vec(out), rhs)
        out
    else
        rhs
    end
    chol = cholesky!(Hermitian(A), check=check)
    if issuccess(chol)
        ldiv!(chol, rhs_mat)
    end
    if N == 1
        copyto!(rhs, vec(rhs_mat))
        return rhs, chol
    end
    return rhs_mat, chol
end

"""
    pseudoinverse_from_roc_svd(backend, U, S, Vt, inv_s_host)

Assemble the pseudoinverse `V * Diagonal(inv_s) * U'` from the compact rocSOLVER
SVD factors.

`rocSOLVER.gesvd!` returns `Vt`, so the final matrix product is expressed with
transpose flags rather than materialized transposes.
"""
function pseudoinverse_from_roc_svd(backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    U::AMDGPU.ROCArray{T,2}, S::AMDGPU.ROCArray{T,1}, Vt::AMDGPU.ROCArray{T,2},
    inv_s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s = Calibration.materialize_build(backend, S, inv_s_host)
    U_scaled = copy(U)
    U_scaled .*= reshape(inv_s, 1, :)
    # Mathematically this is V * U_scaled', but rocSOLVER returns Vt and the
    # arrays are already column-major, so transpose flags are enough here.
    return AMDGPU.rocBLAS.gemm('T', 'T', Vt, U_scaled)
end

function inverse_scaling_and_stats(::Calibration.ExactPseudoInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s_host = similar(s_host)
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = iszero(s_host[i]) ? zero(T) : inv(s_host[i])
    end
    effective_rank = count(!iszero, s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, Calibration.InverseStats(s_host, cond, effective_rank, 0)
end

function inverse_scaling_and_stats(policy::Calibration.TSVDInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.n_trunc >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("TSVD n_trunc must be >= 0"))
    isempty(s_host) && return similar(s_host), Calibration.InverseStats(s_host, T(Inf), 0, 0)
    cutoff = Calibration._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s_host)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s_host = similar(s_host)
    fill!(inv_s_host, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s_host[i] = inv(s_host[i])
    end
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, Calibration.InverseStats(s_host, cond, effective_rank, length(s_host) - effective_rank)
end

function inverse_scaling_and_stats(policy::Calibration.TikhonovInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.lambda >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("Tikhonov lambda must be >= 0"))
    inv_s_host = similar(s_host)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = s_host[i] / (s_host[i]^2 + λ2)
    end
    cutoff = Calibration._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s_host)
    denom = max(isempty(s_host) ? zero(T) : s_host[end], T(policy.lambda))
    cond = (isempty(s_host) || iszero(denom)) ? T(Inf) : s_host[begin] / denom
    return inv_s_host, Calibration.InverseStats(s_host, cond, effective_rank, 0)
end

function roc_inverse_operator(
    backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2},
    policy::Calibration.InversePolicy,
) where {T<:AbstractFloat}
    svd_parts = roc_svd(A)
    inv_s_host, stats = inverse_scaling_and_stats(policy, svd_parts.s_host)
    if isempty(svd_parts.s_host)
        empty_inverse = Calibration.materialize_build(backend, similar(A, T, size(A, 2), size(A, 1)))
        return empty_inverse, stats
    end
    M = pseudoinverse_from_roc_svd(backend, svd_parts.U, svd_parts.S, svd_parts.Vt, inv_s_host)
    return M, stats
end

function Calibration.inverse_operator(backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, ::Calibration.ExactPseudoInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, Calibration.ExactPseudoInverse())
end

function Calibration.inverse_operator(backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::Calibration.TSVDInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, policy)
end

function Calibration.inverse_operator(backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::Calibration.TikhonovInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, policy)
end

"""
    stable_hermitian_right_division(_, rhs, gram)

Solve the right-division `rhs / gram` through a left solve on the transposed
system.

The preferred path is Cholesky on the Hermitian Gram matrix. If that fails, the
implementation falls back to LU so the higher-level algorithm remains robust on
ill-conditioned runtime/calibration cases.
"""
function Tomography.stable_hermitian_right_division(
    _backend::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    rhs::AMDGPU.ROCArray{T,2},
    gram::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    gram_factor = copy(gram)
    rhs_t = permutedims(rhs, (2, 1))
    _, fact = roc_cholesky_solve!(gram_factor, rhs_t; check=false)
    if issuccess(fact)
        return permutedims(rhs_t, (2, 1))
    end
    gram_lu = copy(gram)
    roc_lu_solve!(gram_lu, rhs_t)
    return permutedims(rhs_t, (2, 1))
end

function WavefrontSensors.solve_lift_fallback!(diag::WavefrontSensors.LiFTDiagnostics{T},
    rhs::AMDGPU.ROCArray{T,1}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::WavefrontSensors.LiFTDampingMode) where {T<:AbstractFloat}
    H_mat = dense_copy_to_roc(H)
    svd_parts = roc_svd(H_mat)
    λ = WavefrontSensors.fallback_damping_lambda(damping, T, H_mat)
    work = AMDGPU.ROCArray{T}(undef, length(svd_parts.S))
    mul!(work, transpose(svd_parts.U), residual)
    @. work = ifelse(iszero(svd_parts.S^2 + λ), zero(T), (svd_parts.S * work) / (svd_parts.S^2 + λ))
    mul!(rhs, adjoint(svd_parts.Vt), work)
    diag.regularization = λ
    diag.used_fallback = true
    return rhs
end

roc_factor_matrix(factor::AMDGPU.ROCArray{T,2}) where {T} = factor
roc_factor_matrix(factor::AbstractMatrix{T}) where {T} = dense_copy_to_roc(factor)

"""
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)

Solve the LiFT normal equations on ROCArray inputs.

The main path uses Cholesky on the normal matrix, optionally with diagonal
loading from the damping policy. If repeated factorization attempts fail, the
code falls back to the SVD-based Levenberg-Marquardt solve to preserve the same
robustness guarantees as the CPU implementation.
"""
function WavefrontSensors.solve_normal_system!(diag::WavefrontSensors.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    ::WavefrontSensors.LiFTDampingNone) where {T<:AbstractFloat}
    factor_mat = roc_factor_matrix(factor)
    copy_dense_to_roc!(factor_mat, normal)
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = WavefrontSensors.regularization_load(normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if !issuccess(chol)
            λ *= T(10)
            @views factor_mat[diagind(factor_mat)] .+= λ
            _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
            if !issuccess(chol)
                return WavefrontSensors.solve_lift_fallback!(diag, rhs, H, residual,
                    WavefrontSensors.LiFTLevenbergMarquardt(lambda0=λ))
            end
        end
    end
    diag.regularization = λ
    return rhs
end

function WavefrontSensors.solve_normal_system!(diag::WavefrontSensors.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::WavefrontSensors.LiFTLevenbergMarquardt) where {T<:AbstractFloat}
    factor_mat = roc_factor_matrix(factor)
    copy_dense_to_roc!(factor_mat, normal)
    λ = WavefrontSensors.damping_lambda(damping, normal)
    if λ > zero(T)
        @views factor_mat[diagind(factor_mat)] .+= λ
    end
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), WavefrontSensors.regularization_load(normal))
        copy_dense_to_roc!(factor_mat, normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if λ > T(1e12)
            return WavefrontSensors.solve_lift_fallback!(diag, rhs, H, residual, damping)
        end
    end
    diag.regularization = λ
    return rhs
end

end
