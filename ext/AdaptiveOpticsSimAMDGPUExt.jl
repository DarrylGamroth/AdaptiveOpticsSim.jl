module AdaptiveOpticsSimAMDGPUExt

using AdaptiveOpticsSim
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
AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.ROCArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = :amdgpu
AdaptiveOpticsSim.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
AdaptiveOpticsSim.array_backend_selector(::Type{<:AMDGPU.ROCArray}) = AdaptiveOpticsSim.AMDGPUBackend()
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)
AdaptiveOpticsSim.physical_device_identifier(array::AMDGPU.ROCArray) =
    AMDGPU.device_id(AMDGPU.device(array))
function AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AMDGPU.rocFFT.ROCFFTPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
function AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AbstractFFTs.ScaledPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
AdaptiveOpticsSim.default_build_backend(::AMDGPU.ROCArray) = AdaptiveOpticsSim.GPUArrayBuildBackend(AdaptiveOpticsSim.AMDGPUBackendTag)
AdaptiveOpticsSim.prepare_build_matrix(::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag}, A::AbstractMatrix) = Matrix(A)
AdaptiveOpticsSim.grouped_accumulation_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.PyramidWFS},
) = AdaptiveOpticsSim.GroupedStaged2DPlan()
AdaptiveOpticsSim.grouped_accumulation_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.BioEdgeWFS},
) = AdaptiveOpticsSim.GroupedStaged2DPlan()
function AdaptiveOpticsSim.sh_sensing_execution_plan(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AdaptiveOpticsSim.ShackHartmannWFS,
)
    return AdaptiveOpticsSim.ShackHartmannWFSRocmHostStatsPlan()
end

function AdaptiveOpticsSim.compute_intensity_safe!(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    wfs::AdaptiveOpticsSim.ShackHartmannWFS,
    tel::AdaptiveOpticsSim.Telescope,
    src::AdaptiveOpticsSim.AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int,
)
    propagation = wfs.front_end.propagation
    T = eltype(propagation.intensity)
    field_host = Matrix{Complex{T}}(undef, size(propagation.field)...)
    fill!(field_host, zero(eltype(field_host)))
    # AMDGPU 2.7 routes host conversion of a ROCArray view through a generic
    # gather kernel that is not usable on this path. Copy the dense parents
    # through the direct transfer path, then slice on the host for this
    # deliberately conservative ROCm fallback.
    reflectivity_host = @view Array(AdaptiveOpticsSim.pupil_reflectivity(tel))[xs:xe, ys:ye]
    opd_host = @view Array(tel.state.opd)[xs:xe, ys:ye]
    phasor_host = Array(propagation.phasor)
    opd_to_cycles = T(2) / AdaptiveOpticsSim.wavelength(src)
    amp_scale = sqrt(T(AdaptiveOpticsSim.photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    @views @. field_host[ox+1:ox+sub, oy+1:oy+sub] =
        amp_scale * sqrt(reflectivity_host) * cispi(opd_to_cycles * opd_host)
    @. field_host *= phasor_host
    fft_host = AbstractFFTs.fft(field_host)
    intensity_scale = AdaptiveOpticsSim.sh_fft_intensity_scale(T, size(field_host, 1))
    intensity_host = @. abs2(fft_host) * intensity_scale
    copyto!(propagation.intensity, intensity_host)
    AMDGPU.synchronize()
    return propagation.intensity
end

AdaptiveOpticsSim.detector_execution_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detector},
) = AdaptiveOpticsSim.DetectorHostMirrorPlan()
AdaptiveOpticsSim._detector_value_plan(
    plan::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
) = plan
AdaptiveOpticsSim.can_apply_device_readout_correction(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AdaptiveOpticsSim.FrameReadoutCorrectionModel,
) = false
AdaptiveOpticsSim.counting_output_execution_plan(
    ::Type{<:AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.AbstractCountingDetector},
    ::Type{<:AMDGPU.ROCArray{T,2}},
) where {T<:Integer} = AdaptiveOpticsSim.DetectorHostMirrorPlan()
AdaptiveOpticsSim.reduction_execution_plan(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AMDGPU.ROCArray,
) = AdaptiveOpticsSim.HostMirrorReductionPlan()
AdaptiveOpticsSim.randn_backend_async!(::AdaptiveOpticsSim.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
AdaptiveOpticsSim._randn_backend!(::AdaptiveOpticsSim.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
function AdaptiveOpticsSim._randn_frame_noise!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    out::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    AdaptiveOpticsSim.randn_backend!(rng, out)
    return out
end
function AdaptiveOpticsSim._randn_frame_noise!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    AdaptiveOpticsSim.randn_backend!(rng, cube)
    return cube
end
function AdaptiveOpticsSim._poisson_noise_frame!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    img::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.detector_host_frame!(det, img)
    AdaptiveOpticsSim._poisson_noise!(AdaptiveOpticsSim.ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end
function AdaptiveOpticsSim._poisson_noise_frame!(
    ::AdaptiveOpticsSim.DetectorHostMirrorPlan,
    det::AdaptiveOpticsSim.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.detector_host_cube!(det, cube)
    AdaptiveOpticsSim._poisson_noise!(AdaptiveOpticsSim.ScalarCPUStyle(), rng, host)
    copyto!(cube, host)
    return cube
end
function AdaptiveOpticsSim.randn_phase_noise!(rng::AbstractRNG, out::AMDGPU.ROCArray{T,2}, host::Matrix{T}) where {T<:AbstractFloat}
    if size(host) != size(out)
        host = Matrix{T}(undef, size(out)...)
    end
    randn!(rng, host)
    copyto!(out, host)
    return host
end
function AdaptiveOpticsSim._fill_phase_psd!(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    psd::AMDGPU.ROCArray{T,2},
    freqs::AMDGPU.ROCArray{T,1},
    coeff::T,
    two_pi_sq::T,
    inv_L0_sq::T,
    exponent::T,
    inv_fm_sq::T,
    n::Int,
) where {T<:AbstractFloat}
    host_psd = Matrix{T}(undef, size(psd))
    AdaptiveOpticsSim._fill_phase_psd!(AdaptiveOpticsSim.ScalarCPUStyle(), host_psd,
        Array(freqs), coeff, two_pi_sq, inv_L0_sq, exponent, inv_fm_sq, n)
    copyto!(psd, host_psd)
    return psd
end

# AMDGPU 2.7/GPUCompiler currently fails IR validation for the variable-trip
# KernelAbstractions slope kernels on gfx1030. Keep these backend-specific
# fallbacks explicit so other accelerator backends retain the device kernels.
function AdaptiveOpticsSim._geometric_slopes!(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    AdaptiveOpticsSim._geometric_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), sub, n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

function AdaptiveOpticsSim._edge_geometric_slopes!(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    edge_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    AdaptiveOpticsSim._edge_geometric_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), Array(edge_mask), sub,
        n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

AdaptiveOpticsSim.backend_matmul(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'N', A, B)
AdaptiveOpticsSim.backend_matmul_transpose_right(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
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
    return (; U, S, Vt, s_host=AdaptiveOpticsSim.singular_values_host(S))
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
function pseudoinverse_from_roc_svd(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    U::AMDGPU.ROCArray{T,2}, S::AMDGPU.ROCArray{T,1}, Vt::AMDGPU.ROCArray{T,2},
    inv_s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s = AdaptiveOpticsSim.materialize_build(backend, S, inv_s_host)
    U_scaled = copy(U)
    U_scaled .*= reshape(inv_s, 1, :)
    # Mathematically this is V * U_scaled', but rocSOLVER returns Vt and the
    # arrays are already column-major, so transpose flags are enough here.
    return AMDGPU.rocBLAS.gemm('T', 'T', Vt, U_scaled)
end

function inverse_scaling_and_stats(::AdaptiveOpticsSim.ExactPseudoInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    inv_s_host = similar(s_host)
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = iszero(s_host[i]) ? zero(T) : inv(s_host[i])
    end
    effective_rank = count(!iszero, s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function inverse_scaling_and_stats(policy::AdaptiveOpticsSim.TSVDInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.n_trunc >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("TSVD n_trunc must be >= 0"))
    isempty(s_host) && return similar(s_host), AdaptiveOpticsSim.InverseStats(s_host, T(Inf), 0, 0)
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s_host)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s_host = similar(s_host)
    fill!(inv_s_host, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s_host[i] = inv(s_host[i])
    end
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, length(s_host) - effective_rank)
end

function inverse_scaling_and_stats(policy::AdaptiveOpticsSim.TikhonovInverse, s_host::AbstractVector{T}) where {T<:AbstractFloat}
    policy.lambda >= 0 || throw(AdaptiveOpticsSim.InvalidConfiguration("Tikhonov lambda must be >= 0"))
    inv_s_host = similar(s_host)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = s_host[i] / (s_host[i]^2 + λ2)
    end
    cutoff = AdaptiveOpticsSim._inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s_host)
    denom = max(isempty(s_host) ? zero(T) : s_host[end], T(policy.lambda))
    cond = (isempty(s_host) || iszero(denom)) ? T(Inf) : s_host[begin] / denom
    return inv_s_host, AdaptiveOpticsSim.InverseStats(s_host, cond, effective_rank, 0)
end

function roc_inverse_operator(
    backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2},
    policy::AdaptiveOpticsSim.InversePolicy,
) where {T<:AbstractFloat}
    svd_parts = roc_svd(A)
    inv_s_host, stats = inverse_scaling_and_stats(policy, svd_parts.s_host)
    if isempty(svd_parts.s_host)
        empty_inverse = AdaptiveOpticsSim.materialize_build(backend, similar(A, T, size(A, 2), size(A, 1)))
        return empty_inverse, stats
    end
    M = pseudoinverse_from_roc_svd(backend, svd_parts.U, svd_parts.S, svd_parts.Vt, inv_s_host)
    return M, stats
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, ::AdaptiveOpticsSim.ExactPseudoInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, AdaptiveOpticsSim.ExactPseudoInverse())
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TSVDInverse) where {T<:AbstractFloat}
    return roc_inverse_operator(backend, A, policy)
end

function AdaptiveOpticsSim.inverse_operator(backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
    A::AMDGPU.ROCArray{T,2}, policy::AdaptiveOpticsSim.TikhonovInverse) where {T<:AbstractFloat}
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
function AdaptiveOpticsSim.stable_hermitian_right_division(
    _backend::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag},
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

function AdaptiveOpticsSim.solve_lift_fallback!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T},
    rhs::AMDGPU.ROCArray{T,1}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::AdaptiveOpticsSim.LiFTDampingMode) where {T<:AbstractFloat}
    H_mat = dense_copy_to_roc(H)
    svd_parts = roc_svd(H_mat)
    λ = AdaptiveOpticsSim.fallback_damping_lambda(damping, T, H_mat)
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
function AdaptiveOpticsSim.solve_normal_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    ::AdaptiveOpticsSim.LiFTDampingNone) where {T<:AbstractFloat}
    factor_mat = roc_factor_matrix(factor)
    copy_dense_to_roc!(factor_mat, normal)
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = AdaptiveOpticsSim.regularization_load(normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if !issuccess(chol)
            λ *= T(10)
            @views factor_mat[diagind(factor_mat)] .+= λ
            _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
            if !issuccess(chol)
                return AdaptiveOpticsSim.solve_lift_fallback!(diag, rhs, H, residual,
                    AdaptiveOpticsSim.LiFTLevenbergMarquardt(lambda0=λ))
            end
        end
    end
    diag.regularization = λ
    return rhs
end

function AdaptiveOpticsSim.solve_normal_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T}, rhs::AMDGPU.ROCArray{T,1},
    factor::AbstractMatrix{T}, normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::AdaptiveOpticsSim.LiFTLevenbergMarquardt) where {T<:AbstractFloat}
    factor_mat = roc_factor_matrix(factor)
    copy_dense_to_roc!(factor_mat, normal)
    λ = AdaptiveOpticsSim.damping_lambda(damping, normal)
    if λ > zero(T)
        @views factor_mat[diagind(factor_mat)] .+= λ
    end
    _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), AdaptiveOpticsSim.regularization_load(normal))
        copy_dense_to_roc!(factor_mat, normal)
        @views factor_mat[diagind(factor_mat)] .+= λ
        _, chol = roc_cholesky_solve!(factor_mat, rhs; check=false)
        if λ > T(1e12)
            return AdaptiveOpticsSim.solve_lift_fallback!(diag, rhs, H, residual, damping)
        end
    end
    diag.regularization = λ
    return rhs
end

end
