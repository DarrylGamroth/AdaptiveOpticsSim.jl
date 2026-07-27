module AdaptiveOpticsSimCUDAExt

using AdaptiveOpticsSim
using CUDA
using LinearAlgebra

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.CUDABackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.CuArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = :cuda
AdaptiveOpticsSim.gpu_backend_name(::Type{<:CUDA.CuArray}) = :cuda
AdaptiveOpticsSim.array_backend_selector(::Type{<:CUDA.CuArray}) = AdaptiveOpticsSim.CUDABackend()
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.CUDABackendTag}, value, dims::Vararg{Int}) = CUDA.fill(value, dims...)
AdaptiveOpticsSim.compute_device_identifier(array::CUDA.CuArray) =
    CUDA.deviceid(CUDA.device(array))

struct CUDAPreparedDeviceExecutionContext <:
    AdaptiveOpticsSim._AbstractPreparedDeviceExecutionContext
    device::CUDA.CuDevice
    stream::CUDA.CuStream
    compute_device::AdaptiveOpticsSim.AcceleratorComputeDevice
end

function AdaptiveOpticsSim._prepare_device_execution_context(
    storage::CUDA.CuArray,
)
    device = CUDA.device(storage)
    stream = CUDA.device!(device) do
        CUDA.CuStream()
    end
    return CUDAPreparedDeviceExecutionContext(
        device,
        stream,
        AdaptiveOpticsSim.compute_device(storage),
    )
end

@inline AdaptiveOpticsSim._prepared_device_execution_compute_device(
    context::CUDAPreparedDeviceExecutionContext,
) = context.compute_device

function AdaptiveOpticsSim._with_prepared_device_execution_context(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
) where {F}
    return CUDA.device!(context.device) do
        CUDA.stream!(f, context.stream)
    end
end

@inline function AdaptiveOpticsSim._synchronize_prepared_device_execution_context!(
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.synchronize(context.stream)
    return nothing
end

function AdaptiveOpticsSim.solve_lift_fallback!(
    diag::AdaptiveOpticsSim.LiFTDiagnostics{T},
    rhs::CUDA.AnyCuVector{T},
    H::CUDA.AnyCuMatrix{T},
    residual::CUDA.AnyCuVector{T},
    damping::AdaptiveOpticsSim.LiFTDampingMode,
) where {T<:AbstractFloat}
    # CUSOLVER's SVD requires a dense CuMatrix rather than a wrapped view.
    F = svd(CUDA.CuArray(H); full=false)
    λ = AdaptiveOpticsSim.fallback_damping_lambda(damping, T, H)
    work = CUDA.CuArray{T}(undef, length(F.S))
    mul!(work, transpose(F.U), residual)
    @. work = ifelse(iszero(F.S^2 + λ), zero(T), (F.S * work) / (F.S^2 + λ))
    mul!(rhs, adjoint(F.Vt), work)
    diag.regularization = λ
    diag.used_fallback = true
    return rhs
end

end
