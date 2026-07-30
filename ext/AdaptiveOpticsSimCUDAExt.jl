module AdaptiveOpticsSimCUDAExt

import AdaptiveOpticsSim
import AdaptiveOpticsSim: Backends, WavefrontSensors
using CUDA
using LinearAlgebra

Backends.gpu_backend_loaded(::Type{Backends.CUDABackendTag}) = true
Backends.gpu_backend_array_type(::Type{Backends.CUDABackendTag}) = CUDA.CuArray
Backends.gpu_backend_name(::Type{Backends.CUDABackendTag}) = :cuda
Backends.gpu_backend_name(::Type{<:CUDA.CuArray}) = :cuda
Backends.array_backend_selector(::Type{<:CUDA.CuArray}) = Backends.CUDABackend()
Backends.disable_scalar_backend!(::Type{Backends.CUDABackendTag}) = CUDA.allowscalar(false)
Backends.backend_rand(::Type{Backends.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.rand(T, dims...)
Backends.backend_randn(::Type{Backends.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.randn(T, dims...)
Backends.backend_zeros(::Type{Backends.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.zeros(T, dims...)
Backends.backend_fill(::Type{Backends.CUDABackendTag}, value, dims::Vararg{Int}) = CUDA.fill(value, dims...)
Backends.compute_device_identifier(array::CUDA.CuArray) =
    CUDA.deviceid(CUDA.device(array))

struct CUDAPreparedDeviceExecutionContext <:
    Backends._AbstractPreparedDeviceExecutionContext
    device::CUDA.CuDevice
    stream::CUDA.CuStream
    compute_device::Backends.AcceleratorComputeDevice
end

function Backends._prepare_device_execution_context(
    storage::CUDA.CuArray,
)
    device = CUDA.device(storage)
    stream = CUDA.device!(device) do
        CUDA.CuStream()
    end
    return CUDAPreparedDeviceExecutionContext(
        device,
        stream,
        Backends.compute_device(storage),
    )
end

@inline Backends._prepared_device_execution_compute_device(
    context::CUDAPreparedDeviceExecutionContext,
) = context.compute_device

function Backends._with_prepared_device_execution_context(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
) where {F}
    return CUDA.device!(context.device) do
        CUDA.stream!(f, context.stream)
    end
end

@inline function Backends._synchronize_prepared_device_execution_context!(
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.synchronize(context.stream)
    return nothing
end

function WavefrontSensors.solve_lift_fallback!(
    diag::WavefrontSensors.LiFTDiagnostics{T},
    rhs::CUDA.AnyCuVector{T},
    H::CUDA.AnyCuMatrix{T},
    residual::CUDA.AnyCuVector{T},
    damping::WavefrontSensors.LiFTDampingMode,
) where {T<:AbstractFloat}
    # CUSOLVER's SVD requires a dense CuMatrix rather than a wrapped view.
    F = svd(CUDA.CuArray(H); full=false)
    λ = WavefrontSensors.fallback_damping_lambda(damping, T, H)
    work = CUDA.CuArray{T}(undef, length(F.S))
    mul!(work, transpose(F.U), residual)
    @. work = ifelse(iszero(F.S^2 + λ), zero(T), (F.S * work) / (F.S^2 + λ))
    mul!(rhs, adjoint(F.Vt), work)
    diag.regularization = λ
    diag.used_fallback = true
    return rhs
end

end
