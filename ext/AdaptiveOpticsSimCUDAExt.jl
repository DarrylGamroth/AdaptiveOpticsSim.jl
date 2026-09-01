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

function Backends._prepare_graph_rng(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
    seed::UInt64,
)
    return Backends._prepare_counter_rng(device, seed)
end

function Backends.compute_device_availability(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend,I},
) where {I<:Integer}
    identifier = try
        Int(Backends.compute_device_identifier(device))
    catch
        return Backends.ComputeDeviceUnavailable(:invalid_device_identifier)
    end
    identifier >= 0 || return Backends.ComputeDeviceUnavailable(
        :invalid_device_identifier)
    CUDA.functional() || return Backends.ComputeDeviceUnavailable(
        :backend_runtime_unavailable)
    try
        CUDA.CuDevice(identifier)
    catch
        return Backends.ComputeDeviceUnavailable(:device_unavailable)
    end
    return Backends.ComputeDeviceAvailable()
end

Backends.compute_device_availability(
    ::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
) = Backends.ComputeDeviceUnavailable(:invalid_device_identifier)

@noinline function _require_cuda_device(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend,I},
) where {I<:Integer}
    availability = Backends.compute_device_availability(device)
    Backends.compute_device_is_available(availability) ||
        Backends._throw_compute_device_error(
            :select,
            Backends.compute_device_unavailable_reason(availability),
            device,
            "CUDA cannot address the requested device ordinal",
        )
    return CUDA.CuDevice(Int(Backends.compute_device_identifier(device)))
end

@noinline function _require_cuda_device(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
)
    Backends._throw_compute_device_error(
        :select,
        :invalid_device_identifier,
        device,
        "CUDA device identifiers must be integer ordinals",
    )
end

function Backends._with_compute_device(
    f::F,
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
) where {F}
    return CUDA.device!(f, _require_cuda_device(device))
end

struct CUDAPreparedDeviceExecutionContext{
    D<:Backends.AcceleratorComputeDevice{Backends.CUDABackend},
} <: Backends._AbstractPreparedDeviceExecutionContext
    device::CUDA.CuDevice
    stream::CUDA.CuStream
    compute_device::D
end

struct CUDAPreparedDeviceGraph{Graph,Executable}
    graph::Graph
    executable::Executable
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

function Backends._prepare_device_execution_context(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
)
    runtime_device = _require_cuda_device(device)
    stream = CUDA.device!(runtime_device) do
        CUDA.CuStream()
    end
    return CUDAPreparedDeviceExecutionContext(
        runtime_device,
        stream,
        device,
    )
end

@inline Backends._prepared_device_execution_compute_device(
    context::CUDAPreparedDeviceExecutionContext,
) = context.compute_device

@inline function Backends._with_prepared_device_execution_context(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
) where {F}
    target_context = CUDA.context(context.device)
    CUDA.CUDACore.isvalid(target_context) || throw(
        AdaptiveOpticsSim.InvalidConfiguration(
            "the prepared CUDA context is no longer valid",
        ),
    )
    old_context = CUDA.context!(target_context)
    old_stream = CUDA.stream()
    CUDA.stream!(context.stream)
    try
        return f()
    finally
        CUDA.stream!(old_stream)
        if old_context !== nothing && old_context != target_context &&
                CUDA.CUDACore.isvalid(old_context)
            CUDA.context!(old_context)
        end
    end
end

@inline function Backends._synchronize_prepared_device_execution_context!(
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.synchronize(context.stream)
    return nothing
end

function Backends._capture_prepared_device_graph(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
) where {F}
    graph = CUDA.capture(f)
    executable = CUDA.instantiate(graph)
    return CUDAPreparedDeviceGraph(graph, executable)
end

@inline function Backends._launch_prepared_device_graph!(
    captured::CUDAPreparedDeviceGraph,
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.launch(captured.executable, context.stream)
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
