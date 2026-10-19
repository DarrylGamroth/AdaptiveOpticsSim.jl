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

function Backends._prepare_mul_coefficients(
    output::CUDA.StridedCuVecOrMat{T},
) where {T<:Union{Float32,Float64,ComplexF32,ComplexF64}}
    return Backends._with_compute_device(Backends.compute_device(output)) do
        values = CUDA.CuArray(T[one(T), zero(T)])
        return (alpha=CUDA.CUDACore.CuRefArray(values, 1),
            beta=CUDA.CUDACore.CuRefArray(values, 2))
    end
end

# Passing retained CuRefArray values directly to the typed BLAS wrappers avoids
# constructing transient device scalar allocations during graph capture.
@inline function Backends._enqueue_mul!(
    C::CUDA.StridedCuMatrix{T},
    A::CUDA.StridedCuMatrix{T},
    B::CUDA.StridedCuMatrix{T},
    coefficients::NamedTuple{(:alpha,:beta),Tuple{R,R}},
) where {T<:Union{Float32,Float64,ComplexF32,ComplexF64},
        R<:CUDA.CUDACore.CuRefArray{T}}
    return CUDA.CUBLAS.gemm!('N', 'N', coefficients.alpha, A, B,
        coefficients.beta, C)
end

@inline function Backends._enqueue_mul!(
    C::CUDA.StridedCuMatrix{T},
    A::CUDA.StridedCuMatrix{T},
    B::Transpose{T,<:CUDA.StridedCuMatrix{T}},
    coefficients::NamedTuple{(:alpha,:beta),Tuple{R,R}},
) where {T<:Union{Float32,Float64,ComplexF32,ComplexF64},
        R<:CUDA.CUDACore.CuRefArray{T}}
    return CUDA.CUBLAS.gemm!('N', 'T', coefficients.alpha, A, parent(B),
        coefficients.beta, C)
end

@inline function Backends._enqueue_mul!(
    C::CUDA.StridedCuVector{T},
    A::CUDA.StridedCuMatrix{T},
    B::CUDA.StridedCuVector{T},
    coefficients::NamedTuple{(:alpha,:beta),Tuple{R,R}},
) where {T<:Union{Float32,Float64,ComplexF32,ComplexF64},
        R<:CUDA.CUDACore.CuRefArray{T}}
    return CUDA.CUBLAS.gemv!('N', coefficients.alpha, A, B,
        coefficients.beta, C)
end

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

# The prepared captured stream has one submission owner and no host callbacks.
# Observe completion cooperatively before the public synchronization call.
# CUDA 6.4.2's cooperative wait allocates even for an already-complete stream;
# its public blocking path skips that setup and retains driver and Julia kernel
# error checks. Earlier versions instead need the default completed-stream fast
# path to avoid their unconditional blocking-wait memory-pressure heuristic.
# Other CUDA policies retain semantics but need allocation qualification.
# Keep safepoints while waiting; ordinary Julia GC remains enabled.
@inline function Backends._synchronize_prepared_device_execution_context_blocking!(
    context::CUDAPreparedDeviceExecutionContext,
)
    while !CUDA.isdone(context.stream)
        GC.safepoint()
        yield()
    end
    @static if pkgversion(CUDA) >= v"6.4.2"
        CUDA.synchronize(context.stream; blocking=true)
    else
        CUDA.synchronize(context.stream)
    end
    return nothing
end

@inline function Backends._prepare_device_execution_event(
    ::CUDAPreparedDeviceExecutionContext,
)
    return CUDA.CuEvent(CUDA.EVENT_DISABLE_TIMING)
end

@inline function Backends._record_prepared_device_execution_event!(
    event::CUDA.CuEvent,
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.record(event, context.stream)
    return nothing
end

@inline function Backends._wait_prepared_device_execution_event!(
    event::CUDA.CuEvent,
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.wait(event, context.stream)
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

struct CUDAOutputSnapshot{D,H,V,C,S,T}
    device_values::D
    host_owner::H
    values::V
    source_context::C
    transfer_context::C
    snapshot_graph::S
    transfer_graph::T
    ready::CUDA.CuEvent
    complete::CUDA.CuEvent
end

function _capture_output_copy(destination, source, context)
    return Backends._with_prepared_device_execution_context(context) do
        # Array APIs establish and compile their exact storage bindings cold.
        # Replay only launches retained commands, so CuArray cross-stream
        # bookkeeping cannot introduce a source-side wait on the consumer.
        copyto!(destination, source)
        Backends._synchronize_prepared_device_execution_context!(context)
        Backends._capture_prepared_device_graph(context) do
            copyto!(destination, source)
            nothing
        end
    end
end

function Backends._prepare_output_snapshot(source::CUDA.CuArray{T,N},
    source_context::CUDAPreparedDeviceExecutionContext) where {T,N}
    return Backends._with_prepared_device_execution_context(source_context) do
        device_values = similar(source)
        host_owner = CUDA.CuArray{T,N,CUDA.HostMemory}(undef, size(source))
        # This public array-to-array wrap creates a host view without exposing
        # addresses. The retained host_owner must outlive every view reader.
        values = Base.unsafe_wrap(Array, host_owner)
        fill!(values, zero(T))
        transfer_context = Backends._prepare_device_execution_context(
            source_context.compute_device)
        ready = CUDA.CuEvent(CUDA.EVENT_DISABLE_TIMING)
        complete = CUDA.CuEvent(CUDA.EVENT_DISABLE_TIMING)
        snapshot_graph = _capture_output_copy(device_values, source, source_context)
        transfer_graph = _capture_output_copy(host_owner, device_values, transfer_context)
        CUDAOutputSnapshot(device_values, host_owner, values, source_context,
            transfer_context, snapshot_graph, transfer_graph, ready, complete)
    end
end

@inline function Backends._submit_output_snapshot!(transfer::CUDAOutputSnapshot, source)
    Backends._with_prepared_device_execution_context(transfer.source_context) do
        Backends._launch_prepared_device_graph!(transfer.snapshot_graph, transfer.source_context)
        CUDA.record(transfer.ready, transfer.source_context.stream)
    end
    return nothing
end

@inline function Backends._start_output_snapshot_transfer!(transfer::CUDAOutputSnapshot)
    Backends._with_prepared_device_execution_context(transfer.transfer_context) do
        CUDA.wait(transfer.ready, transfer.transfer_context.stream)
        Backends._launch_prepared_device_graph!(transfer.transfer_graph, transfer.transfer_context)
        CUDA.record(transfer.complete, transfer.transfer_context.stream)
    end
    return nothing
end

@inline function Backends._poll_output_snapshot_transfer!(transfer::CUDAOutputSnapshot)
    return Backends._with_prepared_device_execution_context(transfer.transfer_context) do
        CUDA.isdone(transfer.complete)
    end
end

function Backends._drain_output_snapshot!(transfer::CUDAOutputSnapshot)
    # Cold shutdown/reset only, after joining the consumer. Drain both streams
    # even when submission failed before recording either completion event.
    try
        Backends._with_prepared_device_execution_context(transfer.source_context) do
            Backends._synchronize_prepared_device_execution_context_blocking!(
                transfer.source_context)
        end
    finally
        Backends._with_prepared_device_execution_context(transfer.transfer_context) do
            Backends._synchronize_prepared_device_execution_context_blocking!(
                transfer.transfer_context)
        end
    end
    return nothing
end

end
