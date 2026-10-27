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

# Fixed-sampling sweeps qualify these FFT sizes on the sm86 / 5 MiB device
# class. A broader cache-size rule regressed at side 560: do not interpolate
# this table to other FFT sizes or point counts. Other devices and precisions
# keep the existing limit. See the modulation-batches benchmark results.
function _pyramid_modulation_batch_device_qualified(field)
    device = CUDA.device(field)
    CUDA.capability(device) == v"8.6" || return false
    CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_L2_CACHE_SIZE) == 5 * 1024 * 1024 ||
        return false
    return size(field, 1) == size(field, 2)
end

function WavefrontSensors._pyramid_modulation_batch_limit(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend},
    field::AbstractMatrix{ComplexF32}, point_count::Int, plane_bytes::Int,
)
    default_limit = WavefrontSensors._PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT
    _pyramid_modulation_batch_device_qualified(field) || return default_limit
    side = size(field, 1)
    if (side == 2880 && point_count == 192) ||
            (side == 3072 && point_count == 256)
        return 16
    end
    limit = if (side == 256 && point_count == 32) ||
            (side == 272 && point_count in (24, 30, 32))
        8
    elseif (side == 408 && point_count == 32) ||
            (side == 512 && point_count in (12, 32, 96)) ||
            (side == 528 && point_count == 64)
        2
    else
        return default_limit
    end
    candidate = WavefrontSensors._pyramid_modulation_batch_size(
        point_count, plane_bytes, limit)
    # A prime point count must not collapse to one FFT per launch merely to
    # satisfy the smaller limit.
    return candidate >= 2 ? candidate : default_limit
end

# The finer physical propagation grids benefit from sixteen transforms per
# batch. This raises the main field-stack bound only for measured combinations;
# neighboring sizes, precisions, point counts and device classes retain 256 MiB.
function WavefrontSensors._pyramid_modulation_batch_budget_bytes(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend},
    field::AbstractMatrix{ComplexF32}, point_count::Int, plane_bytes::Int,
)
    default_budget = WavefrontSensors._PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES
    _pyramid_modulation_batch_device_qualified(field) || return default_budget
    side = size(field, 1)
    if (side == 2880 && point_count == 192) ||
            (side == 3072 && point_count == 256)
        return 1536 * 1024 * 1024
    end
    return default_budget
end

# Keep the frequency planes transposed between the forward and inverse axis
# transforms only where the complete captured physical path was measured.
function WavefrontSensors._prepare_pyramid_modulation_fft(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend},
    field_stack::CUDA.CuArray{ComplexF32,3}, point_count::Int,
)
    if point_count == 192 && size(field_stack, 1) == 2880 &&
            _pyramid_modulation_batch_device_qualified(field_stack)
        return WavefrontSensors._prepare_pyramid_transpose_fft(field_stack)
    end
    return WavefrontSensors._prepare_pyramid_native_fft(field_stack)
end

# Synchronous kernel phases retain their current-stream completion boundary.
# CUDA's cooperative wait can enter worker/task bookkeeping for longer kernels;
# the blocking wait avoids that latency and still checks device exceptions.
# Low-level callers that need task cooperation can select blocking=false.
@inline function Backends.synchronize_backend!(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend};
    blocking::Bool=true,
)
    CUDA.synchronize(; blocking)
    return nothing
end

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
    context::CUDA.CuContext
    stream::CUDA.CuStream
    compute_device::D
end

struct CUDAPreparedDeviceGraph{Graph,Executable}
    graph::Graph
    executable::Executable
    # cuBLAS kernels recorded into the graph can retain this scratch address.
    workspace::CUDA.CuArray{UInt8,1,CUDA.DeviceMemory}
end

function Backends._prepare_device_execution_context(
    storage::CUDA.CuArray,
)
    device = CUDA.device(storage)
    stream = CUDA.device!(device) do
        CUDA.CuStream(; flags=CUDA.STREAM_NON_BLOCKING)
    end
    return CUDAPreparedDeviceExecutionContext(
        device,
        CUDA.context(device),
        stream,
        Backends.compute_device(storage),
    )
end

function Backends._prepare_device_execution_context(
    device::Backends.AcceleratorComputeDevice{Backends.CUDABackend},
)
    runtime_device = _require_cuda_device(device)
    stream = CUDA.device!(runtime_device) do
        CUDA.CuStream(; flags=CUDA.STREAM_NON_BLOCKING)
    end
    return CUDAPreparedDeviceExecutionContext(
        runtime_device,
        CUDA.context(runtime_device),
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
    # Retain the exact context that owns the stream. Resolving it from the
    # device again also enters CUDA's context cache on every execution.
    target_context = context.context
    CUDA.CUDACore.isvalid(target_context) || throw(
        AdaptiveOpticsSim.InvalidConfiguration(
            "the prepared CUDA context is no longer valid",
        ),
    )
    return CUDA.context!(target_context) do
        CUDA.stream!(context.stream) do
            f()
        end
    end
end

@inline function Backends._synchronize_prepared_device_execution_context!(
    context::CUDAPreparedDeviceExecutionContext,
)
    CUDA.synchronize(context.stream)
    return nothing
end

function Backends._copy_prepared_device_buffer!(
    context::CUDAPreparedDeviceExecutionContext,
    destination::Array{T,N}, source::CUDA.CuArray{T,N,CUDA.DeviceMemory},
) where {T,N}
    isbitstype(T) || return copyto!(destination, source)
    size(destination) == size(source) || throw(AdaptiveOpticsSim.DimensionMismatchError(
        "prepared CUDA host copies require equal array shapes",
    ))
    # The HIL owner has completed graph work and retains both exact buffers
    # through its final prepared-stream wait. Avoid copyto!'s two extra
    # cooperative waits; pointer conversion still tracks CUDA array use.
    GC.@preserve destination source begin
        unsafe_copyto!(pointer(destination), pointer(source), length(source);
            stream=context.stream, async=true)
    end
    return destination
end

# The prepared captured stream has one submission owner and no host callbacks.
# Observe completion cooperatively before the public synchronization call.
# CUDA 6.4.2's cooperative wait allocates even for an already-complete stream.
# CUDA 6.5 has a completed-stream fast path that avoids pre-wait GC/statistics
# work. Both public paths retain resource bookkeeping and kernel error checks.
# Other CUDA policies retain semantics but need allocation qualification.
# Keep safepoints while waiting; ordinary Julia GC remains enabled.
@inline function Backends._synchronize_prepared_device_execution_context_blocking!(
    context::CUDAPreparedDeviceExecutionContext,
)
    while !CUDA.isdone(context.stream)
        GC.safepoint()
        yield()
    end
    @static if Base.pkgversion(CUDA) >= v"6.5"
        CUDA.synchronize(context.stream; blocking=false)
    else
        CUDA.synchronize(context.stream; blocking=true)
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

function _with_prepared_cublas_workspace(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
    workspace::CUDA.CuArray{UInt8,1,CUDA.DeviceMemory},
) where {F}
    # The caller has activated the prepared stream. Obtain its task-local
    # handle first: CUDA resets the workspace when updating a handle's stream.
    handle = CUDA.CUBLAS.handle()
    UInt(pointer(workspace)) % 256 == 0 || throw(
        AdaptiveOpticsSim.InvalidConfiguration(
            "captured cuBLAS workspace must be 256-byte aligned",
        ),
    )
    GC.@preserve workspace begin
        CUDA.CUBLAS.cublasSetWorkspace_v2(handle, pointer(workspace), sizeof(workspace))
        try
            return f()
        finally
            # Setting the stream restores cuBLAS's default workspace, even for
            # the same stream, without disturbing CUDA's task-local handle cache.
            CUDA.CUBLAS.cublasSetStream_v2(handle, context.stream)
        end
    end
end

function _require_cuda_graph_fixed_storage(graph)
    driver = CUDA.CUDACore
    count = Ref{Csize_t}(0)
    driver.cuGraphGetNodes(graph, C_NULL, count)
    iszero(count[]) && return graph
    nodes = Vector{driver.CUgraphNode}(undef, count[])
    driver.cuGraphGetNodes(graph, nodes, count)
    for node in nodes
        kind = Ref{driver.CUgraphNodeType}()
        driver.cuGraphNodeGetType(node, kind)
        if kind[] == driver.CU_GRAPH_NODE_TYPE_MEM_ALLOC ||
                kind[] == driver.CU_GRAPH_NODE_TYPE_MEM_FREE
            throw(AdaptiveOpticsSim.InvalidConfiguration(
                "captured CUDA graph contains a device allocation or free node",
            ))
        elseif kind[] == driver.CU_GRAPH_NODE_TYPE_GRAPH
            child = Ref{driver.CUgraph}(C_NULL)
            driver.cuGraphChildGraphNodeGetGraph(node, child)
            _require_cuda_graph_fixed_storage(child[])
        end
    end
    return graph
end

function Backends._capture_prepared_device_graph(
    f::F,
    context::CUDAPreparedDeviceExecutionContext,
) where {F}
    # Keep capture-only scratch with the executable, so ordinary stream owners
    # do not pay this cold storage cost. Each executable has independent storage.
    # Size follows cuBLAS guidance: 4 MiB before sm90, 32 MiB from sm90 onward.
    workspace = CUDA.CuArray{UInt8,1,CUDA.DeviceMemory}(
        undef, CUDA.capability(context.device) >= v"9.0" ? 32 << 20 : 4 << 20,
    )
    graph = _with_prepared_cublas_workspace(context, workspace) do
        # CUDA 6.5 otherwise creates a separate capture stream. Keep capture
        # on the stream whose cuBLAS handle owns the retained workspace.
        @static if Base.pkgversion(CUDA) >= v"6.5"
            CUDA.capture(f; stream=context.stream)
        else
            CUDA.capture(f)
        end
    end
    _require_cuda_graph_fixed_storage(graph)
    executable = CUDA.instantiate(graph)
    return CUDAPreparedDeviceGraph(graph, executable, workspace)
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
