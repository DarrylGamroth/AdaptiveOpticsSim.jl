module AdaptiveOpticsSimCUDAExt

using AdaptiveOpticsSim
using CUDA

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.CUDABackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.CuArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = :cuda
AdaptiveOpticsSim.gpu_backend_name(::Type{<:CUDA.CuArray}) = :cuda
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.CUDABackendTag}, value, dims::Vararg{Int}) = CUDA.fill(value, dims...)

struct AO188CUDAStreamState
    high::CUDA.CuStream
    low::CUDA.CuStream
end

function AdaptiveOpticsSim.init_branch_execution_state(::AdaptiveOpticsSim.BackendStreamBranchExecution, ::CUDA.CuArray)
    return AO188CUDAStreamState(
        CUDA.CuStream(; flags=CUDA.STREAM_NON_BLOCKING),
        CUDA.CuStream(; flags=CUDA.STREAM_NON_BLOCKING),
    )
end

function AdaptiveOpticsSim.measure_branches_backend!(::AdaptiveOpticsSim.BackendStreamBranchExecution,
    surrogate::AdaptiveOpticsSim.AO1883kSurrogate, state::AO188CUDAStreamState)
    high_rng, low_rng = AdaptiveOpticsSim._frame_rngs!(surrogate.rng)
    high_task = Threads.@spawn CUDA.stream!(state.high) do
        AdaptiveOpticsSim._measure_high!(surrogate, high_rng)
        CUDA.synchronize(state.high)
    end
    low_task = Threads.@spawn CUDA.stream!(state.low) do
        AdaptiveOpticsSim._measure_low!(surrogate, low_rng)
        CUDA.synchronize(state.low)
    end
    fetch(high_task)
    fetch(low_task)
    return surrogate
end

end
