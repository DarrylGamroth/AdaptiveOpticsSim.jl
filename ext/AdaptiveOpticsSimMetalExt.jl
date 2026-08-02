module AdaptiveOpticsSimMetalExt

import AdaptiveOpticsSim
import AdaptiveOpticsSim: Backends, Plant
using Metal

Backends.gpu_backend_loaded(::Type{Backends.MetalBackendTag}) = true
Backends.gpu_backend_array_type(::Type{Backends.MetalBackendTag}) = Metal.MtlArray
Backends.gpu_backend_name(::Type{<:Metal.MtlArray}) = :metal
Backends.array_backend_selector(::Type{<:Metal.MtlArray}) = Backends.MetalBackend()
Backends.backend_rand(::Type{Backends.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.rand(T, dims...)
Backends.backend_randn(::Type{Backends.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.randn(T, dims...)
Backends.backend_zeros(::Type{Backends.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.zeros(T, dims...)
Backends.backend_fill(::Type{Backends.MetalBackendTag}, value, dims::Vararg{Int}) = Metal.fill(value, dims...)
Backends.compute_device_identifier(array::Metal.MtlArray) =
    UInt(pointer(Metal.device(array)))

function Plant.structural_array_bytes(array::Metal.MtlArray,
    target::Backends.AbstractComputeDevice)
    Plant._require_structural_array_target(array, target)
    return Plant._contiguous_structural_array_bytes(array)
end

end
