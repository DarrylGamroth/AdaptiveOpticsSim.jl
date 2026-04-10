module AdaptiveOpticsSimMetalExt

using AdaptiveOpticsSim
using Metal

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.MetalBackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.MetalBackendTag}) = Metal.MtlArray
AdaptiveOpticsSim.gpu_backend_name(::Type{<:Metal.MtlArray}) = :metal
AdaptiveOpticsSim.array_backend_selector(::Type{<:Metal.MtlArray}) = AdaptiveOpticsSim.MetalBackend()
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.MetalBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = Metal.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.MetalBackendTag}, value, dims::Vararg{Int}) = Metal.fill(value, dims...)

end
