module AdaptiveOpticsSimAMDGPUExt

using AdaptiveOpticsSim
using AMDGPU

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.ROCArray
AdaptiveOpticsSim.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)

end
