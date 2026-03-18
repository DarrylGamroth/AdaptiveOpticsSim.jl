module AdaptiveOpticsSimAMDGPUExt

using AdaptiveOpticsSim
using AMDGPU
using AbstractFFTs

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.ROCArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = :amdgpu
AdaptiveOpticsSim.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AMDGPU.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)
AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AMDGPU.rocFFT.ROCFFTPlan) = (plan * buffer; buffer)
AdaptiveOpticsSim.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AbstractFFTs.ScaledPlan) = (plan * buffer; buffer)
AdaptiveOpticsSim.default_build_backend(::AMDGPU.ROCArray) = AdaptiveOpticsSim.GPUArrayBuildBackend(AdaptiveOpticsSim.AMDGPUBackendTag)
AdaptiveOpticsSim.prepare_build_matrix(::AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.AMDGPUBackendTag}, A::AbstractMatrix) = Matrix(A)

function AdaptiveOpticsSim.solve_lift_system!(diag::AdaptiveOpticsSim.LiFTDiagnostics{T},
    residual::AMDGPU.ROCArray{T,1}, rhs::AMDGPU.ROCArray{T,1}, H::AbstractMatrix{T},
    normal::AbstractMatrix{T}, factor::AbstractMatrix{T},
    effective_mode::AdaptiveOpticsSim.LiFTSolveMode, damping::AdaptiveOpticsSim.LiFTDampingMode) where {T<:AbstractFloat}
    host_residual = Vector(residual)
    host_rhs = Vector(rhs)
    host_H = Matrix(H)
    host_normal = Matrix(normal)
    host_factor = Matrix(factor)
    host_delta = AdaptiveOpticsSim.solve_lift_system!(diag, host_residual, host_rhs, host_H, host_normal, host_factor, effective_mode, damping)
    if effective_mode isa AdaptiveOpticsSim.LiFTSolveQR
        copyto!(residual, host_delta)
        return residual
    end
    copyto!(rhs, host_delta)
    return rhs
end

end
