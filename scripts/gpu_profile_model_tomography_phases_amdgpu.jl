using AdaptiveOpticsSim
using LinearAlgebra
using SparseArrays

try
    using AMDGPU
catch err
    error("gpu_profile_model_tomography_phases_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_profile_model_tomography_phases_contract.jl"))
run_gpu_model_tomography_phase_profile(AMDGPUBackendTag)
