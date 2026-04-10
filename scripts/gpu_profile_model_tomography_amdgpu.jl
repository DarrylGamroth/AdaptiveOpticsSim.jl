using AdaptiveOpticsSim
using Profile

try
    using AMDGPU
catch err
    error("gpu_profile_model_tomography_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_profile_model_tomography_contract.jl"))
run_gpu_model_tomography_profile(AdaptiveOpticsSim.AMDGPUBackendTag)
