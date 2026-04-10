using AdaptiveOpticsSim

try
    using AMDGPU
catch err
    error("gpu_hil_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_hil_contract.jl"))
run_gpu_hil_smoke(AdaptiveOpticsSim.AMDGPUBackendTag)
