using AdaptiveOpticsSim

try
    using AMDGPU
catch err
    error("gpu_sync_audit_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_sync_audit_contract.jl"))
run_gpu_sync_audit(AdaptiveOpticsSim.AMDGPUBackendTag)
