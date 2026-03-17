using AdaptiveOpticsSim

try
    using AMDGPU
catch err
    error("gpu_smoke_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_smoke_contract.jl"))
run_gpu_smoke_matrix(AMDGPUBackendTag)
