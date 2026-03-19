using AdaptiveOpticsSim

try
    using AMDGPU
catch err
    error("gpu_crossover_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "backend_crossover_contract.jl"))
run_backend_crossover_sweep(GPUSweepTarget{AMDGPUBackendTag}())
