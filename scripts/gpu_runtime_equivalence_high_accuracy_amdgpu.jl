using AdaptiveOpticsSim

try
    using AMDGPU
catch err
    error("gpu_runtime_equivalence_high_accuracy_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_runtime_equivalence_contract.jl"))
run_gpu_runtime_equivalence_high_accuracy(AMDGPUBackendTag; branch_mode=SequentialExecution())
