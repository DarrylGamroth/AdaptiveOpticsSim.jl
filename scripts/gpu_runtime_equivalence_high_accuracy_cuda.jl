using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_runtime_equivalence_high_accuracy_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_runtime_equivalence_contract.jl"))
run_gpu_runtime_equivalence_high_accuracy(CUDABackendTag; branch_mode=BackendStreamBranchExecution())
