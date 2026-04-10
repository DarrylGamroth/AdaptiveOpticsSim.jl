using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_runtime_equivalence_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_runtime_equivalence_contract.jl"))
run_gpu_runtime_equivalence(AdaptiveOpticsSim.CUDABackendTag; branch_mode=BackendStreamExecution())
