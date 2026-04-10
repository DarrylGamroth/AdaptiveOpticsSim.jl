using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_crossover_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "backend_crossover_contract.jl"))
run_backend_crossover_sweep(GPUSweepTarget{AdaptiveOpticsSim.CUDABackendTag}())
