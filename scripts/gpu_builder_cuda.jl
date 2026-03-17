using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_builder_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_builder_contract.jl"))
run_gpu_builder_smoke(CUDABackendTag)
