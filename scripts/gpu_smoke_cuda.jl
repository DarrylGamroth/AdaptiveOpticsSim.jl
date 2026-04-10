using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_smoke_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_smoke_contract.jl"))
run_gpu_smoke_matrix(AdaptiveOpticsSim.CUDABackendTag)
