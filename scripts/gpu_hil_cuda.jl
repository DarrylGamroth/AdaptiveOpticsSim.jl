using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_hil_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_hil_contract.jl"))
run_gpu_hil_smoke(CUDABackendTag)
