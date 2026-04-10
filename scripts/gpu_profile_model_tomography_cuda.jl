using AdaptiveOpticsSim
using Profile

try
    using CUDA
catch err
    error("gpu_profile_model_tomography_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_profile_model_tomography_contract.jl"))
run_gpu_model_tomography_profile(AdaptiveOpticsSim.CUDABackendTag)
