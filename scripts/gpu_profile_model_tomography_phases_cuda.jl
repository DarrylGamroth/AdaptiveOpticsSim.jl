using AdaptiveOpticsSim
using LinearAlgebra
using SparseArrays

try
    using CUDA
catch err
    error("gpu_profile_model_tomography_phases_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_profile_model_tomography_phases_contract.jl"))
run_gpu_model_tomography_phase_profile(CUDABackendTag)
