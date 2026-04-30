include("gpu_target_env.jl")
import CUDA
include(joinpath(dirname(@__DIR__), "scripts", "gpu_runtime_equivalence_contract.jl"))
include("runtests_gpu_target_common.jl")
run_gpu_backend_target(AdaptiveOpticsSim.CUDABackendTag)
