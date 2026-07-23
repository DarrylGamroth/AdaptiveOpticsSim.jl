import CUDA
include("runtests_gpu_target_common.jl")
run_gpu_backend_target(AdaptiveOpticsSim.CUDABackendTag)
