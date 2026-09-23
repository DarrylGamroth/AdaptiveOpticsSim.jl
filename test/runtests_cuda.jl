# Keep optional backend loading ahead of the shared harness, which loads
# AdaptiveOpticsSim and its CUDA extension.
import CUDA
include("runtests_gpu_target_common.jl")
run_gpu_backend_target(AdaptiveOpticsSim.Backends.CUDABackendTag)
