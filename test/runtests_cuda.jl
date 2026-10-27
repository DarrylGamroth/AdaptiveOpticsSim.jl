# Keep optional backend loading ahead of the shared harness, which loads
# AdaptiveOpticsSim and its CUDA extension.
import CUDA
include("runtests_gpu_target_common.jl")
run_gpu_backend_target(AdaptiveOpticsSim.Backends.CUDABackendTag)
include("cuda/captured_storage.jl")
include("cuda/prepared_context.jl")
include("cuda/prepared_wait.jl")
include("cuda/optical_quantities.jl")
include("cuda/geometric_slopes.jl")
include("cuda/kernel_completion.jl")
include("cuda/centering_precision.jl")

include("cuda/hil_boundary.jl")
