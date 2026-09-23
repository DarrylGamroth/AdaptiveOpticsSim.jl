# This import must precede the shared harness loading AdaptiveOpticsSim. On
# Julia 1.12, loading AMDGPU after that package can leave extension methods at
# an earlier world age during the optional smoke test.
import AMDGPU
include("runtests_gpu_target_common.jl")
run_gpu_backend_target(AdaptiveOpticsSim.Backends.AMDGPUBackendTag)
