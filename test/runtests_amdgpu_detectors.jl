import AMDGPU
include("runtests_gpu_target_common.jl")
run_gpu_detector_target(AdaptiveOpticsSim.Backends.AMDGPUBackendTag)
