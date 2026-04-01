include("runtests_head.jl")

include("testsets/core_optics.jl")
include("testsets/atmosphere.jl")
include("testsets/control_and_runtime.jl")
include("testsets/detectors_and_wfs.jl")
include("testsets/calibration_and_analysis.jl")
include("testsets/reference_and_tutorials.jl")

include("backend_optional_common.jl")
include("optional_amdgpu_backends.jl")
include("optional_cuda_backends.jl")
