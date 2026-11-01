# Optional AMDGPU selector. CPU graphs here are numerical references, never a
# simulator fallback. Run with one Julia/BLAS thread and the selected AMD target.
include("proper/compiled_paths_helpers.jl")
using AMDGPU
import AdaptiveOpticsSim.Backends
AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const COMPILED_TARGET=Backends.AcceleratorComputeDevice(Backends.AMDGPUBackend(),1)

const COMPILED_ARRAY=ROCArray
const COMPILED_SYNCHRONIZE=AMDGPU.synchronize
include("proper/compiled_paths_gpu.jl")
@test !any(id->id.name=="CUDA",keys(Base.loaded_modules))
