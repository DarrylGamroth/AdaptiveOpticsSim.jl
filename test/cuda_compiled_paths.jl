# Optional CUDA selector. CPU execution is a numerical oracle only.
include("proper/compiled_paths_helpers.jl")
using CUDA
import AdaptiveOpticsSim.Backends
CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const COMPILED_TARGET=Backends.AcceleratorComputeDevice(Backends.CUDABackend(),0)
const COMPILED_ARRAY=CuArray
const COMPILED_SYNCHRONIZE=CUDA.synchronize
include("proper/compiled_paths_gpu.jl")
