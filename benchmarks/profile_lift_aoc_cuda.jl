using CUDA

include("profile_lift_aoc_common.jl")

CUDA.functional() || error("CUDA is not functional on this host")
CUDA.allowscalar(false)

cuda_device_allocated_bytes(f) = CUDA.@allocated f()
cuda_profile_region(f) = (CUDA.@profile external=true f(); nothing)

profile_lift_adapter("cuda", CUDA.CuArray, Backends.CUDABackend(), CUDA.synchronize,
    cuda_device_allocated_bytes, CUDA.used_memory, cuda_profile_region)
