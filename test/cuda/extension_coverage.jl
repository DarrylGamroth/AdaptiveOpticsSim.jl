using CUDA

include(joinpath(@__DIR__, "..", "backend_extension_coverage_common.jl"))

CUDA.functional() ||
    error("CUDA extension coverage requires a functional CUDA device")

@testset "CUDA extension coverage" begin
    run_backend_extension_coverage(Backends.CUDABackendTag)
end
