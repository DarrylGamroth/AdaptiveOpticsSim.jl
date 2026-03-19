include(joinpath(@__DIR__, "backend_benchmark_contract.jl"))

try
    using CUDA
catch err
    error("benchmark_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

CUDA.functional() || error("benchmark_cuda.jl requires a functional CUDA driver/device")

run_backend_benchmark_suite(GPUBenchmarkTarget{CUDABackendTag}())
