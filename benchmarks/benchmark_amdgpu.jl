include(joinpath(@__DIR__, "backend_benchmark_contract.jl"))

try
    using AMDGPU
catch err
    error("benchmark_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

AMDGPU.functional() || error("benchmark_amdgpu.jl requires a functional ROCm installation and GPU")

run_backend_benchmark_suite(GPUBenchmarkTarget{AMDGPUBackendTag}())
