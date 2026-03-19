include(joinpath(@__DIR__, "backend_benchmark_contract.jl"))
run_backend_benchmark_suite(CPUBenchmarkTarget())
