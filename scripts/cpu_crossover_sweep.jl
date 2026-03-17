using AdaptiveOpticsSim

include(joinpath(@__DIR__, "backend_crossover_contract.jl"))
run_backend_crossover_sweep(CPUSweepTarget())
