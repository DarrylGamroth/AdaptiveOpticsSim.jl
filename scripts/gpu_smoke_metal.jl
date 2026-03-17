using AdaptiveOpticsSim

try
    using Metal
catch err
    error("gpu_smoke_metal.jl requires Metal.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_smoke_contract.jl"))
run_gpu_smoke_matrix(MetalBackendTag)
