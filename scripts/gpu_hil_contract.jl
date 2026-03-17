using AdaptiveOpticsSim

function run_gpu_hil_smoke(::Type{B}) where {B<:GPUBackendTag}
    include(joinpath(@__DIR__, "gpu_smoke_contract.jl"))
    include(joinpath(@__DIR__, "gpu_builder_contract.jl"))
    run_gpu_smoke_matrix(B)
    run_gpu_builder_smoke(B)
    println("gpu_hil_smoke complete")
    return nothing
end
