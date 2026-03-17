using AdaptiveOpticsSim

if !isdefined(@__MODULE__, :run_gpu_smoke_matrix)
    include(joinpath(@__DIR__, "gpu_smoke_contract.jl"))
end

if !isdefined(@__MODULE__, :run_gpu_builder_smoke)
    include(joinpath(@__DIR__, "gpu_builder_contract.jl"))
end

function run_gpu_hil_smoke(::Type{B}) where {B<:GPUBackendTag}
    Base.invokelatest(run_gpu_smoke_matrix, B)
    Base.invokelatest(run_gpu_builder_smoke, B)
    println("gpu_hil_smoke complete")
    return nothing
end
