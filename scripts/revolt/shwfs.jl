include(joinpath(@__DIR__, "common.jl"))

revolt_setup(; T::Type{<:AbstractFloat}=Float64, backend=Array, response_mode::Symbol=:default) =
    revolt_setup(:shwfs; T=T, backend=backend, response_mode=response_mode)
