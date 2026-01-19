abstract type AbstractController end

mutable struct DiscreteIntegratorController{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractController
    gain::T
    tau::T
    i_state::V
    dm_state::V
end

function DiscreteIntegratorController(n::Int; gain::Real=0.3, tau::Real=0.02, T=Float64, backend=Array)
    i_state = backend{T}(undef, n)
    dm_state = backend{T}(undef, n)
    fill!(i_state, zero(T))
    fill!(dm_state, zero(T))
    return DiscreteIntegratorController{T, typeof(i_state)}(T(gain), T(tau), i_state, dm_state)
end

function update!(ctrl::DiscreteIntegratorController, slopes::AbstractVector, dt::Real)
    ctrl.i_state .+= ctrl.gain .* slopes .* dt
    ctrl.dm_state .+= (ctrl.i_state .- ctrl.dm_state) .* (dt / ctrl.tau)
    return ctrl.dm_state
end
