#
# Simple discrete-time controller models
#
# The current core controller is a leaky/discrete integrator followed by a
# first-order DM lag state.
#
# For each update:
# 1. integrate the incoming slopes with gain `gain`
# 2. relax the DM state toward that integrated command with time constant `tau`
#
# This is intentionally lightweight and serves as a reusable controller model
# for examples, tutorial transfer-function reasoning, and simple closed-loop
# tests.
#
abstract type AbstractController end

"""
    DiscreteIntegratorController(n; gain=0.3, tau=0.02, T=Float64, backend::AbstractArrayBackend=CPUBackend())

Construct a simple discrete-time integral controller with a first-order DM lag
state.
"""
mutable struct DiscreteIntegratorController{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractController
    gain::T
    tau::T
    i_state::V
    dm_state::V
end

function DiscreteIntegratorController(n::Int; gain::Real=0.3, tau::Real=0.02, T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    i_state = backend{T}(undef, n)
    dm_state = backend{T}(undef, n)
    fill!(i_state, zero(T))
    fill!(dm_state, zero(T))
    return DiscreteIntegratorController{T, typeof(i_state)}(T(gain), T(tau), i_state, dm_state)
end

"""
    update!(ctrl, slopes, dt)

Advance the controller state by one sample period.

`i_state` integrates the incoming slopes, and `dm_state` then applies a
first-order lag toward that integral state using the controller time constant.
"""
function update!(ctrl::DiscreteIntegratorController, slopes::AbstractVector, dt::Real)
    ctrl.i_state .+= ctrl.gain .* slopes .* dt
    ctrl.dm_state .+= (ctrl.i_state .- ctrl.dm_state) .* (dt / ctrl.tau)
    return ctrl.dm_state
end
