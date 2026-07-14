#
# Simple discrete-time controller models
#
# The current core controller is a leaky/discrete integrator followed by a
# first-order DM lag state.
#
# For each update:
# 1. integrate the incoming command-like input with gain `gain`
# 2. relax the DM state toward that integrated command with time constant `tau`
#
# This is intentionally lightweight and serves as a reusable controller model
# for examples, tutorial transfer-function reasoning, and simple closed-loop
# tests.
#
abstract type AbstractController end

"""
    controller_output(ctrl)

Return the current command-like output state of a controller.
"""
function controller_output(::AbstractController)
    throw(InvalidConfiguration("controller_output is not defined for this controller family"))
end

"""
    reset_controller!(ctrl)

Reset any internal controller state for a fresh control run.
"""
function reset_controller!(::AbstractController)
    throw(InvalidConfiguration("reset_controller! is not defined for this controller family"))
end

supports_controller_reset(::AbstractController) = false
runtime_controller_storage(::AbstractController) = nothing
runtime_controller_ownership_roots(controller::AbstractController) =
    (controller,)

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

@inline controller_output(ctrl::DiscreteIntegratorController) = ctrl.dm_state
@inline supports_controller_reset(::DiscreteIntegratorController) = true
@inline runtime_controller_storage(ctrl::DiscreteIntegratorController) =
    (ctrl.i_state, ctrl.dm_state)

function reset_controller!(ctrl::DiscreteIntegratorController)
    fill!(ctrl.i_state, zero(eltype(ctrl.i_state)))
    fill!(ctrl.dm_state, zero(eltype(ctrl.dm_state)))
    return ctrl
end

function DiscreteIntegratorController(n::Int; gain::Real=0.3, tau::Real=0.02, T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    n > 0 || throw(InvalidConfiguration("controller length must be greater than zero"))
    tau_value = T(tau)
    tau_value > zero(T) ||
        throw(InvalidConfiguration("controller tau must be greater than zero at controller precision"))
    backend = _resolve_array_backend(backend)
    i_state = backend{T}(undef, n)
    dm_state = backend{T}(undef, n)
    fill!(i_state, zero(T))
    fill!(dm_state, zero(T))
    return DiscreteIntegratorController{T, typeof(i_state)}(T(gain), tau_value, i_state, dm_state)
end

"""
    update!(ctrl, input, dt)

Advance the controller state by one sample period.

`i_state` integrates the incoming command-like input, and `dm_state` then
applies a first-order lag toward that integral state using the controller time
constant.
"""
function update!(ctrl::DiscreteIntegratorController, input::AbstractVector, dt::Real)
    length(input) == length(ctrl.i_state) ||
        throw(DimensionMismatchError("controller input length must match controller state"))
    T = eltype(ctrl.i_state)
    dt_value = T(dt)
    dt_value >= zero(T) ||
        throw(InvalidConfiguration("controller dt must be >= 0 at controller precision"))
    ctrl.i_state .+= ctrl.gain .* input .* dt_value
    ctrl.dm_state .+= (ctrl.i_state .- ctrl.dm_state) .* (dt_value / ctrl.tau)
    return controller_output(ctrl)
end
