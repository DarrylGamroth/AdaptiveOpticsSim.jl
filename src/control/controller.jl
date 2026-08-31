#
# Simple discrete-time controller models
#
# The current core controller is a discrete integrator followed by a first-order
# DM lag state.
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
    VectorDelayLine(reference, delay_samples)

Preallocated fixed-length delay for vector-valued control signals. The delay
line owns `delay_samples` historical vectors and returns its internal scratch
vector from [`shift_delay!`](@ref).
"""
mutable struct VectorDelayLine{A<:AbstractMatrix,V<:AbstractVector}
    buffer::A
    scratch::V
    head::Int
end

function VectorDelayLine(reference::AbstractVector{T},
    delay_samples::Integer) where {T<:AbstractFloat}
    delay_samples >= 0 ||
        throw(InvalidConfiguration("delay_samples must be >= 0"))
    n_delay = Int(delay_samples)
    buffer = similar(reference, T, length(reference), n_delay)
    fill!(buffer, zero(T))
    scratch = similar(reference, T, length(reference))
    fill!(scratch, zero(T))
    return VectorDelayLine{typeof(buffer),typeof(scratch)}(
        buffer,
        scratch,
        1,
    )
end

"""
    shift_delay!(line, sample)

Insert `sample` and return the vector leaving the delay line. The returned
vector aliases `line.scratch` and remains valid only until the next call.
"""
function shift_delay!(line::VectorDelayLine, sample::AbstractVector)
    length(sample) == length(line.scratch) ||
        throw(DimensionMismatchError(
            "delay-line sample length must match its reference length"))
    n_delay = size(line.buffer, 2)
    if n_delay == 0
        copyto!(line.scratch, sample)
        return line.scratch
    end
    head = line.head
    copyto!(line.scratch, @view(line.buffer[:, head]))
    copyto!(@view(line.buffer[:, head]), sample)
    line.head = head == n_delay ? 1 : head + 1
    return line.scratch
end

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
    DiscreteIntegratorPlan(gain, tau)

Run-immutable coefficients for a discrete integrator followed by a first-order
DM lag. `tau` is expressed in seconds.
"""
struct DiscreteIntegratorPlan{T<:AbstractFloat}
    gain::T
    tau_s::T

    function DiscreteIntegratorPlan(gain::T, tau_s::T) where {T<:AbstractFloat}
        isfinite(gain) ||
            throw(InvalidConfiguration("controller gain must be finite"))
        isfinite(tau_s) && tau_s > zero(T) || throw(InvalidConfiguration(
            "controller tau must be finite and greater than zero",
        ))
        return new{T}(gain, tau_s)
    end
end

"""Persistent mathematical state for one discrete-integrator run."""
mutable struct DiscreteIntegratorState{T<:AbstractFloat,V<:AbstractVector{T}}
    integrated_command::V
    command::V
end

"""
Replaceable scratch for one discrete-integrator execution owner.

Recreating this workspace cannot change the controller trajectory.
"""
struct DiscreteIntegratorWorkspace{T<:AbstractFloat,V<:AbstractVector{T}}
    next_integrated_command::V
    next_command::V
end

function DiscreteIntegratorState(reference::AbstractVector{T}) where {T<:AbstractFloat}
    integrated_command = similar(reference, T, length(reference))
    command = similar(reference, T, length(reference))
    fill!(integrated_command, zero(T))
    fill!(command, zero(T))
    return DiscreteIntegratorState(integrated_command, command)
end

function DiscreteIntegratorWorkspace(reference::AbstractVector{T}) where {T<:AbstractFloat}
    next_integrated_command = similar(reference, T, length(reference))
    next_command = similar(reference, T, length(reference))
    fill!(next_integrated_command, zero(T))
    fill!(next_command, zero(T))
    return DiscreteIntegratorWorkspace(next_integrated_command, next_command)
end

"""
    DiscreteIntegratorController(n; gain=0.3, tau=0.02, T=Float64,
                                 backend::AbstractArrayBackend=CPUBackend())

Prepare a single-writer controller owner with a run-immutable plan, persistent
state, and replaceable workspace.
"""
struct DiscreteIntegratorController{P<:DiscreteIntegratorPlan,
    S<:DiscreteIntegratorState,W<:DiscreteIntegratorWorkspace} <: AbstractController
    plan::P
    state::S
    workspace::W
end

@inline controller_output(ctrl::DiscreteIntegratorController) = ctrl.state.command
@inline supports_controller_reset(::DiscreteIntegratorController) = true
@inline runtime_controller_storage(ctrl::DiscreteIntegratorController) =
    (
        ctrl.state.integrated_command,
        ctrl.state.command,
        ctrl.workspace.next_integrated_command,
        ctrl.workspace.next_command,
    )
@inline runtime_controller_ownership_roots(ctrl::DiscreteIntegratorController) =
    (ctrl.state, ctrl.workspace)

@inline discrete_integrator_plan(ctrl::DiscreteIntegratorController) = ctrl.plan
@inline discrete_integrator_state(ctrl::DiscreteIntegratorController) = ctrl.state
@inline discrete_integrator_workspace(ctrl::DiscreteIntegratorController) = ctrl.workspace

function reset_controller!(
    state::DiscreteIntegratorState,
    workspace::DiscreteIntegratorWorkspace,
)
    fill!(state.integrated_command, zero(eltype(state.integrated_command)))
    fill!(state.command, zero(eltype(state.command)))
    fill!(workspace.next_integrated_command, zero(eltype(workspace.next_integrated_command)))
    fill!(workspace.next_command, zero(eltype(workspace.next_command)))
    return state
end

function reset_controller!(ctrl::DiscreteIntegratorController)
    reset_controller!(ctrl.state, ctrl.workspace)
    return ctrl
end

function DiscreteIntegratorController(n::Int; gain::Real=0.3, tau::Real=0.02,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    n > 0 || throw(InvalidConfiguration("controller length must be greater than zero"))
    backend = _resolve_array_backend(backend)
    reference = backend{T}(undef, n)
    fill!(reference, zero(T))
    plan = DiscreteIntegratorPlan(T(gain), T(tau))
    command = similar(reference)
    fill!(command, zero(T))
    state = DiscreteIntegratorState(reference, command)
    workspace = DiscreteIntegratorWorkspace(reference)
    return DiscreteIntegratorController(plan, state, workspace)
end

function _prepare_controller_update!(
    state::DiscreteIntegratorState{T},
    workspace::DiscreteIntegratorWorkspace{T},
    plan::DiscreteIntegratorPlan{T},
    input::AbstractVector,
    dt::Real,
) where {T}
    axes(input) == axes(state.integrated_command) == axes(state.command) ==
        axes(workspace.next_integrated_command) == axes(workspace.next_command) ||
        throw(DimensionMismatchError(
            "controller input, state, and workspace axes must match",
        ))
    dt_s = T(dt)
    isfinite(dt_s) && dt_s >= zero(T) || throw(InvalidConfiguration(
        "controller dt must be finite and nonnegative at controller precision",
    ))
    lag = dt_s / plan.tau_s
    @. workspace.next_integrated_command =
        state.integrated_command + plan.gain * input * dt_s
    @. workspace.next_command =
        state.command + (workspace.next_integrated_command - state.command) * lag

    return nothing
end

@inline function _commit_controller_update!(
    state::DiscreteIntegratorState,
    workspace::DiscreteIntegratorWorkspace,
)
    copyto!(state.integrated_command, workspace.next_integrated_command)
    copyto!(state.command, workspace.next_command)
    return state.command
end

function update!(
    state::DiscreteIntegratorState,
    workspace::DiscreteIntegratorWorkspace,
    plan::DiscreteIntegratorPlan,
    input::AbstractVector,
    dt::Real,
)
    _prepare_controller_update!(state, workspace, plan, input, dt)
    return _commit_controller_update!(state, workspace)
end

function update!(
    output::AbstractVector,
    state::DiscreteIntegratorState,
    workspace::DiscreteIntegratorWorkspace,
    plan::DiscreteIntegratorPlan,
    input::AbstractVector,
    dt::Real,
)
    axes(output) == axes(state.command) || throw(DimensionMismatchError(
        "controller output and state axes must match",
    ))
    _prepare_controller_update!(state, workspace, plan, input, dt)
    copyto!(output, workspace.next_command)
    _commit_controller_update!(state, workspace)
    return output
end

"""
    update!(ctrl, input, dt)

Advance the controller state by one sample period.

The persistent `integrated_command` integrates the incoming command-like input,
and `command` then applies a first-order lag toward that integral state using
the controller time constant. Replacement workspace is prepared before either
state value commits.
"""
function update!(ctrl::DiscreteIntegratorController, input::AbstractVector, dt::Real)
    return update!(ctrl.state, ctrl.workspace, ctrl.plan, input, dt)
end
