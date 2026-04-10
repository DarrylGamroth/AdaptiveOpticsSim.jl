#
# Closed-loop runtime execution and simulation I/O
#
# This file defines the low-latency control loop used by HIL and deterministic
# runtime benchmarks.
#
# The core runtime step is:
# 1. advance and propagate the atmosphere
# 2. apply the current DM command to the telescope OPD
# 3. measure WFS outputs, optionally through a detector
# 4. reconstruct a new command from the sensed slopes
# 5. apply the signed command back to the DM coefficient vector
#
# `SimulationInterface` and `SimulationReadout` expose that state as the
# controller-facing I/O surface without changing the underlying runtime object.
#
@kernel function apply_command_kernel!(coefs, cmd, sign, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds coefs[i] = sign * cmd[i]
    end
end

abstract type AbstractControlSimulation end
abstract type AbstractExecutionPolicy end
abstract type AbstractRuntimeProfile end

struct SequentialExecution <: AbstractExecutionPolicy end
struct ThreadedExecution <: AbstractExecutionPolicy end
struct BackendStreamExecution <: AbstractExecutionPolicy end

struct ScientificRuntimeProfile <: AbstractRuntimeProfile end
struct HILRuntimeProfile <: AbstractRuntimeProfile end

struct RuntimeLatencyModel
    measurement_delay_frames::Int
    readout_delay_frames::Int
    reconstruction_delay_frames::Int
    dm_delay_frames::Int
end

default_runtime_profile() = ScientificRuntimeProfile()

function RuntimeLatencyModel(;
    measurement_delay_frames::Integer=0,
    readout_delay_frames::Integer=0,
    reconstruction_delay_frames::Integer=0,
    dm_delay_frames::Integer=0,
)
    delays = (
        Int(measurement_delay_frames),
        Int(readout_delay_frames),
        Int(reconstruction_delay_frames),
        Int(dm_delay_frames),
    )
    all(>=(0), delays) || throw(InvalidConfiguration("runtime latency delays must be >= 0"))
    return RuntimeLatencyModel(delays...)
end

default_runtime_latency(::ScientificRuntimeProfile) = RuntimeLatencyModel()
default_runtime_latency(::HILRuntimeProfile) = RuntimeLatencyModel(
    measurement_delay_frames=1,
    readout_delay_frames=1,
    reconstruction_delay_frames=0,
    dm_delay_frames=1,
)

runtime_science_zero_padding(::ScientificRuntimeProfile) = 2
runtime_science_zero_padding(::HILRuntimeProfile) = 0

"""
    VectorDelayLine(ref, delay_frames)

Preallocated fixed-length frame delay for vector-valued control signals.

The delay line stores `delay_frames` historical samples and returns the signal
that exits the tail when a new sample is shifted in.
"""
mutable struct VectorDelayLine{A<:AbstractMatrix,V<:AbstractVector}
    buffer::A
    scratch::V
end

function VectorDelayLine(ref::AbstractVector{T}, delay_frames::Int) where {T<:AbstractFloat}
    delay_frames >= 0 || throw(InvalidConfiguration("delay_frames must be >= 0"))
    buffer = similar(ref, T, length(ref), delay_frames)
    fill!(buffer, zero(T))
    scratch = similar(ref, T, length(ref))
    fill!(scratch, zero(T))
    return VectorDelayLine{typeof(buffer),typeof(scratch)}(buffer, scratch)
end

function shift_delay!(line::VectorDelayLine, sample::AbstractVector)
    n_delay = size(line.buffer, 2)
    if n_delay == 0
        copyto!(line.scratch, sample)
        return line.scratch
    end
    copyto!(line.scratch, @view(line.buffer[:, 1]))
    @inbounds for i in 1:n_delay-1
        copyto!(@view(line.buffer[:, i]), @view(line.buffer[:, i + 1]))
    end
    copyto!(@view(line.buffer[:, n_delay]), sample)
    return line.scratch
end

prepare!(sim::AbstractControlSimulation) = sim
supports_prepared_runtime(::Type) = false
supports_prepared_runtime(sim) = supports_prepared_runtime(typeof(sim))
supports_detector_output(::Type) = false
supports_detector_output(sim) = supports_detector_output(typeof(sim))
supports_stacked_sources(::Type) = false
supports_stacked_sources(sim) = supports_stacked_sources(typeof(sim))
supports_grouped_execution(::Type) = false
supports_grouped_execution(sim) = supports_grouped_execution(typeof(sim))

"""
    supports_prepared_runtime(wfs, src)

Return whether a WFS/source pairing exposes a meaningful `prepare_runtime_wfs!`
precomputation hook for repeated runtime stepping.
"""
supports_prepared_runtime(::AbstractWFS, ::Any) = false

"""
    supports_detector_output(wfs, det)

Return whether the WFS exposes a maintained detector-coupled `measure!`
surface for the supplied detector family.
"""
supports_detector_output(::AbstractWFS, ::AbstractDetector) = false

"""
    supports_stacked_sources(wfs, src)

Return whether the WFS/source pairing has maintained stacked-source execution
support, typically for `Asterism` inputs.
"""
supports_stacked_sources(::AbstractWFS, ::Any) = false

"""
    supports_grouped_execution(wfs, src)

Return whether the WFS/source pairing exposes maintained grouped execution,
such as grouped asterism or spectral accumulation.
"""
supports_grouped_execution(::AbstractWFS, ::Any) = false

init_execution_state(::AbstractExecutionPolicy, ref) = nothing

"""
    ClosedLoopRuntime(simulation, reconstructor; ...)

Create the executable closed-loop control runtime for an AO simulation.

The runtime owns the mutable command vector, optional detector outputs, and the
objects needed to execute repeated `sense!` and `step!` calls without
reconstructing the simulation graph.
"""
mutable struct ClosedLoopRuntime{
    SIM<:AOSimulation,
    TEL,
    A,
    S,
    O,
    W,
    R,
    CV,
    SV,
    WD,
    SD,
    RNG,
    RP<:AbstractRuntimeProfile,
    PR<:RuntimeProductRequirements,
    MD,
    RD,
    CD,
    DD,
    CL,
    T<:AbstractFloat,
} <: AbstractControlSimulation
    simulation::SIM
    tel::TEL
    atm::A
    src::S
    optic::O
    wfs::W
    reconstructor::R
    command::CV
    reconstruct_buffer::CV
    slopes::SV
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    profile::RP
    products::PR
    latency::RuntimeLatencyModel
    measurement_delay::MD
    readout_delay::RD
    reconstruction_delay::CD
    dm_delay::DD
    command_layout::CL
    control_sign::T
    science_zero_padding::Int
    prepared::Bool
end

@inline function Base.getproperty(runtime::ClosedLoopRuntime, name::Symbol)
    if name === :dm
        return getfield(runtime, :optic)
    end
    return getfield(runtime, name)
end

@inline function Base.propertynames(runtime::ClosedLoopRuntime, private::Bool=false)
    names = fieldnames(typeof(runtime))
    return :dm in names ? names : (names..., :dm)
end

"""
    SimulationInterface(runtime)

Allocate a controller-facing I/O view for a runtime.

The interface owns copied command, slope, and frame buffers so external code
can snapshot or stream runtime state without mutating the runtime internals.
"""
mutable struct SimulationInterface{RT,C,S,W,SF}
    runtime::RT
    command::C
    slopes::S
    wfs_frame::W
    science_frame::SF
end

"""
    CompositeSimulationInterface(interfaces...)

Aggregate several compatible simulation interfaces into one combined command,
slope, and frame surface.

This is used for grouped multi-branch execution while preserving per-runtime
state internally.
"""
mutable struct CompositeSimulationInterface{IT,C,S,WF,SF,PR,WS,SS}
    interfaces::IT
    command::C
    slopes::S
    wfs_frames::WF
    science_frames::SF
    products::PR
    wfs_stack::WS
    science_stack::SS
end

"""
    SimulationReadout

Zero-copy view of the current output side of a simulation runtime or interface.

This bundles the command, slopes, optional frames, and export metadata that an
external controller or logger would observe at the current simulation state.
"""
struct SimulationReadout{C,S,W,SF,WM,SM,GW,GS}
    command::C
    slopes::S
    wfs_frame::W
    science_frame::SF
    wfs_metadata::WM
    science_metadata::SM
    grouped_wfs_stack::GW
    grouped_science_stack::GS
end

@inline function apply_command!(::ScalarCPUStyle, coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    @inbounds for i in eachindex(coefs, cmd)
        coefs[i] = sign * cmd[i]
    end
    return coefs
end

@inline function apply_command!(style::AcceleratorStyle, coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    launch_kernel!(style, apply_command_kernel!, coefs, cmd, sign, length(coefs); ndrange=length(coefs))
    return coefs
end

@inline function apply_command!(coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    apply_command!(execution_style(coefs), coefs, cmd, sign)
end
