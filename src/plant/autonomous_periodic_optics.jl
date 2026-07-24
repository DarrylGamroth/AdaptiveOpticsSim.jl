#
# Trigger-relative autonomous periodic optical devices
#
# Fast local waveform generators are command-responsive physical devices, but
# their external command surface carries bounded setpoints rather than one
# sample per waveform point. The first native implementation binds an analytic
# circular Pyramid modulator to the existing normalized cycle quadrature.
#

abstract type AbstractWaveformPhaseReference end

"""
Plant-time phase origin for an autonomous waveform that does not consume a
modeled trigger edge.
"""
struct FreeRunningPhaseReference <: AbstractWaveformPhaseReference
    origin::PlantTimestamp
end

FreeRunningPhaseReference(;
    origin::PlantTimestamp=zero(PlantTimestamp)) =
    FreeRunningPhaseReference(origin)

"""
Phase reference synchronized to realization of a modeled trigger source.

This relationship represents a waveform generator whose own cycle reference
also feeds a trigger fan-out. `initial_reference` defines phase before the
first realized source edge.
"""
struct TriggerSourcePhaseReference <: AbstractWaveformPhaseReference
    source::TriggerSourceID
    initial_reference::PlantTimestamp
end

TriggerSourcePhaseReference(source;
    initial_reference::PlantTimestamp=zero(PlantTimestamp)) =
    TriggerSourcePhaseReference(_as_trigger_source_id(source),
        initial_reference)

"""
Phase reset driven by one delivered trigger consumer.

Link delay, skew, jitter, drops, and duplicates therefore affect this waveform
branch independently. `initial_reference` defines phase before the first
delivered reset.
"""
struct TriggerResetPhaseReference <: AbstractWaveformPhaseReference
    consumer::TriggerConsumerID
    initial_reference::PlantTimestamp
end

TriggerResetPhaseReference(consumer;
    initial_reference::PlantTimestamp=zero(PlantTimestamp)) =
    TriggerResetPhaseReference(_as_trigger_consumer_id(consumer),
        initial_reference)

abstract type AbstractAutonomousOpticFidelity end

"""
One normalized optical quadrature over a complete periodic cycle at one frozen
optical epoch. Detector acquisition owns exposure integration.
"""
struct CycleAveragedModulationFidelity <: AbstractAutonomousOpticFidelity end

"""
    AutonomousPeriodicOpticDefinition(optic, path;
        phase_reference, fidelity=CycleAveragedModulationFidelity())

Bind one independently commanded controllable optic to one optical path and
one immutable waveform-fidelity/phase-reference contract. The definition owns
no waveform cursor, trigger state, command storage, or mutable optical
workspace.
"""
struct AutonomousPeriodicOpticDefinition{
    R<:AbstractWaveformPhaseReference,
    F<:AbstractAutonomousOpticFidelity,
}
    optic::ControllableOpticID
    path::OpticalPathID
    phase_reference::R
    fidelity::F
end

function AutonomousPeriodicOpticDefinition(optic, path;
    phase_reference::AbstractWaveformPhaseReference,
    fidelity::AbstractAutonomousOpticFidelity=
        CycleAveragedModulationFidelity())
    return AutonomousPeriodicOpticDefinition(
        _as_controllable_optic_id(optic),
        _as_optical_path_id(path),
        phase_reference,
        fidelity,
    )
end

"""
    CircularPyramidModulator(; radius_endpoint, frequency_endpoint,
        phase_endpoint, enabled_endpoint)

Cold endpoint-role declaration for a local circular Pyramid modulation
generator. Radius, frequency, phase, and enable state remain independent
absolute scalar command endpoints. The cycle-averaged optical implementation
uses radius and enable state; frequency and phase retain the analytic
trigger-relative device phase without changing full-cycle quadrature.

Additional mode or dither semantics are model-specific extensions rather than
an encoded universal command vector.
"""
struct CircularPyramidModulator
    radius_endpoint::CommandEndpointID
    frequency_endpoint::CommandEndpointID
    phase_endpoint::CommandEndpointID
    enabled_endpoint::CommandEndpointID
end

function CircularPyramidModulator(; radius_endpoint, frequency_endpoint,
    phase_endpoint, enabled_endpoint)
    endpoints = (
        _as_command_endpoint_id(radius_endpoint),
        _as_command_endpoint_id(frequency_endpoint),
        _as_command_endpoint_id(phase_endpoint),
        _as_command_endpoint_id(enabled_endpoint),
    )
    @inbounds for right in 2:length(endpoints), left in 1:(right - 1)
        endpoints[left] == endpoints[right] && throw(PlantDefinitionError(
            :autonomous_periodic_optic, :duplicate_endpoint,
            "circular Pyramid modulation endpoint roles must be distinct"))
    end
    return CircularPyramidModulator(endpoints...)
end

plant_model_definition_style(::Type{CircularPyramidModulator}) =
    ColdPlantModelDefinition()

abstract type AbstractControllableOpticExecutionRole end
struct PupilSurfaceExecutionRole <: AbstractControllableOpticExecutionRole end
struct AutonomousPathExecutionRole <: AbstractControllableOpticExecutionRole end

"""
Execution-role trait for a prepared controllable optic. The default applies a
surface to each selected common pupil input. Path-local autonomous devices opt
in explicitly and are evaluated only through their prepared path binding.
"""
controllable_optic_execution_role(::Any) = PupilSurfaceExecutionRole()

struct PreparedCircularPyramidModulator{
    T<:AbstractFloat,E<:Integer,
}
    radius_endpoint::CommandEndpointID
    frequency_endpoint::CommandEndpointID
    phase_endpoint::CommandEndpointID
    enabled_endpoint::CommandEndpointID
end

controllable_optic_execution_role(::PreparedCircularPyramidModulator) =
    AutonomousPathExecutionRole()

@inline function _autonomous_schema(
    definition::ControllableOpticDefinition, endpoint::CommandEndpointID)
    return command_schema(definition, endpoint)
end

function _require_autonomous_scalar_absolute_schema(
    schema::PlantCommandSchema, unit::Symbol,
    sign_convention::Symbol)
    isempty(command_dimensions(schema)) || throw(PlantPreparationError(
        :autonomous_periodic_optic, :command_shape,
        "autonomous waveform endpoint $(command_endpoint_id(schema)) must be scalar"))
    command_semantics(schema) == AbsoluteCommand || throw(
        PlantPreparationError(:autonomous_periodic_optic,
            :command_semantics,
            "autonomous waveform endpoint $(command_endpoint_id(schema)) must use absolute semantics"))
    command_units(schema) == CommandUnit(unit) || throw(
        PlantPreparationError(:autonomous_periodic_optic, :command_units,
            "autonomous waveform endpoint $(command_endpoint_id(schema)) must use $(repr(unit)) units"))
    command_sign_convention(schema) ==
        CommandSignConvention(sign_convention) || throw(
        PlantPreparationError(:autonomous_periodic_optic,
            :command_sign_convention,
            "autonomous waveform endpoint $(command_endpoint_id(schema)) has the wrong sign convention"))
    command_basis(schema) ==
        CommandBasis(:waveform_setpoint, :pyramid_modulation) || throw(
        PlantPreparationError(:autonomous_periodic_optic, :command_basis,
            "autonomous waveform endpoint $(command_endpoint_id(schema)) must use the pyramid-modulation waveform-setpoint basis"))
    return schema
end

function _require_bounded_autonomous_schema(
    schema::PlantCommandSchema,
    bounds::UniformCommandBounds)
    return schema
end

function _require_bounded_autonomous_schema(
    schema::PlantCommandSchema, ::UnboundedCommandValues)
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :unbounded_setpoint,
        "autonomous waveform endpoint $(command_endpoint_id(schema)) must declare finite bounds"))
end

function _require_nonnegative_autonomous_schema(
    schema::PlantCommandSchema)
    _require_bounded_autonomous_schema(schema, command_bounds(schema))
    bounds = command_bounds(schema)
    bounds.lower >= zero(bounds.lower) || throw(PlantPreparationError(
        :autonomous_periodic_optic, :negative_setpoint_bound,
        "autonomous waveform endpoint $(command_endpoint_id(schema)) must have a nonnegative lower bound"))
    return schema
end

function _require_enabled_autonomous_schema(schema::PlantCommandSchema)
    T = command_numeric_type(schema)
    T <: Integer || throw(PlantPreparationError(
        :autonomous_periodic_optic, :enabled_type,
        "autonomous waveform enable endpoint must use an integer scalar type"))
    _require_bounded_autonomous_schema(schema, command_bounds(schema))
    bounds = command_bounds(schema)
    bounds.lower == zero(T) && bounds.upper == one(T) || throw(
        PlantPreparationError(:autonomous_periodic_optic, :enabled_bounds,
            "autonomous waveform enable endpoint must be bounded to 0:1"))
    return schema
end

function prepare_controllable_optic(model::CircularPyramidModulator,
    definition::ControllableOpticDefinition,
    ::AbstractTelescope, ::AbstractAtmosphere)
    length(command_schemas(definition)) == 4 || throw(
        PlantPreparationError(:autonomous_periodic_optic,
            :command_schema_count,
            "CircularPyramidModulator requires exactly four command schemas"))
    radius = _require_nonnegative_autonomous_schema(
        _require_autonomous_scalar_absolute_schema(
            _autonomous_schema(definition, model.radius_endpoint),
            :lambda_over_d, :nonnegative_radius))
    frequency = _require_nonnegative_autonomous_schema(
        _require_autonomous_scalar_absolute_schema(
            _autonomous_schema(definition, model.frequency_endpoint),
            :hertz, :nonnegative_frequency))
    phase = _require_bounded_autonomous_schema(
        _require_autonomous_scalar_absolute_schema(
            _autonomous_schema(definition, model.phase_endpoint),
            :radian, :counterclockwise_phase),
        command_bounds(_autonomous_schema(definition,
            model.phase_endpoint)))
    enabled = _require_enabled_autonomous_schema(
        _require_autonomous_scalar_absolute_schema(
            _autonomous_schema(definition, model.enabled_endpoint),
            :binary_state, :zero_disabled_one_enabled))
    T = command_numeric_type(radius)
    T <: AbstractFloat || throw(PlantPreparationError(
        :autonomous_periodic_optic, :numeric_type,
        "circular Pyramid radius must use an AbstractFloat scalar type"))
    command_numeric_type(frequency) === T &&
        command_numeric_type(phase) === T || throw(
        PlantPreparationError(:autonomous_periodic_optic, :numeric_type,
            "circular Pyramid radius, frequency, and phase endpoints must use one exact floating-point type"))
    E = command_numeric_type(enabled)
    return PreparedCircularPyramidModulator{T,E}(
        model.radius_endpoint, model.frequency_endpoint,
        model.phase_endpoint, model.enabled_endpoint)
end

@inline function _initial_autonomous_command(endpoint_ids::Tuple,
    initial_commands::Tuple, target::CommandEndpointID)
    @inbounds for index in eachindex(endpoint_ids)
        endpoint_ids[index] == target && return initial_commands[index]
    end
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :prepared_binding,
        "autonomous waveform initial state is missing endpoint $target"))
end

mutable struct CircularPyramidModulatorState{T<:AbstractFloat}
    radius::T
    frequency_hz::T
    phase_offset_rad::T
    enabled::Bool
    phase_anchor_timestamp::PlantTimestamp
    phase_at_anchor_rad::T
    reference_timestamp::PlantTimestamp
    reference_sequence::UInt64
    reference_count::UInt64
    last_evaluation_timestamp::PlantTimestamp
    has_evaluation::Bool
    optical_dirty::Bool
end

mutable struct CircularPyramidModulatorWorkspace{
    T<:AbstractFloat,E<:Integer,
}
    staged_radius::T
    staged_frequency_hz::T
    staged_phase_offset_rad::T
    staged_enabled::E
end

function prepare_controllable_optic_state(
    prepared::PreparedCircularPyramidModulator{T,E},
    ::ControllableOpticDefinition, endpoint_ids::Tuple,
    initial_commands::Tuple) where {T,E}
    radius = T(_initial_autonomous_command(endpoint_ids, initial_commands,
        prepared.radius_endpoint))
    frequency = T(_initial_autonomous_command(endpoint_ids, initial_commands,
        prepared.frequency_endpoint))
    phase = T(_initial_autonomous_command(endpoint_ids, initial_commands,
        prepared.phase_endpoint))
    enabled_value = E(_initial_autonomous_command(endpoint_ids,
        initial_commands, prepared.enabled_endpoint))
    enabled_value in (zero(E), one(E)) || throw(PlantPreparationError(
        :autonomous_periodic_optic, :enabled_state,
        "initial autonomous waveform enable state must be zero or one"))
    origin = zero(PlantTimestamp)
    return CircularPyramidModulatorState(radius, frequency, phase,
        enabled_value == one(E), origin, zero(T), origin, UInt64(0),
        UInt64(0), origin, false, true)
end

function prepare_controllable_optic_workspace(
    ::PreparedCircularPyramidModulator{T,E}) where {T,E}
    return CircularPyramidModulatorWorkspace(
        zero(T), zero(T), zero(T), zero(E))
end

function stage_controllable_optic_command!(
    prepared::PreparedCircularPyramidModulator{T,E},
    ::CircularPyramidModulatorState{T},
    workspace::CircularPyramidModulatorWorkspace{T,E},
    endpoint::CommandEndpointID, command,
    ::PlantTimestamp) where {T,E}
    if endpoint == prepared.radius_endpoint
        workspace.staged_radius = T(command)
    elseif endpoint == prepared.frequency_endpoint
        workspace.staged_frequency_hz = T(command)
    elseif endpoint == prepared.phase_endpoint
        workspace.staged_phase_offset_rad = T(command)
    elseif endpoint == prepared.enabled_endpoint
        workspace.staged_enabled = E(command)
    else
        throw(PlantCommandError(:physical_application, :endpoint_mismatch,
            "circular Pyramid modulator received unknown endpoint $endpoint"))
    end
    return nothing
end

@inline function _unwrapped_waveform_phase(
    state::CircularPyramidModulatorState{T},
    timestamp::PlantTimestamp) where {T}
    elapsed_ns = plant_nanoseconds(timestamp) -
        plant_nanoseconds(state.phase_anchor_timestamp)
    elapsed_seconds = T(elapsed_ns) / T(1_000_000_000)
    return state.phase_at_anchor_rad +
        T(2pi) * state.frequency_hz * elapsed_seconds +
        state.phase_offset_rad
end

@inline function _waveform_phase_without_offset(
    state::CircularPyramidModulatorState{T},
    timestamp::PlantTimestamp) where {T}
    return _unwrapped_waveform_phase(state, timestamp) -
        state.phase_offset_rad
end

function commit_controllable_optic_command!(
    prepared::PreparedCircularPyramidModulator{T,E},
    state::CircularPyramidModulatorState{T},
    workspace::CircularPyramidModulatorWorkspace{T,E},
    endpoint::CommandEndpointID,
    timestamp::PlantTimestamp) where {T,E}
    if endpoint == prepared.radius_endpoint
        state.radius = workspace.staged_radius
        state.optical_dirty = true
    elseif endpoint == prepared.frequency_endpoint
        state.phase_at_anchor_rad =
            _waveform_phase_without_offset(state, timestamp)
        state.phase_anchor_timestamp = timestamp
        state.frequency_hz = workspace.staged_frequency_hz
    elseif endpoint == prepared.phase_endpoint
        state.phase_offset_rad = workspace.staged_phase_offset_rad
    elseif endpoint == prepared.enabled_endpoint
        state.enabled = workspace.staged_enabled == one(E)
        state.optical_dirty = true
    else
        throw(PlantCommandError(:physical_application, :endpoint_mismatch,
            "circular Pyramid modulator committed unknown endpoint $endpoint"))
    end
    return nothing
end

@inline _phase_reference_initial(reference::FreeRunningPhaseReference) =
    reference.origin
@inline _phase_reference_initial(reference::TriggerSourcePhaseReference) =
    reference.initial_reference
@inline _phase_reference_initial(reference::TriggerResetPhaseReference) =
    reference.initial_reference

function initialize_autonomous_periodic_optic!(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState,
    reference::AbstractWaveformPhaseReference)
    origin = _phase_reference_initial(reference)
    state.phase_anchor_timestamp = origin
    state.phase_at_anchor_rad = zero(state.phase_at_anchor_rad)
    state.reference_timestamp = origin
    state.reference_sequence = UInt64(0)
    state.reference_count = UInt64(0)
    state.last_evaluation_timestamp = origin
    state.has_evaluation = false
    state.optical_dirty = true
    return nothing
end

function reset_autonomous_periodic_optic_phase!(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState,
    timestamp::PlantTimestamp, sequence::UInt64)
    state.reference_count != typemax(UInt64) || throw(
        PlantScheduleError(:autonomous_periodic_optic,
            :reference_count_overflow,
            "autonomous waveform phase-reference count exceeds UInt64 range"))
    state.phase_anchor_timestamp = timestamp
    state.phase_at_anchor_rad = zero(state.phase_at_anchor_rad)
    state.reference_timestamp = timestamp
    state.reference_sequence = sequence
    state.reference_count += UInt64(1)
    return nothing
end

function initialize_autonomous_periodic_optic!(
    implementation, state, reference::AbstractWaveformPhaseReference)
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :unsupported_initialization,
        "prepared autonomous optic $(typeof(implementation)) does not initialize phase-reference state $(typeof(reference))"))
end

function reset_autonomous_periodic_optic_phase!(
    implementation, state, timestamp::PlantTimestamp, sequence::UInt64)
    throw(PlantScheduleError(:autonomous_periodic_optic,
        :unsupported_phase_reset,
        "prepared autonomous optic $(typeof(implementation)) does not reset phase at $timestamp for sequence $sequence"))
end

struct PreparedCycleAveragedPyramidModulation{M}
    modulation::M
end

@inline _autonomous_optic_couplings_conflict(::Any, ::Any) = false

@inline function _autonomous_optic_couplings_conflict(
    left::PreparedCycleAveragedPyramidModulation,
    right::PreparedCycleAveragedPyramidModulation)
    return left.modulation === right.modulation
end

@inline function _prepared_pyramid_modulation(
    plan::PreparedPyramidOpticalFormation)
    return plan.front_end.modulation
end

function _prepared_pyramid_modulation(
    plan::PreparedPyramidOpticalBundleFormation)
    modulation = _prepared_pyramid_modulation(first(plan.plans))
    @inbounds for component in Base.tail(plan.plans)
        _prepared_pyramid_modulation(component) === modulation || throw(
            PlantPreparationError(:autonomous_periodic_optic,
                :prepared_binding,
                "Pyramid bundle components do not share one prepared modulation"))
    end
    return modulation
end

@inline function _prepared_pyramid_modulation(
    execution::WFSOpticalPathExecution)
    return _prepared_pyramid_modulation(execution.plan)
end

function _prepared_pyramid_modulation(execution)
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :unsupported_path,
        "cycle-averaged Pyramid modulation requires a prepared Pyramid WFS optical path; got $(typeof(execution))"))
end

@inline function _require_circular_prepared_modulation(
    modulation::PreparedFocalPlaneModulation{<:CircularModulation})
    return modulation
end

function _require_circular_prepared_modulation(modulation)
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :unsupported_waveform,
        "native autonomous Pyramid modulation requires a prepared circular modulation; got $(typeof(modulation))"))
end

function prepare_autonomous_periodic_optic(
    implementation::PreparedCircularPyramidModulator{T},
    path::PreparedPathExecutor,
    ::CycleAveragedModulationFidelity) where {T}
    modulation = _require_circular_prepared_modulation(
        _prepared_pyramid_modulation(path.execution))
    eltype(modulation.amplitude_weights) === T || throw(
        PlantPreparationError(:autonomous_periodic_optic, :numeric_type,
            "Pyramid modulation precision differs from its setpoint endpoints"))
    return PreparedCycleAveragedPyramidModulation(modulation)
end

function prepare_autonomous_periodic_optic(
    implementation, path::PreparedPathExecutor,
    fidelity::AbstractAutonomousOpticFidelity)
    throw(PlantPreparationError(:autonomous_periodic_optic,
        :unsupported_binding,
        "prepared autonomous optic $(typeof(implementation)) does not bind path $(path_id(path.definition)) with fidelity $(typeof(fidelity))"))
end

function evaluate_autonomous_periodic_optic!(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState,
    coupling::PreparedCycleAveragedPyramidModulation,
    timestamp::PlantTimestamp)
    if state.optical_dirty
        update_cycle_averaged_circular_modulation!(
            coupling.modulation, state.radius; enabled=state.enabled)
        state.optical_dirty = false
    end
    state.last_evaluation_timestamp = timestamp
    state.has_evaluation = true
    return nothing
end

function evaluate_autonomous_periodic_optic!(
    implementation, state, coupling, timestamp::PlantTimestamp)
    throw(PlantScheduleError(:autonomous_periodic_optic,
        :unsupported_evaluation,
        "prepared autonomous optic $(typeof(implementation)) does not evaluate $(typeof(coupling)) at $timestamp"))
end

@inline function autonomous_waveform_phase(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState{T},
    timestamp::PlantTimestamp) where {T}
    return mod(_unwrapped_waveform_phase(state, timestamp), T(2pi))
end

@inline function autonomous_waveform_offset(
    implementation::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState{T},
    timestamp::PlantTimestamp) where {T}
    state.enabled || return (zero(T), zero(T))
    phase = autonomous_waveform_phase(implementation, state, timestamp)
    sine, cosine = sincos(phase)
    return state.radius * cosine, state.radius * sine
end

@inline autonomous_waveform_reference_timestamp(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState) = state.reference_timestamp
@inline autonomous_waveform_reference_sequence(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState) = state.reference_sequence
@inline autonomous_waveform_reference_count(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState) = state.reference_count
@inline autonomous_waveform_enabled(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState) = state.enabled
@inline autonomous_waveform_radius(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState) = state.radius

struct _PreparedAutonomousPeriodicOptic{I,C,R,F}
    id::ControllableOpticID
    path::OpticalPathID
    optic_slot::UInt32
    path_slot::UInt32
    implementation::I
    coupling::C
    phase_reference::R
    fidelity::F
end
