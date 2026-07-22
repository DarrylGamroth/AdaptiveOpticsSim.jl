const _PLANT_EVENT_LOOP_COMPONENT = :plant_event_loop

@noinline function _plant_event_loop_error(reason::Symbol,
    message::AbstractString)
    throw(PlantScheduleError(_PLANT_EVENT_LOOP_COMPONENT, reason,
        String(message)))
end

abstract type AbstractAcquisitionStartDefinition end

"""Periodic frame-start recurrence for one detector acquisition owner."""
struct PeriodicAcquisitionStart <: AbstractAcquisitionStartDefinition
    schedule::PeriodicSchedule
    origin::PlantTimestamp
end

PeriodicAcquisitionStart(schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp)) =
    PeriodicAcquisitionStart(schedule, origin)

"""Binding from one delivered trigger consumer to a detector frame start."""
struct TriggeredAcquisitionStart <: AbstractAcquisitionStartDefinition
    consumer::TriggerConsumerID
end

const _PreparedDetectorEventLifecycle = Union{
    PreparedGlobalShutterAcquisition,
    PreparedRollingShutterAcquisition,
    PreparedFrameTransferAcquisition,
}
const _DetectorEventLifecycleState = Union{
    GlobalShutterAcquisitionState,
    RollingShutterAcquisitionState,
    FrameTransferAcquisitionState,
}
const _AcquisitionStartDefinition = Union{
    PeriodicAcquisitionStart,
    TriggeredAcquisitionStart,
}

TriggeredAcquisitionStart(consumer::Symbol) =
    TriggeredAcquisitionStart(_as_trigger_consumer_id(consumer))

"""Periodic optical-sample recurrence for one prepared optical path."""
struct OpticalSampleDefinition
    path::OpticalPathID
    schedule::PeriodicSchedule
    origin::PlantTimestamp
end

OpticalSampleDefinition(path, schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp)) =
    OpticalSampleDefinition(_as_optical_path_id(path), schedule, origin)

"""
Detector lifecycle and start-source declaration for one prepared Gate 2
acquisition owner.
"""
struct DetectorEventDefinition{D<:AbstractDetectorAcquisitionEventDefinition,
    S<:AbstractAcquisitionStartDefinition}
    acquisition::AcquisitionID
    lifecycle::D
    start::S
end

DetectorEventDefinition(acquisition,
    lifecycle::AbstractDetectorAcquisitionEventDefinition,
    start::AbstractAcquisitionStartDefinition) = DetectorEventDefinition(
    _as_acquisition_id(acquisition), lifecycle, start)

@inline _require_optical_sample_definition(value::OpticalSampleDefinition) =
    value
@inline _require_detector_event_definition(value::DetectorEventDefinition) =
    value

function _require_optical_sample_definition(value)
    _plant_event_loop_error(:invalid_definition,
        "optical sample entries must be OpticalSampleDefinition values; got $(typeof(value))")
end

function _require_detector_event_definition(value)
    _plant_event_loop_error(:invalid_definition,
        "detector event entries must be DetectorEventDefinition values; got $(typeof(value))")
end

function _event_definition_tuple(values::Tuple, validator)
    foreach(validator, values)
    return values
end

function _event_definition_tuple(values::NamedTuple, validator)
    tuple_values = Tuple(values)
    foreach(validator, tuple_values)
    return tuple_values
end

function _event_definition_tuple(values::AbstractVector, validator)
    tuple_values = Tuple(values)
    foreach(validator, tuple_values)
    return tuple_values
end

function _event_definition_tuple(values, ::Any)
    _plant_event_loop_error(:invalid_registry,
        "plant event definitions must be finite Tuple, NamedTuple, or " *
        "AbstractVector values; got $(typeof(values))")
end

@inline _require_prepared_trigger_topology(::Nothing) = nothing
@inline _require_prepared_trigger_topology(
    topology::PreparedTriggerTopology) = topology

function _require_prepared_trigger_topology(topology)
    _plant_event_loop_error(:invalid_trigger_topology,
        "trigger_topology must be nothing or PreparedTriggerTopology; got $(typeof(topology))")
end

"""
    PlantEventLoopDefinition(optical_samples, detector_events;
        trigger_topology=nothing)

Cold, finite declaration of periodic path samples and independently periodic or
trigger-driven detector lifecycles. Preparation flattens these values into one
fixed-capacity serial scheduler; the declaration stores no mutable cursor,
detector state, or run-length event list.
"""
struct PlantEventLoopDefinition
    optical_samples::Tuple
    detector_events::Tuple
    trigger_topology::Union{Nothing,PreparedTriggerTopology}
end

function PlantEventLoopDefinition(optical_samples, detector_events;
    trigger_topology=nothing)
    samples = _event_definition_tuple(optical_samples,
        _require_optical_sample_definition)
    events = _event_definition_tuple(detector_events,
        _require_detector_event_definition)
    isempty(samples) && _plant_event_loop_error(:empty_paths,
        "plant event loop requires at least one optical sample definition")
    isempty(events) && _plant_event_loop_error(:empty_acquisitions,
        "plant event loop requires at least one detector event definition")
    topology = _require_prepared_trigger_topology(trigger_topology)
    return PlantEventLoopDefinition(samples, events, topology)
end

struct _NoPreparedTriggerTopology end
struct _NoTriggerTopologyState end
struct _NoTriggerTopologyWorkspace end

mutable struct _PlantEventLoopBinding end

@enum _PlantEventActionKind::UInt8 begin
    _TriggerTopologyAction = 0x01
    _AcquisitionBoundaryAction = 0x02
    _AcquisitionStartAction = 0x03
    _RollingBandOpenAction = 0x04
    _OpticalPathSampleAction = 0x05
    _AcquisitionReadoutAction = 0x06
    _AcquisitionReadinessAction = 0x07
end

struct _PlantEventAction
    kind::_PlantEventActionKind
    owner_slot::UInt32
end

struct _PreparedPlantEventPath
    id::OpticalPathID
    path::PreparedPathExecutor
    rngs::PreparedOwnerRNGs
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    handle::EventGeneratorHandle
end

struct _PreparedPlantEventAcquisition
    id::AcquisitionID
    lifecycle::_PreparedDetectorEventLifecycle
    observation::AbstractArray
    rng::Xoshiro
    start::_AcquisitionStartDefinition
    path_slot::UInt32
    start_handle::EventGeneratorHandle
    boundary_handle::EventGeneratorHandle
    band_open_handle::EventGeneratorHandle
    readout_handle::EventGeneratorHandle
    readiness_handle::EventGeneratorHandle
end

"""Run-immutable, fixed-capacity deterministic plant-event composition."""
struct PreparedPlantEventLoop{A<:AbstractTimedAtmosphere,R,T}
    binding::_PlantEventLoopBinding
    atmosphere::A
    atmosphere_rng::R
    scheduler::PreparedEventScheduler
    actions::Memory{_PlantEventAction}
    paths::Memory{_PreparedPlantEventPath}
    acquisitions::Memory{_PreparedPlantEventAcquisition}
    trigger_topology::T
end

@inline plant_event_path_count(prepared::PreparedPlantEventLoop) =
    length(prepared.paths)
@inline plant_event_acquisition_count(prepared::PreparedPlantEventLoop) =
    length(prepared.acquisitions)
@inline plant_event_generator_count(prepared::PreparedPlantEventLoop) =
    event_generator_count(prepared.scheduler)

@inline function _event_frame_execution(
    provider::PreparedFullOpticalProvider{<:FrameAcquisitionExecution})
    return provider.execution
end

function _event_frame_execution(provider)
    _plant_event_loop_error(:unsupported_acquisition,
        "plant event composition currently requires a full-optical " *
        "FrameAcquisitionExecution provider; got $(typeof(provider))")
end

@inline function _event_frame_execution(owner::PreparedAcquisitionOwner)
    return _event_frame_execution(owner.provider.implementation)
end

@inline function _prepare_detector_event_lifecycle(
    execution::FrameAcquisitionExecution, result::IntensityMap,
    definition::GlobalShutterAcquisitionDefinition)
    return prepare_global_shutter_acquisition(execution.detector, result,
        execution.plan, definition)
end

@inline function _prepare_detector_event_lifecycle(
    execution::FrameAcquisitionExecution, result::IntensityMap,
    definition::RollingShutterAcquisitionDefinition)
    return prepare_rolling_shutter_acquisition(execution.detector, result,
        execution.plan, definition)
end

@inline function _prepare_detector_event_lifecycle(
    execution::FrameAcquisitionExecution, result::IntensityMap,
    definition::FrameTransferAcquisitionDefinition)
    return prepare_frame_transfer_acquisition(execution.detector, result,
        execution.plan, definition)
end

function _sorted_optical_sample_definitions(
    definitions::Tuple)
    values = Any[definitions...]
    sort!(values; by=definition -> definition.path.name)
    @inbounds for index in 2:length(values)
        values[index - 1].path == values[index].path &&
            _plant_event_loop_error(:duplicate_path,
                "optical path $(values[index].path) has more than one sample schedule")
    end
    return values
end

function _sorted_detector_event_definitions(definitions::Tuple)
    values = Any[definitions...]
    sort!(values; by=definition -> definition.acquisition.name)
    @inbounds for index in 2:length(values)
        values[index - 1].acquisition == values[index].acquisition &&
            _plant_event_loop_error(:duplicate_acquisition,
                "acquisition $(values[index].acquisition) has more than one detector lifecycle")
    end
    return values
end

function _prepared_event_path_rngs(plant::PreparedPlant,
    path::PreparedPathExecutor)
    Base.@nospecialize plant path
    paths = getfield(plant, :paths)
    rngs = getfield(getfield(plant, :rngs), :paths)
    @inbounds for index in eachindex(paths)
        paths[index] === path && return rngs[index]
    end
    _plant_event_loop_error(:prepared_binding,
        "event path has no exact prepared RNG owner")
end

function _prepared_event_acquisition_rngs(plant::PreparedPlant,
    owner::PreparedAcquisitionOwner)
    Base.@nospecialize plant owner
    acquisitions = getfield(plant, :acquisitions)
    rngs = getfield(getfield(plant, :rngs), :acquisitions)
    @inbounds for index in eachindex(acquisitions)
        acquisitions[index] === owner && return rngs[index]
    end
    _plant_event_loop_error(:prepared_binding,
        "event acquisition has no exact prepared RNG owner")
end

function _event_prepared_path(plant::PreparedPlant, id::OpticalPathID)
    Base.@nospecialize plant
    for path in getfield(plant, :paths)
        path_id(path.definition) == id && return path
    end
    _plant_event_loop_error(:unknown_path,
        "prepared plant has no optical path $id")
end

function _event_prepared_acquisition(plant::PreparedPlant,
    id::AcquisitionID)
    Base.@nospecialize plant
    for owner in getfield(plant, :acquisitions)
        acquisition_id(owner.definition) == id && return owner
    end
    _plant_event_loop_error(:unknown_acquisition,
        "prepared plant has no acquisition $id")
end

function _event_path_slot(paths, id::OpticalPathID)
    @inbounds for index in eachindex(paths)
        paths[index].id == id && return index
    end
    _plant_event_loop_error(:missing_path_schedule,
        "detector acquisition references path $id without an optical sample schedule")
end

function _event_path_slot_from_definitions(definitions, id::OpticalPathID)
    @inbounds for index in eachindex(definitions)
        definitions[index].path == id && return index
    end
    _plant_event_loop_error(:missing_path_schedule,
        "detector acquisition references path $id without an optical sample schedule")
end

function _trigger_consumer_exists(topology::PreparedTriggerTopology,
    id::TriggerConsumerID)
    @inbounds for consumer in topology.consumers
        consumer.id == id && return true
    end
    return false
end

@inline function _require_start_trigger_topology(
    ::PeriodicAcquisitionStart, ::Nothing)
    return nothing
end

@inline function _require_start_trigger_topology(
    ::PeriodicAcquisitionStart, ::PreparedTriggerTopology)
    return nothing
end

function _require_start_trigger_topology(start::TriggeredAcquisitionStart,
    ::Nothing)
    _plant_event_loop_error(:missing_trigger_topology,
        "triggered acquisition $(start.consumer) requires a prepared trigger topology")
end

function _require_start_trigger_topology(start::TriggeredAcquisitionStart,
    topology::PreparedTriggerTopology)
    _trigger_consumer_exists(topology, start.consumer) ||
        _plant_event_loop_error(:unknown_trigger_consumer,
            "triggered acquisition references unknown consumer $(start.consumer)")
    return nothing
end

function _require_unique_trigger_consumers(definitions)
    @inbounds for right in 2:length(definitions)
        right_start = definitions[right].start
        _require_unique_trigger_consumer(right_start, definitions, right)
    end
    return nothing
end

@inline _require_bound_trigger_consumers(::Any, ::Nothing) = nothing

function _require_bound_trigger_consumers(definitions,
    topology::PreparedTriggerTopology)
    @inbounds for consumer in topology.consumers
        bound = false
        for definition in definitions
            if _start_matches_consumer(definition.start, consumer.id)
                bound = true
                break
            end
        end
        bound || _plant_event_loop_error(:unbound_trigger_consumer,
            "trigger consumer $(consumer.id) has no detector acquisition binding")
    end
    return nothing
end

@inline _require_unique_trigger_consumer(
    ::PeriodicAcquisitionStart, ::Any, ::Int) = nothing

function _require_unique_trigger_consumer(start::TriggeredAcquisitionStart,
    definitions, right::Int)
    @inbounds for left in 1:(right - 1)
        _same_trigger_consumer(definitions[left].start, start) &&
            _plant_event_loop_error(:duplicate_trigger_consumer,
                "trigger consumer $(start.consumer) is bound to more than one acquisition")
    end
    return nothing
end

@inline _same_trigger_consumer(::PeriodicAcquisitionStart,
    ::TriggeredAcquisitionStart) = false
@inline _same_trigger_consumer(left::TriggeredAcquisitionStart,
    right::TriggeredAcquisitionStart) = left.consumer == right.consumer

@inline function _periodic_start_timestamp(start::PeriodicAcquisitionStart)
    return schedule_timestamp(start.schedule, 1, start.origin)
end

@inline _start_generator_definition(start::PeriodicAcquisitionStart,
    ordinal::Integer) = EventGeneratorDefinition(start.schedule,
    ExposureOpenPhase, ordinal; origin=start.origin)

@inline _start_generator_definition(::TriggeredAcquisitionStart,
    ordinal::Integer) = EventGeneratorDefinition(zero(PlantTimestamp),
    ExposureOpenPhase, ordinal; active=false)

@inline function _require_periodic_start_spacing(
    ::AbstractPreparedDetectorAcquisition,
    ::TriggeredAcquisitionStart)
    return nothing
end

function _require_periodic_start_spacing(
    prepared::PreparedGlobalShutterAcquisition,
    start::PeriodicAcquisitionStart)
    occupied = prepared.definition.exposure_duration +
        prepared.definition.readout_duration +
        prepared.definition.readiness_delay
    schedule_period(start.schedule) > occupied ||
        _plant_event_loop_error(:acquisition_period,
            "global-shutter period must be strictly later than acquisition readiness")
    return nothing
end

function _require_periodic_start_spacing(
    prepared::PreparedRollingShutterAcquisition,
    start::PeriodicAcquisitionStart)
    occupied = _rolling_frame_readout_offset(prepared) +
        prepared.definition.readiness_delay
    schedule_period(start.schedule) > occupied ||
        _plant_event_loop_error(:acquisition_period,
            "rolling-shutter period must be strictly later than complete-frame readiness")
    return nothing
end

function _require_periodic_start_spacing(
    prepared::PreparedFrameTransferAcquisition,
    start::PeriodicAcquisitionStart)
    period = schedule_period(start.schedule)
    image_reuse = prepared.definition.exposure_duration +
        prepared.transfer_duration
    period >= image_reuse || _plant_event_loop_error(:acquisition_period,
        "frame-transfer period is shorter than image-area exposure plus transfer")
    period > prepared.definition.readout_duration ||
        _plant_event_loop_error(:acquisition_period,
            "frame-transfer period must place the next storage transfer after prior readout")
    return nothing
end

function _initial_trigger_timestamp(topology::PreparedTriggerTopology)
    state = TriggerTopologyState(topology)
    return realized_trigger_source_timestamp(next_trigger_source(topology,
        state))
end

function _append_event_generator_definitions!(definitions,
    sample_definitions, acquisition_definitions, trigger_topology)
    if trigger_topology !== nothing
        push!(definitions, EventGeneratorDefinition(
            _initial_trigger_timestamp(trigger_topology),
            TriggerUpdatePhase, 1))
    end
    for (index, _) in enumerate(acquisition_definitions)
        push!(definitions, EventGeneratorDefinition(zero(PlantTimestamp),
            IntegrationBoundaryPhase, index; active=false))
    end
    for (index, definition) in enumerate(acquisition_definitions)
        push!(definitions, _start_generator_definition(definition.start,
            2index - 1))
        push!(definitions, EventGeneratorDefinition(zero(PlantTimestamp),
            ExposureOpenPhase, 2index; active=false))
    end
    for (index, definition) in enumerate(sample_definitions)
        push!(definitions, EventGeneratorDefinition(definition.schedule,
            OpticalSamplePhase, index; origin=definition.origin))
    end
    for (index, _) in enumerate(acquisition_definitions)
        push!(definitions, EventGeneratorDefinition(zero(PlantTimestamp),
            ReadoutCompletionPhase, index; active=false))
        push!(definitions, EventGeneratorDefinition(zero(PlantTimestamp),
            AcquisitionReadyPhase, index; active=false))
    end
    return definitions
end

@inline _event_action_for_definition(
    ::Val{TriggerUpdatePhase}, ordinal::UInt32) =
    _PlantEventAction(_TriggerTopologyAction, ordinal)
@inline _event_action_for_definition(
    ::Val{IntegrationBoundaryPhase}, ordinal::UInt32) =
    _PlantEventAction(_AcquisitionBoundaryAction, ordinal)
@inline _event_action_for_definition(
    ::Val{OpticalSamplePhase}, ordinal::UInt32) =
    _PlantEventAction(_OpticalPathSampleAction, ordinal)
@inline _event_action_for_definition(
    ::Val{ReadoutCompletionPhase}, ordinal::UInt32) =
    _PlantEventAction(_AcquisitionReadoutAction, ordinal)
@inline _event_action_for_definition(
    ::Val{AcquisitionReadyPhase}, ordinal::UInt32) =
    _PlantEventAction(_AcquisitionReadinessAction, ordinal)

@inline function _event_action_for_definition(
    ::Val{ExposureOpenPhase}, ordinal::UInt32)
    owner = (ordinal + UInt32(1)) >> 1
    kind = isodd(ordinal) ? _AcquisitionStartAction :
        _RollingBandOpenAction
    return _PlantEventAction(kind, owner)
end

@inline function _event_action_for_definition(
    definition::EventGeneratorDefinition)
    return _event_action_for_definition(Val(definition.phase),
        definition.ordinal)
end

function _prepared_event_actions(scheduler::PreparedEventScheduler)
    actions = Memory{_PlantEventAction}(undef,
        event_generator_count(scheduler))
    @inbounds for index in eachindex(actions)
        actions[index] = _event_action_for_definition(
            scheduler.definitions[index])
    end
    return actions
end

function _prepare_event_paths(plant::PreparedPlant, definitions,
    scheduler::PreparedEventScheduler)
    Base.@nospecialize plant
    paths = Memory{_PreparedPlantEventPath}(undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        path = _event_prepared_path(plant, definition.path)
        rngs = _prepared_event_path_rngs(plant, path)
        _require_rng_owner_binding(rngs, path)
        handle = event_generator_handle(scheduler, OpticalSamplePhase, index)
        paths[index] = _PreparedPlantEventPath(definition.path, path, rngs,
            definition.schedule, definition.origin, handle)
    end
    return paths
end

function _prepare_event_acquisition_parts(plant::PreparedPlant,
    definitions, path_definitions, topology)
    Base.@nospecialize plant
    owners = Any[]
    lifecycles = Any[]
    rngs = Any[]
    path_slots = Int[]
    for definition in definitions
        _require_start_trigger_topology(definition.start, topology)
        owner = _event_prepared_acquisition(plant, definition.acquisition)
        execution = _event_frame_execution(owner)
        lifecycle = _prepare_detector_event_lifecycle(execution,
            owner.path_result, definition.lifecycle)
        _require_periodic_start_spacing(lifecycle, definition.start)
        acquisition_rngs = _prepared_event_acquisition_rngs(plant, owner)
        _require_rng_owner_binding(acquisition_rngs, owner)
        push!(owners, owner)
        push!(lifecycles, lifecycle)
        push!(rngs, rng_stream_state(acquisition_rngs, Val(:detector)))
        push!(path_slots, _event_path_slot_from_definitions(path_definitions,
            acquisition_path_id(owner.definition)))
    end
    return owners, lifecycles, rngs, path_slots
end

function _prepare_event_acquisitions(definitions, owners, lifecycles, rngs,
    path_slots, scheduler::PreparedEventScheduler)
    acquisitions = Memory{_PreparedPlantEventAcquisition}(undef,
        length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        acquisitions[index] = _PreparedPlantEventAcquisition(
            definition.acquisition, lifecycles[index],
            acquisition_observation(owners[index]), rngs[index],
            definition.start, UInt32(path_slots[index]),
            event_generator_handle(scheduler, ExposureOpenPhase,
                2index - 1),
            event_generator_handle(scheduler, IntegrationBoundaryPhase,
                index),
            event_generator_handle(scheduler, ExposureOpenPhase, 2index),
            event_generator_handle(scheduler, ReadoutCompletionPhase,
                index),
            event_generator_handle(scheduler, AcquisitionReadyPhase,
                index))
    end
    return acquisitions
end

function _prepared_trigger_topology(topology::Nothing)
    return _NoPreparedTriggerTopology()
end

@inline _prepared_trigger_topology(topology::PreparedTriggerTopology) =
    topology

"""
    prepare_plant_event_loop(plant, definition)

Bind a Gate 2 prepared plant to a flat deterministic scheduler, exact
owner-derived RNG streams, bounded detector lifecycle state, and an optional
prepared trigger topology. Preparation allocates; repeated stepping does not
grow the registry or materialize future events.
"""
function prepare_plant_event_loop(plant::PreparedPlant,
    definition::PlantEventLoopDefinition)
    Base.@nospecialize plant
    sample_definitions = _sorted_optical_sample_definitions(
        definition.optical_samples)
    acquisition_definitions = _sorted_detector_event_definitions(
        definition.detector_events)
    _require_unique_trigger_consumers(acquisition_definitions)
    _require_bound_trigger_consumers(acquisition_definitions,
        definition.trigger_topology)
    owners, lifecycles, rngs, path_slots =
        _prepare_event_acquisition_parts(plant, acquisition_definitions,
            sample_definitions, definition.trigger_topology)
    generator_definitions = EventGeneratorDefinition[]
    _append_event_generator_definitions!(generator_definitions,
        sample_definitions, acquisition_definitions,
        definition.trigger_topology)
    scheduler = prepare_event_scheduler(generator_definitions;
        capacity=length(generator_definitions))
    actions = _prepared_event_actions(scheduler)
    paths = _prepare_event_paths(plant, sample_definitions, scheduler)
    acquisitions = _prepare_event_acquisitions(acquisition_definitions,
        owners, lifecycles, rngs, path_slots, scheduler)
    atmosphere = _require_selection_atmosphere(
        plant_atmosphere(getfield(plant, :definition)))
    atmosphere_rng = _prepared_atmosphere_rng(atmosphere,
        getfield(getfield(plant, :rngs), :atmosphere))
    return PreparedPlantEventLoop(_PlantEventLoopBinding(), atmosphere,
        atmosphere_rng, scheduler, actions, paths, acquisitions,
        _prepared_trigger_topology(definition.trigger_topology))
end

mutable struct PlantEventLoopState{T}
    binding::_PlantEventLoopBinding
    scheduler::EventSchedulerState
    acquisitions::Memory{_DetectorEventLifecycleState}
    path_sampled::Memory{Bool}
    product_sequences::Memory{UInt64}
    product_ready_timestamps::Memory{PlantTimestamp}
    trigger::T
end

@inline _event_acquisition_state(
    prepared::PreparedGlobalShutterAcquisition) =
    GlobalShutterAcquisitionState(prepared)
@inline _event_acquisition_state(
    prepared::PreparedRollingShutterAcquisition) =
    RollingShutterAcquisitionState(prepared)
@inline _event_acquisition_state(
    prepared::PreparedFrameTransferAcquisition) =
    FrameTransferAcquisitionState(prepared)

@inline _event_trigger_state(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyState()
@inline _event_trigger_state(topology::PreparedTriggerTopology) =
    TriggerTopologyState(topology)

function PlantEventLoopState(prepared::PreparedPlantEventLoop)
    acquisition_states = Memory{_DetectorEventLifecycleState}(undef,
        length(prepared.acquisitions))
    @inbounds for index in eachindex(acquisition_states)
        acquisition_states[index] = _event_acquisition_state(
            prepared.acquisitions[index].lifecycle)
    end
    path_sampled = Memory{Bool}(undef, length(prepared.paths))
    fill!(path_sampled, false)
    product_sequences = Memory{UInt64}(undef,
        length(prepared.acquisitions))
    fill!(product_sequences, UInt64(0))
    product_ready_timestamps = Memory{PlantTimestamp}(undef,
        length(prepared.acquisitions))
    fill!(product_ready_timestamps, zero(PlantTimestamp))
    trigger = _event_trigger_state(prepared.trigger_topology)
    return PlantEventLoopState(prepared.binding,
        EventSchedulerState(prepared.scheduler), acquisition_states,
        path_sampled, product_sequences, product_ready_timestamps, trigger)
end

mutable struct PlantEventLoopWorkspace{T}
    binding::_PlantEventLoopBinding
    scheduler::EventSchedulerWorkspace
    due_paths::Memory{Bool}
    trigger::T
    delivery::Base.RefValue{TriggerDelivery}
end

@inline _event_trigger_workspace(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyWorkspace()
@inline _event_trigger_workspace(topology::PreparedTriggerTopology) =
    TriggerTopologyWorkspace(topology)

function PlantEventLoopWorkspace(prepared::PreparedPlantEventLoop)
    due_paths = Memory{Bool}(undef, length(prepared.paths))
    fill!(due_paths, false)
    return PlantEventLoopWorkspace(prepared.binding,
        EventSchedulerWorkspace(prepared.scheduler), due_paths,
        _event_trigger_workspace(prepared.trigger_topology),
        Ref{TriggerDelivery}())
end

@inline function _require_plant_event_loop_binding(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState)
    state.binding === prepared.binding || _plant_event_loop_error(
        :foreign_state,
        "plant event-loop state belongs to another prepared loop")
    _require_scheduler_binding(prepared.scheduler, state.scheduler)
    length(state.acquisitions) == length(prepared.acquisitions) &&
        length(state.path_sampled) == length(prepared.paths) &&
        length(state.product_sequences) == length(prepared.acquisitions) &&
        length(state.product_ready_timestamps) ==
            length(prepared.acquisitions) ||
        _plant_event_loop_error(:prepared_binding,
            "plant event-loop state capacity changed after preparation")
    return nothing
end

@inline function _require_plant_event_loop_binding(
    prepared::PreparedPlantEventLoop, workspace::PlantEventLoopWorkspace)
    workspace.binding === prepared.binding || _plant_event_loop_error(
        :foreign_workspace,
        "plant event-loop workspace belongs to another prepared loop")
    _require_scheduler_binding(prepared.scheduler, workspace.scheduler)
    length(workspace.due_paths) == length(prepared.paths) ||
        _plant_event_loop_error(:prepared_binding,
            "plant event-loop workspace capacity changed after preparation")
    return nothing
end

function _event_acquisition_slot(prepared::PreparedPlantEventLoop,
    id::AcquisitionID)
    @inbounds for index in eachindex(prepared.acquisitions)
        prepared.acquisitions[index].id == id && return index
    end
    _plant_event_loop_error(:unknown_acquisition,
        "prepared plant event loop has no acquisition $id")
end

@inline _event_acquisition_slot(prepared::PreparedPlantEventLoop,
    name::Symbol) = _event_acquisition_slot(prepared, AcquisitionID(name))

function acquisition_product_sequence(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    return @inbounds state.product_sequences[
        _event_acquisition_slot(prepared, _as_acquisition_id(id))]
end

function acquisition_product_ready_timestamp(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    slot = _event_acquisition_slot(prepared, _as_acquisition_id(id))
    @inbounds iszero(state.product_sequences[slot]) && return nothing
    return @inbounds state.product_ready_timestamps[slot]
end

@inline function _event_action(prepared::PreparedPlantEventLoop,
    claim::EventClaim)
    slot = Int(claim.slot)
    1 <= slot <= length(prepared.actions) ||
        _plant_event_loop_error(:invalid_action,
            "event claim does not map to a prepared plant action")
    return @inbounds prepared.actions[slot]
end

@inline function _event_acquisition_binding(
    prepared::PreparedPlantEventLoop, slot::UInt32)
    index = Int(slot)
    1 <= index <= length(prepared.acquisitions) ||
        _plant_event_loop_error(:invalid_action,
            "event action contains an invalid acquisition slot")
    return @inbounds prepared.acquisitions[index]
end

@inline function _event_acquisition_state(state::PlantEventLoopState,
    slot::UInt32)
    index = Int(slot)
    1 <= index <= length(state.acquisitions) ||
        _plant_event_loop_error(:invalid_action,
            "event action contains an invalid acquisition-state slot")
    return @inbounds state.acquisitions[index]
end

@inline function _event_path_binding(prepared::PreparedPlantEventLoop,
    slot::UInt32)
    index = Int(slot)
    1 <= index <= length(prepared.paths) ||
        _plant_event_loop_error(:invalid_action,
            "event action contains an invalid optical-path slot")
    return @inbounds prepared.paths[index]
end

@inline function _require_inactive_event_generator(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    handle::EventGeneratorHandle, timestamp::PlantTimestamp)
    slot = _require_event_generator_slot(prepared.scheduler, handle)
    cursor = state.scheduler.cursors[slot]
    cursor.status == _InactiveEventGenerator ||
        _plant_event_loop_error(:generator_busy,
            "required prepared event generator is already active")
    definition = prepared.scheduler.definitions[slot]
    key = PlantEventKey(timestamp, definition.phase, definition.ordinal,
        cursor.next_occurrence, _PLANT_TIME_TOKEN)
    _require_forward_event_key(state.scheduler, key)
    return key
end

@inline function _event_generator_due_at(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    handle::EventGeneratorHandle, timestamp::PlantTimestamp)
    slot = _require_event_generator_slot(prepared.scheduler, handle)
    cursor = state.scheduler.cursors[slot]
    return cursor.status == _ScheduledEventGenerator &&
        cursor.next_timestamp == timestamp
end

@inline function _next_periodic_start_timestamp(
    start::PeriodicAcquisitionStart, claim::EventClaim)
    occurrence = _next_event_occurrence(claim.key.occurrence)
    return schedule_timestamp(start.schedule, occurrence, start.origin)
end

@inline _next_periodic_start_timestamp(
    ::TriggeredAcquisitionStart, ::EventClaim) = nothing

@inline function _resolve_start_claim!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    start::PeriodicAcquisitionStart, next_timestamp::PlantTimestamp)
    return reschedule_event!(prepared.scheduler, state.scheduler, claim,
        next_timestamp)
end

@inline function _resolve_start_claim!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    ::TriggeredAcquisitionStart, ::Nothing)
    return deactivate_event_generator!(prepared.scheduler, state.scheduler,
        claim)
end

@inline function _first_detector_boundary_timestamp(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    next_read = next_nondestructive_read_timestamp(prepared, state)
    next_read === nothing && return state.exposure_close
    return min(next_read, state.exposure_close)
end

@inline function _initial_detector_boundary_timestamp(
    prepared::PreparedGlobalShutterAcquisition,
    timestamp::PlantTimestamp)
    isempty(prepared.read_offsets) &&
        return timestamp + prepared.definition.exposure_duration
    first_index = iszero(@inbounds(prepared.read_offsets[1])) ? 2 : 1
    first_index > length(prepared.read_offsets) &&
        return timestamp + prepared.definition.exposure_duration
    return timestamp + min(@inbounds(prepared.read_offsets[first_index]),
        prepared.definition.exposure_duration)
end

@inline function _initial_detector_boundary_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    timestamp::PlantTimestamp)
    return timestamp + prepared.definition.exposure_duration
end

@inline function _initial_detector_boundary_timestamp(
    prepared::PreparedFrameTransferAcquisition,
    timestamp::PlantTimestamp)
    return timestamp + prepared.definition.exposure_duration
end

@inline function _first_detector_boundary_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    return next_rolling_band_close_timestamp(prepared, state)
end

@inline function _first_detector_boundary_timestamp(
    ::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState)
    return state.exposure_close
end

@inline _first_detector_band_open_timestamp(
    ::PreparedGlobalShutterAcquisition,
    ::GlobalShutterAcquisitionState) = nothing
@inline _first_detector_band_open_timestamp(
    ::PreparedFrameTransferAcquisition,
    ::FrameTransferAcquisitionState) = nothing
@inline function _first_detector_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    return next_rolling_band_open_timestamp(prepared, state)
end

@inline _initial_detector_band_open_timestamp(
    ::PreparedGlobalShutterAcquisition, ::PlantTimestamp) = nothing
@inline _initial_detector_band_open_timestamp(
    ::PreparedFrameTransferAcquisition, ::PlantTimestamp) = nothing
@inline function _initial_detector_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_RollingExposureEventMode},
    timestamp::PlantTimestamp)
    prepared.band_count == 1 && return nothing
    return timestamp + prepared.line_duration
end
@inline _initial_detector_band_open_timestamp(
    ::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_GlobalResetEventMode},
    ::PlantTimestamp) = nothing

@inline function _take_initial_detector_snapshot!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    next_read = next_nondestructive_read_timestamp(prepared, state)
    next_read == timestamp || return nothing
    take_nondestructive_read!(prepared, state, timestamp, rng)
    return nothing
end

@inline _take_initial_detector_snapshot!(
    ::PreparedRollingShutterAcquisition,
    ::RollingShutterAcquisitionState, ::PlantTimestamp,
    ::AbstractRNG) = nothing
@inline _take_initial_detector_snapshot!(
    ::PreparedFrameTransferAcquisition,
    ::FrameTransferAcquisitionState, ::PlantTimestamp,
    ::AbstractRNG) = nothing

@inline function _require_event_path_available(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    acquisition::_PreparedPlantEventAcquisition,
    timestamp::PlantTimestamp)
    path_slot = Int(acquisition.path_slot)
    @inbounds state.path_sampled[path_slot] && return nothing
    path = @inbounds prepared.paths[path_slot]
    _event_generator_due_at(prepared, state, path.handle, timestamp) ||
        _plant_event_loop_error(:uninitialized_path,
            "acquisition $(acquisition.id) begins before its first optical sample")
    return nothing
end

function _process_acquisition_start!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    timestamp = claim.key.timestamp
    _require_event_path_available(prepared, state, acquisition, timestamp)
    next_start = _next_periodic_start_timestamp(acquisition.start, claim)
    next_start === nothing || begin
        definition = prepared.scheduler.definitions[Int(claim.slot)]
        next_key = PlantEventKey(next_start, definition.phase,
            definition.ordinal, _next_event_occurrence(claim.key.occurrence),
            _PLANT_TIME_TOKEN)
        _require_forward_event_key(state.scheduler, next_key)
    end

    boundary_timestamp = _initial_detector_boundary_timestamp(
        acquisition.lifecycle, timestamp)
    _require_inactive_event_generator(prepared, state,
        acquisition.boundary_handle, boundary_timestamp)
    band_open_timestamp = _initial_detector_band_open_timestamp(
        acquisition.lifecycle, timestamp)
    band_open_timestamp === nothing || _require_inactive_event_generator(
        prepared, state, acquisition.band_open_handle, band_open_timestamp)

    begin_exposure!(acquisition.lifecycle, acquisition_state, timestamp)
    _take_initial_detector_snapshot!(acquisition.lifecycle,
        acquisition_state, timestamp, acquisition.rng)
    boundary_timestamp == _first_detector_boundary_timestamp(
        acquisition.lifecycle, acquisition_state) ||
        _plant_event_loop_error(:prepared_binding,
            "detector boundary changed after acquisition start")
    band_open_timestamp == _first_detector_band_open_timestamp(
        acquisition.lifecycle, acquisition_state) ||
        _plant_event_loop_error(:prepared_binding,
            "rolling row-band schedule changed after acquisition start")

    _resolve_start_claim!(prepared, state, claim, acquisition.start,
        next_start)
    activate_event_generator!(prepared.scheduler, state.scheduler,
        acquisition.boundary_handle, boundary_timestamp)
    band_open_timestamp === nothing || activate_event_generator!(
        prepared.scheduler, state.scheduler, acquisition.band_open_handle,
        band_open_timestamp)
    return nothing
end

@inline function _integrate_event_acquisition_to!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    state.status == DetectorExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "global-shutter integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    return accumulate_exposure_interval!(prepared, state,
        state.integrated_through, timestamp, rng)
end

@inline function _integrate_event_acquisition_to!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    state.status == DetectorExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "rolling-shutter integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    if state.opened_bands == state.closed_bands
        state.integrated_through = timestamp
        return nothing
    end
    return accumulate_rolling_exposure_interval!(prepared, state,
        state.integrated_through, timestamp, rng)
end

@inline function _integrate_event_acquisition_to!(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    state.image_status == _FrameTransferImageActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "frame-transfer integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    return accumulate_exposure_interval!(prepared, state,
        state.integrated_through, timestamp, rng)
end

@enum _DetectorBoundaryDisposition::UInt8 begin
    _RescheduleDetectorBoundary = 0x01
    _ScheduleDetectorReadout = 0x02
end

struct _DetectorBoundaryResult
    disposition::_DetectorBoundaryDisposition
    timestamp::PlantTimestamp
end

function _process_detector_boundary!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
    next_read = next_nondestructive_read_timestamp(prepared, state)
    if next_read !== nothing && next_read == timestamp
        take_nondestructive_read!(prepared, state, timestamp, rng)
        following_read = next_nondestructive_read_timestamp(prepared, state)
        following = following_read === nothing ? state.exposure_close :
            min(following_read, state.exposure_close)
        return _DetectorBoundaryResult(_RescheduleDetectorBoundary,
            following)
    end
    close_exposure!(prepared, state, timestamp)
    return _DetectorBoundaryResult(_ScheduleDetectorReadout,
        state.readout_complete)
end

function _process_detector_boundary!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
    close_next_rolling_band!(prepared, state, timestamp)
    following = next_rolling_band_close_timestamp(prepared, state)
    following === nothing && return _DetectorBoundaryResult(
        _ScheduleDetectorReadout, state.readout_complete)
    return _DetectorBoundaryResult(_RescheduleDetectorBoundary, following)
end

function _process_detector_boundary!(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    if state.image_status == _FrameTransferImageActive
        _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
        close_exposure!(prepared, state, timestamp)
        return _DetectorBoundaryResult(_RescheduleDetectorBoundary,
            state.transfer_complete)
    end
    complete_frame_transfer!(prepared, state, timestamp)
    return _DetectorBoundaryResult(_ScheduleDetectorReadout,
        state.storage_readout_complete)
end

function _process_acquisition_boundary!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    result = _process_detector_boundary!(acquisition.lifecycle,
        acquisition_state, claim.key.timestamp, acquisition.rng)
    if result.disposition == _RescheduleDetectorBoundary
        reschedule_event!(prepared.scheduler, state.scheduler, claim,
            result.timestamp)
    else
        _require_inactive_event_generator(prepared, state,
            acquisition.readout_handle, result.timestamp)
        deactivate_event_generator!(prepared.scheduler, state.scheduler,
            claim)
        activate_event_generator!(prepared.scheduler, state.scheduler,
            acquisition.readout_handle, result.timestamp)
    end
    return nothing
end

function _process_rolling_band_open!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    _integrate_event_acquisition_to!(acquisition.lifecycle,
        acquisition_state, claim.key.timestamp, acquisition.rng)
    open_next_rolling_band!(acquisition.lifecycle, acquisition_state,
        claim.key.timestamp)
    following = next_rolling_band_open_timestamp(acquisition.lifecycle,
        acquisition_state)
    if following === nothing
        deactivate_event_generator!(prepared.scheduler, state.scheduler,
            claim)
    else
        reschedule_event!(prepared.scheduler, state.scheduler, claim,
            following)
    end
    return nothing
end

@inline _event_requires_readiness(
    ::PreparedGlobalShutterAcquisition) = true
@inline _event_requires_readiness(
    ::PreparedRollingShutterAcquisition) = true
@inline _event_requires_readiness(
    ::PreparedFrameTransferAcquisition) = false

@inline _event_readiness_timestamp(
    ::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState) = state.readiness
@inline _event_readiness_timestamp(
    ::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState) = state.readiness

@inline _event_product_sequence(
    state::GlobalShutterAcquisitionState) = state.sequence
@inline _event_product_sequence(
    state::RollingShutterAcquisitionState) = state.sequence
@inline _event_product_sequence(
    state::FrameTransferAcquisitionState) = state.product_sequence

function _publish_acquisition_product!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, slot::UInt32, output)
    acquisition = _event_acquisition_binding(prepared, slot)
    observation = acquisition.observation
    copyto!(observation, output)
    index = Int(slot)
    sequence = _event_product_sequence(
        _event_acquisition_state(state, slot))
    previous_sequence = @inbounds state.product_sequences[index]
    sequence > previous_sequence ||
        _plant_event_loop_error(:product_sequence,
            "detector product sequence did not advance")
    @inbounds state.product_sequences[index] = sequence
    @inbounds state.product_ready_timestamps[index] =
        state.scheduler.current_timestamp
    return observation
end

function _process_acquisition_readout!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    output = complete_readout!(acquisition.lifecycle, acquisition_state,
        claim.key.timestamp, acquisition.rng)
    _publish_acquisition_product!(prepared, state, action.owner_slot, output)
    if _event_requires_readiness(acquisition.lifecycle)
        ready_timestamp = _event_readiness_timestamp(acquisition.lifecycle,
            acquisition_state)
        _require_inactive_event_generator(prepared, state,
            acquisition.readiness_handle, ready_timestamp)
        deactivate_event_generator!(prepared.scheduler, state.scheduler,
            claim)
        activate_event_generator!(prepared.scheduler, state.scheduler,
            acquisition.readiness_handle, ready_timestamp)
    else
        deactivate_event_generator!(prepared.scheduler, state.scheduler,
            claim)
    end
    return nothing
end

function _process_acquisition_readiness!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    mark_acquisition_ready!(acquisition.lifecycle, acquisition_state,
        claim.key.timestamp)
    deactivate_event_generator!(prepared.scheduler, state.scheduler, claim)
    return nothing
end

function _triggered_acquisition_slot(prepared::PreparedPlantEventLoop,
    consumer::TriggerConsumerID)
    @inbounds for index in eachindex(prepared.acquisitions)
        start = prepared.acquisitions[index].start
        _start_matches_consumer(start, consumer) && return UInt32(index)
    end
    _plant_event_loop_error(:unknown_trigger_consumer,
        "delivered trigger consumer $consumer is not bound to an acquisition")
end

@inline _start_matches_consumer(::PeriodicAcquisitionStart,
    ::TriggerConsumerID) = false
@inline _start_matches_consumer(start::TriggeredAcquisitionStart,
    consumer::TriggerConsumerID) = start.consumer == consumer

@inline function _next_trigger_action_timestamp(
    topology::PreparedTriggerTopology, state::TriggerTopologyState)
    source = next_trigger_source(topology, state)
    delivery_timestamp = next_trigger_delivery_timestamp(topology, state)
    delivery_timestamp === nothing &&
        return realized_trigger_source_timestamp(source)
    return min(realized_trigger_source_timestamp(source),
        delivery_timestamp)
end

function _process_trigger_topology!(
    prepared::PreparedPlantEventLoop{<:Any,<:Any,<:PreparedTriggerTopology},
    state::PlantEventLoopState{<:TriggerTopologyState},
    workspace::PlantEventLoopWorkspace{<:TriggerTopologyWorkspace},
    claim::EventClaim)
    topology = prepared.trigger_topology
    trigger_state = state.trigger
    source = next_trigger_source(topology, trigger_state)
    delivery = next_trigger_delivery(topology, trigger_state)
    source_due = delivery === nothing ||
        realized_trigger_source_timestamp(source) <=
            delivered_trigger_edge(delivery).timestamp
    activated_slot = UInt32(0)
    activation_timestamp = zero(PlantTimestamp)
    if source_due
        realized_trigger_source_timestamp(source) == claim.key.timestamp ||
            _plant_event_loop_error(:trigger_schedule,
                "trigger source does not match its scheduler claim")
        realize_next_trigger_source!(workspace.trigger, topology,
            trigger_state)
    else
        delivered = delivered_trigger_edge(delivery)
        delivered.timestamp == claim.key.timestamp ||
            _plant_event_loop_error(:trigger_schedule,
                "trigger delivery does not match its scheduler claim")
        activated_slot = _triggered_acquisition_slot(prepared,
            trigger_delivery_consumer(delivery))
        acquisition = _event_acquisition_binding(prepared, activated_slot)
        _require_inactive_event_generator(prepared, state,
            acquisition.start_handle, delivered.timestamp)
        pop_next_trigger_delivery!(workspace.delivery, topology,
            trigger_state) || _plant_event_loop_error(:trigger_schedule,
            "due trigger delivery disappeared before removal")
        activation_timestamp = delivered.timestamp
    end
    next_timestamp = _next_trigger_action_timestamp(topology, trigger_state)
    reschedule_event!(prepared.scheduler, state.scheduler, claim,
        next_timestamp)
    if !iszero(activated_slot)
        acquisition = _event_acquisition_binding(prepared, activated_slot)
        activate_event_generator!(prepared.scheduler, state.scheduler,
            acquisition.start_handle, activation_timestamp)
    end
    return nothing
end

function _process_trigger_topology!(
    ::PreparedPlantEventLoop{<:Any,<:Any,<:_NoPreparedTriggerTopology},
    ::PlantEventLoopState{<:_NoTriggerTopologyState},
    ::PlantEventLoopWorkspace{<:_NoTriggerTopologyWorkspace}, ::EventClaim)
    _plant_event_loop_error(:invalid_action,
        "trigger action exists without a prepared trigger topology")
end

function _preflight_event_path(path::_PreparedPlantEventPath,
    atmosphere)
    path.path.atmosphere === atmosphere || _plant_event_loop_error(
        :prepared_binding,
        "event path does not retain the plant atmosphere")
    _require_current_path_binding(path.path)
    _require_rng_owner_binding(path.rngs, path.path)
    return nothing
end

function _preflight_event_acquisition(
    acquisition::_PreparedPlantEventAcquisition,
    state::_DetectorEventLifecycleState)
    _require_event_lifecycle_binding(acquisition.lifecycle, state)
    return nothing
end

@inline function _require_event_lifecycle_binding(
    prepared::_PreparedDetectorEventLifecycle,
    state::_DetectorEventLifecycleState)
    getfield(state, :binding) === getfield(prepared, :binding) ||
        _plant_event_loop_error(:foreign_state,
            "detector lifecycle state belongs to another prepared acquisition")
    det = getfield(prepared, :detector)
    plan = getfield(prepared, :plan)
    det.params === plan.detector_params &&
        det.state === plan.detector_state &&
        det.state.frame === plan.detector_frame &&
        det.state.readout_products === getfield(prepared, :readout_products) &&
        typeof(backend(det)) === typeof(plan.detector_backend) ||
        _plant_event_loop_error(:prepared_binding,
            "detector lifecycle storage changed after event-loop preparation")
    return nothing
end

@inline function _preflight_event_integration_to(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp)
    state.status == DetectorExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "global-shutter integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    _require_detector_event_progress(prepared, state)
    timestamp <= state.exposure_close ||
        _detector_acquisition_event_error(:interval_after_close,
            "integration interval extends beyond exposure close")
    _require_interval_before_next_read(prepared, state, timestamp)
    return nothing
end

@inline function _preflight_event_integration_to(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp)
    state.status == DetectorExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "rolling-shutter integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    state.opened_bands == state.closed_bands && return nothing
    _require_rolling_shutter_progress(prepared, state)
    next_transition = _next_rolling_transition_timestamp(prepared, state)
    next_transition === nothing || timestamp <= next_transition ||
        _detector_acquisition_event_error(:missed_band_transition,
            "rolling integration crosses a pending row-band transition")
    return nothing
end

@inline function _preflight_event_integration_to(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp)
    state.image_status == _FrameTransferImageActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "frame-transfer integration target precedes detector progress")
    state.integrated_through == timestamp && return nothing
    _require_frame_transfer_progress(prepared, state)
    timestamp <= state.exposure_close ||
        _detector_acquisition_event_error(:interval_after_close,
            "frame-transfer integration extends beyond exposure close")
    return nothing
end

function _mark_due_event_paths!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp)
    fill!(workspace.due_paths, false)
    scan_due_events!(workspace.scheduler, prepared.scheduler,
        state.scheduler)
    count = workspace.scheduler.due_count
    @inbounds for index in 1:count
        slot = Int(workspace.scheduler.due_slots[index])
        definition = prepared.scheduler.definitions[slot]
        definition.phase == OpticalSamplePhase || continue
        state.scheduler.cursors[slot].next_timestamp == timestamp || continue
        action = prepared.actions[slot]
        workspace.due_paths[Int(action.owner_slot)] = true
    end
    any(workspace.due_paths) || _plant_event_loop_error(
        :invalid_action, "optical sample phase has no due path")
    return nothing
end

function _preflight_atmosphere_time(atmosphere::AbstractTimedAtmosphere,
    timestamp::PlantTimestamp)
    timeline = atmosphere_timeline(atmosphere)
    T = typeof(timeline.model_time)
    target = plant_time_seconds(timestamp, T)
    timeline.initialized && target < timeline.model_time &&
        _plant_event_loop_error(:atmosphere_time_regression,
            "due optical sample precedes the current atmosphere epoch")
    return target
end

function _preflight_due_path_consumers(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    due_paths::Memory{Bool}, timestamp::PlantTimestamp)
    @inbounds for index in eachindex(prepared.acquisitions)
        acquisition = prepared.acquisitions[index]
        due_paths[Int(acquisition.path_slot)] || continue
        acquisition_state = state.acquisitions[index]
        _preflight_event_acquisition(acquisition, acquisition_state)
        _preflight_event_integration_to(acquisition.lifecycle,
            acquisition_state, timestamp)
    end
    return nothing
end

function _integrate_due_path_consumers!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    due_paths::Memory{Bool}, timestamp::PlantTimestamp)
    @inbounds for index in eachindex(prepared.acquisitions)
        acquisition = prepared.acquisitions[index]
        due_paths[Int(acquisition.path_slot)] || continue
        acquisition_state = state.acquisitions[index]
        _integrate_event_acquisition_to!(acquisition.lifecycle,
            acquisition_state, timestamp, acquisition.rng)
    end
    return nothing
end

function _validate_due_path_materializations!(
    prepared::PreparedPlantEventLoop, due_paths::Memory{Bool}, atmosphere,
    epoch)
    @inbounds for index in eachindex(prepared.paths)
        due_paths[index] || continue
        path = prepared.paths[index].path
        validate_path_materialization(path.materialization, path.input,
            atmosphere, epoch)
    end
    return nothing
end

function _materialize_due_paths!(prepared::PreparedPlantEventLoop,
    due_paths::Memory{Bool}, atmosphere, epoch)
    @inbounds for index in eachindex(prepared.paths)
        due_paths[index] || continue
        binding = prepared.paths[index]
        path = binding.path
        materialize_path_input_rngs!(path.materialization, path.input,
            atmosphere, epoch, binding.rngs)
    end
    return nothing
end

function _execute_due_paths!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, due_paths::Memory{Bool})
    @inbounds for index in eachindex(prepared.paths)
        due_paths[index] || continue
        binding = prepared.paths[index]
        execute_path!(binding.path, binding.rngs)
        state.path_sampled[index] = true
    end
    return nothing
end

function _resolve_due_path_claims!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp)
    while true
        count = scan_due_events!(workspace.scheduler, prepared.scheduler,
            state.scheduler)
        iszero(count) && return nothing
        workspace.scheduler.due_timestamp == timestamp || return nothing
        key = due_event_key(workspace.scheduler, prepared.scheduler,
            state.scheduler, 1)
        key.phase == OpticalSamplePhase || return nothing
        claim = claim_next_event!(workspace.scheduler, prepared.scheduler,
            state.scheduler)
        claim === nothing && _plant_event_loop_error(:invalid_action,
            "due optical path disappeared before claim")
        action = _event_action(prepared, claim)
        path = _event_path_binding(prepared, action.owner_slot)
        reschedule_periodic_event!(prepared.scheduler, state.scheduler,
            claim, path.schedule; origin=path.origin)
    end
end

function _process_optical_path_batch!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp)
    atmosphere = prepared.atmosphere
    _mark_due_event_paths!(prepared, state, workspace, timestamp)
    @inbounds for index in eachindex(prepared.paths)
        workspace.due_paths[index] || continue
        _preflight_event_path(prepared.paths[index], atmosphere)
    end
    target_time = _preflight_atmosphere_time(atmosphere, timestamp)
    _preflight_due_path_consumers(prepared, state, workspace.due_paths,
        timestamp)
    _integrate_due_path_consumers!(prepared, state, workspace.due_paths,
        timestamp)
    epoch = advance_to!(atmosphere, target_time, prepared.atmosphere_rng)
    _validate_due_path_materializations!(prepared, workspace.due_paths,
        atmosphere, epoch)
    _materialize_due_paths!(prepared, workspace.due_paths, atmosphere,
        epoch)
    _execute_due_paths!(prepared, state, workspace.due_paths)
    _resolve_due_path_claims!(prepared, state, workspace, timestamp)
    return nothing
end

function _process_ordinary_event!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    claim::EventClaim, action::_PlantEventAction)
    kind = action.kind
    kind == _TriggerTopologyAction &&
        return _process_trigger_topology!(prepared, state, workspace, claim)
    kind == _AcquisitionBoundaryAction &&
        return _process_acquisition_boundary!(prepared, state, claim, action)
    kind == _AcquisitionStartAction &&
        return _process_acquisition_start!(prepared, state, claim, action)
    kind == _RollingBandOpenAction &&
        return _process_rolling_band_open!(prepared, state, claim, action)
    kind == _AcquisitionReadoutAction &&
        return _process_acquisition_readout!(prepared, state, claim, action)
    kind == _AcquisitionReadinessAction &&
        return _process_acquisition_readiness!(prepared, state, claim, action)
    _plant_event_loop_error(:invalid_action,
        "prepared plant event has an unknown action kind")
end

function next_plant_event_timestamp(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    count = scan_due_events!(workspace.scheduler, prepared.scheduler,
        state.scheduler)
    iszero(count) && return nothing
    return workspace.scheduler.due_timestamp
end

"""
    step_plant_events!(prepared, state, workspace)

Process every event at the next canonical plant timestamp in causal phase and
prepared-ordinal order. All due optical paths at that timestamp are integrated,
advanced, materialized, and formed as one bounded batch. Returns the processed
timestamp, or `nothing` when every generator is inactive.
"""
function step_plant_events!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    timestamp = next_plant_event_timestamp(prepared, state, workspace)
    timestamp === nothing && return nothing
    while true
        next_timestamp = next_plant_event_timestamp(prepared, state,
            workspace)
        next_timestamp === nothing && break
        next_timestamp == timestamp || break
        key = due_event_key(workspace.scheduler, prepared.scheduler,
            state.scheduler, 1)
        if key.phase == OpticalSamplePhase
            _process_optical_path_batch!(prepared, state, workspace,
                timestamp)
            continue
        end
        claim = claim_next_event!(workspace.scheduler, prepared.scheduler,
            state.scheduler)
        claim === nothing && _plant_event_loop_error(:invalid_action,
            "due plant event disappeared before claim")
        action = _event_action(prepared, claim)
        _process_ordinary_event!(prepared, state, workspace, claim, action)
    end
    return timestamp
end

@inline function _checked_event_step_limit(limit::Integer)
    limit >= 0 || _plant_event_loop_error(:invalid_step_limit,
        "max_timestamps must be nonnegative")
    limit <= typemax(Int) || _plant_event_loop_error(:invalid_step_limit,
        "max_timestamps exceeds Int range")
    return Int(limit)
end

@inline _checked_event_step_limit(::Bool) =
    _plant_event_loop_error(:invalid_step_limit,
        "max_timestamps must be an integer count, not Bool")

"""Process scheduled timestamps through an inclusive finite plant horizon."""
function run_plant_events_until!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    stop::PlantTimestamp; max_timestamps::Integer=typemax(Int))
    limit = _checked_event_step_limit(max_timestamps)
    count = 0
    while count < limit
        timestamp = next_plant_event_timestamp(prepared, state, workspace)
        timestamp === nothing && break
        timestamp <= stop || break
        step_plant_events!(prepared, state, workspace)
        count += 1
    end
    return count
end
