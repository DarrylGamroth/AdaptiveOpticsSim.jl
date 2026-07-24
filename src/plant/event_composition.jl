const _PLANT_EVENT_LOOP_COMPONENT = :plant_event_loop

@noinline function _plant_event_loop_error(reason::Symbol,
    message::AbstractString)
    throw(PlantScheduleError(_PLANT_EVENT_LOOP_COMPONENT, reason,
        String(message)))
end

abstract type AbstractAcquisitionStartDefinition end

"""Periodic start recurrence for one acquisition owner."""
struct PeriodicAcquisitionStart <: AbstractAcquisitionStartDefinition
    schedule::PeriodicSchedule
    origin::PlantTimestamp
end

PeriodicAcquisitionStart(schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp)) =
    PeriodicAcquisitionStart(schedule, origin)

"""Binding from one delivered trigger consumer to an acquisition start."""
struct TriggeredAcquisitionStart <: AbstractAcquisitionStartDefinition
    consumer::TriggerConsumerID
end

const _PreparedAcquisitionEventLifecycle = Union{
    PreparedGlobalShutterAcquisition,
    PreparedRollingShutterAcquisition,
    PreparedFrameTransferAcquisition,
    PreparedDirectMeasurementAcquisition,
}
const _AcquisitionEventLifecycleState = Union{
    GlobalShutterAcquisitionState,
    RollingShutterAcquisitionState,
    FrameTransferAcquisitionState,
    DirectMeasurementAcquisitionState,
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
Lifecycle and start-source declaration for one prepared acquisition owner.
"""
struct AcquisitionEventDefinition{D<:AbstractAcquisitionLifecycleDefinition,
    S<:AbstractAcquisitionStartDefinition}
    acquisition::AcquisitionID
    lifecycle::D
    start::S
end

AcquisitionEventDefinition(acquisition,
    lifecycle::AbstractAcquisitionLifecycleDefinition,
    start::AbstractAcquisitionStartDefinition) = AcquisitionEventDefinition(
    _as_acquisition_id(acquisition), lifecycle, start)

@inline _require_optical_sample_definition(value::OpticalSampleDefinition) =
    value
@inline _require_acquisition_event_definition(
    value::AcquisitionEventDefinition) =
    value
@inline _require_autonomous_periodic_optic_definition(
    value::AutonomousPeriodicOpticDefinition) = value

function _require_optical_sample_definition(value)
    _plant_event_loop_error(:invalid_definition,
        "optical sample entries must be OpticalSampleDefinition values; got $(typeof(value))")
end

function _require_acquisition_event_definition(value)
    _plant_event_loop_error(:invalid_definition,
        "acquisition event entries must be AcquisitionEventDefinition values; got $(typeof(value))")
end

function _require_autonomous_periodic_optic_definition(value)
    _plant_event_loop_error(:invalid_definition,
        "autonomous optic entries must be AutonomousPeriodicOpticDefinition values; got $(typeof(value))")
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
    PlantEventLoopDefinition(optical_samples, acquisition_events;
        trigger_topology=nothing, autonomous_optics=())

Cold, finite declaration of periodic path samples and independently periodic or
trigger-driven acquisition lifecycles plus optional trigger-relative
autonomous optical devices. Preparation flattens scheduled values into one
fixed-capacity serial scheduler; the declaration stores no mutable cursor,
acquisition state, waveform state, or run-length event list.
"""
struct PlantEventLoopDefinition
    optical_samples::Tuple
    acquisition_events::Tuple
    trigger_topology::Union{Nothing,PreparedTriggerTopology}
    autonomous_optics::Tuple
end

function PlantEventLoopDefinition(optical_samples, acquisition_events;
    trigger_topology=nothing, autonomous_optics=())
    samples = _event_definition_tuple(optical_samples,
        _require_optical_sample_definition)
    events = _event_definition_tuple(acquisition_events,
        _require_acquisition_event_definition)
    autonomous = _event_definition_tuple(autonomous_optics,
        _require_autonomous_periodic_optic_definition)
    isempty(samples) && _plant_event_loop_error(:empty_paths,
        "plant event loop requires at least one optical sample definition")
    isempty(events) && _plant_event_loop_error(:empty_acquisitions,
        "plant event loop requires at least one acquisition event definition")
    topology = _require_prepared_trigger_topology(trigger_topology)
    return PlantEventLoopDefinition(samples, events, topology, autonomous)
end

struct _NoPreparedTriggerTopology end
struct _NoTriggerTopologyState end
struct _NoTriggerTopologyWorkspace end

mutable struct _PlantEventLoopBinding end

@enum _PlantEventActionKind::UInt8 begin
    _TriggerTopologyAction = 0x01
    _CommandEndpointAction = 0x02
    _AcquisitionBoundaryAction = 0x03
    _AcquisitionStartAction = 0x04
    _RollingBandOpenAction = 0x05
    _OpticalPathSampleAction = 0x06
    _AcquisitionReadoutAction = 0x07
    _AcquisitionReadinessAction = 0x08
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
    requires_full_optical::Bool
end

struct _PreparedPlantEventAcquisition
    id::AcquisitionID
    lifecycle::_PreparedAcquisitionEventLifecycle
    product::Union{AbstractArray,WFSMeasurement}
    sample_provider::Union{Nothing,PreparedLinearReducedOrderEventProvider}
    rng::Xoshiro
    start::_AcquisitionStartDefinition
    path_slot::UInt32
    start_handle::EventGeneratorHandle
    boundary_handle::EventGeneratorHandle
    band_open_handle::EventGeneratorHandle
    readout_handle::EventGeneratorHandle
    readiness_handle::EventGeneratorHandle
end

struct _PreparedPlantEventCommandEndpoint{
    B<:_PreparedPlantCommandEndpoint}
    binding::B
    handle::EventGeneratorHandle
end

"""
Explicit atomic command submission for two or more distinct endpoints.

All members request one common effective plant timestamp. Construction alone
does not admit commands, assign transaction identity, group optical surfaces,
or imply a transport encoding.
"""
struct _PlantCommandTransactionToken end
const _PLANT_COMMAND_TRANSACTION_TOKEN = _PlantCommandTransactionToken()

struct PlantCommandTransaction{C<:Tuple}
    commands::C

    PlantCommandTransaction(commands::C,
        ::_PlantCommandTransactionToken) where {C<:Tuple} =
        new{C}(commands)
end

function PlantCommandTransaction(
    commands::Vararg{PlantCommand,N}) where {N}
    N >= 2 || _command_admission_error(:transaction,
        :invalid_member_count,
        "an atomic command transaction requires at least two commands")
    first_timestamp = command_requested_effective_timestamp(first(commands))
    @inbounds for right in 2:N
        command = commands[right]
        command_requested_effective_timestamp(command) == first_timestamp ||
            _command_admission_error(:transaction,
                :effective_timestamp_mismatch,
                "every atomic transaction member must request the same " *
                "effective plant timestamp")
        endpoint = command_endpoint_id(command)
        for left in 1:(right - 1)
            command_endpoint_id(commands[left]) == endpoint &&
                _command_admission_error(:transaction,
                    :duplicate_endpoint,
                    "atomic transaction contains command endpoint $endpoint " *
                    "more than once")
        end
    end
    return PlantCommandTransaction(commands, _PLANT_COMMAND_TRANSACTION_TOKEN)
end

"""Immediate all-member result of one atomic transaction admission."""
struct PlantCommandTransactionAdmission
    transaction::UInt64
    status::CommandAdmissionStatus
    member_count::UInt32
    scheduled_timestamp::Union{Nothing,PlantTimestamp}
end

@inline command_admission_status(
    admission::PlantCommandTransactionAdmission) = admission.status
@inline command_transaction_member_count(
    admission::PlantCommandTransactionAdmission) =
    admission.member_count
@inline function command_transaction_id(
    admission::PlantCommandTransactionAdmission)
    iszero(admission.transaction) && return nothing
    return admission.transaction
end
@inline command_scheduled_timestamp(
    admission::PlantCommandTransactionAdmission) =
    admission.scheduled_timestamp

struct _CommandTransactionAdmissionPlan
    presentation::CommandPresentationID
    sequence_class::CommandSequenceClass
    scheduled_timestamp::PlantTimestamp
    slot::UInt32
    generation::UInt64
end

struct _CommandTransactionPolicyFailure <: Exception
    error::PlantCommandError
end

"""Run-immutable, fixed-capacity deterministic plant-event composition."""
struct PreparedPlantEventLoop{A<:AbstractTimedAtmosphere,R,T}
    binding::_PlantEventLoopBinding
    atmosphere::A
    atmosphere_rng::R
    scheduler::PreparedEventScheduler
    actions::Memory{_PlantEventAction}
    optics::Memory{PreparedControllableOptic}
    command_endpoints::Memory{_PreparedPlantEventCommandEndpoint}
    paths::Memory{_PreparedPlantEventPath}
    acquisitions::Memory{_PreparedPlantEventAcquisition}
    autonomous_optics::Memory{_PreparedAutonomousPeriodicOptic}
    trigger_topology::T
end

@inline plant_event_path_count(prepared::PreparedPlantEventLoop) =
    length(prepared.paths)
@inline plant_event_acquisition_count(prepared::PreparedPlantEventLoop) =
    length(prepared.acquisitions)
@inline plant_event_generator_count(prepared::PreparedPlantEventLoop) =
    event_generator_count(prepared.scheduler)
@inline plant_event_controllable_optic_count(
    prepared::PreparedPlantEventLoop) = length(prepared.optics)
@inline plant_event_command_endpoint_count(
    prepared::PreparedPlantEventLoop) = length(prepared.command_endpoints)
@inline plant_event_autonomous_optic_count(
    prepared::PreparedPlantEventLoop) = length(prepared.autonomous_optics)

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

function _sorted_acquisition_event_definitions(definitions::Tuple)
    values = Any[definitions...]
    sort!(values; by=definition -> definition.acquisition.name)
    @inbounds for index in 2:length(values)
        values[index - 1].acquisition == values[index].acquisition &&
            _plant_event_loop_error(:duplicate_acquisition,
                "acquisition $(values[index].acquisition) has more than one lifecycle")
    end
    return values
end

function _sorted_autonomous_periodic_optic_definitions(
    definitions::Tuple)
    values = Any[definitions...]
    sort!(values; by=definition -> definition.optic.name)
    @inbounds for index in 2:length(values)
        values[index - 1].optic == values[index].optic &&
            _plant_event_loop_error(:duplicate_autonomous_optic,
                "controllable optic $(values[index].optic) has more than one autonomous waveform binding")
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
        "event owner references path $id without an optical sample schedule")
end

function _event_path_slot_from_definitions(definitions, id::OpticalPathID)
    @inbounds for index in eachindex(definitions)
        definitions[index].path == id && return index
    end
    _plant_event_loop_error(:missing_path_schedule,
        "acquisition references path $id without an optical sample schedule")
end

function _trigger_consumer_exists(topology::PreparedTriggerTopology,
    id::TriggerConsumerID)
    @inbounds for consumer in topology.consumers
        consumer.id == id && return true
    end
    return false
end

function _trigger_source_exists(topology::PreparedTriggerTopology,
    id::TriggerSourceID)
    @inbounds for source in topology.sources
        source.id == id && return true
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

@inline _definition_trigger_consumer(
    definition::AcquisitionEventDefinition) =
    _start_trigger_consumer(definition.start)
@inline _start_trigger_consumer(::PeriodicAcquisitionStart) = nothing
@inline _start_trigger_consumer(start::TriggeredAcquisitionStart) =
    start.consumer

@inline _definition_trigger_consumer(
    definition::AutonomousPeriodicOpticDefinition) =
    _phase_reference_trigger_consumer(definition.phase_reference)
@inline _phase_reference_trigger_consumer(
    ::FreeRunningPhaseReference) = nothing
@inline _phase_reference_trigger_consumer(
    ::TriggerSourcePhaseReference) = nothing
@inline _phase_reference_trigger_consumer(
    reference::TriggerResetPhaseReference) = reference.consumer

function _require_autonomous_trigger_topology(
    definition::AutonomousPeriodicOpticDefinition{
        <:FreeRunningPhaseReference}, ::Any)
    return nothing
end

function _require_autonomous_trigger_topology(
    definition::AutonomousPeriodicOpticDefinition{
        <:TriggerSourcePhaseReference}, ::Nothing)
    _plant_event_loop_error(:missing_trigger_topology,
        "autonomous optic $(definition.optic) uses a trigger-source phase reference without a prepared trigger topology")
end

function _require_autonomous_trigger_topology(
    definition::AutonomousPeriodicOpticDefinition{
        <:TriggerSourcePhaseReference},
    topology::PreparedTriggerTopology)
    reference = definition.phase_reference
    _trigger_source_exists(topology, reference.source) ||
        _plant_event_loop_error(:unknown_trigger_source,
            "autonomous optic $(definition.optic) references unknown trigger source $(reference.source)")
    return nothing
end

function _require_autonomous_trigger_topology(
    definition::AutonomousPeriodicOpticDefinition{
        <:TriggerResetPhaseReference}, ::Nothing)
    _plant_event_loop_error(:missing_trigger_topology,
        "autonomous optic $(definition.optic) uses a delivered-trigger phase reset without a prepared trigger topology")
end

function _require_autonomous_trigger_topology(
    definition::AutonomousPeriodicOpticDefinition{
        <:TriggerResetPhaseReference},
    topology::PreparedTriggerTopology)
    reference = definition.phase_reference
    _trigger_consumer_exists(topology, reference.consumer) ||
        _plant_event_loop_error(:unknown_trigger_consumer,
            "autonomous optic $(definition.optic) references unknown trigger consumer $(reference.consumer)")
    return nothing
end

function _record_unique_trigger_consumer!(
    seen::Set{TriggerConsumerID}, definition)
    consumer = _definition_trigger_consumer(definition)
    consumer === nothing && return nothing
    consumer in seen && _plant_event_loop_error(:duplicate_trigger_consumer,
        "trigger consumer $consumer is bound to more than one event owner")
    push!(seen, consumer)
    return nothing
end

function _require_unique_trigger_consumers(
    acquisition_definitions, autonomous_definitions)
    seen = Set{TriggerConsumerID}()
    foreach(definition -> _record_unique_trigger_consumer!(seen, definition),
        acquisition_definitions)
    foreach(definition -> _record_unique_trigger_consumer!(seen, definition),
        autonomous_definitions)
    return nothing
end

@inline _require_bound_trigger_consumers(
    ::Any, ::Any, ::Nothing) = nothing

function _require_bound_trigger_consumers(
    acquisition_definitions, autonomous_definitions,
    topology::PreparedTriggerTopology)
    @inbounds for consumer in topology.consumers
        bound = false
        for definition in acquisition_definitions
            if _definition_trigger_consumer(definition) == consumer.id
                bound = true
                break
            end
        end
        if !bound
            for definition in autonomous_definitions
                if _definition_trigger_consumer(definition) == consumer.id
                    bound = true
                    break
                end
            end
        end
        bound || _plant_event_loop_error(:unbound_trigger_consumer,
            "trigger consumer $(consumer.id) has no acquisition or autonomous-optic binding")
    end
    return nothing
end

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
    ::AbstractPreparedAcquisitionLifecycle,
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
    prepared::PreparedDirectMeasurementAcquisition,
    start::PeriodicAcquisitionStart)
    occupied = prepared.definition.exposure_duration +
        prepared.definition.readout_duration +
        prepared.definition.readiness_delay
    schedule_period(start.schedule) > occupied ||
        _plant_event_loop_error(:acquisition_period,
            "direct-measurement period must be strictly later than acquisition readiness")
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
    sample_definitions, acquisition_definitions, command_endpoints,
    trigger_topology)
    if trigger_topology !== nothing
        push!(definitions, EventGeneratorDefinition(
            _initial_trigger_timestamp(trigger_topology),
            TriggerUpdatePhase, 1))
    end
    for (index, _) in enumerate(command_endpoints)
        push!(definitions, EventGeneratorDefinition(zero(PlantTimestamp),
            CommandApplicationPhase, index; active=false))
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
    ::Val{CommandApplicationPhase}, ordinal::UInt32) =
    _PlantEventAction(_CommandEndpointAction, ordinal)
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

function _prepare_event_controllable_optics(plant::PreparedPlant)
    source = getfield(plant, :controllable_optics)
    optics = Memory{PreparedControllableOptic}(undef, length(source))
    @inbounds for index in eachindex(source)
        optics[index] = source[index]
    end
    return optics
end

function _prepare_event_command_endpoints(plant::PreparedPlant,
    scheduler::PreparedEventScheduler)
    source = getfield(plant, :command_endpoints)
    endpoints = Memory{_PreparedPlantEventCommandEndpoint}(
        undef, length(source))
    @inbounds for index in eachindex(source)
        endpoints[index] = _PreparedPlantEventCommandEndpoint(source[index],
            event_generator_handle(scheduler, CommandApplicationPhase,
                index))
    end
    return endpoints
end

@inline _provider_requires_full_optical(::FullOpticalProviderStyle) = true
@inline _provider_requires_full_optical(
    ::CommandResponsiveReducedOrderProviderStyle) = false
@inline _provider_requires_full_optical(::SyntheticReplayProviderStyle) =
    false

function _event_path_requires_full_optical(id::OpticalPathID, owners)
    @inbounds for owner in owners
        acquisition_path_id(owner.definition) == id || continue
        _provider_requires_full_optical(acquisition_provider_style(owner)) &&
            return true
    end
    return false
end

@inline _require_linear_reduced_order_provider(
    implementation::PreparedLinearReducedOrderProvider) = implementation

function _require_linear_reduced_order_provider(implementation)
    _plant_event_loop_error(:unsupported_acquisition,
        "event composition currently supports the built-in linear reduced-order provider; got $(typeof(implementation))")
end

@inline _require_direct_event_measurement(
    measurement::WFSMeasurement) = measurement

function _require_direct_event_measurement(measurement)
    _plant_event_loop_error(:unsupported_acquisition,
        "linear reduced-order event acquisition requires a WFSMeasurement product; got $(typeof(measurement))")
end

function _prepare_event_paths(plant::PreparedPlant, definitions, owners,
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
            definition.schedule, definition.origin, handle,
            _event_path_requires_full_optical(definition.path, owners))
    end
    return paths
end

function _event_controllable_optic_slot(optics,
    id::ControllableOpticID)
    @inbounds for index in eachindex(optics)
        controllable_optic_id(optics[index].definition) == id &&
            return index
    end
    _plant_event_loop_error(:unknown_controllable_optic,
        "autonomous waveform references unknown controllable optic $id")
end

@inline function _require_autonomous_execution_role(
    ::AutonomousPathExecutionRole, ::ControllableOpticID)
    return nothing
end

function _require_autonomous_execution_role(
    ::PupilSurfaceExecutionRole, id::ControllableOpticID)
    _plant_event_loop_error(:invalid_autonomous_optic,
        "controllable optic $id uses common-pupil surface execution and cannot be bound as a path-local autonomous waveform")
end

function _require_autonomous_execution_role(
    role, id::ControllableOpticID)
    _plant_event_loop_error(:unsupported_optic_execution_role,
        "controllable optic $id declares unsupported execution role $(typeof(role))")
end

function _definition_binds_autonomous_optic(definitions,
    id::ControllableOpticID)
    @inbounds for definition in definitions
        definition.optic == id && return true
    end
    return false
end

@inline function _require_autonomous_definition_for_role(
    ::PupilSurfaceExecutionRole, ::ControllableOpticID, ::Any)
    return nothing
end

function _require_autonomous_definition_for_role(
    ::AutonomousPathExecutionRole, id::ControllableOpticID, definitions)
    _definition_binds_autonomous_optic(definitions, id) ||
        _plant_event_loop_error(:missing_autonomous_binding,
            "path-local autonomous controllable optic $id has no AutonomousPeriodicOpticDefinition")
    return nothing
end

function _require_autonomous_definition_for_role(
    role, id::ControllableOpticID, ::Any)
    _plant_event_loop_error(:unsupported_optic_execution_role,
        "controllable optic $id declares unsupported execution role $(typeof(role))")
end

function _require_all_autonomous_optics_bound(optics, definitions)
    @inbounds for optic in optics
        id = controllable_optic_id(optic.definition)
        role = controllable_optic_execution_role(optic.implementation)
        _require_autonomous_definition_for_role(role, id, definitions)
    end
    return nothing
end

function _require_unique_autonomous_optic_couplings(bindings)
    @inbounds for right in 2:length(bindings), left in 1:(right - 1)
        _autonomous_optic_couplings_conflict(
            bindings[left].coupling, bindings[right].coupling) || continue
        _plant_event_loop_error(:conflicting_autonomous_coupling,
            "autonomous optics $(bindings[left].id) and " *
            "$(bindings[right].id) target the same exclusive prepared " *
            "optical coupling")
    end
    return nothing
end

function _prepare_event_autonomous_optics(definitions, optics, paths,
    topology)
    bindings = Memory{_PreparedAutonomousPeriodicOptic}(undef,
        length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        _require_autonomous_trigger_topology(definition, topology)
        optic_slot = _event_controllable_optic_slot(optics,
            definition.optic)
        optic = optics[optic_slot]
        _require_autonomous_execution_role(
            controllable_optic_execution_role(optic.implementation),
            definition.optic)
        path_slot = _event_path_slot(paths, definition.path)
        path = paths[path_slot]
        path.requires_full_optical || _plant_event_loop_error(
            :autonomous_path_without_full_optics,
            "autonomous optic $(definition.optic) targets path $(definition.path) without a full-optical acquisition")
        coupling = prepare_autonomous_periodic_optic(
            optic.implementation, path.path, definition.fidelity)
        bindings[index] = _PreparedAutonomousPeriodicOptic(
            definition.optic, definition.path, UInt32(optic_slot),
            UInt32(path_slot), optic.implementation, coupling,
            definition.phase_reference, definition.fidelity)
    end
    _require_all_autonomous_optics_bound(optics, definitions)
    _require_unique_autonomous_optic_couplings(bindings)
    return bindings
end

function _prepare_event_acquisition_lifecycle(
    plant::PreparedPlant,
    owner::PreparedAcquisitionOwner,
    definition::AbstractDetectorAcquisitionLifecycleDefinition,
    ::FullOpticalProviderStyle)
    execution = _event_frame_execution(owner)
    lifecycle = _prepare_detector_event_lifecycle(execution,
        owner.path_result, definition)
    return lifecycle, acquisition_observation(owner), nothing
end

function _prepare_event_acquisition_lifecycle(
    plant::PreparedPlant,
    owner::PreparedAcquisitionOwner,
    definition::DirectMeasurementAcquisitionDefinition,
    ::CommandResponsiveReducedOrderProviderStyle)
    implementation = _require_linear_reduced_order_provider(
        owner.provider.implementation)
    measurement = _require_direct_event_measurement(
        acquisition_measurement(owner))
    lifecycle = prepare_direct_measurement_acquisition(measurement,
        definition)
    sample_provider = prepare_linear_reduced_order_event_provider(
        implementation, getfield(plant, :command_endpoints))
    return lifecycle, measurement, sample_provider
end

function _prepare_event_acquisition_lifecycle(
    ::PreparedPlant, owner::PreparedAcquisitionOwner,
    definition::AbstractAcquisitionLifecycleDefinition,
    ::CommandResponsiveReducedOrderProviderStyle)
    _plant_event_loop_error(:unsupported_acquisition,
        "command-responsive reduced-order acquisition $(acquisition_id(owner.definition)) requires DirectMeasurementAcquisitionDefinition; got $(typeof(definition))")
end

function _prepare_event_acquisition_lifecycle(
    ::PreparedPlant, owner::PreparedAcquisitionOwner,
    definition::AbstractAcquisitionLifecycleDefinition,
    ::SyntheticReplayProviderStyle)
    _plant_event_loop_error(:unsupported_acquisition,
        "scheduled synthetic/replay acquisition lifecycles are not supported")
end

function _prepare_event_acquisition_lifecycle(
    ::PreparedPlant, owner::PreparedAcquisitionOwner,
    definition::DirectMeasurementAcquisitionDefinition,
    ::FullOpticalProviderStyle)
    _plant_event_loop_error(:unsupported_acquisition,
        "full-optical acquisition $(acquisition_id(owner.definition)) requires a detector lifecycle, not DirectMeasurementAcquisitionDefinition")
end

function _prepare_event_acquisition_parts(plant::PreparedPlant,
    definitions, path_definitions, topology)
    Base.@nospecialize plant
    owners = Any[]
    lifecycles = Any[]
    products = Any[]
    sample_providers = Any[]
    rngs = Any[]
    path_slots = Int[]
    for definition in definitions
        _require_start_trigger_topology(definition.start, topology)
        owner = _event_prepared_acquisition(plant, definition.acquisition)
        lifecycle, product, sample_provider =
            _prepare_event_acquisition_lifecycle(plant, owner,
                definition.lifecycle, acquisition_provider_style(owner))
        _require_periodic_start_spacing(lifecycle, definition.start)
        acquisition_rngs = _prepared_event_acquisition_rngs(plant, owner)
        _require_rng_owner_binding(acquisition_rngs, owner)
        push!(owners, owner)
        push!(lifecycles, lifecycle)
        push!(products, product)
        push!(sample_providers, sample_provider)
        push!(rngs, rng_stream_state(acquisition_rngs, Val(:detector)))
        push!(path_slots, _event_path_slot_from_definitions(path_definitions,
            acquisition_path_id(owner.definition)))
    end
    return owners, lifecycles, products, sample_providers, rngs, path_slots
end

function _prepare_event_acquisitions(definitions, lifecycles, products,
    sample_providers, rngs, path_slots,
    scheduler::PreparedEventScheduler)
    acquisitions = Memory{_PreparedPlantEventAcquisition}(undef,
        length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        acquisitions[index] = _PreparedPlantEventAcquisition(
            definition.acquisition, lifecycles[index],
            products[index], sample_providers[index], rngs[index],
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

Bind a prepared plant to a flat deterministic scheduler, exact
owner-derived RNG streams, bounded detector lifecycle state, and an optional
prepared trigger topology. Prepared controllable optics and independently
timed command endpoints join the reserved command phase. Preparation allocates;
repeated stepping does not grow the registry or materialize future events.
"""
function prepare_plant_event_loop(plant::PreparedPlant,
    definition::PlantEventLoopDefinition)
    Base.@nospecialize plant
    sample_definitions = _sorted_optical_sample_definitions(
        definition.optical_samples)
    acquisition_definitions = _sorted_acquisition_event_definitions(
        definition.acquisition_events)
    autonomous_definitions =
        _sorted_autonomous_periodic_optic_definitions(
            definition.autonomous_optics)
    _require_unique_trigger_consumers(acquisition_definitions,
        autonomous_definitions)
    _require_bound_trigger_consumers(acquisition_definitions,
        autonomous_definitions,
        definition.trigger_topology)
    owners, lifecycles, products, sample_providers, rngs, path_slots =
        _prepare_event_acquisition_parts(plant, acquisition_definitions,
            sample_definitions, definition.trigger_topology)
    generator_definitions = EventGeneratorDefinition[]
    _append_event_generator_definitions!(generator_definitions,
        sample_definitions, acquisition_definitions,
        getfield(plant, :command_endpoints),
        definition.trigger_topology)
    scheduler = prepare_event_scheduler(generator_definitions;
        capacity=length(generator_definitions))
    actions = _prepared_event_actions(scheduler)
    optics = _prepare_event_controllable_optics(plant)
    command_endpoints = _prepare_event_command_endpoints(plant, scheduler)
    paths = _prepare_event_paths(plant, sample_definitions, owners,
        scheduler)
    acquisitions = _prepare_event_acquisitions(acquisition_definitions,
        lifecycles, products, sample_providers, rngs, path_slots, scheduler)
    autonomous_optics = _prepare_event_autonomous_optics(
        autonomous_definitions, optics, paths, definition.trigger_topology)
    atmosphere = _require_selection_atmosphere(
        plant_atmosphere(getfield(plant, :definition)))
    atmosphere_rng = _prepared_atmosphere_rng(atmosphere,
        getfield(getfield(plant, :rngs), :atmosphere))
    return PreparedPlantEventLoop(_PlantEventLoopBinding(), atmosphere,
        atmosphere_rng, scheduler, actions, optics, command_endpoints, paths,
        acquisitions, autonomous_optics,
        _prepared_trigger_topology(definition.trigger_topology))
end

mutable struct PlantEventLoopState{T}
    binding::_PlantEventLoopBinding
    scheduler::EventSchedulerState
    command_endpoints::Memory{CommandEndpointState}
    command_applications::Memory{CommandApplicationState}
    command_shadow_transactions::Memory{UInt64}
    controllable_optics::Memory{Any}
    acquisitions::Memory{_AcquisitionEventLifecycleState}
    path_sampled::Memory{Bool}
    product_sequences::Memory{UInt64}
    product_ready_timestamps::Memory{PlantTimestamp}
    command_transaction_sequence::UInt64
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
@inline _event_acquisition_state(
    prepared::PreparedDirectMeasurementAcquisition) =
    DirectMeasurementAcquisitionState(prepared)

@inline _event_trigger_state(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyState()
@inline _event_trigger_state(topology::PreparedTriggerTopology) =
    TriggerTopologyState(topology)

function _event_command_states(prepared::PreparedPlantEventLoop,
    initial_timestamp::PlantTimestamp)
    endpoints = Memory{CommandEndpointState}(undef,
        length(prepared.command_endpoints))
    applications = Memory{CommandApplicationState}(undef,
        length(prepared.command_endpoints))
    @inbounds for index in eachindex(prepared.command_endpoints)
        binding = prepared.command_endpoints[index].binding
        endpoint = binding.endpoint
        endpoint_state = CommandEndpointState(endpoint;
            initial_timestamp)
        endpoints[index] = endpoint_state
        applications[index] = CommandApplicationState(endpoint,
            endpoint_state, binding.initial_command;
            safe_command=binding.safe_command)
    end
    return endpoints, applications
end

function _event_optic_endpoint_initials(
    prepared::PreparedPlantEventLoop,
    optic::PreparedControllableOptic)
    ids = map(optic.endpoint_slots) do slot
        command_endpoint_id(
            prepared.command_endpoints[Int(slot)].binding)
    end
    commands = map(optic.endpoint_slots) do slot
        binding = prepared.command_endpoints[Int(slot)].binding
        _copy_prepared_effective_command(binding.endpoint,
            binding.initial_command, "initial physical command")
    end
    return ids, commands
end

function _event_controllable_optic_states(
    prepared::PreparedPlantEventLoop)
    states = Memory{Any}(undef, length(prepared.optics))
    @inbounds for index in eachindex(prepared.optics)
        optic = prepared.optics[index]
        ids, commands = _event_optic_endpoint_initials(prepared, optic)
        states[index] = prepare_controllable_optic_state(
            optic.implementation, optic.definition, ids, commands)
    end
    return states
end

function _initialize_event_autonomous_optics!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState)
    @inbounds for binding in prepared.autonomous_optics
        optic_state = state.controllable_optics[Int(binding.optic_slot)]
        initialize_autonomous_periodic_optic!(binding.implementation,
            optic_state, binding.phase_reference)
    end
    return nothing
end

function PlantEventLoopState(prepared::PreparedPlantEventLoop)
    scheduler = EventSchedulerState(prepared.scheduler)
    command_endpoints, command_applications = _event_command_states(
        prepared, scheduler_timestamp(scheduler))
    command_shadow_transactions = Memory{UInt64}(undef,
        length(prepared.command_endpoints))
    fill!(command_shadow_transactions, UInt64(0))
    controllable_optics = _event_controllable_optic_states(prepared)
    acquisition_states = Memory{_AcquisitionEventLifecycleState}(undef,
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
    state = PlantEventLoopState(prepared.binding,
        scheduler, command_endpoints, command_applications,
        command_shadow_transactions, controllable_optics,
        acquisition_states, path_sampled,
        product_sequences, product_ready_timestamps, UInt64(0), trigger)
    _initialize_event_autonomous_optics!(prepared, state)
    @inbounds for index in eachindex(prepared.command_endpoints)
        _schedule_event_command_endpoint!(prepared, state, index)
    end
    return state
end

mutable struct PlantEventLoopWorkspace{T}
    binding::_PlantEventLoopBinding
    scheduler::EventSchedulerWorkspace
    command_endpoints::Memory{CommandDispositionWorkspace}
    controllable_optics::Memory{Any}
    command_dispositions::Memory{PlantCommandDisposition}
    command_disposition_count::Int
    transaction_endpoint_slots::Memory{UInt32}
    transaction_admissions::Memory{_CommandTransactionAdmissionPlan}
    transaction_claims::Memory{PlantCommandApplicationClaim}
    transaction_staged::Memory{_StagedCommandApplication}
    transaction_count::Int
    due_paths::Memory{Bool}
    trigger::T
    delivery::Base.RefValue{TriggerDelivery}
end

@inline _event_trigger_workspace(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyWorkspace()
@inline _event_trigger_workspace(topology::PreparedTriggerTopology) =
    TriggerTopologyWorkspace(topology)

function PlantEventLoopWorkspace(prepared::PreparedPlantEventLoop)
    command_endpoints = Memory{CommandDispositionWorkspace}(undef,
        length(prepared.command_endpoints))
    disposition_capacity = 0
    @inbounds for index in eachindex(prepared.command_endpoints)
        endpoint = prepared.command_endpoints[index].binding.endpoint
        command_endpoints[index] = CommandDispositionWorkspace(endpoint)
        disposition_capacity += command_endpoint_capacity(endpoint)
    end
    controllable_optics = Memory{Any}(undef, length(prepared.optics))
    @inbounds for index in eachindex(prepared.optics)
        controllable_optics[index] =
            prepare_controllable_optic_workspace(
                prepared.optics[index].implementation)
    end
    endpoint_count = length(prepared.command_endpoints)
    due_paths = Memory{Bool}(undef, length(prepared.paths))
    fill!(due_paths, false)
    return PlantEventLoopWorkspace(prepared.binding,
        EventSchedulerWorkspace(prepared.scheduler), command_endpoints,
        controllable_optics,
        Memory{PlantCommandDisposition}(undef, disposition_capacity), 0,
        Memory{UInt32}(undef, endpoint_count),
        Memory{_CommandTransactionAdmissionPlan}(undef, endpoint_count),
        Memory{PlantCommandApplicationClaim}(undef, endpoint_count),
        Memory{_StagedCommandApplication}(undef, endpoint_count), 0,
        due_paths,
        _event_trigger_workspace(prepared.trigger_topology),
        Ref{TriggerDelivery}())
end

@inline function _require_plant_event_loop_binding(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState)
    state.binding === prepared.binding || _plant_event_loop_error(
        :foreign_state,
        "plant event-loop state belongs to another prepared loop")
    _require_scheduler_binding(prepared.scheduler, state.scheduler)
    length(state.command_endpoints) == length(prepared.command_endpoints) &&
        length(state.command_applications) ==
            length(prepared.command_endpoints) &&
        length(state.command_shadow_transactions) ==
            length(prepared.command_endpoints) &&
        length(state.controllable_optics) == length(prepared.optics) &&
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
    length(workspace.command_endpoints) ==
        length(prepared.command_endpoints) &&
        length(workspace.controllable_optics) == length(prepared.optics) &&
        length(workspace.transaction_endpoint_slots) ==
            length(prepared.command_endpoints) &&
        length(workspace.transaction_admissions) ==
            length(prepared.command_endpoints) &&
        length(workspace.transaction_claims) ==
            length(prepared.command_endpoints) &&
        length(workspace.transaction_staged) ==
            length(prepared.command_endpoints) ||
        _plant_event_loop_error(:prepared_binding,
            "plant event-loop workspace command capacity changed after preparation")
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

function _event_command_endpoint_slot(prepared::PreparedPlantEventLoop,
    id::CommandEndpointID)
    @inbounds for index in eachindex(prepared.command_endpoints)
        command_endpoint_id(prepared.command_endpoints[index].binding) == id &&
            return index
    end
    _plant_event_loop_error(:unknown_command_endpoint,
        "prepared plant event loop has no command endpoint $id")
end

function effective_command(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    slot = _event_command_endpoint_slot(prepared,
        _as_command_endpoint_id(id))
    return effective_command(state.command_applications[slot])
end

function _event_autonomous_optic_slot(
    prepared::PreparedPlantEventLoop, id::ControllableOpticID)
    @inbounds for index in eachindex(prepared.autonomous_optics)
        prepared.autonomous_optics[index].id == id && return index
    end
    _plant_event_loop_error(:unknown_autonomous_optic,
        "prepared plant event loop has no autonomous optic $id")
end

@inline function _event_autonomous_optic_parts(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState, id)
    slot = _event_autonomous_optic_slot(prepared,
        _as_controllable_optic_id(id))
    binding = @inbounds prepared.autonomous_optics[slot]
    optic_state = @inbounds state.controllable_optics[
        Int(binding.optic_slot)]
    return binding, optic_state
end

function autonomous_waveform_phase(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id, timestamp::PlantTimestamp)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_phase(binding.implementation, optic_state,
        timestamp)
end

function autonomous_waveform_phase(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    return autonomous_waveform_phase(prepared, state, id,
        scheduler_timestamp(state.scheduler))
end

function autonomous_waveform_offset(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id, timestamp::PlantTimestamp)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_offset(binding.implementation, optic_state,
        timestamp)
end

function autonomous_waveform_offset(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    return autonomous_waveform_offset(prepared, state, id,
        scheduler_timestamp(state.scheduler))
end

function autonomous_waveform_reference_timestamp(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_reference_timestamp(
        binding.implementation, optic_state)
end

function autonomous_waveform_reference_sequence(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_reference_sequence(
        binding.implementation, optic_state)
end

function autonomous_waveform_reference_count(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_reference_count(
        binding.implementation, optic_state)
end

function autonomous_waveform_enabled(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_enabled(binding.implementation, optic_state)
end

function autonomous_waveform_radius(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, id)
    _require_plant_event_loop_binding(prepared, state)
    binding, optic_state = _event_autonomous_optic_parts(prepared, state, id)
    return autonomous_waveform_radius(binding.implementation, optic_state)
end

@inline _event_command_endpoint(
    prepared::PreparedPlantEventLoop, slot::Integer) =
    @inbounds prepared.command_endpoints[Int(slot)]
@inline _event_command_endpoint_state(
    state::PlantEventLoopState, slot::Integer) =
    @inbounds state.command_endpoints[Int(slot)]
@inline _event_command_application_state(
    state::PlantEventLoopState, slot::Integer) =
    @inbounds state.command_applications[Int(slot)]
@inline _event_command_workspace(
    workspace::PlantEventLoopWorkspace, slot::Integer) =
    @inbounds workspace.command_endpoints[Int(slot)]

@inline command_disposition_count(
    workspace::PlantEventLoopWorkspace) =
    workspace.command_disposition_count

function command_disposition(workspace::PlantEventLoopWorkspace,
    index::Integer)
    1 <= index <= workspace.command_disposition_count ||
        _plant_event_loop_error(:invalid_disposition_index,
            "command disposition index must be within the current event-loop records")
    return @inbounds workspace.command_dispositions[Int(index)]
end

@inline command_disposition(::PlantEventLoopWorkspace, ::Bool) =
    _plant_event_loop_error(:invalid_disposition_index,
        "command disposition index must be an integer count, not Bool")

function clear_command_dispositions!(
    workspace::PlantEventLoopWorkspace)
    workspace.command_disposition_count = 0
    return workspace
end

@inline function _require_empty_event_command_dispositions(
    workspace::PlantEventLoopWorkspace)
    iszero(workspace.command_disposition_count) ||
        _plant_event_loop_error(:unconsumed_command_dispositions,
            "clear event-loop command dispositions before the next " *
            "command admission")
    return nothing
end

function _append_event_command_dispositions!(
    workspace::PlantEventLoopWorkspace,
    endpoint_workspace::CommandDispositionWorkspace)
    count = command_disposition_count(endpoint_workspace)
    first = workspace.command_disposition_count + 1
    last = workspace.command_disposition_count + count
    last <= length(workspace.command_dispositions) ||
        _plant_event_loop_error(:command_disposition_capacity,
            "plant event-loop command disposition capacity was exceeded")
    @inbounds for source in 1:count
        workspace.command_dispositions[first + source - 1] =
            command_disposition(endpoint_workspace, source)
    end
    workspace.command_disposition_count = last
    clear_command_dispositions!(endpoint_workspace)
    return count
end

function _next_event_command_timestamp(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, slot::Integer)
    binding = _event_command_endpoint(prepared, slot).binding
    endpoint = binding.endpoint
    endpoint_state = _event_command_endpoint_state(state, slot)
    application_state = _event_command_application_state(state, slot)
    key = next_command_order_key(endpoint, endpoint_state)
    command_timestamp = key === nothing ? nothing :
        command_scheduled_timestamp(key)
    silence_timestamp = next_command_silence_timestamp(endpoint,
        endpoint_state, application_state)
    command_timestamp === nothing && return silence_timestamp
    silence_timestamp === nothing && return command_timestamp
    return min(command_timestamp, silence_timestamp)
end

function _schedule_event_command_endpoint!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    slot::Integer)
    event_endpoint = _event_command_endpoint(prepared, slot)
    handle = event_endpoint.handle
    desired = _next_event_command_timestamp(prepared, state, slot)
    cursor = state.scheduler.cursors[Int(handle.slot)]
    if desired === nothing
        cursor.status == _ScheduledEventGenerator &&
            deactivate_event_generator!(prepared.scheduler,
                state.scheduler, handle)
        return nothing
    end
    if cursor.status == _InactiveEventGenerator
        activate_event_generator!(prepared.scheduler, state.scheduler,
            handle, desired)
        return desired
    end
    cursor.status == _ScheduledEventGenerator ||
        _plant_event_loop_error(:command_generator_state,
            "command endpoint generator cannot be changed while claimed")
    cursor.next_timestamp == desired && return desired
    deactivate_event_generator!(prepared.scheduler, state.scheduler, handle)
    activate_event_generator!(prepared.scheduler, state.scheduler, handle,
        desired)
    return desired
end

function _resolve_event_command_claim!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    claim::EventClaim,
    slot::Integer)
    desired = _next_event_command_timestamp(prepared, state, slot)
    if desired === nothing
        deactivate_event_generator!(prepared.scheduler, state.scheduler,
            claim)
    else
        reschedule_event!(prepared.scheduler, state.scheduler, claim,
            desired)
    end
    return desired
end

@inline function _require_routed_command_admission_timestamp(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp)
    scheduler_state = state.scheduler
    current = scheduler_timestamp(scheduler_state)
    timestamp < current &&
        _plant_event_loop_error(:command_admission_time_regression,
            "command admission timestamp precedes the current plant-event " *
            "timestamp")
    timestamp == current && scheduler_state.has_last_key &&
        _plant_event_loop_error(:command_admission_time_elapsed,
            "command admission timestamp has already been processed by the " *
            "plant event loop")
    count = scan_due_events!(workspace.scheduler, prepared.scheduler,
        scheduler_state)
    if !iszero(count)
        next_due = workspace.scheduler.due_timestamp
        next_due < timestamp &&
            _plant_event_loop_error(:command_admission_overtakes_event,
                "command admission timestamp follows the next unprocessed " *
                "plant event at $next_due")
    end
    return nothing
end

"""
Route one command into its exact event-loop-owned endpoint and arm the reserved
command-phase generator. Any admission-time terminal dispositions are copied
into the event-loop's bounded disposition ledger.
"""
function admit_plant_command!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    command::PlantCommand, timestamp::PlantTimestamp)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    _require_empty_event_command_dispositions(workspace)
    _require_routed_command_admission_timestamp(prepared, state, workspace,
        timestamp)
    slot = _event_command_endpoint_slot(prepared,
        command_endpoint_id(command))
    event_endpoint = _event_command_endpoint(prepared, slot)
    endpoint_state = _event_command_endpoint_state(state, slot)
    endpoint_workspace = _event_command_workspace(workspace, slot)
    admission = try
        admit_plant_command!(endpoint_workspace,
            event_endpoint.binding.endpoint,
            endpoint_state, command, timestamp)
    catch
        _append_event_command_dispositions!(workspace, endpoint_workspace)
        _schedule_event_command_endpoint!(prepared, state, slot)
        rethrow()
    end
    _append_event_command_dispositions!(workspace, endpoint_workspace)
    _schedule_event_command_endpoint!(prepared, state, slot)
    return admission
end

@inline function _next_command_transaction_sequence(
    state::PlantEventLoopState)
    state.command_transaction_sequence != typemax(UInt64) ||
        _command_admission_error(:transaction, :transaction_overflow,
            "plant command transaction identity exceeds UInt64 range")
    return state.command_transaction_sequence + UInt64(1)
end

function _transaction_command_for_endpoint(
    transaction::PlantCommandTransaction, id::CommandEndpointID)
    for command in transaction.commands
        command_endpoint_id(command) == id && return command
    end
    return nothing
end

function _prepare_transaction_member_slots!(
    prepared::PreparedPlantEventLoop,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction)
    count = 0
    @inbounds for endpoint_slot in eachindex(prepared.command_endpoints)
        binding = prepared.command_endpoints[endpoint_slot].binding
        command = _transaction_command_for_endpoint(transaction,
            command_endpoint_id(binding))
        command === nothing && continue
        for prior in 1:count
            prior_binding = prepared.command_endpoints[
                Int(workspace.transaction_endpoint_slots[prior])].binding
            prior_binding.optic_slot == binding.optic_slot &&
                _command_admission_error(:transaction,
                    :duplicate_physical_optic,
                    "atomic multi-optic transaction contains more than one " *
                    "endpoint owned by controllable optic " *
                    "$(controllable_optic_id(prepared.optics[
                        Int(binding.optic_slot)].definition))")
        end
        count += 1
        workspace.transaction_endpoint_slots[count] = UInt32(endpoint_slot)
    end
    count == length(transaction.commands) ||
        _command_admission_error(:transaction, :unknown_endpoint,
            "atomic transaction references a command endpoint that is not " *
            "owned by this prepared plant event loop")
    workspace.transaction_count = count
    return count
end

function _require_no_equal_time_command(
    state::CommandEndpointState, timestamp::PlantTimestamp)
    @inbounds for index in 1:state.pending_count
        metadata = state.slots[Int(state.calendar[index])]
        metadata.scheduled_timestamp == timestamp &&
            _command_admission_error(:transaction,
                :equal_time_endpoint_conflict,
                "an atomic transaction member endpoint already has a " *
                "command scheduled at $timestamp")
    end
    return nothing
end

function _preflight_command_transaction_member!(
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    endpoint_workspace::CommandDispositionWorkspace,
    command::PlantCommand,
    timestamp::PlantTimestamp)
    _require_command_endpoint_binding(endpoint, endpoint_state)
    _require_command_endpoint_binding(endpoint, endpoint_workspace)
    _require_operational_command_endpoint(endpoint_state)
    _require_empty_command_dispositions(endpoint_workspace)
    _require_idle_command_endpoint(endpoint_state)
    _require_forward_command_timestamp(endpoint_state, timestamp)
    command_effective_time_policy(command_schema(endpoint)).supersession ==
        PreservePendingCommands || _command_admission_error(:transaction,
        :unsupported_supersession,
        "atomic transaction endpoints must preserve pending commands")
    validate_plant_command(endpoint, command)
    sequence_class = _classify_command_sequence(endpoint, endpoint_state,
        command.sequence)
    sequence_action = _command_sequence_action(
        command_sequence_policy(command_schema(endpoint)), sequence_class)
    if sequence_action != AcceptSequence
        reason = _command_sequence_reason(sequence_class)
        message = "atomic transaction member sequence policy " *
            (sequence_action == FailOnSequence ?
                "requires structural failure" : "rejected the command")
        error = PlantCommandError(:transaction, reason.name, message)
        sequence_action == FailOnSequence &&
            throw(_CommandTransactionPolicyFailure(error))
        throw(error)
    end

    requested = command.requested_effective_timestamp
    scheduled = requested
    policy = command_effective_time_policy(command_schema(endpoint))
    if timestamp < requested
        policy.future == AllowFutureCommand ||
            _command_admission_error(:transaction, :future_command,
                "atomic transaction member rejects future commands")
    elseif requested < timestamp
        if policy.late == FailOnLateCommand
            throw(_CommandTransactionPolicyFailure(PlantCommandError(
                :transaction, :late_command,
                "atomic transaction member late-command policy requires " *
                "structural failure")))
        end
        policy.late == ApplyLateCommandNow ||
            _command_admission_error(:transaction, :late_command,
                "atomic transaction member rejects late commands")
        scheduled = timestamp
    end
    endpoint_state.active_count < endpoint.capacity ||
        _command_admission_error(:transaction, :calendar_capacity,
            "atomic transaction member endpoint is at calendar capacity")
    _require_no_equal_time_command(endpoint_state, scheduled)

    slot = _find_free_command_slot(endpoint_state)
    metadata = endpoint_state.slots[Int(slot)]
    generation = _next_command_counter(metadata.generation,
        :slot_generation_overflow, "command payload-slot generation")
    _stage_command_payload!(endpoint_state.payloads, endpoint.schema,
        command.payload)
    return _CommandTransactionAdmissionPlan(
        _next_command_presentation(endpoint_state), sequence_class,
        scheduled, slot, generation)
end

function _commit_command_transaction_member!(
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    command::PlantCommand,
    timestamp::PlantTimestamp,
    transaction::UInt64,
    member_count::UInt32,
    plan::_CommandTransactionAdmissionPlan)
    slot = Int(plan.slot)
    _commit_staged_command_payload!(endpoint_state.payloads, slot, nothing)
    endpoint_state.slots[slot] = _CommandSlotMetadata(
        plan.presentation.value,
        command.sequence.value,
        command.requested_effective_timestamp,
        plan.scheduled_timestamp,
        timestamp,
        transaction,
        member_count,
        plan.generation,
        _PendingCommandSlot,
    )
    _insert_pending_command_slot!(endpoint, endpoint_state, plan.slot)
    endpoint_state.active_count += 1
    _record_accepted_command_sequence!(endpoint, endpoint_state,
        command.sequence)
    endpoint_state.presentation_sequence = plan.presentation.value
    endpoint_state.current_timestamp = timestamp
    endpoint_state.last_admission_timestamp = timestamp
    endpoint_state.has_admission = true
    return nothing
end

function _terminate_event_command_transaction_admission!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    kind::CommandTerminalKind,
    reason::CommandDispositionReason)
    count = workspace.transaction_count
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        endpoint = binding.endpoint
        endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
        endpoint_workspace = _event_command_workspace(workspace,
            endpoint_slot)
        command = _transaction_command_for_endpoint(transaction,
            command_endpoint_id(binding))
        _finish_terminal_admission!(endpoint_workspace, endpoint,
            endpoint_state, command, timestamp,
            _next_command_presentation(endpoint_state), kind,
            reason, nothing)
        _append_event_command_dispositions!(workspace, endpoint_workspace)
        _schedule_event_command_endpoint!(prepared, state, endpoint_slot)
    end
    workspace.transaction_count = 0
    return PlantCommandTransactionAdmission(UInt64(0),
        CommandTerminatedOnAdmission, UInt32(count), nothing)
end

@noinline function _handle_transaction_preflight_failure!(
    ::PreparedPlantEventLoop,
    ::PlantEventLoopState,
    ::PlantEventLoopWorkspace,
    ::PlantCommandTransaction,
    ::PlantTimestamp,
    error::InterruptException)
    throw(error)
end

@noinline function _transaction_preflight_reason(error)
    return CommandDispositionReason(:transaction_preflight_failure)
end

@noinline function _handle_transaction_preflight_failure!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    failure::_CommandTransactionPolicyFailure)
    error = failure.error
    reason = CommandDispositionReason(Symbol(:atomic_transaction_aborted_,
        error.reason))
    _terminate_event_command_transaction_admission!(prepared, state,
        workspace, transaction, timestamp, FailedCommand, reason)
    throw(error)
end

@noinline function _handle_transaction_preflight_failure!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    error::PlantCommandError)
    reason = CommandDispositionReason(Symbol(:atomic_transaction_aborted_,
        error.reason))
    kind = _validation_terminal_kind(error)
    admission = _terminate_event_command_transaction_admission!(prepared,
        state, workspace, transaction, timestamp, kind, reason)
    kind == FailedCommand && throw(error)
    return admission
end

@noinline function _handle_transaction_preflight_failure!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    error)
    reason = _transaction_preflight_reason(error)
    _terminate_event_command_transaction_admission!(prepared, state,
        workspace, transaction, timestamp, FailedCommand, reason)
    _command_admission_error(:transaction, reason.name,
        "atomic transaction preflight failed unexpectedly " *
        "($(typeof(error)))")
end

"""
Atomically admit a bounded multi-optic command transaction. Every member is
validated and staged before any endpoint calendar is mutated. A normal member
rejection terminates every member with the same explicit transaction-abort
reason; successful admission assigns one transaction identity carried by every
member command.
"""
function admit_plant_command_transaction!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    _require_empty_event_command_dispositions(workspace)
    _require_routed_command_admission_timestamp(prepared, state, workspace,
        timestamp)
    count = _prepare_transaction_member_slots!(prepared, workspace,
        transaction)
    scheduled = zero(PlantTimestamp)
    has_scheduled = false
    failure = nothing
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        command = _transaction_command_for_endpoint(transaction,
            command_endpoint_id(binding))
        plan = try
            _preflight_command_transaction_member!(binding.endpoint,
                _event_command_endpoint_state(state, endpoint_slot),
                _event_command_workspace(workspace, endpoint_slot),
                command, timestamp)
        catch error
            failure = error
            break
        end
        if has_scheduled && plan.scheduled_timestamp != scheduled
            failure = PlantCommandError(:transaction,
                :scheduled_timestamp_mismatch,
                "atomic transaction members resolved to different " *
                "scheduled plant timestamps")
            break
        end
        scheduled = plan.scheduled_timestamp
        has_scheduled = true
        workspace.transaction_admissions[index] = plan
    end
    if failure !== nothing
        return _handle_transaction_preflight_failure!(prepared, state,
            workspace, transaction, timestamp, failure)
    end

    transaction_id = _next_command_transaction_sequence(state)
    member_count = UInt32(count)
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        command = _transaction_command_for_endpoint(transaction,
            command_endpoint_id(binding))
        _commit_command_transaction_member!(binding.endpoint,
            _event_command_endpoint_state(state, endpoint_slot),
            command, timestamp, transaction_id, member_count,
            workspace.transaction_admissions[index])
    end
    state.command_transaction_sequence = transaction_id
    @inbounds for index in 1:count
        _schedule_event_command_endpoint!(prepared, state,
            workspace.transaction_endpoint_slots[index])
    end
    workspace.transaction_count = 0
    status = scheduled <= timestamp ?
        CommandAdmittedReady : CommandAdmittedPending
    return PlantCommandTransactionAdmission(transaction_id, status,
        member_count, scheduled)
end

@inline function _event_controllable_optic_parts(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer)
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    optic_slot = Int(binding.optic_slot)
    optic = @inbounds prepared.optics[optic_slot]
    optic_state = @inbounds state.controllable_optics[optic_slot]
    optic_workspace = @inbounds workspace.controllable_optics[optic_slot]
    return binding, optic, optic_state, optic_workspace
end

@inline function _stage_event_controllable_optic_command!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer)
    binding, optic, optic_state, optic_workspace =
        _event_controllable_optic_parts(prepared, state, workspace,
            endpoint_slot)
    application_state =
        _event_command_application_state(state, endpoint_slot)
    stage_controllable_optic_command!(optic.implementation, optic_state,
        optic_workspace, command_endpoint_id(binding),
        _staged_effective_command(application_state),
        state.scheduler.current_timestamp)
    return nothing
end

@inline function _commit_event_controllable_optic_command!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer)
    binding, optic, optic_state, optic_workspace =
        _event_controllable_optic_parts(prepared, state, workspace,
            endpoint_slot)
    commit_controllable_optic_command!(optic.implementation, optic_state,
        optic_workspace, command_endpoint_id(binding),
        state.scheduler.current_timestamp)
    return nothing
end

@noinline function _fail_composed_command_staging!(
    endpoint_workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error::InterruptException)
    throw(error)
end

@noinline function _fail_composed_command_staging!(
    endpoint_workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error::PlantCommandError)
    _finish_command_application!(endpoint_workspace, endpoint,
        endpoint_state, claim, FailedCommand,
        CommandDispositionReason(error.reason))
    throw(error)
end

@noinline function _fail_composed_command_staging!(
    endpoint_workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error)
    _finish_command_application!(endpoint_workspace, endpoint,
        endpoint_state, claim, FailedCommand,
        CommandDispositionReason(:physical_application_failure))
    _command_admission_error(:physical_application,
        :physical_application_failure,
        "failed to stage a physical controllable-optic command " *
        "($(typeof(error)))")
end

function _apply_event_command_claim_parts!(
    binding::_PreparedPlantCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    endpoint_workspace::CommandDispositionWorkspace,
    optic::PreparedControllableOptic,
    optic_state,
    optic_workspace,
    claim::PlantCommandApplicationClaim,
    timestamp::PlantTimestamp)
    endpoint = binding.endpoint
    staged = try
        result = _stage_claimed_plant_command!(endpoint, endpoint_state,
            application_state, claim)
        if result.decision == _AcceptCommandCandidate
            stage_controllable_optic_command!(optic.implementation,
                optic_state, optic_workspace, command_endpoint_id(binding),
                _staged_effective_command(application_state), timestamp)
        end
        result
    catch error
        _fail_composed_command_staging!(endpoint_workspace, endpoint,
            endpoint_state, claim, error)
    end
    decision = staged.decision
    reason = staged.reason
    if decision == _AcceptCommandCandidate
        commit_controllable_optic_command!(optic.implementation, optic_state,
            optic_workspace, command_endpoint_id(binding), timestamp)
        _commit_staged_application!(application_state, endpoint_state)
        _finish_command_application!(endpoint_workspace, endpoint,
            endpoint_state, claim, AppliedCommand, reason)
    else
        kind = decision == _FailCommandCandidate ?
            FailedCommand : RejectedCommand
        _finish_command_application!(endpoint_workspace, endpoint,
            endpoint_state, claim, kind, reason)
    end
    return staged
end

function _apply_event_command_claim!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer,
    claim::PlantCommandApplicationClaim)
    binding, optic, optic_state, optic_workspace =
        _event_controllable_optic_parts(prepared, state, workspace,
            endpoint_slot)
    endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
    application_state =
        _event_command_application_state(state, endpoint_slot)
    endpoint_workspace = _event_command_workspace(workspace, endpoint_slot)
    staged = try
        _apply_event_command_claim_parts!(binding, endpoint_state,
            application_state, endpoint_workspace, optic, optic_state,
            optic_workspace, claim, state.scheduler.current_timestamp)
    catch
        _append_event_command_dispositions!(workspace, endpoint_workspace)
        rethrow()
    end
    _append_event_command_dispositions!(workspace, endpoint_workspace)
    staged.decision == _FailCommandCandidate && _command_admission_error(
        :application, staged.reason.name,
        "application-stage command policy requires structural failure")
    return nothing
end

function _apply_event_command_silence!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer,
    timestamp::PlantTimestamp)
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    endpoint = binding.endpoint
    endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
    application_state =
        _event_command_application_state(state, endpoint_slot)
    endpoint_workspace = _event_command_workspace(workspace, endpoint_slot)
    policy = command_silence_policy(command_schema(endpoint))
    if policy.action != ApplySafeCommand
        transition = apply_command_silence_transition!(endpoint_workspace,
            endpoint, endpoint_state, application_state, timestamp)
        _append_event_command_dispositions!(workspace, endpoint_workspace)
        return transition
    end

    expected = next_command_silence_timestamp(endpoint, endpoint_state,
        application_state)
    expected == timestamp || _command_admission_error(
        :silence, :unexpected_silence_timestamp,
        "command silence transition expected at $expected; got $timestamp")
    key = next_command_order_key(endpoint, endpoint_state)
    key !== nothing && command_scheduled_timestamp(key) <= timestamp &&
        _command_admission_error(:silence, :commands_due,
            "resolve application-ready commands before an equal-time " *
            "command-silence transition")
    origin = _command_silence_origin_timestamp(policy.age_origin,
        endpoint_state, application_state)
    try
        _stage_safe_command!(application_state.values)
        _stage_event_controllable_optic_command!(prepared, state,
            workspace, endpoint_slot)
    catch error
        _handle_safe_command_staging_error(error)
    end
    _commit_event_controllable_optic_command!(prepared, state, workspace,
        endpoint_slot)
    _commit_application_candidate!(application_state.values)
    endpoint_state.current_timestamp = timestamp
    application_state.last_silence_origin_timestamp = origin
    application_state.has_silence_transition = true
    return PlantCommandSilenceTransition(command_endpoint_id(endpoint),
        policy.action, policy.age_origin, origin, expected, timestamp)
end

function _claim_event_command_transaction!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    first_endpoint_slot::Integer,
    first_claim::PlantCommandApplicationClaim,
    timestamp::PlantTimestamp)
    transaction = _command_transaction_id(first_claim)
    expected = Int(_command_transaction_member_count(first_claim))
    workspace.transaction_count = 1
    workspace.transaction_endpoint_slots[1] = UInt32(first_endpoint_slot)
    workspace.transaction_claims[1] = first_claim
    @inbounds for endpoint_slot in eachindex(prepared.command_endpoints)
        endpoint_slot == first_endpoint_slot && continue
        binding = prepared.command_endpoints[endpoint_slot].binding
        endpoint_state = state.command_endpoints[endpoint_slot]
        member_transaction, member_count =
            _next_command_transaction_metadata(binding.endpoint,
                endpoint_state)
        member_transaction == transaction || continue
        Int(member_count) == expected ||
            _plant_event_loop_error(:command_transaction_invariant,
                "atomic transaction members disagree about member count")
        key = next_command_order_key(binding.endpoint, endpoint_state)
        key !== nothing &&
            command_scheduled_timestamp(key) == timestamp ||
            _plant_event_loop_error(:command_transaction_invariant,
                "atomic transaction member is not ready at its common timestamp")
        member_claim = claim_next_application_ready_command!(
            binding.endpoint, endpoint_state, timestamp)
        member_claim === nothing &&
            _plant_event_loop_error(:command_transaction_invariant,
                "atomic transaction member did not produce a command claim")
        workspace.transaction_count += 1
        index = workspace.transaction_count
        workspace.transaction_endpoint_slots[index] = UInt32(endpoint_slot)
        workspace.transaction_claims[index] = member_claim
        iszero(state.command_shadow_transactions[endpoint_slot]) ||
            _plant_event_loop_error(:command_transaction_invariant,
                "command endpoint already owns an unconsumed transaction event")
        state.command_shadow_transactions[endpoint_slot] = transaction
    end
    workspace.transaction_count == expected ||
        _plant_event_loop_error(:command_transaction_invariant,
            "atomic transaction did not yield every declared member")
    return expected
end

@inline _transaction_application_failure_reason(error::PlantCommandError) =
    CommandDispositionReason(error.reason)
@inline _transaction_application_failure_reason(error::InterruptException) =
    throw(error)

@noinline function _transaction_application_failure_reason(error)
    return CommandDispositionReason(:physical_application_failure)
end

@noinline function _fail_event_command_transaction!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    reason::CommandDispositionReason,
    error)
    count = workspace.transaction_count
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        endpoint_workspace = _event_command_workspace(workspace,
            endpoint_slot)
        _finish_command_application!(endpoint_workspace, binding.endpoint,
            _event_command_endpoint_state(state, endpoint_slot),
            workspace.transaction_claims[index], FailedCommand, reason)
        _append_event_command_dispositions!(workspace, endpoint_workspace)
    end
    workspace.transaction_count = 0
    _command_admission_error(:transaction, reason.name,
        "atomic command transaction application failed ($(typeof(error)))")
end

function _finish_rejected_event_command_transaction!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace)
    count = workspace.transaction_count
    failure_reason = nothing
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        staged = workspace.transaction_staged[index]
        endpoint_workspace = _event_command_workspace(workspace,
            endpoint_slot)
        if staged.decision == _AcceptCommandCandidate
            kind = RejectedCommand
            reason = CommandDispositionReason(:atomic_transaction_aborted)
        else
            kind = staged.decision == _FailCommandCandidate ?
                FailedCommand : RejectedCommand
            reason = staged.reason
            kind == FailedCommand && (failure_reason = reason)
        end
        _finish_command_application!(endpoint_workspace, binding.endpoint,
            _event_command_endpoint_state(state, endpoint_slot),
            workspace.transaction_claims[index], kind, reason)
        _append_event_command_dispositions!(workspace, endpoint_workspace)
    end
    workspace.transaction_count = 0
    failure_reason === nothing || _command_admission_error(
        :transaction, failure_reason.name,
        "atomic transaction member application policy requires structural failure")
    return nothing
end

function _apply_event_command_transaction!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    first_endpoint_slot::Integer,
    first_claim::PlantCommandApplicationClaim,
    timestamp::PlantTimestamp)
    count = _claim_event_command_transaction!(prepared, state, workspace,
        first_endpoint_slot, first_claim, timestamp)
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        staged = try
            _stage_claimed_plant_command!(binding.endpoint,
                _event_command_endpoint_state(state, endpoint_slot),
                _event_command_application_state(state, endpoint_slot),
                workspace.transaction_claims[index])
        catch error
            reason = _transaction_application_failure_reason(error)
            _fail_event_command_transaction!(prepared, state, workspace,
                reason, error)
        end
        workspace.transaction_staged[index] = staged
    end
    @inbounds for index in 1:count
        if workspace.transaction_staged[index].decision !=
                _AcceptCommandCandidate
            return _finish_rejected_event_command_transaction!(prepared,
                state, workspace)
        end
    end
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        try
            _stage_event_controllable_optic_command!(prepared, state,
                workspace, endpoint_slot)
        catch error
            reason = _transaction_application_failure_reason(error)
            _fail_event_command_transaction!(prepared, state, workspace,
                reason, error)
        end
    end
    @inbounds for index in 1:count
        _commit_event_controllable_optic_command!(prepared, state,
            workspace, workspace.transaction_endpoint_slots[index])
    end
    @inbounds for index in 1:count
        endpoint_slot = Int(workspace.transaction_endpoint_slots[index])
        binding = _event_command_endpoint(prepared, endpoint_slot).binding
        endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
        application_state =
            _event_command_application_state(state, endpoint_slot)
        endpoint_workspace = _event_command_workspace(workspace,
            endpoint_slot)
        staged = workspace.transaction_staged[index]
        _commit_staged_application!(application_state, endpoint_state)
        _finish_command_application!(endpoint_workspace, binding.endpoint,
            endpoint_state, workspace.transaction_claims[index],
            AppliedCommand, staged.reason)
        _append_event_command_dispositions!(workspace, endpoint_workspace)
    end
    workspace.transaction_count = 0
    return nothing
end

function _process_command_endpoint!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::EventClaim,
    action::_PlantEventAction)
    endpoint_slot = Int(action.owner_slot)
    timestamp = claimed_event_key(claim).timestamp
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    endpoint = binding.endpoint
    endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
    shadow = state.command_shadow_transactions[endpoint_slot]
    if !iszero(shadow)
        state.command_shadow_transactions[endpoint_slot] = UInt64(0)
        _resolve_event_command_claim!(prepared, state, claim, endpoint_slot)
        return nothing
    end
    key = next_command_order_key(endpoint, endpoint_state)
    if key !== nothing && command_scheduled_timestamp(key) <= timestamp
        command_claim = claim_next_application_ready_command!(endpoint,
            endpoint_state, timestamp)
        command_claim === nothing &&
            _plant_event_loop_error(:missing_command_claim,
                "due command endpoint did not produce an application claim")
        if iszero(_command_transaction_id(command_claim))
            _apply_event_command_claim!(prepared, state, workspace,
                endpoint_slot, command_claim)
        else
            _apply_event_command_transaction!(prepared, state, workspace,
                endpoint_slot, command_claim, timestamp)
        end
    else
        silence = next_command_silence_timestamp(endpoint, endpoint_state,
            _event_command_application_state(state, endpoint_slot))
        silence == timestamp ||
            _plant_event_loop_error(:stale_command_event,
                "command generator has no due command or silence transition")
        _apply_event_command_silence!(prepared, state, workspace,
            endpoint_slot, timestamp)
    end
    _resolve_event_command_claim!(prepared, state, claim, endpoint_slot)
    return nothing
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

@inline function _first_acquisition_boundary_timestamp(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState)
    next_read = next_nondestructive_read_timestamp(prepared, state)
    next_read === nothing && return state.exposure_close
    return min(next_read, state.exposure_close)
end

@inline function _initial_acquisition_boundary_timestamp(
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

@inline function _initial_acquisition_boundary_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    timestamp::PlantTimestamp)
    return timestamp + prepared.definition.exposure_duration
end

@inline function _initial_acquisition_boundary_timestamp(
    prepared::PreparedFrameTransferAcquisition,
    timestamp::PlantTimestamp)
    return timestamp + prepared.definition.exposure_duration
end

@inline function _first_acquisition_boundary_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    return next_rolling_band_close_timestamp(prepared, state)
end

@inline function _first_acquisition_boundary_timestamp(
    ::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState)
    return state.exposure_close
end

@inline _first_acquisition_band_open_timestamp(
    ::PreparedGlobalShutterAcquisition,
    ::GlobalShutterAcquisitionState) = nothing
@inline _first_acquisition_band_open_timestamp(
    ::PreparedFrameTransferAcquisition,
    ::FrameTransferAcquisitionState) = nothing
@inline function _first_acquisition_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState)
    return next_rolling_band_open_timestamp(prepared, state)
end

@inline _initial_acquisition_band_open_timestamp(
    ::PreparedGlobalShutterAcquisition, ::PlantTimestamp) = nothing
@inline _initial_acquisition_band_open_timestamp(
    ::PreparedFrameTransferAcquisition, ::PlantTimestamp) = nothing
@inline function _initial_acquisition_band_open_timestamp(
    prepared::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_RollingExposureEventMode},
    timestamp::PlantTimestamp)
    prepared.band_count == 1 && return nothing
    return timestamp + prepared.line_duration
end
@inline _initial_acquisition_band_open_timestamp(
    ::PreparedRollingShutterAcquisition{
        <:Any,<:Any,<:Any,<:_GlobalResetEventMode},
    ::PlantTimestamp) = nothing

@inline function _take_initial_acquisition_snapshot!(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    next_read = next_nondestructive_read_timestamp(prepared, state)
    next_read == timestamp || return nothing
    take_nondestructive_read!(prepared, state, timestamp, rng)
    return nothing
end

@inline _take_initial_acquisition_snapshot!(
    ::PreparedRollingShutterAcquisition,
    ::RollingShutterAcquisitionState, ::PlantTimestamp,
    ::AbstractRNG) = nothing
@inline _take_initial_acquisition_snapshot!(
    ::PreparedFrameTransferAcquisition,
    ::FrameTransferAcquisitionState, ::PlantTimestamp,
    ::AbstractRNG) = nothing

@inline function _initial_acquisition_boundary_timestamp(
    prepared::PreparedDirectMeasurementAcquisition,
    timestamp::PlantTimestamp)
    return timestamp + prepared.definition.exposure_duration
end

@inline function _first_acquisition_boundary_timestamp(
    ::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState)
    return state.exposure_close
end

@inline _first_acquisition_band_open_timestamp(
    ::PreparedDirectMeasurementAcquisition,
    ::DirectMeasurementAcquisitionState) = nothing

@inline _initial_acquisition_band_open_timestamp(
    ::PreparedDirectMeasurementAcquisition,
    ::PlantTimestamp) = nothing

@inline _take_initial_acquisition_snapshot!(
    ::PreparedDirectMeasurementAcquisition,
    ::DirectMeasurementAcquisitionState, ::PlantTimestamp,
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

    boundary_timestamp = _initial_acquisition_boundary_timestamp(
        acquisition.lifecycle, timestamp)
    _require_inactive_event_generator(prepared, state,
        acquisition.boundary_handle, boundary_timestamp)
    band_open_timestamp = _initial_acquisition_band_open_timestamp(
        acquisition.lifecycle, timestamp)
    band_open_timestamp === nothing || _require_inactive_event_generator(
        prepared, state, acquisition.band_open_handle, band_open_timestamp)

    begin_exposure!(acquisition.lifecycle, acquisition_state, timestamp)
    _take_initial_acquisition_snapshot!(acquisition.lifecycle,
        acquisition_state, timestamp, acquisition.rng)
    boundary_timestamp == _first_acquisition_boundary_timestamp(
        acquisition.lifecycle, acquisition_state) ||
        _plant_event_loop_error(:prepared_binding,
            "acquisition boundary changed after acquisition start")
    band_open_timestamp == _first_acquisition_band_open_timestamp(
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

@inline function _integrate_event_acquisition_to!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp,
    ::AbstractRNG)
    state.status == DirectMeasurementExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "direct-measurement integration target precedes acquisition progress")
    state.integrated_through == timestamp && return nothing
    return accumulate_direct_measurement_interval!(prepared, state,
        state.integrated_through, timestamp)
end

@enum _AcquisitionBoundaryDisposition::UInt8 begin
    _RescheduleAcquisitionBoundary = 0x01
    _ScheduleAcquisitionReadout = 0x02
end

struct _AcquisitionBoundaryResult
    disposition::_AcquisitionBoundaryDisposition
    timestamp::PlantTimestamp
end

function _process_acquisition_lifecycle_boundary!(
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
        return _AcquisitionBoundaryResult(_RescheduleAcquisitionBoundary,
            following)
    end
    close_exposure!(prepared, state, timestamp)
    return _AcquisitionBoundaryResult(_ScheduleAcquisitionReadout,
        state.readout_complete)
end

function _process_acquisition_lifecycle_boundary!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
    close_next_rolling_band!(prepared, state, timestamp)
    following = next_rolling_band_close_timestamp(prepared, state)
    following === nothing && return _AcquisitionBoundaryResult(
        _ScheduleAcquisitionReadout, state.readout_complete)
    return _AcquisitionBoundaryResult(_RescheduleAcquisitionBoundary,
        following)
end

function _process_acquisition_lifecycle_boundary!(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState, timestamp::PlantTimestamp,
    rng::AbstractRNG)
    if state.image_status == _FrameTransferImageActive
        _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
        close_exposure!(prepared, state, timestamp)
        return _AcquisitionBoundaryResult(_RescheduleAcquisitionBoundary,
            state.transfer_complete)
    end
    complete_frame_transfer!(prepared, state, timestamp)
    return _AcquisitionBoundaryResult(_ScheduleAcquisitionReadout,
        state.storage_readout_complete)
end

function _process_acquisition_lifecycle_boundary!(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp,
    rng::AbstractRNG)
    _integrate_event_acquisition_to!(prepared, state, timestamp, rng)
    close_exposure!(prepared, state, timestamp)
    return _AcquisitionBoundaryResult(_ScheduleAcquisitionReadout,
        state.readout_complete)
end

function _process_acquisition_boundary!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    result = _process_acquisition_lifecycle_boundary!(
        acquisition.lifecycle,
        acquisition_state, claim.key.timestamp, acquisition.rng)
    if result.disposition == _RescheduleAcquisitionBoundary
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
@inline _event_requires_readiness(
    ::PreparedDirectMeasurementAcquisition) = true

@inline _event_readiness_timestamp(
    ::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState) = state.readiness
@inline _event_readiness_timestamp(
    ::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState) = state.readiness
@inline _event_readiness_timestamp(
    ::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState) = state.readiness

@inline _event_product_sequence(
    state::GlobalShutterAcquisitionState) = state.sequence
@inline _event_product_sequence(
    state::RollingShutterAcquisitionState) = state.sequence
@inline _event_product_sequence(
    state::FrameTransferAcquisitionState) = state.product_sequence
@inline _event_product_sequence(
    state::DirectMeasurementAcquisitionState) = state.sequence

function _publish_acquisition_product!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, slot::UInt32, output)
    acquisition = _event_acquisition_binding(prepared, slot)
    product = acquisition.product
    product === output ||
        copy_acquisition_product!(product, output)
    index = Int(slot)
    sequence = _event_product_sequence(
        _event_acquisition_state(state, slot))
    previous_sequence = @inbounds state.product_sequences[index]
    sequence > previous_sequence ||
        _plant_event_loop_error(:product_sequence,
            "acquisition product sequence did not advance")
    @inbounds state.product_sequences[index] = sequence
    @inbounds state.product_ready_timestamps[index] =
        state.scheduler.current_timestamp
    return product
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

function _triggered_acquisition_slot_or_zero(
    prepared::PreparedPlantEventLoop,
    consumer::TriggerConsumerID)
    @inbounds for index in eachindex(prepared.acquisitions)
        start = prepared.acquisitions[index].start
        _start_matches_consumer(start, consumer) && return UInt32(index)
    end
    return UInt32(0)
end

@inline _start_matches_consumer(::PeriodicAcquisitionStart,
    ::TriggerConsumerID) = false
@inline _start_matches_consumer(start::TriggeredAcquisitionStart,
    consumer::TriggerConsumerID) = start.consumer == consumer

@inline _phase_reference_matches_consumer(
    ::FreeRunningPhaseReference, ::TriggerConsumerID) = false
@inline _phase_reference_matches_consumer(
    ::TriggerSourcePhaseReference, ::TriggerConsumerID) = false
@inline _phase_reference_matches_consumer(
    reference::TriggerResetPhaseReference,
    consumer::TriggerConsumerID) = reference.consumer == consumer

function _triggered_autonomous_optic_slot_or_zero(
    prepared::PreparedPlantEventLoop,
    consumer::TriggerConsumerID)
    @inbounds for index in eachindex(prepared.autonomous_optics)
        binding = prepared.autonomous_optics[index]
        _phase_reference_matches_consumer(binding.phase_reference,
            consumer) && return UInt32(index)
    end
    return UInt32(0)
end

@inline _phase_reference_matches_source(
    ::FreeRunningPhaseReference, ::TriggerSourceID) = false
@inline _phase_reference_matches_source(
    ::TriggerResetPhaseReference, ::TriggerSourceID) = false
@inline _phase_reference_matches_source(
    reference::TriggerSourcePhaseReference,
    source::TriggerSourceID) = reference.source == source

function _notify_autonomous_trigger_source!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    realization::TriggerSourceRealization)
    nominal = nominal_trigger_edge(realization)
    timestamp = realized_trigger_source_timestamp(realization)
    @inbounds for binding in prepared.autonomous_optics
        _phase_reference_matches_source(binding.phase_reference,
            nominal.source_id) || continue
        optic_state = state.controllable_optics[Int(binding.optic_slot)]
        reset_autonomous_periodic_optic_phase!(binding.implementation,
            optic_state, timestamp, nominal.sequence)
    end
    return nothing
end

function _reset_triggered_autonomous_optic!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    slot::UInt32, delivery::TriggerDelivery)
    binding = @inbounds prepared.autonomous_optics[Int(slot)]
    optic_state = @inbounds state.controllable_optics[
        Int(binding.optic_slot)]
    reset_autonomous_periodic_optic_phase!(binding.implementation,
        optic_state, delivered_trigger_edge(delivery).timestamp,
        nominal_trigger_edge(delivery).sequence)
    return nothing
end

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
    autonomous_slot = UInt32(0)
    activation_timestamp = zero(PlantTimestamp)
    if source_due
        realized_trigger_source_timestamp(source) == claim.key.timestamp ||
            _plant_event_loop_error(:trigger_schedule,
                "trigger source does not match its scheduler claim")
        realization = realize_next_trigger_source!(workspace.trigger,
            topology, trigger_state)
        _notify_autonomous_trigger_source!(prepared, state, realization)
    else
        delivered = delivered_trigger_edge(delivery)
        delivered.timestamp == claim.key.timestamp ||
            _plant_event_loop_error(:trigger_schedule,
                "trigger delivery does not match its scheduler claim")
        activated_slot = _triggered_acquisition_slot_or_zero(prepared,
            trigger_delivery_consumer(delivery))
        autonomous_slot = _triggered_autonomous_optic_slot_or_zero(prepared,
            trigger_delivery_consumer(delivery))
        xor(iszero(activated_slot), iszero(autonomous_slot)) ||
            _plant_event_loop_error(:trigger_binding,
                "delivered trigger consumer $(trigger_delivery_consumer(delivery)) must bind exactly one event owner")
        if !iszero(activated_slot)
            acquisition = _event_acquisition_binding(prepared,
                activated_slot)
            _require_inactive_event_generator(prepared, state,
                acquisition.start_handle, delivered.timestamp)
        end
        pop_next_trigger_delivery!(workspace.delivery, topology,
            trigger_state) || _plant_event_loop_error(:trigger_schedule,
            "due trigger delivery disappeared before removal")
        activation_timestamp = delivered.timestamp
        if !iszero(autonomous_slot)
            _reset_triggered_autonomous_optic!(prepared, state,
                autonomous_slot, delivery)
        end
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
    state::_AcquisitionEventLifecycleState)
    _require_event_lifecycle_binding(acquisition.lifecycle, state)
    return nothing
end

@inline function _require_event_lifecycle_binding(
    prepared::Union{
        PreparedGlobalShutterAcquisition,
        PreparedRollingShutterAcquisition,
        PreparedFrameTransferAcquisition,
    },
    state::Union{
        GlobalShutterAcquisitionState,
        RollingShutterAcquisitionState,
        FrameTransferAcquisitionState,
    })
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

@inline function _require_event_lifecycle_binding(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState)
    _require_direct_measurement_binding(prepared, state)
    measurement = prepared.measurement
    storage = measurement_storage(measurement)
    length(storage) == length(prepared.instantaneous_sample) ==
        length(prepared.integrated_sample) ||
        _plant_event_loop_error(:prepared_binding,
            "direct-measurement lifecycle storage changed after event-loop preparation")
    typeof(backend(storage)) ===
        typeof(backend(prepared.instantaneous_sample)) &&
        plane_device(storage) ==
            plane_device(prepared.instantaneous_sample) ||
        _plant_event_loop_error(:prepared_binding,
            "direct-measurement lifecycle memory domain changed after preparation")
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
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    timestamp::PlantTimestamp)
    state.status == DirectMeasurementExposureActive || return nothing
    state.integrated_through <= timestamp ||
        _plant_event_loop_error(:time_regression,
            "direct-measurement integration target precedes acquisition progress")
    state.integrated_through == timestamp && return nothing
    timestamp <= state.exposure_close ||
        _direct_measurement_event_error(:interval_after_close,
            "direct-measurement integration extends beyond exposure close")
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
        binding = prepared.paths[index]
        binding.requires_full_optical || continue
        path = binding.path
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
        binding.requires_full_optical || continue
        path = binding.path
        materialize_path_input_rngs!(path.materialization, path.input,
            atmosphere, epoch, binding.rngs)
    end
    return nothing
end

function _apply_due_controllable_optics!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    due_paths::Memory{Bool})
    isempty(prepared.optics) && return nothing
    @inbounds for path_index in eachindex(prepared.paths)
        due_paths[path_index] || continue
        path = prepared.paths[path_index]
        path.requires_full_optical || continue
        input = path.path.input
        for optic_index in eachindex(prepared.optics)
            optic = prepared.optics[optic_index]
            _apply_event_controllable_optic_surface!(
                controllable_optic_execution_role(optic.implementation),
                input, optic.implementation,
                state.controllable_optics[optic_index])
        end
    end
    return nothing
end

@inline function _apply_event_controllable_optic_surface!(
    ::PupilSurfaceExecutionRole, input, implementation, state)
    return apply_controllable_optic_surface!(input, implementation, state)
end

@inline function _apply_event_controllable_optic_surface!(
    ::AutonomousPathExecutionRole, input, implementation, state)
    return nothing
end

function _evaluate_due_autonomous_optics!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    due_paths::Memory{Bool}, timestamp::PlantTimestamp)
    @inbounds for binding in prepared.autonomous_optics
        due_paths[Int(binding.path_slot)] || continue
        optic_state = state.controllable_optics[Int(binding.optic_slot)]
        evaluate_autonomous_periodic_optic!(binding.implementation,
            optic_state, binding.coupling, timestamp)
    end
    return nothing
end

function _execute_due_paths!(prepared::PreparedPlantEventLoop,
    due_paths::Memory{Bool})
    @inbounds for index in eachindex(prepared.paths)
        due_paths[index] || continue
        binding = prepared.paths[index]
        binding.requires_full_optical || continue
        execute_path!(binding.path, binding.rngs)
    end
    return nothing
end

@inline function _evaluate_event_sample!(
    ::Nothing, ::AbstractPreparedAcquisitionLifecycle,
    ::PlantTimestamp, command_applications)
    return nothing
end

@inline function _evaluate_event_sample!(
    provider::PreparedLinearReducedOrderEventProvider,
    lifecycle::PreparedDirectMeasurementAcquisition,
    timestamp::PlantTimestamp, command_applications)
    return evaluate_linear_reduced_order_sample!(
        lifecycle.instantaneous_sample, provider, timestamp,
        command_applications)
end

function _evaluate_due_reduced_order_samples!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    due_paths::Memory{Bool}, timestamp::PlantTimestamp)
    @inbounds for acquisition in prepared.acquisitions
        due_paths[Int(acquisition.path_slot)] || continue
        _evaluate_event_sample!(acquisition.sample_provider,
            acquisition.lifecycle, timestamp, state.command_applications)
    end
    return nothing
end

function _mark_due_paths_sampled!(state::PlantEventLoopState,
    due_paths::Memory{Bool})
    @inbounds for index in eachindex(due_paths)
        due_paths[index] && (state.path_sampled[index] = true)
    end
    return nothing
end

function _due_full_optical_path(
    prepared::PreparedPlantEventLoop, due_paths::Memory{Bool})
    @inbounds for index in eachindex(prepared.paths)
        due_paths[index] &&
            prepared.paths[index].requires_full_optical && return true
    end
    return false
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
    _preflight_due_path_consumers(prepared, state, workspace.due_paths,
        timestamp)
    _integrate_due_path_consumers!(prepared, state, workspace.due_paths,
        timestamp)
    if _due_full_optical_path(prepared, workspace.due_paths)
        target_time = _preflight_atmosphere_time(atmosphere, timestamp)
        epoch = advance_to!(atmosphere, target_time,
            prepared.atmosphere_rng)
        _validate_due_path_materializations!(prepared,
            workspace.due_paths, atmosphere, epoch)
        _materialize_due_paths!(prepared, workspace.due_paths, atmosphere,
            epoch)
        _apply_due_controllable_optics!(prepared, state,
            workspace.due_paths)
        _evaluate_due_autonomous_optics!(prepared, state,
            workspace.due_paths, timestamp)
        _execute_due_paths!(prepared, workspace.due_paths)
    end
    _evaluate_due_reduced_order_samples!(prepared, state,
        workspace.due_paths, timestamp)
    _mark_due_paths_sampled!(state, workspace.due_paths)
    _resolve_due_path_claims!(prepared, state, workspace, timestamp)
    return nothing
end

function _process_ordinary_event!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, workspace::PlantEventLoopWorkspace,
    claim::EventClaim, action::_PlantEventAction)
    kind = action.kind
    kind == _TriggerTopologyAction &&
        return _process_trigger_topology!(prepared, state, workspace, claim)
    kind == _CommandEndpointAction &&
        return _process_command_endpoint!(prepared, state, workspace, claim,
            action)
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
