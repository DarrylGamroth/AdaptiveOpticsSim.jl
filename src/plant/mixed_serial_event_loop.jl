#
# Deterministic mixed-target serial event composition
#
# This coordinator closes the serial correctness oracle over one
# global scheduler, atmosphere authority, command authority, and product
# sequence authority while entering each path/acquisition target explicitly.
# It owns no task, queue, pacing policy, affinity policy, or transport.
#

mutable struct _PreparedMixedResourcePlantEventLoopBinding end

struct _PreparedMixedSerialEventCommand
    id::CommandEndpointID
    handle::EventGeneratorHandle
end

struct _PreparedMixedSerialPathSchedule
    id::OpticalPathID
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    handle::EventGeneratorHandle
    acquisition_slots::Memory{UInt32}
end

abstract type _AbstractPreparedMixedSerialEventAcquisition end

struct _PreparedMixedSerialEventAcquisition{
    R,L,P,G,S,
} <: _AbstractPreparedMixedSerialEventAcquisition
    id::AcquisitionID
    resources::R
    lifecycle::L
    product::P
    rng::G
    start::S
    path_slot::UInt32
    start_handle::EventGeneratorHandle
    boundary_handle::EventGeneratorHandle
    band_open_handle::EventGeneratorHandle
    readout_handle::EventGeneratorHandle
    readiness_handle::EventGeneratorHandle
end

"""
Prepared deterministic event composition for one statically partitioned
mixed-resource plant.

The loop retains one scheduler, atmosphere authority, command authority, and
product-sequence authority while executing each complete path/acquisition group
on its exact prepared target. It is caller-driven and owns no task, queue,
pacing policy, affinity policy, or transport.
"""
struct PreparedMixedResourcePlantEventLoop{B,E,T}
    binding::B
    execution::E
    scheduler::PreparedEventScheduler
    actions::Memory{_PlantEventAction}
    commands::Memory{_PreparedMixedSerialEventCommand}
    paths::Memory{_PreparedMixedSerialPathSchedule}
    acquisitions::Memory{_AbstractPreparedMixedSerialEventAcquisition}
    trigger_topology::T
end

@enum _MixedSerialEventLoopPhase::UInt8 begin
    _MixedSerialEventLoopIdle = 0x01
    _MixedSerialEventLoopProcessing = 0x02
    _MixedSerialEventLoopFailed = 0x03
end

"""Single-writer mutable state for a `PreparedMixedResourcePlantEventLoop`."""
mutable struct MixedResourcePlantEventLoopState{B,E,T}
    binding::B
    scheduler::EventSchedulerState
    execution::E
    acquisitions::Memory{_AcquisitionEventLifecycleState}
    path_sampled::Memory{Bool}
    product_sequences::Memory{UInt64}
    product_ready_timestamps::Memory{PlantTimestamp}
    trigger::T
    phase::_MixedSerialEventLoopPhase
end

"""Fixed-capacity scratch storage for a `PreparedMixedResourcePlantEventLoop`."""
mutable struct MixedResourcePlantEventLoopWorkspace{B,E,T}
    binding::B
    scheduler::EventSchedulerWorkspace
    execution::E
    trigger::T
    delivery::Base.RefValue{TriggerDelivery}
    event_timestamp::Base.RefValue{PlantTimestamp}
    command_dispositions::Memory{PlantCommandDisposition}
    command_disposition_count::Int
end

@inline function _mixed_resource_plant_partitions(
    prepared::PreparedMixedResourcePlantEventLoop,
)
    return prepared.execution.partitions
end

@inline plant_event_path_count(
    prepared::PreparedMixedResourcePlantEventLoop) = length(prepared.paths)
@inline plant_event_acquisition_count(
    prepared::PreparedMixedResourcePlantEventLoop) =
    length(prepared.acquisitions)
@inline plant_event_generator_count(
    prepared::PreparedMixedResourcePlantEventLoop) =
    event_generator_count(prepared.scheduler)
@inline plant_event_command_endpoint_count(
    prepared::PreparedMixedResourcePlantEventLoop) =
    length(prepared.commands)
@inline plant_event_controllable_optic_count(
    prepared::PreparedMixedResourcePlantEventLoop) = length(
    controllable_optic_definitions(
        plant_definition(_mixed_resource_plant_partitions(prepared))))
@inline plant_event_sampled_aberration_count(
    prepared::PreparedMixedResourcePlantEventLoop) = length(
    sampled_aberration_definitions(
        plant_definition(_mixed_resource_plant_partitions(prepared))))
@inline plant_event_autonomous_optic_count(
    ::PreparedMixedResourcePlantEventLoop) = 0
@inline atmosphere_authority_target(
    prepared::PreparedMixedResourcePlantEventLoop) =
    atmosphere_authority_target(_mixed_resource_plant_partitions(prepared))
@inline command_authority_target(
    prepared::PreparedMixedResourcePlantEventLoop) =
    command_authority_target(_mixed_resource_plant_partitions(prepared))

@noinline function _mixed_serial_event_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(
        :mixed_resource_plant_event_loop, reason, String(message)))
end

@noinline function _mixed_serial_event_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantScheduleError(
        :mixed_resource_plant_event_loop, reason, String(message)))
end

function _mixed_serial_event_path_slot(
    ids::_FixedPlantRegistry{OpticalPathID},
    id::OpticalPathID,
)
    @inbounds for slot in eachindex(ids)
        ids[slot] == id && return slot
    end
    _mixed_serial_event_preparation_error(
        :missing_path_schedule,
        "acquisition references optical path $id without a sample schedule",
    )
end

function _mixed_serial_event_acquisition_resource(
    partitions::PreparedPlantPartitions,
    id::AcquisitionID,
)
    found_partition = nothing
    found_resources = nothing
    found_slot = 0
    for partition in prepared_partitions(partitions)
        @inbounds for slot in eachindex(prepared_acquisitions(partition))
            resources = prepared_acquisitions(partition)[slot]
            acquisition_id(resources.definition) == id || continue
            isnothing(found_resources) || _mixed_serial_event_preparation_error(
                :duplicate_acquisition,
                "prepared partitions contain more than one acquisition $id",
            )
            found_partition = partition
            found_resources = resources
            found_slot = slot
        end
    end
    isnothing(found_resources) && _mixed_serial_event_preparation_error(
        :unknown_acquisition,
        "prepared partitions contain no acquisition $id",
    )
    return found_partition, found_resources, found_slot
end

@inline function _mixed_serial_frame_execution(
    provider::PreparedFullOpticalProvider{<:FrameAcquisitionExecution},
)
    return provider.execution
end

function _mixed_serial_frame_execution(provider)
    _mixed_serial_event_preparation_error(
        :unsupported_acquisition,
        "mixed serial event composition currently requires a full-optical " *
        "FrameAcquisitionExecution provider; got $(typeof(provider))",
    )
end

@inline function _prepare_mixed_serial_detector_lifecycle(
    execution::FrameAcquisitionExecution,
    result::IntensityMap,
    definition::AbstractDetectorAcquisitionLifecycleDefinition,
)
    return _prepare_detector_event_lifecycle(execution, result, definition)
end

function _prepare_mixed_serial_detector_lifecycle(
    ::FrameAcquisitionExecution,
    result,
    ::AbstractDetectorAcquisitionLifecycleDefinition,
)
    _mixed_serial_event_preparation_error(
        :unsupported_path_result,
        "mixed serial detector lifecycle requires an IntensityMap path " *
        "result; got $(typeof(result))",
    )
end

function _prepare_mixed_serial_event_acquisition(
    partitions::PreparedPlantPartitions,
    definition::AcquisitionEventDefinition,
    path_ids::_FixedPlantRegistry{OpticalPathID},
    scheduler::PreparedEventScheduler,
    topology,
    ordinal::Int,
)
    _require_start_trigger_topology(definition.start, topology)
    partition, resources, resource_slot =
        _mixed_serial_event_acquisition_resource(
            partitions, definition.acquisition)
    path = acquisition_path_id(resources.definition)
    path_slot = _mixed_serial_event_path_slot(path_ids, path)
    partition_target(getfield(partitions, :assignment), path) ==
        compute_device(partition) || _mixed_serial_event_preparation_error(
            :acquisition_placement,
            "acquisition $(definition.acquisition) is not co-located with " *
            "its optical path $path",
        )
    provider = acquisition_provider(resources)
    implementation = getfield(provider, :implementation)
    execution = _mixed_serial_frame_execution(implementation)
    lifecycle = _with_completed_prepared_device_execution_context(
        resources.context) do
        _prepare_mixed_serial_detector_lifecycle(
            execution, resources.path_result, definition.lifecycle)
    end
    _require_periodic_start_spacing(lifecycle, definition.start)
    acquisition_rngs = @inbounds getfield(partition, :rngs).acquisitions[
        resource_slot]
    _require_rng_owner_binding(acquisition_rngs, resources)
    rng = rng_stream_state(acquisition_rngs, Val(:detector))
    return _PreparedMixedSerialEventAcquisition(
        definition.acquisition,
        resources,
        lifecycle,
        acquisition_observation(resources),
        rng,
        definition.start,
        UInt32(path_slot),
        event_generator_handle(
            scheduler, ExposureOpenPhase, 2ordinal - 1),
        event_generator_handle(
            scheduler, IntegrationBoundaryPhase, ordinal),
        event_generator_handle(
            scheduler, ExposureOpenPhase, 2ordinal),
        event_generator_handle(
            scheduler, ReadoutCompletionPhase, ordinal),
        event_generator_handle(
            scheduler, AcquisitionReadyPhase, ordinal),
    )
end

function _prepare_mixed_serial_event_acquisitions(
    partitions::PreparedPlantPartitions,
    definitions,
    path_ids::_FixedPlantRegistry{OpticalPathID},
    scheduler::PreparedEventScheduler,
    topology,
)
    acquisitions = Memory{_AbstractPreparedMixedSerialEventAcquisition}(
        undef, length(definitions))
    @inbounds for ordinal in eachindex(definitions)
        acquisitions[ordinal] = _prepare_mixed_serial_event_acquisition(
            partitions,
            definitions[ordinal],
            path_ids,
            scheduler,
            topology,
            ordinal,
        )
    end
    return acquisitions
end

function _mixed_serial_path_acquisition_slots(
    acquisitions::Memory{_AbstractPreparedMixedSerialEventAcquisition},
    path_slot::Int,
)
    slots = UInt32[]
    @inbounds for acquisition_slot in eachindex(acquisitions)
        Int(acquisitions[acquisition_slot].path_slot) == path_slot || continue
        push!(slots, UInt32(acquisition_slot))
    end
    return Memory{UInt32}(slots)
end

function _prepare_mixed_serial_path_schedules(
    definitions,
    acquisitions::Memory{_AbstractPreparedMixedSerialEventAcquisition},
    scheduler::PreparedEventScheduler,
)
    paths = Memory{_PreparedMixedSerialPathSchedule}(
        undef, length(definitions))
    @inbounds for ordinal in eachindex(definitions)
        definition = definitions[ordinal]
        paths[ordinal] = _PreparedMixedSerialPathSchedule(
            definition.path,
            definition.schedule,
            definition.origin,
            event_generator_handle(
                scheduler, OpticalSamplePhase, ordinal),
            _mixed_serial_path_acquisition_slots(
                acquisitions, ordinal),
        )
    end
    return paths
end

function _prepare_mixed_serial_event_commands(
    execution::_PreparedMixedSerialExecution,
    scheduler::PreparedEventScheduler,
)
    lanes = execution.fanout.lanes
    commands = Memory{_PreparedMixedSerialEventCommand}(
        undef, length(lanes))
    for (ordinal, lane) in enumerate(lanes)
        commands[ordinal] = _PreparedMixedSerialEventCommand(
            command_endpoint_id(lane.authority_endpoint),
            event_generator_handle(
                scheduler, CommandApplicationPhase, ordinal),
        )
    end
    return commands
end

function prepare_plant_event_loop(
    partitions::PreparedPlantPartitions,
    definition::PlantEventLoopDefinition,
)
    isempty(definition.autonomous_optics) ||
        _mixed_serial_event_preparation_error(
            :unsupported_autonomous_optics,
            "mixed serial event composition does not yet support " *
            "autonomous periodic optics",
        )
    sample_definitions = _sorted_optical_sample_definitions(
        definition.optical_samples)
    acquisition_definitions = _sorted_acquisition_event_definitions(
        definition.acquisition_events)
    _require_unique_trigger_consumers(acquisition_definitions, ())
    _require_bound_trigger_consumers(
        acquisition_definitions, (), definition.trigger_topology)
    path_ids = OpticalPathID[
        sample.path for sample in sample_definitions]
    execution = _prepare_mixed_serial_execution(partitions, path_ids)
    generator_definitions = EventGeneratorDefinition[]
    _append_event_generator_definitions!(
        generator_definitions,
        sample_definitions,
        acquisition_definitions,
        execution.fanout.lanes,
        definition.trigger_topology,
    )
    scheduler = prepare_event_scheduler(
        generator_definitions; capacity=length(generator_definitions))
    acquisitions = _prepare_mixed_serial_event_acquisitions(
        partitions,
        acquisition_definitions,
        execution.path_ids,
        scheduler,
        definition.trigger_topology,
    )
    return PreparedMixedResourcePlantEventLoop(
        _PreparedMixedResourcePlantEventLoopBinding(),
        execution,
        scheduler,
        _prepared_event_actions(scheduler),
        _prepare_mixed_serial_event_commands(execution, scheduler),
        _prepare_mixed_serial_path_schedules(
            sample_definitions, acquisitions, scheduler),
        acquisitions,
        _prepared_trigger_topology(definition.trigger_topology),
    )
end

@inline function _mixed_serial_event_acquisition_state(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
)
    return _event_acquisition_state(acquisition.lifecycle)
end

@inline _mixed_serial_event_trigger_state(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyState()
@inline _mixed_serial_event_trigger_state(
    topology::PreparedTriggerTopology) = TriggerTopologyState(topology)

@inline function _next_mixed_serial_command_timestamp(
    ::Tuple{},
    ::Tuple{},
    id::CommandEndpointID,
)
    _mixed_serial_event_error(
        :unknown_command_endpoint,
        "mixed serial event loop contains no command endpoint $id",
    )
end

@inline function _next_mixed_serial_command_timestamp(
    lanes::Tuple,
    lane_states::Tuple,
    id::CommandEndpointID,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    if command_endpoint_id(lane.authority_endpoint) == id
        endpoint = lane.authority_endpoint.endpoint
        key = next_command_order_key(endpoint, lane_state.endpoint)
        command_timestamp = key === nothing ? nothing :
            command_scheduled_timestamp(key)
        silence_timestamp = next_command_silence_timestamp(
            endpoint, lane_state.endpoint, lane_state.application)
        command_timestamp === nothing && return silence_timestamp
        silence_timestamp === nothing && return command_timestamp
        return min(command_timestamp, silence_timestamp)
    end
    return _next_mixed_serial_command_timestamp(
        Base.tail(lanes), Base.tail(lane_states), id)
end

@inline function _next_mixed_serial_command_timestamp(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    id::CommandEndpointID,
)
    return _next_mixed_serial_command_timestamp(
        prepared.execution.fanout.lanes,
        state.execution.fanout.lanes,
        id,
    )
end

function _schedule_mixed_serial_command_endpoint!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    slot::Int,
)
    command = @inbounds prepared.commands[slot]
    desired = _next_mixed_serial_command_timestamp(
        prepared, state, command.id)
    cursor = @inbounds state.scheduler.cursors[Int(command.handle.slot)]
    if desired === nothing
        cursor.status == _ScheduledEventGenerator &&
            deactivate_event_generator!(
                prepared.scheduler, state.scheduler, command.handle)
        return nothing
    end
    if cursor.status == _InactiveEventGenerator
        activate_event_generator!(
            prepared.scheduler, state.scheduler, command.handle, desired)
        return desired
    end
    cursor.status == _ScheduledEventGenerator ||
        _mixed_serial_event_error(
            :command_generator_state,
            "command endpoint generator cannot be changed while claimed",
        )
    cursor.next_timestamp == desired && return desired
    deactivate_event_generator!(
        prepared.scheduler, state.scheduler, command.handle)
    activate_event_generator!(
        prepared.scheduler, state.scheduler, command.handle, desired)
    return desired
end

function _schedule_all_mixed_serial_command_endpoints!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
)
    @inbounds for slot in eachindex(prepared.commands)
        _schedule_mixed_serial_command_endpoint!(
            prepared, state, slot)
    end
    return nothing
end

function MixedResourcePlantEventLoopState(
    prepared::PreparedMixedResourcePlantEventLoop,
)
    scheduler = EventSchedulerState(prepared.scheduler)
    execution = _prepare_mixed_serial_execution_state(prepared.execution)
    acquisitions = Memory{_AcquisitionEventLifecycleState}(
        undef, length(prepared.acquisitions))
    @inbounds for slot in eachindex(acquisitions)
        acquisitions[slot] = _mixed_serial_event_acquisition_state(
            prepared.acquisitions[slot])
    end
    path_sampled = Memory{Bool}(undef, length(prepared.paths))
    fill!(path_sampled, false)
    product_sequences = Memory{UInt64}(
        undef, length(prepared.acquisitions))
    fill!(product_sequences, zero(UInt64))
    product_ready_timestamps = Memory{PlantTimestamp}(
        undef, length(prepared.acquisitions))
    fill!(product_ready_timestamps, zero(PlantTimestamp))
    state = MixedResourcePlantEventLoopState(
        prepared.binding,
        scheduler,
        execution,
        acquisitions,
        path_sampled,
        product_sequences,
        product_ready_timestamps,
        _mixed_serial_event_trigger_state(prepared.trigger_topology),
        _MixedSerialEventLoopIdle,
    )
    _schedule_all_mixed_serial_command_endpoints!(prepared, state)
    return state
end

@inline _mixed_serial_event_trigger_workspace(
    ::_NoPreparedTriggerTopology) = _NoTriggerTopologyWorkspace()
@inline _mixed_serial_event_trigger_workspace(
    topology::PreparedTriggerTopology) = TriggerTopologyWorkspace(topology)

@inline _mixed_serial_disposition_capacity(::Tuple{}) = 0

@inline function _mixed_serial_disposition_capacity(workspaces::Tuple)
    return length(first(workspaces).disposition.dispositions) +
        _mixed_serial_disposition_capacity(Base.tail(workspaces))
end

function MixedResourcePlantEventLoopWorkspace(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
)
    prepared.binding === state.binding ||
        _mixed_serial_event_preparation_error(
            :foreign_state,
            "mixed serial event-loop state belongs to another preparation",
        )
    execution = _prepare_mixed_serial_execution_workspace(
        prepared.execution, state.execution)
    capacity = _mixed_serial_disposition_capacity(execution.fanout.lanes)
    return MixedResourcePlantEventLoopWorkspace(
        prepared.binding,
        EventSchedulerWorkspace(prepared.scheduler),
        execution,
        _mixed_serial_event_trigger_workspace(prepared.trigger_topology),
        Ref{TriggerDelivery}(),
        Ref(zero(PlantTimestamp)),
        Memory{PlantCommandDisposition}(undef, capacity),
        0,
    )
end

@inline function _require_mixed_serial_event_ownership(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    prepared.binding === state.binding === workspace.binding ||
        _mixed_serial_event_error(
            :foreign_owner,
            "mixed serial event-loop plan, state, and workspace belong to " *
            "different owners",
        )
    _require_scheduler_binding(prepared.scheduler, state.scheduler)
    _require_scheduler_binding(prepared.scheduler, workspace.scheduler)
    prepared.execution.binding === state.execution.binding ===
        workspace.execution.binding || _mixed_serial_event_error(
            :prepared_binding,
            "mixed serial event-loop optical execution binding changed",
        )
    length(prepared.paths) == length(state.path_sampled) ==
        length(workspace.execution.due_paths) ||
        _mixed_serial_event_error(
            :prepared_binding,
            "mixed serial event-loop path capacity changed after preparation",
        )
    length(prepared.acquisitions) == length(state.acquisitions) ==
        length(state.product_sequences) ==
        length(state.product_ready_timestamps) ||
        _mixed_serial_event_error(
            :prepared_binding,
            "mixed serial event-loop acquisition capacity changed after " *
            "preparation",
        )
    return nothing
end

@inline function _require_mixed_serial_event_idle(
    state::MixedResourcePlantEventLoopState,
)
    state.phase == _MixedSerialEventLoopIdle ||
        _mixed_serial_event_error(
            state.phase == _MixedSerialEventLoopFailed ?
                :event_loop_failed : :event_loop_busy,
            "mixed serial event loop is not idle ($(state.phase))",
        )
    return nothing
end

@inline function _mixed_serial_command_slot(
    prepared::PreparedMixedResourcePlantEventLoop,
    id::CommandEndpointID,
)
    @inbounds for slot in eachindex(prepared.commands)
        prepared.commands[slot].id == id && return slot
    end
    _mixed_serial_event_error(
        :unknown_command_endpoint,
        "mixed serial event loop contains no command endpoint $id",
    )
end

@inline function _append_mixed_serial_dispositions!(
    workspace::MixedResourcePlantEventLoopWorkspace,
    source::CommandDispositionWorkspace,
)
    count = command_disposition_count(source)
    first_slot = workspace.command_disposition_count + 1
    last_slot = workspace.command_disposition_count + count
    last_slot <= length(workspace.command_dispositions) ||
        _mixed_serial_event_error(
            :command_disposition_capacity,
            "mixed serial command disposition capacity was exceeded",
        )
    @inbounds for source_slot in 1:count
        workspace.command_dispositions[first_slot + source_slot - 1] =
            command_disposition(source, source_slot)
    end
    workspace.command_disposition_count = last_slot
    clear_command_dispositions!(source)
    return nothing
end

@inline function _drain_mixed_serial_fanout_dispositions!(
    ::Tuple{},
    ::MixedResourcePlantEventLoopWorkspace,
)
    return nothing
end

@inline function _drain_mixed_serial_fanout_dispositions!(
    lanes::Tuple,
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    _append_mixed_serial_dispositions!(workspace, first(lanes).disposition)
    return _drain_mixed_serial_fanout_dispositions!(
        Base.tail(lanes), workspace)
end

@inline function _require_empty_mixed_serial_dispositions(
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    iszero(workspace.command_disposition_count) ||
        _mixed_serial_event_error(
            :unconsumed_command_dispositions,
            "clear mixed serial command dispositions before another " *
            "command admission",
        )
    return nothing
end

@inline command_disposition_count(
    workspace::MixedResourcePlantEventLoopWorkspace) =
    workspace.command_disposition_count

function command_disposition(
    workspace::MixedResourcePlantEventLoopWorkspace,
    index::Integer,
)
    1 <= index <= workspace.command_disposition_count ||
        _mixed_serial_event_error(
            :invalid_disposition_index,
            "command disposition index must be within the current records",
        )
    return @inbounds workspace.command_dispositions[Int(index)]
end

@inline command_disposition(
    ::MixedResourcePlantEventLoopWorkspace,
    ::Bool,
) = _mixed_serial_event_error(
    :invalid_disposition_index,
    "command disposition index must be an integer count, not Bool",
)

function clear_command_dispositions!(
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    workspace.command_disposition_count = 0
    return workspace
end

@inline function _require_mixed_serial_command_admission_timestamp(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    current = scheduler_timestamp(state.scheduler)
    timestamp < current && _mixed_serial_event_error(
        :command_admission_time_regression,
        "command admission timestamp precedes the current mixed serial " *
        "event timestamp",
    )
    timestamp == current && state.scheduler.has_last_key &&
        _mixed_serial_event_error(
            :command_admission_time_elapsed,
            "command admission timestamp has already been processed by the " *
            "mixed serial event loop",
        )
    count = scan_due_events!(
        workspace.scheduler, prepared.scheduler, state.scheduler)
    !iszero(count) && workspace.scheduler.due_timestamp < timestamp &&
        _mixed_serial_event_error(
            :command_admission_overtakes_event,
            "command admission timestamp follows the next unprocessed " *
            "plant event at $(workspace.scheduler.due_timestamp)",
        )
    return nothing
end

@inline function _admit_mixed_serial_command_lane!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    command::PlantCommand,
    ::PlantTimestamp,
)
    _mixed_serial_event_error(
        :unknown_command_endpoint,
        "mixed serial command authority contains no endpoint " *
        "$(command_endpoint_id(command))",
    )
end

@inline function _admit_mixed_serial_command_lane!(
    lanes::Tuple,
    lane_states::Tuple,
    lane_workspaces::Tuple,
    command::PlantCommand,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    lane_workspace = first(lane_workspaces)
    if command_endpoint_id(lane.authority_endpoint) ==
            command_endpoint_id(command)
        return admit_plant_command!(
            lane_workspace.disposition,
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            command,
            timestamp,
        )
    end
    return _admit_mixed_serial_command_lane!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(lane_workspaces),
        command,
        timestamp,
    )
end

function admit_plant_command!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    command::PlantCommand,
    timestamp::PlantTimestamp,
)
    _require_mixed_serial_event_ownership(prepared, state, workspace)
    _require_mixed_serial_event_idle(state)
    _require_empty_mixed_serial_dispositions(workspace)
    _require_mixed_serial_command_admission_timestamp(
        prepared, state, workspace, timestamp)
    slot = _mixed_serial_command_slot(
        prepared, command_endpoint_id(command))
    fanout = prepared.execution.fanout
    fanout_state = state.execution.fanout
    fanout_workspace = workspace.execution.fanout
    admission = try
        _with_completed_prepared_device_execution_context(
            fanout.authority.context) do
            _admit_mixed_serial_command_lane!(
                fanout.lanes,
                fanout_state.lanes,
                fanout_workspace.lanes,
                command,
                timestamp,
            )
        end
    catch
        _drain_mixed_serial_fanout_dispositions!(
            fanout_workspace.lanes, workspace)
        _schedule_mixed_serial_command_endpoint!(prepared, state, slot)
        rethrow()
    end
    _drain_mixed_serial_fanout_dispositions!(
        fanout_workspace.lanes, workspace)
    _schedule_mixed_serial_command_endpoint!(prepared, state, slot)
    return admission
end

function admit_plant_command_transaction!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
)
    _require_mixed_serial_event_ownership(prepared, state, workspace)
    _require_mixed_serial_event_idle(state)
    _require_empty_mixed_serial_dispositions(workspace)
    _require_mixed_serial_command_admission_timestamp(
        prepared, state, workspace, timestamp)
    admission = try
        _admit_command_fanout_transaction!(
            prepared.execution.fanout,
            state.execution.fanout,
            workspace.execution.fanout,
            transaction,
            timestamp,
        )
    catch
        _drain_mixed_serial_fanout_dispositions!(
            workspace.execution.fanout.lanes, workspace)
        _schedule_all_mixed_serial_command_endpoints!(prepared, state)
        rethrow()
    end
    _drain_mixed_serial_fanout_dispositions!(
        workspace.execution.fanout.lanes, workspace)
    _schedule_all_mixed_serial_command_endpoints!(prepared, state)
    return admission
end

@inline function _mixed_serial_event_action(
    prepared::PreparedMixedResourcePlantEventLoop,
    claim::EventClaim,
)
    slot = Int(claim.slot)
    1 <= slot <= length(prepared.actions) ||
        _mixed_serial_event_error(
            :invalid_action,
            "event claim does not map to a prepared mixed serial action",
        )
    return @inbounds prepared.actions[slot]
end

@inline function _mixed_serial_event_acquisition(
    prepared::PreparedMixedResourcePlantEventLoop,
    slot::UInt32,
)
    index = Int(slot)
    1 <= index <= length(prepared.acquisitions) ||
        _mixed_serial_event_error(
            :invalid_action,
            "event action contains an invalid acquisition slot",
        )
    return @inbounds prepared.acquisitions[index]
end

@inline function _mixed_serial_event_acquisition_state(
    state::MixedResourcePlantEventLoopState,
    slot::UInt32,
)
    index = Int(slot)
    1 <= index <= length(state.acquisitions) ||
        _mixed_serial_event_error(
            :invalid_action,
            "event action contains an invalid acquisition-state slot",
        )
    return @inbounds state.acquisitions[index]
end

@inline function _mixed_serial_event_path(
    prepared::PreparedMixedResourcePlantEventLoop,
    slot::UInt32,
)
    index = Int(slot)
    1 <= index <= length(prepared.paths) ||
        _mixed_serial_event_error(
            :invalid_action,
            "event action contains an invalid optical-path slot",
        )
    return @inbounds prepared.paths[index]
end

@inline function _require_inactive_mixed_serial_generator(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    handle::EventGeneratorHandle,
    timestamp::PlantTimestamp,
)
    slot = _require_event_generator_slot(prepared.scheduler, handle)
    cursor = @inbounds state.scheduler.cursors[slot]
    cursor.status == _InactiveEventGenerator ||
        _mixed_serial_event_error(
            :generator_busy,
            "required mixed serial event generator is already active",
        )
    definition = @inbounds prepared.scheduler.definitions[slot]
    key = PlantEventKey(
        timestamp,
        definition.phase,
        definition.ordinal,
        cursor.next_occurrence,
        _PLANT_TIME_TOKEN,
    )
    _require_forward_event_key(state.scheduler, key)
    return key
end

@inline function _mixed_serial_generator_due_at(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    handle::EventGeneratorHandle,
    timestamp::PlantTimestamp,
)
    slot = _require_event_generator_slot(prepared.scheduler, handle)
    cursor = @inbounds state.scheduler.cursors[slot]
    return cursor.status == _ScheduledEventGenerator &&
        cursor.next_timestamp == timestamp
end

@inline function _require_mixed_serial_event_path_available(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    timestamp::PlantTimestamp,
)
    path_slot = Int(acquisition.path_slot)
    @inbounds state.path_sampled[path_slot] && return nothing
    path = @inbounds prepared.paths[path_slot]
    _mixed_serial_generator_due_at(
        prepared, state, path.handle, timestamp) ||
        _mixed_serial_event_error(
            :uninitialized_path,
            "acquisition $(acquisition.id) begins before its first optical " *
            "sample",
        )
    return nothing
end

@inline function _resolve_mixed_serial_start_claim!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    ::PeriodicAcquisitionStart,
    next_timestamp::PlantTimestamp,
)
    return reschedule_event!(
        prepared.scheduler, state.scheduler, claim, next_timestamp)
end

@inline function _resolve_mixed_serial_start_claim!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    ::TriggeredAcquisitionStart,
    ::Nothing,
)
    return deactivate_event_generator!(
        prepared.scheduler, state.scheduler, claim)
end

function _begin_mixed_serial_acquisition!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::PlantTimestamp,
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        begin_exposure!(
            acquisition.lifecycle, acquisition_state, timestamp)
        _take_initial_acquisition_snapshot!(
            acquisition.lifecycle,
            acquisition_state,
            timestamp,
            acquisition.rng,
        )
        return nothing
    end
end

function _process_mixed_serial_acquisition_start!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    action::_PlantEventAction,
)
    acquisition = _mixed_serial_event_acquisition(
        prepared, action.owner_slot)
    acquisition_state = _mixed_serial_event_acquisition_state(
        state, action.owner_slot)
    timestamp = claim.key.timestamp
    _require_mixed_serial_event_path_available(
        prepared, state, acquisition, timestamp)
    next_start = _next_periodic_start_timestamp(
        acquisition.start, claim)
    next_start === nothing || begin
        definition = @inbounds prepared.scheduler.definitions[
            Int(claim.slot)]
        next_key = PlantEventKey(
            next_start,
            definition.phase,
            definition.ordinal,
            _next_event_occurrence(claim.key.occurrence),
            _PLANT_TIME_TOKEN,
        )
        _require_forward_event_key(state.scheduler, next_key)
    end
    boundary_timestamp = _initial_acquisition_boundary_timestamp(
        acquisition.lifecycle, timestamp)
    _require_inactive_mixed_serial_generator(
        prepared, state, acquisition.boundary_handle, boundary_timestamp)
    band_open_timestamp = _initial_acquisition_band_open_timestamp(
        acquisition.lifecycle, timestamp)
    band_open_timestamp === nothing ||
        _require_inactive_mixed_serial_generator(
            prepared,
            state,
            acquisition.band_open_handle,
            band_open_timestamp,
        )

    _begin_mixed_serial_acquisition!(
        acquisition, acquisition_state, timestamp)
    boundary_timestamp == _first_acquisition_boundary_timestamp(
        acquisition.lifecycle, acquisition_state) ||
        _mixed_serial_event_error(
            :prepared_binding,
            "acquisition boundary changed after acquisition start",
        )
    band_open_timestamp == _first_acquisition_band_open_timestamp(
        acquisition.lifecycle, acquisition_state) ||
        _mixed_serial_event_error(
            :prepared_binding,
            "rolling row-band schedule changed after acquisition start",
        )

    _resolve_mixed_serial_start_claim!(
        prepared, state, claim, acquisition.start, next_start)
    activate_event_generator!(
        prepared.scheduler,
        state.scheduler,
        acquisition.boundary_handle,
        boundary_timestamp,
    )
    band_open_timestamp === nothing || activate_event_generator!(
        prepared.scheduler,
        state.scheduler,
        acquisition.band_open_handle,
        band_open_timestamp,
    )
    return nothing
end

function _process_mixed_serial_lifecycle_boundary!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::PlantTimestamp,
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        _process_acquisition_lifecycle_boundary!(
            acquisition.lifecycle,
            acquisition_state,
            timestamp,
            acquisition.rng,
        )
    end
end

function _process_mixed_serial_acquisition_boundary!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    action::_PlantEventAction,
)
    acquisition = _mixed_serial_event_acquisition(
        prepared, action.owner_slot)
    acquisition_state = _mixed_serial_event_acquisition_state(
        state, action.owner_slot)
    result = _process_mixed_serial_lifecycle_boundary!(
        acquisition, acquisition_state, claim.key.timestamp)
    if result.disposition == _RescheduleAcquisitionBoundary
        reschedule_event!(
            prepared.scheduler, state.scheduler, claim, result.timestamp)
    else
        _require_inactive_mixed_serial_generator(
            prepared,
            state,
            acquisition.readout_handle,
            result.timestamp,
        )
        deactivate_event_generator!(
            prepared.scheduler, state.scheduler, claim)
        activate_event_generator!(
            prepared.scheduler,
            state.scheduler,
            acquisition.readout_handle,
            result.timestamp,
        )
    end
    return nothing
end

function _open_next_mixed_serial_rolling_band!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::PlantTimestamp,
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        _integrate_event_acquisition_to!(
            acquisition.lifecycle,
            acquisition_state,
            timestamp,
            acquisition.rng,
        )
        open_next_rolling_band!(
            acquisition.lifecycle, acquisition_state, timestamp)
        return next_rolling_band_open_timestamp(
            acquisition.lifecycle, acquisition_state)
    end
end

function _process_mixed_serial_rolling_band_open!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    action::_PlantEventAction,
)
    acquisition = _mixed_serial_event_acquisition(
        prepared, action.owner_slot)
    acquisition_state = _mixed_serial_event_acquisition_state(
        state, action.owner_slot)
    following = _open_next_mixed_serial_rolling_band!(
        acquisition, acquisition_state, claim.key.timestamp)
    if following === nothing
        deactivate_event_generator!(
            prepared.scheduler, state.scheduler, claim)
    else
        reschedule_event!(
            prepared.scheduler, state.scheduler, claim, following)
    end
    return nothing
end

function _complete_mixed_serial_acquisition_readout!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::PlantTimestamp,
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        _complete_event_acquisition_readout!(
            acquisition.product,
            acquisition.lifecycle,
            acquisition_state,
            timestamp,
            acquisition.rng,
        )
        return nothing
    end
end

function _process_mixed_serial_acquisition_readout!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    action::_PlantEventAction,
)
    acquisition = _mixed_serial_event_acquisition(
        prepared, action.owner_slot)
    acquisition_state = _mixed_serial_event_acquisition_state(
        state, action.owner_slot)
    _complete_mixed_serial_acquisition_readout!(
        acquisition, acquisition_state, claim.key.timestamp)
    index = Int(action.owner_slot)
    sequence = _event_product_sequence(acquisition_state)
    previous = @inbounds state.product_sequences[index]
    sequence > previous || _mixed_serial_event_error(
        :product_sequence,
        "acquisition product sequence did not advance",
    )
    @inbounds state.product_sequences[index] = sequence
    @inbounds state.product_ready_timestamps[index] = claim.key.timestamp
    if _event_requires_readiness(acquisition.lifecycle)
        ready_timestamp = _event_readiness_timestamp(
            acquisition.lifecycle, acquisition_state)
        _require_inactive_mixed_serial_generator(
            prepared,
            state,
            acquisition.readiness_handle,
            ready_timestamp,
        )
        deactivate_event_generator!(
            prepared.scheduler, state.scheduler, claim)
        activate_event_generator!(
            prepared.scheduler,
            state.scheduler,
            acquisition.readiness_handle,
            ready_timestamp,
        )
    else
        deactivate_event_generator!(
            prepared.scheduler, state.scheduler, claim)
    end
    return nothing
end

function _mark_mixed_serial_acquisition_ready!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::PlantTimestamp,
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        mark_acquisition_ready!(
            acquisition.lifecycle, acquisition_state, timestamp)
        return nothing
    end
end

function _process_mixed_serial_acquisition_readiness!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    action::_PlantEventAction,
)
    acquisition = _mixed_serial_event_acquisition(
        prepared, action.owner_slot)
    acquisition_state = _mixed_serial_event_acquisition_state(
        state, action.owner_slot)
    _mark_mixed_serial_acquisition_ready!(
        acquisition, acquisition_state, claim.key.timestamp)
    deactivate_event_generator!(
        prepared.scheduler, state.scheduler, claim)
    return nothing
end

@inline function _mixed_serial_command_is_due(
    ::Tuple{},
    ::Tuple{},
    id::CommandEndpointID,
    ::PlantTimestamp,
)
    _mixed_serial_event_error(
        :unknown_command_endpoint,
        "mixed serial command authority contains no endpoint $id",
    )
end

@inline function _mixed_serial_command_is_due(
    lanes::Tuple,
    lane_states::Tuple,
    id::CommandEndpointID,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    if command_endpoint_id(lane.authority_endpoint) == id
        endpoint = lane.authority_endpoint.endpoint
        key = next_command_order_key(endpoint, lane_state.endpoint)
        if key !== nothing && command_scheduled_timestamp(key) <= timestamp
            return true
        end
        silence = next_command_silence_timestamp(
            endpoint, lane_state.endpoint, lane_state.application)
        silence == timestamp || _mixed_serial_event_error(
            :stale_command_event,
            "command generator has no due command or silence transition",
        )
        return false
    end
    return _mixed_serial_command_is_due(
        Base.tail(lanes), Base.tail(lane_states), id, timestamp)
end

function _resolve_mixed_serial_command_claim!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    claim::EventClaim,
    id::CommandEndpointID,
)
    desired = _next_mixed_serial_command_timestamp(prepared, state, id)
    if desired === nothing
        deactivate_event_generator!(
            prepared.scheduler, state.scheduler, claim)
    else
        reschedule_event!(
            prepared.scheduler, state.scheduler, claim, desired)
    end
    return nothing
end

function _process_mixed_serial_command_endpoint!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    claim::EventClaim,
    action::_PlantEventAction,
)
    slot = Int(action.owner_slot)
    1 <= slot <= length(prepared.commands) ||
        _mixed_serial_event_error(
            :invalid_action,
            "command action contains an invalid endpoint slot",
        )
    command = @inbounds prepared.commands[slot]
    timestamp = claim.key.timestamp
    fanout = prepared.execution.fanout
    fanout_state = state.execution.fanout
    fanout_workspace = workspace.execution.fanout
    if _mixed_serial_command_is_due(
        fanout.lanes, fanout_state.lanes, command.id, timestamp)
        _apply_next_command_fanout!(
            fanout,
            fanout_state,
            fanout_workspace,
            command.id,
            timestamp,
        )
    else
        _apply_command_silence_fanout!(
            fanout,
            fanout_state,
            fanout_workspace,
            command.id,
            timestamp,
        )
    end
    _drain_mixed_serial_fanout_dispositions!(
        fanout_workspace.lanes, workspace)
    _resolve_mixed_serial_command_claim!(
        prepared, state, claim, command.id)
    _schedule_all_mixed_serial_command_endpoints!(prepared, state)
    return nothing
end

function _triggered_mixed_serial_acquisition_slot_or_zero(
    prepared::PreparedMixedResourcePlantEventLoop,
    consumer::TriggerConsumerID,
)
    @inbounds for slot in eachindex(prepared.acquisitions)
        start = prepared.acquisitions[slot].start
        _start_matches_consumer(start, consumer) && return UInt32(slot)
    end
    return zero(UInt32)
end

function _process_mixed_serial_trigger_topology!(
    prepared::PreparedMixedResourcePlantEventLoop{
        <:Any,<:Any,<:PreparedTriggerTopology},
    state::MixedResourcePlantEventLoopState{<:Any,<:Any,<:TriggerTopologyState},
    workspace::MixedResourcePlantEventLoopWorkspace{
        <:Any,<:Any,<:TriggerTopologyWorkspace},
    claim::EventClaim,
)
    topology = prepared.trigger_topology
    trigger_state = state.trigger
    source = next_trigger_source(topology, trigger_state)
    delivery = next_trigger_delivery(topology, trigger_state)
    source_due = delivery === nothing ||
        realized_trigger_source_timestamp(source) <=
            delivered_trigger_edge(delivery).timestamp
    activated_slot = zero(UInt32)
    activation_timestamp = zero(PlantTimestamp)
    if source_due
        realized_trigger_source_timestamp(source) == claim.key.timestamp ||
            _mixed_serial_event_error(
                :trigger_schedule,
                "trigger source does not match its scheduler claim",
            )
        realize_next_trigger_source!(
            workspace.trigger, topology, trigger_state)
    else
        delivered = delivered_trigger_edge(delivery)
        delivered.timestamp == claim.key.timestamp ||
            _mixed_serial_event_error(
                :trigger_schedule,
                "trigger delivery does not match its scheduler claim",
            )
        activated_slot =
            _triggered_mixed_serial_acquisition_slot_or_zero(
                prepared, trigger_delivery_consumer(delivery))
        iszero(activated_slot) && _mixed_serial_event_error(
            :trigger_binding,
            "delivered trigger consumer " *
            "$(trigger_delivery_consumer(delivery)) has no acquisition",
        )
        acquisition = _mixed_serial_event_acquisition(
            prepared, activated_slot)
        _require_inactive_mixed_serial_generator(
            prepared,
            state,
            acquisition.start_handle,
            delivered.timestamp,
        )
        pop_next_trigger_delivery!(
            workspace.delivery, topology, trigger_state) ||
            _mixed_serial_event_error(
                :trigger_schedule,
                "due trigger delivery disappeared before removal",
            )
        activation_timestamp = delivered.timestamp
    end
    next_timestamp = _next_trigger_action_timestamp(topology, trigger_state)
    reschedule_event!(
        prepared.scheduler, state.scheduler, claim, next_timestamp)
    if !iszero(activated_slot)
        acquisition = _mixed_serial_event_acquisition(
            prepared, activated_slot)
        activate_event_generator!(
            prepared.scheduler,
            state.scheduler,
            acquisition.start_handle,
            activation_timestamp,
        )
    end
    return nothing
end

function _process_mixed_serial_trigger_topology!(
    ::PreparedMixedResourcePlantEventLoop{
        <:Any,<:Any,<:_NoPreparedTriggerTopology},
    ::MixedResourcePlantEventLoopState{
        <:Any,<:Any,<:_NoTriggerTopologyState},
    ::MixedResourcePlantEventLoopWorkspace{
        <:Any,<:Any,<:_NoTriggerTopologyWorkspace},
    ::EventClaim,
)
    _mixed_serial_event_error(
        :invalid_action,
        "trigger action exists without a prepared trigger topology",
    )
end

function _preflight_mixed_serial_acquisition_integration(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::Base.RefValue{PlantTimestamp},
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _require_event_lifecycle_binding(
            acquisition.lifecycle, acquisition_state)
        _preflight_event_integration_to(
            acquisition.lifecycle, acquisition_state, timestamp[])
        return nothing
    end
end

function _integrate_mixed_serial_acquisition_to!(
    acquisition::_AbstractPreparedMixedSerialEventAcquisition,
    acquisition_state::_AcquisitionEventLifecycleState,
    timestamp::Base.RefValue{PlantTimestamp},
)
    return _with_completed_prepared_device_execution_context(
        acquisition.resources.context) do
        _integrate_event_acquisition_to!(
            acquisition.lifecycle,
            acquisition_state,
            timestamp[],
            acquisition.rng,
        )
        return nothing
    end
end

function _mark_due_mixed_serial_paths!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    due = workspace.execution.due_paths
    fill!(due, false)
    scan_due_events!(
        workspace.scheduler, prepared.scheduler, state.scheduler)
    count = workspace.scheduler.due_count
    @inbounds for due_index in 1:count
        scheduler_slot = Int(workspace.scheduler.due_slots[due_index])
        definition = prepared.scheduler.definitions[scheduler_slot]
        definition.phase == OpticalSamplePhase || continue
        state.scheduler.cursors[scheduler_slot].next_timestamp == timestamp ||
            continue
        action = prepared.actions[scheduler_slot]
        due[Int(action.owner_slot)] = true
    end
    any(due) || _mixed_serial_event_error(
        :invalid_action,
        "optical sample phase has no due optical path",
    )
    return nothing
end

function _preflight_due_mixed_serial_acquisitions!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    due::Memory{Bool},
    timestamp::Base.RefValue{PlantTimestamp},
)
    @inbounds for path_slot in eachindex(prepared.paths)
        due[path_slot] || continue
        path = prepared.paths[path_slot]
        for acquisition_slot in path.acquisition_slots
            index = Int(acquisition_slot)
            _preflight_mixed_serial_acquisition_integration(
                prepared.acquisitions[index],
                state.acquisitions[index],
                timestamp,
            )
        end
    end
    return nothing
end

function _integrate_due_mixed_serial_acquisitions!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    due::Memory{Bool},
    timestamp::Base.RefValue{PlantTimestamp},
)
    @inbounds for path_slot in eachindex(prepared.paths)
        due[path_slot] || continue
        path = prepared.paths[path_slot]
        for acquisition_slot in path.acquisition_slots
            index = Int(acquisition_slot)
            _integrate_mixed_serial_acquisition_to!(
                prepared.acquisitions[index],
                state.acquisitions[index],
                timestamp,
            )
        end
        state.path_sampled[path_slot] = true
    end
    return nothing
end

function _resolve_due_mixed_serial_path_claims!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    while true
        count = scan_due_events!(
            workspace.scheduler, prepared.scheduler, state.scheduler)
        iszero(count) && return nothing
        workspace.scheduler.due_timestamp == timestamp || return nothing
        key = due_event_key(
            workspace.scheduler, prepared.scheduler, state.scheduler, 1)
        key.phase == OpticalSamplePhase || return nothing
        claim = claim_next_event!(
            workspace.scheduler, prepared.scheduler, state.scheduler)
        claim === nothing && _mixed_serial_event_error(
            :invalid_action,
            "due optical path disappeared before claim",
        )
        action = _mixed_serial_event_action(prepared, claim)
        path = _mixed_serial_event_path(prepared, action.owner_slot)
        reschedule_periodic_event!(
            prepared.scheduler,
            state.scheduler,
            claim,
            path.schedule;
            origin=path.origin,
        )
    end
end

function _execute_due_mixed_serial_paths!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    _mark_due_mixed_serial_paths!(
        prepared, state, workspace, timestamp)
    due = workspace.execution.due_paths
    workspace.event_timestamp[] = timestamp
    _preflight_due_mixed_serial_acquisitions!(
        prepared, state, due, workspace.event_timestamp)
    _preflight_selected_mixed_serial_paths!(
        prepared.execution,
        state.execution,
        workspace.execution,
        timestamp,
    )
    _integrate_due_mixed_serial_acquisitions!(
        prepared, state, due, workspace.event_timestamp)
    _execute_preflighted_mixed_serial_paths!(
        prepared.execution,
        state.execution,
        workspace.execution,
        timestamp,
    )
    _resolve_due_mixed_serial_path_claims!(
        prepared, state, workspace, timestamp)
    return nothing
end

function _process_mixed_serial_ordinary_event!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    claim::EventClaim,
    action::_PlantEventAction,
)
    kind = action.kind
    kind == _TriggerTopologyAction &&
        return _process_mixed_serial_trigger_topology!(
            prepared, state, workspace, claim)
    kind == _CommandEndpointAction &&
        return _process_mixed_serial_command_endpoint!(
            prepared, state, workspace, claim, action)
    kind == _AcquisitionBoundaryAction &&
        return _process_mixed_serial_acquisition_boundary!(
            prepared, state, claim, action)
    kind == _AcquisitionStartAction &&
        return _process_mixed_serial_acquisition_start!(
            prepared, state, claim, action)
    kind == _RollingBandOpenAction &&
        return _process_mixed_serial_rolling_band_open!(
            prepared, state, claim, action)
    kind == _AcquisitionReadoutAction &&
        return _process_mixed_serial_acquisition_readout!(
            prepared, state, claim, action)
    kind == _AcquisitionReadinessAction &&
        return _process_mixed_serial_acquisition_readiness!(
            prepared, state, claim, action)
    _mixed_serial_event_error(
        :invalid_action,
        "prepared mixed serial event has an unknown action kind",
    )
end

function next_plant_event_timestamp(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    _require_mixed_serial_event_ownership(prepared, state, workspace)
    _require_mixed_serial_event_idle(state)
    count = scan_due_events!(
        workspace.scheduler, prepared.scheduler, state.scheduler)
    iszero(count) && return nothing
    return workspace.scheduler.due_timestamp
end

function step_plant_events!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
)
    _require_mixed_serial_event_ownership(prepared, state, workspace)
    _require_mixed_serial_event_idle(state)
    count = scan_due_events!(
        workspace.scheduler, prepared.scheduler, state.scheduler)
    iszero(count) && return nothing
    timestamp = workspace.scheduler.due_timestamp
    state.phase = _MixedSerialEventLoopProcessing
    try
        while true
            due_count = scan_due_events!(
                workspace.scheduler, prepared.scheduler, state.scheduler)
            iszero(due_count) && break
            workspace.scheduler.due_timestamp == timestamp || break
            key = due_event_key(
                workspace.scheduler,
                prepared.scheduler,
                state.scheduler,
                1,
            )
            if key.phase == OpticalSamplePhase
                _execute_due_mixed_serial_paths!(
                    prepared, state, workspace, timestamp)
                continue
            end
            claim = claim_next_event!(
                workspace.scheduler, prepared.scheduler, state.scheduler)
            claim === nothing && _mixed_serial_event_error(
                :invalid_action,
                "due mixed serial event disappeared before claim",
            )
            action = _mixed_serial_event_action(prepared, claim)
            _process_mixed_serial_ordinary_event!(
                prepared, state, workspace, claim, action)
        end
    catch
        state.phase = _MixedSerialEventLoopFailed
        rethrow()
    end
    state.phase = _MixedSerialEventLoopIdle
    return timestamp
end

function run_plant_events_until!(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    workspace::MixedResourcePlantEventLoopWorkspace,
    stop::PlantTimestamp;
    max_timestamps::Integer=typemax(Int),
)
    limit = _checked_event_step_limit(max_timestamps)
    count = 0
    while count < limit
        timestamp = next_plant_event_timestamp(
            prepared, state, workspace)
        timestamp === nothing && break
        timestamp <= stop || break
        step_plant_events!(prepared, state, workspace)
        count += 1
    end
    return count
end

function _mixed_serial_event_acquisition_slot(
    prepared::PreparedMixedResourcePlantEventLoop,
    id::AcquisitionID,
)
    @inbounds for slot in eachindex(prepared.acquisitions)
        prepared.acquisitions[slot].id == id && return slot
    end
    _mixed_serial_event_error(
        :unknown_acquisition,
        "mixed serial event loop contains no acquisition $id",
    )
end

function effective_command(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    id,
)
    prepared.binding === state.binding || _mixed_serial_event_error(
        :foreign_state,
        "mixed serial event-loop state belongs to another preparation",
    )
    fanout = prepared.execution.fanout
    fanout_state = state.execution.fanout
    application = command_authority_application_state(
        fanout.authority, fanout_state.authority, id)
    return effective_command(application)
end

function acquisition_products(
    prepared::PreparedMixedResourcePlantEventLoop,
    id,
)
    slot = _mixed_serial_event_acquisition_slot(
        prepared, _as_acquisition_id(id))
    return acquisition_products(
        @inbounds prepared.acquisitions[slot].resources)
end

function acquisition_product_sequence(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    id,
)
    prepared.binding === state.binding || _mixed_serial_event_error(
        :foreign_state,
        "mixed serial event-loop state belongs to another preparation",
    )
    slot = _mixed_serial_event_acquisition_slot(
        prepared, _as_acquisition_id(id))
    return @inbounds state.product_sequences[slot]
end

function acquisition_product_ready_timestamp(
    prepared::PreparedMixedResourcePlantEventLoop,
    state::MixedResourcePlantEventLoopState,
    id,
)
    prepared.binding === state.binding || _mixed_serial_event_error(
        :foreign_state,
        "mixed serial event-loop state belongs to another preparation",
    )
    slot = _mixed_serial_event_acquisition_slot(
        prepared, _as_acquisition_id(id))
    @inbounds iszero(state.product_sequences[slot]) && return nothing
    return @inbounds state.product_ready_timestamps[slot]
end
