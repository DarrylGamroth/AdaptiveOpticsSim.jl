#
# Deterministic mixed-target serial optical execution
#
# This private coordinator is the correctness oracle for one atmosphere
# authority, one command authority, and target-local optical-path replicas. It
# intentionally owns no scheduler, task, queue, pacing policy, or transport.
# Topology-sized registries remain cardinality-neutral; the small ordered stage
# lists within one path are sealed as concrete tuples during preparation.
#

mutable struct _PreparedMixedSerialExecutionBinding end

struct _PreparedMixedSerialControllableOpticApplication{O,C}
    owner::O
    coupling::C
end

struct _PreparedMixedSerialCommandDependency{D}
    lane_slot::UInt32
    destination::D
end

abstract type _AbstractPreparedMixedSerialPath end

struct _PreparedMixedSerialPath{P,R,G,S,O<:Tuple,D<:Tuple} <:
       _AbstractPreparedMixedSerialPath
    id::OpticalPathID
    path::P
    route::R
    rngs::G
    sampled_aberrations::S
    controllable_optics::O
    command_dependencies::D
end

struct _PreparedMixedSerialExecution{B,P,F}
    binding::B
    partitions::P
    fanout::F
    path_ids::_FixedPlantRegistry{OpticalPathID}
    paths::Memory{_AbstractPreparedMixedSerialPath}
end

@enum _MixedSerialExecutionPhase::UInt8 begin
    _MixedSerialExecutionIdle = 0x01
    _MixedSerialExecutionMaterializing = 0x02
    _MixedSerialExecutionExecuting = 0x03
    _MixedSerialExecutionFailed = 0x04
end

mutable struct _MixedSerialExecutionState{B,F}
    binding::B
    fanout::F
    phase::_MixedSerialExecutionPhase
    last_timestamp::PlantTimestamp
    has_timestamp::Bool
end

abstract type _AbstractMixedSerialPathWorkspace end

struct _MixedSerialPathWorkspace{Q} <: _AbstractMixedSerialPathWorkspace
    publication::Base.RefValue{Q}
end

mutable struct _MixedSerialBatchWorkspace{F,S,E}
    fanout::F
    fanout_state::S
    next_epoch_sequence::UInt64
    timestamp::PlantTimestamp
    epoch::E
end

struct _MixedSerialExecutionWorkspace{B,F,C}
    binding::B
    fanout::F
    batch::C
    due_paths::Memory{Bool}
    paths::Memory{_AbstractMixedSerialPathWorkspace}
end

@noinline function _mixed_serial_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(:mixed_serial_execution, reason,
        String(message)))
end

@noinline function _mixed_serial_execution_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantScheduleError(:mixed_serial_execution, reason,
        String(message)))
end

function _mixed_serial_partition_path_slot(
    partition::PreparedTargetPartition,
    id::OpticalPathID,
)
    found = 0
    @inbounds for slot in eachindex(prepared_paths(partition))
        path_id(prepared_paths(partition)[slot].definition) == id || continue
        iszero(found) || _mixed_serial_preparation_error(
            :duplicate_path,
            "target partition contains more than one optical path $id",
        )
        found = slot
    end
    iszero(found) && _mixed_serial_preparation_error(
        :missing_path,
        "target partition does not contain optical path $id",
    )
    return found
end

function _mixed_serial_prepared_optics(
    owners::AbstractVector,
)
    optics = Memory{PreparedTargetLocalControllableOptic}(
        undef, length(owners))
    @inbounds for slot in eachindex(owners)
        optics[slot] = prepared_target_local_controllable_optic(owners[slot])
    end
    return optics
end

@inline function _require_mixed_serial_pupil_surface_role(
    ::PupilSurfaceExecutionRole,
    ::ControllableOpticID,
    ::OpticalPathID,
)
    return nothing
end

function _require_mixed_serial_pupil_surface_role(
    role::AbstractControllableOpticExecutionRole,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    _mixed_serial_preparation_error(
        :unsupported_optic_execution_role,
        "mixed serial execution does not yet prepare $(typeof(role)) " *
        "controllable optic $optic on optical path $path",
    )
end

@inline function _same_mixed_serial_endpoint_owner(
    left::TargetLocalCommandEndpointOwner,
    right::TargetLocalCommandEndpointOwner,
)
    return left.endpoint.binding === right.endpoint.binding &&
        left.state.binding === right.state.binding &&
        left.optic.prepared.binding === right.optic.prepared.binding &&
        left.optic.state.binding === right.optic.state.binding &&
        left.optic.workspace.binding === right.optic.workspace.binding
end

function _mixed_serial_command_dependency(
    fanout::_PreparedCommandFanout,
    endpoint::TargetLocalCommandEndpointOwner,
)
    endpoint_id = command_endpoint_id(endpoint)
    found_lane = 0
    found_destination = nothing
    for (lane_slot, lane) in enumerate(fanout.lanes)
        command_endpoint_id(lane.authority_endpoint) == endpoint_id ||
            continue
        iszero(found_lane) || _mixed_serial_preparation_error(
            :duplicate_command_endpoint,
            "command fanout contains more than one lane for $endpoint_id",
        )
        for route in lane.routes
            destination = _effective_command_route_destination(route)
            _same_mixed_serial_endpoint_owner(destination, endpoint) ||
                continue
            isnothing(found_destination) ||
                _mixed_serial_preparation_error(
                    :duplicate_command_replica,
                    "command fanout contains more than one route to " *
                    "$(compute_device(endpoint)) for $endpoint_id",
                )
            found_destination = destination
        end
        found_lane = lane_slot
    end
    iszero(found_lane) && _mixed_serial_preparation_error(
        :missing_command_endpoint,
        "command fanout contains no lane for $endpoint_id",
    )
    isnothing(found_destination) && _mixed_serial_preparation_error(
        :missing_command_replica,
        "command fanout contains no exact route to $(compute_device(endpoint)) " *
        "for $endpoint_id",
    )
    found_lane <= typemax(UInt32) || _mixed_serial_preparation_error(
        :capacity,
        "command-fanout lane count exceeds UInt32 capacity",
    )
    return _PreparedMixedSerialCommandDependency(
        UInt32(found_lane), found_destination)
end

function _append_mixed_serial_command_dependencies!(
    dependencies::Vector{Any},
    fanout::_PreparedCommandFanout,
    owner::TargetLocalControllableOpticOwner,
)
    for endpoint in owner.prepared.endpoints
        endpoint_owner = target_local_command_endpoint_owner(
            owner, command_endpoint_id(endpoint))
        push!(dependencies,
            _mixed_serial_command_dependency(fanout, endpoint_owner))
    end
    return nothing
end

function _prepare_mixed_serial_controllable_optics(
    partition::PreparedTargetPartition,
    path::PreparedTargetLocalPathResources,
    fanout::_PreparedCommandFanout,
)
    owners = target_local_controllable_optic_owners(partition)
    optics = _mixed_serial_prepared_optics(owners)
    bindings = _prepare_controllable_optic_path_bindings(
        optics, prepared_paths(partition))
    id = path_id(path.definition)
    applications = Any[]
    dependencies = Any[]
    binding_range = prepared_controllable_optic_binding_range(bindings, id)
    sizehint!(applications, length(binding_range))
    @inbounds for binding in binding_range
        slot = prepared_controllable_optic_slot(bindings, binding)
        owner = owners[slot]
        implementation = owner.prepared.implementation
        _require_mixed_serial_pupil_surface_role(
            controllable_optic_execution_role(implementation),
            controllable_optic_id(owner),
            id,
        )
        coupling = prepare_controllable_optic_path_coupling(
            owner.prepared, path)
        push!(applications,
            _PreparedMixedSerialControllableOpticApplication(
                owner, coupling))
        _append_mixed_serial_command_dependencies!(
            dependencies, fanout, owner)
    end
    return Tuple(applications), Tuple(dependencies)
end

function _prepare_mixed_serial_path(
    partitions::PreparedPlantPartitions,
    fanout::_PreparedCommandFanout,
    id::OpticalPathID,
)
    target = partition_target(getfield(partitions, :assignment), id)
    partition = prepared_partition(partitions, target)
    slot = _mixed_serial_partition_path_slot(partition, id)
    path = @inbounds prepared_paths(partition)[slot]
    rngs = @inbounds getfield(partition, :rngs).paths[slot]
    _require_rng_owner_binding(rngs, path)

    sampled = prepared_sampled_aberrations(partition)
    sampled_bindings = _prepare_sampled_aberration_path_bindings(
        sampled, prepared_paths(partition))
    sampled_plan = _prepare_sampled_aberration_path_plan(
        sampled, sampled_bindings, id)
    applications, dependencies =
        _prepare_mixed_serial_controllable_optics(
            partition, path, fanout)
    route = prepare_pupil_opd_publication_route(partitions, id)
    pupil_opd_publication_path_id(route) == id ||
        _mixed_serial_preparation_error(
            :prepared_binding,
            "pupil-OPD publication route does not retain optical path $id",
        )
    return _PreparedMixedSerialPath(
        id,
        path,
        route,
        rngs,
        sampled_plan,
        applications,
        dependencies,
    )
end

function _prepare_mixed_serial_execution(
    partitions::PreparedPlantPartitions,
)
    ids = _canonical_partition_path_ids(plant_definition(partitions))
    isempty(ids) && _mixed_serial_preparation_error(
        :empty_topology,
        "mixed serial execution requires at least one optical path",
    )
    fanout = _prepare_command_fanout(partitions)
    paths = Memory{_AbstractPreparedMixedSerialPath}(undef, length(ids))
    @inbounds for slot in eachindex(ids)
        paths[slot] = _prepare_mixed_serial_path(
            partitions, fanout, ids[slot])
    end
    return _PreparedMixedSerialExecution(
        _PreparedMixedSerialExecutionBinding(),
        partitions,
        fanout,
        ids,
        paths,
    )
end

function _prepare_mixed_serial_execution_state(
    prepared::_PreparedMixedSerialExecution,
)
    return _MixedSerialExecutionState(
        prepared.binding,
        _prepare_command_fanout_state(prepared.fanout),
        _MixedSerialExecutionIdle,
        zero(PlantTimestamp),
        false,
    )
end

function _prepare_mixed_serial_execution_workspace(
    prepared::_PreparedMixedSerialExecution,
    state::_MixedSerialExecutionState,
)
    prepared.binding === state.binding ||
        _mixed_serial_preparation_error(
            :foreign_state,
            "mixed serial execution state belongs to another preparation",
        )
    path_workspaces = Memory{_AbstractMixedSerialPathWorkspace}(
        undef, length(prepared.paths))
    @inbounds for slot in eachindex(prepared.paths)
        output = prepare_pupil_opd_publication_output(
            prepared.paths[slot].route)
        path_workspaces[slot] = _MixedSerialPathWorkspace(output)
    end
    atmosphere = prepared_atmosphere(
        prepared_atmosphere_authority(prepared.partitions))
    timeline = atmosphere_timeline(atmosphere)
    epoch = AtmosphereEpoch(
        atmosphere_identity(atmosphere),
        timeline.model_time,
        timeline.sequence,
    )
    return _MixedSerialExecutionWorkspace(
        prepared.binding,
        _prepare_command_fanout_workspace(
            prepared.fanout, state.fanout),
        _MixedSerialBatchWorkspace(
            prepared.fanout,
            state.fanout,
            UInt64(0),
            zero(PlantTimestamp),
            epoch,
        ),
        Memory{Bool}(undef, length(prepared.paths)),
        path_workspaces,
    )
end

@inline _mixed_serial_path_id(id::OpticalPathID) = id
@inline _mixed_serial_path_id(id::Symbol) = OpticalPathID(id)

function _mixed_serial_path_id(id)
    _mixed_serial_execution_error(
        :invalid_path_id,
        "due optical path identity must be an OpticalPathID or Symbol; got " *
        "$(typeof(id))",
    )
end

function _mixed_serial_path_slot(
    path_ids::_FixedPlantRegistry{OpticalPathID},
    id::OpticalPathID,
)
    @inbounds for slot in eachindex(path_ids)
        path_ids[slot] == id && return slot
    end
    _mixed_serial_execution_error(
        :unknown_path,
        "mixed serial execution contains no optical path $id",
    )
end

function _select_mixed_serial_due_paths!(
    due::Memory{Bool},
    path_ids::_FixedPlantRegistry{OpticalPathID},
    ids::Union{Tuple,AbstractVector},
)
    fill!(due, false)
    count = 0
    for value in ids
        id = _mixed_serial_path_id(value)
        slot = _mixed_serial_path_slot(path_ids, id)
        due[slot] && _mixed_serial_execution_error(
            :duplicate_due_path,
            "optical path $id occurs more than once in the due set",
        )
        due[slot] = true
        count += 1
    end
    iszero(count) && _mixed_serial_execution_error(
        :empty_due_paths,
        "mixed serial execution requires at least one due optical path",
    )
    return nothing
end

function _select_mixed_serial_due_paths!(
    ::Memory{Bool},
    ::_FixedPlantRegistry{OpticalPathID},
    ids,
)
    _mixed_serial_execution_error(
        :invalid_due_paths,
        "due optical paths must be a Tuple or AbstractVector; got " *
        "$(typeof(ids))",
    )
end

@inline function _require_mixed_serial_binding(
    prepared::_PreparedMixedSerialExecution,
    state::_MixedSerialExecutionState,
    workspace::_MixedSerialExecutionWorkspace,
)
    prepared.binding === state.binding === workspace.binding ||
        _mixed_serial_execution_error(
            :foreign_owner,
            "mixed serial preparation, state, and workspace belong to " *
            "different owners",
        )
    length(prepared.path_ids) == length(prepared.paths) ==
        length(workspace.paths) == length(workspace.due_paths) ||
        _mixed_serial_execution_error(
            :prepared_binding,
            "mixed serial workspace topology changed after preparation",
        )
    workspace.batch.fanout.binding === prepared.fanout.binding &&
        workspace.batch.fanout_state.binding === state.fanout.binding ||
        _mixed_serial_execution_error(
            :prepared_binding,
            "mixed serial batch workspace no longer matches its fanout",
        )
    state.phase == _MixedSerialExecutionIdle ||
        _mixed_serial_execution_error(
            state.phase == _MixedSerialExecutionFailed ?
                :executor_failed : :executor_busy,
            "mixed serial execution is not idle ($(state.phase))",
        )
    _require_command_fanout_binding(
        prepared.fanout, state.fanout, workspace.fanout)
    _require_command_fanout_workspace_idle(workspace.fanout.lanes)
    return nothing
end

@inline function _require_mixed_serial_timestamp(
    state::_MixedSerialExecutionState,
    timestamp::PlantTimestamp,
)
    state.has_timestamp && timestamp <= state.last_timestamp &&
        _mixed_serial_execution_error(
            :nonmonotonic_batch_time,
            "optical-path batch timestamp $timestamp must follow " *
            "$(state.last_timestamp)",
        )
    return nothing
end

function _preflight_mixed_serial_atmosphere_time(
    atmosphere::AbstractTimedAtmosphere,
    timestamp::PlantTimestamp,
)
    timeline = atmosphere_timeline(atmosphere)
    T = typeof(timeline.model_time)
    target = plant_time_seconds(timestamp, T)
    timeline.initialized && target < timeline.model_time &&
        _mixed_serial_execution_error(
            :atmosphere_time_regression,
            "due optical sample precedes the current atmosphere epoch",
        )
    advances = !timeline.initialized || target != timeline.model_time
    advances && timeline.sequence == typemax(UInt64) &&
        _mixed_serial_execution_error(
            :atmosphere_sequence_exhausted,
            "atmosphere epoch sequence is exhausted",
        )
    return target, advances ? timeline.sequence + one(UInt64) :
        timeline.sequence
end

function _preflight_mixed_serial_path(
    path::_PreparedMixedSerialPath,
    batch::_MixedSerialBatchWorkspace,
)
    route = path.route
    route.state.phase == _PupilOPDRouteIdle ||
        _mixed_serial_execution_error(
            :path_input_route_busy,
            "pupil-OPD route for $(path.id) is not idle",
        )
    route.state.generation != typemax(UInt64) ||
        _mixed_serial_execution_error(
            :path_input_generation_exhausted,
            "pupil-OPD route generation for $(path.id) is exhausted",
        )
    batch.next_epoch_sequence > route.state.last_epoch_sequence ||
        _mixed_serial_execution_error(
            :path_input_epoch_not_new,
            "next atmosphere epoch is not newer than the last publication " *
            "for $(path.id)",
        )
    _validate_pupil_opd_route_binding(route)
    _require_rng_owner_binding(path.rngs, path.path)
    _require_mixed_serial_command_dependencies!(
        path.command_dependencies,
        batch.fanout.lanes,
        batch.fanout_state.lanes,
        batch.fanout_state.authority.publication_sequences,
        batch.timestamp,
    )
    return nothing
end

@noinline function _require_mixed_serial_command_dependency(
    lane::_PreparedCommandFanoutLane,
    lane_state::_CommandFanoutLaneState,
    dependency::_PreparedMixedSerialCommandDependency,
    publication_sequence::UInt64,
    timestamp::PlantTimestamp,
)
    destination = dependency.destination
    endpoint = lane.authority_endpoint.endpoint
    command_endpoint_id(endpoint) == command_endpoint_id(destination) ||
        _mixed_serial_execution_error(
            :command_dependency_binding,
            "prepared command dependency no longer matches its fanout lane",
        )
    _require_command_endpoint_binding(endpoint, lane_state.endpoint)
    _require_command_application_binding(
        endpoint, lane_state.endpoint, lane_state.application)
    _require_operational_command_endpoint(lane_state.endpoint)
    _require_idle_command_endpoint(lane_state.endpoint)
    destination.state.has_staged_publication &&
        _mixed_serial_execution_error(
            :command_replica_staged,
            "target-local command endpoint $(command_endpoint_id(destination)) " *
            "retains a staged publication",
        )
    destination.optic.state.has_staged_publication &&
        _mixed_serial_execution_error(
            :controllable_optic_staged,
            "target-local controllable optic " *
            "$(controllable_optic_id(destination)) retains a staged command",
        )

    key = next_command_order_key(endpoint, lane_state.endpoint)
    key !== nothing && command_scheduled_timestamp(key) <= timestamp &&
        _mixed_serial_execution_error(
            :command_not_applied,
            "command endpoint $(command_endpoint_id(destination)) has an " *
            "unapplied command due no later than optical sample $timestamp",
        )
    silence = next_command_silence_timestamp(
        endpoint, lane_state.endpoint, lane_state.application)
    silence !== nothing && silence <= timestamp &&
        _mixed_serial_execution_error(
            :command_silence_not_applied,
            "command endpoint $(command_endpoint_id(destination)) has an " *
            "unapplied silence transition due no later than optical sample " *
            "$timestamp",
        )

    local_state = destination.state
    if iszero(publication_sequence)
        local_state.has_publication && _mixed_serial_execution_error(
            :command_replica_ahead,
            "target-local command endpoint $(command_endpoint_id(destination)) " *
            "has a publication before its authority",
        )
        iszero(local_state.last_sequence) ||
            _mixed_serial_execution_error(
                :command_replica_sequence,
                "target-local command endpoint " *
                "$(command_endpoint_id(destination)) has a nonzero baseline " *
                "publication sequence",
            )
        return nothing
    end

    local_state.has_publication || _mixed_serial_execution_error(
        :command_replica_missing,
        "target-local command endpoint $(command_endpoint_id(destination)) " *
        "has not received authority publication $publication_sequence",
    )
    local_state.last_sequence == publication_sequence ||
        _mixed_serial_execution_error(
            :command_replica_sequence,
            "target-local command endpoint $(command_endpoint_id(destination)) " *
            "is at sequence $(local_state.last_sequence); authority is at " *
            "$publication_sequence",
        )
    applied = last_command_application_timestamp(lane_state.application)
    local_state.last_timestamp == applied ||
        _mixed_serial_execution_error(
            :command_replica_timestamp,
            "target-local command endpoint $(command_endpoint_id(destination)) " *
            "does not match its authority application timestamp",
        )
    applied <= timestamp || _mixed_serial_execution_error(
        :command_from_future,
        "command endpoint $(command_endpoint_id(destination)) was applied " *
        "after optical sample $timestamp",
    )
    return nothing
end

@inline _require_mixed_serial_command_dependencies!(
    ::Tuple{}, ::Tuple, ::Tuple, ::Memory{UInt64}, ::PlantTimestamp) = nothing

@inline function _require_mixed_serial_command_dependencies!(
    dependencies::Tuple,
    lanes::Tuple,
    lane_states::Tuple,
    sequences::Memory{UInt64},
    timestamp::PlantTimestamp,
)
    dependency = first(dependencies)
    slot = Int(dependency.lane_slot)
    @inbounds _require_mixed_serial_command_dependency(
        lanes[slot],
        lane_states[slot],
        dependency,
        sequences[slot],
        timestamp,
    )
    return _require_mixed_serial_command_dependencies!(
        Base.tail(dependencies), lanes, lane_states, sequences, timestamp)
end

function _preflight_mixed_serial_due_paths(
    prepared::_PreparedMixedSerialExecution,
    due::Memory{Bool},
    batch::_MixedSerialBatchWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        due[slot] || continue
        _preflight_mixed_serial_path(
            prepared.paths[slot], batch)::Nothing
    end
    return nothing
end

function _materialize_mixed_serial_path_input!(
    path::_PreparedMixedSerialPath,
    workspace::_MixedSerialPathWorkspace,
    batch::_MixedSerialBatchWorkspace,
)
    status = materialize_pupil_opd_publication!(
        workspace.publication, path.route, batch.epoch, batch.timestamp)
    status == PupilOPDPublicationSucceeded ||
        _mixed_serial_execution_error(
            :path_input_materialization,
            "pupil-OPD materialization for $(path.id) failed ($status)",
        )
    return nothing
end

function _submit_mixed_serial_path_input!(
    path::_PreparedMixedSerialPath,
    workspace::_MixedSerialPathWorkspace,
)
    status = submit_pupil_opd_publication!(
        path.route, workspace.publication[])
    status == PupilOPDPublicationSucceeded ||
        _mixed_serial_execution_error(
            :path_input_submission,
            "pupil-OPD submission for $(path.id) failed ($status)",
        )
    return nothing
end

function _complete_mixed_serial_path_input!(
    path::_PreparedMixedSerialPath,
    workspace::_MixedSerialPathWorkspace,
)
    status = _complete_pupil_opd_publication!(
        path.route, workspace.publication[])
    status == PupilOPDPublicationSucceeded ||
        _mixed_serial_execution_error(
            :path_input_completion,
            "pupil-OPD completion for $(path.id) failed ($status)",
        )
    return nothing
end

function _apply_mixed_serial_path_input!(
    path::_PreparedMixedSerialPath,
    workspace::_MixedSerialPathWorkspace,
)
    status = apply_pupil_opd_publication!(
        path.route, workspace.publication[])
    status == PupilOPDPublicationSucceeded ||
        _mixed_serial_execution_error(
            :path_input_application,
            "pupil-OPD application for $(path.id) failed ($status)",
        )
    return nothing
end

@inline _apply_mixed_serial_controllable_optics!(input, ::Tuple{}) = input

@inline function _apply_mixed_serial_controllable_optics!(
    input,
    applications::Tuple,
)
    application = first(applications)
    owner = application.owner
    apply_controllable_optic_surface!(
        input,
        owner.prepared.implementation,
        owner.state.physical,
        application.coupling,
    )
    return _apply_mixed_serial_controllable_optics!(
        input, Base.tail(applications))
end

function _execute_mixed_serial_path!(
    path::PreparedTargetLocalPathResources,
    rngs::PreparedOwnerRNGs,
    sampled::_PreparedSampledAberrationPathPlan,
    optics::Tuple,
)
    _require_rng_owner_binding(rngs, path)
    return _with_completed_prepared_device_execution_context(
        path.context) do
        _apply_sampled_aberration_path_plan!(path.input, sampled)
        _apply_mixed_serial_controllable_optics!(path.input, optics)
        execute_path_rngs!(path.result, path.input, path.execution, rngs)
        return nothing
    end
end

function _execute_mixed_serial_path!(
    path::_PreparedMixedSerialPath,
)
    return _execute_mixed_serial_path!(
        path.path,
        path.rngs,
        path.sampled_aberrations,
        path.controllable_optics,
    )
end

function _reclaim_mixed_serial_path_input!(
    path::_PreparedMixedSerialPath,
    workspace::_MixedSerialPathWorkspace,
)
    status = reclaim_pupil_opd_publication!(
        path.route, workspace.publication[])
    status == PupilOPDPublicationSucceeded ||
        _mixed_serial_execution_error(
            :path_input_reclamation,
            "pupil-OPD reclamation for $(path.id) failed ($status)",
        )
    return nothing
end

function _materialize_mixed_serial_due_paths!(
    prepared::_PreparedMixedSerialExecution,
    workspace::_MixedSerialExecutionWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        workspace.due_paths[slot] || continue
        _materialize_mixed_serial_path_input!(
            prepared.paths[slot], workspace.paths[slot],
            workspace.batch)::Nothing
    end
    return nothing
end

function _submit_mixed_serial_due_paths!(
    prepared::_PreparedMixedSerialExecution,
    workspace::_MixedSerialExecutionWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        workspace.due_paths[slot] || continue
        _submit_mixed_serial_path_input!(
            prepared.paths[slot], workspace.paths[slot])::Nothing
    end
    return nothing
end

function _complete_mixed_serial_due_paths!(
    prepared::_PreparedMixedSerialExecution,
    workspace::_MixedSerialExecutionWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        workspace.due_paths[slot] || continue
        _complete_mixed_serial_path_input!(
            prepared.paths[slot], workspace.paths[slot])::Nothing
    end
    return nothing
end

function _apply_mixed_serial_due_path_inputs!(
    prepared::_PreparedMixedSerialExecution,
    workspace::_MixedSerialExecutionWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        workspace.due_paths[slot] || continue
        _apply_mixed_serial_path_input!(
            prepared.paths[slot], workspace.paths[slot])::Nothing
    end
    return nothing
end

function _execute_and_reclaim_mixed_serial_due_paths!(
    prepared::_PreparedMixedSerialExecution,
    workspace::_MixedSerialExecutionWorkspace,
)
    @inbounds for slot in eachindex(prepared.paths)
        workspace.due_paths[slot] || continue
        _execute_mixed_serial_path!(prepared.paths[slot])::Nothing
        _reclaim_mixed_serial_path_input!(
            prepared.paths[slot], workspace.paths[slot])::Nothing
    end
    return nothing
end

function _execute_mixed_serial_paths!(
    prepared::_PreparedMixedSerialExecution,
    state::_MixedSerialExecutionState,
    workspace::_MixedSerialExecutionWorkspace,
    due_paths,
    timestamp::PlantTimestamp,
)
    _require_mixed_serial_binding(prepared, state, workspace)
    _require_mixed_serial_timestamp(state, timestamp)
    try
        _select_mixed_serial_due_paths!(
            workspace.due_paths, prepared.path_ids, due_paths)
    catch
        fill!(workspace.due_paths, false)
        rethrow()
    end

    authority = prepared_atmosphere_authority(prepared.partitions)
    atmosphere = prepared_atmosphere(authority)
    target_time, next_epoch_sequence = try
        target, sequence =
            _preflight_mixed_serial_atmosphere_time(atmosphere, timestamp)
        workspace.batch.next_epoch_sequence = sequence
        workspace.batch.timestamp = timestamp
        _preflight_mixed_serial_due_paths(
            prepared,
            workspace.due_paths,
            workspace.batch,
        )
        (target, sequence)
    catch
        fill!(workspace.due_paths, false)
        rethrow()
    end

    state.phase = _MixedSerialExecutionMaterializing
    try
        atmosphere_rng = _prepared_atmosphere_rng(
            atmosphere, getfield(authority, :rngs))
        epoch = _with_completed_prepared_device_execution_context(
            authority.context) do
            advance_to!(atmosphere, target_time, atmosphere_rng)
        end
        epoch_sequence(epoch) == next_epoch_sequence ||
            _mixed_serial_execution_error(
                :atmosphere_epoch,
                "advanced atmosphere did not publish the preflighted epoch",
            )
        workspace.batch.epoch = epoch
        _materialize_mixed_serial_due_paths!(
            prepared, workspace)
        _submit_mixed_serial_due_paths!(prepared, workspace)
        _complete_mixed_serial_due_paths!(prepared, workspace)
        _apply_mixed_serial_due_path_inputs!(prepared, workspace)
        state.phase = _MixedSerialExecutionExecuting
        _execute_and_reclaim_mixed_serial_due_paths!(prepared, workspace)
    catch
        state.phase = _MixedSerialExecutionFailed
        rethrow()
    end

    fill!(workspace.due_paths, false)
    state.last_timestamp = timestamp
    state.has_timestamp = true
    state.phase = _MixedSerialExecutionIdle
    return nothing
end

function _execute_mixed_serial_paths!(
    ::_PreparedMixedSerialExecution,
    ::_MixedSerialExecutionState,
    ::_MixedSerialExecutionWorkspace,
    due_paths,
    timestamp,
)
    _mixed_serial_execution_error(
        :invalid_timestamp,
        "mixed serial optical execution requires a PlantTimestamp; got " *
        "$(typeof(timestamp))",
    )
end
