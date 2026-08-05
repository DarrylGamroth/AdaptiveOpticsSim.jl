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
mutable struct _PlantEventLoopStateBinding end
mutable struct _PlantEventLoopWorkspaceBinding end

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

"""
One contiguous run of path bindings that share a prepared geometric coupling.
The run preserves canonical binding order and never combines command state.
"""
struct _PreparedControllableOpticPathCouplingGroup
    first_binding::Int
    binding_count::Int
    representative_coupling_slot::UInt32
end

struct _PreparedPathExecutionGroupAcquisition
    slot::UInt32
    id::AcquisitionID
end

@enum _OpticalPathBatchPhase::UInt8 begin
    _OpticalPathBatchIdle = 0x00
    _OpticalPathBatchMaterializing = 0x01
    _OpticalPathBatchExecuting = 0x02
    _OpticalPathBatchFinalizing = 0x03
    _OpticalPathBatchAbandoned = 0x04
end

@enum _OpticalPathBatchGroupStatus::UInt8 begin
    _OpticalPathBatchGroupNotDue = 0x00
    _OpticalPathBatchGroupAwaitingMaterialization = 0x01
    _OpticalPathBatchGroupMaterializing = 0x02
    _OpticalPathBatchGroupReady = 0x03
    _OpticalPathBatchGroupExecuting = 0x04
    _OpticalPathBatchGroupComplete = 0x05
end

mutable struct _OpticalPathBatchWorkspace{E<:AtmosphereEpoch}
    binding::_PlantEventLoopWorkspaceBinding
    state_binding::_PlantEventLoopStateBinding
    generation::UInt64
    scheduler_revision::UInt64
    timestamp::PlantTimestamp
    epoch::E
    has_epoch::Bool
    phase::_OpticalPathBatchPhase
    due_group_slots::Memory{UInt32}
    due_group_count::Int
    group_status::Memory{_OpticalPathBatchGroupStatus}
end

struct _OpticalPathBatchClaimToken end
const _OPTICAL_PATH_BATCH_CLAIM_TOKEN = _OpticalPathBatchClaimToken()

"""
Opaque ownership claim for one same-timestamp optical-path batch.

The claim identifies a prepared event loop, exact state and workspace owners,
one scheduler revision, one atmosphere epoch, and one bounded due-group set.
It does not retain atmosphere layer storage.
"""
struct OpticalPathBatchClaim
    plant_binding_id::UInt64
    state_binding_id::UInt64
    workspace_binding_id::UInt64
    generation::UInt64
    scheduler_revision::UInt64
    timestamp::PlantTimestamp
    epoch_sequence::UInt64
    has_epoch::Bool

    function OpticalPathBatchClaim(
        plant_binding_id::UInt64,
        state_binding_id::UInt64,
        workspace_binding_id::UInt64,
        generation::UInt64,
        scheduler_revision::UInt64,
        timestamp::PlantTimestamp,
        epoch_sequence::UInt64,
        has_epoch::Bool,
        ::_OpticalPathBatchClaimToken,
    )
        return new(
            plant_binding_id,
            state_binding_id,
            workspace_binding_id,
            generation,
            scheduler_revision,
            timestamp,
            epoch_sequence,
            has_epoch,
        )
    end
end

"""
Execution-policy seam for one due optical-path batch.

Core supplies the deterministic serial implementation. HIL runtimes may
subtype this interface to coordinate externally owned execution workers while
using the same bounded batch lifecycle.
"""
abstract type AbstractOpticalPathBatchExecutor end

"""Canonical deterministic serial optical-path batch executor."""
struct SerialOpticalPathBatchExecutor <: AbstractOpticalPathBatchExecutor end

"""
Immutable compute-domain requirements for one prepared path execution group.

The semantic array `backend` and concrete `compute_device` are intentionally
separate. `requires_full_optical` distinguishes groups that must materialize
and execute their prepared optical path from reduced-order or replay groups.
This contract describes requirements only; it assigns no worker, stream,
placement, or scheduling policy.
"""
struct PathExecutionRequirements
    backend::AbstractArrayBackend
    compute_device::AbstractComputeDevice
    requires_full_optical::Bool

    function PathExecutionRequirements(
        backend::B,
        compute_device::D,
        requires_full_optical::Bool,
    ) where {B<:AbstractArrayBackend,D<:AbstractComputeDevice}
        typeof(compute_device_backend(compute_device)) === B || throw(
            PlantPreparationError(
                :path_execution_group,
                :compute_device,
                "path execution backend and compute-device family differ",
            ))
        return new(backend, compute_device, requires_full_optical)
    end
end

"""Cold compatibility result for one prepared group and exact target device."""
abstract type AbstractPathExecutionTargetSupport end

"""The prepared group already owns all state required on the exact target."""
struct SupportedPathExecutionTarget <: AbstractPathExecutionTargetSupport end

"""The prepared group cannot execute on the target for the structured `reason`."""
struct UnsupportedPathExecutionTarget <:
    AbstractPathExecutionTargetSupport
    reason::Symbol
end

function path_execution_target_supported(
    ::SupportedPathExecutionTarget,
)
    return true
end

function path_execution_target_supported(
    ::UnsupportedPathExecutionTarget,
)
    return false
end
@inline path_execution_target_rejection_reason(
    ::SupportedPathExecutionTarget,
) = nothing
@inline path_execution_target_rejection_reason(
    support::UnsupportedPathExecutionTarget,
) = support.reason

"""
Report whether an already prepared path group can execute on `target`.

Prepared groups are bound to one exact device. A same-family alternate device
therefore requires target-local preparation rather than an implicit runtime
move. Resource-local preparation extends this cold contract; warmed execution
never selects or migrates a target.
"""
@inline function path_execution_target_support(
    requirements::PathExecutionRequirements,
    target::AbstractComputeDevice,
)
    typeof(compute_device_backend(target)) === typeof(requirements.backend) ||
        return UnsupportedPathExecutionTarget(:backend_mismatch)
    target == requirements.compute_device ||
        return UnsupportedPathExecutionTarget(:requires_repreparation)
    return SupportedPathExecutionTarget()
end

"""
Run-immutable owner of one prepared direction-dependent optical path and every
compatible acquisition consumer scheduled on that path.

The group owns path-local mutable products, workspaces, and RNG streams through
`path`; repeated execution therefore has exactly one writer. `id` is the stable
declared path identity and does not depend on resource or owner order. `ordinal`
and acquisition slots are run-local bounded references. The group owns no task,
queue, CPU affinity, or transport policy.
"""
struct PreparedPathExecutionGroup
    ordinal::UInt32
    id::OpticalPathID
    requirements::PathExecutionRequirements
    path::PreparedPathExecutor
    # Retain the exact atmosphere as a prepared heterogeneous owner. This
    # avoids reboxing the immutable atmosphere when the path preflight crosses
    # its model-specific function barrier.
    atmosphere::AbstractTimedAtmosphere
    # Retain the heterogeneous input box created during preparation. The
    # event-loop barrier reuses it across all visible optics instead of
    # repeatedly boxing an immutable optical-product wrapper.
    input::Union{AbstractOpticalProduct,Tuple}
    rngs::PreparedOwnerRNGs
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    handle::EventGeneratorHandle
    sampled_aberration_binding_start::Int
    sampled_aberration_binding_stop::Int
    optic_binding_start::Int
    optic_binding_stop::Int
    optic_couplings::Memory{AbstractPupilSurfacePathCoupling}
    optic_coupling_groups::Memory{
        _PreparedControllableOpticPathCouplingGroup}
    acquisitions::Memory{_PreparedPathExecutionGroupAcquisition}
    autonomous_optic_slots::Memory{UInt32}
end

struct _PreparedPlantEventAcquisition
    id::AcquisitionID
    # The event loop deliberately retains a topology-bounded heterogeneous
    # registry. Store immutable lifecycle implementations behind their
    # abstract owner seam so they are boxed once during preparation instead of
    # being reboxed from an inline union on every event.
    lifecycle::AbstractPreparedAcquisitionLifecycle
    products::AcquisitionProducts
    product::Union{AbstractArray,WFSMeasurement}
    sample_provider::Union{Nothing,PreparedLinearReducedOrderEventProvider}
    rng::Xoshiro
    start::AbstractAcquisitionStartDefinition
    path_slot::UInt32
    start_handle::EventGeneratorHandle
    boundary_handle::EventGeneratorHandle
    band_open_handle::EventGeneratorHandle
    readout_handle::EventGeneratorHandle
    readiness_handle::EventGeneratorHandle
end

struct _PreparedPlantEventCommandEndpoint
    binding::_PreparedPlantCommandEndpoint
    endpoint::PreparedCommandEndpoint
    id::CommandEndpointID
    handle::EventGeneratorHandle
end

mutable struct _PreparedDevicePathBatchOwnerBinding end
struct _PreparedDevicePathBatchOwnerToken end
const _PREPARED_DEVICE_PATH_BATCH_OWNER_TOKEN =
    _PreparedDevicePathBatchOwnerToken()

"""
Run-immutable submission owner for one compatible, single-device set of
prepared path execution groups.

`implementation` is one small model-specific prepared pipeline behind a
topology-bounded concrete wrapper. It retains numerical plans, workspace,
backend execution context, and an explicit completion boundary, but no task,
queue, pacing, affinity, or transport policy.
"""
struct PreparedDevicePathBatchOwner
    binding::_PreparedDevicePathBatchOwnerBinding
    event_loop_binding::_PlantEventLoopBinding
    device::AbstractComputeDevice
    group_slots::Memory{UInt32}
    implementation::Any

    function PreparedDevicePathBatchOwner(
        event_loop_binding::_PlantEventLoopBinding,
        device::AbstractComputeDevice,
        group_slots::Memory{UInt32},
        implementation,
        ::_PreparedDevicePathBatchOwnerToken,
    )
        return new(
            _PreparedDevicePathBatchOwnerBinding(),
            event_loop_binding,
            device,
            group_slots,
            implementation,
        )
    end
end

"""Run-immutable, fixed-capacity deterministic plant-event composition."""
struct PreparedPlantEventLoop{
    C<:AbstractComputeDevice,
    X,
    A<:AbstractTimedAtmosphere,
    R,
    O<:_PreparedControllableOpticRegistry,
    T,
}
    binding::_PlantEventLoopBinding
    target::C
    context::X
    atmosphere::A
    atmosphere_rng::R
    scheduler::PreparedEventScheduler
    actions::Memory{_PlantEventAction}
    optics::O
    optic_path_bindings::PreparedControllableOpticPathBindings
    sampled_aberrations::Memory{PreparedSampledAberration}
    sampled_aberration_path_bindings::
        PreparedSampledAberrationPathBindings
    command_endpoints::Memory{_PreparedPlantEventCommandEndpoint}
    path_groups::Memory{PreparedPathExecutionGroup}
    device_path_batch_owners::Memory{PreparedDevicePathBatchOwner}
    path_device_batch_owner_slots::Memory{UInt32}
    acquisitions::Memory{_PreparedPlantEventAcquisition}
    autonomous_optics::Memory{_PreparedAutonomousPeriodicOptic}
    trigger_topology::T
end

@inline compute_device(prepared::PreparedPlantEventLoop) = prepared.target

@inline plant_event_path_count(prepared::PreparedPlantEventLoop) =
    length(prepared.path_groups)
@inline path_execution_group_count(prepared::PreparedPlantEventLoop) =
    length(prepared.path_groups)

function path_execution_group(prepared::PreparedPlantEventLoop,
    ordinal::Integer)
    checkbounds(prepared.path_groups, ordinal)
    return @inbounds prepared.path_groups[ordinal]
end

@inline path_execution_group_ordinal(group::PreparedPathExecutionGroup) =
    Int(group.ordinal)
@inline path_execution_group_path_id(group::PreparedPathExecutionGroup) =
    group.id
@inline path_execution_group_requirements(
    group::PreparedPathExecutionGroup,
) = group.requirements
@inline path_execution_target_support(
    group::PreparedPathExecutionGroup,
    target::AbstractComputeDevice,
) = path_execution_target_support(group.requirements, target)
@inline path_execution_backend(requirements::PathExecutionRequirements) =
    requirements.backend
@inline path_execution_compute_device(
    requirements::PathExecutionRequirements,
) = requirements.compute_device
@inline path_execution_requires_full_optical(
    requirements::PathExecutionRequirements,
) = requirements.requires_full_optical
@inline path_execution_group_acquisition_count(
    group::PreparedPathExecutionGroup) = length(group.acquisitions)

function path_execution_group_acquisition_id(
    group::PreparedPathExecutionGroup, index::Integer)
    checkbounds(group.acquisitions, index)
    return @inbounds group.acquisitions[index].id
end

"""Return the exact canonical plant timestamp owned by a batch claim."""
@inline optical_path_batch_timestamp(claim::OpticalPathBatchClaim) =
    claim.timestamp

@inline plant_event_acquisition_count(prepared::PreparedPlantEventLoop) =
    length(prepared.acquisitions)
@inline plant_event_generator_count(prepared::PreparedPlantEventLoop) =
    event_generator_count(prepared.scheduler)
@inline plant_event_controllable_optic_count(
    prepared::PreparedPlantEventLoop) = length(prepared.optics)
@inline plant_event_sampled_aberration_count(
    prepared::PreparedPlantEventLoop) =
    length(prepared.sampled_aberrations)
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
    _require_prepared_acquisition(execution.acquisition, result)
    return _prepare_global_shutter_acquisition(execution.acquisition,
        definition)
end

@inline function _prepare_detector_event_lifecycle(
    execution::FrameAcquisitionExecution, result::IntensityMap,
    definition::RollingShutterAcquisitionDefinition)
    _require_prepared_acquisition(execution.acquisition, result)
    return _prepare_rolling_shutter_acquisition(execution.acquisition,
        definition)
end

@inline function _prepare_detector_event_lifecycle(
    execution::FrameAcquisitionExecution, result::IntensityMap,
    definition::FrameTransferAcquisitionDefinition)
    _require_prepared_acquisition(execution.acquisition, result)
    return _prepare_frame_transfer_acquisition(execution.acquisition,
        definition)
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
        path_id(path) == id && return path
    end
    _plant_event_loop_error(:unknown_path,
        "prepared plant has no optical path $id")
end

function _event_prepared_acquisition(plant::PreparedPlant,
    id::AcquisitionID)
    Base.@nospecialize plant
    for owner in getfield(plant, :acquisitions)
        acquisition_id(owner) == id && return owner
    end
    _plant_event_loop_error(:unknown_acquisition,
        "prepared plant has no acquisition $id")
end

function _event_path_group_slot(groups, id::OpticalPathID)
    @inbounds for index in eachindex(groups)
        groups[index].id == id && return index
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

function _phase_reference_trigger_consumer(
    reference::AbstractWaveformPhaseReference)
    _plant_event_loop_error(
        :unsupported_phase_reference,
        "plant event composition does not support phase reference $(typeof(reference)); use a declared free-running, trigger-source, or trigger-reset relationship",
    )
end

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
    return getfield(plant, :controllable_optics)
end

function _prepare_event_sampled_aberrations(plant::PreparedPlant)
    source = getfield(plant, :sampled_aberrations)
    aberrations = Memory{PreparedSampledAberration}(
        undef, length(source))
    @inbounds for index in eachindex(source)
        aberrations[index] = source[index]
    end
    return aberrations
end

function _prepare_event_command_endpoints(plant::PreparedPlant,
    scheduler::PreparedEventScheduler)
    source = getfield(plant, :command_endpoints)
    endpoints = Memory{_PreparedPlantEventCommandEndpoint}(
        undef, length(source))
    @inbounds for index in eachindex(source)
        endpoints[index] = _PreparedPlantEventCommandEndpoint(
            source[index], source[index].endpoint,
            command_endpoint_id(source[index]),
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
    has_acquisition = false
    @inbounds for owner in owners
        acquisition_path_id(owner) == id || continue
        has_acquisition = true
        _provider_requires_full_optical(acquisition_provider_style(owner)) &&
            return true
    end
    # An explicitly scheduled path with no acquisition consumer is itself the
    # demand for a device-ready optical product. Reduced-order and replay-only
    # consumers remain the only reason to bypass that path's optical work.
    return !has_acquisition
end

function _require_prepared_event_path_coupling(
    coupling::AbstractPupilSurfacePathCoupling,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    Base.ismutabletype(typeof(coupling)) && _plant_event_loop_error(
        :mutable_optic_path_coupling,
        "controllable optic $optic returned mutable path coupling " *
        "$(typeof(coupling)) for $path")
    return coupling
end

function _require_prepared_event_path_coupling(
    coupling,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    _plant_event_loop_error(:invalid_optic_path_coupling,
        "controllable optic $optic must return an " *
        "AbstractPupilSurfacePathCoupling for $path; got " *
        "$(typeof(coupling))")
end

function _prepare_event_path_optic_coupling(
    ::PupilSurfaceExecutionRole,
    optic::PreparedControllableOptic,
    definition::ControllableOpticDefinition,
    path::PreparedPathExecutor,
)
    coupling = prepare_controllable_optic_path_coupling(
        optic.implementation, definition, path)
    return _require_prepared_event_path_coupling(
        coupling,
        controllable_optic_id(definition),
        path_id(path),
    )
end

@inline function _prepare_event_path_optic_coupling(
    ::AutonomousPathExecutionRole,
    optic::PreparedControllableOptic,
    definition::ControllableOpticDefinition,
    path::PreparedPathExecutor,
)
    return _prepare_event_autonomous_path_coupling(
        controllable_optic_placement(definition),
        controllable_optic_id(definition),
        path_id(path),
    )
end

@inline function _prepare_event_autonomous_path_coupling(
    ::FocalPlanePlacement,
    ::ControllableOpticID,
    ::OpticalPathID,
)
    return _NoPupilSurfacePathCoupling()
end

function _prepare_event_autonomous_path_coupling(
    placement::AbstractOpticalPlacement,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    _plant_event_loop_error(:invalid_optic_placement,
        "path-local autonomous controllable optic $optic requires " *
        "FocalPlanePlacement on path $path; got $(typeof(placement))")
end

function _prepare_event_path_optic_coupling(
    role::AbstractControllableOpticExecutionRole,
    optic::PreparedControllableOptic,
    definition::ControllableOpticDefinition,
    path::PreparedPathExecutor,
)
    _plant_event_loop_error(:unsupported_optic_execution_role,
        "controllable optic $(controllable_optic_id(definition)) on " *
        "path $(path_id(path)) declares unsupported execution " *
        "role $(typeof(role))")
end

@inline _coupling_group_member(
    ::PupilSurfaceExecutionRole) = true
@inline _coupling_group_member(
    ::AutonomousPathExecutionRole) = false
@inline _coupling_group_member(
    ::AbstractControllableOpticExecutionRole) = false

@inline function _append_event_path_coupling_group!(
    groups::Vector{_PreparedControllableOpticPathCouplingGroup},
    first_binding::Int,
    binding_count::Int,
    representative_coupling_slot::Int,
)
    iszero(binding_count) && return nothing
    representative_coupling_slot <= typemax(UInt32) ||
        _plant_event_loop_error(:capacity,
            "path-local controllable-optic coupling count exceeds UInt32 " *
            "capacity")
    push!(groups, _PreparedControllableOpticPathCouplingGroup(
        first_binding,
        binding_count,
        UInt32(representative_coupling_slot),
    ))
    return nothing
end

function _prepare_event_path_optic_couplings(
    plant::PreparedPlant,
    path::PreparedPathExecutor,
    bindings::PreparedControllableOpticPathBindings,
    binding_range::UnitRange{Int},
    requires_full_optical::Bool,
)
    coupling_count = requires_full_optical ? length(binding_range) : 0
    couplings = Memory{AbstractPupilSurfacePathCoupling}(
        undef, coupling_count)
    groups = _PreparedControllableOpticPathCouplingGroup[]
    requires_full_optical || return couplings,
        Memory{_PreparedControllableOpticPathCouplingGroup}(undef, 0)

    optics = getfield(plant, :controllable_optics)
    group_first_binding = 0
    group_binding_count = 0
    group_representative_slot = 0
    previous_coupling = nothing
    @inbounds for (coupling_slot, binding) in
        enumerate(binding_range)
        optic_slot = prepared_controllable_optic_slot(bindings, binding)
        optic = optics[optic_slot]
        definition =
            _prepared_controllable_optic_definition(optics, optic)
        role = controllable_optic_execution_role(optic.implementation)
        coupling = _prepare_event_path_optic_coupling(
            role, optic, definition, path)
        couplings[coupling_slot] = coupling
        if !_coupling_group_member(role)
            _append_event_path_coupling_group!(
                groups,
                group_first_binding,
                group_binding_count,
                group_representative_slot,
            )
            group_binding_count = 0
            previous_coupling = nothing
            continue
        end
        if iszero(group_binding_count)
            group_first_binding = binding
            group_binding_count = 1
            group_representative_slot = coupling_slot
        elseif _same_pupil_footprint_coupling(
            previous_coupling, coupling)
            group_binding_count += 1
        else
            _append_event_path_coupling_group!(
                groups,
                group_first_binding,
                group_binding_count,
                group_representative_slot,
            )
            group_first_binding = binding
            group_binding_count = 1
            group_representative_slot = coupling_slot
        end
        previous_coupling = coupling
    end
    _append_event_path_coupling_group!(
        groups,
        group_first_binding,
        group_binding_count,
        group_representative_slot,
    )
    group_memory =
        Memory{_PreparedControllableOpticPathCouplingGroup}(
            undef, length(groups))
    copyto!(group_memory, groups)
    return couplings, group_memory
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

function _prepare_path_execution_group_acquisitions(
    acquisition_definitions,
    acquisition_path_slots,
    path_slot::Int,
)
    count = 0
    @inbounds for acquisition_slot in eachindex(acquisition_path_slots)
        acquisition_path_slots[acquisition_slot] == path_slot &&
            (count += 1)
    end
    acquisitions = Memory{_PreparedPathExecutionGroupAcquisition}(
        undef, count)
    member = 1
    @inbounds for acquisition_slot in eachindex(acquisition_path_slots)
        acquisition_path_slots[acquisition_slot] == path_slot || continue
        acquisitions[member] = _PreparedPathExecutionGroupAcquisition(
            UInt32(acquisition_slot),
            acquisition_definitions[acquisition_slot].acquisition,
        )
        member += 1
    end
    return acquisitions
end

function _prepare_path_execution_group_autonomous_optic_slots(
    definitions,
    path::OpticalPathID,
)
    count = 0
    @inbounds for definition in definitions
        definition.path == path && (count += 1)
    end
    slots = Memory{UInt32}(undef, count)
    member = 1
    @inbounds for index in eachindex(definitions)
        definitions[index].path == path || continue
        index <= typemax(UInt32) || _plant_event_loop_error(
            :capacity,
            "autonomous-optic registry exceeds UInt32 path-group capacity",
        )
        slots[member] = UInt32(index)
        member += 1
    end
    return slots
end

function _prepare_path_execution_groups(
    plant::PreparedPlant, definitions, owners,
    acquisition_definitions, acquisition_path_slots,
    autonomous_definitions,
    scheduler::PreparedEventScheduler)
    Base.@nospecialize plant
    groups = Memory{PreparedPathExecutionGroup}(undef, length(definitions))
    bindings = getfield(plant, :controllable_optic_path_bindings)
    sampled_bindings =
        getfield(plant, :sampled_aberration_path_bindings)
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        path = _event_prepared_path(plant, definition.path)
        rngs = _prepared_event_path_rngs(plant, path)
        _require_rng_owner_binding(rngs, path)
        handle = event_generator_handle(scheduler, OpticalSamplePhase, index)
        binding_range = prepared_controllable_optic_binding_range(
            bindings, definition.path)
        sampled_binding_range =
            prepared_sampled_aberration_binding_range(
                sampled_bindings, definition.path)
        requires_full_optical =
            _event_path_requires_full_optical(definition.path, owners)
        optic_couplings, optic_coupling_groups =
            _prepare_event_path_optic_couplings(
                plant,
                path,
                bindings,
                binding_range,
                requires_full_optical,
            )
        acquisitions =
            _prepare_path_execution_group_acquisitions(
                acquisition_definitions,
                acquisition_path_slots,
                index,
            )
        autonomous_optic_slots =
            _prepare_path_execution_group_autonomous_optic_slots(
                autonomous_definitions,
                definition.path,
            )
        groups[index] = PreparedPathExecutionGroup(
            UInt32(index),
            definition.path,
            PathExecutionRequirements(
                getfield(path.key, :backend),
                getfield(path.key, :device),
                requires_full_optical,
            ),
            path,
            path.atmosphere,
            path.input,
            rngs,
            definition.schedule,
            definition.origin,
            handle,
            first(sampled_binding_range), last(sampled_binding_range),
            first(binding_range), last(binding_range),
            optic_couplings, optic_coupling_groups,
            acquisitions, autonomous_optic_slots)
    end
    return groups
end

function _event_controllable_optic_slot(optics,
    id::ControllableOpticID)
    @inbounds for index in eachindex(optics)
        definition = _prepared_controllable_optic_definition(optics, index)
        controllable_optic_id(definition) == id &&
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
    @inbounds for index in eachindex(optics)
        optic = optics[index]
        definition = _prepared_controllable_optic_definition(optics, optic)
        id = controllable_optic_id(definition)
        role = controllable_optic_execution_role(optic.implementation)
        _require_autonomous_definition_for_role(role, id, definitions)
    end
    return nothing
end

function _require_autonomous_path_visibility(
    ::AllPathVisibility,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    _plant_event_loop_error(:invalid_autonomous_visibility,
        "path-local autonomous controllable optic $optic must select only " *
        "its coupled path $path")
end

function _require_autonomous_path_visibility(
    visibility::SelectedPathVisibility,
    optic::ControllableOpticID,
    path::OpticalPathID,
)
    paths = selected_path_ids(visibility)
    length(paths) == 1 && only(paths) == path && return nothing
    _plant_event_loop_error(:invalid_autonomous_visibility,
        "path-local autonomous controllable optic $optic must select only " *
        "its coupled path $path; got $(paths)")
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

function _prepare_event_autonomous_optics(definitions, optics, path_groups,
    topology)
    bindings = Memory{_PreparedAutonomousPeriodicOptic}(undef,
        length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        _require_autonomous_trigger_topology(definition, topology)
        optic_slot = _event_controllable_optic_slot(optics,
            definition.optic)
        optic = optics[optic_slot]
        optic_definition =
            _prepared_controllable_optic_definition(optics, optic)
        _require_autonomous_execution_role(
            controllable_optic_execution_role(optic.implementation),
            definition.optic)
        _require_autonomous_path_visibility(
            controllable_optic_visibility(optic_definition), definition.optic,
            definition.path)
        path_slot = _event_path_group_slot(path_groups, definition.path)
        group = path_groups[path_slot]
        group.requirements.requires_full_optical || _plant_event_loop_error(
            :autonomous_path_without_full_optics,
            "autonomous optic $(definition.optic) targets path " *
            "$(definition.path) without full-optical execution demand")
        coupling = prepare_autonomous_periodic_optic(
            optic.implementation, group.path, definition.fidelity)
        bindings[index] = _PreparedAutonomousPeriodicOptic(
            definition.optic, definition.path, UInt32(optic_slot),
            UInt32(path_slot), optic.implementation, coupling,
            definition.phase_reference, definition.fidelity)
    end
    _require_all_autonomous_optics_bound(optics, definitions)
    _require_unique_autonomous_optic_couplings(bindings)
    return bindings
end

function _event_path_visible_command_endpoint_slots(
    plant::PreparedPlant,
    path::OpticalPathID,
)
    bindings = getfield(plant, :controllable_optic_path_bindings)
    optics = getfield(plant, :controllable_optics)
    binding_range =
        prepared_controllable_optic_binding_range(bindings, path)
    endpoint_slots = UInt32[]
    @inbounds for binding in binding_range
        optic_slot =
            prepared_controllable_optic_slot(bindings, binding)
        append!(endpoint_slots, optics[optic_slot].endpoint_slots)
    end
    sort!(endpoint_slots)
    @inbounds for index in 2:length(endpoint_slots)
        endpoint_slots[index - 1] == endpoint_slots[index] &&
            _plant_event_loop_error(
                :duplicate_visible_command_endpoint,
                "optical path $path resolves command endpoint slot " *
                "$(endpoint_slots[index]) more than once",
            )
    end
    return Tuple(endpoint_slots)
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
    visible_endpoint_slots =
        _event_path_visible_command_endpoint_slots(
            plant,
            acquisition_path_id(owner),
        )
    sample_provider = prepare_linear_reduced_order_event_provider(
        implementation,
        getfield(plant, :command_endpoints),
        visible_endpoint_slots,
    )
    return lifecycle, measurement, sample_provider
end

function _prepare_event_acquisition_lifecycle(
    ::PreparedPlant, owner::PreparedAcquisitionOwner,
    definition::AbstractAcquisitionLifecycleDefinition,
    ::CommandResponsiveReducedOrderProviderStyle)
    _plant_event_loop_error(:unsupported_acquisition,
        "command-responsive reduced-order acquisition $(acquisition_id(owner)) requires DirectMeasurementAcquisitionDefinition; got $(typeof(definition))")
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
        "full-optical acquisition $(acquisition_id(owner)) requires a detector lifecycle, not DirectMeasurementAcquisitionDefinition")
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
            acquisition_path_id(owner)))
    end
    return owners, lifecycles, products, sample_providers, rngs, path_slots
end

function _prepare_event_acquisitions(definitions, owners, lifecycles, products,
    sample_providers, rngs, path_slots,
    scheduler::PreparedEventScheduler)
    acquisitions = Memory{_PreparedPlantEventAcquisition}(undef,
        length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        acquisitions[index] = _PreparedPlantEventAcquisition(
            definition.acquisition, lifecycles[index],
            acquisition_products(owners[index]),
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

function _prepare_plant_event_loop_in_context(
    plant::PreparedPlant,
    definition::PlantEventLoopDefinition,
    device_batch_selection::Val,
)
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
    optic_path_bindings =
        getfield(plant, :controllable_optic_path_bindings)
    sampled_aberrations = _prepare_event_sampled_aberrations(plant)
    sampled_aberration_path_bindings =
        getfield(plant, :sampled_aberration_path_bindings)
    command_endpoints = _prepare_event_command_endpoints(plant, scheduler)
    path_groups = _prepare_path_execution_groups(
        plant, sample_definitions, owners, acquisition_definitions,
        path_slots, autonomous_definitions, scheduler)
    acquisitions = _prepare_event_acquisitions(
        acquisition_definitions, owners, lifecycles, products,
        sample_providers, rngs, path_slots, scheduler)
    autonomous_optics = _prepare_event_autonomous_optics(
        autonomous_definitions, optics, path_groups,
        definition.trigger_topology)
    atmosphere = _require_selection_atmosphere(prepared_atmosphere(plant))
    atmosphere_rng = _prepared_atmosphere_rng(atmosphere,
        getfield(getfield(plant, :rngs), :atmosphere))
    binding = _PlantEventLoopBinding()
    device_batch_owners, path_device_batch_owner_slots =
        _prepare_device_path_batch_owners(
            device_batch_selection,
            binding,
            path_groups,
            atmosphere,
        )
    prepared = PreparedPlantEventLoop(binding, compute_device(plant),
        getfield(plant, :context), atmosphere,
        atmosphere_rng, scheduler, actions, optics, optic_path_bindings,
        sampled_aberrations, sampled_aberration_path_bindings,
        command_endpoints, path_groups, device_batch_owners,
        path_device_batch_owner_slots, acquisitions, autonomous_optics,
        _prepared_trigger_topology(definition.trigger_topology))
    return _require_exact_prepared_event_loop_target(
        prepared, compute_device(plant))
end

function _prepare_plant_event_loop(
    plant::PreparedPlant,
    definition::PlantEventLoopDefinition,
    device_batch_selection::Val,
)
    return _with_completed_prepared_device_execution_context(
        getfield(plant, :context)) do
        _prepare_plant_event_loop_in_context(
            plant, definition, device_batch_selection)
    end
end

"""
    prepare_plant_event_loop(plant, definition)

Bind a prepared plant to a flat deterministic scheduler, exact owner-derived
RNG streams, bounded detector lifecycle state, explicitly compatible
single-accelerator path batches, and an optional prepared trigger topology.
Prepared controllable optics and independently timed command endpoints join
the reserved command phase. Preparation allocates; repeated stepping does not
grow the registry or materialize future events. Core creates no execution
task, queue, pacing, affinity, or transport policy.
"""
function prepare_plant_event_loop(
    plant::PreparedPlant,
    definition::PlantEventLoopDefinition,
)
    return _prepare_plant_event_loop(
        plant,
        definition,
        Val(:accelerator),
    )
end

mutable struct PlantEventLoopState{T,O<:_ControllableOpticStateRegistry}
    binding::_PlantEventLoopBinding
    state_binding::_PlantEventLoopStateBinding
    scheduler::EventSchedulerState
    command_endpoints::Memory{CommandEndpointState}
    command_applications::Memory{CommandApplicationState}
    command_shadow_transactions::Memory{UInt64}
    controllable_optics::O
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
    ids = Tuple(map(optic.endpoint_slots) do slot
        command_endpoint_id(
            prepared.command_endpoints[Int(slot)].binding)
    end)
    commands = Tuple(map(optic.endpoint_slots) do slot
        binding = prepared.command_endpoints[Int(slot)].binding
        _copy_prepared_effective_command(binding.endpoint,
            binding.initial_command, "initial physical command",
            prepared.target)
    end)
    return ids, commands
end

function _event_controllable_optic_state(
    prepared::PreparedPlantEventLoop,
    optic::PreparedControllableOptic,
)
    ids, commands = _event_optic_endpoint_initials(prepared, optic)
    definition =
        _prepared_controllable_optic_definition(prepared.optics, optic)
    return prepare_controllable_optic_state(
        optic.implementation, definition, ids, commands)
end

function _event_controllable_optic_state_family(
    prepared::PreparedPlantEventLoop,
    family::_PreparedControllableOpticFamily,
)
    first_optic = first(family.values)
    first_state = _event_controllable_optic_state(prepared, first_optic)
    S = typeof(first_state)
    states = Vector{S}(undef, length(family.values))
    states[1] = first_state
    @inbounds for index in 2:length(states)
        state = _event_controllable_optic_state(
            prepared, family.values[index])
        typeof(state) === S || throw(PlantPreparationError(
            :controllable_optic,
            :state_family_type,
            "one exact prepared controllable-optic family returned " *
            "different physical-state types",
        ))
        states[index] = state
    end
    fixed = FixedSizeVectorDefault{S}(states)
    return _ControllableOpticStateFamily{S,typeof(fixed)}(fixed)
end

@inline _event_controllable_optic_state_families(
    prepared::PreparedPlantEventLoop,
    ::Tuple{},
) = ()

function _event_controllable_optic_state_families(
    prepared::PreparedPlantEventLoop,
    families::Tuple,
)
    first = _event_controllable_optic_state_family(prepared, families[1])
    rest = _event_controllable_optic_state_families(
        prepared, Base.tail(families))
    return (first, rest...)
end

function _event_controllable_optic_states(
    prepared::PreparedPlantEventLoop,
)
    optics = prepared.optics
    groups = _event_controllable_optic_state_families(
        prepared, optics.groups)
    return _ControllableOpticStateRegistry(groups, optics.slots)
end

function _initialize_event_autonomous_optics!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState)
    @inbounds for binding in prepared.autonomous_optics
        _initialize_event_autonomous_optic!(
            prepared,
            state,
            binding.optic_slot,
            binding.implementation,
            binding.phase_reference,
        )
    end
    return nothing
end

@inline function _initialize_event_autonomous_optic_family!(
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    implementation,
    phase_reference,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _initialize_event_autonomous_optic_family!(
    state_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    implementation,
    phase_reference,
)
    if family_slot == 1
        optic_state = @inbounds state_families[1].values[member_slot]
        initialize_autonomous_periodic_optic!(
            implementation, optic_state, phase_reference)
        return nothing
    end
    return _initialize_event_autonomous_optic_family!(
        Base.tail(state_families),
        family_slot - 1,
        member_slot,
        implementation,
        phase_reference,
    )
end

@inline function _initialize_event_autonomous_optic!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    optic_slot::UInt32,
    implementation,
    phase_reference,
)
    slot = @inbounds prepared.optics.slots[Int(optic_slot)]
    return _initialize_event_autonomous_optic_family!(
        state.controllable_optics.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        implementation,
        phase_reference,
    )
end

function _plant_event_loop_state(prepared::PreparedPlantEventLoop)
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
    path_sampled = Memory{Bool}(undef, length(prepared.path_groups))
    fill!(path_sampled, false)
    product_sequences = Memory{UInt64}(undef,
        length(prepared.acquisitions))
    fill!(product_sequences, UInt64(0))
    product_ready_timestamps = Memory{PlantTimestamp}(undef,
        length(prepared.acquisitions))
    fill!(product_ready_timestamps, zero(PlantTimestamp))
    trigger = _event_trigger_state(prepared.trigger_topology)
    state = PlantEventLoopState(prepared.binding,
        _PlantEventLoopStateBinding(),
        scheduler, command_endpoints, command_applications,
        command_shadow_transactions, controllable_optics,
        acquisition_states, path_sampled,
        product_sequences, product_ready_timestamps, UInt64(0), trigger)
    _initialize_event_autonomous_optics!(prepared, state)
    @inbounds for index in eachindex(prepared.command_endpoints)
        _schedule_event_command_endpoint!(prepared, state, index)
    end
    return _require_exact_plant_event_loop_state_target(
        prepared, state, prepared.target)
end

function PlantEventLoopState(prepared::PreparedPlantEventLoop)
    return _with_completed_prepared_device_execution_context(
        prepared.context) do
        _plant_event_loop_state(prepared)
    end
end

mutable struct PlantEventLoopWorkspace{
    T,
    B,
    O<:_ControllableOpticWorkspaceRegistry,
}
    binding::_PlantEventLoopBinding
    scheduler::EventSchedulerWorkspace
    command_endpoints::Memory{CommandDispositionWorkspace}
    controllable_optics::O
    command_dispositions::Memory{PlantCommandDisposition}
    command_disposition_count::Int
    transaction_endpoint_slots::Memory{UInt32}
    transaction_admissions::Memory{_CommandTransactionAdmissionPlan}
    transaction_claims::Memory{PlantCommandApplicationClaim}
    transaction_staged::Memory{_StagedCommandApplication}
    transaction_count::Int
    due_paths::Memory{Bool}
    optical_path_batch::B
    trigger::T
    delivery::Base.RefValue{TriggerDelivery}
end

@inline _event_trigger_workspace(::_NoPreparedTriggerTopology) =
    _NoTriggerTopologyWorkspace()
@inline _event_trigger_workspace(topology::PreparedTriggerTopology) =
    TriggerTopologyWorkspace(topology)

function _optical_path_batch_workspace(
    prepared::PreparedPlantEventLoop,
)
    group_count = length(prepared.path_groups)
    status = Memory{_OpticalPathBatchGroupStatus}(undef, group_count)
    fill!(status, _OpticalPathBatchGroupNotDue)
    atmosphere = prepared.atmosphere
    timeline = atmosphere_timeline(atmosphere)
    epoch = AtmosphereEpoch(
        atmosphere_identity(atmosphere),
        timeline.model_time,
        timeline.sequence,
    )
    return _OpticalPathBatchWorkspace(
        _PlantEventLoopWorkspaceBinding(),
        _PlantEventLoopStateBinding(),
        UInt64(0),
        UInt64(0),
        zero(PlantTimestamp),
        epoch,
        false,
        _OpticalPathBatchIdle,
        Memory{UInt32}(undef, group_count),
        0,
        status,
    )
end

function _event_controllable_optic_workspace_family(
    family::_PreparedControllableOpticFamily,
)
    first_workspace = prepare_controllable_optic_workspace(
        first(family.values).implementation)
    W = typeof(first_workspace)
    workspaces = Vector{W}(undef, length(family.values))
    workspaces[1] = first_workspace
    @inbounds for index in 2:length(workspaces)
        workspace = prepare_controllable_optic_workspace(
            family.values[index].implementation)
        typeof(workspace) === W || throw(PlantPreparationError(
            :controllable_optic,
            :workspace_family_type,
            "one exact prepared controllable-optic family returned " *
            "different physical-workspace types",
        ))
        workspaces[index] = workspace
    end
    fixed = FixedSizeVectorDefault{W}(workspaces)
    return _ControllableOpticWorkspaceFamily{W,typeof(fixed)}(fixed)
end

@inline _event_controllable_optic_workspace_families(::Tuple{}) = ()

function _event_controllable_optic_workspace_families(families::Tuple)
    first = _event_controllable_optic_workspace_family(families[1])
    rest = _event_controllable_optic_workspace_families(Base.tail(families))
    return (first, rest...)
end

function _event_controllable_optic_workspaces(
    prepared::PreparedPlantEventLoop,
)
    optics = prepared.optics
    groups = _event_controllable_optic_workspace_families(optics.groups)
    return _ControllableOpticWorkspaceRegistry(groups, optics.slots)
end

function _plant_event_loop_workspace(prepared::PreparedPlantEventLoop)
    command_endpoints = Memory{CommandDispositionWorkspace}(undef,
        length(prepared.command_endpoints))
    disposition_capacity = 0
    @inbounds for index in eachindex(prepared.command_endpoints)
        endpoint = prepared.command_endpoints[index].binding.endpoint
        command_endpoints[index] = CommandDispositionWorkspace(endpoint)
        disposition_capacity += command_endpoint_capacity(endpoint)
    end
    controllable_optics = _event_controllable_optic_workspaces(prepared)
    endpoint_count = length(prepared.command_endpoints)
    due_paths = Memory{Bool}(undef, length(prepared.path_groups))
    fill!(due_paths, false)
    optical_path_batch = _optical_path_batch_workspace(prepared)
    workspace = PlantEventLoopWorkspace(prepared.binding,
        EventSchedulerWorkspace(prepared.scheduler), command_endpoints,
        controllable_optics,
        Memory{PlantCommandDisposition}(undef, disposition_capacity), 0,
        Memory{UInt32}(undef, endpoint_count),
        Memory{_CommandTransactionAdmissionPlan}(undef, endpoint_count),
        Memory{PlantCommandApplicationClaim}(undef, endpoint_count),
        Memory{_StagedCommandApplication}(undef, endpoint_count), 0,
        due_paths, optical_path_batch,
        _event_trigger_workspace(prepared.trigger_topology),
        Ref{TriggerDelivery}())
    return _require_exact_plant_event_loop_workspace_target(
        prepared, workspace, prepared.target)
end

function PlantEventLoopWorkspace(prepared::PreparedPlantEventLoop)
    return _with_completed_prepared_device_execution_context(
        prepared.context) do
        _plant_event_loop_workspace(prepared)
    end
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
        state.controllable_optics.slots === prepared.optics.slots &&
        length(state.acquisitions) == length(prepared.acquisitions) &&
        length(state.path_sampled) == length(prepared.path_groups) &&
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
        workspace.controllable_optics.slots === prepared.optics.slots &&
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
    length(workspace.due_paths) == length(prepared.path_groups) ||
        _plant_event_loop_error(:prepared_binding,
            "plant event-loop workspace capacity changed after preparation")
    length(prepared.path_device_batch_owner_slots) ==
        length(prepared.path_groups) ||
        _plant_event_loop_error(:prepared_binding,
            "prepared device path-batch membership changed after preparation")
    batch = workspace.optical_path_batch
    length(batch.due_group_slots) == length(prepared.path_groups) &&
        length(batch.group_status) == length(prepared.path_groups) ||
        _plant_event_loop_error(:prepared_binding,
            "optical-path batch workspace capacity changed after preparation")
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

"""
Return the exact run-bound products for one prepared event-loop acquisition.

The returned wrapper and its storage are the same objects used by event-loop
execution; this cold binding seam does not copy or transfer ownership.
"""
function acquisition_products(
    prepared::PreparedPlantEventLoop,
    id)
    slot = _event_acquisition_slot(
        prepared, _as_acquisition_id(id))
    return @inbounds prepared.acquisitions[slot].products
end

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
        prepared.command_endpoints[index].id == id &&
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
    endpoint = _event_command_endpoint(prepared, slot).endpoint
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
    _require_idle_optical_path_batch(workspace)
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
            event_endpoint.endpoint,
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

"""
    fail_pending_plant_commands!(
        prepared, state, workspace, endpoint_id; reason=:endpoint_failure)

Boundedly terminate every unclaimed pending command for one event-loop-owned
endpoint at its current canonical coordinate. An ingress admission may have
advanced the endpoint beyond the scheduler's last processed plant event, so
the failure coordinate is the later of those two coordinates. Terminal
dispositions are copied into the event loop's bounded ledger and the
endpoint's command generator is rescheduled without changing its effective
command or physical optic state.
"""
function fail_pending_plant_commands!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_id::CommandEndpointID;
    reason=CommandDispositionReason(:endpoint_failure))
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    _require_command_failure_batch_state(workspace)
    _require_empty_event_command_dispositions(workspace)
    slot = _event_command_endpoint_slot(prepared, endpoint_id)
    event_endpoint = _event_command_endpoint(prepared, slot)
    endpoint_state = _event_command_endpoint_state(state, slot)
    endpoint_workspace = _event_command_workspace(workspace, slot)
    timestamp = max(
        scheduler_timestamp(state.scheduler),
        command_endpoint_timestamp(endpoint_state))
    count = fail_pending_plant_commands!(
        endpoint_workspace,
        event_endpoint.binding.endpoint,
        endpoint_state,
        timestamp;
        reason)
    _append_event_command_dispositions!(workspace, endpoint_workspace)
    _schedule_event_command_endpoint!(prepared, state, slot)
    return count
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
            if prior_binding.optic_slot == binding.optic_slot
                definition = _prepared_controllable_optic_definition(
                    prepared.optics, Int(binding.optic_slot))
                _command_admission_error(:transaction,
                    :duplicate_physical_optic,
                    "atomic multi-optic transaction contains more than one " *
                    "endpoint owned by controllable optic " *
                    "$(controllable_optic_id(definition))")
            end
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
    _require_idle_optical_path_batch(workspace)
    _require_empty_event_command_dispositions(workspace)
    _require_routed_command_admission_timestamp(prepared, state, workspace,
        timestamp)
    count = _prepare_transaction_member_slots!(prepared, workspace,
        transaction)
    transaction_id = try
        _next_command_transaction_sequence(
            state.command_transaction_sequence)
    catch
        workspace.transaction_count = 0
        rethrow()
    end
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

@inline function _stage_event_controllable_optic_family_command!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    application_state::CommandApplicationState,
    timestamp::PlantTimestamp,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _stage_event_controllable_optic_family_command!(
    prepared_families::Tuple,
    state_families::Tuple,
    workspace_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    application_state::CommandApplicationState,
    timestamp::PlantTimestamp,
)
    if family_slot == 1
        optic = @inbounds prepared_families[1].values[member_slot]
        optic_state = @inbounds state_families[1].values[member_slot]
        optic_workspace = @inbounds workspace_families[1].values[member_slot]
        stage_controllable_optic_command!(
            optic.implementation,
            optic_state,
            optic_workspace,
            command_endpoint_id(binding),
            _staged_effective_command(application_state),
            timestamp,
        )
        return nothing
    end
    return _stage_event_controllable_optic_family_command!(
        Base.tail(prepared_families),
        Base.tail(state_families),
        Base.tail(workspace_families),
        family_slot - 1,
        member_slot,
        binding,
        application_state,
        timestamp,
    )
end

@inline function _stage_event_controllable_optic_command!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer)
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    application_state =
        _event_command_application_state(state, endpoint_slot)
    slot = @inbounds prepared.optics.slots[Int(binding.optic_slot)]
    return _stage_event_controllable_optic_family_command!(
        prepared.optics.groups,
        state.controllable_optics.groups,
        workspace.controllable_optics.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        binding,
        application_state,
        state.scheduler.current_timestamp,
    )
end

@inline function _commit_event_controllable_optic_family_command!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    timestamp::PlantTimestamp,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _commit_event_controllable_optic_family_command!(
    prepared_families::Tuple,
    state_families::Tuple,
    workspace_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    timestamp::PlantTimestamp,
)
    if family_slot == 1
        optic = @inbounds prepared_families[1].values[member_slot]
        optic_state = @inbounds state_families[1].values[member_slot]
        optic_workspace = @inbounds workspace_families[1].values[member_slot]
        commit_controllable_optic_command!(
            optic.implementation,
            optic_state,
            optic_workspace,
            command_endpoint_id(binding),
            timestamp,
        )
        return nothing
    end
    return _commit_event_controllable_optic_family_command!(
        Base.tail(prepared_families),
        Base.tail(state_families),
        Base.tail(workspace_families),
        family_slot - 1,
        member_slot,
        binding,
        timestamp,
    )
end

@inline function _commit_event_controllable_optic_command!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    endpoint_slot::Integer)
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    slot = @inbounds prepared.optics.slots[Int(binding.optic_slot)]
    return _commit_event_controllable_optic_family_command!(
        prepared.optics.groups,
        state.controllable_optics.groups,
        workspace.controllable_optics.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        binding,
        state.scheduler.current_timestamp,
    )
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
    binding = _event_command_endpoint(prepared, endpoint_slot).binding
    endpoint_state = _event_command_endpoint_state(state, endpoint_slot)
    application_state =
        _event_command_application_state(state, endpoint_slot)
    endpoint_workspace = _event_command_workspace(workspace, endpoint_slot)
    slot = @inbounds prepared.optics.slots[Int(binding.optic_slot)]
    staged = try
        _apply_event_command_claim_family!(
            prepared.optics.groups,
            state.controllable_optics.groups,
            workspace.controllable_optics.groups,
            Int(slot.family_slot),
            Int(slot.member_slot),
            binding,
            endpoint_state,
            application_state,
            endpoint_workspace,
            claim,
            state.scheduler.current_timestamp,
        )
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

@inline function _apply_event_command_claim_family!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    endpoint_workspace::CommandDispositionWorkspace,
    claim::PlantCommandApplicationClaim,
    timestamp::PlantTimestamp,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _apply_event_command_claim_family!(
    prepared_families::Tuple,
    state_families::Tuple,
    workspace_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    binding::_PreparedPlantCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    endpoint_workspace::CommandDispositionWorkspace,
    claim::PlantCommandApplicationClaim,
    timestamp::PlantTimestamp,
)
    if family_slot == 1
        optic = @inbounds prepared_families[1].values[member_slot]
        optic_state = @inbounds state_families[1].values[member_slot]
        optic_workspace = @inbounds workspace_families[1].values[member_slot]
        return _apply_event_command_claim_parts!(
            binding,
            endpoint_state,
            application_state,
            endpoint_workspace,
            optic,
            optic_state,
            optic_workspace,
            claim,
            timestamp,
        )
    end
    return _apply_event_command_claim_family!(
        Base.tail(prepared_families),
        Base.tail(state_families),
        Base.tail(workspace_families),
        family_slot - 1,
        member_slot,
        binding,
        endpoint_state,
        application_state,
        endpoint_workspace,
        claim,
        timestamp,
    )
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
    endpoint = _event_command_endpoint(prepared, endpoint_slot).endpoint
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

@inline function _event_path_group(prepared::PreparedPlantEventLoop,
    slot::UInt32)
    index = Int(slot)
    1 <= index <= length(prepared.path_groups) ||
        _plant_event_loop_error(:invalid_action,
            "event action contains an invalid path-execution-group slot")
    return @inbounds prepared.path_groups[index]
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
    group = @inbounds prepared.path_groups[path_slot]
    _event_generator_due_at(prepared, state, group.handle, timestamp) ||
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

Base.@noinline function _complete_event_acquisition_readout!(
    product,
    lifecycle::AbstractPreparedAcquisitionLifecycle,
    acquisition_state::AbstractAcquisitionLifecycleState,
    timestamp::PlantTimestamp,
    rng::AbstractRNG,
)
    output = complete_readout!(
        lifecycle, acquisition_state, timestamp, rng)
    product === output ||
        copy_acquisition_product!(product, output)
    return nothing
end

function _process_acquisition_readout!(prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState, claim::EventClaim,
    action::_PlantEventAction)
    acquisition = _event_acquisition_binding(prepared, action.owner_slot)
    acquisition_state = _event_acquisition_state(state, action.owner_slot)
    _complete_event_acquisition_readout!(
        acquisition.product,
        acquisition.lifecycle,
        acquisition_state,
        claim.key.timestamp,
        acquisition.rng,
    )
    index = Int(action.owner_slot)
    sequence = _event_product_sequence(acquisition_state)
    previous_sequence = @inbounds state.product_sequences[index]
    sequence > previous_sequence ||
        _plant_event_loop_error(:product_sequence,
            "acquisition product sequence did not advance")
    @inbounds state.product_sequences[index] = sequence
    @inbounds state.product_ready_timestamps[index] =
        state.scheduler.current_timestamp
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
        _reset_event_autonomous_optic!(
            prepared,
            state,
            binding.optic_slot,
            binding.implementation,
            timestamp,
            nominal.sequence,
        )
    end
    return nothing
end

function _reset_triggered_autonomous_optic!(
    prepared::PreparedPlantEventLoop, state::PlantEventLoopState,
    slot::UInt32, delivery::TriggerDelivery)
    binding = @inbounds prepared.autonomous_optics[Int(slot)]
    return _reset_event_autonomous_optic!(
        prepared,
        state,
        binding.optic_slot,
        binding.implementation,
        delivered_trigger_edge(delivery).timestamp,
        nominal_trigger_edge(delivery).sequence,
    )
end

@inline function _reset_event_autonomous_optic_family!(
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    implementation,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _reset_event_autonomous_optic_family!(
    state_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    implementation,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    if family_slot == 1
        optic_state = @inbounds state_families[1].values[member_slot]
        reset_autonomous_periodic_optic_phase!(
            implementation, optic_state, timestamp, sequence)
        return nothing
    end
    return _reset_event_autonomous_optic_family!(
        Base.tail(state_families),
        family_slot - 1,
        member_slot,
        implementation,
        timestamp,
        sequence,
    )
end

@inline function _reset_event_autonomous_optic!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    optic_slot::UInt32,
    implementation,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    slot = @inbounds prepared.optics.slots[Int(optic_slot)]
    return _reset_event_autonomous_optic_family!(
        state.controllable_optics.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        implementation,
        timestamp,
        sequence,
    )
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
    prepared::PreparedPlantEventLoop{
        <:Any,<:Any,<:Any,<:Any,<:Any,<:PreparedTriggerTopology},
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
    ::PreparedPlantEventLoop{
        <:Any,<:Any,<:Any,<:Any,<:Any,<:_NoPreparedTriggerTopology},
    ::PlantEventLoopState{<:_NoTriggerTopologyState},
    ::PlantEventLoopWorkspace{<:_NoTriggerTopologyWorkspace}, ::EventClaim)
    _plant_event_loop_error(:invalid_action,
        "trigger action exists without a prepared trigger topology")
end

function _preflight_event_path(group::PreparedPathExecutionGroup)
    _preflight_event_path(
        group.path, group.rngs, group.atmosphere)
    return nothing
end

Base.@noinline function _preflight_event_path(
    path::PreparedPathExecutor,
    rngs::PreparedOwnerRNGs,
    atmosphere,
)
    path.atmosphere === atmosphere || _plant_event_loop_error(
        :prepared_binding,
        "event path does not retain the plant atmosphere")
    _require_current_path_binding(path)
    _require_rng_owner_binding(rngs, path)
    return nothing
end

function _preflight_event_acquisition(
    acquisition::_PreparedPlantEventAcquisition,
    state::_AcquisitionEventLifecycleState)
    _require_event_lifecycle_binding(acquisition.lifecycle, state)
    return nothing
end

@inline _require_event_lifecycle_binding(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState,
) = _require_detector_event_binding(prepared, state)

@inline _require_event_lifecycle_binding(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState,
) = _require_rolling_shutter_event_binding(prepared, state)

@inline _require_event_lifecycle_binding(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState,
) = _require_frame_transfer_event_binding(prepared, state)

function _require_event_lifecycle_binding(
    ::Union{
        PreparedGlobalShutterAcquisition,
        PreparedRollingShutterAcquisition,
        PreparedFrameTransferAcquisition,
    },
    ::Union{
        GlobalShutterAcquisitionState,
        RollingShutterAcquisitionState,
        FrameTransferAcquisitionState,
    },
)
    _plant_event_loop_error(
        :prepared_binding,
        "detector lifecycle state kind changed after preparation",
    )
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
        compute_device(storage) ==
            compute_device(prepared.instantaneous_sample) ||
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

@inline function _require_idle_optical_path_batch(
    workspace::PlantEventLoopWorkspace,
)
    workspace.optical_path_batch.phase == _OpticalPathBatchIdle ||
        _plant_event_loop_error(
            :optical_path_batch_active,
            "complete the active optical-path batch before another event-loop operation",
        )
    return nothing
end

@inline function _require_command_failure_batch_state(
    workspace::PlantEventLoopWorkspace)
    phase = workspace.optical_path_batch.phase
    phase in (
        _OpticalPathBatchIdle,
        _OpticalPathBatchAbandoned,
    ) || _plant_event_loop_error(
        :optical_path_batch_active,
        "complete or explicitly abandon the active optical-path batch before failing pending commands",
    )
    return nothing
end

@inline function _next_optical_path_batch_generation(
    generation::UInt64,
)
    generation != typemax(UInt64) || _plant_event_loop_error(
        :optical_path_batch_generation_overflow,
        "optical-path batch generation exceeds UInt64 range",
    )
    return generation + UInt64(1)
end

@inline _plant_event_loop_binding_id(prepared::PreparedPlantEventLoop) =
    UInt64(objectid(prepared.binding))
@inline _plant_event_loop_state_binding_id(state::PlantEventLoopState) =
    UInt64(objectid(state.state_binding))
@inline _plant_event_loop_workspace_binding_id(
    workspace::PlantEventLoopWorkspace) =
    UInt64(objectid(workspace.optical_path_batch.binding))

@inline function _require_current_optical_path_batch(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    batch = workspace.optical_path_batch
    claim.plant_binding_id == _plant_event_loop_binding_id(prepared) ||
        _plant_event_loop_error(
            :foreign_optical_path_batch_claim,
            "optical-path batch claim belongs to another prepared event loop",
        )
    claim.workspace_binding_id ==
        _plant_event_loop_workspace_binding_id(workspace) ||
        _plant_event_loop_error(
            :foreign_optical_path_batch_workspace,
            "optical-path batch claim belongs to another event-loop workspace",
        )
    claim.state_binding_id == _plant_event_loop_state_binding_id(state) &&
        batch.state_binding === state.state_binding ||
        _plant_event_loop_error(
            :foreign_optical_path_batch_state,
            "optical-path batch claim belongs to another event-loop state",
        )
    batch.phase != _OpticalPathBatchIdle &&
        batch.phase != _OpticalPathBatchAbandoned &&
        batch.generation == claim.generation &&
        batch.scheduler_revision == claim.scheduler_revision &&
        batch.timestamp == claim.timestamp &&
        batch.has_epoch == claim.has_epoch &&
        (!batch.has_epoch ||
            batch.epoch.sequence == claim.epoch_sequence) ||
        _plant_event_loop_error(
            :stale_optical_path_batch_claim,
            "optical-path batch claim is stale or no longer active",
        )
    state.scheduler.revision == claim.scheduler_revision ||
        _plant_event_loop_error(
            :optical_path_batch_scheduler_changed,
            "event scheduler changed while an optical-path batch was active",
        )
    return batch
end

"""
Return the batch's atmosphere epoch identity token, or `nothing` when no due
group uses a full optical provider. The active state and workspace are required
because the compact claim stores only checked owner and epoch identities; the
token remains in bounded workspace storage and does not retain atmosphere
layers.
"""
function optical_path_batch_epoch(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    return batch.has_epoch ? batch.epoch : nothing
end

@inline function _require_optical_path_batch_phase(
    batch::_OpticalPathBatchWorkspace,
    expected::_OpticalPathBatchPhase,
    operation::AbstractString,
)
    batch.phase == expected || _plant_event_loop_error(
        :invalid_optical_path_batch_phase,
        "$operation is invalid in optical-path batch phase $(batch.phase)",
    )
    return nothing
end

@inline function _require_due_path_execution_group(
    prepared::PreparedPlantEventLoop,
    batch::_OpticalPathBatchWorkspace,
    ordinal::Integer,
)
    1 <= ordinal <= length(prepared.path_groups) ||
        _plant_event_loop_error(
            :invalid_path_execution_group,
            "path execution-group ordinal must be within the prepared registry",
        )
    slot = Int(ordinal)
    @inbounds batch.group_status[slot] !=
        _OpticalPathBatchGroupNotDue || _plant_event_loop_error(
            :path_execution_group_not_due,
            "path execution group $ordinal is not due in the active optical-path batch",
        )
    return slot
end

@inline _require_due_path_execution_group(
    ::PreparedPlantEventLoop,
    ::_OpticalPathBatchWorkspace,
    ::Bool,
) = _plant_event_loop_error(
    :invalid_path_execution_group,
    "path execution-group ordinal must be an integer count, not Bool",
)

"""
    optical_path_batch_due_group_count(prepared, state, workspace, claim)

Return the fixed number of path execution groups due in the active batch.
"""
function optical_path_batch_due_group_count(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    return batch.due_group_count
end

"""
    optical_path_batch_due_group_ordinal(
        prepared, state, workspace, claim, index)

Return the canonical prepared-group ordinal at one position in the active
batch's bounded due-group registry.
"""
function optical_path_batch_due_group_ordinal(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    index::Integer,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    1 <= index <= batch.due_group_count || _plant_event_loop_error(
        :invalid_due_path_execution_group,
        "due path execution-group index must be within the active batch",
    )
    return Int(@inbounds batch.due_group_slots[Int(index)])
end

@inline optical_path_batch_due_group_ordinal(
    ::PreparedPlantEventLoop,
    ::PlantEventLoopState,
    ::PlantEventLoopWorkspace,
    ::OpticalPathBatchClaim,
    ::Bool,
) = _plant_event_loop_error(
    :invalid_due_path_execution_group,
    "due path execution-group index must be an integer count, not Bool",
)

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
    @inbounds for path_slot in eachindex(prepared.path_groups)
        due_paths[path_slot] || continue
        group = prepared.path_groups[path_slot]
        for member in group.acquisitions
            index = Int(member.slot)
            acquisition = prepared.acquisitions[index]
            acquisition_state = state.acquisitions[index]
            _preflight_event_acquisition(acquisition, acquisition_state)
            _preflight_event_integration_to(acquisition.lifecycle,
                acquisition_state, timestamp)
        end
    end
    return nothing
end

function _validate_due_path_materializations!(
    prepared::PreparedPlantEventLoop, due_paths::Memory{Bool}, atmosphere,
    epoch)
    @inbounds for index in eachindex(prepared.path_groups)
        due_paths[index] || continue
        group = prepared.path_groups[index]
        group.requirements.requires_full_optical || continue
        _validate_due_path_materialization!(
            group.path,
            atmosphere,
            epoch,
        )
    end
    return nothing
end

Base.@noinline function _validate_due_path_materialization!(
    path::PreparedPathExecutor,
    atmosphere,
    epoch,
)
    validate_path_materialization(
        path.materialization,
        path.input,
        atmosphere,
        epoch,
    )
    return nothing
end

Base.@noinline function _materialize_due_path!(
    path::PreparedPathExecutor,
    atmosphere,
    epoch,
    rngs::PreparedOwnerRNGs,
)
    materialize_path_input_rngs!(
        path.materialization,
        path.input,
        atmosphere,
        epoch,
        rngs,
    )
    return nothing
end

@inline function _require_next_optical_path_batch(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    count = scan_due_events!(
        workspace.scheduler, prepared.scheduler, state.scheduler)
    iszero(count) && _plant_event_loop_error(
        :optical_path_batch_not_due,
        "no plant event is due for optical-path batch execution",
    )
    workspace.scheduler.due_timestamp == timestamp ||
        _plant_event_loop_error(
            :optical_path_batch_timestamp_mismatch,
            "next plant event is due at $(workspace.scheduler.due_timestamp), not $timestamp",
        )
    key = due_event_key(
        workspace.scheduler, prepared.scheduler, state.scheduler, 1)
    key.phase == OpticalSamplePhase || _plant_event_loop_error(
        :optical_path_batch_not_due,
        "process earlier causal phases before beginning the optical-path batch at $timestamp",
    )
    return nothing
end

"""
    begin_optical_path_batch!(prepared, state, workspace, timestamp)

Open the exact due optical-path batch at `timestamp`. The operation performs a
mutation-free preflight, advances the shared atmosphere at most once, validates
every due full-optical materializer against that publication, and returns an
opaque bounded claim. It does not materialize or execute a path group.
"""
function begin_optical_path_batch!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    result = _with_completed_prepared_device_execution_context(
        prepared.context) do
        _begin_optical_path_batch_in_context!(
            prepared, state, workspace, timestamp)
    end
    return result::OpticalPathBatchClaim
end

function _begin_optical_path_batch_in_context!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    _require_idle_optical_path_batch(workspace)
    _require_next_optical_path_batch(
        prepared, state, workspace, timestamp)

    atmosphere = prepared.atmosphere
    _mark_due_event_paths!(prepared, state, workspace, timestamp)
    @inbounds for index in eachindex(prepared.path_groups)
        workspace.due_paths[index] || continue
        _preflight_event_path(prepared.path_groups[index])
    end
    _preflight_due_path_consumers(
        prepared, state, workspace.due_paths, timestamp)

    batch = workspace.optical_path_batch
    generation = _next_optical_path_batch_generation(batch.generation)
    has_epoch = _due_full_optical_path(prepared, workspace.due_paths)
    epoch = batch.epoch
    if has_epoch
        target_time = _preflight_atmosphere_time(atmosphere, timestamp)
        epoch = advance_to!(
            atmosphere, target_time, prepared.atmosphere_rng)
        _validate_due_path_materializations!(
            prepared, workspace.due_paths, atmosphere, epoch)
    end

    fill!(batch.group_status, _OpticalPathBatchGroupNotDue)
    due_group_count = 0
    @inbounds for index in eachindex(workspace.due_paths)
        workspace.due_paths[index] || continue
        due_group_count += 1
        batch.due_group_slots[due_group_count] = UInt32(index)
        batch.group_status[index] =
            _OpticalPathBatchGroupAwaitingMaterialization
    end
    batch.state_binding = state.state_binding
    batch.generation = generation
    batch.scheduler_revision = state.scheduler.revision
    batch.timestamp = timestamp
    batch.epoch = epoch
    batch.has_epoch = has_epoch
    batch.due_group_count = due_group_count
    batch.phase = _OpticalPathBatchMaterializing
    return OpticalPathBatchClaim(
        _plant_event_loop_binding_id(prepared),
        _plant_event_loop_state_binding_id(state),
        _plant_event_loop_workspace_binding_id(workspace),
        generation,
        batch.scheduler_revision,
        timestamp,
        has_epoch ? epoch.sequence : UInt64(0),
        has_epoch,
        _OPTICAL_PATH_BATCH_CLAIM_TOKEN,
    )
end

"""
    materialize_path_execution_group!(
        prepared, state, workspace, claim, ordinal)

Materialize one due group's path-local input exactly once against the batch
epoch. Reduced-order and replay groups cross the same bounded lifecycle without
reading the atmosphere. The atmosphere writer must remain held until every due
group has crossed this phase. A model exception after the group enters
materialization is fail-stop for this run because its product or RNG state may
already be partially mutated; core does not attempt rollback or retry.
"""
function materialize_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
)
    return _with_prepared_device_execution_context(
        prepared.context) do
        _materialize_path_execution_group_in_context!(
            prepared, state, workspace, claim, ordinal, Val(true))
    end
end

@inline _synchronize_independent_path_execution_group!(
    ::Val{false}, context) = nothing

@inline _synchronize_independent_path_execution_group!(
    ::Val{true}, context) =
    _synchronize_prepared_device_execution_context!(context)

function _materialize_path_execution_group_in_context!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
    synchronize::Val,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch, _OpticalPathBatchMaterializing,
        "path execution-group materialization")
    slot = _require_due_path_execution_group(prepared, batch, ordinal)
    _require_independent_path_execution_group(prepared, slot)
    status = @inbounds batch.group_status[slot]
    status == _OpticalPathBatchGroupAwaitingMaterialization || begin
        reason = status == _OpticalPathBatchGroupMaterializing ?
            :path_execution_group_materialization_active :
            :duplicate_path_execution_group_materialization
        _plant_event_loop_error(
            reason,
            "path execution group $slot has already entered materialization",
        )
    end
    @inbounds batch.group_status[slot] =
        _OpticalPathBatchGroupMaterializing
    group = @inbounds prepared.path_groups[slot]
    if group.requirements.requires_full_optical
        _validate_epoch_identity(
            atmosphere_identity(prepared.atmosphere),
            prepared.atmosphere,
            batch.epoch,
        )
        _materialize_due_path!(
            group.path, prepared.atmosphere, batch.epoch, group.rngs)
    end
    _synchronize_independent_path_execution_group!(
        synchronize, prepared.context)
    @inbounds batch.group_status[slot] = _OpticalPathBatchGroupReady
    return nothing
end

"""
    seal_optical_path_batch_materialization!(
        prepared, state, workspace, claim)

Close the atmosphere-read phase only after every due group has completed its
materialization transition. Once sealed, group execution no longer reads
mutable atmosphere layers.
"""
function seal_optical_path_batch_materialization!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch, _OpticalPathBatchMaterializing,
        "optical-path materialization sealing")
    @inbounds for index in 1:batch.due_group_count
        slot = Int(batch.due_group_slots[index])
        batch.group_status[slot] == _OpticalPathBatchGroupReady ||
            _plant_event_loop_error(
                :incomplete_optical_path_batch_materialization,
                "path execution group $slot has not completed materialization",
            )
    end
    batch.phase = _OpticalPathBatchExecuting
    return nothing
end

Base.@noinline function _apply_due_path_sampled_aberrations!(
    path::PreparedPathExecutor,
    aberrations,
    bindings::PreparedSampledAberrationPathBindings,
    binding_range::UnitRange{Int},
)
    _apply_sampled_aberration_bindings_noreturn!(
        path.input,
        aberrations,
        bindings,
        binding_range,
    )
    return nothing
end

Base.@noinline function _apply_due_path_controllable_optics!(
    input,
    binding_start::Int,
    binding_stop::Int,
    couplings::Memory{AbstractPupilSurfacePathCoupling},
    optics::_PreparedControllableOpticRegistry,
    states::_ControllableOpticStateRegistry,
    bindings::PreparedControllableOpticPathBindings,
)
    Base.@nospecialize input
    coupling_slot = 1
    @inbounds for binding in binding_start:binding_stop
        optic_index = prepared_controllable_optic_slot(bindings, binding)
        _apply_prepared_event_controllable_optic_slot!(
            input,
            optics,
            states,
            optic_index,
            couplings[coupling_slot],
        )
        coupling_slot += 1
    end
    return nothing
end

@inline function _apply_prepared_event_controllable_optic_family!(
    input,
    ::Tuple{},
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    coupling::AbstractPupilSurfacePathCoupling,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _apply_prepared_event_controllable_optic_family!(
    input,
    prepared_families::Tuple,
    state_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    coupling::AbstractPupilSurfacePathCoupling,
)
    if family_slot == 1
        optic = @inbounds prepared_families[1].values[member_slot]
        state = @inbounds state_families[1].values[member_slot]
        return _apply_prepared_event_controllable_optic_surface!(
            input, optic, state, coupling)
    end
    return _apply_prepared_event_controllable_optic_family!(
        input,
        Base.tail(prepared_families),
        Base.tail(state_families),
        family_slot - 1,
        member_slot,
        coupling,
    )
end

@inline function _apply_prepared_event_controllable_optic_slot!(
    input,
    optics::_PreparedControllableOpticRegistry,
    states::_ControllableOpticStateRegistry,
    optic_index::Int,
    coupling::AbstractPupilSurfacePathCoupling,
)
    slot = @inbounds optics.slots[optic_index]
    return _apply_prepared_event_controllable_optic_family!(
        input,
        optics.groups,
        states.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        coupling,
    )
end

@inline function _apply_prepared_event_controllable_optic_surface!(
    input,
    optic::PreparedControllableOptic,
    state,
    coupling::AbstractPupilSurfacePathCoupling,
)
    implementation = optic.implementation
    _apply_event_controllable_optic_surface!(
        controllable_optic_execution_role(implementation),
        input,
        implementation,
        state,
        coupling,
    )
    return nothing
end

@inline function _apply_event_controllable_optic_surface!(
    ::PupilSurfaceExecutionRole, input, implementation, state,
    coupling::AbstractPupilSurfacePathCoupling)
    apply_controllable_optic_surface!(
        input, implementation, state, coupling)
    return nothing
end

@inline function _apply_event_controllable_optic_surface!(
    ::AutonomousPathExecutionRole, input, implementation, state,
    ::_NoPupilSurfacePathCoupling)
    return nothing
end

Base.@noinline function _execute_due_path!(
    path::PreparedPathExecutor,
    rngs::PreparedOwnerRNGs,
)
    _execute_path_in_context!(path, rngs)
    return nothing
end

Base.@noinline function _stage_materialized_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    group::PreparedPathExecutionGroup,
    timestamp::PlantTimestamp,
)
    for member in group.acquisitions
        index = Int(member.slot)
        acquisition = @inbounds prepared.acquisitions[index]
        acquisition_state = @inbounds state.acquisitions[index]
        _integrate_event_acquisition_to!(
            acquisition.lifecycle,
            acquisition_state,
            timestamp,
            acquisition.rng,
        )
    end

    if group.requirements.requires_full_optical
        sampled_binding_range = (
            group.sampled_aberration_binding_start:
            group.sampled_aberration_binding_stop
        )
        _apply_due_path_sampled_aberrations!(
            group.path,
            prepared.sampled_aberrations,
            prepared.sampled_aberration_path_bindings,
            sampled_binding_range,
        )
        _apply_due_path_controllable_optics!(
            group.input,
            group.optic_binding_start,
            group.optic_binding_stop,
            group.optic_couplings,
            prepared.optics,
            state.controllable_optics,
            prepared.optic_path_bindings,
        )
        @inbounds for binding_slot in group.autonomous_optic_slots
            binding = prepared.autonomous_optics[Int(binding_slot)]
            _evaluate_event_autonomous_optic!(
                prepared,
                state,
                binding.optic_slot,
                binding.implementation,
                binding.coupling,
                timestamp,
            )
        end
    end
    return nothing
end

@inline function _evaluate_event_autonomous_optic_family!(
    ::Tuple{},
    family_slot::Int,
    member_slot::Int,
    implementation,
    coupling,
    timestamp::PlantTimestamp,
)
    return _prepared_controllable_optic_slot_error(
        family_slot, member_slot)
end

@inline function _evaluate_event_autonomous_optic_family!(
    state_families::Tuple,
    family_slot::Int,
    member_slot::Int,
    implementation,
    coupling,
    timestamp::PlantTimestamp,
)
    if family_slot == 1
        optic_state = @inbounds state_families[1].values[member_slot]
        evaluate_autonomous_periodic_optic!(
            implementation, optic_state, coupling, timestamp)
        return nothing
    end
    return _evaluate_event_autonomous_optic_family!(
        Base.tail(state_families),
        family_slot - 1,
        member_slot,
        implementation,
        coupling,
        timestamp,
    )
end

@inline function _evaluate_event_autonomous_optic!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    optic_slot::UInt32,
    implementation,
    coupling,
    timestamp::PlantTimestamp,
)
    slot = @inbounds prepared.optics.slots[Int(optic_slot)]
    return _evaluate_event_autonomous_optic_family!(
        state.controllable_optics.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
        implementation,
        coupling,
        timestamp,
    )
end

Base.@noinline function _finish_materialized_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    group::PreparedPathExecutionGroup,
    timestamp::PlantTimestamp,
)
    for member in group.acquisitions
        acquisition = @inbounds prepared.acquisitions[Int(member.slot)]
        _evaluate_event_sample!(
            acquisition.sample_provider,
            acquisition.lifecycle,
            timestamp,
            state.command_applications,
        )
    end
    return nothing
end

Base.@noinline function _execute_materialized_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    group::PreparedPathExecutionGroup,
    timestamp::PlantTimestamp,
)
    _stage_materialized_path_execution_group!(
        prepared, state, group, timestamp)
    group.requirements.requires_full_optical &&
        _execute_due_path!(group.path, group.rngs)
    _finish_materialized_path_execution_group!(
        prepared, state, group, timestamp)
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

function _due_full_optical_path(
    prepared::PreparedPlantEventLoop, due_paths::Memory{Bool})
    @inbounds for index in eachindex(prepared.path_groups)
        due_paths[index] &&
            prepared.path_groups[index].requirements.requires_full_optical &&
            return true
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
        group = _event_path_group(prepared, action.owner_slot)
        reschedule_periodic_event!(prepared.scheduler, state.scheduler,
            claim, group.schedule; origin=group.origin)
    end
end

"""
    execute_path_execution_group!(
        prepared, state, workspace, claim, ordinal)

Execute one materialized due path group exactly once. The call mutates only
that group's path-local products, acquisition owners, RNG streams, and
path-local autonomous-optic state while borrowing held effective commands
read-only. A model exception after execution begins is fail-stop for this run;
core does not retry partially mutated products, acquisition state, or RNG
streams.
"""
function execute_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
)
    return _with_prepared_device_execution_context(
        prepared.context) do
        _execute_path_execution_group_in_context!(
            prepared, state, workspace, claim, ordinal, Val(true))
    end
end

function _execute_path_execution_group_in_context!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
    synchronize::Val,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch, _OpticalPathBatchExecuting,
        "path execution-group execution")
    slot = _require_due_path_execution_group(prepared, batch, ordinal)
    _require_independent_path_execution_group(prepared, slot)
    status = @inbounds batch.group_status[slot]
    status == _OpticalPathBatchGroupReady || begin
        reason = status == _OpticalPathBatchGroupComplete ?
            :duplicate_path_execution_group :
            status == _OpticalPathBatchGroupExecuting ?
                :path_execution_group_active :
                :path_execution_group_not_materialized
        _plant_event_loop_error(
            reason,
            "path execution group $slot is not ready for execution",
        )
    end

    @inbounds batch.group_status[slot] = _OpticalPathBatchGroupExecuting
    group = @inbounds prepared.path_groups[slot]
    _execute_materialized_path_execution_group!(
        prepared, state, group, claim.timestamp)
    _synchronize_independent_path_execution_group!(
        synchronize, prepared.context)
    @inbounds state.path_sampled[slot] = true
    @inbounds batch.group_status[slot] = _OpticalPathBatchGroupComplete
    return nothing
end

"""
    complete_optical_path_batch!(prepared, state, workspace, claim)

Resolve every due optical-sample generator and reclaim the bounded batch
workspace after all due groups complete. Incomplete or copied stale claims fail
before scheduler state changes.
"""
function complete_optical_path_batch!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch, _OpticalPathBatchExecuting,
        "optical-path batch completion")
    @inbounds for index in 1:batch.due_group_count
        slot = Int(batch.due_group_slots[index])
        batch.group_status[slot] == _OpticalPathBatchGroupComplete ||
            _plant_event_loop_error(
                :incomplete_optical_path_batch_execution,
                "path execution group $slot has not completed execution",
            )
    end

    batch.phase = _OpticalPathBatchFinalizing
    _resolve_due_path_claims!(
        prepared, state, workspace, claim.timestamp)
    fill!(workspace.due_paths, false)
    fill!(batch.group_status, _OpticalPathBatchGroupNotDue)
    batch.due_group_count = 0
    batch.has_epoch = false
    batch.scheduler_revision = state.scheduler.revision
    batch.phase = _OpticalPathBatchIdle
    return claim.timestamp
end

"""
    abandon_optical_path_batch!(prepared, state, workspace, claim)

Relinquish the coordinator-local ownership of a failed optical-path batch
without resolving its scheduler claims or attempting to undo any partially
mutated plant, product, acquisition, atmosphere, or RNG state. The event loop
is permanently fail-stopped after abandonment: only bounded terminal failure
of pending commands remains permitted, and further admission or plant
execution requires a fresh prepared run.
"""
function abandon_optical_path_batch!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    fill!(workspace.due_paths, false)
    fill!(batch.group_status, _OpticalPathBatchGroupNotDue)
    batch.due_group_count = 0
    batch.has_epoch = false
    batch.phase = _OpticalPathBatchAbandoned
    return claim.timestamp
end

function execute_optical_path_batch!(
    executor::AbstractOpticalPathBatchExecutor,
    ::PreparedPlantEventLoop,
    ::PlantEventLoopState,
    ::PlantEventLoopWorkspace,
    ::PlantTimestamp,
)
    _plant_event_loop_error(
        :unsupported_optical_path_batch_executor,
        "optical-path batch executor $(typeof(executor)) does not implement execute_optical_path_batch!",
    )
end

"""
Run one optical-path batch through the canonical serial order using the same
public materialization, execution, and completion lifecycle available to HIL
execution policies.
"""
function execute_optical_path_batch!(
    executor::SerialOpticalPathBatchExecutor,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    result = _with_completed_prepared_device_execution_context(
        prepared.context) do
        _execute_serial_optical_path_batch_in_context!(
            executor, prepared, state, workspace, timestamp)
    end
    return result::PlantTimestamp
end

function _execute_serial_optical_path_batch_in_context!(
    ::SerialOpticalPathBatchExecutor,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    claim = _begin_optical_path_batch_in_context!(
        prepared, state, workspace, timestamp)
    try
        count = optical_path_batch_due_group_count(
            prepared, state, workspace, claim)
        @inbounds for index in 1:count
            ordinal = optical_path_batch_due_group_ordinal(
                prepared, state, workspace, claim, index)
            _materialize_due_path_execution_group!(
                prepared, state, workspace, claim, ordinal)
        end
        seal_optical_path_batch_materialization!(
            prepared, state, workspace, claim)
        @inbounds for index in 1:count
            ordinal = optical_path_batch_due_group_ordinal(
                prepared, state, workspace, claim, index)
            _execute_due_path_execution_group!(
                prepared, state, workspace, claim, ordinal)
        end
    catch
        abandon_optical_path_batch!(
            prepared, state, workspace, claim)
        rethrow()
    end
    return complete_optical_path_batch!(
        prepared, state, workspace, claim)
end

@inline function _execute_optical_path_batch_in_context!(
    executor::AbstractOpticalPathBatchExecutor,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    return execute_optical_path_batch!(
        executor, prepared, state, workspace, timestamp)
end

@inline function _execute_optical_path_batch_in_context!(
    executor::SerialOpticalPathBatchExecutor,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    return _execute_serial_optical_path_batch_in_context!(
        executor, prepared, state, workspace, timestamp)
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
    _require_idle_optical_path_batch(workspace)
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
    return step_plant_events!(
        prepared,
        state,
        workspace,
        SerialOpticalPathBatchExecutor(),
    )
end

"""
    step_plant_events!(prepared, state, workspace, batch_executor)

Process the next canonical timestamp using `batch_executor` only for the
optical-sample phase. The executor must complete the bounded batch before
returning; all other causal phases remain owned by the serial event
coordinator.
"""
function step_plant_events!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    batch_executor::AbstractOpticalPathBatchExecutor,
)
    result = _with_completed_prepared_device_execution_context(
        prepared.context) do
        _step_plant_events_in_context!(
            prepared, state, workspace, batch_executor)
    end
    return result::Union{Nothing,PlantTimestamp}
end

function _step_plant_events_in_context!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    batch_executor::AbstractOpticalPathBatchExecutor,
)
    _require_plant_event_loop_binding(prepared, state)
    _require_plant_event_loop_binding(prepared, workspace)
    _require_idle_optical_path_batch(workspace)
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
            _execute_optical_path_batch_in_context!(
                batch_executor, prepared, state, workspace, timestamp)
            _require_idle_optical_path_batch(workspace)
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
    return run_plant_events_until!(
        prepared,
        state,
        workspace,
        stop,
        SerialOpticalPathBatchExecutor();
        max_timestamps,
    )
end

"""
Process a finite plant horizon while delegating every optical-sample batch to
one explicit execution policy.
"""
function run_plant_events_until!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    stop::PlantTimestamp,
    batch_executor::AbstractOpticalPathBatchExecutor;
    max_timestamps::Integer=typemax(Int),
)
    result = _with_completed_prepared_device_execution_context(
        prepared.context) do
        _run_plant_events_until_in_context!(
            prepared, state, workspace, stop, batch_executor;
            max_timestamps)
    end
    return result::Int
end

function _run_plant_events_until_in_context!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    stop::PlantTimestamp,
    batch_executor::AbstractOpticalPathBatchExecutor;
    max_timestamps::Integer=typemax(Int),
)
    limit = _checked_event_step_limit(max_timestamps)
    count = 0
    while count < limit
        timestamp = next_plant_event_timestamp(prepared, state, workspace)
        timestamp === nothing && break
        timestamp <= stop || break
        _step_plant_events_in_context!(
            prepared, state, workspace, batch_executor)
        count += 1
    end
    return count
end
