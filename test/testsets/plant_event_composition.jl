const AOSEvents = AdaptiveOpticsSim

function event_test_intensity_map(values::AbstractMatrix{T}) where
    {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.75e-6)))
    return IntensityMap(metadata, values)
end

function run_rolling_event_frame!(prepared, state, rng,
    start::PlantTimestamp)
    begin_exposure!(prepared, state, start)
    while true
        next_open = next_rolling_band_open_timestamp(prepared, state)
        next_close = next_rolling_band_close_timestamp(prepared, state)
        next_open === nothing && next_close === nothing && break
        timestamp = next_open === nothing ? next_close :
            next_close === nothing ? next_open : min(next_open, next_close)
        if integrated_through_timestamp(state) < timestamp &&
                rolling_opened_band_count(state) >
                    rolling_closed_band_count(state)
            accumulate_rolling_exposure_interval!(prepared, state,
                integrated_through_timestamp(state), timestamp, rng)
        end
        next_close == timestamp &&
            close_next_rolling_band!(prepared, state, timestamp)
        next_open == timestamp &&
            open_next_rolling_band!(prepared, state, timestamp)
    end
    output = complete_readout!(prepared, state,
        readout_complete_timestamp(state), rng)
    mark_acquisition_ready!(prepared, state,
        acquisition_readiness_timestamp(state))
    return output
end

function run_frame_transfer_event!(prepared, state, rng,
    start::PlantTimestamp)
    close = start + prepared.definition.exposure_duration
    transfer = close + prepared.transfer_duration
    readout = transfer + prepared.definition.readout_duration
    begin_exposure!(prepared, state, start)
    accumulate_exposure_interval!(prepared, state, start, close, rng)
    close_exposure!(prepared, state, close)
    complete_frame_transfer!(prepared, state, transfer)
    return complete_readout!(prepared, state, readout, rng)
end

function event_test_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

mutable struct EventCompositionConcurrencyProbe
    expected::Int
    enabled::Threads.Atomic{Bool}
    arrivals::Threads.Atomic{Int}
    release::Threads.Atomic{Bool}
    thread_ids::Memory{Int}
    timeout_ns::UInt64
end

function EventCompositionConcurrencyProbe(
    expected::Integer;
    timeout_seconds::Real=5,
)
    expected > 1 || throw(ArgumentError(
        "concurrency probe requires at least two arrivals"))
    timeout_seconds > 0 || throw(ArgumentError(
        "concurrency probe timeout must be positive"))
    thread_ids = Memory{Int}(undef, Int(expected))
    fill!(thread_ids, 0)
    return EventCompositionConcurrencyProbe(
        Int(expected),
        Threads.Atomic{Bool}(true),
        Threads.Atomic{Int}(0),
        Threads.Atomic{Bool}(false),
        thread_ids,
        UInt64(round(Int, timeout_seconds * 1.0e9)),
    )
end

@inline event_composition_execution_probe!(::Nothing) = nothing

struct EventCompositionFailureProbe
    error::Exception
end

event_composition_execution_probe!(
    probe::EventCompositionFailureProbe) = throw(probe.error)

function event_composition_execution_probe!(
    probe::EventCompositionConcurrencyProbe,
)
    probe.enabled[] || return nothing
    slot = Threads.atomic_add!(probe.arrivals, 1) + 1
    slot <= probe.expected || return nothing
    @inbounds probe.thread_ids[slot] = Threads.threadid()
    if slot == probe.expected
        probe.enabled[] = false
        probe.release[] = true
        return nothing
    end
    start = time_ns()
    while !probe.release[]
        GC.safepoint()
        time_ns() - start <= probe.timeout_ns || error(
            "path groups did not overlap within the concurrency-probe timeout")
    end
    return nothing
end

struct EventCompositionPathModel{R,P}
    zero_padding::Int
    revision::R
    execution_probe::P
end

EventCompositionPathModel(zero_padding::Int, revision) =
    EventCompositionPathModel(zero_padding, revision, nothing)

struct EventCompositionPathExecution{E,C,P}
    imaging::E
    executions::C
    execution_probe::P
end

abstract type EventCompositionSensorKind end
struct EventGlobalCMOS <: EventCompositionSensorKind end
struct EventRollingCMOS <: EventCompositionSensorKind end
struct EventCCD <: EventCompositionSensorKind end
struct EventFrameTransferEMCCD <: EventCompositionSensorKind end
struct EventHgCdTeRamp <: EventCompositionSensorKind end

struct EventCompositionAcquisitionModel{
    T<:AbstractFloat,K<:EventCompositionSensorKind}
    exposure::T
    kind::K
end

Plant.plant_model_definition_style(
    ::Type{<:EventCompositionPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:EventCompositionAcquisitionModel}) = ColdPlantModelDefinition()

function Plant.validate_path_execution_binding(
    execution::EventCompositionPathExecution, input, result)
    return Plant.validate_path_execution_binding(
        execution.imaging, input, result)
end

function Plant.execute_path!(result, input,
    execution::EventCompositionPathExecution)
    Plant.validate_path_execution_binding(execution, input,
        result)
    event_composition_execution_probe!(execution.execution_probe)
    execution.executions[] += 1
    return Plant.execute_path!(result, input,
        execution.imaging)
end

function Plant.prepare_path_executor(
    model::EventCompositionPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    execution = EventCompositionPathExecution(
        imaging, Ref(0), model.execution_probe)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        execution;
        materialization=prepare_pupil_opd_materialization(atmosphere,
            telescope, source, pupil),
        optical_model=(kind=:event_composition_direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=model.revision,
    )
end

@inline event_composition_sensor(::EventGlobalCMOS, ::Type{T}) where {T} =
    CMOSSensor(timing_model=GlobalShutter(), T=T)
@inline event_composition_sensor(::EventRollingCMOS, ::Type{T}) where {T} =
    CMOSSensor(timing_model=RollingShutter(T(0.005)), T=T)
@inline event_composition_sensor(::EventCCD, ::Type{T}) where {T} =
    CCDSensor(T=T)
@inline event_composition_sensor(::EventFrameTransferEMCCD,
    ::Type{T}) where {T} = EMCCDSensor(
    acquisition_mode=FrameTransferAcquisition(transfer_time=T(0.02), T=T),
    T=T)
@inline event_composition_sensor(::EventHgCdTeRamp, ::Type{T}) where {T} =
    HgCdTeSensor(sampling_mode=UpTheRampSampling(3),
        read_time=zero(T), T=T)

function Plant.prepare_acquisition_provider(
    model::EventCompositionAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), gain=one(T),
        response_model=NullFrameResponse(),
        sensor=event_composition_sensor(model.kind, T), T=T,
        backend=path.key.backend)
    execution = FrameAcquisitionExecution(detector, path.result)
    metadata = (kind=:event_composition_frame,
        units=:detected_electrons, geometry=path.result.metadata,
        detector=detector_export_metadata(detector))
    products = AcquisitionProducts(execution.observation; metadata)
    return prepare_full_optical_provider(execution, products)
end

function event_composition_fixture(; reverse_order::Bool=false,
    unbound_trigger_consumer::Bool=false,
    faulted_trigger_fanout::Bool=false,
    missing_path_schedule::Bool=false,
    aligned_optical_samples::Bool=false,
    execution_probe=nothing)
    T = Float64
    telescope = Telescope(resolution=4, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[1], wind_speed=T[7], wind_direction=T[35],
        altitude=T[0], layer_ids=(:ground,), T=T)
    science_source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(80), coordinates=(T(0), T(0)), T=T)
    ngs_source = Source(band=:custom, wavelength=T(0.7e-6),
        photon_irradiance=T(65), coordinates=(T(2), T(35)), T=T)
    lgs_source = LGSSource(wavelength=T(589e-9),
        photon_irradiance=T(70), coordinates=(T(-3), T(80)),
        altitude=T(90_000), T=T)
    paths = (
        OpticalPathDefinition(:science, science_source,
            EventCompositionPathModel(1, UInt(1), execution_probe)),
        OpticalPathDefinition(:ngs, ngs_source,
            EventCompositionPathModel(1, UInt(2), execution_probe)),
        OpticalPathDefinition(:lgs, lgs_source,
            EventCompositionPathModel(1, UInt(3), execution_probe)),
    )
    acquisitions = (
        AcquisitionDefinition(:science_cmos, :science,
            EventCompositionAcquisitionModel(T(0.2), EventGlobalCMOS())),
        AcquisitionDefinition(:science_ccd, :science,
            EventCompositionAcquisitionModel(T(0.3), EventCCD())),
        AcquisitionDefinition(:science_rolling, :science,
            EventCompositionAcquisitionModel(T(0.15), EventRollingCMOS())),
        AcquisitionDefinition(:ngs_saphira, :ngs,
            EventCompositionAcquisitionModel(T(0.2), EventHgCdTeRamp())),
        AcquisitionDefinition(:lgs_emccd, :lgs,
            EventCompositionAcquisitionModel(T(0.2),
                EventFrameTransferEMCCD())),
    )
    plant = prepare_plant(PlantDefinition(; telescope, atmosphere,
        paths, acquisitions);
        run_seed=0x7900)

    trigger_faults = faulted_trigger_fanout ? TriggerFaultTrace(
        TriggerFaultTraceEntry(2, :common_trigger_drop;
            action=DropTriggerEdge),
        TriggerFaultTraceEntry(3, :common_trigger_phase_error;
            phase_step=PlantTimeOffset(20_000_000),
            jitter=PlantTimeOffset(10_000_000),
            timestamp_label_offset=PlantTimeOffset(7_000_000)),
        TriggerFaultTraceEntry(4, :common_trigger_duplicate;
            action=DuplicateTriggerEdge,
            duplicate_delay=PlantDuration(250_000_000)),
    ) : TriggerFaultTrace()
    trigger_source = TriggerSourceDefinition(:ngs_camera_trigger,
        PeriodicSchedule(
            period_ns=faulted_trigger_fanout ? 600_000_000 : 400_000_000,
            phase_ns=50_000_000); faults=trigger_faults)
    trigger_consumer = TriggerConsumerDefinition(:ngs_camera,
        TriggerSourceID(:ngs_camera_trigger))
    trigger_consumers = faulted_trigger_fanout ? (
        trigger_consumer,
        TriggerConsumerDefinition(:science_camera,
            TriggerSourceID(:ngs_camera_trigger)),
    ) : unbound_trigger_consumer ? (
        trigger_consumer,
        TriggerConsumerDefinition(:unbound_camera,
            TriggerSourceID(:ngs_camera_trigger)),
    ) : (trigger_consumer,)
    trigger_topology = prepare_trigger_topology((trigger_source,), (),
        trigger_consumers;
        in_flight_capacity=(unbound_trigger_consumer ||
            faulted_trigger_fanout) ? 8 : 1)
    samples = (
        OpticalSampleDefinition(:science,
            PeriodicSchedule(period_ns=100_000_000, phase_ns=0)),
        OpticalSampleDefinition(:ngs,
            PeriodicSchedule(period_ns=125_000_000,
                phase_ns=aligned_optical_samples ? 0 : 50_000_000)),
        OpticalSampleDefinition(:lgs,
            PeriodicSchedule(period_ns=200_000_000,
                phase_ns=aligned_optical_samples ? 0 : 25_000_000)),
    )
    events = (
        AcquisitionEventDefinition(:science_cmos,
            GlobalShutterAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(20_000_000),
                readiness_delay=PlantDuration(10_000_000)),
            faulted_trigger_fanout ?
                TriggeredAcquisitionStart(:science_camera) :
                PeriodicAcquisitionStart(PeriodicSchedule(
                    period_ns=500_000_000,
                    phase_ns=aligned_optical_samples ? 1 : 0))),
        AcquisitionEventDefinition(:science_ccd,
            GlobalShutterAcquisitionDefinition(PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=600_000_000, phase_ns=50_000_000))),
        AcquisitionEventDefinition(:science_rolling,
            RollingShutterAcquisitionDefinition(PlantDuration(150_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=500_000_000, phase_ns=250_000_000))),
        AcquisitionEventDefinition(:ngs_saphira,
            GlobalShutterAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(10_000_000)),
            TriggeredAcquisitionStart(:ngs_camera)),
        AcquisitionEventDefinition(:lgs_emccd,
            FrameTransferAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=350_000_000, phase_ns=25_000_000))),
    )
    selected_samples = missing_path_schedule ? Base.front(samples) : samples
    ordered_samples = reverse_order ? reverse(selected_samples) :
        selected_samples
    ordered_events = reverse_order ? reverse(events) : events
    definition = PlantEventLoopDefinition(ordered_samples, ordered_events;
        trigger_topology)
    prepared = prepare_plant_event_loop(plant, definition)
    return plant, prepared, PlantEventLoopState(prepared),
        PlantEventLoopWorkspace(prepared)
end

struct ReverseOpticalPathBatchExecutor <:
    Plant.AbstractOpticalPathBatchExecutor end

function Plant.execute_optical_path_batch!(
    ::ReverseOpticalPathBatchExecutor,
    prepared::Plant.PreparedPlantEventLoop,
    state::Plant.PlantEventLoopState,
    workspace::Plant.PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    claim = Plant.begin_optical_path_batch!(
        prepared, state, workspace, timestamp)
    count = Plant.optical_path_batch_due_group_count(
        prepared, state, workspace, claim)
    @inbounds for index in count:-1:1
        ordinal = Plant.optical_path_batch_due_group_ordinal(
            prepared, state, workspace, claim, index)
        Plant.materialize_path_execution_group!(
            prepared, state, workspace, claim, ordinal)
    end
    Plant.seal_optical_path_batch_materialization!(
        prepared, state, workspace, claim)
    @inbounds for index in count:-1:1
        ordinal = Plant.optical_path_batch_due_group_ordinal(
            prepared, state, workspace, claim, index)
        Plant.execute_path_execution_group!(
            prepared, state, workspace, claim, ordinal)
    end
    return Plant.complete_optical_path_batch!(
        prepared, state, workspace, claim)
end

@enum EventCompositionWorkerPhase::UInt8 begin
    EventCompositionMaterializationPhase = 0x01
    EventCompositionExecutionPhase = 0x02
end

struct EventCompositionWorkerCommand
    phase::EventCompositionWorkerPhase
    claim::Plant.OpticalPathBatchClaim
    measure_allocations::Bool
end

struct EventCompositionWorkerCompletion
    phase::EventCompositionWorkerPhase
    exception::Any
    allocated_bytes::Int
    thread_id::Int
    task_id::UInt
end

"""
Test-only fixed-owner executor used to prove that core's independently callable
group seams can be driven concurrently. The `Channel`s are synchronization for
this test harness; they are not the bounded SPSC transport proposed for the HIL
runtime.
"""
mutable struct FixedOwnerOpticalPathBatchExecutor{P,S,W} <:
    Plant.AbstractOpticalPathBatchExecutor
    prepared::P
    state::S
    workspace::W
    commands::Memory{
        Channel{Union{Nothing,EventCompositionWorkerCommand}}}
    completions::Memory{Channel{EventCompositionWorkerCompletion}}
    tasks::Memory{Task}
    due_group_ordinals::Memory{Int}
    task_ids::Memory{UInt}
    thread_ids::Memory{Int}
    batch_count::Int
    forward_phase_count::Int
    reverse_phase_count::Int
    materialization_count::Int
    execution_count::Int
    measure_allocations::Bool
    measured_group_call_count::Int
    measured_group_allocation_bytes::Int
    maximum_group_allocation_bytes::Int
    closed::Bool
end

function execute_event_composition_worker_command!(
    ordinal::Int,
    prepared,
    state,
    workspace,
    command::EventCompositionWorkerCommand,
)
    if command.phase == EventCompositionMaterializationPhase
        Plant.materialize_path_execution_group!(
            prepared, state, workspace, command.claim, ordinal)
    else
        Plant.execute_path_execution_group!(
            prepared, state, workspace, command.claim, ordinal)
    end
    return nothing
end

function event_composition_group_worker!(
    ordinal::Int,
    prepared,
    state,
    workspace,
    commands::Channel{Union{Nothing,EventCompositionWorkerCommand}},
    completions::Channel{EventCompositionWorkerCompletion},
)
    task_id = objectid(current_task())
    while true
        command = take!(commands)
        command === nothing && break
        allocated_bytes = 0
        exception = try
            if command.measure_allocations
                allocated_bytes = @allocated(
                    execute_event_composition_worker_command!(
                        ordinal, prepared, state, workspace, command))
            else
                execute_event_composition_worker_command!(
                    ordinal, prepared, state, workspace, command)
            end
            nothing
        catch caught
            caught
        end
        put!(completions, EventCompositionWorkerCompletion(
            command.phase,
            exception,
            allocated_bytes,
            Threads.threadid(),
            task_id,
        ))
    end
    return nothing
end

function FixedOwnerOpticalPathBatchExecutor(
    prepared::Plant.PreparedPlantEventLoop,
    state::Plant.PlantEventLoopState,
    workspace::Plant.PlantEventLoopWorkspace,
)
    count = path_execution_group_count(prepared)
    commands = Memory{
        Channel{Union{Nothing,EventCompositionWorkerCommand}}}(undef, count)
    completions = Memory{
        Channel{EventCompositionWorkerCompletion}}(undef, count)
    tasks = Memory{Task}(undef, count)
    @inbounds for ordinal in 1:count
        command = Channel{
            Union{Nothing,EventCompositionWorkerCommand}}(1)
        completion = Channel{EventCompositionWorkerCompletion}(1)
        commands[ordinal] = command
        completions[ordinal] = completion
        tasks[ordinal] = Threads.@spawn event_composition_group_worker!(
            ordinal,
            prepared,
            state,
            workspace,
            command,
            completion,
        )
    end
    task_ids = Memory{UInt}(undef, count)
    fill!(task_ids, UInt(0))
    thread_ids = Memory{Int}(undef, count)
    fill!(thread_ids, 0)
    return FixedOwnerOpticalPathBatchExecutor(
        prepared,
        state,
        workspace,
        commands,
        completions,
        tasks,
        Memory{Int}(undef, count),
        task_ids,
        thread_ids,
        0,
        0,
        0,
        0,
        0,
        false,
        0,
        0,
        0,
        false,
    )
end

function require_event_composition_executor_owners(
    executor::FixedOwnerOpticalPathBatchExecutor,
    prepared,
    state,
    workspace,
)
    prepared === executor.prepared ||
        error("fixed-owner executor received a foreign prepared plant")
    state === executor.state ||
        error("fixed-owner executor received a foreign state owner")
    workspace === executor.workspace ||
        error("fixed-owner executor received a foreign workspace owner")
    executor.closed &&
        error("fixed-owner executor cannot be used after close")
    return nothing
end

@inline function event_composition_dispatch_index(
    count::Int,
    index::Int,
    reverse_order::Bool,
)
    return reverse_order ? count - index + 1 : index
end

function dispatch_event_composition_group_phase!(
    executor::FixedOwnerOpticalPathBatchExecutor,
    claim::Plant.OpticalPathBatchClaim,
    due_count::Int,
    phase::EventCompositionWorkerPhase,
    reverse_order::Bool,
)
    @inbounds for index in 1:due_count
        due_index = event_composition_dispatch_index(
            due_count, index, reverse_order)
        ordinal = executor.due_group_ordinals[due_index]
        put!(executor.commands[ordinal],
            EventCompositionWorkerCommand(
                phase, claim, executor.measure_allocations))
    end
    @inbounds for index in 1:due_count
        ordinal = executor.due_group_ordinals[index]
        completion = take!(executor.completions[ordinal])
        completion.phase == phase ||
            error("fixed-owner worker completed the wrong phase")
        known_task_id = executor.task_ids[ordinal]
        if iszero(known_task_id)
            executor.task_ids[ordinal] = completion.task_id
        else
            completion.task_id == known_task_id ||
                error("path group changed task owner")
        end
        executor.thread_ids[ordinal] = completion.thread_id
        completion.exception === nothing ||
            throw(completion.exception)
        if executor.measure_allocations
            executor.measured_group_call_count += 1
            executor.measured_group_allocation_bytes +=
                completion.allocated_bytes
            executor.maximum_group_allocation_bytes = max(
                executor.maximum_group_allocation_bytes,
                completion.allocated_bytes,
            )
        end
    end
    if reverse_order
        executor.reverse_phase_count += 1
    else
        executor.forward_phase_count += 1
    end
    if phase == EventCompositionMaterializationPhase
        executor.materialization_count += due_count
    else
        executor.execution_count += due_count
    end
    return nothing
end

function Plant.execute_optical_path_batch!(
    executor::FixedOwnerOpticalPathBatchExecutor,
    prepared::Plant.PreparedPlantEventLoop,
    state::Plant.PlantEventLoopState,
    workspace::Plant.PlantEventLoopWorkspace,
    timestamp::PlantTimestamp,
)
    require_event_composition_executor_owners(
        executor, prepared, state, workspace)
    claim = Plant.begin_optical_path_batch!(
        prepared, state, workspace, timestamp)
    due_count = Plant.optical_path_batch_due_group_count(
        prepared, state, workspace, claim)
    @inbounds for index in 1:due_count
        executor.due_group_ordinals[index] =
            Plant.optical_path_batch_due_group_ordinal(
                prepared, state, workspace, claim, index)
    end

    reverse_materialization = iseven(executor.batch_count)
    dispatch_event_composition_group_phase!(
        executor,
        claim,
        due_count,
        EventCompositionMaterializationPhase,
        reverse_materialization,
    )
    Plant.seal_optical_path_batch_materialization!(
        prepared, state, workspace, claim)
    dispatch_event_composition_group_phase!(
        executor,
        claim,
        due_count,
        EventCompositionExecutionPhase,
        !reverse_materialization,
    )
    completed = Plant.complete_optical_path_batch!(
        prepared, state, workspace, claim)
    executor.batch_count += 1
    return completed
end

function close_event_composition_executor!(
    executor::FixedOwnerOpticalPathBatchExecutor,
)
    executor.closed && return executor
    @inbounds for command in executor.commands
        put!(command, nothing)
    end
    @inbounds for task in executor.tasks
        wait(task)
    end
    executor.closed = true
    return executor
end

@inline function run_event_composition_window!(prepared, state, workspace,
    stop::PlantTimestamp)
    return run_plant_events_until!(prepared, state, workspace, stop)
end

function event_composition_storage_signature(prepared, state, workspace)
    return (
        Base.summarysize(prepared.scheduler.definitions),
        Base.summarysize(prepared.actions),
        Base.summarysize(prepared.path_groups),
        Base.summarysize(prepared.acquisitions),
        Base.summarysize(state.scheduler.cursors),
        Base.summarysize(state.acquisitions),
        Base.summarysize(state.path_sampled),
        Base.summarysize(state.product_sequences),
        Base.summarysize(state.product_ready_timestamps),
        Base.summarysize(workspace.scheduler.due_slots),
        Base.summarysize(workspace.due_paths),
        Base.summarysize(workspace.optical_path_batch.due_group_slots),
        Base.summarysize(workspace.optical_path_batch.group_status),
    )
end

function optical_path_batch_allocation_signature!(
    prepared,
    state,
    workspace,
    timestamp::PlantTimestamp,
)
    local claim
    begin_bytes = @allocated claim = Plant.begin_optical_path_batch!(
        prepared, state, workspace, timestamp)
    count = Plant.optical_path_batch_due_group_count(
        prepared, state, workspace, claim)
    materialization_bytes = 0
    @inbounds for index in 1:count
        ordinal = Plant.optical_path_batch_due_group_ordinal(
            prepared, state, workspace, claim, index)
        materialization_bytes += @allocated(
            Plant.materialize_path_execution_group!(
                prepared, state, workspace, claim, ordinal))
    end
    seal_bytes = @allocated(
        Plant.seal_optical_path_batch_materialization!(
            prepared, state, workspace, claim))
    execution_bytes = 0
    @inbounds for index in 1:count
        ordinal = Plant.optical_path_batch_due_group_ordinal(
            prepared, state, workspace, claim, index)
        execution_bytes += @allocated(
            Plant.execute_path_execution_group!(
                prepared, state, workspace, claim, ordinal))
    end
    completion_bytes = @allocated Plant.complete_optical_path_batch!(
        prepared, state, workspace, claim)
    return (
        begin_phase=begin_bytes,
        materialization=materialization_bytes,
        seal=seal_bytes,
        execution=execution_bytes,
        completion=completion_bytes,
    )
end

@testset "Rolling-shutter row-band event lifecycle" begin
    values = reshape(collect(1.0:24.0), 6, 4)
    map = event_test_intensity_map(values)
    exposure = PlantDuration(1_000_000_000)
    line = PlantDuration(200_000_000)
    definition = RollingShutterAcquisitionDefinition(exposure;
        readiness_delay=PlantDuration(17))

    detector = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2)))
    mtf_before = detector_mtf(detector, 0.15, -0.2)
    prepared = prepare_rolling_shutter_acquisition(detector, map,
        definition)
    state = RollingShutterAcquisitionState(prepared)

    tiny_line_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(1.0e-12;
            row_group_size=2)))
    tiny_line_error = event_test_error() do
        prepare_rolling_shutter_acquisition(tiny_line_detector, map,
            definition)
    end
    @test tiny_line_error isa DetectorAcquisitionError
    @test tiny_line_error.reason == :unrepresentable_line_duration

    @test rolling_band_count(prepared) == 3
    @test rolling_band_rows(prepared, 1) == 1:2
    @test rolling_band_rows(prepared, 2) == 3:4
    @test rolling_band_rows(prepared, 3) == 5:6
    begin_exposure!(prepared, state, PlantTimestamp(100))
    @test rolling_opened_band_count(state) == 1
    @test rolling_closed_band_count(state) == 0
    @test next_rolling_band_open_timestamp(prepared, state) ==
        PlantTimestamp(100) + line
    @test next_rolling_band_close_timestamp(prepared, state) ==
        PlantTimestamp(100) + exposure
    @test readout_complete_timestamp(state) ==
        PlantTimestamp(100) + exposure + 3line

    first_open = next_rolling_band_open_timestamp(prepared, state)
    accumulate_rolling_exposure_interval!(prepared, state,
        PlantTimestamp(100), first_open, Xoshiro(700))
    open_next_rolling_band!(prepared, state, first_open)
    @test rolling_opened_band_count(state) == 2

    # Restart on a fresh owner for the complete constant-rate oracle.
    oracle_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2)))
    oracle_prepared = prepare_rolling_shutter_acquisition(oracle_detector,
        map, definition)
    oracle_state = RollingShutterAcquisitionState(oracle_prepared)
    output = run_rolling_event_frame!(oracle_prepared, oracle_state,
        Xoshiro(701), PlantTimestamp(0))
    @test output == values
    @test rolling_opened_band_count(oracle_state) == 3
    @test rolling_closed_band_count(oracle_state) == 3
    @test detector_acquisition_status(oracle_state) ==
        DetectorAcquisitionReady
    @test detector_mtf(oracle_detector, 0.15, -0.2) == mtf_before

    reset_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2, exposure_mode=GlobalResetExposure())))
    reset_prepared = prepare_rolling_shutter_acquisition(reset_detector,
        event_test_intensity_map(ones(6, 2)), definition)
    reset_state = RollingShutterAcquisitionState(reset_prepared)
    reset_output = run_rolling_event_frame!(reset_prepared, reset_state,
        Xoshiro(702), PlantTimestamp(0))
    @test reset_output == [
        1.0 1.0
        1.0 1.0
        1.2 1.2
        1.2 1.2
        1.4 1.4
        1.4 1.4
    ]
    @test rolling_band_open_timestamp(reset_prepared, reset_state, 3) ==
        PlantTimestamp(0)
    @test rolling_band_close_timestamp(reset_prepared, reset_state, 3) ==
        PlantTimestamp(1_400_000_000)

    if coverage_instrumented()
        @test_skip "rolling event allocation gate disabled under coverage instrumentation"
    else
        next_start = acquisition_readiness_timestamp(oracle_state) +
            PlantDuration(1)
        run_rolling_event_frame!(oracle_prepared, oracle_state,
            Xoshiro(703), next_start)
        following_start = acquisition_readiness_timestamp(oracle_state) +
            PlantDuration(1)
        allocation_rng = Xoshiro(704)
        @test @allocated(run_rolling_event_frame!(oracle_prepared,
            oracle_state, allocation_rng, following_start)) == 0
    end
end

@testset "Frame-transfer storage and overlapping acquisition" begin
    values = ones(3, 3)
    map = event_test_intensity_map(values)
    detector = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
        gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=0.1)))
    prepared = prepare_frame_transfer_acquisition(detector, map,
        FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(800_000_000)))
    state = FrameTransferAcquisitionState(prepared)
    rng = Xoshiro(710)
    mtf_before = detector_mtf(detector, 0.2, 0.1)

    tiny_transfer_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=1.0e-12)))
    tiny_transfer_error = event_test_error() do
        prepare_frame_transfer_acquisition(tiny_transfer_detector, map,
            FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000)))
    end
    @test tiny_transfer_error isa DetectorAcquisitionError
    @test tiny_transfer_error.reason == :unrepresentable_transfer_duration

    @test frame_transfer_storage_capacity(prepared) == 1
    @test !Base.mightalias(prepared.storage_frame, detector.state.frame)
    @test !Base.mightalias(prepared.storage_frame,
        detector.state.accum_buffer)
    @test !Base.mightalias(prepared.storage_frame, output_frame(detector))

    first_start = PlantTimestamp(0)
    first_close = PlantTimestamp(1_000_000_000)
    first_transfer = PlantTimestamp(1_100_000_000)
    first_readout = PlantTimestamp(1_900_000_000)
    begin_exposure!(prepared, state, first_start)
    accumulate_exposure_interval!(prepared, state, first_start,
        first_close, rng)
    close_exposure!(prepared, state, first_close)
    complete_frame_transfer!(prepared, state, first_transfer)
    @test frame_transfer_image_ready(state)
    @test frame_transfer_readout_pending(state)
    @test frame_transfer_storage_sequence(state) == 1

    regression = event_test_error() do
        begin_exposure!(prepared, state, first_transfer - PlantDuration(1))
    end
    @test regression isa DetectorAcquisitionError
    @test regression.reason == :time_regression

    slow_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=0.1)))
    slow_prepared = prepare_frame_transfer_acquisition(slow_detector, map,
        FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(1_500_000_000)))
    slow_state = FrameTransferAcquisitionState(slow_prepared)
    begin_exposure!(slow_prepared, slow_state, first_start)
    accumulate_exposure_interval!(slow_prepared, slow_state, first_start,
        first_close, Xoshiro(711))
    close_exposure!(slow_prepared, slow_state, first_close)
    complete_frame_transfer!(slow_prepared, slow_state, first_transfer)
    capacity_error = event_test_error() do
        begin_exposure!(slow_prepared, slow_state, first_transfer)
    end
    @test capacity_error isa DetectorAcquisitionError
    @test capacity_error.reason == :storage_capacity
    @test detector_acquisition_sequence(slow_state) == 1
    @test frame_transfer_storage_sequence(slow_state) == 1

    fill!(values, 2.0)
    second_start = first_transfer
    second_close = PlantTimestamp(2_100_000_000)
    second_transfer = PlantTimestamp(2_200_000_000)
    begin_exposure!(prepared, state, second_start)
    accumulate_exposure_interval!(prepared, state, second_start,
        first_readout, rng)
    first_output = copy(complete_readout!(prepared, state, first_readout,
        rng))
    @test first_output == ones(3, 3)
    @test frame_transfer_image_sequence(state) == 2
    @test frame_transfer_product_sequence(state) == 1

    accumulate_exposure_interval!(prepared, state, first_readout,
        second_close, rng)
    close_exposure!(prepared, state, second_close)
    complete_frame_transfer!(prepared, state, second_transfer)
    second_output = copy(complete_readout!(prepared, state,
        PlantTimestamp(3_000_000_000), rng))
    @test second_output == fill(2.0, 3, 3)
    @test frame_transfer_product_sequence(state) == 2
    @test frame_transfer_storage_empty(state)
    @test detector_mtf(detector, 0.2, 0.1) == mtf_before

    warm_start = PlantTimestamp(3_100_000_000)
    run_frame_transfer_event!(prepared, state, rng, warm_start)
    allocation_start = PlantTimestamp(5_100_000_000)
    if coverage_instrumented()
        @test_skip "frame-transfer allocation gate disabled under coverage instrumentation"
    else
        @test @allocated(run_frame_transfer_event!(prepared, state, rng,
            allocation_start)) == 0
    end
end

@testset "Multi-rate plant event composition" begin
    plant, prepared, state, workspace = event_composition_fixture()
    @test plant_event_path_count(prepared) == 3
    @test path_execution_group_count(prepared) == 3
    @test plant_event_acquisition_count(prepared) == 5
    for id in AcquisitionID.((
        :lgs_emccd,
        :ngs_saphira,
        :science_ccd,
        :science_cmos,
        :science_rolling,
    ))
        @test acquisition_products(prepared, id) ===
            acquisition_products(prepared_acquisition(plant, id))
    end
    unknown_products = event_test_error() do
        acquisition_products(prepared, :unknown_acquisition)
    end
    @test unknown_products isa PlantScheduleError
    @test unknown_products.reason == :unknown_acquisition
    @test plant_event_generator_count(prepared) == 1 + 3 + 5 * 5
    @test next_plant_event_timestamp(prepared, state, workspace) ==
        PlantTimestamp(0)

    @test prepared.path_groups isa Memory{PreparedPathExecutionGroup}
    expected_group_paths = OpticalPathID.((:lgs, :ngs, :science))
    expected_group_acquisitions = (
        AcquisitionID.((:lgs_emccd,)),
        AcquisitionID.((:ngs_saphira,)),
        AcquisitionID.(
            (:science_ccd, :science_cmos, :science_rolling)),
    )
    for ordinal in 1:path_execution_group_count(prepared)
        group = @inferred path_execution_group(prepared, ordinal)
        @test path_execution_group_ordinal(group) == ordinal
        @test path_execution_group_path_id(group) ==
            expected_group_paths[ordinal]
        requirements =
            @inferred path_execution_group_requirements(group)
        @test path_execution_backend(requirements) == CPUBackend()
        @test path_execution_compute_device(requirements) ==
            AdaptiveOpticsSim.Backends.HostComputeDevice()
        target_support =
            @inferred AbstractPathExecutionTargetSupport path_execution_target_support(
                group,
                AdaptiveOpticsSim.Backends.HostComputeDevice(),
            )
        @test typeof(target_support) === SupportedPathExecutionTarget
        @test path_execution_target_supported(target_support)
        @test isnothing(
            path_execution_target_rejection_reason(target_support))
        @test path_execution_requires_full_optical(requirements)
        @test path_execution_group_acquisition_count(group) ==
            length(expected_group_acquisitions[ordinal])
        @test Tuple(path_execution_group_acquisition_id(group, index)
            for index in
                1:path_execution_group_acquisition_count(group)) ==
            expected_group_acquisitions[ordinal]
        @test @inferred(path_execution_group_acquisition_id(group, 1)) ==
            first(expected_group_acquisitions[ordinal])
    end
    first_group = path_execution_group(prepared, 1)
    first_requirements = path_execution_group_requirements(first_group)
    path_execution_backend(first_requirements)
    path_execution_compute_device(first_requirements)
    path_execution_target_support(
        first_requirements,
        AdaptiveOpticsSim.Backends.HostComputeDevice(),
    )
    path_execution_requires_full_optical(first_requirements)
    if coverage_instrumented()
        @test_skip "path-requirements allocation gate disabled under coverage instrumentation"
    else
        @test @allocated(path_execution_group_requirements(
            first_group)) == 0
        @test @allocated(path_execution_backend(
            first_requirements)) == 0
        @test @allocated(path_execution_compute_device(
            first_requirements)) == 0
        @test @allocated(path_execution_target_support(
            first_requirements,
            AdaptiveOpticsSim.Backends.HostComputeDevice())) == 0
        @test @allocated(path_execution_requires_full_optical(
            first_requirements)) == 0
    end
    @test_throws BoundsError path_execution_group(prepared, 0)
    @test_throws BoundsError path_execution_group(prepared, 4)
    @test_throws BoundsError path_execution_group_acquisition_id(
        path_execution_group(prepared, 1), 2)
    @test_throws PlantPreparationError PathExecutionRequirements(
        CUDABackend(),
        AdaptiveOpticsSim.Backends.HostComputeDevice(),
        true,
    )
    cuda_device_0 = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 0)
    cuda_device_1 = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 1)
    backend_mismatch =
        @inferred AbstractPathExecutionTargetSupport path_execution_target_support(
            first_requirements, cuda_device_0)
    @test typeof(backend_mismatch) === UnsupportedPathExecutionTarget
    @test !path_execution_target_supported(backend_mismatch)
    @test path_execution_target_rejection_reason(backend_mismatch) ==
        :backend_mismatch
    cuda_requirements = PathExecutionRequirements(
        CUDABackend(), cuda_device_0, true)
    requires_repreparation =
        @inferred AbstractPathExecutionTargetSupport path_execution_target_support(
            cuda_requirements, cuda_device_1)
    @test typeof(requires_repreparation) === UnsupportedPathExecutionTarget
    @test !path_execution_target_supported(requires_repreparation)
    @test path_execution_target_rejection_reason(requires_repreparation) ==
        :requires_repreparation
    for name in (
        :PreparedPathExecutionGroup,
        :PathExecutionRequirements,
        :AbstractPathExecutionTargetSupport,
        :SupportedPathExecutionTarget,
        :UnsupportedPathExecutionTarget,
        :path_execution_group_count,
        :path_execution_group,
        :path_execution_group_ordinal,
        :path_execution_group_path_id,
        :path_execution_group_requirements,
        :path_execution_backend,
        :path_execution_compute_device,
        :path_execution_target_support,
        :path_execution_target_supported,
        :path_execution_target_rejection_reason,
        :path_execution_requires_full_optical,
        :path_execution_group_acquisition_count,
        :path_execution_group_acquisition_id,
        :OpticalPathBatchClaim,
        :AbstractOpticalPathBatchExecutor,
        :SerialOpticalPathBatchExecutor,
        :optical_path_batch_timestamp,
        :optical_path_batch_epoch,
        :optical_path_batch_due_group_count,
        :optical_path_batch_due_group_ordinal,
        :begin_optical_path_batch!,
        :materialize_path_execution_group!,
        :seal_optical_path_batch_materialization!,
        :execute_path_execution_group!,
        :complete_optical_path_batch!,
        :abandon_optical_path_batch!,
        :execute_optical_path_batch!,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end

    @testset "Optical path batch lifecycle" begin
        batch_plant, batch_prepared, batch_state, batch_workspace =
            event_composition_fixture(aligned_optical_samples=true)
        claim = @inferred Plant.begin_optical_path_batch!(
            batch_prepared,
            batch_state,
            batch_workspace,
            PlantTimestamp(0),
        )
        @test claim isa Plant.OpticalPathBatchClaim
        @test isbitstype(typeof(claim))
        @test Plant.optical_path_batch_timestamp(claim) ==
            PlantTimestamp(0)
        @test Plant.optical_path_batch_epoch(
            batch_prepared, batch_state, batch_workspace, claim) ==
            current_epoch(batch_prepared.atmosphere)
        @test @inferred(Plant.optical_path_batch_due_group_count(
            batch_prepared, batch_state, batch_workspace, claim)) == 3
        @test Tuple(Plant.optical_path_batch_due_group_ordinal(
            batch_prepared, batch_state, batch_workspace, claim, index)
            for index in 1:3) == (1, 2, 3)

        active_error = event_test_error() do
            Plant.begin_optical_path_batch!(
                batch_prepared,
                batch_state,
                batch_workspace,
                PlantTimestamp(0),
            )
        end
        @test active_error isa PlantScheduleError
        @test active_error.reason == :optical_path_batch_active
        active_timestamp_error = event_test_error() do
            next_plant_event_timestamp(
                batch_prepared, batch_state, batch_workspace)
        end
        @test active_timestamp_error isa PlantScheduleError
        @test active_timestamp_error.reason == :optical_path_batch_active

        foreign_state = PlantEventLoopState(batch_prepared)
        foreign_state_error = event_test_error() do
            Plant.optical_path_batch_due_group_count(
                batch_prepared, foreign_state, batch_workspace, claim)
        end
        @test foreign_state_error isa PlantScheduleError
        @test foreign_state_error.reason ==
            :foreign_optical_path_batch_state
        foreign_workspace = PlantEventLoopWorkspace(batch_prepared)
        foreign_workspace_error = event_test_error() do
            Plant.optical_path_batch_due_group_count(
                batch_prepared, batch_state, foreign_workspace, claim)
        end
        @test foreign_workspace_error isa PlantScheduleError
        @test foreign_workspace_error.reason ==
            :foreign_optical_path_batch_workspace

        invalid_due_error = event_test_error() do
            Plant.optical_path_batch_due_group_ordinal(
                batch_prepared, batch_state, batch_workspace, claim, 4)
        end
        @test invalid_due_error isa PlantScheduleError
        @test invalid_due_error.reason ==
            :invalid_due_path_execution_group
        invalid_group_error = event_test_error() do
            Plant.materialize_path_execution_group!(
                batch_prepared, batch_state, batch_workspace, claim, 4)
        end
        @test invalid_group_error isa PlantScheduleError
        @test invalid_group_error.reason == :invalid_path_execution_group

        early_execute_error = event_test_error() do
            Plant.execute_path_execution_group!(
                batch_prepared, batch_state, batch_workspace, claim, 1)
        end
        @test early_execute_error isa PlantScheduleError
        @test early_execute_error.reason ==
            :invalid_optical_path_batch_phase
        early_seal_error = event_test_error() do
            Plant.seal_optical_path_batch_materialization!(
                batch_prepared, batch_state, batch_workspace, claim)
        end
        @test early_seal_error isa PlantScheduleError
        @test early_seal_error.reason ==
            :incomplete_optical_path_batch_materialization

        for ordinal in 3:-1:1
            Plant.materialize_path_execution_group!(
                batch_prepared, batch_state, batch_workspace, claim,
                ordinal)
            if ordinal == 3
                duplicate_materialization_error = event_test_error() do
                    Plant.materialize_path_execution_group!(
                        batch_prepared, batch_state, batch_workspace, claim,
                        ordinal)
                end
                @test duplicate_materialization_error isa PlantScheduleError
                @test duplicate_materialization_error.reason ==
                    :duplicate_path_execution_group_materialization
            end
        end
        Plant.seal_optical_path_batch_materialization!(
            batch_prepared, batch_state, batch_workspace, claim)
        incomplete_error = event_test_error() do
            Plant.complete_optical_path_batch!(
                batch_prepared, batch_state, batch_workspace, claim)
        end
        @test incomplete_error isa PlantScheduleError
        @test incomplete_error.reason ==
            :incomplete_optical_path_batch_execution

        sealed_epoch = current_epoch(batch_prepared.atmosphere)
        advanced_epoch = advance_to!(
            batch_prepared.atmosphere,
            0.001;
            rng=Xoshiro(0x7a0f),
        )
        @test epoch_sequence(advanced_epoch) >
            epoch_sequence(sealed_epoch)
        @test Plant.optical_path_batch_epoch(
            batch_prepared, batch_state, batch_workspace, claim) ==
            sealed_epoch

        for ordinal in 3:-1:1
            Plant.execute_path_execution_group!(
                batch_prepared, batch_state, batch_workspace, claim,
                ordinal)
            if ordinal == 3
                duplicate_execution_error = event_test_error() do
                    Plant.execute_path_execution_group!(
                        batch_prepared, batch_state, batch_workspace, claim,
                        ordinal)
                end
                @test duplicate_execution_error isa PlantScheduleError
                @test duplicate_execution_error.reason ==
                    :duplicate_path_execution_group
            end
        end
        @test Plant.complete_optical_path_batch!(
            batch_prepared, batch_state, batch_workspace, claim) ==
            PlantTimestamp(0)
        @test prepared_path(batch_plant, :science).execution.executions[] ==
            1
        @test prepared_path(batch_plant, :ngs).execution.executions[] == 1
        @test prepared_path(batch_plant, :lgs).execution.executions[] == 1
        stale_error = event_test_error() do
            Plant.optical_path_batch_due_group_count(
                batch_prepared, batch_state, batch_workspace, claim)
        end
        @test stale_error isa PlantScheduleError
        @test stale_error.reason == :stale_optical_path_batch_claim

        _, abandoned_prepared, abandoned_state, abandoned_workspace =
            event_composition_fixture(aligned_optical_samples=true)
        abandoned_claim = Plant.begin_optical_path_batch!(
            abandoned_prepared,
            abandoned_state,
            abandoned_workspace,
            PlantTimestamp(0),
        )
        @test @inferred(Plant.abandon_optical_path_batch!(
            abandoned_prepared,
            abandoned_state,
            abandoned_workspace,
            abandoned_claim,
        )) == PlantTimestamp(0)
        abandoned_claim_error = event_test_error() do
            Plant.optical_path_batch_due_group_count(
                abandoned_prepared,
                abandoned_state,
                abandoned_workspace,
                abandoned_claim,
            )
        end
        @test abandoned_claim_error isa PlantScheduleError
        @test abandoned_claim_error.reason ==
            :stale_optical_path_batch_claim
        abandoned_step_error = event_test_error() do
            step_plant_events!(
                abandoned_prepared,
                abandoned_state,
                abandoned_workspace,
            )
        end
        @test abandoned_step_error isa PlantScheduleError
        @test abandoned_step_error.reason == :optical_path_batch_active
        abandoned_timestamp_error = event_test_error() do
            next_plant_event_timestamp(
                abandoned_prepared,
                abandoned_state,
                abandoned_workspace,
            )
        end
        @test abandoned_timestamp_error isa PlantScheduleError
        @test abandoned_timestamp_error.reason ==
            :optical_path_batch_active
        repeated_abandonment_error = event_test_error() do
            Plant.abandon_optical_path_batch!(
                abandoned_prepared,
                abandoned_state,
                abandoned_workspace,
                abandoned_claim,
            )
        end
        @test repeated_abandonment_error isa PlantScheduleError
        @test repeated_abandonment_error.reason ==
            :stale_optical_path_batch_claim

        serial_failure = ErrorException(
            "intentional serial optical-path failure")
        _, serial_failure_prepared, serial_failure_state,
            serial_failure_workspace = event_composition_fixture(
                aligned_optical_samples=true,
                execution_probe=EventCompositionFailureProbe(
                    serial_failure),
            )
        observed_serial_failure = event_test_error() do
            step_plant_events!(
                serial_failure_prepared,
                serial_failure_state,
                serial_failure_workspace,
            )
        end
        @test observed_serial_failure === serial_failure
        @test serial_failure_workspace.optical_path_batch.phase ==
            Plant._OpticalPathBatchAbandoned
        serial_failure_step_error = event_test_error() do
            step_plant_events!(
                serial_failure_prepared,
                serial_failure_state,
                serial_failure_workspace,
            )
        end
        @test serial_failure_step_error isa PlantScheduleError
        @test serial_failure_step_error.reason ==
            :optical_path_batch_active

        _, stale_prepared, stale_state, stale_workspace =
            event_composition_fixture(aligned_optical_samples=true)
        stale_epoch_claim = Plant.begin_optical_path_batch!(
            stale_prepared,
            stale_state,
            stale_workspace,
            PlantTimestamp(0),
        )
        foreign_claim_error = event_test_error() do
            Plant.optical_path_batch_due_group_count(
                stale_prepared,
                stale_state,
                stale_workspace,
                claim,
            )
        end
        @test foreign_claim_error isa PlantScheduleError
        @test foreign_claim_error.reason ==
            :foreign_optical_path_batch_claim
        advance_to!(
            stale_prepared.atmosphere,
            0.001;
            rng=Xoshiro(0x7a10),
        )
        stale_epoch_error = event_test_error() do
            Plant.materialize_path_execution_group!(
                stale_prepared,
                stale_state,
                stale_workspace,
                stale_epoch_claim,
                1,
            )
        end
        @test stale_epoch_error isa AtmosphereEpochError
        failed_materialization_retry = event_test_error() do
            Plant.materialize_path_execution_group!(
                stale_prepared,
                stale_state,
                stale_workspace,
                stale_epoch_claim,
                1,
            )
        end
        @test failed_materialization_retry isa PlantScheduleError
        @test failed_materialization_retry.reason ==
            :path_execution_group_materialization_active

        serial_plant, serial_prepared, serial_state, serial_workspace =
            event_composition_fixture(aligned_optical_samples=true)
        reverse_plant, reverse_prepared, reverse_state, reverse_workspace =
            event_composition_fixture(aligned_optical_samples=true)
        serial_timestamp = step_plant_events!(
            serial_prepared,
            serial_state,
            serial_workspace,
        )
        reverse_timestamp = step_plant_events!(
            reverse_prepared,
            reverse_state,
            reverse_workspace,
            ReverseOpticalPathBatchExecutor(),
        )
        @test reverse_timestamp == serial_timestamp == PlantTimestamp(0)
        @test epoch_time(current_epoch(reverse_prepared.atmosphere)) ==
            epoch_time(current_epoch(serial_prepared.atmosphere))
        @test epoch_sequence(current_epoch(reverse_prepared.atmosphere)) ==
            epoch_sequence(current_epoch(serial_prepared.atmosphere))
        for id in (:science, :ngs, :lgs)
            serial_path = prepared_path(serial_plant, id)
            reverse_path = prepared_path(reverse_plant, id)
            @test reverse_path.result.values == serial_path.result.values
            @test reverse_path.input.opd == serial_path.input.opd
            @test reverse_path.execution.executions[] ==
                serial_path.execution.executions[] == 1
            serial_group = path_execution_group(
                serial_prepared,
                findfirst(group ->
                    path_execution_group_path_id(group) ==
                        OpticalPathID(id),
                    serial_prepared.path_groups),
            )
            reverse_group = path_execution_group(
                reverse_prepared,
                findfirst(group ->
                    path_execution_group_path_id(group) ==
                        OpticalPathID(id),
                    reverse_prepared.path_groups),
            )
            @test copy(Plant.rng_stream_state(
                reverse_group.rngs, Val(:provider))) ==
                copy(Plant.rng_stream_state(
                    serial_group.rngs, Val(:provider)))
        end

        @testset "Fixed-owner grouped CPU execution" begin
            if Threads.nthreads() < 4
                @test_skip "grouped CPU overlap requires one coordinator and three worker threads"
            elseif BLAS.get_num_threads() != 1
                @test_skip "grouped CPU overlap requires the declared one-thread BLAS policy"
            else
                budget = Plant.grouped_cpu_execution_budget(
                    cpu_context_count=Threads.nthreads(),
                    julia_thread_count=Threads.nthreads(),
                    outer_owner_count=3,
                    group_julia_thread_count=1,
                    fft_thread_count=1,
                    blas_thread_count=1,
                )
                environment = Plant.CPUExecutionEnvironment(
                    available_cpu_context_count=Threads.nthreads(),
                    fft_thread_count=1,
                )
                @test Plant.validate_cpu_execution_budget(
                    budget, environment) === budget

                serial_grouped_plant,
                    serial_grouped_prepared,
                    serial_grouped_state,
                    serial_grouped_workspace =
                    event_composition_fixture(
                        aligned_optical_samples=true)
                concurrency_probe = EventCompositionConcurrencyProbe(3)
                grouped_plant,
                    grouped_prepared,
                    grouped_state,
                    grouped_workspace =
                    event_composition_fixture(
                        aligned_optical_samples=true,
                        execution_probe=concurrency_probe)
                executor = FixedOwnerOpticalPathBatchExecutor(
                    grouped_prepared,
                    grouped_state,
                    grouped_workspace,
                )
                horizon = PlantTimestamp(1_500_000_000)
                warmup_horizon = PlantTimestamp(500_000_000)
                local serial_timestamp_count
                local grouped_timestamp_count
                try
                    serial_warmup_count = run_plant_events_until!(
                        serial_grouped_prepared,
                        serial_grouped_state,
                        serial_grouped_workspace,
                        warmup_horizon,
                    )
                    grouped_warmup_count = run_plant_events_until!(
                        grouped_prepared,
                        grouped_state,
                        grouped_workspace,
                        warmup_horizon,
                        executor,
                    )
                    @test grouped_warmup_count == serial_warmup_count
                    executor.measure_allocations = true
                    serial_timestamp_count = run_plant_events_until!(
                        serial_grouped_prepared,
                        serial_grouped_state,
                        serial_grouped_workspace,
                        horizon,
                    )
                    grouped_timestamp_count = run_plant_events_until!(
                        grouped_prepared,
                        grouped_state,
                        grouped_workspace,
                        horizon,
                        executor,
                    )

                    @test grouped_timestamp_count ==
                        serial_timestamp_count
                    @test scheduler_timestamp(
                        grouped_state.scheduler) ==
                        scheduler_timestamp(
                            serial_grouped_state.scheduler)
                    @test grouped_state.scheduler.revision ==
                        serial_grouped_state.scheduler.revision
                    @test grouped_state.scheduler.cursors ==
                        serial_grouped_state.scheduler.cursors
                    @test grouped_state.path_sampled ==
                        serial_grouped_state.path_sampled
                    @test grouped_state.product_sequences ==
                        serial_grouped_state.product_sequences
                    @test grouped_state.product_ready_timestamps ==
                        serial_grouped_state.product_ready_timestamps

                    grouped_epoch =
                        current_epoch(grouped_prepared.atmosphere)
                    serial_epoch =
                        current_epoch(serial_grouped_prepared.atmosphere)
                    @test epoch_time(grouped_epoch) ==
                        epoch_time(serial_epoch)
                    @test epoch_sequence(grouped_epoch) ==
                        epoch_sequence(serial_epoch)
                    for index in eachindex(
                        grouped_prepared.atmosphere_rng.streams)
                        @test copy(grouped_prepared.atmosphere_rng.
                            streams[index].state) ==
                            copy(serial_grouped_prepared.atmosphere_rng.
                                streams[index].state)
                    end

                    for id in (:science, :ngs, :lgs)
                        serial_path =
                            prepared_path(serial_grouped_plant, id)
                        grouped_path = prepared_path(grouped_plant, id)
                        @test grouped_path.input.opd ==
                            serial_path.input.opd
                        @test grouped_path.result.values ==
                            serial_path.result.values
                        @test grouped_path.execution.executions[] ==
                            serial_path.execution.executions[]
                    end

                    for ordinal in
                            1:path_execution_group_count(grouped_prepared)
                        grouped_group = path_execution_group(
                            grouped_prepared, ordinal)
                        serial_group = path_execution_group(
                            serial_grouped_prepared, ordinal)
                        @test path_execution_group_path_id(
                            grouped_group) ==
                            path_execution_group_path_id(serial_group)
                        @test copy(Plant.rng_stream_state(
                            grouped_group.rngs, Val(:provider))) ==
                            copy(Plant.rng_stream_state(
                                serial_group.rngs, Val(:provider)))
                    end

                    acquisition_ids = (
                        :science_cmos,
                        :science_ccd,
                        :science_rolling,
                        :ngs_saphira,
                        :lgs_emccd,
                    )
                    for id in acquisition_ids
                        @test acquisition_product_sequence(
                            grouped_prepared, grouped_state, id) ==
                            acquisition_product_sequence(
                                serial_grouped_prepared,
                                serial_grouped_state,
                                id,
                            )
                        @test acquisition_product_ready_timestamp(
                            grouped_prepared, grouped_state, id) ==
                            acquisition_product_ready_timestamp(
                                serial_grouped_prepared,
                                serial_grouped_state,
                                id,
                            )
                        @test acquisition_observation(
                            prepared_acquisition(grouped_plant, id)) ==
                            acquisition_observation(prepared_acquisition(
                                serial_grouped_plant, id))
                    end
                    for index in eachindex(grouped_prepared.acquisitions)
                        @test copy(
                            grouped_prepared.acquisitions[index].rng) ==
                            copy(serial_grouped_prepared.
                                acquisitions[index].rng)
                    end

                    @test executor.batch_count > 1
                    @test executor.forward_phase_count > 0
                    @test executor.reverse_phase_count > 0
                    @test executor.materialization_count ==
                        executor.execution_count
                    @test executor.execution_count ==
                        sum(prepared_path(grouped_plant, id).
                            execution.executions[] for
                            id in (:science, :ngs, :lgs))
                    @test executor.measured_group_call_count > 0
                    if coverage_instrumented()
                        @test_skip "grouped worker allocation gate disabled under coverage instrumentation"
                    else
                        @test executor.maximum_group_allocation_bytes <= 512
                        @test executor.measured_group_allocation_bytes <=
                            512 * executor.measured_group_call_count
                    end
                    @test length(executor.tasks) == 3
                    @test all(task -> !istaskdone(task), executor.tasks)
                    @test all(id -> !iszero(id), executor.task_ids)
                    @test length(unique(collect(executor.task_ids))) == 3
                    @test concurrency_probe.arrivals[] == 3
                    @test all(
                        id -> !iszero(id), concurrency_probe.thread_ids)
                    @test length(unique(collect(
                        concurrency_probe.thread_ids))) == 3
                finally
                    close_event_composition_executor!(executor)
                end
                @test executor.closed
                @test all(istaskdone, executor.tasks)
            end
        end

        if coverage_instrumented()
            @test_skip "optical-path batch allocation gate disabled under coverage instrumentation"
        else
            _, allocation_prepared, allocation_state,
                allocation_workspace = event_composition_fixture(
                    aligned_optical_samples=true)
            signature = optical_path_batch_allocation_signature!(
                allocation_prepared,
                allocation_state,
                allocation_workspace,
                PlantTimestamp(0),
            )
            # Julia 1.12 warmed contract for three full-optical groups. Cold
            # preparation, compilation, and exception formatting are outside
            # the window; current heterogeneous model barriers remain inside.
            # The batch phase transitions and reclamation themselves are
            # allocation-free, and the complete bounded orchestration stays
            # within the existing 2 KiB per-timestamp ceiling.
            @test signature.begin_phase <= 1_024
            @test signature.materialization <= 512
            @test signature.seal == 0
            @test signature.execution <= 512
            @test signature.completion == 0
            @test sum(signature) <= 2_048

            _, abandonment_prepared, abandonment_state,
                abandonment_workspace = event_composition_fixture(
                    aligned_optical_samples=true)
            abandonment_claim = Plant.begin_optical_path_batch!(
                abandonment_prepared,
                abandonment_state,
                abandonment_workspace,
                PlantTimestamp(0),
            )
            @test @allocated(Plant.abandon_optical_path_batch!(
                abandonment_prepared,
                abandonment_state,
                abandonment_workspace,
                abandonment_claim,
            )) == 0
        end
    end

    missing_path_error = event_test_error() do
        event_composition_fixture(missing_path_schedule=true)
    end
    @test missing_path_error isa PlantScheduleError
    @test missing_path_error.reason == :missing_path_schedule

    unbound_error = event_test_error() do
        event_composition_fixture(unbound_trigger_consumer=true)
    end
    @test unbound_error isa PlantScheduleError
    @test unbound_error.component == :plant_event_loop
    @test unbound_error.reason == :unbound_trigger_consumer

    _, atomic_prepared, atomic_state, atomic_workspace =
        event_composition_fixture()
    run_plant_events_until!(atomic_prepared, atomic_state,
        atomic_workspace, PlantTimestamp(299_999_999))
    ccd_slot = findfirst(acquisition ->
        acquisition.id == AcquisitionID(:science_ccd),
        atomic_prepared.acquisitions)
    rolling_slot = findfirst(acquisition ->
        acquisition.id == AcquisitionID(:science_rolling),
        atomic_prepared.acquisitions)
    ccd_state = atomic_state.acquisitions[ccd_slot]
    ccd_progress = integrated_through_timestamp(ccd_state)
    rolling_lifecycle = atomic_prepared.acquisitions[
        rolling_slot].lifecycle
    rolling_lifecycle.detector.state.readout_ready = true
    preflight_error = event_test_error() do
        step_plant_events!(atomic_prepared, atomic_state, atomic_workspace)
    end
    @test preflight_error isa DetectorAcquisitionError
    @test preflight_error.reason == :detector_state_changed
    @test integrated_through_timestamp(ccd_state) == ccd_progress

    horizon = PlantTimestamp(1_500_000_000)
    timestamp_count = run_plant_events_until!(prepared, state, workspace,
        horizon)
    @test timestamp_count > 0
    @test acquisition_product_sequence(prepared, state,
        :science_cmos) == 3
    @test acquisition_product_sequence(prepared, state,
        :science_ccd) == 2
    @test acquisition_product_sequence(prepared, state,
        :science_rolling) == 3
    @test acquisition_product_sequence(prepared, state,
        :ngs_saphira) == 4
    @test acquisition_product_sequence(prepared, state,
        :lgs_emccd) == 3
    @test acquisition_product_ready_timestamp(prepared, state,
        :science_cmos) == PlantTimestamp(1_220_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :science_ccd) == PlantTimestamp(950_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :ngs_saphira) == PlantTimestamp(1_460_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :lgs_emccd) == PlantTimestamp(1_245_000_000)

    @test prepared_path(plant, :science).execution.executions[] == 16
    @test prepared_path(plant, :ngs).execution.executions[] == 12
    @test prepared_path(plant, :lgs).execution.executions[] == 8
    science_cmos = prepared_acquisition(plant, :science_cmos)
    science_ccd = prepared_acquisition(plant, :science_ccd)
    @test science_cmos.path_result === science_ccd.path_result
    @test all(isfinite, acquisition_observation(science_cmos))
    @test all(isfinite, acquisition_observation(science_ccd))

    reordered_plant, reordered, reordered_state, reordered_workspace =
        event_composition_fixture(reverse_order=true)
    @test path_execution_group_count(reordered) ==
        path_execution_group_count(prepared)
    for ordinal in 1:path_execution_group_count(prepared)
        canonical_group = path_execution_group(prepared, ordinal)
        reordered_group = path_execution_group(reordered, ordinal)
        @test path_execution_group_ordinal(reordered_group) ==
            path_execution_group_ordinal(canonical_group)
        @test path_execution_group_path_id(reordered_group) ==
            path_execution_group_path_id(canonical_group)
        @test Tuple(path_execution_group_acquisition_id(reordered_group,
            index) for index in
                1:path_execution_group_acquisition_count(
                    reordered_group)) ==
            Tuple(path_execution_group_acquisition_id(canonical_group,
                index) for index in
                    1:path_execution_group_acquisition_count(
                        canonical_group))
    end
    @test run_plant_events_until!(reordered, reordered_state,
        reordered_workspace, horizon) == timestamp_count
    for id in (:science_cmos, :science_ccd, :science_rolling,
            :ngs_saphira, :lgs_emccd)
        @test acquisition_product_sequence(reordered, reordered_state, id) ==
            acquisition_product_sequence(prepared, state, id)
        @test acquisition_product_ready_timestamp(reordered,
            reordered_state, id) == acquisition_product_ready_timestamp(
            prepared, state, id)
        @test acquisition_observation(prepared_acquisition(reordered_plant,
            id)) == acquisition_observation(prepared_acquisition(plant, id))
    end


    storage_before = event_composition_storage_signature(prepared, state,
        workspace)
    run_event_composition_window!(prepared, state, workspace,
        PlantTimestamp(2_000_000_000))
    if coverage_instrumented()
        @test_skip "plant event-loop allocation gate disabled under coverage instrumentation"
    else
        processed = Ref(0)
        allocated = @allocated processed[] =
            run_event_composition_window!(prepared, state, workspace,
                PlantTimestamp(2_500_000_000))
        # The heterogeneous orchestration barrier deliberately avoids
        # specializing on the complete instrument tuple. Keep its bounded
        # dynamic-dispatch cost below 2 KiB per processed timestamp while the
        # detector and optical kernels remain allocation-free after warmup.
        @test allocated <= 2_048 * processed[]
    end
    long_run_timestamps = run_event_composition_window!(prepared, state,
        workspace, PlantTimestamp(25_000_000_000))
    @test 0 < long_run_timestamps < 2_000
    @test event_composition_storage_signature(prepared, state, workspace) ==
        storage_before
    @test plant_event_generator_count(prepared) == 29

    faulted_plant, faulted, faulted_state, faulted_workspace =
        event_composition_fixture(faulted_trigger_fanout=true)
    faulted_timestamp_count = run_plant_events_until!(faulted,
        faulted_state, faulted_workspace, PlantTimestamp(3_000_000_000))
    @test faulted_timestamp_count > 0
    @test acquisition_product_sequence(faulted, faulted_state,
        :science_cmos) == 5
    @test acquisition_product_sequence(faulted, faulted_state,
        :ngs_saphira) == 5
    @test acquisition_product_ready_timestamp(faulted, faulted_state,
        :science_cmos) == PlantTimestamp(2_690_000_000)
    @test acquisition_product_ready_timestamp(faulted, faulted_state,
        :ngs_saphira) == PlantTimestamp(2_680_000_000)
    @test all(isfinite, acquisition_observation(prepared_acquisition(
        faulted_plant, :science_cmos)))
    @test all(isfinite, acquisition_observation(prepared_acquisition(
        faulted_plant, :ngs_saphira)))
end
