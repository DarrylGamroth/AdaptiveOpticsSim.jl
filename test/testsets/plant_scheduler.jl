function assert_event_scheduler_error(f, reason::Symbol)
    try
        f()
        @test false
    catch error
        @test error isa PlantScheduleError
        if error isa PlantScheduleError
            @test error.component == :event_scheduler
            @test error.reason == reason
            @test !isempty(error.msg)
        end
    end
end

function collect_scheduler_trace(definitions)
    scheduler = prepare_event_scheduler(definitions;
        capacity=length(definitions))
    state = EventSchedulerState(scheduler)
    workspace = EventSchedulerWorkspace(scheduler)
    trace = PlantEventKey[]
    while true
        claim = claim_next_event!(workspace, scheduler, state)
        claim === nothing && break
        push!(trace, claimed_event_key(claim))
        deactivate_event_generator!(scheduler, state, claim)
    end
    return trace
end

function scheduler_hot_cycle!(scheduler, state, workspace, schedule)
    claim = claim_next_event!(workspace, scheduler, state)
    claim === nothing && return zero(PlantTimestamp)
    key = reschedule_periodic_event!(scheduler, state, claim, schedule)
    return key.timestamp
end

@testset "Fixed-capacity event scheduler preparation" begin
    periodic = PeriodicSchedule(period_ns=500_000, phase_ns=125_000)
    periodic_definition = EventGeneratorDefinition(
        periodic, OpticalSamplePhase, 7;
        origin=PlantTimestamp(1_000_000),
    )
    inactive_definition = EventGeneratorDefinition(
        PlantTimestamp(2_000_000), AcquisitionReadyPhase, 3;
        occurrence=4, active=false,
    )

    @test isbitstype(EventGeneratorDefinition)
    @test isbitstype(EventGeneratorCursor)
    @test isbitstype(EventGeneratorHandle)
    @test isbitstype(EventClaim)
    @test first_event_timestamp(periodic_definition) ==
        PlantTimestamp(1_125_000)
    @test first_event_occurrence(periodic_definition) == 1
    @test event_generator_phase(periodic_definition) == OpticalSamplePhase
    @test event_generator_ordinal(periodic_definition) == 7
    @test initially_active(periodic_definition)
    @test !initially_active(inactive_definition)
    @test first_event_occurrence(inactive_definition) == 4

    scheduler = prepare_event_scheduler(
        (inactive_definition, periodic_definition); capacity=4)
    @test event_generator_count(scheduler) == 2
    @test event_scheduler_capacity(scheduler) == 4
    @test getfield(scheduler, :definitions) isa
        Memory{EventGeneratorDefinition}
    @test length(EventSchedulerState(scheduler).cursors) == 2
    @test length(EventSchedulerWorkspace(scheduler).due_slots) == 4
    @test event_generator_handle(scheduler, OpticalSamplePhase, 7) isa
        EventGeneratorHandle

    vector_definitions = [periodic_definition, inactive_definition]
    vector_scheduler = prepare_event_scheduler(vector_definitions;
        capacity=4)
    vector_definitions[1] = inactive_definition
    @test event_generator_count(vector_scheduler) == 2
    @test event_generator_phase(getfield(
        vector_scheduler, :definitions)[1]) == OpticalSamplePhase

    same_ordinal_other_phase = EventGeneratorDefinition(
        PlantTimestamp(1), ExposureOpenPhase, 7)
    @test event_generator_count(prepare_event_scheduler(
        (periodic_definition, same_ordinal_other_phase))) == 2

    duplicate = EventGeneratorDefinition(
        PlantTimestamp(9), OpticalSamplePhase, 7)
    assert_event_scheduler_error(
        () -> prepare_event_scheduler((periodic_definition, duplicate)),
        :duplicate_ordinal,
    )
    assert_event_scheduler_error(
        () -> prepare_event_scheduler(
            (periodic_definition, inactive_definition); capacity=1),
        :capacity_overflow,
    )
    assert_event_scheduler_error(
        () -> prepare_event_scheduler((periodic_definition,); capacity=-1),
        :invalid_capacity,
    )
    assert_event_scheduler_error(
        () -> prepare_event_scheduler((periodic_definition,); capacity=true),
        :invalid_capacity,
    )
    assert_event_scheduler_error(
        () -> prepare_event_scheduler((periodic_definition, :invalid)),
        :invalid_definition,
    )
    assert_event_scheduler_error(
        () -> prepare_event_scheduler((definition=periodic_definition,)),
        :invalid_registry,
    )
    assert_event_scheduler_error(
        () -> event_generator_handle(scheduler, OpticalSamplePhase, 99),
        :unknown_generator,
    )
    assert_event_scheduler_error(
        () -> EventSchedulerState(scheduler;
            initial_timestamp=PlantTimestamp(1_125_001)),
        :time_regression,
    )

    for name in (
        :EventGeneratorDefinition,
        :PreparedEventScheduler,
        :EventSchedulerState,
        :EventSchedulerWorkspace,
        :EventClaim,
        :prepare_event_scheduler,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
end

@testset "Deterministic scheduler scan and equal-time ordering" begin
    definitions = (
        EventGeneratorDefinition(
            PlantTimestamp(40), ReadoutCompletionPhase, 1),
        EventGeneratorDefinition(
            PlantTimestamp(10), OpticalSamplePhase, 2),
        EventGeneratorDefinition(
            PlantTimestamp(10), ExposureOpenPhase, 3),
        EventGeneratorDefinition(
            PlantTimestamp(10), OpticalSamplePhase, 1),
    )
    scheduler = prepare_event_scheduler(reverse(definitions); capacity=8)
    state = EventSchedulerState(scheduler)
    workspace = EventSchedulerWorkspace(scheduler)

    @test scan_due_events!(workspace, scheduler, state) == 3
    @test due_event_count(workspace, scheduler, state) == 3
    @test due_event_timestamp(workspace, scheduler, state) ==
        PlantTimestamp(10)
    scanned = ntuple(index -> due_event_key(
        workspace, scheduler, state, index), 3)
    @test Tuple(key.phase for key in scanned) == (
        ExposureOpenPhase,
        OpticalSamplePhase,
        OpticalSamplePhase,
    )
    @test Tuple(key.ordinal for key in scanned) ==
        (UInt32(3), UInt32(1), UInt32(2))

    forward_trace = collect_scheduler_trace(definitions)
    reverse_trace = collect_scheduler_trace(reverse(definitions))
    @test forward_trace == reverse_trace
    @test issorted(forward_trace)
    @test Tuple(key.timestamp for key in forward_trace) == (
        PlantTimestamp(10),
        PlantTimestamp(10),
        PlantTimestamp(10),
        PlantTimestamp(40),
    )

    first_claim = claim_next_event!(workspace, scheduler, state)
    @test claimed_event_key(first_claim) == scanned[1]
    @test event_generator_handle(first_claim) == event_generator_handle(
        scheduler, ExposureOpenPhase, 3)
    @test scheduler_timestamp(state) == PlantTimestamp(10)
    reschedule_event!(scheduler, state, first_claim, PlantTimestamp(50))
    second_claim = claim_next_event!(workspace, scheduler, state)
    @test claimed_event_key(second_claim) == scanned[2]
    deactivate_event_generator!(scheduler, state, second_claim)

    duplicate_definitions = (
        EventGeneratorDefinition(
            PlantTimestamp(100), OpticalSamplePhase, 1),
        EventGeneratorDefinition(
            PlantTimestamp(100), OpticalSamplePhase, 2),
    )
    duplicate_scheduler = prepare_event_scheduler(duplicate_definitions)
    duplicate_state = EventSchedulerState(duplicate_scheduler)
    duplicate_workspace = EventSchedulerWorkspace(duplicate_scheduler)
    claim_1 = claim_next_event!(duplicate_workspace, duplicate_scheduler,
        duplicate_state)
    @test claimed_event_key(claim_1).ordinal == 1
    same_time_key = reschedule_event!(duplicate_scheduler, duplicate_state,
        claim_1, PlantTimestamp(100))
    @test same_time_key.occurrence == 2
    claim_2 = claim_next_event!(duplicate_workspace, duplicate_scheduler,
        duplicate_state)
    @test claimed_event_key(claim_2) == same_time_key
    deactivate_event_generator!(duplicate_scheduler, duplicate_state,
        claim_2)
    claim_3 = claim_next_event!(duplicate_workspace, duplicate_scheduler,
        duplicate_state)
    @test claimed_event_key(claim_3).ordinal == 2

    long_period = PeriodicSchedule(period_ns=1_000_000_000_000,
        phase_ns=999_999_999_999)
    long_scheduler = prepare_event_scheduler((EventGeneratorDefinition(
        long_period, AtmosphereEvolutionPhase, 1),))
    long_state = EventSchedulerState(long_scheduler)
    long_workspace = EventSchedulerWorkspace(long_scheduler)
    long_claim = claim_next_event!(long_workspace, long_scheduler,
        long_state)
    @test claimed_event_key(long_claim).timestamp ==
        PlantTimestamp(999_999_999_999)
    next_long = reschedule_periodic_event!(long_scheduler, long_state,
        long_claim, long_period)
    @test next_long.timestamp == PlantTimestamp(1_999_999_999_999)
end

@testset "Scheduler activation, deactivation, and structural failures" begin
    active = EventGeneratorDefinition(
        PlantTimestamp(10), OpticalSamplePhase, 1)
    inactive = EventGeneratorDefinition(
        PlantTimestamp(20), AcquisitionReadyPhase, 1; active=false)
    scheduler = prepare_event_scheduler((active, inactive); capacity=2)
    state = EventSchedulerState(scheduler)
    workspace = EventSchedulerWorkspace(scheduler)
    active_handle = event_generator_handle(
        scheduler, OpticalSamplePhase, 1)
    inactive_handle = event_generator_handle(
        scheduler, AcquisitionReadyPhase, 1)

    @test deactivate_event_generator!(scheduler, state, active_handle) ===
        active_handle
    @test scan_due_events!(workspace, scheduler, state) == 0
    @test isnothing(due_event_timestamp(workspace, scheduler, state))
    @test isnothing(claim_next_event!(workspace, scheduler, state))
    @test activate_event_generator!(scheduler, state, inactive_handle,
        PlantTimestamp(30)).timestamp == PlantTimestamp(30)
    assert_event_scheduler_error(
        () -> activate_event_generator!(scheduler, state, inactive_handle,
            PlantTimestamp(31)),
        :invalid_transition,
    )

    @test scan_due_events!(workspace, scheduler, state) == 1
    assert_event_scheduler_error(
        () -> due_event_key(workspace, scheduler, state, 0),
        :invalid_due_index,
    )
    assert_event_scheduler_error(
        () -> due_event_key(workspace, scheduler, state, true),
        :invalid_due_index,
    )
    claim = claim_next_event!(workspace, scheduler, state)
    assert_event_scheduler_error(
        () -> scan_due_events!(workspace, scheduler, state),
        :outstanding_claim,
    )
    assert_event_scheduler_error(
        () -> activate_event_generator!(scheduler, state, active_handle,
            PlantTimestamp(40)),
        :outstanding_claim,
    )
    assert_event_scheduler_error(
        () -> reschedule_event!(scheduler, state, claim,
            PlantTimestamp(29)),
        :time_regression,
    )
    next_key = reschedule_event!(scheduler, state, claim,
        PlantTimestamp(40))
    @test next_key.occurrence == 2
    assert_event_scheduler_error(
        () -> reschedule_event!(scheduler, state, claim,
            PlantTimestamp(50)),
        :stale_claim,
    )
    assert_event_scheduler_error(
        () -> due_event_count(workspace, scheduler, state),
        :stale_due_scan,
    )

    claim_2 = claim_next_event!(workspace, scheduler, state)
    returned_handle = deactivate_event_generator!(scheduler, state, claim_2)
    @test activate_event_generator!(scheduler, state, returned_handle,
        PlantTimestamp(40)).occurrence == 3
    claim_3 = claim_next_event!(workspace, scheduler, state)
    @test claimed_event_key(claim_3).occurrence == 3
    deactivate_event_generator!(scheduler, state, claim_3)

    assert_event_scheduler_error(
        () -> activate_event_generator!(scheduler, state, active_handle,
            PlantTimestamp(20)),
        :time_regression,
    )

    other_scheduler = prepare_event_scheduler((active,))
    other_state = EventSchedulerState(other_scheduler)
    other_workspace = EventSchedulerWorkspace(other_scheduler)
    other_handle = event_generator_handle(
        other_scheduler, OpticalSamplePhase, 1)
    other_claim = claim_next_event!(other_workspace, other_scheduler,
        other_state)
    assert_event_scheduler_error(
        () -> scan_due_events!(workspace, scheduler, other_state),
        :foreign_state,
    )
    assert_event_scheduler_error(
        () -> scan_due_events!(other_workspace, scheduler, state),
        :foreign_workspace,
    )
    assert_event_scheduler_error(
        () -> activate_event_generator!(scheduler, state, other_handle,
            PlantTimestamp(100)),
        :foreign_generator,
    )
    assert_event_scheduler_error(
        () -> reschedule_event!(scheduler, state, other_claim,
            PlantTimestamp(100)),
        :foreign_claim,
    )

    overflow_definition = EventGeneratorDefinition(
        PlantTimestamp(1), OpticalSamplePhase, 1;
        occurrence=typemax(UInt64),
    )
    overflow_scheduler = prepare_event_scheduler((overflow_definition,))
    overflow_state = EventSchedulerState(overflow_scheduler)
    overflow_workspace = EventSchedulerWorkspace(overflow_scheduler)
    overflow_claim = claim_next_event!(overflow_workspace,
        overflow_scheduler, overflow_state)
    assert_event_scheduler_error(
        () -> reschedule_event!(overflow_scheduler, overflow_state,
            overflow_claim, PlantTimestamp(2)),
        :occurrence_overflow,
    )
end

@testset "Scheduler long-run storage, inference, and allocation" begin
    schedule = PeriodicSchedule(period_ns=17, phase_ns=3)
    scheduler = prepare_event_scheduler((EventGeneratorDefinition(
        schedule, OpticalSamplePhase, 1),); capacity=16)
    state = EventSchedulerState(scheduler)
    workspace = EventSchedulerWorkspace(scheduler)
    registry_size = Base.summarysize(getfield(scheduler, :definitions))
    cursor_size = Base.summarysize(state.cursors)
    workspace_size = Base.summarysize(workspace.due_slots)

    @test @inferred(scheduler_hot_cycle!(
        scheduler, state, workspace, schedule)) == PlantTimestamp(20)
    for _ in 1:100_000
        scheduler_hot_cycle!(scheduler, state, workspace, schedule)
    end
    @test scheduler_timestamp(state) ==
        schedule_timestamp(schedule, 100_001)
    @test event_generator_count(scheduler) == 1
    @test event_scheduler_capacity(scheduler) == 16
    @test Base.summarysize(getfield(scheduler, :definitions)) == registry_size
    @test Base.summarysize(state.cursors) == cursor_size
    @test Base.summarysize(workspace.due_slots) == workspace_size

    scheduler_hot_cycle!(scheduler, state, workspace, schedule)
    if coverage_instrumented()
        @test_skip "allocation gate disabled under coverage instrumentation"
    else
        @test @allocated(scheduler_hot_cycle!(
            scheduler, state, workspace, schedule)) == 0
    end
end
