const AOSTrigger = Plant

function assert_trigger_topology_error(f, reason::Symbol)
    try
        f()
        @test false
    catch error
        @test error isa PlantScheduleError
        if error isa PlantScheduleError
            @test error.component == :trigger_topology
            @test error.reason == reason
            @test !isempty(error.msg)
        end
    end
end

function collect_trigger_deliveries(topology, source_count)
    state = AOSTrigger.TriggerTopologyState(topology)
    workspace = AOSTrigger.TriggerTopologyWorkspace(topology)
    deliveries = AOSTrigger.TriggerDelivery[]
    realized = 0
    while realized < source_count ||
            AOSTrigger.pending_trigger_delivery_count(state) > 0
        next_delivery = AOSTrigger.next_trigger_delivery(topology, state)
        if realized < source_count
            next_source = AOSTrigger.next_trigger_source(topology, state)
            if next_delivery === nothing ||
                    AOSTrigger.realized_trigger_source_timestamp(next_source) <=
                        AOSTrigger.delivered_trigger_edge(next_delivery).timestamp
                AOSTrigger.realize_next_trigger_source!(workspace, topology,
                    state)
                realized += 1
                continue
            end
        end
        push!(deliveries,
            AOSTrigger.pop_next_trigger_delivery!(topology, state))
    end
    return deliveries
end

function trigger_hot_cycle!(workspace, topology, state, first_output,
    second_output)
    realization = AOSTrigger.realize_next_trigger_source!(workspace,
        topology, state)
    AOSTrigger.pop_next_trigger_delivery!(first_output, topology, state)
    AOSTrigger.pop_next_trigger_delivery!(second_output, topology, state)
    return plant_nanoseconds(
        AOSTrigger.realized_trigger_source_timestamp(realization)) +
        plant_nanoseconds(AOSTrigger.delivered_trigger_edge(
            first_output[]).timestamp) +
        plant_nanoseconds(AOSTrigger.delivered_trigger_edge(
            second_output[]).timestamp)
end

@testset "Trigger identities, traces, and prepared topology" begin
    source_id = AOSTrigger.TriggerSourceID(:camera_clock)
    link_id = AOSTrigger.TriggerLinkID(:fanout_a)
    consumer_id = AOSTrigger.TriggerConsumerID(:wfs_a)
    fault_id = AOSTrigger.TriggerFaultID(:source_step)
    @test source_id == AOSTrigger.TriggerSourceID(:camera_clock)
    @test isequal(source_id, AOSTrigger.TriggerSourceID(:camera_clock))
    @test source_id != AOSTrigger.TriggerSourceID(:other_clock)
    @test length(Set((source_id,
        AOSTrigger.TriggerSourceID(:camera_clock)))) == 1
    @test sprint(show, source_id) == "TriggerSourceID(:camera_clock)"
    @test sprint(show, link_id) == "TriggerLinkID(:fanout_a)"
    @test sprint(show, consumer_id) == "TriggerConsumerID(:wfs_a)"
    @test sprint(show, fault_id) == "TriggerFaultID(:source_step)"

    late_entry = AOSTrigger.TriggerFaultTraceEntry(3, :late;
        jitter=AOSTrigger.PlantTimeOffset(2))
    early_entry = AOSTrigger.TriggerFaultTraceEntry(2, fault_id;
        phase_step=AOSTrigger.PlantTimeOffset(5))
    trace = AOSTrigger.TriggerFaultTrace(late_entry, early_entry)
    @test typeof(trace) === typeof(AOSTrigger.TriggerFaultTrace(early_entry))
    source = AOSTrigger.TriggerSourceDefinition(source_id,
        PeriodicSchedule(period_ns=100); faults=trace)
    link = AOSTrigger.TriggerLinkDefinition(link_id, source_id;
        propagation_delay=PlantDuration(10))
    consumer = AOSTrigger.TriggerConsumerDefinition(consumer_id, link_id)
    topology = AOSTrigger.prepare_trigger_topology((source,), (link,),
        (consumer,); in_flight_capacity=1)

    @test !Base.ismutabletype(AOSTrigger.TriggerSourceID)
    @test !Base.ismutabletype(AOSTrigger.TriggerLinkID)
    @test !Base.ismutabletype(AOSTrigger.TriggerConsumerID)
    @test !Base.ismutabletype(AOSTrigger.TriggerFaultID)
    @test !Base.ismutabletype(AOSTrigger.TriggerFaultTraceEntry)
    @test !Base.ismutabletype(typeof(trace))
    @test isbitstype(AOSTrigger.TriggerSourceHandle)
    @test !Base.ismutabletype(AOSTrigger.NominalTriggerEdge)
    @test isbitstype(AOSTrigger.DeliveredTriggerEdge)
    @test isbitstype(AOSTrigger.ReportedTriggerTimestamp)
    @test !Base.ismutabletype(AOSTrigger.TriggerDelivery)
    @test AOSTrigger.trigger_source_id(source) == source_id
    @test AOSTrigger.trigger_link_id(link) == link_id
    @test AOSTrigger.trigger_consumer_id(consumer) == consumer_id
    @test AOSTrigger.trigger_parent_id(link) == source_id
    @test AOSTrigger.trigger_parent_id(consumer) == link_id
    @test AOSTrigger.trigger_source_count(topology) == 1
    @test AOSTrigger.trigger_link_count(topology) == 1
    @test AOSTrigger.trigger_consumer_count(topology) == 1
    @test AOSTrigger.trigger_in_flight_capacity(topology) == 1
    @test AOSTrigger.required_trigger_in_flight_capacity(topology) == 1
    @test AOSTrigger.trigger_source_handle(topology, :camera_clock) isa
        AOSTrigger.TriggerSourceHandle
    assert_trigger_topology_error(
        () -> AOSTrigger.trigger_source_handle(topology, :missing),
        :unknown_source,
    )
    @test getfield(topology, :sources) isa
        Memory{AOSTrigger._PreparedTriggerSource}
    @test getfield(topology, :links) isa
        Memory{AOSTrigger._PreparedTriggerLink}
    @test getfield(topology, :consumers) isa
        Memory{AOSTrigger._PreparedTriggerConsumer}

    @test Base.ispublic(Plant, :PreparedTriggerTopology)
    @test !Base.isexported(Plant, :PreparedTriggerTopology)
    @test Base.isexported(Plant, :TriggerSourceDefinition)
    @test Base.isexported(Plant, :prepare_trigger_topology)
    @test !Base.isexported(AdaptiveOpticsSim, :TriggerSourceDefinition)

    state = AOSTrigger.TriggerTopologyState(topology)
    workspace = AOSTrigger.TriggerTopologyWorkspace(topology)
    AOSTrigger.realize_next_trigger_source!(workspace, topology, state)
    @test AOSTrigger.next_trigger_delivery_timestamp(topology, state) ==
        PlantTimestamp(10)
    output = Ref{AOSTrigger.TriggerDelivery}()
    @test AOSTrigger.pop_next_trigger_delivery!(output, topology, state)
    @test AOSTrigger.trigger_delivery_consumer(output[]) == consumer_id
    @test !AOSTrigger.pop_next_trigger_delivery!(output, topology, state)
    assert_trigger_topology_error(
        () -> AOSTrigger.contains_trigger_fault(topology,
            AOSTrigger.trigger_delivery_faults(output[]), :missing),
        :unknown_fault,
    )

    named_topology = AOSTrigger.prepare_trigger_topology(
        (camera_clock=source,), (fanout_a=link,), (wfs_a=consumer,);
        in_flight_capacity=1)
    @test AOSTrigger.trigger_source_count(named_topology) == 1
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology(
            (wrong_name=source,), (fanout_a=link,), (wfs_a=consumer,);
            in_flight_capacity=1),
        :identity_mismatch,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.TriggerSourceDefinition("invalid",
            PeriodicSchedule(period_ns=1)),
        :invalid_id,
    )
end

@testset "Equal-time trigger source and delivery ordering" begin
    source_a = AOSTrigger.TriggerSourceDefinition(:source_a,
        PeriodicSchedule(period_ns=100))
    source_b = AOSTrigger.TriggerSourceDefinition(:source_b,
        PeriodicSchedule(period_ns=100))
    consumer_a = AOSTrigger.TriggerConsumerDefinition(:camera_a,
        AOSTrigger.TriggerSourceID(:source_a))
    consumer_b = AOSTrigger.TriggerConsumerDefinition(:camera_b,
        AOSTrigger.TriggerSourceID(:source_b))
    topology = AOSTrigger.prepare_trigger_topology((source_b, source_a), (),
        (consumer_b, consumer_a); in_flight_capacity=2)
    state = AOSTrigger.TriggerTopologyState(topology)
    workspace = AOSTrigger.TriggerTopologyWorkspace(topology)
    first_source = AOSTrigger.realize_next_trigger_source!(workspace,
        topology, state)
    second_source = AOSTrigger.realize_next_trigger_source!(workspace,
        topology, state)
    @test AOSTrigger.nominal_trigger_edge(first_source).source_id ==
        AOSTrigger.TriggerSourceID(:source_a)
    @test AOSTrigger.nominal_trigger_edge(second_source).source_id ==
        AOSTrigger.TriggerSourceID(:source_b)
    @test AOSTrigger.trigger_delivery_consumer(
        AOSTrigger.pop_next_trigger_delivery!(topology, state)) ==
        AOSTrigger.TriggerConsumerID(:camera_a)
    @test AOSTrigger.trigger_delivery_consumer(
        AOSTrigger.pop_next_trigger_delivery!(topology, state)) ==
        AOSTrigger.TriggerConsumerID(:camera_b)

    duplicate = AOSTrigger.TriggerFaultTraceEntry(1, :link_duplicate;
        action=AOSTrigger.DuplicateTriggerEdge)
    source = AOSTrigger.TriggerSourceDefinition(:clock,
        PeriodicSchedule(period_ns=100))
    link = AOSTrigger.TriggerLinkDefinition(:fanout,
        AOSTrigger.TriggerSourceID(:clock);
        faults=AOSTrigger.TriggerFaultTrace(duplicate))
    consumer = AOSTrigger.TriggerConsumerDefinition(:camera,
        AOSTrigger.TriggerLinkID(:fanout))
    duplicate_topology = AOSTrigger.prepare_trigger_topology((source,),
        (link,), (consumer,); in_flight_capacity=2)
    duplicate_state = AOSTrigger.TriggerTopologyState(duplicate_topology)
    duplicate_workspace = AOSTrigger.TriggerTopologyWorkspace(
        duplicate_topology)
    AOSTrigger.realize_next_trigger_source!(duplicate_workspace,
        duplicate_topology, duplicate_state)
    primary = AOSTrigger.pop_next_trigger_delivery!(duplicate_topology,
        duplicate_state)
    repeated = AOSTrigger.pop_next_trigger_delivery!(duplicate_topology,
        duplicate_state)
    @test AOSTrigger.delivered_trigger_edge(primary).timestamp ==
        AOSTrigger.delivered_trigger_edge(repeated).timestamp
    @test AOSTrigger.delivered_trigger_edge(primary).occurrence == 1
    @test AOSTrigger.delivered_trigger_edge(repeated).occurrence == 2
    @test !AOSTrigger.delivered_trigger_edge(primary).duplicate
    @test AOSTrigger.delivered_trigger_edge(repeated).duplicate
    @test AOSTrigger.contains_trigger_fault(duplicate_topology,
        AOSTrigger.trigger_delivery_faults(repeated), :link_duplicate)
end

@testset "Correlated source and isolated branch trigger faults" begin
    source_fault = AOSTrigger.TriggerFaultTraceEntry(2,
        :common_phase_duplicate;
        jitter=AOSTrigger.PlantTimeOffset(2),
        phase_step=AOSTrigger.PlantTimeOffset(5),
        timestamp_label_offset=AOSTrigger.PlantTimeOffset(100),
        action=AOSTrigger.DuplicateTriggerEdge,
        duplicate_delay=PlantDuration(1))
    branch_drop = AOSTrigger.TriggerFaultTraceEntry(2, :branch_b_drop;
        action=AOSTrigger.DropTriggerEdge)
    branch_jitter = AOSTrigger.TriggerFaultTraceEntry(3, :branch_a_jitter;
        jitter=AOSTrigger.PlantTimeOffset(5))
    source = AOSTrigger.TriggerSourceDefinition(:clock,
        PeriodicSchedule(period_ns=100);
        faults=AOSTrigger.TriggerFaultTrace(source_fault))
    link_a = AOSTrigger.TriggerLinkDefinition(:branch_a,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(10),
        timestamp_label_offset=AOSTrigger.PlantTimeOffset(2),
        faults=AOSTrigger.TriggerFaultTrace(branch_jitter))
    link_b = AOSTrigger.TriggerLinkDefinition(:branch_b,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(20),
        faults=AOSTrigger.TriggerFaultTrace(branch_drop))
    nested = AOSTrigger.TriggerLinkDefinition(:nested,
        AOSTrigger.TriggerLinkID(:branch_b);
        propagation_delay=PlantDuration(3))
    consumer_a = AOSTrigger.TriggerConsumerDefinition(:camera_a,
        AOSTrigger.TriggerLinkID(:branch_a))
    consumer_b = AOSTrigger.TriggerConsumerDefinition(:camera_b,
        AOSTrigger.TriggerLinkID(:branch_b))
    consumer_nested = AOSTrigger.TriggerConsumerDefinition(:modulator,
        AOSTrigger.TriggerLinkID(:nested))

    forward = AOSTrigger.prepare_trigger_topology(
        (source,), (link_a, link_b, nested),
        (consumer_a, consumer_b, consumer_nested);
        in_flight_capacity=6)
    reversed = AOSTrigger.prepare_trigger_topology(
        (source,), (nested, link_b, link_a),
        (consumer_nested, consumer_b, consumer_a);
        in_flight_capacity=6)
    forward_deliveries = collect_trigger_deliveries(forward, 3)
    reverse_deliveries = collect_trigger_deliveries(reversed, 3)
    @test forward_deliveries == reverse_deliveries

    sequence_1 = filter(delivery ->
        AOSTrigger.nominal_trigger_edge(delivery).sequence == 1,
        forward_deliveries)
    @test Tuple(AOSTrigger.trigger_delivery_consumer(d).name for d in
        sequence_1) == (:camera_a, :camera_b, :modulator)
    @test Tuple(plant_nanoseconds(
        AOSTrigger.delivered_trigger_edge(d).timestamp) for d in sequence_1) ==
        (10, 20, 23)
    @test Tuple(plant_nanoseconds(
        AOSTrigger.reported_trigger_timestamp(d).timestamp)
        for d in sequence_1) == (12, 20, 23)

    sequence_2 = filter(delivery ->
        AOSTrigger.nominal_trigger_edge(delivery).sequence == 2,
        forward_deliveries)
    @test length(sequence_2) == 2
    @test all(AOSTrigger.trigger_delivery_consumer(d).name == :camera_a
        for d in sequence_2)
    @test Tuple(plant_nanoseconds(
        AOSTrigger.delivered_trigger_edge(d).timestamp) for d in sequence_2) ==
        (117, 118)
    @test Tuple(AOSTrigger.delivered_trigger_edge(d).occurrence
        for d in sequence_2) == (UInt64(1), UInt64(2))
    @test Tuple(AOSTrigger.delivered_trigger_edge(d).duplicate
        for d in sequence_2) == (false, true)
    @test Tuple(plant_nanoseconds(
        AOSTrigger.reported_trigger_timestamp(d).timestamp)
        for d in sequence_2) == (219, 220)
    @test all(AOSTrigger.contains_trigger_fault(forward,
        AOSTrigger.trigger_delivery_faults(d), :common_phase_duplicate)
        for d in sequence_2)
    @test all(!AOSTrigger.contains_trigger_fault(forward,
        AOSTrigger.trigger_delivery_faults(d), :branch_b_drop)
        for d in sequence_2)

    sequence_3 = filter(delivery ->
        AOSTrigger.nominal_trigger_edge(delivery).sequence == 3,
        forward_deliveries)
    @test Tuple(AOSTrigger.trigger_delivery_consumer(d).name for d in
        sequence_3) == (:camera_a, :camera_b, :modulator)
    @test Tuple(plant_nanoseconds(
        AOSTrigger.delivered_trigger_edge(d).timestamp) for d in sequence_3) ==
        (220, 225, 228)
    camera_a = only(filter(d ->
        AOSTrigger.trigger_delivery_consumer(d).name == :camera_a,
        sequence_3))
    camera_b = only(filter(d ->
        AOSTrigger.trigger_delivery_consumer(d).name == :camera_b,
        sequence_3))
    @test AOSTrigger.contains_trigger_fault(forward,
        AOSTrigger.trigger_delivery_faults(camera_a), :branch_a_jitter)
    @test !AOSTrigger.contains_trigger_fault(forward,
        AOSTrigger.trigger_delivery_faults(camera_b), :branch_a_jitter)

    state = AOSTrigger.TriggerTopologyState(forward)
    workspace = AOSTrigger.TriggerTopologyWorkspace(forward)
    AOSTrigger.realize_next_trigger_source!(workspace, forward, state)
    while AOSTrigger.pending_trigger_delivery_count(state) > 0
        AOSTrigger.pop_next_trigger_delivery!(forward, state)
    end
    AOSTrigger.realize_next_trigger_source!(workspace, forward, state)
    @test AOSTrigger.trigger_fault_observation_count(workspace) == 2
    observations = ntuple(index ->
        AOSTrigger.trigger_fault_observation(workspace, index), 2)
    @test Set(AOSTrigger.trigger_fault_id(forward, observation).name
        for observation in observations) ==
        Set((:common_phase_duplicate, :branch_b_drop))
    @test Set(AOSTrigger.trigger_fault_location(forward, observation).name
        for observation in observations) == Set((:clock, :branch_b))
    assert_trigger_topology_error(
        () -> AOSTrigger.trigger_fault_id(reversed, first(observations)),
        :foreign_observation,
    )

    other_source = AOSTrigger.TriggerSourceDefinition(:other_clock,
        PeriodicSchedule(period_ns=100);
        faults=AOSTrigger.TriggerFaultTrace(
            AOSTrigger.TriggerFaultTraceEntry(1, :different_fault)))
    other_consumer = AOSTrigger.TriggerConsumerDefinition(:other_camera,
        AOSTrigger.TriggerSourceID(:other_clock))
    other = AOSTrigger.prepare_trigger_topology((other_source,), (),
        (other_consumer,); in_flight_capacity=1)
    assert_trigger_topology_error(
        () -> AOSTrigger.contains_trigger_fault(other,
            AOSTrigger.trigger_delivery_faults(first(sequence_2)),
            :different_fault),
        :foreign_fault_set,
    )
end

@testset "Trigger chronology, bounds, and structural failures" begin
    source = AOSTrigger.TriggerSourceDefinition(:clock,
        PeriodicSchedule(period_ns=10))
    delayed = AOSTrigger.TriggerLinkDefinition(:delayed,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(20))
    consumer = AOSTrigger.TriggerConsumerDefinition(:camera,
        AOSTrigger.TriggerLinkID(:delayed))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (delayed,),
            (consumer,); in_flight_capacity=2),
        :insufficient_in_flight_capacity,
    )
    topology = AOSTrigger.prepare_trigger_topology((source,), (delayed,),
        (consumer,); in_flight_capacity=3)
    @test AOSTrigger.required_trigger_in_flight_capacity(topology) == 3
    state = AOSTrigger.TriggerTopologyState(topology)
    workspace = AOSTrigger.TriggerTopologyWorkspace(topology)
    AOSTrigger.realize_next_trigger_source!(workspace, topology, state)
    assert_trigger_topology_error(
        () -> AOSTrigger.pop_next_trigger_delivery!(topology, state),
        :source_due,
    )
    AOSTrigger.realize_next_trigger_source!(workspace, topology, state)
    AOSTrigger.realize_next_trigger_source!(workspace, topology, state)
    @test AOSTrigger.pending_trigger_delivery_count(state) == 3
    assert_trigger_topology_error(
        () -> AOSTrigger.realize_next_trigger_source!(workspace, topology,
            state),
        :delivery_due,
    )
    @test plant_nanoseconds(AOSTrigger.delivered_trigger_edge(
        AOSTrigger.pop_next_trigger_delivery!(topology, state)).timestamp) == 20

    overtaking = AOSTrigger.TriggerFaultTraceEntry(2, :early;
        jitter=AOSTrigger.PlantTimeOffset(-100))
    bad_link = AOSTrigger.TriggerLinkDefinition(:bad,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(100),
        faults=AOSTrigger.TriggerFaultTrace(overtaking))
    bad_consumer = AOSTrigger.TriggerConsumerDefinition(:bad_camera,
        AOSTrigger.TriggerLinkID(:bad))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (bad_link,),
            (bad_consumer,); in_flight_capacity=20),
        :ordinary_edge_overtaking,
    )

    duplicate_trace = AOSTrigger.TriggerFaultTrace(
        AOSTrigger.TriggerFaultTraceEntry(2, :a),
        AOSTrigger.TriggerFaultTraceEntry(2, :b),
    )
    duplicate_source = AOSTrigger.TriggerSourceDefinition(:duplicate,
        PeriodicSchedule(period_ns=100); faults=duplicate_trace)
    direct_consumer = AOSTrigger.TriggerConsumerDefinition(:direct,
        AOSTrigger.TriggerSourceID(:duplicate))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((duplicate_source,), (),
            (direct_consumer,); in_flight_capacity=1),
        :duplicate_fault_sequence,
    )

    duplicate_id = AOSTrigger.TriggerSourceDefinition(:clock,
        PeriodicSchedule(period_ns=20))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source, duplicate_id), (),
            (AOSTrigger.TriggerConsumerDefinition(:direct,
                AOSTrigger.TriggerSourceID(:clock)),);
            in_flight_capacity=1),
        :duplicate_id,
    )
    missing_link = AOSTrigger.TriggerLinkDefinition(:missing,
        AOSTrigger.TriggerLinkID(:absent))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (missing_link,),
            (consumer,); in_flight_capacity=3),
        :unknown_parent,
    )
    cycle_a = AOSTrigger.TriggerLinkDefinition(:cycle_a,
        AOSTrigger.TriggerLinkID(:cycle_b))
    cycle_b = AOSTrigger.TriggerLinkDefinition(:cycle_b,
        AOSTrigger.TriggerLinkID(:cycle_a))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,),
            (cycle_a, cycle_b),
            (AOSTrigger.TriggerConsumerDefinition(:cycle_consumer,
                AOSTrigger.TriggerLinkID(:cycle_a)),);
            in_flight_capacity=3),
        :cyclic_topology,
    )
    orphan = AOSTrigger.TriggerLinkDefinition(:orphan,
        AOSTrigger.TriggerSourceID(:clock))
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,),
            (delayed, orphan), (consumer,); in_flight_capacity=3),
        :unconsumed_link,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((), (),
            (direct_consumer,); in_flight_capacity=1),
        :empty_sources,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (), ();
            in_flight_capacity=1),
        :empty_consumers,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (),
            (AOSTrigger.TriggerConsumerDefinition(:direct,
                AOSTrigger.TriggerSourceID(:clock)),);
            in_flight_capacity=true),
        :invalid_capacity,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.prepare_trigger_topology((source,), (),
            (AOSTrigger.TriggerConsumerDefinition(:direct,
                AOSTrigger.TriggerSourceID(:clock)),);
            in_flight_capacity="one"),
        :invalid_capacity,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.TriggerFaultTraceEntry(1, :invalid;
            duplicate_delay=PlantDuration(1)),
        :invalid_fault_trace,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.TriggerFaultTraceEntry(0, :invalid),
        :invalid_sequence,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger.TriggerFaultTraceEntry("one", :invalid),
        :invalid_sequence,
    )

    direct = AOSTrigger.prepare_trigger_topology((source,), (),
        (AOSTrigger.TriggerConsumerDefinition(:direct,
            AOSTrigger.TriggerSourceID(:clock)),);
        in_flight_capacity=1)
    corrupted_state = AOSTrigger.TriggerTopologyState(direct)
    corrupted_state.last_ordinary_delivery[1] = PlantTimestamp(1_000)
    corrupted_state.has_last_ordinary_delivery[1] = true
    assert_trigger_topology_error(
        () -> AOSTrigger.realize_next_trigger_source!(
            AOSTrigger.TriggerTopologyWorkspace(direct), direct,
            corrupted_state),
        :ordinary_edge_overtaking,
    )
    assert_trigger_topology_error(
        () -> AOSTrigger._trigger_fault_slot(
            AOSTrigger.TriggerFaultID[], AOSTrigger.TriggerFaultID(:missing)),
        :unknown_fault,
    )
end

@testset "Trigger long-run storage, inference, and allocation" begin
    source = AOSTrigger.TriggerSourceDefinition(:clock,
        PeriodicSchedule(period_ns=100))
    fast_link = AOSTrigger.TriggerLinkDefinition(:fast,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(10))
    slow_link = AOSTrigger.TriggerLinkDefinition(:slow,
        AOSTrigger.TriggerSourceID(:clock);
        propagation_delay=PlantDuration(20))
    fast_consumer = AOSTrigger.TriggerConsumerDefinition(:fast_camera,
        AOSTrigger.TriggerLinkID(:fast))
    slow_consumer = AOSTrigger.TriggerConsumerDefinition(:slow_camera,
        AOSTrigger.TriggerLinkID(:slow))
    topology = AOSTrigger.prepare_trigger_topology((source,),
        (slow_link, fast_link), (slow_consumer, fast_consumer);
        in_flight_capacity=2)
    state = AOSTrigger.TriggerTopologyState(topology)
    workspace = AOSTrigger.TriggerTopologyWorkspace(topology)
    first_output = Ref{AOSTrigger.TriggerDelivery}()
    second_output = Ref{AOSTrigger.TriggerDelivery}()

    @test @inferred(trigger_hot_cycle!(workspace, topology, state,
        first_output, second_output)) == 30
    trigger_hot_cycle!(workspace, topology, state, first_output,
        second_output)
    if coverage_instrumented()
        @test_skip "trigger allocation assertion is disabled under coverage"
    else
        @test @allocated(trigger_hot_cycle!(workspace, topology, state,
            first_output, second_output)) == 0
    end
    pending_storage = Base.summarysize(state.pending) +
        Base.summarysize(state.pending_active)
    for _ in 1:100_000
        trigger_hot_cycle!(workspace, topology, state, first_output,
            second_output)
    end
    @test Base.summarysize(state.pending) +
        Base.summarysize(state.pending_active) ==
        pending_storage
    @test length(state.pending) == 2
    @test length(workspace.propagation) == 3
    @test length(workspace.staged_deliveries) == 2
    @test AOSTrigger.pending_trigger_delivery_count(state) == 0
end
