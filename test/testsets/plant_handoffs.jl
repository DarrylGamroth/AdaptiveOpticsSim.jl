function assert_handoff_preparation_error(f, reason::Symbol)
    error = captured_handoff_preparation_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === :handoff
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function handoff_test_buffers(capacity::Int=2)
    sources = ntuple(_ -> zeros(Float64, 2, 3), capacity)
    destinations = ntuple(
        _ -> HandoffTestArray(fill(-1.0, 2, 3)), capacity)
    return sources, destinations
end

function handoff_test_slot_output(handoff)
    I = typeof(handoff_preparation_identity(handoff))
    return Ref{HandoffSlotReference{I}}()
end

function handoff_allocation_sample!(handoff, contract)
    reference_output = handoff_test_slot_output(handoff)
    reference_bytes = @allocated begin
        @assert try_next_free_handoff_slot!(reference_output, handoff) ==
            HandoffTransitionSucceeded
    end
    reference = reference_output[]
    source_bytes = @allocated producer_handoff_payload(handoff, reference)
    source = producer_handoff_payload(handoff, reference)
    fill!(source, 7.0)
    publication = HandoffTestPublication(contract.identity, UInt64(99))
    submit_bytes = @allocated begin
        @assert submit_handoff!(handoff, reference, publication) ==
            HandoffTransitionSucceeded
    end
    complete_bytes = @allocated begin
        @assert complete_handoff!(handoff, reference) ==
            HandoffTransitionSucceeded
    end
    payload_output = Ref{typeof(getfield(handoff, :destination_slots)[1])}()
    publication_output = Ref{typeof(publication)}()
    borrow_bytes = @allocated begin
        @assert try_borrow_completed_handoff!(
            payload_output, publication_output, handoff, reference) ==
            HandoffTransitionSucceeded
    end
    reclaim_bytes = @allocated begin
        @assert reclaim_handoff!(handoff, reference) ==
            HandoffTransitionSucceeded
    end
    return (
        reference=reference_bytes,
        source=source_bytes,
        submit=submit_bytes,
        complete=complete_bytes,
        borrow=borrow_bytes,
        reclaim=reclaim_bytes,
    )
end

@testset "Typed bounded cross-domain handoffs" begin
    reset_handoff_test_transfer_controls!()
    partitions = handoff_test_partitions()
    contract_identity = HandoffTestContractIdentity()
    contract = HandoffTestContract(
        contract_identity, (Base.OneTo(2), Base.OneTo(3)))
    sources, destinations = handoff_test_buffers()
    handoff = prepare_cross_domain_handoff(
        partitions, contract, sources, destinations)

    @test handoff isa PreparedCrossDomainHandoff
    @test handoff_contract(handoff) === contract
    @test getfield(handoff_preparation_identity(handoff), :partitions) ===
        partitions
    @test handoff_source_target(handoff) == HostComputeDevice()
    @test handoff_destination_target(handoff) == HANDOFF_TEST_ACCELERATOR
    @test handoff_capacity(handoff) == 2
    @test handoff_payload_bytes(handoff) == 6 * sizeof(Float64)
    @test_throws ErrorException destinations[1][1]
    @test_throws ErrorException setfield!(
        handoff_preparation_identity(handoff), :partitions, partitions)
    @test :PreparedCrossDomainHandoff in names(Plant)
    @test :PreparedCrossDomainHandoff ∉ names(AdaptiveOpticsSim)
    @test isconcretetype(typeof(handoff))

    reference_output = handoff_test_slot_output(handoff)
    @test @inferred(try_next_free_handoff_slot!(
        reference_output, handoff)) == HandoffTransitionSucceeded
    reference = reference_output[]
    @test reference.slot == UInt32(1)
    @test reference.generation == UInt64(1)
    @test reference.identity === handoff_preparation_identity(handoff)
    @test @inferred(handoff_slot_status(handoff, reference)) ==
        HandoffSlotFree
    @test @inferred(producer_handoff_payload(handoff, reference)) === sources[1]

    payload_output = Ref{typeof(destinations[1])}()
    publication = HandoffTestPublication(contract_identity, UInt64(1))
    publication_output = Ref{typeof(publication)}()
    @test try_borrow_completed_handoff!(
        payload_output, publication_output, handoff, reference) ==
        HandoffPayloadNotVisible
    @test reclaim_handoff!(handoff, reference) == HandoffSlotNotTerminal

    foreign_publication = HandoffTestPublication(
        HandoffTestContractIdentity(), UInt64(1))
    @test @inferred(submit_handoff!(
        handoff, reference, foreign_publication)) ==
        HandoffPublicationRejected
    @test submit_handoff!(handoff, reference, :wrong_type) ==
        HandoffPublicationRejected
    @test handoff_slot_status(handoff, reference) == HandoffSlotFree

    copyto!(producer_handoff_payload(handoff, reference),
        reshape(collect(1.0:6.0), 2, 3))
    HANDOFF_TEST_PENDING_OBSERVATIONS[] = UInt8(1)
    @test @inferred(submit_handoff!(handoff, reference, publication)) ==
        HandoffTransitionSucceeded
    @test handoff_slot_status(handoff, reference) ==
        HandoffTransferSubmitted
    @test submit_handoff!(handoff, reference, publication) ==
        HandoffSlotNotFree
    @test try_borrow_completed_handoff!(
        payload_output, publication_output, handoff, reference) ==
        HandoffPayloadNotVisible
    @test @inferred(try_complete_handoff!(handoff, reference)) ==
        HandoffTransferPending
    @test handoff_slot_status(handoff, reference) ==
        HandoffTransferSubmitted
    @test try_complete_handoff!(handoff, reference) ==
        HandoffTransitionSucceeded
    @test handoff_slot_status(handoff, reference) ==
        HandoffTransferCompleted
    @test try_complete_handoff!(handoff, reference) ==
        DuplicateHandoffCompletion
    @test @inferred(try_borrow_completed_handoff!(
        payload_output, publication_output, handoff, reference)) ==
        HandoffTransitionSucceeded
    @test payload_output[] === destinations[1]
    @test publication_output[] === publication
    @test payload_output[].storage == reshape(collect(1.0:6.0), 2, 3)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
    @test HANDOFF_TEST_CONTEXT_ENTRIES[] >= UInt64(3)

    @test @inferred(reclaim_handoff!(handoff, reference)) ==
        HandoffTransitionSucceeded
    @test submit_handoff!(handoff, reference, publication) ==
        StaleHandoffSlot
    @test_throws PlantPreparationError handoff_slot_status(handoff, reference)
    @test try_next_free_handoff_slot!(reference_output, handoff) ==
        HandoffTransitionSucceeded
    reused_reference = reference_output[]
    @test reused_reference.slot == reference.slot
    @test reused_reference.generation == reference.generation + one(UInt64)

    blocking_sources, blocking_destinations = handoff_test_buffers(1)
    blocking_handoff = prepare_cross_domain_handoff(
        partitions, contract, blocking_sources, blocking_destinations)
    blocking_output = handoff_test_slot_output(blocking_handoff)
    @test try_next_free_handoff_slot!(blocking_output, blocking_handoff) ==
        HandoffTransitionSucceeded
    blocking_reference = blocking_output[]
    fill!(producer_handoff_payload(
        blocking_handoff, blocking_reference), 8.0)
    HANDOFF_TEST_PENDING_OBSERVATIONS[] = typemax(UInt8)
    @test submit_handoff!(
        blocking_handoff, blocking_reference, publication) ==
        HandoffTransitionSucceeded
    @test @inferred(complete_handoff!(
        blocking_handoff, blocking_reference)) ==
        HandoffTransitionSucceeded
    @test handoff_slot_status(blocking_handoff, blocking_reference) ==
        HandoffTransferCompleted
    @test reclaim_handoff!(blocking_handoff, blocking_reference) ==
        HandoffTransitionSucceeded
    HANDOFF_TEST_PENDING_OBSERVATIONS[] = zero(UInt8)

    second_sources, second_destinations = handoff_test_buffers()
    second_handoff = prepare_cross_domain_handoff(
        partitions, contract, second_sources, second_destinations)
    foreign_output = handoff_test_slot_output(second_handoff)
    @test try_next_free_handoff_slot!(foreign_output, second_handoff) ==
        HandoffTransitionSucceeded
    @test submit_handoff!(handoff, foreign_output[], publication) ==
        ForeignHandoffPreparation
    invalid_slot = HandoffSlotReference(
        handoff_preparation_identity(handoff), UInt32(0), UInt64(1))
    invalid_generation = HandoffSlotReference(
        handoff_preparation_identity(handoff), UInt32(1), UInt64(0))
    @test submit_handoff!(handoff, invalid_slot, publication) ==
        InvalidHandoffSlot
    @test submit_handoff!(handoff, invalid_generation, publication) ==
        InvalidHandoffGeneration
end

@testset "Handoff capacity and transfer failures" begin
    reset_handoff_test_transfer_controls!()
    partitions = handoff_test_partitions()
    contract_identity = HandoffTestContractIdentity()
    contract = HandoffTestContract(
        contract_identity, (Base.OneTo(2), Base.OneTo(3)))
    sources, destinations = handoff_test_buffers()
    handoff = prepare_cross_domain_handoff(
        partitions, contract, sources, destinations)
    first_output = handoff_test_slot_output(handoff)
    second_output = handoff_test_slot_output(handoff)
    full_output = handoff_test_slot_output(handoff)
    publication = HandoffTestPublication(contract_identity, UInt64(2))

    @test try_next_free_handoff_slot!(first_output, handoff) ==
        HandoffTransitionSucceeded
    first_reference = first_output[]
    fill!(producer_handoff_payload(handoff, first_reference), 2.0)
    @test submit_handoff!(handoff, first_reference, publication) ==
        HandoffTransitionSucceeded
    @test try_next_free_handoff_slot!(second_output, handoff) ==
        HandoffTransitionSucceeded
    second_reference = second_output[]
    @test second_reference.slot != first_reference.slot
    fill!(producer_handoff_payload(handoff, second_reference), 3.0)
    @test submit_handoff!(handoff, second_reference, publication) ==
        HandoffTransitionSucceeded
    @test try_next_free_handoff_slot!(full_output, handoff) ==
        HandoffSlotsExhausted
    for reference in (first_reference, second_reference)
        @test try_complete_handoff!(handoff, reference) ==
            HandoffTransitionSucceeded
        @test reclaim_handoff!(handoff, reference) ==
            HandoffTransitionSucceeded
    end

    failure_output = handoff_test_slot_output(handoff)
    @test try_next_free_handoff_slot!(failure_output, handoff) ==
        HandoffTransitionSucceeded
    failure_reference = failure_output[]
    fill!(producer_handoff_payload(handoff, failure_reference), 4.0)
    HANDOFF_TEST_FAIL_SUBMISSION[] = true
    @test submit_handoff!(handoff, failure_reference, publication) ==
        HandoffSubmissionFailed
    HANDOFF_TEST_FAIL_SUBMISSION[] = false
    @test handoff_slot_status(handoff, failure_reference) ==
        HandoffTransferFailed
    @test handoff_slot_failure_reason(handoff, failure_reference) ==
        :submission_failed
    @test try_complete_handoff!(handoff, failure_reference) ==
        HandoffAlreadyFailed
    failed_payload_output = Ref{typeof(destinations[1])}()
    failed_publication_output = Ref{typeof(publication)}()
    @test try_borrow_completed_handoff!(
        failed_payload_output,
        failed_publication_output,
        handoff,
        failure_reference,
    ) == HandoffAlreadyFailed
    @test reclaim_handoff!(handoff, failure_reference) ==
        HandoffTransitionSucceeded
    @test reclaim_handoff!(handoff, failure_reference) == StaleHandoffSlot

    completion_output = handoff_test_slot_output(handoff)
    @test try_next_free_handoff_slot!(completion_output, handoff) ==
        HandoffTransitionSucceeded
    completion_reference = completion_output[]
    fill!(producer_handoff_payload(handoff, completion_reference), 5.0)
    @test submit_handoff!(handoff, completion_reference, publication) ==
        HandoffTransitionSucceeded
    HANDOFF_TEST_FAIL_COMPLETION[] = true
    @test try_complete_handoff!(handoff, completion_reference) ==
        HandoffCompletionFailed
    HANDOFF_TEST_FAIL_COMPLETION[] = false
    @test handoff_slot_failure_reason(handoff, completion_reference) ==
        :completion_failed
    @test reclaim_handoff!(handoff, completion_reference) ==
        HandoffTransitionSucceeded

    blocking_failure_sources, blocking_failure_destinations =
        handoff_test_buffers(1)
    blocking_failure_handoff = prepare_cross_domain_handoff(
        partitions,
        contract,
        blocking_failure_sources,
        blocking_failure_destinations,
    )
    blocking_failure_output =
        handoff_test_slot_output(blocking_failure_handoff)
    @test try_next_free_handoff_slot!(
        blocking_failure_output, blocking_failure_handoff) ==
        HandoffTransitionSucceeded
    blocking_failure_reference = blocking_failure_output[]
    @test submit_handoff!(
        blocking_failure_handoff,
        blocking_failure_reference,
        publication,
    ) == HandoffTransitionSucceeded
    HANDOFF_TEST_FAIL_COMPLETION[] = true
    @test complete_handoff!(
        blocking_failure_handoff, blocking_failure_reference) ==
        HandoffCompletionFailed
    HANDOFF_TEST_FAIL_COMPLETION[] = false
    @test handoff_slot_failure_reason(
        blocking_failure_handoff, blocking_failure_reference) ==
        :synchronous_completion_failed
    @test reclaim_handoff!(
        blocking_failure_handoff, blocking_failure_reference) ==
        HandoffTransitionSucceeded

    uncertain_submission_output = handoff_test_slot_output(handoff)
    @test try_next_free_handoff_slot!(
        uncertain_submission_output, handoff) == HandoffTransitionSucceeded
    uncertain_submission_reference = uncertain_submission_output[]
    HANDOFF_TEST_THROW_SUBMISSION[] = true
    @test_throws ErrorException submit_handoff!(
        handoff, uncertain_submission_reference, publication)
    HANDOFF_TEST_THROW_SUBMISSION[] = false
    @test handoff_slot_status(handoff, uncertain_submission_reference) ==
        HandoffTransferUncertain
    @test handoff_slot_failure_reason(
        handoff, uncertain_submission_reference) == :submission_exception
    @test try_complete_handoff!(handoff, uncertain_submission_reference) ==
        HandoffAlreadyUncertain
    @test try_borrow_completed_handoff!(
        failed_payload_output,
        failed_publication_output,
        handoff,
        uncertain_submission_reference,
    ) == HandoffAlreadyUncertain
    @test reclaim_handoff!(handoff, uncertain_submission_reference) ==
        HandoffSlotNotTerminal

    uncertain_sources, uncertain_destinations = handoff_test_buffers(1)
    uncertain_completion_handoff = prepare_cross_domain_handoff(
        partitions, contract, uncertain_sources, uncertain_destinations)
    uncertain_completion_output =
        handoff_test_slot_output(uncertain_completion_handoff)
    @test try_next_free_handoff_slot!(
        uncertain_completion_output, uncertain_completion_handoff) ==
        HandoffTransitionSucceeded
    uncertain_completion_reference = uncertain_completion_output[]
    @test submit_handoff!(
        uncertain_completion_handoff,
        uncertain_completion_reference,
        publication,
    ) == HandoffTransitionSucceeded
    HANDOFF_TEST_THROW_COMPLETION[] = true
    @test_throws ErrorException try_complete_handoff!(
        uncertain_completion_handoff, uncertain_completion_reference)
    HANDOFF_TEST_THROW_COMPLETION[] = false
    @test handoff_slot_status(
        uncertain_completion_handoff, uncertain_completion_reference) ==
        HandoffTransferUncertain
    @test handoff_slot_failure_reason(
        uncertain_completion_handoff, uncertain_completion_reference) ==
        :completion_exception
    @test reclaim_handoff!(
        uncertain_completion_handoff, uncertain_completion_reference) ==
        HandoffSlotNotTerminal

    invalid_sources, invalid_destinations = handoff_test_buffers(1)
    invalid_submission_handoff = prepare_cross_domain_handoff(
        partitions, contract, invalid_sources, invalid_destinations)
    invalid_submission_output =
        handoff_test_slot_output(invalid_submission_handoff)
    @test try_next_free_handoff_slot!(
        invalid_submission_output, invalid_submission_handoff) ==
        HandoffTransitionSucceeded
    invalid_submission_reference = invalid_submission_output[]
    HANDOFF_TEST_INVALID_SUBMISSION[] = true
    @test_throws InvalidConfiguration submit_handoff!(
        invalid_submission_handoff,
        invalid_submission_reference,
        publication,
    )
    HANDOFF_TEST_INVALID_SUBMISSION[] = false
    @test handoff_slot_status(
        invalid_submission_handoff, invalid_submission_reference) ==
        HandoffTransferUncertain
    @test handoff_slot_failure_reason(
        invalid_submission_handoff, invalid_submission_reference) ==
        :invalid_submission_status
    @test try_borrow_completed_handoff!(
        failed_payload_output,
        failed_publication_output,
        invalid_submission_handoff,
        invalid_submission_reference,
    ) == HandoffAlreadyUncertain
    @test reclaim_handoff!(
        invalid_submission_handoff, invalid_submission_reference) ==
        HandoffSlotNotTerminal

    invalid_completion_sources, invalid_completion_destinations =
        handoff_test_buffers(1)
    invalid_completion_handoff = prepare_cross_domain_handoff(
        partitions,
        contract,
        invalid_completion_sources,
        invalid_completion_destinations,
    )
    invalid_completion_output =
        handoff_test_slot_output(invalid_completion_handoff)
    @test try_next_free_handoff_slot!(
        invalid_completion_output, invalid_completion_handoff) ==
        HandoffTransitionSucceeded
    invalid_completion_reference = invalid_completion_output[]
    @test submit_handoff!(
        invalid_completion_handoff,
        invalid_completion_reference,
        publication,
    ) == HandoffTransitionSucceeded
    HANDOFF_TEST_INVALID_COMPLETION[] = true
    @test_throws InvalidConfiguration try_complete_handoff!(
        invalid_completion_handoff, invalid_completion_reference)
    HANDOFF_TEST_INVALID_COMPLETION[] = false
    @test handoff_slot_status(
        invalid_completion_handoff, invalid_completion_reference) ==
        HandoffTransferUncertain
    @test handoff_slot_failure_reason(
        invalid_completion_handoff, invalid_completion_reference) ==
        :invalid_completion_status
    @test try_borrow_completed_handoff!(
        failed_payload_output,
        failed_publication_output,
        invalid_completion_handoff,
        invalid_completion_reference,
    ) == HandoffAlreadyUncertain
    @test reclaim_handoff!(
        invalid_completion_handoff, invalid_completion_reference) ==
        HandoffSlotNotTerminal

    generation_output = handoff_test_slot_output(handoff)
    @test try_next_free_handoff_slot!(generation_output, handoff) ==
        HandoffTransitionSucceeded
    generation_reference = generation_output[]
    generation_state = getfield(handoff, :states)[Int(generation_reference.slot)]
    generation_state.generation = typemax(UInt64)
    generation_state.status = HandoffTransferCompleted
    terminal_reference = HandoffSlotReference(
        handoff_preparation_identity(handoff),
        generation_reference.slot,
        typemax(UInt64),
    )
    @test reclaim_handoff!(handoff, terminal_reference) ==
        HandoffGenerationExhausted
    @test handoff_slot_status(handoff, terminal_reference) ==
        HandoffTransferCompleted
end

@testset "Bidirectional fake-device handoff" begin
    reset_handoff_test_transfer_controls!()
    partitions = handoff_test_partitions()
    identity = HandoffTestContractIdentity()
    contract = HandoffTestContract(
        identity, (Base.OneTo(2), Base.OneTo(3)))
    source = HandoffTestArray(reshape(collect(11.0:16.0), 2, 3))
    destination = zeros(Float64, 2, 3)
    handoff = prepare_cross_domain_handoff(
        partitions, contract, (source,), (destination,))
    reference_output = handoff_test_slot_output(handoff)
    @test try_next_free_handoff_slot!(reference_output, handoff) ==
        HandoffTransitionSucceeded
    reference = reference_output[]
    publication = HandoffTestPublication(identity, UInt64(3))
    @test submit_handoff!(handoff, reference, publication) ==
        HandoffTransitionSucceeded
    @test try_complete_handoff!(handoff, reference) ==
        HandoffTransitionSucceeded
    payload_output = Ref{typeof(destination)}()
    publication_output = Ref{typeof(publication)}()
    @test try_borrow_completed_handoff!(
        payload_output, publication_output, handoff, reference) ==
        HandoffTransitionSucceeded
    @test payload_output[] === destination
    @test destination == source.storage
    @test reclaim_handoff!(handoff, reference) ==
        HandoffTransitionSucceeded
end

@testset "Handoff preparation rejects invalid storage" begin
    partitions = handoff_test_partitions()
    identity = HandoffTestContractIdentity()
    contract = HandoffTestContract(
        identity, (Base.OneTo(2), Base.OneTo(3)))
    sources, destinations = handoff_test_buffers()

    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions, contract, (), (destinations[1],)),
        :invalid_capacity,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions, contract, sources, (destinations[1],)),
        :capacity_mismatch,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions, contract, sources, sources),
        :same_target,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (zeros(Float32, 2, 3),),
            (HandoffTestArray(zeros(Float32, 2, 3)),),
        ),
        :payload_eltype,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (zeros(Float64, 3, 2),),
            (HandoffTestArray(zeros(Float64, 3, 2)),),
        ),
        :payload_axes,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (sources[1], sources[1]),
            destinations,
        ),
        :aliased_source_slots,
    )
    shared = zeros(Float64, 2, 3)
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (shared,),
            (HandoffTestArray(shared),),
        ),
        :aliased_transfer_slots,
    )
    mixed_sources = (
        HandoffTestArray(zeros(Float64, 2, 3), HANDOFF_TEST_ACCELERATOR),
        HandoffTestArray(
            zeros(Float64, 2, 3), HANDOFF_TEST_SECOND_ACCELERATOR),
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            mixed_sources,
            (zeros(Float64, 2, 3), zeros(Float64, 2, 3)),
        ),
        :wrong_target,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (HandoffTestArray(zeros(Float64, 2, 3),
                HANDOFF_TEST_ACCELERATOR),),
            (HandoffTestArray(zeros(Float64, 2, 3),
                HANDOFF_TEST_SECOND_ACCELERATOR),),
        ),
        :unsupported_target_pair,
    )
    assert_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions, :not_a_contract, sources, destinations),
        :invalid_contract,
    )
    unsupported_destination = HandoffUnsupportedTestArray(
        zeros(Float64, 2, 3), HANDOFF_TEST_ACCELERATOR)
    unsupported_error = captured_handoff_preparation_error(
        () -> prepare_cross_domain_handoff(
            partitions,
            contract,
            (zeros(Float64, 2, 3),),
            (unsupported_destination,),
        ))
    @test unsupported_error isa Backends.ComputeDeviceError
    if unsupported_error isa Backends.ComputeDeviceError
        @test unsupported_error.operation === :prepare_transfer
        @test unsupported_error.reason === :unsupported_transfer_pair
        @test unsupported_error.device == HANDOFF_TEST_ACCELERATOR
    end
end

@testset "Warmed handoff allocation contract" begin
    reset_handoff_test_transfer_controls!()
    partitions = handoff_test_partitions()
    identity = HandoffTestContractIdentity()
    contract = HandoffTestContract(
        identity, (Base.OneTo(2), Base.OneTo(3)))
    sources, destinations = handoff_test_buffers(1)
    handoff = prepare_cross_domain_handoff(
        partitions, contract, sources, destinations)
    handoff_allocation_sample!(handoff, contract)
    allocations = @inferred handoff_allocation_sample!(handoff, contract)
    if !coverage_instrumented()
        @test allocations == (
            reference=0,
            source=0,
            submit=0,
            complete=0,
            borrow=0,
            reclaim=0,
        )
    end
end
