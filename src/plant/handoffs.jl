#
# Typed, fixed-capacity cross-device handoffs
#
# Core owns exact-device transfer semantics and the per-slot transfer state.
# It deliberately does not own a queue, producer lease, wait strategy,
# backpressure policy, task, or Agent lifecycle. One external owner must
# serialize every lifecycle operation. Before submission that owner grants
# exclusive mutable access to a source slot; after successful submission the
# handoff owns that slot until terminal reclamation.
#

"""
Typed semantic contract for one prepared handoff payload.

The publication type is immutable and stored inline. Callers must not mutate
any semantic state referenced by a publication between successful submission
and reclamation.
"""
abstract type AbstractHandoffPayloadContract{M} end

function handoff_payload_eltype(contract::AbstractHandoffPayloadContract)
    throw(PlantPreparationError(
        :handoff_contract,
        :unsupported_payload_eltype,
        "handoff contract $(typeof(contract)) must define " *
        "Plant.handoff_payload_eltype",
    ))
end

function handoff_payload_axes(contract::AbstractHandoffPayloadContract)
    throw(PlantPreparationError(
        :handoff_contract,
        :unsupported_payload_axes,
        "handoff contract $(typeof(contract)) must define " *
        "Plant.handoff_payload_axes",
    ))
end

function validate_handoff_publication(
    contract::AbstractHandoffPayloadContract,
    publication,
)
    throw(PlantPreparationError(
        :handoff_publication,
        :unsupported_publication,
        "handoff contract $(typeof(contract)) must define " *
        "Plant.validate_handoff_publication for $(typeof(publication))",
    ))
end

mutable struct _HandoffPreparationIdentityToken end
const _HANDOFF_PREPARATION_IDENTITY_TOKEN =
    _HandoffPreparationIdentityToken()

"""
Opaque run-local identity of one prepared cross-device handoff.

The identity is bound to the exact `PreparedPlantPartitions` object used at
preparation.  It is not a stable digest, placement-plan identity, payload
lease, or transport identifier.
"""
mutable struct HandoffPreparationIdentity{P<:PreparedPlantPartitions}
    const partitions::P

    function HandoffPreparationIdentity(
        token::_HandoffPreparationIdentityToken,
        partitions::P,
    ) where {P<:PreparedPlantPartitions}
        token === _HANDOFF_PREPARATION_IDENTITY_TOKEN || throw(
            ArgumentError("invalid internal handoff-preparation token"))
        return new{P}(partitions)
    end
end

"""Current lifecycle state of one prepared handoff slot."""
@enum HandoffSlotStatus::UInt8 begin
    HandoffSlotFree = 0x01
    HandoffTransferSubmitted = 0x02
    HandoffTransferCompleted = 0x03
    HandoffTransferFailed = 0x04
    HandoffTransferUncertain = 0x05
end

"""Nonthrowing result of one expected handoff ownership transition."""
@enum HandoffStatus::UInt8 begin
    HandoffTransitionSucceeded = 0x01
    HandoffTransferPending = 0x02
    HandoffSlotsExhausted = 0x03
    ForeignHandoffPreparation = 0x04
    InvalidHandoffSlot = 0x05
    InvalidHandoffGeneration = 0x06
    StaleHandoffSlot = 0x07
    HandoffSlotNotFree = 0x08
    HandoffPublicationRejected = 0x09
    HandoffSubmissionFailed = 0x0a
    HandoffNotSubmitted = 0x0b
    DuplicateHandoffCompletion = 0x0c
    HandoffCompletionFailed = 0x0d
    HandoffAlreadyFailed = 0x0e
    HandoffPayloadNotVisible = 0x0f
    HandoffSlotNotTerminal = 0x10
    HandoffGenerationExhausted = 0x11
    HandoffAlreadyUncertain = 0x12
end

"""Generation-bearing descriptor for one slot in one prepared handoff."""
struct HandoffSlotReference{I<:HandoffPreparationIdentity}
    identity::I
    slot::UInt32
    generation::UInt64
end

mutable struct _HandoffSlotState
    generation::UInt64
    status::HandoffSlotStatus
    failure_reason::Symbol
end

"""
One typed, fixed-capacity, caller-driven transfer between distinct exact
compute devices.

The source and destination arrays remain caller-provided storage, but the
prepared handoff has exclusive ownership of every destination slot.  Callers
must access destination storage only through a successful completed-slot
borrow.  Slot state and publication metadata are host-resident and concrete.
One external owner must serialize slot selection, submission, completion
observation, borrowing, and reclamation. The object is not thread-safe. HIL
supplies descriptor/lease handoff around that owner rather than invoking these
mutations concurrently from producer and consumer Agents.
"""
struct PreparedCrossDomainHandoff{
    M,
    C<:AbstractHandoffPayloadContract{M},
    I<:HandoffPreparationIdentity,
    ST<:AbstractComputeDevice,
    DT<:AbstractComputeDevice,
    S<:AbstractArray,
    D<:AbstractArray,
    SC,
    DC,
    B,
}
    identity::I
    contract::C
    source_target::ST
    destination_target::DT
    source_context::SC
    destination_context::DC
    source_slots::Memory{S}
    destination_slots::Memory{D}
    backend_states::Memory{B}
    publications::Memory{M}
    states::Memory{_HandoffSlotState}
    logical_payload_bytes::UInt64
end

@inline handoff_preparation_identity(
    handoff::PreparedCrossDomainHandoff) = handoff.identity
@inline handoff_contract(handoff::PreparedCrossDomainHandoff) =
    handoff.contract
@inline handoff_source_target(handoff::PreparedCrossDomainHandoff) =
    handoff.source_target
@inline handoff_destination_target(handoff::PreparedCrossDomainHandoff) =
    handoff.destination_target
@inline handoff_capacity(handoff::PreparedCrossDomainHandoff) =
    length(handoff.states)
@inline handoff_payload_bytes(handoff::PreparedCrossDomainHandoff) =
    handoff.logical_payload_bytes

@noinline function _handoff_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(:handoff, reason, String(message)))
end

function _copy_handoff_slots(slots, label::Symbol)
    count = length(slots)
    count > 0 || _handoff_preparation_error(
        :invalid_capacity,
        "$label handoff slots must not be empty",
    )
    count <= typemax(UInt32) || _handoff_preparation_error(
        :capacity,
        "$label handoff-slot count exceeds UInt32 capacity",
    )
    iteration = iterate(slots)
    isnothing(iteration) && _handoff_preparation_error(
        :invalid_capacity,
        "$label handoff slots reported positive length but iterate as empty",
    )
    first_slot, iteration_state = iteration
    first_slot = _require_handoff_array_slot(first_slot, label)
    S = typeof(first_slot)
    copied = Memory{S}(undef, count)
    copied[1] = first_slot
    ordinal = 2
    while true
        iteration = iterate(slots, iteration_state)
        isnothing(iteration) && break
        slot, iteration_state = iteration
        ordinal <= count || _handoff_preparation_error(
            :invalid_capacity,
            "$label handoff slot iterator yielded more entries than length",
        )
        typeof(slot) === S || _handoff_preparation_error(
            :heterogeneous_storage,
            "$label handoff slots must have one concrete array type; got " *
            "$S and $(typeof(slot))",
        )
        @inbounds copied[ordinal] = slot
        ordinal += 1
    end
    ordinal == count + 1 || _handoff_preparation_error(
        :invalid_capacity,
        "$label handoff slot iterator yielded fewer entries than length",
    )
    return copied
end

@inline _require_handoff_array_slot(slot::AbstractArray, ::Symbol) = slot

function _require_handoff_array_slot(slot, label::Symbol)
    _handoff_preparation_error(
        :invalid_storage,
        "$label handoff slots must contain arrays; got $(typeof(slot))",
    )
end

@inline _require_handoff_payload_type(::Type{T}) where {T} = T

function _require_handoff_payload_type(value)
    _handoff_preparation_error(
        :invalid_payload_eltype,
        "handoff_payload_eltype must return a type; got $(typeof(value))",
    )
end

@inline _require_handoff_payload_axis(::AbstractUnitRange) = nothing

function _require_handoff_payload_axis(axis)
    _handoff_preparation_error(
        :invalid_payload_axes,
        "every handoff payload axis must be an AbstractUnitRange; got " *
        "$(typeof(axis))",
    )
end

function _require_handoff_payload_axes(expected_axes::Tuple)
    foreach(_require_handoff_payload_axis, expected_axes)
    return expected_axes
end

function _require_handoff_payload_axes(value)
    _handoff_preparation_error(
        :invalid_payload_axes,
        "handoff_payload_axes must return a tuple; got $(typeof(value))",
    )
end

function _require_handoff_contract(
    contract::C,
    ::Type{M},
) where {M,C<:AbstractHandoffPayloadContract{M}}
    isconcretetype(C) || _handoff_preparation_error(
        :abstract_contract,
        "handoff payload contract must be concrete",
    )
    ismutabletype(C) && _handoff_preparation_error(
        :mutable_contract,
        "handoff payload contract must be immutable",
    )
    isconcretetype(M) || _handoff_preparation_error(
        :abstract_publication,
        "handoff publication metadata type must be concrete",
    )
    ismutabletype(M) && _handoff_preparation_error(
        :mutable_publication,
        "handoff publication metadata type must be immutable",
    )
    Base.allocatedinline(M) || _handoff_preparation_error(
        :boxed_publication,
        "handoff publication metadata must use Julia's inline array representation",
    )
    T = _require_handoff_payload_type(handoff_payload_eltype(contract))
    isconcretetype(T) || _handoff_preparation_error(
        :abstract_payload_eltype,
        "handoff payload element type must be concrete",
    )
    isbitstype(T) || _handoff_preparation_error(
        :nonisbits_payload_eltype,
        "handoff payload element type must be isbits; got $T",
    )
    expected_axes = _require_handoff_payload_axes(
        handoff_payload_axes(contract))
    return T, expected_axes
end

function _require_handoff_slot_storage(
    slots::Memory,
    target::AbstractComputeDevice,
    expected_eltype::Type,
    expected_axes::Tuple,
    label::Symbol,
)
    @inbounds for ordinal in eachindex(slots)
        slot = slots[ordinal]
        eltype(slot) === expected_eltype || _handoff_preparation_error(
            :payload_eltype,
            "$label handoff slot $ordinal has element type $(eltype(slot)); " *
            "expected $expected_eltype",
        )
        axes(slot) == expected_axes || _handoff_preparation_error(
            :payload_axes,
            "$label handoff slot $ordinal has axes $(axes(slot)); expected " *
            "$expected_axes",
        )
        compute_device(slot) == target || _handoff_preparation_error(
            :wrong_target,
            "$label handoff slot $ordinal does not occupy exact target $target",
        )
    end
    return nothing
end

function _require_nonaliasing_handoff_slots(
    source_slots::Memory,
    destination_slots::Memory,
)
    @inbounds for left in eachindex(source_slots)
        for right in (left + 1):lastindex(source_slots)
            Base.mightalias(source_slots[left], source_slots[right]) &&
                _handoff_preparation_error(
                    :aliased_source_slots,
                    "source handoff slots $left and $right may alias",
                )
        end
    end
    @inbounds for left in eachindex(destination_slots)
        for right in (left + 1):lastindex(destination_slots)
            Base.mightalias(destination_slots[left], destination_slots[right]) &&
                _handoff_preparation_error(
                    :aliased_destination_slots,
                    "destination handoff slots $left and $right may alias",
                )
        end
    end
    @inbounds for source in eachindex(source_slots)
        for destination in eachindex(destination_slots)
            Base.mightalias(source_slots[source],
                destination_slots[destination]) &&
                _handoff_preparation_error(
                    :aliased_transfer_slots,
                    "source handoff slot $source and destination slot " *
                    "$destination may alias",
                )
        end
    end
    return nothing
end

function _checked_handoff_payload_bytes(
    ::Type{T},
    expected_axes::Tuple,
) where {T}
    count = UInt128(1)
    for axis in expected_axes
        count *= UInt128(length(axis))
        count <= typemax(UInt64) || _handoff_preparation_error(
            :payload_size_overflow,
            "handoff payload element count exceeds UInt64 capacity",
        )
    end
    bytes = count * UInt128(sizeof(T))
    bytes <= typemax(UInt64) || _handoff_preparation_error(
        :payload_size_overflow,
        "handoff logical payload byte count exceeds UInt64 capacity",
    )
    return UInt64(bytes)
end

function _prepare_handoff_backend_states(
    source_slots::Memory,
    destination_slots::Memory,
    source_context,
    destination_context,
)
    first_state = _prepare_array_transfer(
        source_slots[1],
        destination_slots[1],
        source_context,
        destination_context,
    )
    B = typeof(first_state)
    isconcretetype(B) || _handoff_preparation_error(
        :abstract_backend_state,
        "prepared handoff backend state must be concrete",
    )
    states = Memory{B}(undef, length(source_slots))
    states[1] = first_state
    @inbounds for ordinal in 2:length(source_slots)
        state = _prepare_array_transfer(
            source_slots[ordinal],
            destination_slots[ordinal],
            source_context,
            destination_context,
        )
        typeof(state) === B || _handoff_preparation_error(
            :heterogeneous_backend_state,
            "all prepared handoff slots must use one concrete backend-state type",
        )
        states[ordinal] = state
    end
    return states
end

@inline _require_supported_handoff_target_pair(
    ::HostComputeDevice,
    ::AcceleratorComputeDevice,
) = nothing

@inline _require_supported_handoff_target_pair(
    ::AcceleratorComputeDevice,
    ::HostComputeDevice,
) = nothing

function _require_supported_handoff_target_pair(
    source::AbstractComputeDevice,
    destination::AbstractComputeDevice,
)
    _handoff_preparation_error(
        :unsupported_target_pair,
        "Gate 9A cross-domain handoffs support host-to-accelerator or " *
        "accelerator-to-host transfers; got $source to $destination",
    )
end

"""
    prepare_cross_domain_handoff(partitions, contract, source_slots,
        destination_slots)

Prepare fixed paired payload slots for one explicit transfer between distinct
exact compute devices.  The caller retains the arrays and must provide
exclusive source-slot ownership before submission.  Preparation transfers
exclusive destination-slot ownership to the handoff; callers may read a
destination only after a successful completed-slot borrow.  Backend support
is explicitly dispatched for the concrete source/destination storage pair;
there is no generic cross-backend copy fallback.
"""
function prepare_cross_domain_handoff(
    partitions::PreparedPlantPartitions,
    contract::C,
    source_slot_inputs,
    destination_slot_inputs,
) where {M,C<:AbstractHandoffPayloadContract{M}}
    source_slots = _copy_handoff_slots(source_slot_inputs, :source)
    destination_slots = _copy_handoff_slots(
        destination_slot_inputs, :destination)
    length(source_slots) == length(destination_slots) ||
        _handoff_preparation_error(
            :capacity_mismatch,
            "source and destination handoff capacities must be equal",
        )
    expected_eltype, expected_axes = _require_handoff_contract(contract, M)
    source_target = compute_device(source_slots[1])
    destination_target = compute_device(destination_slots[1])
    source_target != destination_target || _handoff_preparation_error(
        :same_target,
        "cross-domain handoff source and destination must occupy distinct " *
        "exact compute devices",
    )
    _require_supported_handoff_target_pair(
        source_target, destination_target)
    _require_handoff_slot_storage(
        source_slots,
        source_target,
        expected_eltype,
        expected_axes,
        :source,
    )
    _require_handoff_slot_storage(
        destination_slots,
        destination_target,
        expected_eltype,
        expected_axes,
        :destination,
    )
    _require_nonaliasing_handoff_slots(source_slots, destination_slots)
    source_context = _prepare_device_execution_context(source_target)
    destination_context = _prepare_device_execution_context(destination_target)
    _prepared_device_execution_compute_device(source_context) == source_target ||
        _handoff_preparation_error(
            :source_context_target,
            "prepared source context belongs to another exact compute device",
        )
    _prepared_device_execution_compute_device(destination_context) ==
        destination_target || _handoff_preparation_error(
            :destination_context_target,
            "prepared destination context belongs to another exact compute device",
        )
    backend_states = _prepare_handoff_backend_states(
        source_slots,
        destination_slots,
        source_context,
        destination_context,
    )
    identity = HandoffPreparationIdentity(
        _HANDOFF_PREPARATION_IDENTITY_TOKEN, partitions)
    publications = Memory{M}(undef, length(source_slots))
    states = Memory{_HandoffSlotState}(undef, length(source_slots))
    @inbounds for ordinal in eachindex(states)
        states[ordinal] = _HandoffSlotState(
            one(UInt64), HandoffSlotFree, :none)
    end
    return PreparedCrossDomainHandoff(
        identity,
        contract,
        source_target,
        destination_target,
        source_context,
        destination_context,
        source_slots,
        destination_slots,
        backend_states,
        publications,
        states,
        _checked_handoff_payload_bytes(expected_eltype, expected_axes),
    )
end

function prepare_cross_domain_handoff(
    ::PreparedPlantPartitions,
    contract,
    source_slots,
    destination_slots,
)
    _handoff_preparation_error(
        :invalid_contract,
        "handoff contract must subtype AbstractHandoffPayloadContract; got " *
        "$(typeof(contract))",
    )
end

@inline function _handoff_reference_status(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    reference.identity === handoff.identity ||
        return ForeignHandoffPreparation, 0
    iszero(reference.slot) && return InvalidHandoffSlot, 0
    UInt64(reference.slot) <= UInt64(length(handoff.states)) ||
        return InvalidHandoffSlot, 0
    iszero(reference.generation) && return InvalidHandoffGeneration, 0
    ordinal = Int(reference.slot)
    state = @inbounds handoff.states[ordinal]
    reference.generation == state.generation ||
        return StaleHandoffSlot, 0
    return HandoffTransitionSucceeded, ordinal
end

@noinline function _throw_handoff_reference_error(status::HandoffStatus)
    _handoff_preparation_error(
        :invalid_slot_reference,
        "handoff slot reference is not usable: $status",
    )
end

"""
Write the first currently free generation-bearing slot reference to `output`.

This does not claim the slot. The externally serialized owner must ensure one
producer has exclusive pre-submit access until `submit_handoff!` succeeds.
"""
function try_next_free_handoff_slot!(
    output::Base.RefValue{HandoffSlotReference{I}},
    handoff::PreparedCrossDomainHandoff{M,C,I},
) where {M,C,I}
    @inbounds for ordinal in eachindex(handoff.states)
        state = handoff.states[ordinal]
        state.status == HandoffSlotFree || continue
        output[] = HandoffSlotReference(
            handoff.identity,
            UInt32(ordinal),
            state.generation,
        )
        return HandoffTransitionSucceeded
    end
    return HandoffSlotsExhausted
end

"""Return caller-writable source storage for one current free slot."""
function producer_handoff_payload(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded ||
        _throw_handoff_reference_error(status)
    state = @inbounds handoff.states[ordinal]
    state.status == HandoffSlotFree || _handoff_preparation_error(
        :wrong_slot_owner,
        "producer access requires a free handoff slot",
    )
    return @inbounds handoff.source_slots[ordinal]
end

@inline function _validate_handoff_publication(
    contract::AbstractHandoffPayloadContract,
    publication,
)
    try
        result = validate_handoff_publication(contract, publication)
        isnothing(result) || return false
        return true
    catch error
        return _handoff_publication_validation_failure(error)
    end
end

@inline _handoff_publication_validation_failure(
    ::PlantPreparationError) = false

@noinline function _handoff_publication_validation_failure(error)
    throw(error)
end

"""Launch one explicit backend transfer without asserting completion."""
function submit_handoff!(
    handoff::PreparedCrossDomainHandoff{M},
    reference::HandoffSlotReference,
    publication::M,
) where {M}
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded || return status
    state = @inbounds handoff.states[ordinal]
    state.status == HandoffSlotFree || return HandoffSlotNotFree
    _validate_handoff_publication(handoff.contract, publication) ||
        return HandoffPublicationRejected
    @inbounds handoff.publications[ordinal] = publication
    state.status = HandoffTransferUncertain
    state.failure_reason = :submission_in_progress
    submission = try
        @inbounds _submit_prepared_array_transfer!(
            handoff.backend_states[ordinal],
            handoff.destination_slots[ordinal],
            handoff.source_slots[ordinal],
        )
    catch
        state.failure_reason = :submission_exception
        rethrow()
    end
    if submission == _PreparedArrayTransferSubmitted
        state.status = HandoffTransferSubmitted
        state.failure_reason = :none
        return HandoffTransitionSucceeded
    end
    if submission == _PreparedArrayTransferSubmissionFailed
        state.status = HandoffTransferFailed
        state.failure_reason = :submission_failed
        return HandoffSubmissionFailed
    end
    state.failure_reason = :invalid_submission_status
    throw(InvalidConfiguration(
        "prepared array-transfer submission returned unsupported status " *
        "$(repr(submission)); slot ownership is now uncertain"))
end

function submit_handoff!(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
    publication,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded || return status
    state = @inbounds handoff.states[ordinal]
    state.status == HandoffSlotFree || return HandoffSlotNotFree
    return HandoffPublicationRejected
end

"""
Observe one submitted transfer once without polling or waiting in Core.

A pending observation leaves the slot submitted.  Only a completed observation
makes the destination payload visible to the consumer.
"""
function try_complete_handoff!(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded || return status
    state = @inbounds handoff.states[ordinal]
    state.status == HandoffTransferCompleted &&
        return DuplicateHandoffCompletion
    state.status == HandoffTransferFailed && return HandoffAlreadyFailed
    state.status == HandoffTransferUncertain && return HandoffAlreadyUncertain
    state.status == HandoffTransferSubmitted || return HandoffNotSubmitted
    state.status = HandoffTransferUncertain
    state.failure_reason = :completion_observation_in_progress
    completion = try
        @inbounds _observe_prepared_array_transfer_completion!(
            handoff.backend_states[ordinal])
    catch
        state.failure_reason = :completion_exception
        rethrow()
    end
    if completion == _PreparedArrayTransferPending
        state.status = HandoffTransferSubmitted
        state.failure_reason = :none
        return HandoffTransferPending
    end
    if completion == _PreparedArrayTransferCompleted
        state.status = HandoffTransferCompleted
        state.failure_reason = :none
        return HandoffTransitionSucceeded
    end
    if completion == _PreparedArrayTransferCompletionFailed
        state.status = HandoffTransferFailed
        state.failure_reason = :completion_failed
        return HandoffCompletionFailed
    end
    state.failure_reason = :invalid_completion_status
    throw(InvalidConfiguration(
        "prepared array-transfer completion observation returned unsupported " *
        "status $(repr(completion)); slot ownership is now uncertain"))
end

"""
Borrow the completed destination and its typed publication metadata.

The handoff retains slot ownership until `reclaim_handoff!`; the external
single-consumer contract must prevent concurrent use and reclamation.
"""
function try_borrow_completed_handoff!(
    payload_output::Base.RefValue{D},
    publication_output::Base.RefValue{M},
    handoff::PreparedCrossDomainHandoff{M,C,I,ST,DT,S,D},
    reference::HandoffSlotReference,
) where {M,C,I,ST,DT,S,D}
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded || return status
    state = @inbounds handoff.states[ordinal]
    state.status == HandoffTransferFailed && return HandoffAlreadyFailed
    state.status == HandoffTransferUncertain && return HandoffAlreadyUncertain
    state.status == HandoffTransferCompleted ||
        return HandoffPayloadNotVisible
    @inbounds begin
        payload_output[] = handoff.destination_slots[ordinal]
        publication_output[] = handoff.publications[ordinal]
    end
    return HandoffTransitionSucceeded
end

"""Reclaim one completed or failed slot and advance its generation."""
function reclaim_handoff!(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded || return status
    state = @inbounds handoff.states[ordinal]
    terminal = state.status == HandoffTransferCompleted ||
        state.status == HandoffTransferFailed
    terminal || return HandoffSlotNotTerminal
    state.generation == typemax(UInt64) &&
        return HandoffGenerationExhausted
    state.generation += one(UInt64)
    state.status = HandoffSlotFree
    state.failure_reason = :none
    return HandoffTransitionSucceeded
end

function handoff_slot_status(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded ||
        _throw_handoff_reference_error(status)
    return @inbounds handoff.states[ordinal].status
end

function handoff_slot_failure_reason(
    handoff::PreparedCrossDomainHandoff,
    reference::HandoffSlotReference,
)
    status, ordinal = _handoff_reference_status(handoff, reference)
    status == HandoffTransitionSucceeded ||
        _throw_handoff_reference_error(status)
    return @inbounds handoff.states[ordinal].failure_reason
end
