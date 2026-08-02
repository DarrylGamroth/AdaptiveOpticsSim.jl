struct HandoffTestBackend <: AbstractArrayBackend end

const HANDOFF_TEST_ACCELERATOR =
    AcceleratorComputeDevice(HandoffTestBackend(), UInt32(1))
const HANDOFF_TEST_SECOND_ACCELERATOR =
    AcceleratorComputeDevice(HandoffTestBackend(), UInt32(2))

struct HandoffTestArray{T,N,A<:Array{T,N},D} <: DenseArray{T,N}
    storage::A
    device::D
end

HandoffTestArray(storage::Array) =
    HandoffTestArray(storage, HANDOFF_TEST_ACCELERATOR)

Base.size(array::HandoffTestArray) = size(array.storage)
Base.axes(array::HandoffTestArray) = axes(array.storage)
Base.IndexStyle(::Type{<:HandoffTestArray}) = IndexLinear()
Base.getindex(::HandoffTestArray, ::Int) =
    error("scalar indexing is disabled for HandoffTestArray")
Base.setindex!(::HandoffTestArray, value, ::Int) =
    error("scalar indexing is disabled for HandoffTestArray")
Base.similar(array::HandoffTestArray, ::Type{T}, dimensions::Dims{N}) where {
    T,N} = HandoffTestArray(Array{T}(undef, dimensions), array.device)
Base.similar(array::HandoffTestArray, dimensions::Dims{N}) where {N} =
    similar(array, eltype(array), dimensions)
Base.mightalias(left::HandoffTestArray, right::HandoffTestArray) =
    Base.mightalias(left.storage, right.storage)
Base.mightalias(left::HandoffTestArray, right::Array) =
    Base.mightalias(left.storage, right)
Base.mightalias(left::Array, right::HandoffTestArray) =
    Base.mightalias(left, right.storage)

Backends.array_backend_selector(::Type{<:HandoffTestArray}) =
    HandoffTestBackend()
Backends.backend(::HandoffTestArray) = HandoffTestBackend()
Backends.compute_device(array::HandoffTestArray) = array.device
Backends.execution_style(::HandoffTestArray) = ScalarCPUStyle()
Backends.compute_device_availability(
    ::AcceleratorComputeDevice{HandoffTestBackend}) = ComputeDeviceAvailable()

struct HandoffTestContext{D} <:
    Backends._AbstractPreparedDeviceExecutionContext
    target::D
end

Backends._prepare_device_execution_context(
    target::AcceleratorComputeDevice{HandoffTestBackend}) =
    HandoffTestContext(target)
Backends._prepared_device_execution_compute_device(
    context::HandoffTestContext) = context.target

const HANDOFF_TEST_ACTIVE_DEVICE = Ref{UInt32}(0)
const HANDOFF_TEST_CONTEXT_ENTRIES = Ref{UInt64}(0)

function Backends._with_prepared_device_execution_context(
    f::F,
    context::HandoffTestContext,
) where {F}
    old_device = HANDOFF_TEST_ACTIVE_DEVICE[]
    HANDOFF_TEST_ACTIVE_DEVICE[] =
        UInt32(compute_device_identifier(context.target))
    HANDOFF_TEST_CONTEXT_ENTRIES[] += one(UInt64)
    try
        return f()
    finally
        HANDOFF_TEST_ACTIVE_DEVICE[] = old_device
    end
end

Backends._synchronize_prepared_device_execution_context!(
    ::HandoffTestContext) = nothing

mutable struct HandoffTestHostToDeviceTransfer{C}
    context::C
    remaining_pending::UInt8
    submitted::Bool
end

mutable struct HandoffTestDeviceToHostTransfer{C}
    context::C
    remaining_pending::UInt8
    submitted::Bool
end

const HandoffTestPreparedTransfer = Union{
    HandoffTestHostToDeviceTransfer,
    HandoffTestDeviceToHostTransfer,
}

const HANDOFF_TEST_PENDING_OBSERVATIONS = Ref{UInt8}(0)
const HANDOFF_TEST_FAIL_SUBMISSION = Ref(false)
const HANDOFF_TEST_FAIL_COMPLETION = Ref(false)
const HANDOFF_TEST_THROW_SUBMISSION = Ref(false)
const HANDOFF_TEST_THROW_COMPLETION = Ref(false)
const HANDOFF_TEST_INVALID_SUBMISSION = Ref(false)
const HANDOFF_TEST_INVALID_COMPLETION = Ref(false)

function reset_handoff_test_transfer_controls!()
    HANDOFF_TEST_PENDING_OBSERVATIONS[] = zero(UInt8)
    HANDOFF_TEST_FAIL_SUBMISSION[] = false
    HANDOFF_TEST_FAIL_COMPLETION[] = false
    HANDOFF_TEST_THROW_SUBMISSION[] = false
    HANDOFF_TEST_THROW_COMPLETION[] = false
    HANDOFF_TEST_INVALID_SUBMISSION[] = false
    HANDOFF_TEST_INVALID_COMPLETION[] = false
    HANDOFF_TEST_ACTIVE_DEVICE[] = zero(UInt32)
    HANDOFF_TEST_CONTEXT_ENTRIES[] = zero(UInt64)
    return nothing
end

function Backends._prepare_array_transfer(
    ::Array,
    ::HandoffTestArray,
    ::Backends._HostPreparedDeviceExecutionContext,
    destination_context::HandoffTestContext,
)
    return HandoffTestHostToDeviceTransfer(
        destination_context, zero(UInt8), false)
end

function Backends._prepare_array_transfer(
    ::HandoffTestArray,
    ::Array,
    source_context::HandoffTestContext,
    ::Backends._HostPreparedDeviceExecutionContext,
)
    return HandoffTestDeviceToHostTransfer(
        source_context, zero(UInt8), false)
end

function Backends._submit_prepared_array_transfer!(
    transfer::HandoffTestHostToDeviceTransfer,
    destination::HandoffTestArray,
    source::Array,
)
    return Backends._with_prepared_device_execution_context(
        transfer.context) do
        HANDOFF_TEST_ACTIVE_DEVICE[] ==
            UInt32(compute_device_identifier(destination.device)) ||
            error("fake destination context is not active")
        if HANDOFF_TEST_FAIL_SUBMISSION[]
            return Backends._PreparedArrayTransferSubmissionFailed
        end
        copyto!(destination.storage, source)
        HANDOFF_TEST_THROW_SUBMISSION[] &&
            error("injected uncertain handoff submission")
        HANDOFF_TEST_INVALID_SUBMISSION[] && return :invalid_submission
        transfer.remaining_pending = HANDOFF_TEST_PENDING_OBSERVATIONS[]
        transfer.submitted = true
        return Backends._PreparedArrayTransferSubmitted
    end
end

function Backends._submit_prepared_array_transfer!(
    transfer::HandoffTestDeviceToHostTransfer,
    destination::Array,
    source::HandoffTestArray,
)
    return Backends._with_prepared_device_execution_context(
        transfer.context) do
        HANDOFF_TEST_ACTIVE_DEVICE[] ==
            UInt32(compute_device_identifier(source.device)) ||
            error("fake source context is not active")
        if HANDOFF_TEST_FAIL_SUBMISSION[]
            return Backends._PreparedArrayTransferSubmissionFailed
        end
        copyto!(destination, source.storage)
        HANDOFF_TEST_THROW_SUBMISSION[] &&
            error("injected uncertain handoff submission")
        HANDOFF_TEST_INVALID_SUBMISSION[] && return :invalid_submission
        transfer.remaining_pending = HANDOFF_TEST_PENDING_OBSERVATIONS[]
        transfer.submitted = true
        return Backends._PreparedArrayTransferSubmitted
    end
end

function Backends._observe_prepared_array_transfer_completion!(
    transfer::HandoffTestPreparedTransfer,
)
    return Backends._with_prepared_device_execution_context(
        transfer.context) do
        transfer.submitted || error("fake transfer was not submitted")
        HANDOFF_TEST_FAIL_COMPLETION[] && begin
            transfer.submitted = false
            return Backends._PreparedArrayTransferCompletionFailed
        end
        HANDOFF_TEST_THROW_COMPLETION[] &&
            error("injected uncertain handoff completion")
        HANDOFF_TEST_INVALID_COMPLETION[] && return :invalid_completion
        if !iszero(transfer.remaining_pending)
            transfer.remaining_pending -= one(UInt8)
            return Backends._PreparedArrayTransferPending
        end
        transfer.submitted = false
        return Backends._PreparedArrayTransferCompleted
    end
end

struct HandoffUnsupportedTestArray{T,N,A<:Array{T,N},D} <:
    AbstractArray{T,N}
    storage::A
    device::D
end

Base.size(array::HandoffUnsupportedTestArray) = size(array.storage)
Base.axes(array::HandoffUnsupportedTestArray) = axes(array.storage)
Base.IndexStyle(::Type{<:HandoffUnsupportedTestArray}) = IndexLinear()
Base.getindex(::HandoffUnsupportedTestArray, ::Int) =
    error("scalar indexing is disabled for HandoffUnsupportedTestArray")
Base.setindex!(::HandoffUnsupportedTestArray, value, ::Int) =
    error("scalar indexing is disabled for HandoffUnsupportedTestArray")
Base.mightalias(left::Array, right::HandoffUnsupportedTestArray) =
    Base.mightalias(left, right.storage)
Base.mightalias(left::HandoffUnsupportedTestArray, right::Array) =
    Base.mightalias(left.storage, right)
Backends.array_backend_selector(::Type{<:HandoffUnsupportedTestArray}) =
    HandoffTestBackend()
Backends.backend(::HandoffUnsupportedTestArray) = HandoffTestBackend()
Backends.compute_device(array::HandoffUnsupportedTestArray) = array.device
Backends.execution_style(::HandoffUnsupportedTestArray) = ScalarCPUStyle()

mutable struct HandoffTestContractIdentity end

struct HandoffTestPublication{I}
    identity::I
    sequence::UInt64
end

struct HandoffTestContract{I,A} <:
    AbstractHandoffPayloadContract{HandoffTestPublication{I}}
    identity::I
    axes::A
end

Plant.handoff_payload_eltype(::HandoffTestContract) = Float64
Plant.handoff_payload_axes(contract::HandoffTestContract) = contract.axes

function Plant.validate_handoff_publication(
    contract::HandoffTestContract,
    publication::HandoffTestPublication,
)
    publication.identity === contract.identity || throw(
        PlantPreparationError(
            :handoff_publication,
            :foreign_contract,
            "test publication belongs to another handoff contract",
        ))
    iszero(publication.sequence) && throw(PlantPreparationError(
        :handoff_publication,
        :invalid_sequence,
        "test publication sequence must be positive",
    ))
    return nothing
end

function handoff_test_partitions()
    definition = PlantDefinition(;
        telescope=TelescopeDefinition(
            resolution=4,
            diameter=4.0,
            central_obstruction=0.0,
            revision=1,
        ),
        atmosphere=KolmogorovAtmosphereDefinition(r0=0.2, L0=25.0),
    )
    assignment = resolve_plant_partition_assignment(
        definition, HostComputeDevice())
    return prepare_plant_partitions(
        definition,
        assignment;
        run_seed=0x206,
    )
end

function captured_handoff_preparation_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end
