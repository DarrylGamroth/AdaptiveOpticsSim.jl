const _SNAPSHOT_FREE = UInt8(0)
const _SNAPSHOT_CAPTURING = UInt8(1)
const _SNAPSHOT_PENDING = UInt8(2)
const _SNAPSHOT_TRANSFERRING = UInt8(3)
const _SNAPSHOT_READY = UInt8(4)
const _SNAPSHOT_FAILED = UInt8(5)
const _SNAPSHOT_CLOSED = UInt8(6)

mutable struct _GraphOutputSnapshotState
    phase::Base.Threads.Atomic{UInt8}
    sequence::UInt64
    generation::UInt64
    pts::Int64
    model_time_ns::Int64
end

"""
One optional, capacity-one snapshot of an exact graph output. The graph owner
alone submits; one observation consumer polls, reads, and releases. The owner
retains all transfer storage, including backend-owned host views. Busy consumers
drop observations without delaying a graph step.
"""
struct PreparedGraphOutputSnapshot{G,A,T}
    graph::G
    source::A
    transfer::T
    state::_GraphOutputSnapshotState
end

"""
    prepare_graph_output_snapshot(graph, output::Union{Symbol,Val})

Bind one completed graph output to a prepared host snapshot. Preparation is cold
and requires no pending graph step. CPU and CUDA arrays are supported; other
accelerators are rejected. CUDA prepares two captured array copies: a device
snapshot ordered before the next plant step, and a separate mapped-host copy.
"""
function prepare_graph_output_snapshot(graph::PreparedAlgorithmGraph, output)
    (graph_failed(graph) || graph_step_pending(graph)) && throw(AlgorithmGraphError(
        "output snapshot preparation requires a healthy completed graph"))
    source = graph_output(graph, output)
    (eltype(source) <: AbstractFloat && isbitstype(eltype(source))) ||
        throw(AlgorithmGraphError("output snapshots require isbits floating-point arrays"))
    transfer = Backends._prepare_output_snapshot(source, graph.context)
    state = _GraphOutputSnapshotState(Base.Threads.Atomic{UInt8}(_SNAPSHOT_FREE),
        UInt64(0), UInt64(0), Int64(0), Int64(0))
    return PreparedGraphOutputSnapshot(graph, source, transfer, state)
end

"""
    try_snapshot_graph_output!(snapshot; generation=UInt64(0), pts=Int64(0),
                               model_time_ns=Int64(0)) -> Bool

The graph's single producer makes one slot reservation attempt after a completed
frame. Sequence zero is permitted for pre-admission warmup and remains zero in
the published metadata. A busy, failed, or closed observer returns false. Metadata and
the graph sequence are published together before the consumer can start its
transfer. Submission does not wait for the consumer or GPU completion. The
producer must submit before the next graph step; it must not reset the graph
concurrently with either snapshot participant.
"""
function try_snapshot_graph_output!(snapshot::PreparedGraphOutputSnapshot;
    generation::UInt64=UInt64(0), pts::Int64=Int64(0), model_time_ns::Int64=Int64(0))
    state = snapshot.state
    Base.Threads.atomic_cas!(state.phase, _SNAPSHOT_FREE, _SNAPSHOT_CAPTURING) ==
        _SNAPSHOT_FREE || return false
    try
        graph = snapshot.graph
        (graph_failed(graph) || graph_step_pending(graph)) && throw(AlgorithmGraphError(
            "output snapshot submission requires a healthy completed graph"))
        Backends._submit_output_snapshot!(snapshot.transfer, snapshot.source)
        state.sequence = graph_step_sequence(graph)
        state.generation = generation
        state.pts = pts
        state.model_time_ns = model_time_ns
        state.phase[] = _SNAPSHOT_PENDING
        return true
    catch
        # Ownership is uncertain until a cold drain succeeds. Never reuse the
        # arrays or event generation after an exceptional submission.
        state.phase[] = _SNAPSHOT_FAILED
        rethrow()
    end
end

"""
    poll_graph_output_snapshot!(snapshot) -> Bool

The single consumer starts any pending host transfer and makes one nonblocking
completion query. True grants read ownership of the host values and metadata
until `release_graph_output_snapshot!`. Repeated polls while ready return true.
False means no sample is ready, including a failed or closed observer.
"""
function poll_graph_output_snapshot!(snapshot::PreparedGraphOutputSnapshot)
    state = snapshot.state
    phase = state.phase[]
    phase == _SNAPSHOT_READY && return true
    phase in (_SNAPSHOT_PENDING, _SNAPSHOT_TRANSFERRING) || return false
    try
        if phase == _SNAPSHOT_PENDING
            state.phase[] = _SNAPSHOT_TRANSFERRING
            Backends._start_output_snapshot_transfer!(snapshot.transfer)
        end
        Backends._poll_output_snapshot_transfer!(snapshot.transfer) || return false
        state.phase[] = _SNAPSHOT_READY
        return true
    catch
        state.phase[] = _SNAPSHOT_FAILED
        rethrow()
    end
end

"""
Return retained host storage for cold preparation of consumers. During execution,
read it only after a successful poll and before release. Do not resize, replace,
or write it, and retain the snapshot for as long as this storage is referenced.
"""
@inline snapshot_values(snapshot::PreparedGraphOutputSnapshot) =
    Backends._output_snapshot_values(snapshot.transfer)

"Return the metadata associated with the consumer's currently ready snapshot."
function snapshot_metadata(snapshot::PreparedGraphOutputSnapshot)
    state = snapshot.state
    state.phase[] == _SNAPSHOT_READY || throw(AlgorithmGraphError(
        "snapshot metadata requires consumer read ownership"))
    return (; state.sequence, state.generation, state.pts, state.model_time_ns)
end

"Release consumer read ownership after all readers have finished using values."
function release_graph_output_snapshot!(snapshot::PreparedGraphOutputSnapshot)
    Base.Threads.atomic_cas!(snapshot.state.phase, _SNAPSHOT_READY, _SNAPSHOT_FREE) ==
        _SNAPSHOT_READY || throw(AlgorithmGraphError("no ready snapshot to release"))
    return nothing
end

"""
Drain and discard a snapshot before resetting the graph. Cold, potentially
blocking: stop producer admission and join the consumer first. Successful drain
invalidates old metadata; a failed drain retains ownership in failed state.
"""
function reset_graph_output_snapshot!(snapshot::PreparedGraphOutputSnapshot)
    state = snapshot.state
    state.phase[] == _SNAPSHOT_CLOSED && throw(AlgorithmGraphError("snapshot is closed"))
    state.phase[] = _SNAPSHOT_FAILED
    Backends._drain_output_snapshot!(snapshot.transfer)
    state.sequence = UInt64(0)
    state.generation = UInt64(0)
    state.pts = Int64(0)
    state.model_time_ns = Int64(0)
    state.phase[] = _SNAPSHOT_FREE
    return snapshot
end

"Drain retained work after stopping both participants; close is idempotent."
function Base.close(snapshot::PreparedGraphOutputSnapshot)
    snapshot.state.phase[] == _SNAPSHOT_CLOSED && return nothing
    snapshot.state.phase[] = _SNAPSHOT_FAILED
    Backends._drain_output_snapshot!(snapshot.transfer)
    snapshot.state.phase[] = _SNAPSHOT_CLOSED
    return nothing
end
