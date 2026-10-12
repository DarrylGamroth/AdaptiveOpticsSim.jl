# Backend-owned arrays and fixed-address transfer work for optional graph output
# observation. The producer and consumer each have one retained context; neither
# completion polling nor publication may introduce a host completion wait.
struct _HostOutputSnapshot{A}
    values::A
end

function _prepare_output_snapshot(source::AbstractArray{T,N},
    ::_HostPreparedDeviceExecutionContext) where {T,N}
    values = FixedSizeArray{T,N}(undef, size(source))
    fill!(values, zero(T))
    return _HostOutputSnapshot(values)
end

function _prepare_output_snapshot(source::AbstractArray,
    context::_AbstractPreparedDeviceExecutionContext)
    throw(InvalidConfiguration(
        "optional output snapshots are unsupported for $(typeof(context))"))
end

@inline _output_snapshot_values(transfer) = transfer.values
@inline function _submit_output_snapshot!(transfer::_HostOutputSnapshot, source)
    copyto!(transfer.values, source)
    return nothing
end
@inline _start_output_snapshot_transfer!(::_HostOutputSnapshot) = nothing
@inline _poll_output_snapshot_transfer!(::_HostOutputSnapshot) = true
@inline _drain_output_snapshot!(::_HostOutputSnapshot) = nothing
