#
# Prepared single-device Plant path batches
#
# This layer binds already prepared Gate 6 path groups to the explicit
# atmosphere and direct-imaging batch capabilities. It deliberately does not
# create a task, queue, clock, affinity policy, transport, or resource planner.
#

struct _DirectImagingDevicePathBatchKey{
    B<:AbstractArrayBackend,
    D<:AbstractComputeDevice,
    I<:OpticalPlaneMetadata,
    P,
    S,
    U,
    M,
    F,
}
    backend::B
    device::D
    input_metadata::I
    product_contract::P
    source_type::Type{S}
    input_type::Type{U}
    materialization_type::Type{M}
    propagation_type::Type{F}
    zero_padding::Int
    schedule::PeriodicSchedule
    origin::PlantTimestamp
end

@inline function Base.:(==)(
    left::_DirectImagingDevicePathBatchKey,
    right::_DirectImagingDevicePathBatchKey,
)
    return typeof(left.backend) === typeof(right.backend) &&
        left.device == right.device &&
        left.input_metadata == right.input_metadata &&
        left.product_contract == right.product_contract &&
        left.source_type === right.source_type &&
        left.input_type === right.input_type &&
        left.materialization_type === right.materialization_type &&
        left.propagation_type === right.propagation_type &&
        left.zero_padding == right.zero_padding &&
        left.schedule == right.schedule &&
        left.origin == right.origin
end

@inline Base.isequal(
    left::_DirectImagingDevicePathBatchKey,
    right::_DirectImagingDevicePathBatchKey,
) = left == right

struct _PreparedDirectImagingDevicePathBatch{
    K<:_DirectImagingDevicePathBatchKey,
    C,
    A<:PreparedAtmosphereDirectionBatch,
    O<:PreparedDirectImagingBatch,
    V<:Memory,
    I<:Memory,
    R<:Memory,
}
    key::K
    context::C
    atmosphere_batch::A
    optical_batch::O
    atmosphere_outputs::V
    path_inputs::I
    path_results::R
end

@inline _device_path_batch_selected(
    ::Val{:accelerator},
    ::AcceleratorComputeDevice,
) = true
@inline _device_path_batch_selected(
    ::Val{:accelerator},
    ::HostComputeDevice,
) = false
@inline _device_path_batch_selected(
    ::Val{:all},
    ::AbstractComputeDevice,
) = true
@inline _device_path_batch_selected(
    ::Val{:none},
    ::AbstractComputeDevice,
) = false

function _device_path_batch_selected(
    selection::Val,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :device_path_batch,
        :invalid_selection,
        "unsupported internal device path-batch selection $(selection)",
    ))
end

@inline function _device_path_batch_candidate(
    selection::Val,
    group::PreparedPathExecutionGroup,
)
    requirements = group.requirements
    requirements.requires_full_optical || return nothing
    _device_path_batch_selected(
        selection, requirements.compute_device) || return nothing
    path = group.path
    return _device_path_batch_candidate(
        selection,
        group,
        path,
        path.input,
        path.result,
        path.source,
        path.materialization,
        path.execution,
        source_radiometry(path.source),
        requirements.compute_device,
        path.telescope,
        path.atmosphere,
    )
end

@inline _device_path_batch_candidate(
    ::Val,
    ::PreparedPathExecutionGroup,
    ::PreparedPathExecutor,
    args...,
) = nothing

function _device_path_batch_candidate(
    ::Val,
    group::PreparedPathExecutionGroup,
    path::PreparedPathExecutor,
    input::PupilFunction,
    result::IntensityMap,
    source::Source,
    materialization::PreparedPupilOPDMaterialization,
    execution::PreparedDirectImaging{<:PupilFunction},
    ::PhysicalPhotonIrradianceSource,
    device::AbstractComputeDevice,
    ::Telescope,
    ::AbstractTimedAtmosphere,
)
    validate_path_execution_binding(execution, input, result)
    execution.input === input || throw(PlantPreparationError(
        :device_path_batch,
        :prepared_binding,
        "direct-imaging device-batch candidate does not retain its exact path input",
    ))
    execution.output === result || throw(PlantPreparationError(
        :device_path_batch,
        :prepared_binding,
        "direct-imaging device-batch candidate does not retain its exact path result",
    ))
    n1, n2 = input.metadata.dimensions
    n1 == n2 || throw(PlantPreparationError(
        :device_path_batch,
        :grid,
        "direct-imaging device-batch candidate requires a square pupil grid",
    ))
    padded1, padded2 = execution.field.metadata.dimensions
    padded1 == padded2 && padded1 % n1 == 0 || throw(
        PlantPreparationError(
            :device_path_batch,
            :fft_shape,
            "direct-imaging device-batch candidate has an incompatible FFT shape",
        ),
    )
    result.metadata.dimensions == (padded1, padded2) || throw(
        PlantPreparationError(
            :device_path_batch,
            :product_contract,
            "direct-imaging device-batch candidate output and FFT dimensions differ",
        ),
    )
    typeof(backend(input)) === typeof(group.requirements.backend) ||
        throw(PlantPreparationError(
            :device_path_batch,
            :backend,
            "direct-imaging device-batch candidate backend changed after preparation",
        ))
    compute_device(input.opd) == device &&
        compute_device(result.values) == device ||
        throw(PlantPreparationError(
            :device_path_batch,
            :device,
            "direct-imaging device-batch candidate storage occupies a " *
            "different compute device",
        ))
    propagation = execution.plan.propagation
    return _DirectImagingDevicePathBatchKey(
        group.requirements.backend,
        device,
        input.metadata,
        _direct_imaging_batch_product_contract(result.metadata),
        typeof(source),
        typeof(input),
        typeof(materialization),
        typeof(propagation),
        padded1 ÷ n1,
        group.schedule,
        group.origin,
    )
end

function _matching_device_path_batch_key(keys, candidate)
    @inbounds for index in eachindex(keys)
        isequal(keys[index], candidate) && return index
    end
    return nothing
end

function _copy_device_path_batch_slots(slots)
    copied = Memory{UInt32}(undef, length(slots))
    copyto!(copied, slots)
    return copied
end

function _device_path_batch_stack_views(
    stack::AbstractArray{T,3},
) where {T}
    count = size(stack, 3)
    first_view = @view stack[:, :, 1]
    V = typeof(first_view)
    views = Memory{V}(undef, count)
    views[1] = first_view
    @inbounds for index in 2:count
        view = @view stack[:, :, index]
        typeof(view) === V || throw(PlantPreparationError(
            :device_path_batch,
            :storage,
            "atmosphere batch slices have heterogeneous storage types",
        ))
        views[index] = view
    end
    return views
end

function _require_device_path_batch_product_binding(
    result::IntensityMap,
    batch_product::IntensityMap,
    device::AbstractComputeDevice,
)
    result.metadata == batch_product.metadata || throw(
        PlantPreparationError(
            :device_path_batch,
            :product_contract,
            "batch-backed and path-local direct-imaging products have different metadata",
        ),
    )
    size(result.values) == size(batch_product.values) ||
        throw(PlantPreparationError(
            :device_path_batch,
            :product_contract,
            "batch-backed and path-local direct-imaging products have different dimensions",
        ))
    eltype(result.values) === eltype(batch_product.values) ||
        throw(PlantPreparationError(
            :device_path_batch,
            :product_contract,
            "batch-backed and path-local direct-imaging products use " *
            "different numeric types",
        ))
    compute_device(result.values) == device &&
        compute_device(batch_product.values) == device ||
        throw(PlantPreparationError(
            :device_path_batch,
            :device,
            "direct-imaging batch products occupy a different compute device",
        ))
    return nothing
end

function _collect_device_path_batch_members(
    groups::Memory{PreparedPathExecutionGroup},
    slots::Memory{UInt32},
    key::_DirectImagingDevicePathBatchKey,
    atmosphere::AbstractTimedAtmosphere,
)
    first_group = @inbounds groups[Int(first(slots))]
    first_input = first_group.path.input
    first_result = first_group.path.result
    first_source = first_group.path.source
    I = typeof(first_input)
    R = typeof(first_result)
    S = typeof(first_source)
    inputs = Memory{I}(undef, length(slots))
    results = Memory{R}(undef, length(slots))
    sources = Memory{S}(undef, length(slots))
    @inbounds for member in eachindex(slots)
        slot = Int(slots[member])
        group = groups[slot]
        candidate = _device_path_batch_candidate(Val(:all), group)
        candidate == key || throw(PlantPreparationError(
            :device_path_batch,
            :compatibility,
            "prepared path group $slot changed device-batch compatibility",
        ))
        path = group.path
        path.atmosphere === atmosphere || throw(PlantPreparationError(
            :device_path_batch,
            :prepared_binding,
            "prepared path group $slot belongs to another atmosphere",
        ))
        typeof(path.input) === I &&
            typeof(path.result) === R &&
            typeof(path.source) === S ||
            throw(PlantPreparationError(
                :device_path_batch,
                :compatibility,
                "prepared path group $slot has heterogeneous batch storage",
            ))
        inputs[member] = path.input
        results[member] = path.result
        sources[member] = path.source
    end
    return inputs, results, sources
end

function _prepare_direct_imaging_device_path_batch_owner(
    event_loop_binding::_PlantEventLoopBinding,
    groups::Memory{PreparedPathExecutionGroup},
    candidate_slots,
    key::_DirectImagingDevicePathBatchKey,
    atmosphere::AbstractTimedAtmosphere,
)
    length(candidate_slots) >= 2 || throw(PlantPreparationError(
        :device_path_batch,
        :member_count,
        "a prepared device path batch requires at least two path groups",
    ))
    slots = _copy_device_path_batch_slots(candidate_slots)
    first_group = @inbounds groups[Int(first(slots))]
    context = _prepare_device_execution_context(first_group.path.input.opd)
    _prepared_device_execution_compute_device(context) == key.device ||
        throw(PlantPreparationError(
            :device_path_batch,
            :device,
            "prepared backend execution context belongs to another compute device",
        ))

    implementation = _with_prepared_device_execution_context(context) do
        inputs, results, sources = _collect_device_path_batch_members(
            groups, slots, key, atmosphere)
        optical_batch = prepare_direct_imaging_batch(
            inputs,
            sources;
            zero_padding=key.zero_padding,
        )
        n1, n2 = first(inputs).metadata.dimensions
        T = eltype(first(inputs).opd)
        atmosphere_output = similar(
            first(inputs).opd,
            T,
            n1,
            n2,
            length(slots),
        )
        atmosphere_batch = prepare_atmosphere_direction_batch(
            atmosphere,
            first_group.path.telescope,
            sources,
            atmosphere_output,
        )
        atmosphere_outputs =
            _device_path_batch_stack_views(atmosphere_output)
        @inbounds for member in eachindex(slots)
            _require_device_path_batch_product_binding(
                results[member],
                optical_batch.output[member],
                key.device,
            )
        end
        return _PreparedDirectImagingDevicePathBatch(
            key,
            context,
            atmosphere_batch,
            optical_batch,
            atmosphere_outputs,
            inputs,
            results,
        )
    end
    return PreparedDevicePathBatchOwner(
        event_loop_binding,
        key.device,
        slots,
        implementation,
        _PREPARED_DEVICE_PATH_BATCH_OWNER_TOKEN,
    )
end

function _prepare_device_path_batch_owners(
    selection::Val,
    event_loop_binding::_PlantEventLoopBinding,
    groups::Memory{PreparedPathExecutionGroup},
    atmosphere::AbstractTimedAtmosphere,
)
    keys = Any[]
    candidate_groups = Vector{Vector{UInt32}}()
    @inbounds for slot in eachindex(groups)
        candidate = _device_path_batch_candidate(selection, groups[slot])
        candidate === nothing && continue
        key_slot = _matching_device_path_batch_key(keys, candidate)
        if key_slot === nothing
            push!(keys, candidate)
            push!(candidate_groups, UInt32[UInt32(slot)])
        else
            push!(candidate_groups[key_slot], UInt32(slot))
        end
    end

    owner_count = count(group -> length(group) >= 2, candidate_groups)
    owners = Memory{PreparedDevicePathBatchOwner}(undef, owner_count)
    path_owner_slots = Memory{UInt32}(undef, length(groups))
    fill!(path_owner_slots, UInt32(0))
    owner_slot = 0
    @inbounds for key_slot in eachindex(keys)
        slots = candidate_groups[key_slot]
        length(slots) >= 2 || continue
        owner_slot += 1
        owner = _prepare_direct_imaging_device_path_batch_owner(
            event_loop_binding,
            groups,
            slots,
            keys[key_slot],
            atmosphere,
        )
        owners[owner_slot] = owner
        for group_slot in owner.group_slots
            index = Int(group_slot)
            iszero(path_owner_slots[index]) || throw(
                PlantPreparationError(
                    :device_path_batch,
                    :duplicate_owner,
                    "path execution group $index belongs to more than one " *
                    "device batch owner",
                ),
            )
            path_owner_slots[index] = UInt32(owner_slot)
        end
    end
    return owners, path_owner_slots
end

"""Return the number of prepared single-device path-batch owners."""
@inline device_path_batch_owner_count(
    prepared::PreparedPlantEventLoop,
) = length(prepared.device_path_batch_owners)

"""Return one exact prepared device path-batch owner."""
function device_path_batch_owner(
    prepared::PreparedPlantEventLoop,
    ordinal::Integer,
)
    checkbounds(prepared.device_path_batch_owners, ordinal)
    return @inbounds prepared.device_path_batch_owners[Int(ordinal)]
end

"""Return the exact compute device retained by a prepared path-batch owner."""
@inline device_path_batch_compute_device(
    owner::PreparedDevicePathBatchOwner,
) = owner.device

"""Return the array backend retained by a prepared device path-batch owner."""
@inline device_path_batch_backend(
    owner::PreparedDevicePathBatchOwner,
) = compute_device_backend(owner.device)

"""Return the number of path execution groups retained by an owner."""
@inline device_path_batch_group_count(
    owner::PreparedDevicePathBatchOwner,
) = length(owner.group_slots)

"""Return the prepared event-loop group ordinal for one owner member."""
function device_path_batch_group_ordinal(
    owner::PreparedDevicePathBatchOwner,
    index::Integer,
)
    checkbounds(owner.group_slots, index)
    return Int(@inbounds owner.group_slots[Int(index)])
end

"""
Return the device path-batch owner ordinal for a path execution group.

Returns `nothing` when the group remains independently executable.
"""
function path_execution_group_device_batch_owner_ordinal(
    prepared::PreparedPlantEventLoop,
    group_ordinal::Integer,
)
    checkbounds(prepared.path_device_batch_owner_slots, group_ordinal)
    slot = @inbounds prepared.path_device_batch_owner_slots[
        Int(group_ordinal)]
    return iszero(slot) ? nothing : Int(slot)
end

function _require_device_path_batch_owner(
    prepared::PreparedPlantEventLoop,
    owner::PreparedDevicePathBatchOwner,
)
    owner.event_loop_binding === prepared.binding ||
        _plant_event_loop_error(
            :foreign_device_path_batch_owner,
            "device path-batch owner belongs to another prepared event loop",
        )
    found = 0
    @inbounds for index in eachindex(prepared.device_path_batch_owners)
        candidate = prepared.device_path_batch_owners[index]
        candidate.binding === owner.binding || continue
        iszero(found) || _plant_event_loop_error(
            :duplicate_device_path_batch_owner,
            "prepared event loop retains one device path-batch owner more than once",
        )
        found = index
    end
    iszero(found) && _plant_event_loop_error(
        :foreign_device_path_batch_owner,
        "device path-batch owner is not retained by this prepared event loop",
    )
    return found
end

function _validate_device_path_batch_implementation(
    implementation::_PreparedDirectImagingDevicePathBatch,
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    owner_slot::Int,
)
    count = length(owner.group_slots)
    count >= 2 || _plant_event_loop_error(
        :prepared_binding,
        "prepared device path-batch membership is no longer valid",
    )
    implementation.key.device == owner.device &&
        _prepared_device_execution_compute_device(
            implementation.context) == owner.device ||
        _plant_event_loop_error(
            :device_path_batch_device_mismatch,
            "prepared device path-batch execution context occupies another compute device",
        )
    typeof(implementation.key.backend) ===
        typeof(compute_device_backend(owner.device)) ||
        _plant_event_loop_error(
            :device_path_batch_backend_mismatch,
            "prepared device path-batch backend and compute device differ",
        )
    length(implementation.atmosphere_outputs) == count &&
        length(implementation.path_inputs) == count &&
        length(implementation.path_results) == count &&
        length(implementation.optical_batch.inputs) == count &&
        length(implementation.optical_batch.output) == count ||
        _plant_event_loop_error(
            :prepared_binding,
            "prepared device path-batch storage membership changed",
        )
    validate_direct_imaging_batch(implementation.optical_batch)
    @inbounds for member in 1:count
        group_slot = Int(owner.group_slots[member])
        prepared.path_device_batch_owner_slots[group_slot] ==
            UInt32(owner_slot) ||
            _plant_event_loop_error(
                :prepared_binding,
                "path execution group $group_slot changed device-batch ownership",
            )
        group = prepared.path_groups[group_slot]
        group.requirements.compute_device == owner.device ||
            _plant_event_loop_error(
                :device_path_batch_device_mismatch,
                "path execution group $group_slot occupies another compute device",
            )
        group.input === implementation.path_inputs[member] &&
            implementation.atmosphere_batch.params.sources[member] ==
                implementation.optical_batch.sources[member] &&
            implementation.optical_batch.inputs[member] ===
                implementation.path_inputs[member] ||
            _plant_event_loop_error(
                :prepared_binding,
                "path execution group $group_slot changed device-batch storage binding",
            )
        compute_device(
            implementation.atmosphere_outputs[member]) == owner.device ||
            _plant_event_loop_error(
                :device_path_batch_device_mismatch,
                "atmosphere batch output occupies another compute device",
            )
        _require_device_path_batch_product_binding(
            implementation.path_results[member],
            implementation.optical_batch.output[member],
            owner.device,
        )
    end
    return implementation
end

@inline function _validate_device_path_batch_owner_implementation(
    implementation::_PreparedDirectImagingDevicePathBatch,
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
)
    owner_slot = _require_device_path_batch_owner(prepared, owner)
    _validate_device_path_batch_implementation(
        implementation, owner, prepared, owner_slot)
    return owner_slot
end

@inline function _require_independent_path_execution_group(
    prepared::PreparedPlantEventLoop,
    slot::Int,
)
    owner_slot = @inbounds prepared.path_device_batch_owner_slots[slot]
    iszero(owner_slot) || _plant_event_loop_error(
        :device_path_batch_owner_required,
        "path execution group $slot must be submitted by prepared device " *
        "batch owner $(Int(owner_slot))",
    )
    return nothing
end

function _require_device_path_batch_group_status(
    owner::PreparedDevicePathBatchOwner,
    batch::_OpticalPathBatchWorkspace,
    expected::_OpticalPathBatchGroupStatus,
    operation::AbstractString,
)
    @inbounds for group_slot in owner.group_slots
        slot = Int(group_slot)
        batch.group_status[slot] == expected || _plant_event_loop_error(
            :invalid_device_path_batch_group_status,
            "path execution group $slot is not ready for $operation",
        )
    end
    return nothing
end

function _mark_device_path_batch_group_status!(
    owner::PreparedDevicePathBatchOwner,
    batch::_OpticalPathBatchWorkspace,
    status::_OpticalPathBatchGroupStatus,
)
    @inbounds for group_slot in owner.group_slots
        batch.group_status[Int(group_slot)] = status
    end
    return nothing
end

function _validate_device_path_batch_materialization(
    implementation::_PreparedDirectImagingDevicePathBatch,
    prepared::PreparedPlantEventLoop,
    epoch::AtmosphereEpoch,
)
    validate_atmosphere_direction_batch(
        implementation.atmosphere_batch,
        prepared.atmosphere,
        epoch,
    )
    return nothing
end

Base.@noinline function _materialize_device_path_batch!(
    implementation::_PreparedDirectImagingDevicePathBatch,
    prepared::PreparedPlantEventLoop,
    epoch::AtmosphereEpoch,
)
    _with_prepared_device_execution_context(implementation.context) do
        render_atmosphere_directions!(
            implementation.atmosphere_batch,
            prepared.atmosphere,
            epoch,
        )
        @inbounds for member in eachindex(implementation.path_inputs)
            copyto!(
                implementation.path_inputs[member].opd,
                implementation.atmosphere_outputs[member],
            )
        end
        _synchronize_prepared_device_execution_context!(
            implementation.context)
    end
    return nothing
end

"""
    materialize_device_path_batch!(
        owner, prepared, state, workspace, claim)

Materialize every due member of one exact prepared device owner from the
active Gate 6 atmosphere epoch. Return establishes completion of all layer
reads and same-device path-input handoffs, so the atmosphere materialization
barrier may subsequently be sealed.
"""
function materialize_device_path_batch!(
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch,
        _OpticalPathBatchMaterializing,
        "device path-batch materialization",
    )
    return _materialize_prepared_device_path_batch!(
        owner.implementation,
        owner,
        prepared,
        batch,
    )
end

Base.@noinline function _materialize_prepared_device_path_batch!(
    implementation::_PreparedDirectImagingDevicePathBatch,
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    batch::_OpticalPathBatchWorkspace,
)
    _validate_device_path_batch_owner_implementation(
        implementation, owner, prepared)
    _require_device_path_batch_group_status(
        owner,
        batch,
        _OpticalPathBatchGroupAwaitingMaterialization,
        "device path-batch materialization",
    )
    batch.has_epoch || _plant_event_loop_error(
        :missing_optical_path_batch_epoch,
        "full-optical device path batch has no atmosphere epoch",
    )
    _validate_device_path_batch_materialization(
        implementation, prepared, batch.epoch)

    _mark_device_path_batch_group_status!(
        owner, batch, _OpticalPathBatchGroupMaterializing)
    _materialize_device_path_batch!(
        implementation, prepared, batch.epoch)
    _mark_device_path_batch_group_status!(
        owner, batch, _OpticalPathBatchGroupReady)
    return nothing
end

Base.@noinline function _execute_device_path_batch!(
    implementation::_PreparedDirectImagingDevicePathBatch,
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    timestamp::PlantTimestamp,
)
    _with_prepared_device_execution_context(implementation.context) do
        @inbounds for group_slot in owner.group_slots
            group = prepared.path_groups[Int(group_slot)]
            _stage_materialized_path_execution_group!(
                prepared, state, group, timestamp)
        end
        form_direct_image!(implementation.optical_batch)
        @inbounds for member in eachindex(implementation.path_results)
            copyto!(
                implementation.path_results[member].values,
                implementation.optical_batch.output[member].values,
            )
        end
        @inbounds for group_slot in owner.group_slots
            group = prepared.path_groups[Int(group_slot)]
            _finish_materialized_path_execution_group!(
                prepared, state, group, timestamp)
        end
        _synchronize_prepared_device_execution_context!(
            implementation.context)
    end
    return nothing
end

"""
    execute_device_path_batch!(
        owner, prepared, state, workspace, claim)

Apply every member's path-local acquisition integration, sampled aberrations,
controllable optics, and autonomous optics in canonical order, execute one
prepared optical batch, and complete the explicit same-device handoff to the
original path products. Return establishes the owner's backend completion
boundary before any member is marked complete.
"""
function execute_device_path_batch!(
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
)
    batch = _require_current_optical_path_batch(
        prepared, state, workspace, claim)
    _require_optical_path_batch_phase(
        batch,
        _OpticalPathBatchExecuting,
        "device path-batch execution",
    )
    return _execute_prepared_device_path_batch!(
        owner.implementation,
        owner,
        prepared,
        state,
        batch,
    )
end

Base.@noinline function _execute_prepared_device_path_batch!(
    implementation::_PreparedDirectImagingDevicePathBatch,
    owner::PreparedDevicePathBatchOwner,
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    batch::_OpticalPathBatchWorkspace,
)
    _validate_device_path_batch_owner_implementation(
        implementation, owner, prepared)
    _require_device_path_batch_group_status(
        owner,
        batch,
        _OpticalPathBatchGroupReady,
        "device path-batch execution",
    )

    _mark_device_path_batch_group_status!(
        owner, batch, _OpticalPathBatchGroupExecuting)
    _execute_device_path_batch!(
        implementation,
        owner,
        prepared,
        state,
        batch.timestamp,
    )
    @inbounds for group_slot in owner.group_slots
        slot = Int(group_slot)
        state.path_sampled[slot] = true
        batch.group_status[slot] = _OpticalPathBatchGroupComplete
    end
    return nothing
end

@inline function _materialize_due_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
)
    owner_slot = @inbounds prepared.path_device_batch_owner_slots[
        Int(ordinal)]
    iszero(owner_slot) && return materialize_path_execution_group!(
        prepared, state, workspace, claim, ordinal)
    owner = @inbounds prepared.device_path_batch_owners[Int(owner_slot)]
    Int(first(owner.group_slots)) == ordinal &&
        return materialize_device_path_batch!(
            owner, prepared, state, workspace, claim)
    return nothing
end

@inline function _execute_due_path_execution_group!(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    claim::OpticalPathBatchClaim,
    ordinal::Integer,
)
    owner_slot = @inbounds prepared.path_device_batch_owner_slots[
        Int(ordinal)]
    iszero(owner_slot) && return execute_path_execution_group!(
        prepared, state, workspace, claim, ordinal)
    owner = @inbounds prepared.device_path_batch_owners[Int(owner_slot)]
    Int(first(owner.group_slots)) == ordinal &&
        return execute_device_path_batch!(
            owner, prepared, state, workspace, claim)
    return nothing
end
