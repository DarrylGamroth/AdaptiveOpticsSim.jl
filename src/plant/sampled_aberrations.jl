#
# Prepared native sampled OPD aberrations
#
# Declarations retain caller-owned NCPA/OPDMap storage. Preparation creates
# run-owned copies and exact path-local pupil-footprint couplings. Repeated
# execution performs bounded offset lookup only.
#

@inline _require_sampled_aberration_registration(
    registration::PupilRelayRegistration) = registration

struct _PreparedSampledAberrationToken end
const _PREPARED_SAMPLED_ABERRATION_TOKEN =
    _PreparedSampledAberrationToken()

"""
Run-owned native sampled-aberration OPD and its immutable declaration.
"""
struct PreparedSampledAberration{
    D<:SampledAberrationDefinition,
    A<:AbstractMatrix,
}
    definition::D
    opd::A

    function PreparedSampledAberration(
        ::_PreparedSampledAberrationToken,
        definition::D,
        opd::A,
    ) where {D<:SampledAberrationDefinition,A<:AbstractMatrix}
        return new{D,A}(definition, opd)
    end
end

@inline sampled_aberration_id(
    aberration::PreparedSampledAberration) =
    sampled_aberration_id(aberration.definition)
@inline sampled_aberration_opd(
    aberration::PreparedSampledAberration) = aberration.opd

function _validate_prepared_sampled_aberration_storage(
    definition::SampledAberrationDefinition,
    opd::AbstractMatrix,
)
    metadata = sampled_aberration_metadata(definition)
    size(opd) == metadata.dimensions || throw(PlantPreparationError(
        :sampled_aberration, :surface_dimensions,
        "prepared sampled-aberration OPD dimensions $(size(opd)) do not " *
        "match declared dimensions $(metadata.dimensions)"))
    eltype(opd) === metadata.numeric_type || throw(PlantPreparationError(
        :sampled_aberration, :surface_numeric_type,
        "prepared sampled-aberration OPD numeric type $(eltype(opd)) does " *
        "not match declared type $(metadata.numeric_type)"))
    typeof(backend(opd)) === typeof(metadata.backend) || throw(
        PlantPreparationError(:sampled_aberration, :surface_backend,
            "prepared sampled-aberration OPD changed array backend"))
    compute_device(opd) == metadata.device || throw(PlantPreparationError(
        :sampled_aberration, :surface_device,
        "prepared sampled-aberration OPD changed compute device"))
    return opd
end

function _prepare_sampled_aberration(
    definition::SampledAberrationDefinition)
    source = surface_opd(sampled_aberration_surface(definition))
    copied = similar(source)
    copyto!(copied, source)
    _validate_prepared_sampled_aberration_storage(definition, copied)
    return PreparedSampledAberration(
        _PREPARED_SAMPLED_ABERRATION_TOKEN, definition, copied)
end

function _canonical_sampled_aberration_definitions(
    definition::PlantDefinition)
    aberrations = collect(sampled_aberration_definitions(definition))
    sort!(aberrations; by=aberration ->
        String(sampled_aberration_id(aberration).name))
    return aberrations
end

function _prepare_sampled_aberrations(definition::PlantDefinition)
    declarations = _canonical_sampled_aberration_definitions(definition)
    prepared = Memory{PreparedSampledAberration}(
        undef, length(declarations))
    @inbounds for index in eachindex(declarations)
        prepared[index] =
            _prepare_sampled_aberration(declarations[index])
    end
    return prepared
end

struct _PreparedSampledAberrationPathBindingsToken end
const _PREPARED_SAMPLED_ABERRATION_PATH_BINDINGS_TOKEN =
    _PreparedSampledAberrationPathBindingsToken()

"""
Canonical bounded path-to-sampled-aberration bindings.

Paths use stable `OpticalPathID` order. Within a path, an optional replacement
surface is first, followed by additive surfaces in placement and stable
`SampledAberrationID` order. More than one visible replacement is rejected
because replacement surfaces do not commute.
"""
struct PreparedSampledAberrationPathBindings
    path_ids::Memory{OpticalPathID}
    binding_offsets::Memory{Int}
    aberration_slots::Memory{UInt32}
    couplings::Memory{AbstractPupilSurfacePathCoupling}

    function PreparedSampledAberrationPathBindings(
        ::_PreparedSampledAberrationPathBindingsToken,
        path_ids::Memory{OpticalPathID},
        binding_offsets::Memory{Int},
        aberration_slots::Memory{UInt32},
        couplings::Memory{AbstractPupilSurfacePathCoupling},
    )
        return new(path_ids, binding_offsets, aberration_slots, couplings)
    end
end

@inline _sampled_aberration_application_rank(::DMReplace) = UInt8(0)
@inline _sampled_aberration_application_rank(::DMAdditive) = UInt8(1)

@inline function _sampled_aberration_slot_isless(
    left::UInt32,
    right::UInt32,
    aberrations::AbstractVector,
)
    left_aberration = aberrations[Int(left)]
    right_aberration = aberrations[Int(right)]
    left_definition = left_aberration.definition
    right_definition = right_aberration.definition
    left_rank = _sampled_aberration_application_rank(
        sampled_aberration_application(left_definition))
    right_rank = _sampled_aberration_application_rank(
        sampled_aberration_application(right_definition))
    left_rank == right_rank || return left_rank < right_rank
    left_placement = sampled_aberration_placement(left_definition)
    right_placement = sampled_aberration_placement(right_definition)
    if !_same_optical_placement(left_placement, right_placement)
        return _optical_placement_isless(left_placement, right_placement)
    end
    return String(sampled_aberration_id(left_definition).name) <
        String(sampled_aberration_id(right_definition).name)
end

function _visible_sampled_aberration_slots(
    aberrations::AbstractVector,
    path::OpticalPathID,
)
    slots = UInt32[]
    replacement = nothing
    sizehint!(slots, length(aberrations))
    @inbounds for slot in eachindex(aberrations)
        aberration = aberrations[slot]
        definition = aberration.definition
        _visible_on_path(sampled_aberration_visibility(definition), path) ||
            continue
        slot <= typemax(UInt32) || throw(PlantPreparationError(
            :sampled_aberration, :capacity,
            "prepared sampled-aberration count exceeds UInt32 capacity"))
        application = sampled_aberration_application(definition)
        replacement = _record_sampled_aberration_replacement(
            replacement, application, sampled_aberration_id(definition),
            path)
        push!(slots, UInt32(slot))
    end
    sort!(slots; lt=(left, right) ->
        _sampled_aberration_slot_isless(left, right, aberrations))
    return slots
end

@inline _record_sampled_aberration_replacement(
    replacement, ::DMAdditive, ::SampledAberrationID,
    ::OpticalPathID) = replacement

@inline _record_sampled_aberration_replacement(
    ::Nothing, ::DMReplace, id::SampledAberrationID,
    ::OpticalPathID) = id

function _record_sampled_aberration_replacement(
    previous::SampledAberrationID,
    ::DMReplace,
    current::SampledAberrationID,
    path::OpticalPathID,
)
    throw(PlantPreparationError(:sampled_aberration,
        :ambiguous_replacement_order,
        "optical path $path sees noncommuting replacement aberrations " *
        "$previous and $current; at most one DMReplace surface may be " *
        "visible on a path"))
end

function _rethrow_sampled_aberration_coupling_error(
    error::PlantPreparationError,
    aberration::SampledAberrationID,
    path::OpticalPathID,
)
    throw(PlantPreparationError(:sampled_aberration, error.reason,
        "sampled aberration $aberration on $path: $(error.msg)"))
end

function _rethrow_sampled_aberration_coupling_error(
    error,
    ::SampledAberrationID,
    ::OpticalPathID,
)
    throw(error)
end

function _prepare_sampled_aberration_path_coupling(
    aberration::PreparedSampledAberration,
    path::PreparedPathExecutor,
)
    definition = aberration.definition
    id = sampled_aberration_id(definition)
    path_identity = path_id(path.definition)
    try
        return prepare_sampled_pupil_footprint_coupling(
            sampled_aberration_metadata(definition),
            aberration.opd,
            path,
            sampled_aberration_placement(definition);
            registration=sampled_aberration_registration(definition),
        )
    catch error
        _rethrow_sampled_aberration_coupling_error(
            error, id, path_identity)
    end
end

function _prepare_sampled_aberration_path_bindings(
    aberrations::AbstractVector,
    paths::AbstractVector,
)
    canonical_path_slots = _canonical_prepared_path_slots(paths)
    path_ids = OpticalPathID[]
    binding_offsets = Int[1]
    aberration_slots = UInt32[]
    coupling_values = AbstractPupilSurfacePathCoupling[]
    sizehint!(path_ids, length(paths))
    sizehint!(binding_offsets, length(paths) + 1)
    sizehint!(aberration_slots, length(paths) * length(aberrations))
    sizehint!(coupling_values, length(paths) * length(aberrations))
    @inbounds for path_slot in canonical_path_slots
        path = paths[path_slot]
        id = path_id(path.definition)
        push!(path_ids, id)
        visible_slots =
            _visible_sampled_aberration_slots(aberrations, id)
        for aberration_slot in visible_slots
            aberration = aberrations[Int(aberration_slot)]
            coupling =
                _prepare_sampled_aberration_path_coupling(aberration, path)
            push!(aberration_slots, aberration_slot)
            push!(coupling_values, coupling)
        end
        push!(binding_offsets, length(aberration_slots) + 1)
    end
    couplings = Memory{AbstractPupilSurfacePathCoupling}(
        undef, length(coupling_values))
    copyto!(couplings, coupling_values)
    return PreparedSampledAberrationPathBindings(
        _PREPARED_SAMPLED_ABERRATION_PATH_BINDINGS_TOKEN,
        _optic_binding_memory(path_ids),
        _optic_binding_memory(binding_offsets),
        _optic_binding_memory(aberration_slots),
        couplings,
    )
end

@inline prepared_sampled_aberration_path_count(
    bindings::PreparedSampledAberrationPathBindings) =
    length(bindings.path_ids)
@inline prepared_sampled_aberration_binding_count(
    bindings::PreparedSampledAberrationPathBindings) =
    length(bindings.aberration_slots)
@inline prepared_sampled_aberration_path_id(
    bindings::PreparedSampledAberrationPathBindings,
    ordinal::Integer,
) = bindings.path_ids[ordinal]

function _prepared_sampled_aberration_path_ordinal(
    bindings::PreparedSampledAberrationPathBindings,
    path::OpticalPathID,
)
    @inbounds for ordinal in eachindex(bindings.path_ids)
        bindings.path_ids[ordinal] == path && return ordinal
    end
    throw(PlantPreparationError(:path, :unknown_id,
        "prepared sampled-aberration bindings have no optical path $path"))
end

@inline function prepared_sampled_aberration_binding_range(
    bindings::PreparedSampledAberrationPathBindings,
    path,
)
    ordinal = _prepared_sampled_aberration_path_ordinal(
        bindings, _as_optical_path_id(path))
    first_binding = @inbounds bindings.binding_offsets[ordinal]
    last_binding = @inbounds bindings.binding_offsets[ordinal + 1] - 1
    return first_binding:last_binding
end

@inline prepared_sampled_aberration_slot(
    bindings::PreparedSampledAberrationPathBindings,
    binding::Integer,
) = Int(bindings.aberration_slots[binding])

@inline prepared_sampled_aberration_path_coupling(
    bindings::PreparedSampledAberrationPathBindings,
    binding::Integer,
) = bindings.couplings[binding]

struct _PreparedSampledAberrationApplication{
    A<:AbstractMatrix,
    C<:AbstractPupilSurfacePathCoupling,
    M<:DMApplyMode,
}
    opd::A
    coupling::C
    application::M
end

struct _PreparedSampledAberrationPathPlan{B<:Tuple}
    bindings::B
end

@inline function _prepared_sampled_aberration_application(
    aberration::PreparedSampledAberration,
    coupling::C,
) where {C<:AbstractPupilSurfacePathCoupling}
    definition = aberration.definition
    return _PreparedSampledAberrationApplication(
        aberration.opd,
        coupling,
        sampled_aberration_application(definition),
    )
end

function _prepare_sampled_aberration_path_plan(
    aberrations::AbstractVector,
    bindings::PreparedSampledAberrationPathBindings,
    path::OpticalPathID,
)
    applications = Any[]
    binding_range =
        prepared_sampled_aberration_binding_range(bindings, path)
    sizehint!(applications, length(binding_range))
    @inbounds for binding in binding_range
        slot = prepared_sampled_aberration_slot(bindings, binding)
        push!(applications, _prepared_sampled_aberration_application(
            aberrations[slot],
            prepared_sampled_aberration_path_coupling(bindings, binding),
        ))
    end
    return _PreparedSampledAberrationPathPlan(Tuple(applications))
end

@inline _apply_sampled_aberration_applications!(
    input::PupilFunction, ::Tuple{}) = input

@inline function _apply_sampled_aberration_applications!(
    input::PupilFunction,
    applications::Tuple,
)
    application = first(applications)
    apply_sampled_pupil_surface!(
        input,
        application.opd,
        application.coupling,
        application.application,
    )
    return _apply_sampled_aberration_applications!(
        input, Base.tail(applications))
end

@inline function _apply_sampled_aberration_path_plan!(
    input,
    ::_PreparedSampledAberrationPathPlan{Tuple{}},
)
    return input
end

@inline function _apply_sampled_aberration_path_plan!(
    input::PupilFunction,
    ::_PreparedSampledAberrationPathPlan{Tuple{}},
)
    return input
end

@inline function _apply_sampled_aberration_path_plan!(
    input::PupilFunction,
    plan::_PreparedSampledAberrationPathPlan,
)
    return _apply_sampled_aberration_applications!(input, plan.bindings)
end

function _apply_sampled_aberration_bindings!(
    input::PupilFunction,
    aberrations,
    bindings::PreparedSampledAberrationPathBindings,
    binding_range::UnitRange{Int},
)
    _apply_sampled_aberration_bindings_noreturn!(
        input,
        aberrations,
        bindings,
        binding_range,
    )
    return input
end

function _apply_sampled_aberration_bindings_noreturn!(
    input::PupilFunction,
    aberrations,
    bindings::PreparedSampledAberrationPathBindings,
    binding_range::UnitRange{Int},
)
    @inbounds for binding in binding_range
        slot = prepared_sampled_aberration_slot(bindings, binding)
        aberration = aberrations[slot]
        _apply_prepared_sampled_aberration!(
            input,
            aberration,
            prepared_sampled_aberration_path_coupling(bindings, binding),
        )
    end
    return nothing
end

Base.@noinline function _apply_prepared_sampled_aberration!(
    input::PupilFunction,
    aberration::PreparedSampledAberration{D,A},
    coupling::AbstractPupilSurfacePathCoupling,
) where {D,A}
    apply_sampled_pupil_surface!(
        input,
        aberration.opd,
        coupling,
        sampled_aberration_application(aberration.definition),
    )
    return nothing
end

function _apply_sampled_aberration_bindings!(
    input,
    ::Any,
    ::PreparedSampledAberrationPathBindings,
    ::UnitRange{Int},
)
    throw(PlantPreparationError(:sampled_aberration,
        :unsupported_path_input,
        "sampled aberrations require a PupilFunction path input; got " *
        "$(typeof(input))"))
end

@inline function apply_sampled_aberrations!(
    path::PreparedPathExecutor,
    aberrations,
    bindings::PreparedSampledAberrationPathBindings,
)
    range = prepared_sampled_aberration_binding_range(
        bindings, path_id(path.definition))
    isempty(range) && return path
    _apply_sampled_aberration_bindings!(
        path.input, aberrations, bindings, range)
    return path
end
