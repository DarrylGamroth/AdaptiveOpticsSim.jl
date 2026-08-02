#
# Authority-owned pupil-OPD publication routes
#
# A route is a single-writer lifecycle owner. One external owner must serialize
# every operation; these values are deliberately not thread-safe. Direct routes
# materialize into the exact target-local pupil. Remote routes materialize into
# authority-owned pupils and transfer only their OPD arrays through the bounded
# handoff primitive.
#

mutable struct _PupilOPDPublicationRouteIdentityToken end
const _PUPIL_OPD_PUBLICATION_ROUTE_IDENTITY_TOKEN =
    _PupilOPDPublicationRouteIdentityToken()

"""Opaque run-local identity of one exact prepared pupil-OPD route."""
mutable struct _PupilOPDPublicationRouteIdentity
    function _PupilOPDPublicationRouteIdentity(
        token::_PupilOPDPublicationRouteIdentityToken,
    )
        token === _PUPIL_OPD_PUBLICATION_ROUTE_IDENTITY_TOKEN || throw(
            ArgumentError("invalid internal pupil-OPD route token"))
        return new()
    end
end

mutable struct _MaterializedPupilOPDPublicationToken end
const _MATERIALIZED_PUPIL_OPD_PUBLICATION_TOKEN =
    _MaterializedPupilOPDPublicationToken()

"""
    MaterializedPupilOPDPublication

Immutable metadata for one atmosphere-derived pupil-OPD publication. The
record binds the opaque prepared-route identity, atmosphere epoch, plant
timestamp, optical path, route generation, and remote handoff slot when one
exists. It retains neither prepared plant resources nor atmospheric layers.
"""
struct MaterializedPupilOPDPublication{A,E}
    route_identity::_PupilOPDPublicationRouteIdentity
    atmosphere_identity::A
    epoch::E
    timestamp::PlantTimestamp
    path::OpticalPathID
    route_generation::UInt64
    handoff_slot::UInt32
    handoff_generation::UInt64

    function MaterializedPupilOPDPublication(
        token::_MaterializedPupilOPDPublicationToken,
        route_identity::_PupilOPDPublicationRouteIdentity,
        atmosphere_identity::A,
        epoch::E,
        timestamp::PlantTimestamp,
        path::OpticalPathID,
        route_generation::UInt64,
        handoff_slot::UInt32,
        handoff_generation::UInt64,
    ) where {A,E}
        token === _MATERIALIZED_PUPIL_OPD_PUBLICATION_TOKEN || throw(
            ArgumentError("invalid internal pupil-OPD publication token"))
        return new{A,E}(
            route_identity,
            atmosphere_identity,
            epoch,
            timestamp,
            path,
            route_generation,
            handoff_slot,
            handoff_generation,
        )
    end
end

"""Expected result of one pupil-OPD publication-route transition."""
@enum PupilOPDPublicationStatus::UInt8 begin
    PupilOPDPublicationSucceeded = 0x01
    PupilOPDPublicationTransferPending = 0x02
    PupilOPDPublicationRouteBusy = 0x03
    PupilOPDPublicationRejected = 0x04
    PupilOPDPublicationNotMaterialized = 0x05
    PupilOPDPublicationNotSubmitted = 0x06
    PupilOPDPublicationNotCompleted = 0x07
    PupilOPDPublicationTransferFailed = 0x08
    PupilOPDPublicationNotApplied = 0x09
    PupilOPDPublicationAlreadyApplied = 0x0a
    PupilOPDPublicationUncertain = 0x0b
    PupilOPDPublicationGenerationExhausted = 0x0c
end

@enum _PupilOPDPublicationPhase::UInt8 begin
    _PupilOPDRouteIdle = 0x01
    _PupilOPDRouteMaterialized = 0x02
    _PupilOPDRouteSubmitted = 0x03
    _PupilOPDRouteCompleted = 0x04
    _PupilOPDRouteApplied = 0x05
    _PupilOPDRouteFailed = 0x06
    _PupilOPDRouteUncertain = 0x07
end

"""
Prepared single-writer route for one exact target-local pupil-OPD input.

Concrete direct and remote routes implement the lifecycle. Core creates no
task, queue, wait strategy, event-loop integration, or scheduling policy.
"""
abstract type PreparedPupilOPDPublicationRoute end

mutable struct _DirectPupilOPDPublicationState{Q}
    publication::Base.RefValue{Q}
    generation::UInt64
    last_epoch_sequence::UInt64
    phase::_PupilOPDPublicationPhase
end

mutable struct _RemotePupilOPDPublicationState{Q,F,D}
    reference::Base.RefValue{F}
    publication::Base.RefValue{Q}
    borrowed_publication::Base.RefValue{Q}
    borrowed_payload::Base.RefValue{D}
    generation::UInt64
    last_epoch_sequence::UInt64
    phase::_PupilOPDPublicationPhase
end

struct PreparedDirectPupilOPDPublicationRoute{
    P,A,T,R,M,Q,
} <: PreparedPupilOPDPublicationRoute
    identity::_PupilOPDPublicationRouteIdentity
    partitions::P
    authority::A
    partition::T
    path::R
    materialization::M
    state::_DirectPupilOPDPublicationState{Q}
end

struct _PupilOPDHandoffContract{
    Q,A,E,T,AX<:Tuple,
} <: AbstractHandoffPayloadContract{Q}
    route_identity::_PupilOPDPublicationRouteIdentity
    atmosphere_identity::A
    path::OpticalPathID
    epoch_type::Type{E}
    payload_eltype::Type{T}
    payload_axes::AX
end

struct PreparedRemotePupilOPDPublicationRoute{
    P,A,T,R,U,M,H,Q,F,D,
} <: PreparedPupilOPDPublicationRoute
    identity::_PupilOPDPublicationRouteIdentity
    partitions::P
    authority::A
    partition::T
    path::R
    source_pupil::U
    materialization::M
    handoff::H
    state::_RemotePupilOPDPublicationState{Q,F,D}
end

@inline pupil_opd_publication_route_identity(
    route::PreparedPupilOPDPublicationRoute) = route.identity
@inline pupil_opd_publication_authority_target(
    route::PreparedPupilOPDPublicationRoute) = compute_device(route.authority)
@inline pupil_opd_publication_route_target(
    route::PreparedPupilOPDPublicationRoute) = compute_device(route.partition)
@inline pupil_opd_publication_path_id(
    route::PreparedPupilOPDPublicationRoute) = path_id(route.path.definition)
@inline pupil_opd_publication_handoff_capacity(
    ::PreparedDirectPupilOPDPublicationRoute) = 0
@inline pupil_opd_publication_handoff_capacity(
    ::PreparedRemotePupilOPDPublicationRoute) = 1
@inline function prepare_pupil_opd_publication_output(
    ::PreparedDirectPupilOPDPublicationRoute{P,A,T,R,M,Q},
) where {P,A,T,R,M,Q}
    return Ref{Q}()
end

@inline function prepare_pupil_opd_publication_output(
    ::PreparedRemotePupilOPDPublicationRoute{P,A,T,R,U,M,H,Q},
) where {P,A,T,R,U,M,H,Q}
    return Ref{Q}()
end

@inline pupil_opd_publication_route_identity(
    publication::MaterializedPupilOPDPublication) =
    publication.route_identity
@inline pupil_opd_publication_atmosphere_identity(
    publication::MaterializedPupilOPDPublication) =
    publication.atmosphere_identity
@inline pupil_opd_publication_epoch(
    publication::MaterializedPupilOPDPublication) = publication.epoch
@inline pupil_opd_publication_timestamp(
    publication::MaterializedPupilOPDPublication) = publication.timestamp
@inline pupil_opd_publication_path_id(
    publication::MaterializedPupilOPDPublication) = publication.path

@noinline function _pupil_opd_publication_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(
        :pupil_opd_publication, reason, String(message)))
end

@inline _require_pupil_opd_path_input(input::PupilFunction) = input

function _require_pupil_opd_path_input(input)
    _pupil_opd_publication_error(
        :unsupported_input,
        "pupil-OPD publication requires one native PupilFunction input; " *
        "got $(typeof(input))",
    )
end

"""
    prepare_pupil_opd_publication_materialization(
        model, atmosphere, telescope, source, pupil)

Qualified model opt-in seam for authority-owned pupil-OPD publication. A
supported optical-path model must return one native
`PreparedPupilOPDMaterialization` bound to the exact supplied objects.
"""
function prepare_pupil_opd_publication_materialization(
    model,
    atmosphere::AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AbstractSource,
    pupil::PupilFunction,
)
    _pupil_opd_publication_error(
        :unsupported_model,
        "path model $(typeof(model)) does not opt in to native pupil-OPD " *
        "publication",
    )
end

function _require_pupil_opd_publication_materialization(
    materialization::PreparedPupilOPDMaterialization,
    atmosphere::AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AbstractSource,
    pupil::PupilFunction,
)
    validate_path_materialization_binding(
        materialization, pupil, atmosphere, source)
    _validate_path_materialization_telescope(materialization, telescope)
    return materialization
end

function _require_pupil_opd_publication_materialization(
    materialization,
    ::AbstractTimedAtmosphere,
    ::Telescope,
    ::AbstractSource,
    ::PupilFunction,
)
    _pupil_opd_publication_error(
        :invalid_materialization,
        "pupil-OPD publication opt-in must return " *
        "PreparedPupilOPDMaterialization; got $(typeof(materialization))",
    )
end

function _pupil_opd_publication_path(
    partitions::PreparedPlantPartitions,
    id::OpticalPathID,
)
    found_partition = nothing
    found_path = nothing
    for partition in prepared_partitions(partitions)
        for path in prepared_paths(partition)
            path_id(path.definition) == id || continue
            isnothing(found_path) || _pupil_opd_publication_error(
                :duplicate_path,
                "prepared partitions contain more than one path $id",
            )
            found_partition = partition
            found_path = path
        end
    end
    isnothing(found_path) && _pupil_opd_publication_error(
        :unknown_path,
        "prepared partitions contain no path $id",
    )
    return found_partition, found_path
end

@inline _pupil_opd_publication_path(
    partitions::PreparedPlantPartitions,
    name::Symbol,
) = _pupil_opd_publication_path(partitions, OpticalPathID(name))

function _pupil_opd_publication_path(
    ::PreparedPlantPartitions,
    id,
)
    _pupil_opd_publication_error(
        :invalid_path_id,
        "pupil-OPD publication path must be an OpticalPathID or Symbol; " *
        "got $(typeof(id))",
    )
end

@inline function _pupil_opd_publication_epoch_type(authority)
    timeline = atmosphere_timeline(prepared_atmosphere(authority))
    identity = atmosphere_authority_identity(authority)
    return AtmosphereEpoch{typeof(timeline.model_time),typeof(identity)}
end

function _new_pupil_opd_publication(
    route::PreparedPupilOPDPublicationRoute,
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
    reference::HandoffSlotReference,
)
    return MaterializedPupilOPDPublication(
        _MATERIALIZED_PUPIL_OPD_PUBLICATION_TOKEN,
        route.identity,
        atmosphere_authority_identity(route.authority),
        epoch,
        timestamp,
        path_id(route.path.definition),
        route.state.generation,
        reference.slot,
        reference.generation,
    )
end

function _new_pupil_opd_publication(
    route::PreparedDirectPupilOPDPublicationRoute,
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
)
    return MaterializedPupilOPDPublication(
        _MATERIALIZED_PUPIL_OPD_PUBLICATION_TOKEN,
        route.identity,
        atmosphere_authority_identity(route.authority),
        epoch,
        timestamp,
        path_id(route.path.definition),
        route.state.generation,
        UInt32(0),
        UInt64(0),
    )
end

@inline handoff_payload_eltype(contract::_PupilOPDHandoffContract) =
    contract.payload_eltype
@inline handoff_payload_axes(contract::_PupilOPDHandoffContract) =
    contract.payload_axes

function validate_handoff_publication(
    contract::_PupilOPDHandoffContract{Q},
    publication::Q,
) where {Q}
    publication.route_identity === contract.route_identity ||
        _pupil_opd_publication_error(
            :wrong_route, "publication belongs to another pupil-OPD route")
    publication.atmosphere_identity === contract.atmosphere_identity ||
        _pupil_opd_publication_error(
            :wrong_atmosphere,
            "publication belongs to another atmosphere authority",
        )
    typeof(publication.epoch) === contract.epoch_type ||
        _pupil_opd_publication_error(
            :wrong_epoch_type,
            "publication epoch type does not match the prepared route",
        )
    publication.epoch.identity === contract.atmosphere_identity ||
        _pupil_opd_publication_error(
            :wrong_epoch,
            "publication epoch belongs to another atmosphere",
        )
    publication.path == contract.path || _pupil_opd_publication_error(
        :wrong_path, "publication belongs to another optical path")
    iszero(publication.handoff_slot) && _pupil_opd_publication_error(
        :invalid_slot, "remote publication handoff slot must be nonzero")
    iszero(publication.handoff_generation) &&
        _pupil_opd_publication_error(
            :invalid_generation,
            "remote publication handoff generation must be nonzero",
        )
    return nothing
end

function _validate_pupil_opd_route_binding(
    route::PreparedPupilOPDPublicationRoute,
)
    prepared_atmosphere_authority(route.partitions) === route.authority ||
        _pupil_opd_publication_error(
            :prepared_binding, "route no longer matches its atmosphere authority")
    atmosphere_authority_binding(route.partition) ===
        atmosphere_authority_binding(route.authority) ||
        _pupil_opd_publication_error(
            :prepared_binding,
            "target partition does not retain the route's authority binding",
        )
    target = compute_device(route.partition)
    getfield(route.path.key, :device) == target ||
        _pupil_opd_publication_error(
            :prepared_binding, "path-result key no longer matches its target")
    _prepared_device_execution_compute_device(route.path.context) == target ||
        _pupil_opd_publication_error(
            :prepared_binding, "path context no longer matches its target")
    validate_telescope_target(route.path.telescope, target)
    input = _require_pupil_opd_path_input(path_input(route.path))
    _validate_path_input(input)
    input.metadata.device == target &&
        compute_device(input.support) == target &&
        compute_device(input.amplitude) == target &&
        compute_device(input.opd) == target ||
        _pupil_opd_publication_error(
            :prepared_binding, "path pupil no longer occupies its target")
    revision = aperture_revision(route.path.telescope)
    revision == getfield(getfield(route.path.key, :revisions), :telescope) ||
        _pupil_opd_publication_error(
            :revision, "path telescope revision changed after preparation")
    _require_path_input_revisions(input, revision)
    validate_path_execution_binding(
        route.path.execution, input, route.path.result)
    validate_path_execution_target(route.path.execution, target)
    _validate_pupil_opd_publication_materialization_binding(route)
    return nothing
end

function _validate_pupil_opd_publication_materialization_binding(
    route::PreparedDirectPupilOPDPublicationRoute,
)
    atmosphere = prepared_atmosphere(route.authority)
    telescope = prepared_telescope(route.authority)
    source = route.path.source
    input = _require_pupil_opd_path_input(path_input(route.path))
    validate_path_materialization_binding(
        route.materialization, input, atmosphere, source)
    _validate_path_materialization_telescope(
        route.materialization, telescope)
    return nothing
end

function _validate_pupil_opd_publication_materialization_binding(
    route::PreparedRemotePupilOPDPublicationRoute,
)
    atmosphere = prepared_atmosphere(route.authority)
    telescope = prepared_telescope(route.authority)
    source = route.path.source
    validate_path_materialization_binding(
        route.materialization, route.source_pupil, atmosphere, source)
    _validate_path_materialization_telescope(
        route.materialization, telescope)
    return nothing
end

function _validate_pupil_opd_epoch(
    route::PreparedPupilOPDPublicationRoute,
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
)
    _validate_pupil_opd_route_binding(route)
    atmosphere = prepared_atmosphere(route.authority)
    _validate_epoch_identity(
        atmosphere_authority_identity(route.authority), atmosphere, epoch)
    epoch_sequence(epoch) > route.state.last_epoch_sequence ||
        _pupil_opd_publication_error(
            :epoch_sequence,
            "pupil-OPD route epochs must advance monotonically",
        )
    T = typeof(epoch_time(epoch))
    isequal(epoch_time(epoch), plant_time_seconds(timestamp, T)) ||
        _pupil_opd_publication_error(
            :timestamp,
            "plant timestamp does not identify the supplied atmosphere epoch",
        )
    return nothing
end

function _validate_active_pupil_opd_publication(
    route::PreparedPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    publication.route_identity === route.identity || return false
    publication.atmosphere_identity ===
        atmosphere_authority_identity(route.authority) || return false
    publication.epoch.identity === publication.atmosphere_identity ||
        return false
    publication.path == path_id(route.path.definition) || return false
    publication.route_generation == route.state.generation || return false
    active = route.state.publication[]
    return publication.epoch == active.epoch &&
        publication.timestamp == active.timestamp &&
        publication.handoff_slot == active.handoff_slot &&
        publication.handoff_generation == active.handoff_generation
end

function _prepare_direct_pupil_opd_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedAtmosphereAuthority,
    partition::PreparedTargetPartition,
    path::PreparedTargetLocalPathResources,
)
    input = _require_pupil_opd_path_input(path_input(path))
    source = path.source
    atmosphere = prepared_atmosphere(authority)
    materialization = prepare_pupil_opd_publication_materialization(
        path_model(path.definition),
        atmosphere,
        prepared_telescope(authority),
        source,
        input,
    )
    materialization = _require_pupil_opd_publication_materialization(
        materialization, atmosphere, prepared_telescope(authority), source,
        input)
    identity = _PupilOPDPublicationRouteIdentity(
        _PUPIL_OPD_PUBLICATION_ROUTE_IDENTITY_TOKEN)
    E = _pupil_opd_publication_epoch_type(authority)
    Q = MaterializedPupilOPDPublication{
        typeof(atmosphere_authority_identity(authority)),E,
    }
    return PreparedDirectPupilOPDPublicationRoute(
        identity,
        partitions,
        authority,
        partition,
        path,
        materialization,
        _DirectPupilOPDPublicationState(
            Ref{Q}(), UInt64(1), UInt64(0), _PupilOPDRouteIdle),
    )
end

function _prepare_pupil_opd_source_pupil(model, authority, source,
    ::Type{T}) where {T}
    atmosphere = prepared_atmosphere(authority)
    telescope = prepared_telescope(authority)
    return _with_completed_prepared_device_execution_context(
        authority.context) do
        pupil = PupilFunction(telescope; T=T)
        materialization = prepare_pupil_opd_publication_materialization(
            model, atmosphere, telescope, source, pupil)
        materialization = _require_pupil_opd_publication_materialization(
            materialization, atmosphere, telescope, source, pupil)
        return pupil, materialization
    end
end

function _prepare_pupil_opd_destination_slot(path)
    input = _require_pupil_opd_path_input(path_input(path))
    return _with_completed_prepared_device_execution_context(
        path.context) do
        return similar(input.opd)
    end
end

function _prepare_remote_pupil_opd_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedAtmosphereAuthority,
    partition::PreparedTargetPartition,
    path::PreparedTargetLocalPathResources,
)
    input = _require_pupil_opd_path_input(path_input(path))
    source = path.source
    source_pupil, materialization = _prepare_pupil_opd_source_pupil(
        path_model(path.definition), authority, source, eltype(input.opd))
    destination_slot = _prepare_pupil_opd_destination_slot(path)
    identity = _PupilOPDPublicationRouteIdentity(
        _PUPIL_OPD_PUBLICATION_ROUTE_IDENTITY_TOKEN)
    E = _pupil_opd_publication_epoch_type(authority)
    Q = MaterializedPupilOPDPublication{
        typeof(atmosphere_authority_identity(authority)),E,
    }
    payload_axes = axes(input.opd)
    contract = _PupilOPDHandoffContract{
        Q,typeof(atmosphere_authority_identity(authority)),E,
        eltype(input.opd),typeof(payload_axes),
    }(
        identity,
        atmosphere_authority_identity(authority),
        path_id(path.definition),
        E,
        eltype(input.opd),
        payload_axes,
    )
    handoff = prepare_cross_domain_handoff(
        partitions, contract, (source_pupil.opd,), (destination_slot,))
    F = HandoffSlotReference{typeof(handoff.identity)}
    D = typeof(destination_slot)
    return PreparedRemotePupilOPDPublicationRoute(
        identity,
        partitions,
        authority,
        partition,
        path,
        source_pupil,
        materialization,
        handoff,
        _RemotePupilOPDPublicationState(
            Ref{F}(),
            Ref{Q}(),
            Ref{Q}(),
            Ref{D}(),
            UInt64(1),
            UInt64(0),
            _PupilOPDRouteIdle,
        ),
    )
end

"""
    prepare_pupil_opd_publication_route(partitions, path)

Prepare one exact native pupil-OPD route. A same-target route materializes
directly and allocates no handoff. A remote route owns exactly one fixed paired
handoff slot because the route cannot materialize another epoch until its
current publication has been applied and reclaimed.
"""
function prepare_pupil_opd_publication_route(
    partitions::PreparedPlantPartitions,
    id,
)
    partition, path = _pupil_opd_publication_path(partitions, id)
    authority = prepared_atmosphere_authority(partitions)
    if compute_device(partition) == compute_device(authority)
        return _prepare_direct_pupil_opd_publication_route(
            partitions, authority, partition, path)
    end
    return _prepare_remote_pupil_opd_publication_route(
        partitions, authority, partition, path)
end

function materialize_pupil_opd_publication!(
    output::Base.RefValue{Q},
    route::PreparedDirectPupilOPDPublicationRoute{P,A,T,R,M,Q},
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
) where {P,A,T,R,M,Q}
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteIdle ||
        return PupilOPDPublicationRouteBusy
    _validate_pupil_opd_epoch(route, epoch, timestamp)
    validate_path_materialization(
        route.materialization,
        path_input(route.path),
        prepared_atmosphere(route.authority),
        epoch,
    )
    route.state.phase = _PupilOPDRouteUncertain
    _with_completed_prepared_device_execution_context(route.authority.context) do
        materialize_path_input!(
            route.materialization,
            path_input(route.path),
            prepared_atmosphere(route.authority),
            epoch,
        )
    end
    publication = _new_pupil_opd_publication(route, epoch, timestamp)
    route.state.publication[] = publication
    output[] = publication
    route.state.last_epoch_sequence = epoch_sequence(epoch)
    route.state.phase = _PupilOPDRouteApplied
    return PupilOPDPublicationSucceeded
end

function materialize_pupil_opd_publication!(
    output::Base.RefValue{Q},
    route::PreparedRemotePupilOPDPublicationRoute{P,A,T,R,U,M,H,Q},
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
) where {P,A,T,R,U,M,H,Q}
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteIdle ||
        return PupilOPDPublicationRouteBusy
    _validate_pupil_opd_epoch(route, epoch, timestamp)
    status = try_next_free_handoff_slot!(
        route.state.reference, route.handoff)
    if status != HandoffTransitionSucceeded
        route.state.phase = _PupilOPDRouteUncertain
        return PupilOPDPublicationUncertain
    end
    # The primitive has selected a free slot but does not claim it until
    # submission. The route nevertheless becomes fail-stop here: after any
    # later exception it cannot prove whether preparation reached a state that
    # is safe to repeat, so external recovery must dispose of the whole route.
    route.state.phase = _PupilOPDRouteUncertain
    reference = route.state.reference[]
    reference.slot == UInt32(1) || begin
        route.state.phase = _PupilOPDRouteUncertain
        return PupilOPDPublicationUncertain
    end
    pupil = route.source_pupil
    materialization = route.materialization
    producer_handoff_payload(route.handoff, reference) === pupil.opd ||
        _pupil_opd_publication_error(
            :prepared_binding,
            "remote route source pupil lost its exact handoff-slot binding",
        )
    validate_path_materialization(
        materialization,
        pupil,
        prepared_atmosphere(route.authority),
        epoch,
    )
    _with_completed_prepared_device_execution_context(route.authority.context) do
        materialize_path_input!(materialization, pupil,
            prepared_atmosphere(route.authority), epoch)
    end
    publication = _new_pupil_opd_publication(
        route, epoch, timestamp, reference)
    route.state.publication[] = publication
    output[] = publication
    route.state.last_epoch_sequence = epoch_sequence(epoch)
    route.state.phase = _PupilOPDRouteMaterialized
    return PupilOPDPublicationSucceeded
end

function submit_pupil_opd_publication!(
    route::PreparedDirectPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteApplied ||
        return PupilOPDPublicationNotMaterialized
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    return PupilOPDPublicationSucceeded
end

function submit_pupil_opd_publication!(
    route::PreparedRemotePupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteMaterialized ||
        return PupilOPDPublicationNotMaterialized
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    route.state.phase = _PupilOPDRouteUncertain
    status = submit_handoff!(
        route.handoff, route.state.reference[], publication)
    if status == HandoffTransitionSucceeded
        route.state.phase = _PupilOPDRouteSubmitted
        return PupilOPDPublicationSucceeded
    end
    if status == HandoffSubmissionFailed
        route.state.phase = _PupilOPDRouteFailed
        return PupilOPDPublicationTransferFailed
    end
    route.state.phase = _PupilOPDRouteUncertain
    return PupilOPDPublicationUncertain
end

function try_complete_pupil_opd_publication!(
    route::PreparedDirectPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteApplied ||
        return PupilOPDPublicationNotSubmitted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    return PupilOPDPublicationSucceeded
end

function try_complete_pupil_opd_publication!(
    route::PreparedRemotePupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteSubmitted ||
        return PupilOPDPublicationNotSubmitted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    route.state.phase = _PupilOPDRouteUncertain
    status = try_complete_handoff!(
        route.handoff, route.state.reference[])
    if status == HandoffTransitionSucceeded
        route.state.phase = _PupilOPDRouteCompleted
        return PupilOPDPublicationSucceeded
    end
    if status == HandoffTransferPending
        route.state.phase = _PupilOPDRouteSubmitted
        return PupilOPDPublicationTransferPending
    end
    if status == HandoffCompletionFailed
        route.state.phase = _PupilOPDRouteFailed
        return PupilOPDPublicationTransferFailed
    end
    route.state.phase = _PupilOPDRouteUncertain
    return PupilOPDPublicationUncertain
end

# Deterministic serial-oracle seam. Concurrent owners use the nonblocking
# `try_complete_pupil_opd_publication!` transition and their own wait strategy.
function _complete_pupil_opd_publication!(
    route::PreparedDirectPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteApplied ||
        return PupilOPDPublicationNotSubmitted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    return PupilOPDPublicationSucceeded
end

function _complete_pupil_opd_publication!(
    route::PreparedRemotePupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteSubmitted ||
        return PupilOPDPublicationNotSubmitted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    route.state.phase = _PupilOPDRouteUncertain
    status = complete_handoff!(route.handoff, route.state.reference[])
    if status == HandoffTransitionSucceeded
        route.state.phase = _PupilOPDRouteCompleted
        return PupilOPDPublicationSucceeded
    end
    if status == HandoffCompletionFailed
        route.state.phase = _PupilOPDRouteFailed
        return PupilOPDPublicationTransferFailed
    end
    route.state.phase = _PupilOPDRouteUncertain
    return PupilOPDPublicationUncertain
end

function apply_pupil_opd_publication!(
    route::PreparedDirectPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteApplied ||
        return PupilOPDPublicationNotCompleted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    return PupilOPDPublicationSucceeded
end

function apply_pupil_opd_publication!(
    route::PreparedRemotePupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    if route.state.phase == _PupilOPDRouteApplied
        _validate_active_pupil_opd_publication(route, publication) ||
            return PupilOPDPublicationRejected
        return PupilOPDPublicationAlreadyApplied
    end
    route.state.phase == _PupilOPDRouteCompleted ||
        return PupilOPDPublicationNotCompleted
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    # The remote atmosphere epoch may no longer be current here. Revalidate
    # only the exact target-local resources and execution binding immediately
    # before any completed payload can be borrowed or copied.
    route.state.phase = _PupilOPDRouteUncertain
    _validate_pupil_opd_route_binding(route)
    status = try_borrow_completed_handoff!(
        route.state.borrowed_payload,
        route.state.borrowed_publication,
        route.handoff,
        route.state.reference[],
    )
    status == HandoffTransitionSucceeded || begin
        return PupilOPDPublicationUncertain
    end
    borrowed_publication = route.state.borrowed_publication[]
    _validate_active_pupil_opd_publication(
        route, borrowed_publication) || begin
        route.state.phase = _PupilOPDRouteUncertain
        return PupilOPDPublicationUncertain
    end
    borrowed_publication.epoch == publication.epoch &&
        borrowed_publication.timestamp == publication.timestamp || begin
        route.state.phase = _PupilOPDRouteUncertain
        return PupilOPDPublicationUncertain
    end
    destination = _require_pupil_opd_path_input(path_input(route.path)).opd
    staging = route.state.borrowed_payload[]
    if !(eltype(staging) === eltype(destination) &&
            axes(staging) == axes(destination))
        route.state.phase = _PupilOPDRouteUncertain
        _pupil_opd_publication_error(
            :payload_contract,
            "completed pupil-OPD staging does not match the target-local pupil",
        )
    end
    if compute_device(staging) != compute_device(destination)
        route.state.phase = _PupilOPDRouteUncertain
        _pupil_opd_publication_error(
            :wrong_target,
            "completed pupil-OPD staging is not on the target-local device",
        )
    end
    _with_completed_prepared_device_execution_context(route.path.context) do
        copyto!(destination, staging)
    end
    route.state.phase = _PupilOPDRouteApplied
    return PupilOPDPublicationSucceeded
end

function reclaim_pupil_opd_publication!(
    route::PreparedDirectPupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    route.state.phase == _PupilOPDRouteApplied ||
        return PupilOPDPublicationNotApplied
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    route.state.generation == typemax(UInt64) &&
        return PupilOPDPublicationGenerationExhausted
    route.state.generation += one(UInt64)
    route.state.phase = _PupilOPDRouteIdle
    return PupilOPDPublicationSucceeded
end

function reclaim_pupil_opd_publication!(
    route::PreparedRemotePupilOPDPublicationRoute,
    publication::MaterializedPupilOPDPublication,
)
    route.state.phase == _PupilOPDRouteUncertain &&
        return PupilOPDPublicationUncertain
    terminal = route.state.phase == _PupilOPDRouteApplied ||
        route.state.phase == _PupilOPDRouteFailed
    terminal || return PupilOPDPublicationNotApplied
    _validate_active_pupil_opd_publication(route, publication) ||
        return PupilOPDPublicationRejected
    route.state.generation == typemax(UInt64) &&
        return PupilOPDPublicationGenerationExhausted
    status = reclaim_handoff!(route.handoff, route.state.reference[])
    status == HandoffGenerationExhausted &&
        return PupilOPDPublicationGenerationExhausted
    status == HandoffTransitionSucceeded || begin
        route.state.phase = _PupilOPDRouteUncertain
        return PupilOPDPublicationUncertain
    end
    route.state.generation += one(UInt64)
    route.state.phase = _PupilOPDRouteIdle
    return PupilOPDPublicationSucceeded
end

# Expected RTC-facing validation is status based. An unrelated value must not
# escape as a MethodError or enter any route/handoff transition.
@inline submit_pupil_opd_publication!(
    ::PreparedPupilOPDPublicationRoute, publication) =
    PupilOPDPublicationRejected
@inline try_complete_pupil_opd_publication!(
    ::PreparedPupilOPDPublicationRoute, publication) =
    PupilOPDPublicationRejected
@inline apply_pupil_opd_publication!(
    ::PreparedPupilOPDPublicationRoute, publication) =
    PupilOPDPublicationRejected
@inline reclaim_pupil_opd_publication!(
    ::PreparedPupilOPDPublicationRoute, publication) =
    PupilOPDPublicationRejected
