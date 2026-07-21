#
# Stable prepared RNG ownership
#
# Stateful streams are derived during plant preparation from explicit owner
# identities. The derivation deliberately does not use Julia's `hash`, whose
# process seed is not a replay contract.
#

"""Versioned domain for stable prepared RNG-stream derivation."""
struct RNGDerivationVersion
    value::UInt32

    function RNGDerivationVersion(value::UInt32)
        iszero(value) && throw(InvalidConfiguration(
            "RNG derivation version must be positive"))
        return new(value)
    end
end

function RNGDerivationVersion(value::Integer)
    value > 0 || throw(InvalidConfiguration(
        "RNG derivation version must be positive"))
    value <= typemax(UInt32) || throw(InvalidConfiguration(
        "RNG derivation version exceeds UInt32"))
    return RNGDerivationVersion(UInt32(value))
end

Base.:(==)(left::RNGDerivationVersion, right::RNGDerivationVersion) =
    left.value == right.value
Base.isequal(left::RNGDerivationVersion, right::RNGDerivationVersion) =
    isequal(left.value, right.value)
Base.hash(version::RNGDerivationVersion, seed::UInt) =
    hash(version.value, hash(RNGDerivationVersion, seed))

function Base.show(io::IO, version::RNGDerivationVersion)
    print(io, "RNGDerivationVersion(", version.value, ")")
end

const _DEFAULT_RNG_DERIVATION_VERSION = RNGDerivationVersion(UInt32(1))
const _RNG_DERIVATION_ALGORITHM = :fnv1a_splitmix64_v1
const _RNG_STREAM_ALGORITHM = :xoshiro

@inline _prepare_rng_derivation_version(version::RNGDerivationVersion) =
    version

function _prepare_rng_derivation_version(version::Integer)
    version > 0 || throw(PlantPreparationError(
        :rng, :invalid_derivation_version,
        "RNG derivation version must be positive"))
    version <= typemax(UInt32) || throw(PlantPreparationError(
        :rng, :invalid_derivation_version,
        "RNG derivation version exceeds UInt32"))
    return RNGDerivationVersion(UInt32(version))
end

function _prepare_rng_derivation_version(version)
    throw(PlantPreparationError(:rng, :invalid_derivation_version,
        "rng_derivation_version must be an Integer or RNGDerivationVersion; got $(typeof(version))"))
end

"""
    RNGOwnerIdentity(category, component, role)

Stable run-local identity of one stateful random stream. All three fields are
declared semantic names; execution order and resource placement are excluded.
"""
struct RNGOwnerIdentity
    category::Symbol
    component::Symbol
    role::Symbol

    function RNGOwnerIdentity(category::Symbol, component::Symbol,
        role::Symbol)
        isempty(String(category)) && throw(InvalidConfiguration(
            "RNG owner category must not be empty"))
        isempty(String(component)) && throw(InvalidConfiguration(
            "RNG owner component must not be empty"))
        isempty(String(role)) && throw(InvalidConfiguration(
            "RNG owner role must not be empty"))
        return new(category, component, role)
    end
end

Base.:(==)(left::RNGOwnerIdentity, right::RNGOwnerIdentity) =
    left.category == right.category &&
    left.component == right.component &&
    left.role == right.role
Base.isequal(left::RNGOwnerIdentity, right::RNGOwnerIdentity) =
    isequal(left.category, right.category) &&
    isequal(left.component, right.component) &&
    isequal(left.role, right.role)
Base.hash(identity::RNGOwnerIdentity, seed::UInt) = hash(identity.role,
    hash(identity.component, hash(identity.category,
        hash(RNGOwnerIdentity, seed))))

function Base.show(io::IO, identity::RNGOwnerIdentity)
    print(io, "RNGOwnerIdentity(", repr(identity.category), ", ",
        repr(identity.component), ", ", repr(identity.role), ")")
end

struct RNGOwnerBinding{O}
    identity::RNGOwnerIdentity
    owner::O
end

mutable struct RNGOwnerToken end

struct PreparedRNGStream{R<:AbstractRNG}
    identity::RNGOwnerIdentity
    derived_seed::UInt64
    state::R
end

struct PreparedOwnerRNGs{B<:RNGOwnerToken,S<:NamedTuple}
    binding::B
    streams::S
end

struct PreparedAtmosphereRNGs{B,O<:Tuple,S<:Tuple}
    binding::B
    owner_bindings::O
    streams::S
end

struct PreparedPlantRNGs{A,P<:Tuple,Q<:Tuple}
    run_seed::UInt64
    derivation_version::RNGDerivationVersion
    atmosphere::A
    paths::P
    acquisitions::Q
end

@inline rng_stream(group::PreparedOwnerRNGs, ::Val{role}) where {role} =
    getproperty(group.streams, role)
@inline rng_stream_state(group::PreparedOwnerRNGs, role::Val) =
    rng_stream(group, role).state

@inline function _rng_derivation_byte(state::UInt64, byte::UInt8)
    return (state ⊻ UInt64(byte)) * UInt64(0x00000100000001b3)
end

function _rng_derivation_uint64(state::UInt64, value::UInt64)
    @inbounds for shift in 0:8:56
        state = _rng_derivation_byte(state,
            UInt8((value >> shift) & UInt64(0xff)))
    end
    return state
end

function _rng_derivation_symbol(state::UInt64, value::Symbol)
    bytes = codeunits(String(value))
    state = _rng_derivation_uint64(state, UInt64(length(bytes)))
    @inbounds for byte in bytes
        state = _rng_derivation_byte(state, byte)
    end
    return state
end

function _derive_rng_stream_seed(run_seed::UInt64,
    version::RNGDerivationVersion, identity::RNGOwnerIdentity)
    state = UInt64(0xcbf29ce484222325)
    state = _rng_derivation_symbol(state, :AdaptiveOpticsSim)
    state = _rng_derivation_uint64(state, run_seed)
    state = _rng_derivation_uint64(state, UInt64(version.value))
    state = _rng_derivation_symbol(state, identity.category)
    state = _rng_derivation_symbol(state, identity.component)
    state = _rng_derivation_symbol(state, identity.role)
    return splitmix64(state)
end

function _prepare_run_seed(seed::Integer)
    seed >= 0 || throw(PlantPreparationError(:rng, :invalid_seed,
        "plant run seed must be non-negative"))
    seed <= typemax(UInt64) || throw(PlantPreparationError(
        :rng, :invalid_seed, "plant run seed exceeds UInt64"))
    return UInt64(seed)
end

function _prepare_run_seed(seed)
    throw(PlantPreparationError(:rng, :invalid_seed,
        "plant run seed must be an Integer; got $(typeof(seed))"))
end

@inline function _prepare_rng_stream(binding::RNGOwnerBinding,
    run_seed::UInt64, version::RNGDerivationVersion)
    seed = _derive_rng_stream_seed(run_seed, version, binding.identity)
    return PreparedRNGStream(binding.identity, seed, runtime_rng(seed))
end

@inline _rng_owner_binding_token(owner) = owner

@inline _atmosphere_rng_binding_token(
    atmosphere::AbstractTimedAtmosphere) = atmosphere_identity(atmosphere)
@inline _atmosphere_rng_binding_token(atmosphere::AbstractAtmosphere) =
    atmosphere
@inline _atmosphere_rng_stream_binding_token(
    atmosphere::AbstractTimedAtmosphere) = atmosphere_identity(atmosphere)
@inline _atmosphere_rng_stream_binding_token(
    layer::MovingAtmosphereLayer) = layer.state
@inline _atmosphere_rng_stream_binding_token(
    layer::InfiniteAtmosphereLayer) = layer.state
@inline _atmosphere_rng_stream_binding_token(owner) = owner

@inline _require_rng_owner_role(role::Symbol, component::Symbol) = role

function _require_rng_owner_role(role, component::Symbol)
    throw(PlantPreparationError(:rng, :invalid_owner_id,
        "$component additional RNG owner roles must be Symbols; got $(typeof(role))"))
end

function _require_additional_rng_roles(roles::Tuple, component::Symbol)
    @inbounds for role in roles
        resolved = _require_rng_owner_role(role, component)
        isempty(String(resolved)) && throw(PlantPreparationError(:rng,
            :invalid_owner_id,
            "$component additional RNG owner roles must not be empty"))
    end
    return roles
end

function _require_additional_rng_roles(roles, component::Symbol)
    throw(PlantPreparationError(:rng, :invalid_owner_id,
        "$component additional RNG owner roles must be a Tuple; got $(typeof(roles))"))
end

"""Qualified extension seam for extra stateful path/provider RNG roles."""
@inline additional_path_rng_owner_roles(execution) = ()

"""Qualified extension seam for extra stateful acquisition/device RNG roles."""
@inline additional_acquisition_rng_owner_roles(execution) = ()

@inline function _path_rng_owner_roles(execution)
    additional = _require_additional_rng_roles(
        additional_path_rng_owner_roles(execution), :path)
    return (:provider, additional...)
end

@inline function _acquisition_rng_owner_roles(execution)
    additional = _require_additional_rng_roles(
        additional_acquisition_rng_owner_roles(execution), :acquisition)
    return (:detector, additional...)
end

@inline function _rng_owner_bindings(owner, category::Symbol,
    component::Symbol, roles::Tuple)
    token = _rng_owner_binding_token(owner)
    return map(role -> RNGOwnerBinding(
        RNGOwnerIdentity(category, component, role), token), roles)
end

@inline function _atmosphere_rng_owner_bindings(
    atmosphere::AbstractAtmosphere)
    return (RNGOwnerBinding(
        RNGOwnerIdentity(:atmosphere, :plant, :state), atmosphere),)
end

function _atmosphere_layer_rng_binding(::Nothing, layer)
    throw(PlantPreparationError(:rng, :missing_owner_id,
        "every multilayer atmosphere layer requires an explicit layer_ids entry for plant preparation"))
end

@inline function _atmosphere_layer_rng_binding(id::AtmosphereLayerID,
    layer)
    return RNGOwnerBinding(
        RNGOwnerIdentity(:atmosphere, id.name, :layer_state), layer)
end

function _multilayer_rng_owner_bindings(atmosphere)
    layers = atmosphere.layers
    ids = atmosphere.params.layer_ids
    length(ids) == length(layers) || throw(PlantPreparationError(
        :rng, :owner_topology,
        "atmosphere layer identities do not match current layer topology"))
    return ntuple(index ->
        _atmosphere_layer_rng_binding(ids[index], layers[index]),
        length(layers))
end

@inline _atmosphere_rng_owner_bindings(atmosphere::MultiLayerAtmosphere) =
    _multilayer_rng_owner_bindings(atmosphere)
@inline _atmosphere_rng_owner_bindings(
    atmosphere::InfiniteMultiLayerAtmosphere) =
    _multilayer_rng_owner_bindings(atmosphere)

@inline _flatten_rng_binding_groups(::Tuple{}) = ()

@inline function _flatten_rng_binding_groups(groups::Tuple)
    return (first(groups)...,
        _flatten_rng_binding_groups(Base.tail(groups))...)
end

function _require_unique_rng_owner_identities(bindings::Tuple)
    @inbounds for index in eachindex(bindings)
        identity = bindings[index].identity
        for previous in firstindex(bindings):(index - 1)
            identity == bindings[previous].identity && throw(
                PlantPreparationError(:rng, :duplicate_owner_id,
                    "duplicate RNG owner identity $identity"))
        end
    end
    return bindings
end

function _require_unique_rng_stream_seeds(bindings::Tuple,
    run_seed::UInt64, version::RNGDerivationVersion)
    @inbounds for index in eachindex(bindings)
        seed = _derive_rng_stream_seed(run_seed, version,
            bindings[index].identity)
        for previous in firstindex(bindings):(index - 1)
            previous_seed = _derive_rng_stream_seed(run_seed, version,
                bindings[previous].identity)
            seed == previous_seed && throw(PlantPreparationError(
                :rng, :derived_seed_collision,
                "distinct RNG owner identities produced the same derived seed"))
        end
    end
    return bindings
end

function _prepare_owner_rngs(bindings::Tuple, run_seed::UInt64,
    version::RNGDerivationVersion)
    roles = map(binding -> binding.identity.role, bindings)
    streams = map(binding ->
        _prepare_rng_stream(binding, run_seed, version), bindings)
    return PreparedOwnerRNGs(first(bindings).owner,
        NamedTuple{roles}(streams))
end

@inline _prepare_owner_rng_groups(::Tuple{}, ::UInt64,
    ::RNGDerivationVersion) = ()

@inline function _prepare_owner_rng_groups(bindings::Tuple,
    run_seed::UInt64, version::RNGDerivationVersion)
    return (
        _prepare_owner_rngs(first(bindings), run_seed, version),
        _prepare_owner_rng_groups(Base.tail(bindings), run_seed,
            version)...,
    )
end

function _prepare_atmosphere_rngs(atmosphere, bindings::Tuple,
    run_seed::UInt64, version::RNGDerivationVersion)
    streams = map(binding ->
        _prepare_rng_stream(binding, run_seed, version), bindings)
    owner_bindings = map(binding ->
        _atmosphere_rng_stream_binding_token(binding.owner), bindings)
    return PreparedAtmosphereRNGs(
        _atmosphere_rng_binding_token(atmosphere), owner_bindings, streams)
end

function _require_rng_owner_binding(group::PreparedOwnerRNGs, owner)
    group.binding === _rng_owner_binding_token(owner) || throw(
        PlantPreparationError(:rng,
        :prepared_binding,
        "prepared RNG streams do not retain their exact state owner"))
    return group
end

function _require_single_atmosphere_rng(
    rngs::PreparedAtmosphereRNGs, atmosphere)
    rngs.binding === _atmosphere_rng_binding_token(atmosphere) || throw(
        PlantPreparationError(:rng,
        :prepared_binding,
        "prepared atmosphere RNGs do not retain the plant atmosphere"))
    length(rngs.streams) == 1 || throw(PlantPreparationError(:rng,
        :owner_topology,
        "single-owner atmosphere requires exactly one prepared RNG stream"))
    length(rngs.owner_bindings) == 1 || throw(PlantPreparationError(:rng,
        :owner_topology,
        "single-owner atmosphere requires exactly one RNG owner binding"))
    first(rngs.owner_bindings) ===
        _atmosphere_rng_stream_binding_token(atmosphere) || throw(
        PlantPreparationError(:rng,
        :prepared_binding,
        "prepared atmosphere RNG stream does not retain the atmosphere"))
    return first(rngs.streams).state
end

function _require_multilayer_atmosphere_rngs(
    rngs::PreparedAtmosphereRNGs, atmosphere)
    rngs.binding === _atmosphere_rng_binding_token(atmosphere) || throw(
        PlantPreparationError(:rng,
        :prepared_binding,
        "prepared atmosphere RNGs do not retain the plant atmosphere"))
    layers = atmosphere.layers
    length(rngs.streams) == length(layers) || throw(PlantPreparationError(
        :rng, :owner_topology,
        "prepared atmosphere RNG count does not match current layers"))
    length(rngs.owner_bindings) == length(layers) || throw(
        PlantPreparationError(:rng, :owner_topology,
            "prepared atmosphere RNG bindings do not match current layers"))
    @inbounds for index in eachindex(layers)
        rngs.owner_bindings[index] ===
            _atmosphere_rng_stream_binding_token(layers[index]) || throw(
            PlantPreparationError(:rng, :prepared_binding,
                "prepared atmosphere layer RNG binding changed"))
    end
    return rngs.streams
end

# Preserve the concrete RNG's complete Random.jl behavior for ordinary and
# extension-defined single-owner atmospheres. Built-in multilayer models use
# the internal router so each layer receives its own stream.
@inline _prepared_atmosphere_rng(atmosphere::AbstractAtmosphere,
    rngs::PreparedAtmosphereRNGs) =
    _require_single_atmosphere_rng(rngs, atmosphere)

@inline function _prepared_atmosphere_rng(
    atmosphere::MultiLayerAtmosphere,
    rngs::PreparedAtmosphereRNGs,
)
    _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    return rngs
end

@inline function _prepared_atmosphere_rng(
    atmosphere::InfiniteMultiLayerAtmosphere,
    rngs::PreparedAtmosphereRNGs,
)
    _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    return rngs
end

@inline function validate_atmosphere_rng_binding(
    rngs::PreparedAtmosphereRNGs, atmosphere::AbstractAtmosphere)
    _require_single_atmosphere_rng(rngs, atmosphere)
    return rngs
end

@inline function validate_atmosphere_rng_binding(
    rngs::PreparedAtmosphereRNGs, atmosphere::MultiLayerAtmosphere)
    _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    return rngs
end

@inline function validate_atmosphere_rng_binding(
    rngs::PreparedAtmosphereRNGs,
    atmosphere::InfiniteMultiLayerAtmosphere)
    _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    return rngs
end

function initialize_atmosphere!(atmosphere::MultiLayerAtmosphere,
    rngs::PreparedAtmosphereRNGs)
    streams = _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    @inbounds for index in eachindex(atmosphere.layers)
        ensure_initialized!(atmosphere.layers[index], streams[index].state)
    end
    return atmosphere
end

function evolve_atmosphere!(atmosphere::MultiLayerAtmosphere,
    duration::Real, rngs::PreparedAtmosphereRNGs)
    _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    @inbounds for layer in atmosphere.layers
        evolve_layer!(layer, duration)
    end
    return atmosphere
end

@inline evolve_initial_atmosphere!(atmosphere::MultiLayerAtmosphere,
    duration::Real, rngs::PreparedAtmosphereRNGs) =
    evolve_atmosphere!(atmosphere, duration, rngs)

function initialize_atmosphere!(atmosphere::InfiniteMultiLayerAtmosphere,
    rngs::PreparedAtmosphereRNGs)
    streams = _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    @inbounds for index in eachindex(atmosphere.layers)
        ensure_initialized!(atmosphere.layers[index].screen,
            streams[index].state)
    end
    return atmosphere
end

function evolve_atmosphere!(atmosphere::InfiniteMultiLayerAtmosphere,
    duration::Real, rngs::PreparedAtmosphereRNGs)
    streams = _require_multilayer_atmosphere_rngs(rngs, atmosphere)
    @inbounds for index in eachindex(atmosphere.layers)
        evolve_layer!(atmosphere.layers[index], duration,
            streams[index].state)
    end
    return atmosphere
end

@inline evolve_initial_atmosphere!(
    atmosphere::InfiniteMultiLayerAtmosphere,
    duration::Real,
    rngs::PreparedAtmosphereRNGs,
) = evolve_atmosphere!(atmosphere, duration, rngs)

@inline function advance_by!(atmosphere::MultiLayerAtmosphere,
    duration::Real, rngs::PreparedAtmosphereRNGs)
    return _advance_by_with_rng!(atmosphere, duration, rngs)
end

@inline function advance_by!(atmosphere::InfiniteMultiLayerAtmosphere,
    duration::Real, rngs::PreparedAtmosphereRNGs)
    return _advance_by_with_rng!(atmosphere, duration, rngs)
end

@inline function advance_to!(atmosphere::MultiLayerAtmosphere,
    target_time::Real, rngs::PreparedAtmosphereRNGs)
    return _advance_to_with_rng!(atmosphere, target_time, rngs)
end

@inline function advance_to!(atmosphere::InfiniteMultiLayerAtmosphere,
    target_time::Real, rngs::PreparedAtmosphereRNGs)
    return _advance_to_with_rng!(atmosphere, target_time, rngs)
end

@inline _rng_group_streams(group::PreparedOwnerRNGs) = values(group.streams)

@inline _flatten_rng_group_streams(::Tuple{}) = ()

@inline function _flatten_rng_group_streams(groups::Tuple)
    return (_rng_group_streams(first(groups))...,
        _flatten_rng_group_streams(Base.tail(groups))...)
end

@inline function _rng_stream_metadata(stream::PreparedRNGStream)
    identity = stream.identity
    return (
        identity=(category=identity.category,
            component=identity.component, role=identity.role),
        derived_seed=stream.derived_seed,
    )
end

@inline function _rng_metadata_order(record)
    identity = record.identity
    return (String(identity.category), String(identity.component),
        String(identity.role))
end

function _rng_replay_metadata(rngs::PreparedPlantRNGs)
    streams = (
        rngs.atmosphere.streams...,
        _flatten_rng_group_streams(rngs.paths)...,
        _flatten_rng_group_streams(rngs.acquisitions)...,
    )
    owners = collect(map(_rng_stream_metadata, streams))
    sort!(owners; by=_rng_metadata_order)
    return (
        run_seed=rngs.run_seed,
        derivation_version=rngs.derivation_version.value,
        derivation_algorithm=_RNG_DERIVATION_ALGORITHM,
        stream_algorithm=_RNG_STREAM_ALGORITHM,
        owners=Tuple(owners),
    )
end
