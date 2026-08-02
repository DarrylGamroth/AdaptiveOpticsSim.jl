#
# Native deformable-mirror plant model
#
# DeformableMirrorModel is cold configuration. Preparation retains immutable
# influence/operator data, while each event-loop state and workspace owns a
# separate live DeformableMirror runtime. Command staging forms a complete
# surface in the workspace runtime; commit publishes it with an O(1) swap.
#

@inline function _deformable_mirror_definition_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantDefinitionError(
        :controllable_optic,
        reason,
        String(message),
    ))
end

@inline function _deformable_mirror_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(
        :controllable_optic,
        reason,
        String(message),
    ))
end

@inline _deformable_mirror_error_message(error::Exception) =
    sprint(showerror, error)
@inline _deformable_mirror_error_message(error) = repr(error)

@noinline function _throw_deformable_mirror_definition_failure(
    error::InterruptException,
    ::Symbol,
)
    throw(error)
end

@noinline function _throw_deformable_mirror_definition_failure(
    error,
    reason::Symbol,
)
    _deformable_mirror_definition_error(
        reason,
        _deformable_mirror_error_message(error),
    )
end

function _resolve_deformable_mirror_topology(
    n_act::Integer,
    ::Nothing,
    ::Type{T},
) where {T<:AbstractFloat}
    n_act > 0 || _deformable_mirror_definition_error(
        :invalid_deformable_mirror_topology,
        "deformable-mirror n_act must be positive",
    )
    resolved = ActuatorGridTopology(Int(n_act); T)
    return resolved::ActuatorGridTopology{T}
end

function _resolve_deformable_mirror_topology(
    ::Bool,
    ::Nothing,
    ::Type{<:AbstractFloat},
)
    _deformable_mirror_definition_error(
        :invalid_deformable_mirror_topology,
        "deformable-mirror n_act must be an integer, not Bool",
    )
end

@inline function _resolve_deformable_mirror_topology(
    ::Nothing,
    topology::AbstractDMTopology,
    ::Type{<:AbstractFloat},
)
    return topology
end

function _resolve_deformable_mirror_topology(
    n_act,
    topology,
    ::Type{<:AbstractFloat},
)
    _deformable_mirror_definition_error(
        :invalid_deformable_mirror_topology,
        "specify exactly one of n_act or topology; got n_act=$(repr(n_act)) " *
        "and topology $(typeof(topology))",
    )
end

@inline function _resolve_deformable_mirror_influence(
    ::Nothing,
    ::Nothing,
    ::Nothing,
    ::Type{T},
) where {T<:AbstractFloat}
    return GaussianInfluenceWidth{T}(T(0.2))
end

function _resolve_deformable_mirror_influence(
    width::Real,
    ::Nothing,
    ::Nothing,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = try
        T(width)
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_influence,
        )
    end
    isfinite(converted) && converted > zero(T) ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_influence,
            "Gaussian influence width must be finite and positive",
        )
    return GaussianInfluenceWidth{T}(converted)
end

function _resolve_deformable_mirror_influence(
    ::Nothing,
    coupling::Real,
    ::Nothing,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = try
        T(coupling)
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_influence,
        )
    end
    isfinite(converted) && zero(T) < converted < one(T) ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_influence,
            "Gaussian mechanical coupling must be finite and lie in (0, 1)",
        )
    return GaussianMechanicalCoupling{T}(converted)
end

@inline function _resolve_deformable_mirror_influence(
    ::Nothing,
    ::Nothing,
    model::GaussianInfluenceWidth,
    ::Type{T},
) where {T<:AbstractFloat}
    width = try
        T(model.width)
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_influence,
        )
    end
    isfinite(width) && width > zero(T) ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_influence,
            "Gaussian influence width must be finite and positive",
        )
    return model
end

@inline function _resolve_deformable_mirror_influence(
    ::Nothing,
    ::Nothing,
    model::GaussianMechanicalCoupling,
    ::Type{T},
) where {T<:AbstractFloat}
    coupling = try
        T(model.coupling)
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_influence,
        )
    end
    isfinite(coupling) && zero(T) < coupling < one(T) ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_influence,
            "Gaussian mechanical coupling must be finite and lie in (0, 1)",
        )
    return model
end

@inline function _resolve_deformable_mirror_influence(
    ::Nothing,
    ::Nothing,
    model::AbstractDMInfluenceModel,
    ::Type{<:AbstractFloat},
)
    return model
end

function _resolve_deformable_mirror_influence(
    influence_width,
    mechanical_coupling,
    influence_model,
    ::Type{<:AbstractFloat},
)
    _deformable_mirror_definition_error(
        :invalid_deformable_mirror_influence,
        "specify at most one of influence_width, mechanical_coupling, or " *
        "influence_model; got $(typeof(influence_width)), " *
        "$(typeof(mechanical_coupling)), and $(typeof(influence_model))",
    )
end

@inline function _resolve_deformable_mirror_actuator_model(
    ::Nothing,
)
    return LinearStaticActuators()
end

function _resolve_deformable_mirror_actuator_model(
    model::AbstractDMActuatorModel,
)
    return try
        validate_dm_actuator_model(model)
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_actuator_model,
        )
    end
end

function _resolve_deformable_mirror_actuator_model(model)
    _deformable_mirror_definition_error(
        :invalid_deformable_mirror_actuator_model,
        "actuator_model must be an AbstractDMActuatorModel or nothing; got " *
        "$(typeof(model))",
    )
end

function _convert_deformable_mirror_misregistration(
    misregistration::Misregistration,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = try
        Misregistration(
            shift_x=T(misregistration.shift_x),
            shift_y=T(misregistration.shift_y),
            rotation_deg=T(rad2deg(misregistration.rotation_rad)),
            anamorphosis_angle=T(
                rad2deg(misregistration.anamorphosis_angle_rad),
            ),
            tangential_scaling=T(misregistration.tangential_scaling),
            radial_scaling=T(misregistration.radial_scaling),
            T=T,
        )
    catch error
        _throw_deformable_mirror_definition_failure(
            error,
            :invalid_deformable_mirror_misregistration,
        )
    end
    values = (
        converted.shift_x,
        converted.shift_y,
        converted.rotation_rad,
        converted.anamorphosis_angle_rad,
        converted.tangential_scaling,
        converted.radial_scaling,
        converted.transform...,
    )
    all(isfinite, values) &&
        converted.tangential_scaling > zero(T) &&
        converted.radial_scaling > zero(T) ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_misregistration,
            "actuator misregistration must be finite with positive scale " *
            "factors in numeric type $(T)",
        )
    return converted
end

function _convert_deformable_mirror_misregistration(
    misregistration,
    ::Type{<:AbstractFloat},
)
    _deformable_mirror_definition_error(
        :invalid_deformable_mirror_misregistration,
        "misregistration must be a Misregistration; got " *
        "$(typeof(misregistration))",
    )
end

@inline function _require_deformable_mirror_registration(
    registration::PupilRelayRegistration,
)
    return registration
end

function _require_deformable_mirror_registration(registration)
    _deformable_mirror_definition_error(
        :invalid_pupil_relay_registration,
        "pupil_relay_registration must be a PupilRelayRegistration; got " *
        "$(typeof(registration))",
    )
end

"""
    DeformableMirrorModel(; n_act=nothing, topology=nothing,
        influence_width=nothing, mechanical_coupling=nothing,
        influence_model=nothing, actuator_model=nothing,
        misregistration=Misregistration(T=T),
        pupil_relay_registration=PupilRelayRegistration(T=T),
        T=Float64)

Cold native deformable-mirror model for a `ControllableOpticDefinition`.
Exactly one of `n_act` or `topology` is required. Influence selection follows
the standalone `DeformableMirror` vocabulary, but no live command, surface, or
workspace storage is created until plant preparation.
"""
struct DeformableMirrorModel{
    T<:AbstractFloat,
    TP<:AbstractDMTopology,
    I<:AbstractDMInfluenceModel,
    A<:AbstractDMActuatorModel,
    M<:Misregistration{T},
    R<:PupilRelayRegistration,
}
    topology::TP
    influence_model::I
    actuator_model::A
    misregistration::M
    pupil_relay_registration::R
end

struct _DefaultDeformableMirrorMisregistration end
struct _DefaultDeformableMirrorPupilRelayRegistration end

const _DEFAULT_DEFORMABLE_MIRROR_MISREGISTRATION =
    _DefaultDeformableMirrorMisregistration()
const _DEFAULT_DEFORMABLE_MIRROR_PUPIL_RELAY_REGISTRATION =
    _DefaultDeformableMirrorPupilRelayRegistration()

@inline _resolve_deformable_mirror_misregistration(
    ::_DefaultDeformableMirrorMisregistration,
    ::Type{T},
) where {T<:AbstractFloat} = Misregistration(T=T)

@inline _resolve_deformable_mirror_misregistration(
    misregistration,
    ::Type{T},
) where {T<:AbstractFloat} =
    _convert_deformable_mirror_misregistration(misregistration, T)

@inline _resolve_deformable_mirror_registration(
    ::_DefaultDeformableMirrorPupilRelayRegistration,
    ::Type{T},
) where {T<:AbstractFloat} = PupilRelayRegistration(T=T)

@inline _resolve_deformable_mirror_registration(
    registration,
    ::Type{<:AbstractFloat},
) = _require_deformable_mirror_registration(registration)

function DeformableMirrorModel(;
    n_act=nothing,
    topology=nothing,
    influence_width=nothing,
    mechanical_coupling=nothing,
    influence_model=nothing,
    actuator_model=nothing,
    T::Type{<:AbstractFloat}=Float64,
    misregistration=_DEFAULT_DEFORMABLE_MIRROR_MISREGISTRATION,
    pupil_relay_registration=
        _DEFAULT_DEFORMABLE_MIRROR_PUPIL_RELAY_REGISTRATION,
)
    isconcretetype(T) || _deformable_mirror_definition_error(
        :invalid_deformable_mirror_numeric_type,
        "deformable-mirror numeric type must be concrete; got $(T)",
    )
    resolved_topology =
        _resolve_deformable_mirror_topology(n_act, topology, T)
    topology_command_count(resolved_topology) > 0 ||
        _deformable_mirror_definition_error(
            :invalid_deformable_mirror_topology,
            "deformable-mirror topology must contain at least one active " *
            "actuator",
        )
    resolved_influence = _resolve_deformable_mirror_influence(
        influence_width,
        mechanical_coupling,
        influence_model,
        T,
    )
    resolved_actuator =
        _resolve_deformable_mirror_actuator_model(actuator_model)
    resolved_misregistration =
        _resolve_deformable_mirror_misregistration(misregistration, T)
    resolved_registration = _resolve_deformable_mirror_registration(
        pupil_relay_registration, T)
    return DeformableMirrorModel{
        T,
        typeof(resolved_topology),
        typeof(resolved_influence),
        typeof(resolved_actuator),
        typeof(resolved_misregistration),
        typeof(resolved_registration),
    }(
        resolved_topology,
        resolved_influence,
        resolved_actuator,
        resolved_misregistration,
        resolved_registration,
    )
end

plant_model_definition_style(::Type{<:DeformableMirrorModel}) =
    ColdPlantModelDefinition()

@inline topology(model::DeformableMirrorModel) = model.topology
@inline influence_model(model::DeformableMirrorModel) =
    model.influence_model
@inline actuator_model(model::DeformableMirrorModel) = model.actuator_model
@inline n_actuators(model::DeformableMirrorModel) =
    topology_command_count(model.topology)

struct _PreparedPlantDeformableMirror{
    T<:AbstractFloat,
    P<:DeformableMirrorParams{T},
    M<:AbstractMatrix{T},
    X,
    Y,
    SM<:OpticalPlaneMetadata,
    R<:PupilRelayRegistration,
}
    endpoint::CommandEndpointID
    params::P
    modes::M
    separable_x::X
    separable_y_t::Y
    surface_metadata::SM
    pupil_relay_registration::R
end

mutable struct _PlantDeformableMirrorState{D<:DeformableMirror}
    active::D
end

mutable struct _PlantDeformableMirrorWorkspace{D<:DeformableMirror}
    staged::D
end

function _require_deformable_mirror_schema(
    model::DeformableMirrorModel{T},
    definition::ControllableOpticDefinition,
) where {T<:AbstractFloat}
    schemas = command_schemas(definition)
    length(schemas) == 1 || _deformable_mirror_preparation_error(
        :deformable_mirror_command_schema,
        "native deformable mirror $(controllable_optic_id(definition)) " *
        "requires exactly one actuator command endpoint",
    )
    schema = only(schemas)
    command_numeric_type(schema) === T ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_numeric_type,
            "native deformable mirror $(controllable_optic_id(definition)) " *
            "uses $(T), but endpoint $(command_endpoint_id(schema)) uses " *
            "$(command_numeric_type(schema))",
        )
    expected_dimensions = (n_actuators(model),)
    command_dimensions(schema) == expected_dimensions ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_dimensions,
            "native deformable mirror $(controllable_optic_id(definition)) " *
            "requires command dimensions $(expected_dimensions); got " *
            "$(command_dimensions(schema))",
        )
    command_units(schema) == CommandUnit(:metre) ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_units,
            "native deformable-mirror actuator commands must use metres",
        )
    command_sign_convention(schema) ==
        CommandSignConvention(:positive_surface_increases_opd) ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_sign_convention,
            "native deformable-mirror actuator commands must use " *
            ":positive_surface_increases_opd",
        )
    command_basis(schema).kind === :actuator ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_basis,
            "native deformable-mirror command basis kind must be :actuator",
        )
    return schema
end

@noinline function _throw_deformable_mirror_preparation_failure(
    error::AdaptiveOpticsSimError,
    id::ControllableOpticID,
)
    _deformable_mirror_preparation_error(
        :deformable_mirror_preparation,
        "failed to prepare native deformable mirror $id: " *
        _deformable_mirror_error_message(error),
    )
end

@noinline function _throw_deformable_mirror_preparation_failure(
    error::InterruptException,
    ::ControllableOpticID,
)
    throw(error)
end

@noinline function _throw_deformable_mirror_preparation_failure(
    error,
    id::ControllableOpticID,
)
    _deformable_mirror_preparation_error(
        :deformable_mirror_preparation,
        "failed to prepare native deformable mirror $id: " *
        _deformable_mirror_error_message(error),
    )
end

@inline function _snapshot_deformable_mirror_topology(
    topology::TP,
) where {TP<:AbstractDMTopology}
    # A topology is structurally immutable but owns mutable coordinate, mask,
    # index, and possibly metadata storage. Plant preparation must detach that
    # storage from the caller-owned cold model before publishing runtime params.
    return deepcopy(topology)::TP
end

function _prepare_deformable_mirror_controllable_optic(
    model::DeformableMirrorModel{T},
    definition::ControllableOpticDefinition,
    telescope::Telescope,
) where {T<:AbstractFloat}
    schema = _require_deformable_mirror_schema(model, definition)
    eltype(pupil_reflectivity(telescope)) === T ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_numeric_type,
            "native deformable mirror $(controllable_optic_id(definition)) " *
            "uses $(T), but the telescope pupil uses " *
            "$(eltype(pupil_reflectivity(telescope)))",
        )
    dm = try
        topology_snapshot =
            _snapshot_deformable_mirror_topology(model.topology)
        DeformableMirror(
            telescope;
            topology=topology_snapshot,
            influence_model=model.influence_model,
            actuator_model=model.actuator_model,
            T,
            misregistration=model.misregistration,
            backend=backend(telescope),
        )
    catch error
        _throw_deformable_mirror_preparation_failure(
            error,
            controllable_optic_id(definition),
        )
    end
    surface_metadata = OpticalPlaneMetadata(
        PupilPlane(),
        surface_opd(dm);
        coordinate_domain=MetricCoordinates(),
        sampling=ntuple(
            index -> T(telescope.aperture.sampling_m[index]),
            Val(2),
        ),
        origin=ntuple(
            index -> T(telescope.aperture.origin_m[index]),
            Val(2),
        ),
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
    )
    return _PreparedPlantDeformableMirror(
        command_endpoint_id(schema),
        dm.params,
        dm.state.modes,
        dm.state.separable_x,
        dm.state.separable_y_t,
        surface_metadata,
        model.pupil_relay_registration,
    )
end

@inline function prepare_controllable_optic(
    model::DeformableMirrorModel,
    definition::ControllableOpticDefinition,
    telescope::Telescope,
    ::AbstractAtmosphere,
)
    return _prepare_deformable_mirror_controllable_optic(
        model, definition, telescope)
end

@inline function prepare_target_local_controllable_optic(
    model::DeformableMirrorModel,
    definition::ControllableOpticDefinition,
    telescope::Telescope,
    ::AbstractTimedAtmosphereDefinition,
    ::AbstractComputeDevice,
)
    return _prepare_deformable_mirror_controllable_optic(
        model, definition, telescope)
end

@inline function _deformable_mirror_separable_runtime(
    ::Nothing,
    ::Nothing,
    actuator_commands,
    opd,
)
    return nothing, nothing
end

function _deformable_mirror_separable_runtime(
    separable_x::AbstractMatrix,
    separable_y_t::AbstractMatrix,
    actuator_commands,
    opd,
)
    command_grid = reshape(
        actuator_commands,
        size(separable_x, 2),
        size(separable_y_t, 1),
    )
    temporary = similar(
        opd,
        eltype(opd),
        size(separable_x, 1),
        size(separable_x, 2),
    )
    fill!(temporary, zero(eltype(temporary)))
    return command_grid, temporary
end

function _deformable_mirror_separable_runtime(
    separable_x,
    separable_y_t,
    ::Any,
    ::Any,
)
    _deformable_mirror_preparation_error(
        :deformable_mirror_preparation,
        "native deformable-mirror separable factors are incomplete: " *
        "$(typeof(separable_x)) and $(typeof(separable_y_t))",
    )
end

function _require_deformable_mirror_command_storage(
    prepared::_PreparedPlantDeformableMirror{T},
    command::AbstractVector{T},
) where {T<:AbstractFloat}
    length(command) == size(prepared.modes, 2) ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_dimensions,
            "native deformable-mirror command length $(length(command)) " *
            "does not match $(size(prepared.modes, 2)) actuators",
        )
    typeof(backend(command)) ===
        typeof(prepared.surface_metadata.backend) ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_backend,
            "native deformable-mirror command storage and sampled surface " *
            "use different array backends",
        )
    compute_device(command) == prepared.surface_metadata.device ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_command_device,
            "native deformable-mirror command storage and sampled surface " *
            "occupy different compute devices",
        )
    return command
end

function _require_deformable_mirror_command_storage(
    ::_PreparedPlantDeformableMirror{T},
    command,
) where {T<:AbstractFloat}
    _deformable_mirror_preparation_error(
        :deformable_mirror_command_storage,
        "native deformable-mirror commands must use AbstractVector{$T}; got " *
        "$(typeof(command))",
    )
end

function _allocate_deformable_mirror_runtime(
    prepared::_PreparedPlantDeformableMirror{T},
    command,
) where {T<:AbstractFloat}
    source = _require_deformable_mirror_command_storage(prepared, command)
    opd = allocate_array(
        prepared.surface_metadata.backend,
        T,
        prepared.surface_metadata.dimensions...,
    )
    fill!(opd, zero(T))
    opd_vec = reshape(opd, :)
    coefs = allocate_array(
        prepared.surface_metadata.backend,
        T,
        size(prepared.modes, 2),
    )
    actuator_coefs = similar(coefs)
    copyto!(coefs, source)
    fill!(actuator_coefs, zero(T))
    coefs_grid, separable_tmp =
        _deformable_mirror_separable_runtime(
            prepared.separable_x,
            prepared.separable_y_t,
            actuator_coefs,
            opd,
        )
    state = DeformableMirrorState{
        T,
        typeof(opd),
        typeof(prepared.modes),
        typeof(opd_vec),
        typeof(coefs_grid),
    }(
        opd,
        opd_vec,
        prepared.modes,
        coefs,
        actuator_coefs,
        coefs_grid,
        prepared.separable_x,
        prepared.separable_y_t,
        separable_tmp,
    )
    dm = DeformableMirror{
        typeof(prepared.params),
        typeof(state),
        typeof(prepared.surface_metadata.backend),
    }(prepared.params, state)
    update_surface!(dm)
    return dm
end

function prepare_controllable_optic_state(
    prepared::_PreparedPlantDeformableMirror,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    length(endpoint_ids) == 1 &&
        length(initial_commands) == 1 &&
        only(endpoint_ids) == prepared.endpoint ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_endpoint_binding,
            "native deformable-mirror endpoint binding changed after " *
            "preparation",
        )
    return _PlantDeformableMirrorState(
        _allocate_deformable_mirror_runtime(
            prepared,
            only(initial_commands),
        ),
    )
end

function prepare_controllable_optic_workspace(
    prepared::_PreparedPlantDeformableMirror{T},
) where {T<:AbstractFloat}
    command = allocate_array(
        prepared.surface_metadata.backend,
        T,
        size(prepared.modes, 2),
    )
    fill!(command, zero(T))
    return _PlantDeformableMirrorWorkspace(
        _allocate_deformable_mirror_runtime(prepared, command),
    )
end

function stage_controllable_optic_command!(
    prepared::_PreparedPlantDeformableMirror{T},
    ::_PlantDeformableMirrorState,
    workspace::_PlantDeformableMirrorWorkspace,
    endpoint::CommandEndpointID,
    effective_command::AbstractVector{T},
    ::PlantTimestamp,
) where {T<:AbstractFloat}
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application,
        :deformable_mirror_endpoint_mismatch,
        "native deformable mirror received command for $endpoint; expected " *
        "$(prepared.endpoint)",
    ))
    _require_deformable_mirror_command_storage(prepared, effective_command)
    set_command!(workspace.staged, effective_command)
    update_surface!(workspace.staged)
    return nothing
end

function stage_controllable_optic_command!(
    prepared::_PreparedPlantDeformableMirror{T},
    ::_PlantDeformableMirrorState,
    ::_PlantDeformableMirrorWorkspace,
    ::CommandEndpointID,
    effective_command,
    ::PlantTimestamp,
) where {T<:AbstractFloat}
    throw(PlantCommandError(
        :physical_application,
        :deformable_mirror_command_storage,
        "native deformable-mirror effective command must use " *
        "AbstractVector{$T}; got $(typeof(effective_command))",
    ))
end

@inline function commit_controllable_optic_command!(
    ::_PreparedPlantDeformableMirror,
    state::_PlantDeformableMirrorState{D},
    workspace::_PlantDeformableMirrorWorkspace{D},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {D<:DeformableMirror}
    previous = state.active
    state.active = workspace.staged
    workspace.staged = previous
    return nothing
end

function _deformable_mirror_surface_metadata(
    prepared::_PreparedPlantDeformableMirror{T},
    destination::PupilFunction,
) where {T<:AbstractFloat}
    eltype(destination.opd) === T ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_surface_numeric_type,
            "native deformable-mirror surface type $(T) does not match path " *
            "OPD type $(eltype(destination.opd))",
        )
    typeof(backend(destination.opd)) ===
        typeof(prepared.surface_metadata.backend) ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_surface_backend,
            "native deformable-mirror surface and path pupil use different " *
            "array backends",
        )
    compute_device(destination.opd) == prepared.surface_metadata.device ||
        _deformable_mirror_preparation_error(
            :deformable_mirror_surface_device,
            "native deformable-mirror surface and path pupil occupy different " *
            "compute devices",
        )
    return prepared.surface_metadata
end

function _deformable_mirror_surface_metadata(
    ::_PreparedPlantDeformableMirror,
    destination,
)
    _deformable_mirror_preparation_error(
        :unsupported_path_input,
        "native deformable-mirror path coupling requires PupilFunction input; " *
        "got $(typeof(destination))",
    )
end

function prepare_controllable_optic_path_coupling(
    prepared::_PreparedPlantDeformableMirror,
    definition::ControllableOpticDefinition,
    path,
)
    destination = path_input(path)
    metadata =
        _deformable_mirror_surface_metadata(prepared, destination)
    return prepare_sampled_pupil_footprint_coupling(
        metadata,
        destination.opd,
        path,
        controllable_optic_placement(definition);
        registration=prepared.pupil_relay_registration,
    )
end

@inline function apply_controllable_optic_surface!(
    input::PupilFunction,
    ::_PreparedPlantDeformableMirror,
    state::_PlantDeformableMirrorState,
    coupling::Union{
        PreparedIdentityPupilFootprintCoupling,
        PreparedPupilFootprintCoupling,
    },
)
    return apply_sampled_pupil_surface!(
        input,
        surface_opd(state.active),
        coupling,
    )
end
