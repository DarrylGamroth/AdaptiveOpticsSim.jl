# Focused effective-command publication and target-local owner qualification.

struct CommandPublicationTestTelescope <: AbstractTelescope end
struct CommandPublicationTestAtmosphereDefinition <:
    AbstractTimedAtmosphereDefinition end
struct CommandPublicationTestModel end
struct CommandPublicationForeignBackend <: AbstractArrayBackend end

Plant.plant_model_definition_style(::Type{CommandPublicationTestModel}) =
    ColdPlantModelDefinition()

struct CommandPublicationTestPrepared{D<:AbstractComputeDevice}
    endpoint::CommandEndpointID
    target::D
end

mutable struct CommandPublicationTestState{A<:AbstractVector}
    active::A
end

mutable struct CommandPublicationTestWorkspace{A<:AbstractVector}
    staging::A
end

function Plant.prepare_target_local_controllable_optic(
    ::CommandPublicationTestModel,
    definition::ControllableOpticDefinition,
    ::CommandPublicationTestTelescope,
    ::CommandPublicationTestAtmosphereDefinition,
    target::AbstractComputeDevice,
)
    return CommandPublicationTestPrepared(
        first(command_endpoint_ids(definition)), target)
end

function Plant.prepare_controllable_optic_state(
    prepared::CommandPublicationTestPrepared,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    return CommandPublicationTestState(first(initial_commands))
end

function Plant.prepare_controllable_optic_workspace(
    prepared::CommandPublicationTestPrepared,
)
    staging = allocate_device_array(prepared.target, Float64, 2)
    fill!(staging, 0.0)
    return CommandPublicationTestWorkspace(staging)
end

function Plant.stage_controllable_optic_command!(
    prepared::CommandPublicationTestPrepared,
    ::CommandPublicationTestState,
    workspace::CommandPublicationTestWorkspace,
    endpoint::CommandEndpointID,
    command::AbstractVector{Float64},
    ::PlantTimestamp,
)
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "publication-test optic received another endpoint"))
    copyto!(workspace.staging, command)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::CommandPublicationTestPrepared,
    state::CommandPublicationTestState{A},
    workspace::CommandPublicationTestWorkspace{A},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {A}
    previous = state.active
    state.active = workspace.staging
    workspace.staging = previous
    return nothing
end

function Plant.validate_controllable_optic_target(
    prepared::CommandPublicationTestPrepared,
    target::AbstractComputeDevice,
)
    prepared.target == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "publication-test preparation occupies another target"))
    return prepared
end

function Plant.validate_controllable_optic_state_target(
    ::CommandPublicationTestPrepared,
    state::CommandPublicationTestState,
    target::AbstractComputeDevice,
)
    compute_device(state.active) == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "publication-test state occupies another target"))
    return state
end

function Plant.validate_controllable_optic_workspace_target(
    ::CommandPublicationTestPrepared,
    workspace::CommandPublicationTestWorkspace,
    target::AbstractComputeDevice,
)
    compute_device(workspace.staging) == target || throw(
        PlantPreparationError(:controllable_optic, :wrong_device,
            "publication-test workspace occupies another target"))
    return workspace
end

Plant.structural_resource_fact(
    ::CommandPublicationTestPrepared,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(id, target, 0, 0)

Plant.structural_resource_fact(
    state::CommandPublicationTestState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(
    id, target, structural_array_bytes(state.active, target), 0)

Plant.structural_resource_fact(
    workspace::CommandPublicationTestWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(
    id, target, 0, structural_array_bytes(workspace.staging, target))

struct SharedWorkspacePublicationTestModel end

Plant.plant_model_definition_style(
    ::Type{SharedWorkspacePublicationTestModel}) = ColdPlantModelDefinition()

struct SharedWorkspacePublicationTestPrepared{E,D<:AbstractComputeDevice}
    endpoints::E
    target::D
end

mutable struct SharedWorkspacePublicationTestState{A<:AbstractVector}
    active::A
end

mutable struct SharedWorkspacePublicationTestWorkspace{A<:AbstractVector}
    staging::A
end

function Plant.prepare_target_local_controllable_optic(
    ::SharedWorkspacePublicationTestModel,
    definition::ControllableOpticDefinition,
    ::CommandPublicationTestTelescope,
    ::CommandPublicationTestAtmosphereDefinition,
    target::AbstractComputeDevice,
)
    return SharedWorkspacePublicationTestPrepared(
        Tuple(command_endpoint_ids(definition)), target)
end

function Plant.prepare_controllable_optic_state(
    ::SharedWorkspacePublicationTestPrepared,
    ::ControllableOpticDefinition,
    ::Tuple,
    initial_commands::Tuple,
)
    return SharedWorkspacePublicationTestState(first(initial_commands))
end


function Plant.prepare_controllable_optic_workspace(
    prepared::SharedWorkspacePublicationTestPrepared,
)
    staging = allocate_device_array(prepared.target, Float64, 2)
    fill!(staging, 0.0)
    return SharedWorkspacePublicationTestWorkspace(staging)
end

function Plant.stage_controllable_optic_command!(
    prepared::SharedWorkspacePublicationTestPrepared,
    ::SharedWorkspacePublicationTestState,
    workspace::SharedWorkspacePublicationTestWorkspace,
    endpoint::CommandEndpointID,
    command::AbstractVector{Float64},
    ::PlantTimestamp,
)
    endpoint in prepared.endpoints || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "shared-workspace test optic received another endpoint"))
    copyto!(workspace.staging, command)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::SharedWorkspacePublicationTestPrepared,
    state::SharedWorkspacePublicationTestState{A},
    workspace::SharedWorkspacePublicationTestWorkspace{A},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {A}
    previous = state.active
    state.active = workspace.staging
    workspace.staging = previous
    return nothing
end

function Plant.validate_controllable_optic_target(
    prepared::SharedWorkspacePublicationTestPrepared,
    target::AbstractComputeDevice,
)
    prepared.target == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "shared-workspace test preparation occupies another target"))
    return prepared
end

function Plant.validate_controllable_optic_state_target(
    ::SharedWorkspacePublicationTestPrepared,
    state::SharedWorkspacePublicationTestState,
    target::AbstractComputeDevice,
)
    compute_device(state.active) == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "shared-workspace test state occupies another target"))
    return state
end

function Plant.validate_controllable_optic_workspace_target(
    ::SharedWorkspacePublicationTestPrepared,
    workspace::SharedWorkspacePublicationTestWorkspace,
    target::AbstractComputeDevice,
)
    compute_device(workspace.staging) == target || throw(
        PlantPreparationError(:controllable_optic, :wrong_device,
            "shared-workspace test workspace occupies another target"))
    return workspace
end

Plant.structural_resource_fact(
    ::SharedWorkspacePublicationTestPrepared,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(id, target, 0, 0)

Plant.structural_resource_fact(
    state::SharedWorkspacePublicationTestState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(
    id, target, structural_array_bytes(state.active, target), 0)

Plant.structural_resource_fact(
    workspace::SharedWorkspacePublicationTestWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
) = KnownStructuralResourceFact(
    id, target, 0, structural_array_bytes(workspace.staging, target))

function command_publication_test_schema(
    ::Type{T}=Float64;
    dimensions=(2,),
    id=:publication_schema,
    version=1,
    endpoint=:publication_endpoint,
    units=:metre,
    sign_convention=:positive_surface_increases_opd,
    basis=CommandBasis(:actuator, :publication_basis),
    basis_revision=1,
    semantics=AbsoluteCommand,
    bounds=UniformCommandBounds(T(-10), T(10)),
    value_policy=CommandValuePolicy(range_stage=EnforceOnApplication),
    sequence_policy=CommandSequencePolicy(),
    effective_time_policy=CommandEffectiveTimePolicy(),
    silence_policy=CommandSilencePolicy(),
) where {T}
    return PlantCommandSchema(
        T,
        dimensions;
        id,
        version,
        endpoint,
        units,
        sign_convention,
        basis,
        basis_revision,
        semantics,
        bounds,
        value_policy,
        sequence_policy,
        effective_time_policy,
        silence_policy,
    )
end

function command_publication_test_definition(schemas)
    return ControllableOpticDefinition(
        :publication_optic,
        CommandPublicationTestModel(),
        schemas;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
end

function command_publication_test_owner(;
    schemas=(command_publication_test_schema(),),
    configurations=nothing,
)
    definition = command_publication_test_definition(schemas)
    resolved_configurations = configurations === nothing ? Tuple(
        CommandEndpointConfiguration(
            command_endpoint_id(schema), fill(0.0, command_dimensions(schema));
            capacity=2,
        ) for schema in schemas) : configurations
    authority = CommandAuthorityIdentity(
        Plant._COMMAND_AUTHORITY_IDENTITY_TOKEN)
    owner = Plant._prepare_target_local_controllable_optic_owner(
        definition,
        CommandPublicationTestTelescope(),
        CommandPublicationTestAtmosphereDefinition(),
        HostComputeDevice(),
        Plant._sorted_command_endpoint_configurations(
            resolved_configurations),
        authority,
    )
    return (; owner, authority, definition, configurations=resolved_configurations)
end

function command_publication_test_publication(
    fixture;
    schema=only(command_schemas(fixture.definition)),
    target=HostComputeDevice(),
    optic=controllable_optic_id(fixture.definition),
    timestamp=PlantTimestamp(1),
    sequence=1,
    authority=fixture.authority,
)
    return Plant._effective_command_publication(
        authority, target, optic, schema, timestamp, sequence)
end

function command_publication_test_error(f, reason::Symbol)
    error = try
        f()
        nothing
    catch caught
        caught
    end
    @test error isa PlantCommandError
    if error isa PlantCommandError
        @test error.stage === :effective_command_publication
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function command_publication_test_fingerprint(fixture)
    endpoint = CommandEndpointID(:publication_endpoint)
    state = target_local_command_endpoint_state(fixture.owner, endpoint)
    physical = target_local_controllable_optic_state(
        fixture.owner).physical
    return (
        effective=copy(effective_command(state)),
        physical=copy(physical.active),
        has_publication=has_effective_command_publication(state),
        timestamp=last_effective_command_publication_timestamp(state),
        sequence=last_effective_command_publication_sequence(state),
    )
end

function command_publication_test_stage_bytes(
    endpoint_owner,
    publication,
    payload,
)
    return @allocated Plant._stage_effective_command_publication!(
        endpoint_owner, publication, payload)
end

function command_publication_test_commit_bytes(endpoint_owner)
    return @allocated Plant._commit_effective_command_publication!(
        endpoint_owner)
end

@testset "Publication identities and sequences are internally sealed" begin
    fixture = command_publication_test_owner()
    schema = only(command_schemas(fixture.definition))
    publication = command_publication_test_publication(fixture)
    endpoint = only(target_local_command_endpoints(
        prepared_target_local_controllable_optic(fixture.owner)))
    endpoint_state = target_local_command_endpoint_state(
        fixture.owner, command_endpoint_id(endpoint))

    @test_throws ArgumentError CommandAuthorityIdentity(
        Plant._CommandAuthorityIdentityToken())
    @test_throws ArgumentError EffectiveCommandPublication(
        Plant._EffectiveCommandPublicationToken(),
        fixture.authority,
        HostComputeDevice(),
        controllable_optic_id(fixture.definition),
        schema,
        PlantTimestamp(1),
        EffectiveCommandPublicationSequence(1),
    )
    @test_throws ArgumentError EffectiveCommandPublicationSequence(
        UInt64(1), Plant._EffectiveCommandPublicationSequenceToken())
    @test_throws ArgumentError TargetLocalCommandEndpointOwner(
        Plant._TargetLocalCommandEndpointOwnerToken(),
        fixture.owner,
        endpoint,
        endpoint_state,
        UInt32(1),
    )
    @test publication isa EffectiveCommandPublication
    for invalid in (0, -1, true, 1.0)
        @test_throws PlantCommandError EffectiveCommandPublicationSequence(
            invalid)
    end
    @test_throws PlantCommandError EffectiveCommandPublicationSequence(
        big(typemax(UInt64)) + 1)
end

@testset "Scalar target-local effective values remain scalar" begin
    schema = command_publication_test_schema(dimensions=())
    target = HostComputeDevice()
    values = Plant._prepare_target_local_effective_values(
        schema, 1.5, target)
    @test Plant._active_effective_command(values) === 1.5
    @test Plant._staged_effective_command(values) === 1.5
    Plant._copy_staged_effective_command!(values, 2.5)
    Plant._commit_staged_effective_command!(values)
    @test Plant._active_effective_command(values) === 2.5
    id = StructuralResourceOwnerID(
        :effective_command_replica, :scalar_publication_endpoint)
    binding = Plant._TargetLocalCommandEndpointBinding()
    endpoint = PreparedTargetLocalCommandEndpoint(
        binding,
        CommandAuthorityIdentity(Plant._COMMAND_AUTHORITY_IDENTITY_TOKEN),
        target,
        ControllableOpticID(:scalar_publication_optic),
        schema,
    )
    origin = zero(PlantTimestamp)
    state = TargetLocalCommandEndpointState(
        binding, values, origin, origin, UInt64(0), UInt64(0), false, false)
    @test validate_target_local_command_endpoint_target(
        endpoint, state, target) === state
    fact = structural_resource_fact(state, id, target)
    @test structural_resource_known(fact)
    @test structural_resident_bytes(fact) == 0
    @test structural_workspace_bytes(fact) == 0
end

@testset "Effective-command publication values are target-local and sealed" begin
    caller_initial = [1.0, 2.0]
    schema = command_publication_test_schema()
    configurations = (
        CommandEndpointConfiguration(
            :publication_endpoint, caller_initial; capacity=2),)
    fixture = command_publication_test_owner(; schemas=(schema,), configurations)
    state = target_local_command_endpoint_state(
        fixture.owner, CommandEndpointID(:publication_endpoint))
    values = getfield(state, :values)
    physical = target_local_controllable_optic_state(fixture.owner).physical

    caller_initial .= 9.0
    @test effective_command(state) == [1.0, 2.0]
    @test physical.active == [1.0, 2.0]
    @test effective_command(state) !== caller_initial
    @test effective_command(state) !== getfield(values, :staging)
    @test effective_command(state) !== physical.active

    second = command_publication_test_owner()
    second_state = target_local_command_endpoint_state(
        second.owner, CommandEndpointID(:publication_endpoint))
    @test effective_command(second_state) !== effective_command(state)

    payload = [3.0, 4.0]
    publication = command_publication_test_publication(fixture)
    Plant._stage_effective_command_publication!(
        fixture.owner, publication, payload)
    payload .= 8.0
    Plant._commit_effective_command_publication!(
        fixture.owner, CommandEndpointID(:publication_endpoint))
    @test effective_command(state) == [3.0, 4.0]
    @test physical.active == [3.0, 4.0]
    @test effective_command(state) !== payload
    @test command_authority_identity(publication) === fixture.authority
    @test command_endpoint_id(publication) ==
        CommandEndpointID(:publication_endpoint)
    @test controllable_optic_id(publication) ==
        ControllableOpticID(:publication_optic)
    @test effective_command_publication_timestamp(publication) ==
        PlantTimestamp(1)
    @test effective_command_publication_sequence(publication) ==
        EffectiveCommandPublicationSequence(1)
end

@testset "Already-effective incremental publications are not reapplied" begin
    schema = command_publication_test_schema(;
        semantics=IncrementalCommand)
    configurations = (
        CommandEndpointConfiguration(
            :publication_endpoint, [1.0, 1.0]; capacity=2),)
    fixture = command_publication_test_owner(; schemas=(schema,), configurations)
    publication = command_publication_test_publication(fixture)
    Plant._stage_effective_command_publication!(
        fixture.owner, publication, [2.0, 3.0])
    Plant._commit_effective_command_publication!(
        fixture.owner, CommandEndpointID(:publication_endpoint))
    endpoint_state = target_local_command_endpoint_state(
        fixture.owner, CommandEndpointID(:publication_endpoint))
    physical_state = target_local_controllable_optic_state(
        fixture.owner).physical
    @test effective_command(endpoint_state) == [2.0, 3.0]
    @test physical_state.active == [2.0, 3.0]
end

@testset "An optic has one physical staging owner" begin
    alpha = command_publication_test_schema(
        id=:shared_alpha_schema,
        endpoint=:shared_alpha,
        basis=CommandBasis(:actuator, :shared_alpha_basis),
    )
    beta = command_publication_test_schema(
        id=:shared_beta_schema,
        endpoint=:shared_beta,
        basis=CommandBasis(:actuator, :shared_beta_basis),
    )
    definition = ControllableOpticDefinition(
        :shared_workspace_optic,
        SharedWorkspacePublicationTestModel(),
        (beta, alpha);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    configurations = Plant._sorted_command_endpoint_configurations((
        CommandEndpointConfiguration(:shared_beta, [0.0, 0.0]; capacity=2),
        CommandEndpointConfiguration(:shared_alpha, [0.0, 0.0]; capacity=2),
    ))
    authority = CommandAuthorityIdentity(
        Plant._COMMAND_AUTHORITY_IDENTITY_TOKEN)
    owner = Plant._prepare_target_local_controllable_optic_owner(
        definition,
        CommandPublicationTestTelescope(),
        CommandPublicationTestAtmosphereDefinition(),
        HostComputeDevice(),
        configurations,
        authority,
    )
    fixture = (; owner, authority, definition)
    alpha_publication = command_publication_test_publication(
        fixture; schema=alpha, optic=ControllableOpticID(:shared_workspace_optic))
    beta_publication = command_publication_test_publication(
        fixture; schema=beta, optic=ControllableOpticID(:shared_workspace_optic))
    alpha_owner = target_local_command_endpoint_owner(
        owner, CommandEndpointID(:shared_alpha))
    beta_owner = target_local_command_endpoint_owner(
        owner, CommandEndpointID(:shared_beta))

    Plant._stage_effective_command_publication!(
        alpha_owner, alpha_publication, [1.0, 1.0])
    command_publication_test_error(:stage_pending) do
        Plant._stage_effective_command_publication!(
            alpha_owner, alpha_publication, [3.0, 3.0])
    end
    command_publication_test_error(:optic_stage_pending) do
        Plant._stage_effective_command_publication!(
            beta_owner, beta_publication, [2.0, 2.0])
    end
    command_publication_test_error(:no_staged_publication) do
        Plant._commit_effective_command_publication!(beta_owner)
    end
    Plant._commit_effective_command_publication!(alpha_owner)
    physical = target_local_controllable_optic_state(owner).physical
    alpha_state = target_local_command_endpoint_state(
        owner, CommandEndpointID(:shared_alpha))
    beta_state = target_local_command_endpoint_state(
        owner, CommandEndpointID(:shared_beta))
    @test effective_command(alpha_state) !== effective_command(beta_state)
    @test physical.active == [1.0, 1.0]
    @test effective_command(alpha_state) == [1.0, 1.0]
    @test effective_command(beta_state) == [0.0, 0.0]

    Plant._stage_effective_command_publication!(
        beta_owner, beta_publication, [2.0, 2.0])
    Plant._commit_effective_command_publication!(beta_owner)
    @test physical.active == [2.0, 2.0]
    @test effective_command(alpha_state) == [1.0, 1.0]
    @test effective_command(beta_state) == [2.0, 2.0]
end

@testset "Publication validation precedes active mutation" begin
    fixture = command_publication_test_owner()
    endpoint = CommandEndpointID(:publication_endpoint)
    first = command_publication_test_publication(
        fixture; timestamp=PlantTimestamp(10), sequence=2)
    Plant._stage_effective_command_publication!(fixture.owner, first, [1.0, 2.0])
    Plant._commit_effective_command_publication!(fixture.owner, endpoint)
    baseline = command_publication_test_fingerprint(fixture)
    foreign_authority = CommandAuthorityIdentity(
        Plant._COMMAND_AUTHORITY_IDENTITY_TOKEN)
    foreign_target = AcceleratorComputeDevice(
        CommandPublicationForeignBackend(), UInt32(1))

    cases = (
        (:foreign_authority, command_publication_test_publication(
            fixture; authority=foreign_authority, sequence=3), [2.0, 3.0]),
        (:wrong_target, command_publication_test_publication(
            fixture; target=foreign_target, sequence=3), [2.0, 3.0]),
        (:foreign_optic, command_publication_test_publication(
            fixture; optic=ControllableOpticID(:foreign_optic), sequence=3),
            [2.0, 3.0]),
        (:foreign_endpoint, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                endpoint=:foreign_endpoint), sequence=3), [2.0, 3.0]),
        (:schema_id, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                id=:foreign_schema), sequence=3), [2.0, 3.0]),
        (:schema_version, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                version=2), sequence=3), [2.0, 3.0]),
        (:numeric_type, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(Float32),
            sequence=3), [2.0, 3.0]),
        (:dimensions, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                dimensions=(3,)), sequence=3), [2.0, 3.0]),
        (:units, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                units=:volt), sequence=3), [2.0, 3.0]),
        (:sign_convention, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                sign_convention=:opposite_sign), sequence=3), [2.0, 3.0]),
        (:basis, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                basis=CommandBasis(:modal, :foreign_basis)), sequence=3),
            [2.0, 3.0]),
        (:basis_revision, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                basis_revision=2), sequence=3), [2.0, 3.0]),
        (:semantics, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                semantics=IncrementalCommand), sequence=3), [2.0, 3.0]),
        (:bounds, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                bounds=UniformCommandBounds(-20.0, 20.0)), sequence=3),
            [2.0, 3.0]),
        (:value_policy, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                value_policy=CommandValuePolicy(
                    nonfinite=FailOnInvalidCommand,
                    range_stage=EnforceOnApplication)), sequence=3),
            [2.0, 3.0]),
        (:sequence_policy, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                sequence_policy=CommandSequencePolicy(
                    duplicate=FailOnSequence)), sequence=3), [2.0, 3.0]),
        (:effective_time_policy, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                effective_time_policy=CommandEffectiveTimePolicy(
                    future=RejectFutureCommand)), sequence=3), [2.0, 3.0]),
        (:silence_policy, command_publication_test_publication(
            fixture; schema=command_publication_test_schema(
                silence_policy=CommandSilencePolicy(
                    ApplySafeCommand,
                    AgeFromApplication;
                    timeout=PlantDuration(1),
                )), sequence=3), [2.0, 3.0]),
        (:dimensions, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=3), [2.0]),
        (:numeric_type, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=3),
            Float64[2.0 3.0]),
        (:nonfinite, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=3), [NaN, 3.0]),
        (:out_of_range, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=3), [11.0, 3.0]),
        (:duplicate_sequence, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=2), [2.0, 3.0]),
        (:stale_sequence, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(10), sequence=1), [2.0, 3.0]),
        (:regressing_timestamp, command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(9), sequence=3), [2.0, 3.0]),
    )

    for (reason, publication, payload) in cases
        command_publication_test_error(reason) do
            Plant._stage_effective_command_publication!(
                fixture.owner, publication, payload)
        end
        @test command_publication_test_fingerprint(fixture) == baseline
    end
end

@testset "Unequal-rate consumers hold canonical publications" begin
    fixture = command_publication_test_owner()
    endpoint = CommandEndpointID(:publication_endpoint)
    inputs = [
        (PlantTimestamp(20), 3, [3.0, 3.0]),
        (PlantTimestamp(10), 2, [2.0, 2.0]),
        (PlantTimestamp(10), 1, [1.0, 1.0]),
    ]
    sort!(inputs; by=input -> (plant_nanoseconds(input[1]), input[2]))
    cursor = 1
    observed = Vector{Vector{Float64}}()
    for sample_timestamp in PlantTimestamp.(Int64[5, 10, 15, 20])
        while cursor <= length(inputs) && inputs[cursor][1] <= sample_timestamp
            timestamp, sequence, payload = inputs[cursor]
            publication = command_publication_test_publication(
                fixture; timestamp, sequence)
            Plant._stage_effective_command_publication!(
                fixture.owner, publication, payload)
            Plant._commit_effective_command_publication!(fixture.owner, endpoint)
            cursor += 1
        end
        state = target_local_command_endpoint_state(fixture.owner, endpoint)
        push!(observed, copy(effective_command(state)))
    end
    @test observed == [[0.0, 0.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0]]

    alpha = command_publication_test_schema(endpoint=:alpha_endpoint)
    zeta = command_publication_test_schema(
        id=:zeta_schema, endpoint=:zeta_endpoint,
        basis=CommandBasis(:actuator, :zeta_basis))
    configurations = (
        CommandEndpointConfiguration(:zeta_endpoint, [4.0, 5.0]; capacity=2),
        CommandEndpointConfiguration(:alpha_endpoint, [1.0, 2.0]; capacity=2),
    )
    reordered = command_publication_test_owner(
        ; schemas=(zeta, alpha), configurations)
    prepared = prepared_target_local_controllable_optic(reordered.owner)
    @test Tuple(command_endpoint_id(endpoint) for endpoint in
        target_local_command_endpoints(prepared)) ==
        (CommandEndpointID(:alpha_endpoint), CommandEndpointID(:zeta_endpoint))
end

@testset "Publication stage and commit are fixed-capacity" begin
    fixture = command_publication_test_owner()
    endpoint = CommandEndpointID(:publication_endpoint)
    endpoint_owner = target_local_command_endpoint_owner(
        fixture.owner, endpoint)
    @test isconcretetype(typeof(endpoint_owner))
    @test all(isconcretetype, fieldtypes(typeof(endpoint_owner)))
    warm = command_publication_test_publication(
        fixture; timestamp=PlantTimestamp(1), sequence=1)
    @test @inferred(Plant._stage_effective_command_publication!(
        endpoint_owner, warm, [1.0, 2.0])) === nothing
    @test @inferred(Plant._commit_effective_command_publication!(
        endpoint_owner)) === nothing

    measured = command_publication_test_publication(
        fixture; timestamp=PlantTimestamp(2), sequence=2)
    payload = [2.0, 3.0]
    Plant._stage_effective_command_publication!(
        endpoint_owner, measured, payload)
    Plant._commit_effective_command_publication!(endpoint_owner)
    if !coverage_instrumented()
        next_publication = command_publication_test_publication(
            fixture; timestamp=PlantTimestamp(3), sequence=3)
        stage_bytes = command_publication_test_stage_bytes(
            endpoint_owner, next_publication, payload)
        commit_bytes = command_publication_test_commit_bytes(endpoint_owner)
        @test stage_bytes == 0
        @test commit_bytes == 0
    end

    state = target_local_command_endpoint_state(fixture.owner, endpoint)
    id = StructuralResourceOwnerID(
        :effective_command_replica, :publication_endpoint)
    fact = structural_resource_fact(state, id, HostComputeDevice())
    @test structural_resource_known(fact)
    @test structural_resident_bytes(fact) == UInt64(2 * sizeof(Float64))
    @test structural_workspace_bytes(fact) == UInt64(2 * sizeof(Float64))
end
