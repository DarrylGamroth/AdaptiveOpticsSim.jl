struct CommandSchemaTestAtmosphere <: AdaptiveOpticsSim.AbstractAtmosphere end

struct CommandSchemaTestOpticModel
    name::Symbol
end

AdaptiveOpticsSim.plant_model_definition_style(
    ::Type{CommandSchemaTestOpticModel},
) = ColdPlantModelDefinition()

function _captured_command_schema_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function _assert_command_error(f, stage::Symbol, reason::Symbol)
    error = _captured_command_schema_error(f)
    @test error isa PlantCommandError
    if error isa PlantCommandError
        @test error.stage === stage
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function _schema_fixture(::Type{T}=Float32;
    id=:woofer_schema,
    version=1,
    endpoint=:woofer_command,
    dimensions=(3,),
    units=:metre,
    sign_convention=:positive_surface_increases_opd,
    basis=CommandBasis(:actuator, :woofer_actuators),
    basis_revision=1,
    semantics=AbsoluteCommand,
    bounds=UniformCommandBounds(T(-1), T(1)),
    value_policy=CommandValuePolicy(),
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

@testset "Immutable versioned plant-command schemas" begin
    for name in (
        :PlantCommandError,
        :PlantCommandSchemaID,
        :PlantCommandSchemaVersion,
        :CommandBasisRevision,
        :CommandUnit,
        :CommandSignConvention,
        :CommandBasis,
        :CommandValueSemantics,
        :AbsoluteCommand,
        :IncrementalCommand,
        :UnboundedCommandValues,
        :UniformCommandBounds,
        :CommandValuePolicy,
        :CommandSequencePolicy,
        :CommandEffectiveTimePolicy,
        :CommandSilencePolicy,
        :PlantCommandSchema,
        :command_schema_id,
        :command_schema_version,
        :command_endpoint_id,
        :command_numeric_type,
        :command_dimensions,
        :command_schemas,
        :command_schema,
        :plant_command_schema,
        :validate_plant_command_payload,
    )
        @test Base.isexported(AdaptiveOpticsSim, name)
    end

    value_policy = CommandValuePolicy(
        RejectInvalidCommand,
        ClipInvalidCommand,
        ValidateOnPresentation,
    )
    sequence_policy = CommandSequencePolicy(
        RejectSequence,
        FailOnSequence,
        AcceptSequence,
        AcceptSequence,
    )
    effective_time_policy = CommandEffectiveTimePolicy(
        AllowFutureCommand,
        ApplyLateCommandNow,
        SupersedeOlderPendingCommands,
    )
    silence_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(1_000_000),
    )
    schema = @inferred _schema_fixture(
        Float32;
        id=PlantCommandSchemaID(:woofer_schema),
        version=PlantCommandSchemaVersion(3),
        endpoint=CommandEndpointID(:woofer_command),
        dimensions=(4,),
        units=CommandUnit(:metre),
        sign_convention=CommandSignConvention(
            :positive_surface_increases_opd,
        ),
        basis=CommandBasis(:actuator, :woofer_actuators),
        basis_revision=CommandBasisRevision(7),
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(-2f0, 2f0),
        value_policy,
        sequence_policy,
        effective_time_policy,
        silence_policy,
    )

    @test command_schema_id(schema) == PlantCommandSchemaID(:woofer_schema)
    @test command_schema_version(schema) == PlantCommandSchemaVersion(3)
    @test command_endpoint_id(schema) == CommandEndpointID(:woofer_command)
    @test command_numeric_type(schema) === Float32
    @test command_dimensions(schema) === (4,)
    @test command_units(schema) == CommandUnit(:metre)
    @test command_sign_convention(schema) == CommandSignConvention(
        :positive_surface_increases_opd,
    )
    @test command_basis(schema) == CommandBasis(
        :actuator,
        :woofer_actuators,
    )
    @test command_basis_revision(schema) == CommandBasisRevision(7)
    @test command_semantics(schema) == AbsoluteCommand
    @test command_bounds(schema) == UniformCommandBounds(-2f0, 2f0)
    @test command_value_policy(schema) === value_policy
    @test command_sequence_policy(schema) === sequence_policy
    @test command_effective_time_policy(schema) === effective_time_policy
    @test command_silence_policy(schema) === silence_policy
    @test !Base.ismutable(schema)
    @test isbitstype(typeof(value_policy))
    @test isbitstype(typeof(sequence_policy))
    @test isbitstype(typeof(effective_time_policy))
    @test isbitstype(typeof(silence_policy))

    @test PlantCommandSchemaID(:woofer_schema) ==
        PlantCommandSchemaID(:woofer_schema)
    @test PlantCommandSchemaID(:woofer_schema) !=
        PlantCommandSchemaID(:tweeter_schema)
    @test length(Set((PlantCommandSchemaID(:same),
        PlantCommandSchemaID(:same)))) == 1
    @test length(Set((PlantCommandSchemaID(:same),
        CommandEndpointID(:same)))) == 2
    @test sprint(show, PlantCommandSchemaID(:woofer_schema)) ==
        "PlantCommandSchemaID(:woofer_schema)"
    @test sprint(show, PlantCommandSchemaVersion(3)) ==
        "PlantCommandSchemaVersion(3)"
    @test sprint(show, CommandBasisRevision(7)) ==
        "CommandBasisRevision(7)"
    @test sprint(show, CommandBasis(:modal, :kl_modes)) ==
        "CommandBasis(:modal, :kl_modes)"
    @test all(!Base.ismutable(value) for value in (
        PlantCommandSchemaID(:schema),
        PlantCommandSchemaVersion(1),
        CommandBasisRevision(1),
        CommandUnit(:meter),
        CommandSignConvention(:positive_opd),
        CommandBasis(:modal, :kl_modes),
        value_policy,
        sequence_policy,
        effective_time_policy,
        silence_policy,
    ))

    for absent in (
        :state,
        :workspace,
        :sequence,
        :calendar,
        :applied_value,
        :session,
        :source_clock,
        :mapping,
        :lease,
        :queue,
        :transport,
        :outcome,
    )
        @test !hasproperty(schema, absent)
    end

    scalar = @inferred _schema_fixture(
        Float64;
        id=:focus_schema,
        endpoint=:focus_command,
        dimensions=(),
        units=:metre,
        sign_convention=:positive_focus_increases_opd,
        basis=CommandBasis(:rigid_body, :focus),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    @test @inferred(validate_plant_command_payload(scalar, 0.25)) === 0.25
    _assert_command_error(
        () -> validate_plant_command_payload(scalar, 0.25f0),
        :payload,
        :element_type,
    )
    _assert_command_error(
        () -> validate_plant_command_payload(scalar, [0.25]),
        :payload,
        :shape,
    )

    payload = Float32[-2, -1, 0, 1]
    @test @inferred(validate_plant_command_payload(schema, payload)) === payload
    matrix = reshape(Float32[-2, -1, 0, 1, 0, 0, 0, 0], 4, 2)
    column = @view matrix[:, 1]
    @test validate_plant_command_payload(schema, column) === column
    @test matrix[:, 1] == Float32[-2, -1, 0, 1]

    _assert_command_error(
        () -> validate_plant_command_payload(schema, Float64[0, 0, 0, 0]),
        :payload,
        :element_type,
    )
    _assert_command_error(
        () -> validate_plant_command_payload(schema, Float32[0, 0, 0]),
        :payload,
        :shape,
    )
    _assert_command_error(
        () -> validate_plant_command_payload(schema, "not numeric"),
        :payload,
        :payload_type,
    )
    _assert_command_error(
        () -> validate_plant_command_payload(schema,
            Float32[0, NaN, 0, 0]),
        :payload,
        :nonfinite_rejected,
    )

    rejected = _schema_fixture(
        Float32;
        bounds=UniformCommandBounds(-1f0, 1f0),
        value_policy=CommandValuePolicy(
            RejectInvalidCommand,
            RejectInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    _assert_command_error(
        () -> validate_plant_command_payload(rejected,
            Float32[0, 2, 0]),
        :payload,
        :out_of_range_rejected,
    )

    failed = _schema_fixture(
        Float32;
        bounds=UniformCommandBounds(-1f0, 1f0),
        value_policy=CommandValuePolicy(
            FailOnInvalidCommand,
            FailOnInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    _assert_command_error(
        () -> validate_plant_command_payload(failed,
            Float32[0, 2, 0]),
        :payload,
        :out_of_range_failure,
    )

    clipped_payload = Float32[0, 2, 0]
    clipped = _schema_fixture(
        Float32;
        bounds=UniformCommandBounds(-1f0, 1f0),
        value_policy=CommandValuePolicy(
            RejectInvalidCommand,
            ClipInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    @test validate_plant_command_payload(clipped, clipped_payload) ===
        clipped_payload
    @test clipped_payload == Float32[0, 2, 0]

    application_checked = _schema_fixture(
        Float32;
        bounds=UniformCommandBounds(-1f0, 1f0),
        value_policy=CommandValuePolicy(
            RejectInvalidCommand,
            RejectInvalidCommand,
            EnforceOnApplication,
        ),
    )
    @test validate_plant_command_payload(
        application_checked,
        Float32[0, 2, 0],
    ) == Float32[0, 2, 0]

    unbounded = _schema_fixture(
        Float32;
        bounds=UnboundedCommandValues(),
    )
    @test validate_plant_command_payload(
        unbounded,
        Float32[-1f20, 0, 1f20],
    ) == Float32[-1f20, 0, 1f20]

    allocation_payload = zeros(Float32, 3)
    validate_plant_command_payload(rejected, allocation_payload)
    if coverage_instrumented()
        @test_skip "command payload allocation assertion disabled under coverage"
    else
        @test @allocated(validate_plant_command_payload(
            rejected,
            allocation_payload,
        )) == 0
    end
end

@testset "Command-schema validation failures" begin
    for (constructor, reason) in (
        (() -> PlantCommandSchemaID(Symbol("")), :empty_id),
        (() -> PlantCommandSchemaID("schema"), :invalid_id),
        (() -> PlantCommandSchemaVersion(0), :invalid_version),
        (() -> PlantCommandSchemaVersion(-1), :invalid_version),
        (() -> PlantCommandSchemaVersion(true), :invalid_version),
        (() -> PlantCommandSchemaVersion(typemax(UInt64)),
            :invalid_version),
        (() -> CommandBasisRevision(0), :invalid_basis_revision),
        (() -> CommandBasisRevision(1.0), :invalid_basis_revision),
        (() -> CommandUnit(Symbol("")), :invalid_unit),
        (() -> CommandUnit("metre"), :invalid_unit),
        (() -> CommandSignConvention(Symbol("")),
            :invalid_sign_convention),
        (() -> CommandBasis(Symbol(""), :basis), :invalid_basis),
        (() -> CommandBasis(:actuator, Symbol("")), :invalid_basis),
        (() -> CommandBasis("actuator", :basis), :invalid_basis),
        (() -> UniformCommandBounds(Float32(NaN), 1f0),
            :nonfinite_bounds),
        (() -> UniformCommandBounds(2f0, 1f0), :invalid_bounds),
        (() -> UniformCommandBounds(-1f0, 1.0), :bounds_type),
        (() -> CommandValuePolicy(
            ClipInvalidCommand,
            RejectInvalidCommand,
            ValidateOnPresentation,
        ), :invalid_nonfinite_policy),
        (() -> CommandSequencePolicy(
            AcceptSequence,
            RejectSequence,
            RejectSequence,
            AcceptSequence,
        ), :invalid_duplicate_policy),
        (() -> CommandSequencePolicy(
            RejectSequence,
            AcceptSequence,
            RejectSequence,
            AcceptSequence,
        ), :invalid_stale_policy),
        (() -> CommandSilencePolicy(
            HoldLastCommand,
            AgeFromApplication;
            timeout=PlantDuration(1),
        ), :invalid_silence_timeout),
        (() -> CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromAdmission,
        ), :missing_silence_timeout),
        (() -> CommandSilencePolicy(
            FailOnCommandSilence,
            AgeFromApplication;
            timeout=PlantDuration(0),
        ), :invalid_silence_timeout),
    )
        _assert_command_error(constructor, :schema, reason)
    end

    for (type, dimensions, reason) in (
        (Bool, (1,), :invalid_numeric_type),
        (ComplexF32, (1,), :invalid_numeric_type),
        (Real, (1,), :invalid_numeric_type),
        (Float32, (0,), :invalid_dimensions),
        (Float32, (-1,), :invalid_dimensions),
        (Float32, (true,), :invalid_dimensions),
        (Float32, [1], :invalid_dimensions),
    )
        _assert_command_error(
            () -> _schema_fixture(
                type;
                dimensions,
                bounds=UnboundedCommandValues(),
            ),
            :schema,
            reason,
        )
    end

    _assert_command_error(
        () -> _schema_fixture(Float32; endpoint=1),
        :schema,
        :invalid_endpoint,
    )
    _assert_command_error(
        () -> _schema_fixture(
            Float32;
            bounds=UniformCommandBounds(-1.0, 1.0),
        ),
        :schema,
        :bounds_type,
    )
    _assert_command_error(
        () -> _schema_fixture(
            Float32;
            bounds=UnboundedCommandValues(),
            value_policy=CommandValuePolicy(
                RejectInvalidCommand,
                ClipInvalidCommand,
                ValidateOnPresentation,
            ),
        ),
        :schema,
        :invalid_unbounded_policy,
    )
    _assert_command_error(
        () -> _schema_fixture(
            Float32;
            bounds=UnboundedCommandValues(),
            value_policy=CommandValuePolicy(
                RejectInvalidCommand,
                RejectInvalidCommand,
                EnforceOnApplication,
            ),
        ),
        :schema,
        :invalid_unbounded_policy,
    )
end

@testset "Command schemas attach to independent optic endpoints" begin
    telescope = Telescope(
        resolution=8,
        diameter=8.0,
        central_obstruction=0.0,
    )
    atmosphere = CommandSchemaTestAtmosphere()
    woofer_schema = _schema_fixture(
        Float32;
        id=:woofer_v1,
        endpoint=:woofer_command,
        dimensions=(6,),
    )
    tweeter_schema = _schema_fixture(
        Float32;
        id=:tweeter_v1,
        endpoint=:tweeter_command,
        dimensions=(12,),
        semantics=IncrementalCommand,
    )
    woofer = ControllableOpticDefinition(
        :woofer,
        CommandSchemaTestOpticModel(:woofer);
        command_schemas=(woofer_command=woofer_schema,),
    )
    tweeter = ControllableOpticDefinition(
        :tweeter,
        CommandSchemaTestOpticModel(:tweeter),
        (tweeter_schema,),
    )
    plant = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        controllable_optics=(woofer=woofer, tweeter=tweeter),
    )

    @test command_schemas(woofer) === (woofer_schema,)
    @test command_endpoint_ids(woofer) ==
        (CommandEndpointID(:woofer_command),)
    @test command_schema(woofer, :woofer_command) === woofer_schema
    @test command_schema(plant, :tweeter_command) === tweeter_schema
    @test plant_command_schema(plant, :woofer_v1) === woofer_schema
    @test command_endpoint_owner(plant, :tweeter_command) === tweeter

    bad_key = (wrong_endpoint=woofer_schema,)
    error = try
        ControllableOpticDefinition(
            :bad_key,
            CommandSchemaTestOpticModel(:bad_key);
            command_schemas=bad_key,
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantDefinitionError
    if error isa PlantDefinitionError
        @test error.component === :command_schema
        @test error.reason === :identity_mismatch
    end

    duplicate_schema_id = _schema_fixture(
        Float32;
        id=:woofer_v1,
        endpoint=:other_command,
    )
    other = ControllableOpticDefinition(
        :other,
        CommandSchemaTestOpticModel(:other),
        (duplicate_schema_id,),
    )
    error = try
        PlantDefinition(
            telescope,
            atmosphere,
            (woofer, other),
            (),
            (),
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantDefinitionError
    if error isa PlantDefinitionError
        @test error.component === :command_schema
        @test error.reason === :duplicate_id
    end
end
