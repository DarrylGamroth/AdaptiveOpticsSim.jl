struct PlantTopologyTestAtmosphere <: AdaptiveOpticsSim.AbstractAtmosphere end

struct PlantTopologyTestOpticalModel
    label::Symbol
end

struct PlantTopologyTestAcquisitionModel
    label::Symbol
end

struct PlantTopologyTestControllableOpticModel
    label::Symbol
    conjugate::Symbol
end

struct PlantTopologyTestInvalidDefinitionStyle end

Plant.plant_model_definition_style(
    ::Type{PlantTopologyTestOpticalModel},
) = ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{PlantTopologyTestAcquisitionModel},
) = ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{PlantTopologyTestControllableOpticModel},
) = ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{PlantTopologyTestInvalidDefinitionStyle},
) = nothing

@inline _topology_endpoint_name(endpoint::Symbol) = endpoint
@inline _topology_endpoint_name(endpoint::CommandEndpointID) = endpoint.name

function topology_command_schema(endpoint;
    id=Symbol(_topology_endpoint_name(endpoint), :_schema))
    return PlantCommandSchema(
        Float32,
        (1,);
        id,
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :test_actuator),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(-1f0, 1f0),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function captured_plant_definition_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_plant_definition_error(f, component::Symbol, reason::Symbol)
    error = captured_plant_definition_error(f)
    @test error isa PlantDefinitionError
    if error isa PlantDefinitionError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

@testset "Immutable plant topology" begin
    for name in (
        :ControllableOpticID,
        :CommandEndpointID,
        :ControllableOpticDefinition,
    )
        @test Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
    for name in (
        :controllable_optic_id,
        :controllable_optic_model,
        :command_schemas,
        :command_endpoint_ids,
        :command_schema,
        :plant_command_schema,
        :controllable_optic_definitions,
        :controllable_optic_definition,
        :command_endpoint_owner,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end

    telescope = Telescope(
        resolution=8,
        diameter=8.0,
        central_obstruction=0.0,
    )
    atmosphere = PlantTopologyTestAtmosphere()
    ngs = Source(band=:I, magnitude=0.0)
    science_source = Source(
        band=:custom,
        wavelength=1.65e-6,
        photon_irradiance=2.0,
    )

    ngs_path = OpticalPathDefinition(
        :ngs,
        ngs,
        PlantTopologyTestOpticalModel(:ngs_relay),
    )
    science_path = OpticalPathDefinition(
        OpticalPathID(:science),
        science_source;
        optical_model=PlantTopologyTestOpticalModel(:science_relay),
    )
    fast_camera = AcquisitionDefinition(
        :fast_camera,
        :science,
        PlantTopologyTestAcquisitionModel(:fast_detector),
    )
    slow_camera = AcquisitionDefinition(
        AcquisitionID(:slow_camera),
        OpticalPathID(:science);
        acquisition_model=PlantTopologyTestAcquisitionModel(:slow_detector),
    )
    ngs_wfs = AcquisitionDefinition(
        :ngs_wfs,
        :ngs,
        PlantTopologyTestAcquisitionModel(:wfs),
    )
    woofer_feedback = AcquisitionDefinition(
        :woofer_feedback,
        :science,
        PlantTopologyTestAcquisitionModel(:device_feedback),
    )
    woofer = ControllableOpticDefinition(
        :woofer,
        PlantTopologyTestControllableOpticModel(:woofer, :pupil);
        command_schemas=(
            woofer_command=topology_command_schema(:woofer_command),
        ),
    )
    tweeter_schema = topology_command_schema(
        CommandEndpointID(:tweeter_command),
    )
    tweeter = ControllableOpticDefinition(
        ControllableOpticID(:tweeter),
        PlantTopologyTestControllableOpticModel(:tweeter, :pupil),
        (tweeter_schema,),
    )
    inferred_tweeter_schema = topology_command_schema(
        CommandEndpointID(:inferred_tweeter_command),
    )
    inferred_tweeter = @inferred ControllableOpticDefinition(
        ControllableOpticID(:inferred_tweeter),
        PlantTopologyTestControllableOpticModel(:inferred_tweeter, :pupil),
        (inferred_tweeter_schema,),
    )
    @test @inferred(command_endpoint_ids(inferred_tweeter)) ===
        (CommandEndpointID(:inferred_tweeter_command),)
    segmented = ControllableOpticDefinition(
        :segmented,
        PlantTopologyTestControllableOpticModel(:segmented, :pupil),
        (
            topology_command_schema(:segment_a),
            topology_command_schema(:segment_b),
        ),
    )

    @test OpticalPathID(:science) == path_id(science_path)
    @test AcquisitionID(:fast_camera) == acquisition_id(fast_camera)
    @test ControllableOpticID(:woofer) == controllable_optic_id(woofer)
    @test command_endpoint_ids(woofer) ==
        (CommandEndpointID(:woofer_command),)
    @test command_endpoint_ids(segmented) ==
        (CommandEndpointID(:segment_a), CommandEndpointID(:segment_b))
    @test command_schema(woofer, :woofer_command) ===
        only(command_schemas(woofer))
    @test isequal(OpticalPathID(:science), OpticalPathID(:science))
    @test !isequal(OpticalPathID(:science), OpticalPathID(:ngs))
    @test isequal(AcquisitionID(:fast_camera), AcquisitionID(:fast_camera))
    @test !isequal(AcquisitionID(:fast_camera), AcquisitionID(:slow_camera))
    @test isequal(ControllableOpticID(:woofer),
        ControllableOpticID(:woofer))
    @test !isequal(ControllableOpticID(:woofer),
        ControllableOpticID(:tweeter))
    @test isequal(CommandEndpointID(:woofer_command),
        CommandEndpointID(:woofer_command))
    @test acquisition_path_id(fast_camera) == OpticalPathID(:science)
    @test path_source(science_path) === science_source
    @test path_model(science_path).label === :science_relay
    @test acquisition_model(fast_camera).label === :fast_detector
    @test controllable_optic_model(woofer).label === :woofer
    @test controllable_optic_model(woofer).conjugate === :pupil
    @test sprint(show, OpticalPathID(:science)) ==
        "OpticalPathID(:science)"
    @test sprint(show, AcquisitionID(:fast_camera)) ==
        "AcquisitionID(:fast_camera)"
    @test sprint(show, ControllableOpticID(:woofer)) ==
        "ControllableOpticID(:woofer)"
    @test sprint(show, CommandEndpointID(:woofer_command)) ==
        "CommandEndpointID(:woofer_command)"
    @test length(Set((OpticalPathID(:science), OpticalPathID(:science)))) == 1
    @test length(Set((OpticalPathID(:science), AcquisitionID(:science)))) == 2
    @test length(Set((
        OpticalPathID(:shared),
        AcquisitionID(:shared),
        ControllableOpticID(:shared),
        CommandEndpointID(:shared),
    ))) == 4

    plant = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        controllable_optics=(
            woofer=woofer,
            tweeter=tweeter,
            segmented=segmented,
        ),
        paths=(ngs=ngs_path, science=science_path),
        acquisitions=(
            fast_camera=fast_camera,
            slow_camera=slow_camera,
            ngs_wfs=ngs_wfs,
            woofer_feedback=woofer_feedback,
        ),
    )

    @test plant_telescope(plant) === telescope
    @test plant_atmosphere(plant) === atmosphere
    @test controllable_optic_definitions(plant) ===
        (woofer, tweeter, segmented)
    @test path_definitions(plant) === (ngs_path, science_path)
    @test acquisition_definitions(plant) ===
        (fast_camera, slow_camera, ngs_wfs, woofer_feedback)
    @test controllable_optic_definition(plant, :woofer) === woofer
    @test controllable_optic_definition(
        plant,
        ControllableOpticID(:tweeter),
    ) === tweeter
    @test command_endpoint_owner(plant, :woofer_command) === woofer
    @test command_schema(plant, :tweeter_command) === tweeter_schema
    @test plant_command_schema(plant, :tweeter_command_schema) ===
        tweeter_schema
    @test command_endpoint_owner(
        plant,
        CommandEndpointID(:segment_b),
    ) === segmented
    @test path_definition(plant, :science) === science_path
    @test path_definition(plant, OpticalPathID(:ngs)) === ngs_path
    @test acquisition_definition(plant, :fast_camera) === fast_camera
    @test acquisition_definition(
        plant,
        AcquisitionID(:slow_camera),
    ) === slow_camera

    # Declarations are organization, not identity. Reordering optics, paths,
    # and acquisitions preserves every explicit reference, and two independent
    # acquisitions may reference one reusable optical path.
    reordered = PlantDefinition(
        telescope,
        atmosphere,
        (segmented, tweeter, woofer),
        (science_path, ngs_path),
        (woofer_feedback, ngs_wfs, slow_camera, fast_camera),
    )
    @test controllable_optic_definition(reordered, :woofer) === woofer
    @test command_endpoint_owner(reordered, :tweeter_command) === tweeter
    @test command_schema(reordered, :tweeter_command) === tweeter_schema
    @test plant_command_schema(reordered, :tweeter_command_schema) ===
        tweeter_schema
    @test path_definition(reordered, :science) === science_path
    @test acquisition_path_id(
        acquisition_definition(reordered, :fast_camera),
    ) == OpticalPathID(:science)
    @test acquisition_path_id(
        acquisition_definition(reordered, :slow_camera),
    ) == OpticalPathID(:science)

    named_positional = PlantDefinition(
        telescope,
        atmosphere,
        (woofer=woofer, tweeter=tweeter, segmented=segmented),
        (ngs=ngs_path, science=science_path),
        (fast_camera=fast_camera, slow_camera=slow_camera, ngs_wfs=ngs_wfs,
            woofer_feedback=woofer_feedback),
    )
    @test controllable_optic_definitions(named_positional) ===
        (woofer, tweeter, segmented)
    @test path_definitions(named_positional) === (ngs_path, science_path)
    @test acquisition_definitions(named_positional) ===
        (fast_camera, slow_camera, ngs_wfs, woofer_feedback)

    empty_plant = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
    )
    @test isempty(path_definitions(empty_plant))
    @test isempty(acquisition_definitions(empty_plant))
    @test isempty(controllable_optic_definitions(empty_plant))

    for value in (
        OpticalPathID(:science),
        AcquisitionID(:fast_camera),
        ControllableOpticID(:woofer),
        CommandEndpointID(:woofer_command),
        ngs_path,
        fast_camera,
        woofer,
        plant,
    )
        @test !Base.ismutable(value)
    end
    for absent in (
        :schedule,
        :workspace,
        :state,
        :rng,
        :queue,
        :transport,
        :hil,
    )
        @test !hasproperty(ngs_path, absent)
        @test !hasproperty(fast_camera, absent)
        @test !hasproperty(woofer, absent)
        @test !hasproperty(plant, absent)
    end

    assert_plant_definition_error(
        () -> OpticalPathID(Symbol("")),
        :path,
        :empty_id,
    )
    assert_plant_definition_error(
        () -> AcquisitionID(Symbol("")),
        :acquisition,
        :empty_id,
    )
    assert_plant_definition_error(
        () -> ControllableOpticID(Symbol("")),
        :controllable_optic,
        :empty_id,
    )
    assert_plant_definition_error(
        () -> CommandEndpointID(Symbol("")),
        :command_endpoint,
        :empty_id,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(1, ngs, PlantTopologyTestOpticalModel(:x)),
        :path,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> AcquisitionDefinition(:camera, 1,
            PlantTopologyTestAcquisitionModel(:x)),
        :path,
        :invalid_id,
    )
    generic_command_schema = topology_command_schema(
        :command;
        id=:generic_command_schema,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            1,
            PlantTopologyTestControllableOpticModel(:x, :pupil),
            (generic_command_schema,),
        ),
        :controllable_optic,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            PlantTopologyTestControllableOpticModel(:x, :pupil),
            (1,),
        ),
        :command_schema,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(:bad, nothing,
            PlantTopologyTestOpticalModel(:x)),
        :path,
        :missing_source,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(:bad, 1,
            PlantTopologyTestOpticalModel(:x)),
        :path,
        :invalid_source,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(:bad, ngs, nothing),
        :path,
        :missing_model,
    )
    assert_plant_definition_error(
        () -> AcquisitionDefinition(:bad, :ngs, nothing),
        :acquisition,
        :missing_model,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            nothing,
            (generic_command_schema,),
        ),
        :controllable_optic,
        :missing_model,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(:bad, ngs, Ref(:live_state)),
        :path,
        :unsupported_model_definition,
    )
    assert_plant_definition_error(
        () -> AcquisitionDefinition(:bad, :ngs, Ref(:live_state)),
        :acquisition,
        :unsupported_model_definition,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            Ref(:live_state),
            (generic_command_schema,),
        ),
        :controllable_optic,
        :unsupported_model_definition,
    )
    live_detector = Detector()
    @test Base.ismutable(live_detector.state)
    assert_plant_definition_error(
        () -> AcquisitionDefinition(:bad, :ngs, live_detector),
        :acquisition,
        :unsupported_model_definition,
    )
    assert_plant_definition_error(
        () -> OpticalPathDefinition(
            :bad,
            ngs,
            PlantTopologyTestInvalidDefinitionStyle(),
        ),
        :path,
        :invalid_model_definition_style,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            PlantTopologyTestInvalidDefinitionStyle(),
            (generic_command_schema,),
        ),
        :controllable_optic,
        :invalid_model_definition_style,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            PlantTopologyTestControllableOpticModel(:bad, :pupil),
            (),
        ),
        :controllable_optic,
        :missing_schema,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            PlantTopologyTestControllableOpticModel(:bad, :pupil),
            [generic_command_schema],
        ),
        :command_schema,
        :invalid_container,
    )
    duplicate_endpoint_a = topology_command_schema(
        :duplicate_endpoint;
        id=:duplicate_endpoint_a,
    )
    duplicate_endpoint_b = topology_command_schema(
        :duplicate_endpoint;
        id=:duplicate_endpoint_b,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad,
            PlantTopologyTestControllableOpticModel(:bad, :pupil),
            (duplicate_endpoint_a, duplicate_endpoint_b),
        ),
        :command_endpoint,
        :duplicate_id,
    )
    duplicate_schema_a = topology_command_schema(
        :schema_endpoint_a;
        id=:duplicate_schema,
    )
    duplicate_schema_b = topology_command_schema(
        :schema_endpoint_b;
        id=:duplicate_schema,
    )
    assert_plant_definition_error(
        () -> ControllableOpticDefinition(
            :bad_schema,
            PlantTopologyTestControllableOpticModel(:bad, :pupil),
            (duplicate_schema_a, duplicate_schema_b),
        ),
        :command_schema,
        :duplicate_id,
    )

    duplicate_woofer = ControllableOpticDefinition(
        :woofer,
        PlantTopologyTestControllableOpticModel(:duplicate, :pupil),
        (topology_command_schema(:duplicate_command),),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
            (woofer, duplicate_woofer),
            (),
            (),
        ),
        :controllable_optic,
        :duplicate_id,
    )
    shared_endpoint_optic = ControllableOpticDefinition(
        :shared_endpoint_optic,
        PlantTopologyTestControllableOpticModel(:shared, :pupil),
        (topology_command_schema(
            :woofer_command;
            id=:shared_endpoint_schema,
        ),),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
            (woofer, shared_endpoint_optic),
            (),
            (),
        ),
        :command_endpoint,
        :duplicate_owner,
    )

    duplicate_ngs_path = OpticalPathDefinition(
        :ngs,
        ngs,
        PlantTopologyTestOpticalModel(:duplicate),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
            (),
            (ngs_path, duplicate_ngs_path),
            (),
        ),
        :path,
        :duplicate_id,
    )
    duplicate_fast_camera = AcquisitionDefinition(
        :fast_camera,
        :ngs,
        PlantTopologyTestAcquisitionModel(:duplicate),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
            (),
            (ngs_path, science_path),
            (fast_camera, duplicate_fast_camera),
        ),
        :acquisition,
        :duplicate_id,
    )
    unknown_path_acquisition = AcquisitionDefinition(
        :unknown_camera,
        :unknown_path,
        PlantTopologyTestAcquisitionModel(:camera),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
            (),
            (ngs_path,),
            (unknown_path_acquisition,),
        ),
        :acquisition,
        :unknown_path,
    )

    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope=telescope,
            atmosphere=atmosphere,
            controllable_optics=(wrong_name=woofer,),
        ),
        :controllable_optic,
        :identity_mismatch,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope=telescope,
            atmosphere=atmosphere,
            paths=(wrong_name=ngs_path,),
        ),
        :path,
        :identity_mismatch,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope=telescope,
            atmosphere=atmosphere,
            paths=(ngs=ngs_path,),
            acquisitions=(wrong_name=ngs_wfs,),
        ),
        :acquisition,
        :identity_mismatch,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (), [ngs_path], ()),
        :path,
        :invalid_container,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (), (ngs_path,),
            [ngs_wfs]),
        :acquisition,
        :invalid_container,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, [woofer], (), ()),
        :controllable_optic,
        :invalid_container,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (1,), (), ()),
        :controllable_optic,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (), (1,), ()),
        :path,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (), (ngs_path,), (1,)),
        :acquisition,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(nothing, atmosphere, (), (), ()),
        :plant,
        :invalid_telescope,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, nothing, (), (), ()),
        :plant,
        :invalid_atmosphere,
    )
    assert_plant_definition_error(
        () -> path_definition(plant, :unknown),
        :path,
        :unknown_id,
    )
    assert_plant_definition_error(
        () -> path_definition(plant, 1),
        :path,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> acquisition_definition(plant, :unknown),
        :acquisition,
        :unknown_id,
    )
    assert_plant_definition_error(
        () -> acquisition_definition(plant, 1),
        :acquisition,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> controllable_optic_definition(plant, :unknown),
        :controllable_optic,
        :unknown_id,
    )
    assert_plant_definition_error(
        () -> controllable_optic_definition(plant, 1),
        :controllable_optic,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> command_endpoint_owner(plant, :unknown),
        :command_endpoint,
        :unknown_id,
    )
    assert_plant_definition_error(
        () -> command_endpoint_owner(plant, 1),
        :command_endpoint,
        :invalid_id,
    )
    assert_plant_definition_error(
        () -> plant_command_schema(plant, :unknown),
        :command_schema,
        :unknown_id,
    )
    assert_plant_definition_error(
        () -> plant_command_schema(plant, 1),
        :command_schema,
        :invalid_id,
    )

    preparation_error = captured_plant_definition_error() do
        prepare_plant(PlantDefinition(
            telescope=telescope,
            atmosphere=atmosphere,
            controllable_optics=(woofer,),
        ); run_seed=0x1234)
    end
    @test preparation_error isa PlantPreparationError
    if preparation_error isa PlantPreparationError
        @test preparation_error.component === :command_endpoint
        @test preparation_error.reason === :configuration_count
    end

    unsupported_model_error = captured_plant_definition_error() do
        prepare_plant(PlantDefinition(
            telescope=telescope,
            atmosphere=atmosphere,
            controllable_optics=(woofer,),
        ); run_seed=0x1234,
            command_endpoints=(
                CommandEndpointConfiguration(:woofer_command,
                    Float32[0]; capacity=1),
            ))
    end
    @test unsupported_model_error isa PlantPreparationError
    if unsupported_model_error isa PlantPreparationError
        @test unsupported_model_error.component === :controllable_optic
        @test unsupported_model_error.reason === :unsupported_model
    end
end
