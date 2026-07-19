struct PlantTopologyTestAtmosphere <: AdaptiveOpticsSim.AbstractAtmosphere end

struct PlantTopologyTestOpticalModel
    label::Symbol
end

struct PlantTopologyTestAcquisitionModel
    label::Symbol
end

struct PlantTopologyTestInvalidDefinitionStyle end

AdaptiveOpticsSim.plant_model_definition_style(
    ::Type{PlantTopologyTestOpticalModel},
) = ColdPlantModelDefinition()

AdaptiveOpticsSim.plant_model_definition_style(
    ::Type{PlantTopologyTestAcquisitionModel},
) = ColdPlantModelDefinition()

AdaptiveOpticsSim.plant_model_definition_style(
    ::Type{PlantTopologyTestInvalidDefinitionStyle},
) = nothing

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

    @test OpticalPathID(:science) == path_id(science_path)
    @test AcquisitionID(:fast_camera) == acquisition_id(fast_camera)
    @test isequal(OpticalPathID(:science), OpticalPathID(:science))
    @test !isequal(OpticalPathID(:science), OpticalPathID(:ngs))
    @test isequal(AcquisitionID(:fast_camera), AcquisitionID(:fast_camera))
    @test !isequal(AcquisitionID(:fast_camera), AcquisitionID(:slow_camera))
    @test acquisition_path_id(fast_camera) == OpticalPathID(:science)
    @test path_source(science_path) === science_source
    @test path_model(science_path).label === :science_relay
    @test acquisition_model(fast_camera).label === :fast_detector
    @test sprint(show, OpticalPathID(:science)) ==
        "OpticalPathID(:science)"
    @test sprint(show, AcquisitionID(:fast_camera)) ==
        "AcquisitionID(:fast_camera)"
    @test length(Set((OpticalPathID(:science), OpticalPathID(:science)))) == 1
    @test length(Set((OpticalPathID(:science), AcquisitionID(:science)))) == 2

    plant = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        paths=(ngs=ngs_path, science=science_path),
        acquisitions=(
            fast_camera=fast_camera,
            slow_camera=slow_camera,
            ngs_wfs=ngs_wfs,
        ),
    )

    @test plant_telescope(plant) === telescope
    @test plant_atmosphere(plant) === atmosphere
    @test path_definitions(plant) === (ngs_path, science_path)
    @test acquisition_definitions(plant) ===
        (fast_camera, slow_camera, ngs_wfs)
    @test path_definition(plant, :science) === science_path
    @test path_definition(plant, OpticalPathID(:ngs)) === ngs_path
    @test acquisition_definition(plant, :fast_camera) === fast_camera
    @test acquisition_definition(
        plant,
        AcquisitionID(:slow_camera),
    ) === slow_camera

    # Declarations are organization, not identity. Reordering paths and
    # acquisitions preserves every explicit reference, and two independent
    # acquisitions may reference one reusable optical path.
    reordered = PlantDefinition(
        telescope,
        atmosphere,
        (science_path, ngs_path),
        (ngs_wfs, slow_camera, fast_camera),
    )
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
        (ngs=ngs_path, science=science_path),
        (fast_camera=fast_camera, slow_camera=slow_camera, ngs_wfs=ngs_wfs),
    )
    @test path_definitions(named_positional) === (ngs_path, science_path)
    @test acquisition_definitions(named_positional) ===
        (fast_camera, slow_camera, ngs_wfs)

    empty_plant = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
    )
    @test isempty(path_definitions(empty_plant))
    @test isempty(acquisition_definitions(empty_plant))

    for value in (
        OpticalPathID(:science),
        AcquisitionID(:fast_camera),
        ngs_path,
        fast_camera,
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
        () -> OpticalPathDefinition(:bad, ngs, Ref(:live_state)),
        :path,
        :unsupported_model_definition,
    )
    assert_plant_definition_error(
        () -> AcquisitionDefinition(:bad, :ngs, Ref(:live_state)),
        :acquisition,
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

    duplicate_ngs_path = OpticalPathDefinition(
        :ngs,
        ngs,
        PlantTopologyTestOpticalModel(:duplicate),
    )
    assert_plant_definition_error(
        () -> PlantDefinition(
            telescope,
            atmosphere,
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
        () -> PlantDefinition(telescope, atmosphere, [ngs_path], ()),
        :path,
        :invalid_container,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (ngs_path,), [ngs_wfs]),
        :acquisition,
        :invalid_container,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (1,), ()),
        :path,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, atmosphere, (ngs_path,), (1,)),
        :acquisition,
        :invalid_definition,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(nothing, atmosphere, (), ()),
        :plant,
        :invalid_telescope,
    )
    assert_plant_definition_error(
        () -> PlantDefinition(telescope, nothing, (), ()),
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
end
