struct MCAOMOAOPathModel end

struct MCAOMOAOUnconvertibleReal <: Real end

struct MCAOMOAOFrameAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

struct MCAOMOAOZeroPupilMaterialization{P}
    destination::P
end

Plant.plant_model_definition_style(::Type{MCAOMOAOPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:MCAOMOAOFrameAcquisitionModel},
) = ColdPlantModelDefinition()

function Plant.validate_path_materialization_binding(
    materialization::MCAOMOAOZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AbstractSource,
)
    materialization.destination === input || throw(
        PlantPreparationError(
            :path,
            :prepared_binding,
            "MCAO/MOAO zero-pupil materialization belongs to another path",
        ),
    )
    return nothing
end

function Plant.validate_path_materialization_target(
    materialization::MCAOMOAOZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "MCAO/MOAO zero-pupil materialization belongs to another path"))
    Plant._require_exact_plant_product_target(
        materialization.destination, target,
        "MCAO/MOAO zero-pupil materialization destination")
    return materialization
end

function Plant.validate_path_materialization(
    materialization::MCAOMOAOZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(
            :path,
            :prepared_binding,
            "MCAO/MOAO zero-pupil materialization belongs to another path",
        ),
    )
    return input
end

function Plant.materialize_path_input!(
    materialization::MCAOMOAOZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(
            :path,
            :prepared_binding,
            "MCAO/MOAO zero-pupil materialization belongs to another path",
        ),
    )
    fill!(input.opd, zero(eltype(input.opd)))
    return input
end

function Plant.prepare_path_executor(
    ::MCAOMOAOPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source; zero_padding=1)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        context=context,
        materialization=MCAOMOAOZeroPupilMaterialization(pupil),
        optical_model=:native_deformable_mirror_mcao_moao_test,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::MCAOMOAOFrameAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(
        exposure_duration=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=path.key.backend,
    )
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(
        execution.observation;
        metadata=(
            kind=:native_deformable_mirror_mcao_moao_test,
            units=:detected_electrons,
            geometry=path.result.metadata,
        ),
    )
    return prepare_full_optical_provider(execution, products)
end

function mcao_moao_command_schema(
    endpoint::Symbol;
    T::Type{<:AbstractFloat}=Float64,
    dimensions=(1,),
    units=:metre,
    sign_convention=:positive_surface_increases_opd,
    basis=CommandBasis(:actuator, endpoint),
)
    return PlantCommandSchema(
        T,
        dimensions;
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units,
        sign_convention,
        basis,
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(T(-1e-3), T(1e-3)),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function mcao_moao_dm_definition(
    id::Symbol,
    placement,
    visibility,
    resolution::Integer;
    T::Type{<:AbstractFloat}=Float64,
)
    endpoint = Symbol(id, :_command)
    topology = SampledActuatorTopology(zeros(T, 2, 1); T=T)
    influence = DenseInfluenceMatrix(
        fill(one(T), Int(resolution)^2, 1),
    )
    model = DeformableMirrorModel(;
        topology,
        influence_model=influence,
        T=T,
    )
    return ControllableOpticDefinition(
        id,
        model,
        (mcao_moao_command_schema(endpoint; T=T),);
        placement,
        visibility,
    )
end

function mcao_moao_atmosphere(telescope::Telescope)
    T = eltype(pupil_reflectivity(telescope))
    return MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction=T[0],
        altitude=T[0],
        layer_ids=(:ground,),
        T=T,
        backend=backend(telescope),
    )
end

function mcao_moao_full_optical_fixture()
    T = Float64
    resolution = 5
    telescope = Telescope(
        resolution=resolution,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = mcao_moao_atmosphere(telescope)
    sources = (
        lgs=LGSSource(
            altitude=T(90_000),
            photon_irradiance=T(100),
            T=T,
        ),
        ngs=Source(
            band=:custom,
            wavelength=T(0.7e-6),
            photon_irradiance=T(100),
            T=T,
        ),
        science_a=Source(
            band=:custom,
            wavelength=T(0.8e-6),
            photon_irradiance=T(100),
            T=T,
        ),
        science_b=Source(
            band=:custom,
            wavelength=T(0.9e-6),
            photon_irradiance=T(100),
            T=T,
        ),
    )
    paths = map(keys(sources), values(sources)) do id, source
        OpticalPathDefinition(
            id,
            source,
            MCAOMOAOPathModel(),
        )
    end
    acquisitions = map(keys(sources)) do id
        AcquisitionDefinition(
            Symbol(id, :_camera),
            id,
            MCAOMOAOFrameAcquisitionModel(T(0.02)),
        )
    end
    optics = (
        mcao_moao_dm_definition(
            :science_b_moao,
            AtmosphericConjugatePlacement(T(5_000)),
            SelectedPathVisibility(:science_b),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :common_high,
            AtmosphericConjugatePlacement(T(10_000)),
            AllPathVisibility(),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :science_a_moao,
            AtmosphericConjugatePlacement(T(5_000)),
            SelectedPathVisibility(:science_a),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :common_ground,
            PupilPlanePlacement(),
            AllPathVisibility(),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :common_mid,
            AtmosphericConjugatePlacement(T(5_000)),
            AllPathVisibility(),
            resolution;
            T=T,
        ),
    )
    initial_commands = (
        common_ground=T(1e-9),
        common_mid=T(2e-9),
        common_high=T(3e-9),
        science_a_moao=T(4e-9),
        science_b_moao=T(8e-9),
    )
    configurations = map(optics) do optic
        id = controllable_optic_id(optic).name
        value = getproperty(initial_commands, id)
        CommandEndpointConfiguration(
            only(command_endpoint_ids(optic)),
            T[value];
            capacity=32,
        )
    end
    plant = prepare_plant(
        PlantDefinition(;
            telescope=plant_test_telescope_definition(telescope),
            atmosphere=plant_test_atmosphere_definition(atmosphere),
            controllable_optics=optics,
            paths,
            acquisitions,
        ), PLANT_TEST_HOST_TARGET;
        run_seed=0x8601,
        command_endpoints=configurations,
    )
    optical_samples = map(keys(sources)) do id
        OpticalSampleDefinition(
            id,
            PeriodicSchedule(period_ns=100_000_000),
        )
    end
    acquisition_events = map(keys(sources)) do id
        AcquisitionEventDefinition(
            Symbol(id, :_camera),
            GlobalShutterAcquisitionDefinition(
                PlantDuration(20_000_000),
            ),
            PeriodicAcquisitionStart(
                PeriodicSchedule(period_ns=1_000_000_000),
            ),
        )
    end
    event_loop = prepare_plant_event_loop(
        plant,
        PlantEventLoopDefinition(
            optical_samples,
            acquisition_events,
        ),
    )
    return (
        plant=plant,
        event_loop=event_loop,
        state=PlantEventLoopState(event_loop),
        workspace=PlantEventLoopWorkspace(event_loop),
        initial_commands=initial_commands,
    )
end

function mcao_moao_path_opd(fixture, path::Symbol)
    return path_input(prepared_path(fixture.plant, path)).opd
end

function mcao_moao_schema(fixture, optic::Symbol)
    definition = first(
        filter(
            value -> controllable_optic_id(value) ==
                ControllableOpticID(optic),
            controllable_optic_definitions(
                getfield(fixture.plant, :definition),
            ),
        ),
    )
    return only(command_schemas(definition))
end

function mcao_moao_optic_schema(optics::Tuple, optic::Symbol)
    id = ControllableOpticID(optic)
    @inbounds for definition in optics
        controllable_optic_id(definition) == id &&
            return only(command_schemas(definition))
    end
    error("test fixture has no controllable optic $id")
end

function mcao_moao_reduced_order_response(
    schema::PlantCommandSchema;
    gain::Real=1,
)
    T = command_numeric_type(schema)
    return Plant.ReducedOrderCommandResponse(
        command_endpoint_id(schema),
        reshape(T[gain], 1, 1);
        units=command_units(schema),
        sign_convention=command_sign_convention(schema),
        basis=command_basis(schema),
        basis_revision=command_basis_revision(schema),
    )
end

function mcao_moao_reduced_order_model(
    schemas::Tuple,
    ;
    T::Type{<:AbstractFloat} =
        isempty(schemas) ? Float64 : command_numeric_type(first(schemas)),
)
    responses = map(mcao_moao_reduced_order_response, schemas)
    return Plant.LinearReducedOrderAcquisitionModel(
        Plant.HarmonicDisturbanceModel(
            T[0],
            T[0];
            offsets=T[0],
        ),
        ones(T, 1, 1),
        ones(T, 1, 1),
        responses;
        measurement_kind=:modal_residual,
        measurement_units=:metre,
        residual_kind=:modal_wavefront_error,
        residual_units=:metre,
        calibration_revision=1,
        operating_envelope=(
            maximum_absolute_residual_m=T(1e-3),
            sample_period_ns=1_000_000,
        ),
        omitted_effects=(
            :diffraction,
            :detector_noise,
            :nonlinear_wavefront_sensor_response,
        ),
    )
end

function mcao_moao_reduced_order_fixture(;
    science_a_response_optics=(:common, :science_a_moao),
    science_b_response_optics=(:common, :science_b_moao),
)
    T = Float64
    resolution = 4
    telescope = Telescope(
        resolution=resolution,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = mcao_moao_atmosphere(telescope)
    sources = (
        science_a=Source(
            band=:custom,
            wavelength=T(0.8e-6),
            photon_irradiance=T(100),
            T=T,
        ),
        science_b=Source(
            band=:custom,
            wavelength=T(0.9e-6),
            photon_irradiance=T(100),
            T=T,
        ),
        uncorrected=Source(
            band=:custom,
            wavelength=T(1.0e-6),
            photon_irradiance=T(100),
            T=T,
        ),
    )
    paths = map(keys(sources), values(sources)) do id, source
        OpticalPathDefinition(id, source, MCAOMOAOPathModel())
    end
    placement = AtmosphericConjugatePlacement(T(5_000))
    optics = (
        mcao_moao_dm_definition(
            :science_b_moao,
            placement,
            SelectedPathVisibility(:science_b),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :common,
            placement,
            SelectedPathVisibility(:science_a, :science_b),
            resolution;
            T=T,
        ),
        mcao_moao_dm_definition(
            :science_a_moao,
            placement,
            SelectedPathVisibility(:science_a),
            resolution;
            T=T,
        ),
    )
    schemas_a = map(science_a_response_optics) do id
        mcao_moao_optic_schema(optics, id)
    end
    schemas_b = map(science_b_response_optics) do id
        mcao_moao_optic_schema(optics, id)
    end
    acquisitions = (
        AcquisitionDefinition(
            :reduced_science_a,
            :science_a,
            mcao_moao_reduced_order_model(schemas_a),
        ),
        AcquisitionDefinition(
            :reduced_science_b,
            :science_b,
            mcao_moao_reduced_order_model(schemas_b),
        ),
        AcquisitionDefinition(
            :reduced_uncorrected,
            :uncorrected,
            mcao_moao_reduced_order_model((); T=T),
        ),
    )
    initial_commands = (
        common=T(1e-9),
        science_a_moao=T(2e-9),
        science_b_moao=T(4e-9),
    )
    configurations = map(optics) do optic
        id = controllable_optic_id(optic).name
        CommandEndpointConfiguration(
            only(command_endpoint_ids(optic)),
            T[getproperty(initial_commands, id)];
            capacity=16,
        )
    end
    plant = prepare_plant(
        PlantDefinition(;
            telescope=plant_test_telescope_definition(telescope),
            atmosphere=plant_test_atmosphere_definition(atmosphere),
            controllable_optics=optics,
            paths,
            acquisitions,
        ), PLANT_TEST_HOST_TARGET;
        run_seed=0x8602,
        command_endpoints=configurations,
    )
    event_loop = prepare_plant_event_loop(
        plant,
        PlantEventLoopDefinition(
            (
                OpticalSampleDefinition(
                    :science_a,
                    PeriodicSchedule(period_ns=1_000_000),
                ),
                OpticalSampleDefinition(
                    :science_b,
                    PeriodicSchedule(period_ns=1_000_000),
                ),
                OpticalSampleDefinition(
                    :uncorrected,
                    PeriodicSchedule(period_ns=1_000_000),
                ),
            ),
            (
                AcquisitionEventDefinition(
                    :reduced_science_a,
                    DirectMeasurementAcquisitionDefinition(
                        PlantDuration(1_000_000),
                    ),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=2_000_000),
                    ),
                ),
                AcquisitionEventDefinition(
                    :reduced_science_b,
                    DirectMeasurementAcquisitionDefinition(
                        PlantDuration(1_000_000),
                    ),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=2_000_000),
                    ),
                ),
                AcquisitionEventDefinition(
                    :reduced_uncorrected,
                    DirectMeasurementAcquisitionDefinition(
                        PlantDuration(1_000_000),
                    ),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=2_000_000),
                    ),
                ),
            ),
        ),
    )
    return (
        plant=plant,
        event_loop=event_loop,
        state=PlantEventLoopState(event_loop),
        workspace=PlantEventLoopWorkspace(event_loop),
        optics=optics,
    )
end

function mcao_moao_reduced_measurement(fixture, acquisition::Symbol)
    owner = prepared_acquisition(fixture.plant, acquisition)
    return copy(measurement_storage(acquisition_measurement(owner)))
end

function mcao_moao_native_hot_path_allocations(
    implementation,
    state,
    workspace,
    endpoint::CommandEndpointID,
    command,
    pupil::PupilFunction,
    coupling,
)
    timestamp = PlantTimestamp(200_000_000)

    stage_controllable_optic_command!(
        implementation,
        state,
        workspace,
        endpoint,
        command,
        timestamp,
    )
    commit_controllable_optic_command!(
        implementation,
        state,
        workspace,
        endpoint,
        timestamp,
    )
    stage_bytes = @allocated stage_controllable_optic_command!(
        implementation,
        state,
        workspace,
        endpoint,
        command,
        timestamp,
    )
    commit_bytes = @allocated commit_controllable_optic_command!(
        implementation,
        state,
        workspace,
        endpoint,
        timestamp,
    )

    fill!(pupil.opd, zero(eltype(pupil.opd)))
    apply_controllable_optic_surface!(
        pupil,
        implementation,
        state,
        coupling,
    )
    fill!(pupil.opd, zero(eltype(pupil.opd)))
    apply_bytes = @allocated apply_controllable_optic_surface!(
        pupil,
        implementation,
        state,
        coupling,
    )
    return (; stage_bytes, commit_bytes, apply_bytes)
end

function mcao_moao_native_runtime_binding(
    fixture,
    optic_id::Symbol,
    path_id_value::Symbol,
)
    optic_id_value = ControllableOpticID(optic_id)
    optic_slot = findfirst(fixture.event_loop.optics) do optic
        controllable_optic_id(optic.definition) == optic_id_value
    end
    optic_slot === nothing &&
        error("test fixture has no prepared optic $optic_id_value")
    event_group_slot = findfirst(fixture.event_loop.path_groups) do group
        group.id == OpticalPathID(path_id_value)
    end
    event_group_slot === nothing &&
        error("test fixture has no path execution group $path_id_value")
    event_group = fixture.event_loop.path_groups[event_group_slot]
    binding_slot = nothing
    bindings = fixture.event_loop.optic_path_bindings
    for (coupling_slot, binding) in enumerate(
        event_group.optic_binding_start:event_group.optic_binding_stop,
    )
        prepared_controllable_optic_slot(bindings, binding) ==
            optic_slot || continue
        binding_slot = coupling_slot
        break
    end
    binding_slot === nothing && error(
        "test fixture optic $optic_id_value is hidden on $path_id_value",
    )
    optic = fixture.event_loop.optics[optic_slot]
    endpoint = only(command_endpoint_ids(optic.definition))
    return (
        implementation=optic.implementation,
        state=fixture.state.controllable_optics[optic_slot],
        workspace=fixture.workspace.controllable_optics[optic_slot],
        endpoint=endpoint,
        pupil=event_group.path.input,
        coupling=event_group.optic_couplings[binding_slot],
    )
end

@testset "Native deformable-mirror plant model contract" begin
    model = @inferred DeformableMirrorModel(n_act=2)
    @test n_actuators(model) == 4
    @test topology(model) isa ActuatorGridTopology
    @test influence_model(model) isa GaussianInfluenceWidth
    @test actuator_model(model) isa LinearStaticActuators
    @test plant_model_definition_style(typeof(model)) isa
        ColdPlantModelDefinition
    @test !Base.ismutable(model)
    @test Base.isexported(Plant, :DeformableMirrorModel)
    @test !Base.isexported(AdaptiveOpticsSim, :DeformableMirrorModel)

    width_model =
        @inferred DeformableMirrorModel(n_act=2, influence_width=0.25)
    @test influence_model(width_model).width == 0.25
    coupling_model =
        @inferred DeformableMirrorModel(n_act=2, mechanical_coupling=0.25)
    @test influence_model(coupling_model).coupling == 0.25
    explicit_width = GaussianInfluenceWidth(0.3)
    explicit_coupling = GaussianMechanicalCoupling(0.3)
    @test influence_model(DeformableMirrorModel(
        n_act=2,
        influence_model=explicit_width,
    )) === explicit_width
    @test influence_model(DeformableMirrorModel(
        n_act=2,
        influence_model=explicit_coupling,
    )) === explicit_coupling

    for operation in (
        () -> DeformableMirrorModel(),
        () -> DeformableMirrorModel(n_act=true),
        () -> DeformableMirrorModel(n_act=0),
        () -> DeformableMirrorModel(
            n_act=2,
            topology=ActuatorGridTopology(2),
        ),
        () -> DeformableMirrorModel(n_act=2, influence_width=0),
        () -> DeformableMirrorModel(
            n_act=2,
            influence_width=MCAOMOAOUnconvertibleReal(),
        ),
        () -> DeformableMirrorModel(n_act=2, mechanical_coupling=0),
        () -> DeformableMirrorModel(
            n_act=2,
            mechanical_coupling=MCAOMOAOUnconvertibleReal(),
        ),
        () -> DeformableMirrorModel(
            n_act=2,
            influence_width=0.2,
            mechanical_coupling=0.1,
        ),
        () -> DeformableMirrorModel(
            n_act=2,
            influence_model=GaussianInfluenceWidth(-0.2),
        ),
        () -> DeformableMirrorModel(
            n_act=2,
            influence_model=GaussianMechanicalCoupling(1.0),
        ),
        () -> DeformableMirrorModel(
            n_act=2,
            actuator_model=ClippedActuators(1.0, -1.0),
        ),
        () -> DeformableMirrorModel(n_act=2, actuator_model=:linear),
        () -> DeformableMirrorModel(n_act=2, misregistration=nothing),
        () -> DeformableMirrorModel(
            n_act=2,
            pupil_relay_registration=nothing,
        ),
        () -> DeformableMirrorModel(n_act=2, T=AbstractFloat),
    )
        error = try
            operation()
            nothing
        catch caught
            caught
        end
        @test error isa PlantDefinitionError
        @test error.component === :controllable_optic
    end

    @test _deformable_mirror_error_message(:sentinel) == ":sentinel"
    definition_error = try
        _throw_deformable_mirror_definition_failure(
            :sentinel,
            :test_definition_failure,
        )
        nothing
    catch caught
        caught
    end
    @test definition_error isa PlantDefinitionError
    @test definition_error.reason === :test_definition_failure
    @test_throws InterruptException _throw_deformable_mirror_definition_failure(
        InterruptException(),
        :test_definition_interrupt,
    )
end

@testset "Native DM preparation snapshots cold topology" begin
    T = Float64
    telescope = Telescope(
        resolution=4,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = mcao_moao_atmosphere(telescope)
    metadata_values = T[1]
    cold_topology = SampledActuatorTopology(
        zeros(T, 2, 1);
        metadata=(calibration=metadata_values,),
        T=T,
    )
    model = DeformableMirrorModel(
        topology=cold_topology,
        influence_model=DenseInfluenceMatrix(
            ones(T, telescope.params.resolution^2, 1),
        ),
        T=T,
    )
    definition = ControllableOpticDefinition(
        :snapshot_native_dm,
        model,
        (mcao_moao_command_schema(:snapshot_native_dm),);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )

    prepared = prepare_controllable_optic(
        model,
        definition,
        telescope,
        atmosphere,
    )
    prepared_topology = prepared.params.topology
    @test prepared_topology !== cold_topology
    @test prepared_topology.coords !== cold_topology.coords
    @test prepared_topology.active_coords !== cold_topology.active_coords
    @test prepared_topology.valid_actuators !==
        cold_topology.valid_actuators
    @test prepared_topology.active_indices !== cold_topology.active_indices
    @test prepared_topology.metadata.calibration !== metadata_values

    cold_topology.coords[1, 1] = T(0.25)
    cold_topology.active_coords[2, 1] = T(0.5)
    cold_topology.valid_actuators[1] = false
    cold_topology.active_indices[1] = 2
    metadata_values[1] = T(2)
    @test iszero(prepared_topology.coords[1, 1])
    @test iszero(prepared_topology.active_coords[2, 1])
    @test prepared_topology.valid_actuators == Bool[true]
    @test prepared_topology.active_indices == [1]
    @test prepared_topology.metadata.calibration == T[1]
end

@testset "Native DM preparation rejects incompatible schemas and models" begin
    T = Float64
    resolution = 4
    telescope = Telescope(
        resolution=resolution,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = mcao_moao_atmosphere(telescope)
    topology = SampledActuatorTopology(zeros(T, 2, 1); T=T)
    model = DeformableMirrorModel(
        topology=topology,
        influence_model=DenseInfluenceMatrix(
            ones(T, resolution^2, 1),
        ),
        T=T,
    )
    base = mcao_moao_command_schema(:native_dm)
    incompatible = (
        (
            :deformable_mirror_command_numeric_type,
            (
                mcao_moao_command_schema(
                    :native_dm;
                    T=Float32,
                ),
            ),
        ),
        (
            :deformable_mirror_command_dimensions,
            (
                mcao_moao_command_schema(
                    :native_dm;
                    dimensions=(2,),
                ),
            ),
        ),
        (
            :deformable_mirror_command_units,
            (
                mcao_moao_command_schema(
                    :native_dm;
                    units=:radian,
                ),
            ),
        ),
        (
            :deformable_mirror_command_sign_convention,
            (
                mcao_moao_command_schema(
                    :native_dm;
                    sign_convention=:positive_command_increases_residual,
                ),
            ),
        ),
        (
            :deformable_mirror_command_basis,
            (
                mcao_moao_command_schema(
                    :native_dm;
                    basis=CommandBasis(:modal, :native_dm),
                ),
            ),
        ),
        (
            :deformable_mirror_command_schema,
            (
                base,
                mcao_moao_command_schema(:second_native_dm),
            ),
        ),
    )
    for (reason, schemas) in incompatible
        definition = ControllableOpticDefinition(
            :native_dm,
            model,
            schemas;
            placement=PupilPlanePlacement(),
            visibility=AllPathVisibility(),
        )
        error = try
            prepare_controllable_optic(
                model,
                definition,
                telescope,
                atmosphere,
            )
            nothing
        catch caught
            caught
        end
        @test error isa PlantPreparationError
        @test error.component === :controllable_optic
        @test error.reason === reason
    end

    invalid_model = DeformableMirrorModel(
        topology=topology,
        influence_model=DenseInfluenceMatrix(
            ones(T, resolution^2 - 1, 1),
        ),
        T=T,
    )
    invalid_definition = ControllableOpticDefinition(
        :invalid_native_dm,
        invalid_model,
        (base,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    error = try
        prepare_controllable_optic(
            invalid_model,
            invalid_definition,
            telescope,
            atmosphere,
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    @test error.component === :controllable_optic
    @test error.reason === :deformable_mirror_preparation
end

@testset "Native DM separable runtime and structured boundary errors" begin
    T = Float64
    telescope = Telescope(
        resolution=4,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = mcao_moao_atmosphere(telescope)
    model = DeformableMirrorModel(
        n_act=2,
        mechanical_coupling=T(0.25),
        T=T,
    )
    schema = mcao_moao_command_schema(
        :separable_native_dm;
        T=T,
        dimensions=(4,),
    )
    definition = ControllableOpticDefinition(
        :separable_native_dm,
        model,
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    prepared = prepare_controllable_optic(
        model,
        definition,
        telescope,
        atmosphere,
    )
    endpoint = command_endpoint_id(schema)
    initial_command = zeros(T, 4)
    state = prepare_controllable_optic_state(
        prepared,
        definition,
        (endpoint,),
        (initial_command,),
    )
    workspace = prepare_controllable_optic_workspace(prepared)
    @test state.active.state.coefs_grid !== nothing
    @test state.active.state.separable_tmp !== nothing
    @test workspace.staged.state.coefs_grid !== nothing
    @test workspace.staged.state.separable_tmp !== nothing

    command = T[1e-9, 2e-9, 3e-9, 4e-9]
    stage_controllable_optic_command!(
        prepared,
        state,
        workspace,
        endpoint,
        command,
        PlantTimestamp(1),
    )
    commit_controllable_optic_command!(
        prepared,
        state,
        workspace,
        endpoint,
        PlantTimestamp(1),
    )
    @test any(!iszero, surface_opd(state.active))

    storage_error = try
        prepare_controllable_optic_state(
            prepared,
            definition,
            (endpoint,),
            (zeros(Float32, 4),),
        )
        nothing
    catch caught
        caught
    end
    @test storage_error isa PlantPreparationError
    @test storage_error.reason === :deformable_mirror_command_storage

    stage_error = try
        stage_controllable_optic_command!(
            prepared,
            state,
            workspace,
            endpoint,
            zeros(Float32, 4),
            PlantTimestamp(2),
        )
        nothing
    catch caught
        caught
    end
    @test stage_error isa PlantCommandError
    @test stage_error.stage === :physical_application
    @test stage_error.reason === :deformable_mirror_command_storage

    incomplete_error = try
        _deformable_mirror_separable_runtime(
            zeros(T, 4, 2),
            nothing,
            zeros(T, 4),
            zeros(T, 4, 4),
        )
        nothing
    catch caught
        caught
    end
    @test incomplete_error isa PlantPreparationError
    @test incomplete_error.reason === :deformable_mirror_preparation

    input_error = try
        _deformable_mirror_surface_metadata(prepared, zeros(T, 4, 4))
        nothing
    catch caught
        caught
    end
    @test input_error isa PlantPreparationError
    @test input_error.reason === :unsupported_path_input

    preparation_error = try
        _throw_deformable_mirror_preparation_failure(
            ArgumentError("test preparation failure"),
            ControllableOpticID(:separable_native_dm),
        )
        nothing
    catch caught
        caught
    end
    @test preparation_error isa PlantPreparationError
    @test preparation_error.reason === :deformable_mirror_preparation
    @test_throws InterruptException _throw_deformable_mirror_preparation_failure(
        InterruptException(),
        ControllableOpticID(:separable_native_dm),
    )
end

@testset "Native DM staged publication and allocation contract" begin
    fixture = mcao_moao_full_optical_fixture()
    binding = mcao_moao_native_runtime_binding(
        fixture,
        :common_mid,
        :science_a,
    )
    command = [7e-9]
    @test binding.state.active.params === binding.workspace.staged.params
    @test binding.state.active.state.modes ===
        binding.workspace.staged.state.modes
    @test binding.state.active.state.coefs !==
        binding.workspace.staged.state.coefs
    @test surface_opd(binding.state.active) !==
        surface_opd(binding.workspace.staged)
    @test binding.implementation.surface_metadata.sampling ==
        binding.pupil.metadata.sampling
    @test binding.implementation.surface_metadata.origin ==
        binding.pupil.metadata.origin
    storage_ids = Set((
        objectid(surface_opd(binding.state.active)),
        objectid(surface_opd(binding.workspace.staged)),
    ))
    bytes = mcao_moao_native_hot_path_allocations(
        binding.implementation,
        binding.state,
        binding.workspace,
        binding.endpoint,
        command,
        binding.pupil,
        binding.coupling,
    )
    @test all(isapprox.(
        binding.pupil.opd,
        7e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test Set((
        objectid(surface_opd(binding.state.active)),
        objectid(surface_opd(binding.workspace.staged)),
    )) == storage_ids
    if coverage_instrumented()
        @test bytes.stage_bytes >= 0
        @test bytes.commit_bytes >= 0
        @test bytes.apply_bytes >= 0
    else
        @test bytes.stage_bytes == 0
        @test bytes.commit_bytes == 0
        @test bytes.apply_bytes == 0
    end
end

@testset "Reduced-order responses match exact path visibility" begin
    fixture = mcao_moao_reduced_order_fixture()
    run_plant_events_until!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(1_000_000),
    )
    @test mcao_moao_reduced_measurement(
        fixture,
        :reduced_science_a,
    ) ≈ [3e-9]
    @test mcao_moao_reduced_measurement(
        fixture,
        :reduced_science_b,
    ) ≈ [5e-9]
    @test mcao_moao_reduced_measurement(
        fixture,
        :reduced_uncorrected,
    ) == [0.0]
    @test all(
        path -> all(iszero, path_result(path).values),
        prepared_paths(fixture.plant),
    )

    schema_a = mcao_moao_optic_schema(
        fixture.optics,
        :science_a_moao,
    )
    admission = admit_plant_command!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        PlantCommand(
            schema_a,
            1,
            PlantTimestamp(2_000_000),
            [3e-9],
        ),
        PlantTimestamp(1_000_001),
    )
    @test command_admission_status(admission) ==
        CommandAdmittedPending
    run_plant_events_until!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(3_000_000),
    )
    @test mcao_moao_reduced_measurement(
        fixture,
        :reduced_science_a,
    ) ≈ [4e-9]
    @test mcao_moao_reduced_measurement(
        fixture,
        :reduced_science_b,
    ) ≈ [5e-9]

    for response_optics in (
        (:common,),
        (:common, :science_b_moao),
    )
        error = try
            mcao_moao_reduced_order_fixture(
                science_a_response_optics=response_optics,
            )
            nothing
        catch caught
            caught
        end
        @test error isa PlantPreparationError
        @test error.component === :reduced_order
        @test error.reason === :path_visibility
    end
end

@testset "Common MCAO surfaces and target-local MOAO isolation" begin
    fixture = mcao_moao_full_optical_fixture()
    @test step_plant_events!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
    ) == PlantTimestamp(0)

    common = 6e-9
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :lgs),
        common;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :ngs),
        common;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_a),
        10e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_b),
        14e-9;
        rtol=0,
        atol=16eps(Float64),
    ))

    schema = mcao_moao_schema(fixture, :science_a_moao)
    admission = admit_plant_command!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        PlantCommand(
            schema,
            1,
            PlantTimestamp(100_000_000),
            [5e-9],
        ),
        PlantTimestamp(1),
    )
    @test command_admission_status(admission) ==
        CommandAdmittedPending
    @test step_plant_events!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
    ) == PlantTimestamp(20_000_000)
    @test step_plant_events!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
    ) == PlantTimestamp(100_000_000)
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_a),
        11e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_b),
        14e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :ngs),
        common;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :lgs),
        common;
        rtol=0,
        atol=16eps(Float64),
    ))

    clear_command_dispositions!(fixture.workspace)
    common_schema = mcao_moao_schema(fixture, :common_mid)
    common_admission = admit_plant_command!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        PlantCommand(
            common_schema,
            1,
            PlantTimestamp(200_000_000),
            [6e-9],
        ),
        PlantTimestamp(100_000_001),
    )
    @test command_admission_status(common_admission) ==
        CommandAdmittedPending
    @test step_plant_events!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
    ) == PlantTimestamp(200_000_000)
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :ngs),
        10e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_a),
        15e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_b),
        18e-9;
        rtol=0,
        atol=16eps(Float64),
    ))

    clear_command_dispositions!(fixture.workspace)
    science_b_schema = mcao_moao_schema(
        fixture,
        :science_b_moao,
    )
    transaction = Plant.PlantCommandTransaction(
        PlantCommand(
            schema,
            2,
            PlantTimestamp(300_000_000),
            [7e-9],
        ),
        PlantCommand(
            science_b_schema,
            1,
            PlantTimestamp(300_000_000),
            [9e-9],
        ),
    )
    transaction_admission = Plant.admit_plant_command_transaction!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
        transaction,
        PlantTimestamp(200_000_001),
    )
    @test command_admission_status(transaction_admission) ==
        CommandAdmittedPending
    @test Plant.command_transaction_member_count(
        transaction_admission,
    ) == 2
    @test step_plant_events!(
        fixture.event_loop,
        fixture.state,
        fixture.workspace,
    ) == PlantTimestamp(300_000_000)
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_a),
        17e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :science_b),
        19e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
    @test all(isapprox.(
        mcao_moao_path_opd(fixture, :lgs),
        10e-9;
        rtol=0,
        atol=16eps(Float64),
    ))
end
