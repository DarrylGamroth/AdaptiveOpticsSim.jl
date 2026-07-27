struct DeviceBatchTestPathModel
    zero_padding::Int
end

struct DeviceBatchTestAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

function device_batch_test_command_schema(::Type{T}) where {
    T<:AbstractFloat}
    return PlantCommandSchema(
        T,
        (4,);
        id=:device_batch_dm_schema,
        version=1,
        endpoint=:device_batch_dm_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :device_batch_dm),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function device_batch_test_aberration_metadata(
    prototype::PupilFunction,
    surface::AbstractMatrix,
)
    return OpticalPlaneMetadata(
        PupilPlane(),
        surface;
        coordinate_domain=MetricCoordinates(),
        sampling=prototype.metadata.sampling,
        origin=prototype.metadata.origin,
        orientation=prototype.metadata.orientation,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
    )
end

function device_batch_test_physical_definitions(
    telescope::Telescope,
    backend::AdaptiveOpticsSim.AbstractArrayBackend,
    ::Type{T},
    ;
    selected_path::Symbol=:alpha,
) where {T<:AbstractFloat}
    prototype = PupilFunction(telescope; T, backend)
    static_opd = similar(prototype.opd)
    fill!(static_opd, T(2e-9))
    static_metadata =
        device_batch_test_aberration_metadata(prototype, static_opd)
    common_aberration = SampledAberrationDefinition(
        :device_batch_static,
        OPDMap(static_opd),
        static_metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMAdditive(),
    )
    ncpa_opd = similar(prototype.opd)
    ncpa_host = Matrix{T}(undef, size(ncpa_opd))
    scale = T(1e-9) / T(size(ncpa_host, 1))
    @inbounds for column in axes(ncpa_host, 2),
        row in axes(ncpa_host, 1)
        ncpa_host[row, column] = scale * T(row - column)
    end
    copyto!(ncpa_opd, ncpa_host)
    ncpa_metadata =
        device_batch_test_aberration_metadata(prototype, ncpa_opd)
    alpha_ncpa = SampledAberrationDefinition(
        :device_batch_alpha_ncpa,
        NCPA(ncpa_opd, nothing, nothing),
        ncpa_metadata;
        placement=PupilPlanePlacement(),
        visibility=SelectedPathVisibility(selected_path),
        application=DMAdditive(),
    )
    schema = device_batch_test_command_schema(T)
    optic = ControllableOpticDefinition(
        :device_batch_dm,
        DeformableMirrorModel(
            n_act=2,
            influence_width=T(0.3),
            T=T,
        ),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    initial_command = T[1e-9, 0, -5e-10, 0]
    configuration = CommandEndpointConfiguration(
        :device_batch_dm_command,
        initial_command;
        capacity=4,
        backend,
    )
    return (
        sampled_aberrations=(common_aberration, alpha_ncpa),
        optic,
        schema,
        configuration,
    )
end

Plant.plant_model_definition_style(
    ::Type{DeviceBatchTestPathModel},
) = ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{<:DeviceBatchTestAcquisitionModel},
) = ColdPlantModelDefinition()

function Plant.prepare_path_executor(
    model::DeviceBatchTestPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    execution = prepare_direct_imaging(
        pupil,
        source;
        zero_padding=model.zero_padding,
    )
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(execution),
        execution;
        materialization=prepare_pupil_opd_materialization(
            atmosphere,
            telescope,
            source,
            pupil,
        ),
        optical_model=(
            kind=:device_batch_test_direct_imaging,
            zero_padding=model.zero_padding,
        ),
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::DeviceBatchTestAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    result = path_result(path)
    T = eltype(result.values)
    detector = Detector(
        integration_time=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=path_result_key(path).backend,
    )
    execution = FrameAcquisitionExecution(detector, result)
    products = AcquisitionProducts(
        execution.observation;
        metadata=(
            kind=:device_batch_test_frame,
            units=:detected_electrons,
            geometry=result.metadata,
            detector=detector_export_metadata(detector),
        ),
    )
    return prepare_full_optical_provider(execution, products)
end

@inline function prepare_device_batch_test_event_loop(
    plant,
    definition,
    ::Val{:public},
)
    return prepare_plant_event_loop(plant, definition)
end

@inline function prepare_device_batch_test_event_loop(
    plant,
    definition,
    ::Val{:all},
)
    return Plant._prepare_plant_event_loop(plant, definition, Val(:all))
end

@inline function prepare_device_batch_test_event_loop(
    plant,
    definition,
    ::Val{:none},
)
    return Plant._prepare_plant_event_loop(plant, definition, Val(:none))
end

function device_batch_test_fixture(;
    backend::AdaptiveOpticsSim.AbstractArrayBackend=CPUBackend(),
    selection::Val=Val(:public),
    T::Type{<:AbstractFloat}=Float64,
    include_beta::Bool=true,
    include_unequal_rate::Bool=true,
    include_phase_offset::Bool=true,
    include_origin_offset::Bool=true,
    include_lgs::Bool=true,
    include_physical_state::Bool=true,
)
    telescope = Telescope(
        resolution=8,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[0.65, 0.35],
        wind_speed=T[7, 11],
        wind_direction=T[20, 125],
        altitude=T[0, 5_000],
        layer_ids=(:ground, :high),
        T=T,
        backend=backend,
    )
    alpha_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(60),
        coordinates=(zero(T), zero(T)),
        T=T,
    )
    beta_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(45),
        coordinates=(T(4), T(35)),
        T=T,
    )
    unequal_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(35),
        coordinates=(T(7), T(110)),
        T=T,
    )
    phase_offset_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(30),
        coordinates=(T(6), T(250)),
        T=T,
    )
    origin_offset_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(25),
        coordinates=(T(5), T(290)),
        T=T,
    )
    lgs_source = LGSSource(
        wavelength=T(589e-9),
        photon_irradiance=T(50),
        coordinates=(T(9), T(210)),
        altitude=T(90_000),
        T=T,
    )

    path_definitions = OpticalPathDefinition[
        OpticalPathDefinition(
            :alpha,
            alpha_source,
            DeviceBatchTestPathModel(2),
        ),
    ]
    include_beta && push!(
        path_definitions,
        OpticalPathDefinition(
            :beta,
            beta_source,
            DeviceBatchTestPathModel(2),
        ),
    )
    include_unequal_rate && push!(
        path_definitions,
        OpticalPathDefinition(
            :unequal,
            unequal_source,
            DeviceBatchTestPathModel(2),
        ),
    )
    include_phase_offset && push!(
        path_definitions,
        OpticalPathDefinition(
            :phase_offset,
            phase_offset_source,
            DeviceBatchTestPathModel(2),
        ),
    )
    include_origin_offset && push!(
        path_definitions,
        OpticalPathDefinition(
            :origin_offset,
            origin_offset_source,
            DeviceBatchTestPathModel(2),
        ),
    )
    include_lgs && push!(
        path_definitions,
        OpticalPathDefinition(
            :lgs,
            lgs_source,
            DeviceBatchTestPathModel(2),
        ),
    )

    acquisition_definitions = AcquisitionDefinition[]
    sample_definitions = OpticalSampleDefinition[]
    event_definitions = AcquisitionEventDefinition[]
    for definition in path_definitions
        id = path_id(definition).name
        camera_id = Symbol(id, :_camera)
        sample_period = id === :unequal ? 150_000_000 : 100_000_000
        sample_phase = id === :phase_offset ? 50_000_000 : 0
        sample_origin = id === :origin_offset ?
            PlantTimestamp(50_000_000) : zero(PlantTimestamp)
        acquisition_phase =
            id in (:phase_offset, :origin_offset) ? 60_000_000 : 10_000_000
        push!(
            sample_definitions,
            OpticalSampleDefinition(
                id,
                PeriodicSchedule(
                    period_ns=sample_period,
                    phase_ns=sample_phase,
                ),
                origin=sample_origin,
            ),
        )
        push!(
            acquisition_definitions,
            AcquisitionDefinition(
                camera_id,
                id,
                DeviceBatchTestAcquisitionModel(T(0.08)),
            ),
        )
        push!(
            event_definitions,
            AcquisitionEventDefinition(
                camera_id,
                GlobalShutterAcquisitionDefinition(
                    PlantDuration(80_000_000);
                    readout_duration=PlantDuration(10_000_000),
                    readiness_delay=PlantDuration(5_000_000),
                ),
                PeriodicAcquisitionStart(
                    PeriodicSchedule(
                        period_ns=200_000_000,
                        phase_ns=acquisition_phase,
                    ),
                ),
            ),
        )
    end

    physical = include_physical_state ?
        device_batch_test_physical_definitions(telescope, backend, T) :
        nothing
    sampled_aberrations = isnothing(physical) ?
        () : physical.sampled_aberrations
    controllable_optics = isnothing(physical) ? () : (physical.optic,)
    command_endpoints = isnothing(physical) ? () : (physical.configuration,)
    definition = PlantDefinition(
        ;
        telescope,
        atmosphere,
        sampled_aberrations,
        controllable_optics,
        paths=Tuple(path_definitions),
        acquisitions=Tuple(acquisition_definitions),
    )
    plant = prepare_plant(
        definition;
        run_seed=0x7_400,
        command_endpoints,
    )
    event_definition = PlantEventLoopDefinition(
        Tuple(sample_definitions),
        Tuple(event_definitions),
    )
    prepared = prepare_device_batch_test_event_loop(
        plant,
        event_definition,
        selection,
    )
    return (
        plant,
        prepared,
        state=PlantEventLoopState(prepared),
        workspace=PlantEventLoopWorkspace(prepared),
        event_definition,
        path_ids=Tuple(path_id(definition).name for definition in
            path_definitions),
        acquisition_ids=Tuple(
            acquisition_id(definition).name for definition in
            acquisition_definitions
        ),
        command_schema=isnothing(physical) ? nothing : physical.schema,
    )
end

function submit_device_batch_test_command!(
    fixture;
    sequence::Integer=1,
    effective_timestamp::PlantTimestamp=PlantTimestamp(200_000_000),
)
    schema = fixture.command_schema
    isnothing(schema) && error(
        "device-batch test fixture has no physical command state")
    current = effective_command(
        fixture.prepared,
        fixture.state,
        :device_batch_dm_command,
    )
    T = eltype(current)
    payload = similar(current)
    copyto!(payload, T[-5e-10, 2e-9, 0, 1e-9])
    return admit_plant_command!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantCommand(schema, sequence, effective_timestamp, payload),
        PlantTimestamp(0),
    )
end

function compare_device_batch_test_command_state(first, second)
    first_effective = effective_command(
        first.prepared,
        first.state,
        :device_batch_dm_command,
    )
    second_effective = effective_command(
        second.prepared,
        second.state,
        :device_batch_dm_command,
    )
    @test Array(second_effective) == Array(first_effective)
    @test last_command_application_timestamp(
        second.state.command_applications[1],
    ) == last_command_application_timestamp(
        first.state.command_applications[1],
    )

    first_surface =
        surface_opd(only(first.state.controllable_optics).active)
    second_surface =
        surface_opd(only(second.state.controllable_optics).active)
    @test Array(second_surface) == Array(first_surface)

    @test command_disposition_count(second.workspace) ==
        command_disposition_count(first.workspace)
    for index in 1:command_disposition_count(first.workspace)
        first_disposition = command_disposition(first.workspace, index)
        second_disposition = command_disposition(second.workspace, index)
        @test command_sequence(second_disposition) ==
            command_sequence(first_disposition)
        @test command_terminal_kind(second_disposition) ==
            command_terminal_kind(first_disposition)
        @test command_disposition_reason(second_disposition) ==
            command_disposition_reason(first_disposition)
        @test command_requested_effective_timestamp(second_disposition) ==
            command_requested_effective_timestamp(first_disposition)
        @test command_terminal_timestamp(second_disposition) ==
            command_terminal_timestamp(first_disposition)
        @test command_lateness(second_disposition) ==
            command_lateness(first_disposition)
    end
    return nothing
end
