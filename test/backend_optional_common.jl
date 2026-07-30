backend_package_name(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag}) = "CUDA"
backend_package_name(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag}) = "AMDGPU"

backend_label(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag}) = "CUDA"
backend_label(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag}) = "AMDGPU"

backend_full_smoke_env(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_CUDA"
backend_full_smoke_env(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_AMDGPU"

if !isdefined(@__MODULE__, :ContractRateModel)
    include(joinpath(@__DIR__, "wfs_stage_contract_fixtures.jl"))
end

backend_selector(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag}) = AdaptiveOpticsSim.Backends.CUDABackend()
backend_selector(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag}) = AdaptiveOpticsSim.Backends.AMDGPUBackend()

function run_optional_command_application_checks(::Type{B}, BackendArray) where {
    B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    schema = PlantCommandSchema(
        T,
        (3,);
        id=:optional_command_application_schema,
        version=1,
        endpoint=:optional_command_application,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :optional_command_application),
        basis_revision=1,
        semantics=IncrementalCommand,
        bounds=UniformCommandBounds(T(-1), T(1)),
        value_policy=CommandValuePolicy(
            out_of_range=ClipInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    endpoint = prepare_command_endpoint(schema;
        capacity=1, ordinal=1, backend=selector)
    endpoint_state = Plant.CommandEndpointState(endpoint)
    initial = BackendArray(zeros(T, 3))
    safe = BackendArray(fill(T(-0.25), 3))
    application_state = Plant.CommandApplicationState(
        endpoint, endpoint_state, initial; safe_command=safe)
    workspace = Plant.CommandDispositionWorkspace(endpoint)
    payload = BackendArray(T[2, 0.25, -2])
    admit_plant_command!(workspace, endpoint, endpoint_state,
        PlantCommand(schema, 1, PlantTimestamp(1), payload),
        PlantTimestamp(1))
    claim = claim_next_application_ready_command!(endpoint, endpoint_state,
        PlantTimestamp(1))
    apply_claimed_plant_command!(workspace, endpoint, endpoint_state,
        application_state, claim)
    @test effective_command(application_state) isa BackendArray
    @test Array(effective_command(application_state)) == T[1, 0.25, -1]
    clear_command_dispositions!(workspace)
    transition = apply_command_silence_transition!(workspace, endpoint,
        endpoint_state, application_state, PlantTimestamp(11))
    @test command_silence_action(transition) == ApplySafeCommand
    @test Array(effective_command(application_state)) == fill(T(-0.25), 3)
    return nothing
end

struct OptionalControllerRoutingModel end
struct OptionalPreparedControllerRoutingModel end

Plant.plant_model_definition_style(
    ::Type{OptionalControllerRoutingModel}) = ColdPlantModelDefinition()

Plant.prepare_controllable_optic(
    ::OptionalControllerRoutingModel,
    ::ControllableOpticDefinition,
    ::Telescope,
    ::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
) = OptionalPreparedControllerRoutingModel()

function run_optional_controller_routing_checks(::Type{B},
    BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T, backend=selector)
    atmosphere = OptionalStaticAtmosphere(telescope; T, backend=selector)
    schema = PlantCommandSchema(
        T,
        (2,);
        id=:optional_controller_routing_schema,
        version=1,
        endpoint=:optional_controller_routing,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :optional_controller_routing),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
    definition = PlantDefinition(; telescope, atmosphere,
        controllable_optics=(
            ControllableOpticDefinition(:optional_controller_routing,
                OptionalControllerRoutingModel(), (schema,);
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility()),))
    plant = prepare_plant(definition;
        run_seed=0x6101,
        command_endpoints=(
            CommandEndpointConfiguration(:optional_controller_routing,
                zeros(T, 2); capacity=2, backend=selector),))

    flat_output = BackendArray(T[0, 1, 2, 0])
    routed_output = @view flat_output[2:3]
    routing = prepare_controller_output_routing(
        plant,
        (command=routed_output,),
        ControllerOutputRoute(:command, :optional_controller_routing),
    )
    @test controller_output_payload(
        routing, Val(:optional_controller_routing)) === routed_output
    @test Array(controller_output_payload(
        routing, :optional_controller_routing)) == T[1, 2]

    error = try
        prepare_controller_output_routing(
            plant,
            (command=zeros(T, 2),),
            ControllerOutputRoute(:command, :optional_controller_routing),
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    error isa PlantPreparationError &&
        @test error.reason === :backend
    return nothing
end

function run_optional_cycle_averaged_modulation_checks(::Type{B},
    BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    policy = CircularModulation(T(2);
        samples=5, phase_offset=T(0.3), T=T)
    cpu = AdaptiveOpticsSim.Optics.prepare_focal_plane_modulation(
        policy, 8, zeros(T, 8, 8), T)
    device = AdaptiveOpticsSim.Optics.prepare_focal_plane_modulation(
        policy, 8, BackendArray(zeros(T, 8, 8)), T)
    weights = copy(device.amplitude_weights)

    AdaptiveOpticsSim.Optics.update_cycle_averaged_circular_modulation!(
        cpu, T(1.5))
    AdaptiveOpticsSim.Optics.update_cycle_averaged_circular_modulation!(
        device, T(1.5))
    @test device.phases isa BackendArray
    @test device.amplitude_weights == weights
    @test isapprox(sum(abs2, device.amplitude_weights), one(T);
        rtol=8eps(T))
    @test isapprox(Array(device.phases), cpu.phases;
        rtol=8eps(T), atol=8eps(T))

    AdaptiveOpticsSim.Optics.update_cycle_averaged_circular_modulation!(
        device, T(1.5); enabled=false)
    @test Array(device.phases) == ones(Complex{T}, 8, 8, 5)
    @test device.amplitude_weights == weights
    return nothing
end

struct OptionalStaticAtmosphere{A,B<:AbstractArrayBackend} <: AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere
    screen::A
end

struct OptionalPreparedDirectPathModel
    zero_padding::Int
end

struct OptionalPreparedIlluminationPathModel{T<:AbstractFloat}
    photon_rate::T
end

struct OptionalPreparedFrameAcquisitionModel{T<:AbstractFloat}
    exposure::T
end

struct OptionalIlluminationIdentityExecution{P}
    product::P
end

Plant.plant_model_definition_style(
    ::Type{OptionalPreparedDirectPathModel}) = ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{<:OptionalPreparedIlluminationPathModel}) =
    ColdPlantModelDefinition()

Plant.plant_model_definition_style(
    ::Type{<:OptionalPreparedFrameAcquisitionModel}) =
    ColdPlantModelDefinition()

function Plant.validate_path_execution_binding(
    execution::OptionalIlluminationIdentityExecution, input, result)
    execution.product === input && input === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "optional illumination identity binding changed"))
    return nothing
end

function Plant.execute_path!(result, input,
    execution::OptionalIlluminationIdentityExecution)
    Plant.validate_path_execution_binding(execution, input,
        result)
    return result
end

function Plant.prepare_path_executor(
    model::OptionalPreparedDirectPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    return Plant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        materialization=optional_path_materialization(atmosphere, telescope,
            source, pupil),
        optical_model=(kind=:direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_path_executor(
    model::OptionalPreparedIlluminationPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    values = similar(pupil_reflectivity(telescope), T,
        telescope.params.resolution, telescope.params.resolution)
    fill!(values, zero(T))
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=MetricCoordinates(),
        sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(source))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    destination = IntensityMap(metadata, values)
    entry = prepare_illumination_entry(
        UniformIntensityIllumination(model.photon_rate;
            combination=SingleIllumination()),
        destination, DetectorInputIlluminationEntry();
        visibility=(downstream_path=path_id(definition),
            starts_at=:detector_input))
    execution = OptionalIlluminationIdentityExecution(destination)
    return Plant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        destination,
        destination,
        execution;
        materialization=entry,
        optical_model=(kind=:uniform_detector_illumination,
            photon_rate=T(model.photon_rate)),
        propagation_model=:detector_input_identity,
        model_revisions=UInt(1),
    )
end

@inline optional_path_materialization(
    ::OptionalStaticAtmosphere,
    ::Telescope,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
    ::PupilFunction,
) = Plant.AtmosphereIndependentPath()

@inline optional_path_materialization(
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
) = prepare_pupil_opd_materialization(atmosphere, telescope, source, pupil)

function Plant.prepare_acquisition_provider(
    model::OptionalPreparedFrameAcquisitionModel,
    definition::AcquisitionDefinition,
    path::Plant.PreparedPathExecutor,
)
    Plant.require_path_result(path)
    result = path_result(path)
    T = eltype(result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), response_model=NullFrameResponse(),
        T=T, backend=path_result_key(path).backend)
    execution = Plant.FrameAcquisitionExecution(detector, result)
    metadata = (kind=:detector_frame, units=:detected_electrons,
        geometry=result.metadata,
        detector=detector_export_metadata(detector))
    products = Plant.AcquisitionProducts(execution.observation;
        metadata)
    return prepare_full_optical_provider(execution, products)
end

AdaptiveOpticsSim.Backends.backend(::OptionalStaticAtmosphere{<:Any,B}) where {B} = B()
AdaptiveOpticsSim.Atmospheres.advance!(atm::OptionalStaticAtmosphere, tel::Telescope, rng::AbstractRNG) = atm
AdaptiveOpticsSim.Atmospheres.advance!(atm::OptionalStaticAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng()) = atm

function optional_detector_calibration_signature_allocation_bytes(
    det, seed::UInt)
    AdaptiveOpticsSim.detector_calibration_signature(det, seed)
    return @allocated AdaptiveOpticsSim.detector_calibration_signature(
        det, seed)
end

function AdaptiveOpticsSim.Atmospheres.propagate!(atm::OptionalStaticAtmosphere,
    pupil::PupilFunction)
    copyto!(pupil.opd, atm.screen)
    return pupil
end

function OptionalStaticAtmosphere(tel::Telescope; T::Type{<:AbstractFloat}=Float32, backend::AbstractArrayBackend=backend(tel))
    selector = AdaptiveOpticsSim.Backends.require_same_backend(tel, AdaptiveOpticsSim.Backends._resolve_backend_selector(backend))
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(selector)
    host = zeros(T, tel.params.resolution, tel.params.resolution)
    host .*= Array(pupil_mask(tel))
    screen = array_backend{T}(undef, size(host)...)
    copyto!(screen, host)
    return OptionalStaticAtmosphere{typeof(screen),typeof(selector)}(screen)
end

function run_optional_prepared_plant_checks(::Type{B},
    BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T, backend=selector)
    atmosphere = MultiLayerAtmosphere(telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8, 4],
        wind_direction=T[0, 90],
        altitude=T[0, 5_000],
        layer_ids=(:ground, :high),
        T=T,
        backend=selector,
    )
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(3), T=T)
    sampled_prototype =
        PupilFunction(telescope; T=T, backend=selector)
    sampled_opd = similar(sampled_prototype.opd)
    fill!(sampled_opd, T(0.125))
    sampled_metadata = OpticalPlaneMetadata(
        PupilPlane(), sampled_opd;
        coordinate_domain=MetricCoordinates(),
        sampling=sampled_prototype.metadata.sampling,
        origin=sampled_prototype.metadata.origin,
        centering=sampled_prototype.metadata.centering,
        orientation=sampled_prototype.metadata.orientation,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
    )
    sampled_aberration = SampledAberrationDefinition(
        :optional_science_static,
        OPDMap(sampled_opd),
        sampled_metadata;
        placement=PupilPlanePlacement(),
        visibility=SelectedPathVisibility(:science),
        application=DMAdditive(),
    )
    path_definition = OpticalPathDefinition(:science, source,
        OptionalPreparedDirectPathModel(2))
    illumination_path_definition = OpticalPathDefinition(:illumination,
        source, OptionalPreparedIlluminationPathModel(T(8)))
    fast_definition = AcquisitionDefinition(:fast_science, :science,
        OptionalPreparedFrameAcquisitionModel(T(0.25)))
    slow_definition = AcquisitionDefinition(:slow_science, :science,
        OptionalPreparedFrameAcquisitionModel(T(0.75)))
    illumination_definition = AcquisitionDefinition(:illumination_frame,
        :illumination, OptionalPreparedFrameAcquisitionModel(T(0.125)))
    definition = PlantDefinition(; telescope, atmosphere,
        sampled_aberrations=(sampled_aberration,),
        paths=(science=path_definition,
            illumination=illumination_path_definition),
        acquisitions=(fast_science=fast_definition,
            slow_science=slow_definition,
            illumination_frame=illumination_definition))

    plant = prepare_plant(definition; run_seed=0x6100)
    prepared_sampled = Plant.prepared_sampled_aberration(
        plant, :optional_science_static)
    prepared_sampled_opd =
        Plant.sampled_aberration_opd(prepared_sampled)
    @test prepared_sampled_opd isa BackendArray
    @test compute_device(prepared_sampled_opd) ==
        sampled_metadata.device
    fill!(sampled_opd, zero(T))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(prepared_sampled_opd))
    @test all(==(T(0.125)), Array(prepared_sampled_opd))
    path = prepared_path(plant, :science)
    illumination_path = prepared_path(plant, :illumination)
    fast = prepared_acquisition(plant, :fast_science)
    slow = prepared_acquisition(plant, :slow_science)
    illumination = prepared_acquisition(plant, :illumination_frame)
    @test path_result(path).values isa BackendArray
    @test path_input(illumination_path).values isa BackendArray
    @test path_input(illumination_path) === path_result(illumination_path)
    @test acquisition_observation(fast) isa BackendArray
    @test acquisition_observation(slow) isa BackendArray
    @test acquisition_observation(illumination) isa BackendArray
    @test path_result_key(path).device ==
        compute_device(acquisition_observation(fast))
    @test acquisition_observation(fast) !== acquisition_observation(slow)

    pupil = path_input(path)
    surface_host = zeros(T, size(pupil.opd))
    surface_host[4, 4] = one(T)
    surface = similar(pupil.opd)
    copyto!(surface, surface_host)
    identity_coupling = prepare_sampled_pupil_footprint_coupling(
        pupil.metadata,
        surface,
        path,
        PupilPlanePlacement(),
    )
    @test identity_coupling isa
        PreparedIdentityPupilFootprintCoupling
    fill!(pupil.opd, zero(T))
    @test apply_sampled_pupil_surface!(
        pupil, surface, identity_coupling) === pupil
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(pupil.opd))
    @test Array(pupil.opd) == surface_host
    fill!(pupil.opd, T(7))
    @test apply_sampled_pupil_surface!(
        pupil, surface, identity_coupling, DMReplace()) === pupil
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(pupil.opd))
    @test Array(pupil.opd) == surface_host

    half_sample = (
        pupil.metadata.sampling[1] / T(2),
        pupil.metadata.sampling[2] / T(2),
    )
    transformed_coupling = prepare_sampled_pupil_footprint_coupling(
        pupil.metadata,
        surface,
        path,
        PupilPlanePlacement();
        registration=PupilRelayRegistration(
            decenter_m=half_sample,
            T=T,
        ),
    )
    @test transformed_coupling isa PreparedPupilFootprintCoupling
    fill!(pupil.opd, zero(T))
    @test apply_sampled_pupil_surface!(
        pupil, surface, transformed_coupling) === pupil
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(pupil.opd))
    transformed_expected = zeros(T, size(surface_host))
    transformed_expected[3:4, 3:4] .= T(0.25)
    @test Array(pupil.opd) ≈ transformed_expected rtol=zero(T) atol=8eps(T)
    fill!(pupil.opd, T(7))
    @test apply_sampled_pupil_surface!(
        pupil, surface, transformed_coupling, DMReplace()) === pupil
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(pupil.opd))
    @test Array(pupil.opd) ≈ transformed_expected rtol=zero(T) atol=8eps(T)

    selection = prepare_acquisition_selection(plant,
        (:slow_science, :illumination_frame, :fast_science))
    @test @inferred(execute_acquisition_selection_at!(selection,
        T(1e-3))) === selection
    fast_products = acquisition_products(fast)
    slow_products = acquisition_products(slow)
    @test epoch_time(current_epoch(atmosphere)) == T(1e-3)
    @test path_input(path).opd isa BackendArray
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(path_result(path).values))
    @test fast_products === acquisition_products(fast)
    @test slow_products === acquisition_products(slow)
    fast_host = Array(acquisition_observation(fast))
    slow_host = Array(acquisition_observation(slow))
    illumination_host = Array(acquisition_observation(illumination))
    @test all(isfinite, fast_host)
    @test sum(fast_host) > zero(T)
    @test slow_host ≈ T(3) .* fast_host rtol=T(3e-5) atol=T(3e-5)
    @test all(==(one(T)), illumination_host)
    @test all(==(T(8)), Array(path_result(illumination_path).values))
    @test illumination_entry_boundary(
        Plant.path_materialization(illumination_path)) isa
        DetectorInputIlluminationEntry

    synthetic_metadata = acquisition_product_metadata(fast)
    copied_destination = Plant.AcquisitionProducts(
        similar(acquisition_observation(fast));
        metadata=synthetic_metadata)
    copied_source = Plant.AcquisitionProducts(
        similar(acquisition_observation(fast));
        metadata=synthetic_metadata)
    fill!(copied_source.observation, T(9))
    copied_provider = prepare_copied_synthetic_provider(
        copied_destination, copied_source)
    @test Plant.execute_acquisition_provider!(copied_provider,
        path_result(path), Xoshiro(0x6110)) === copied_destination
    replay_destination = Plant.AcquisitionProducts(
        similar(acquisition_observation(fast));
        metadata=synthetic_metadata)
    replay_first = Plant.AcquisitionProducts(
        similar(acquisition_observation(fast));
        metadata=synthetic_metadata)
    replay_second = Plant.AcquisitionProducts(
        similar(acquisition_observation(fast));
        metadata=synthetic_metadata)
    fill!(replay_first.observation, T(10))
    fill!(replay_second.observation, T(11))
    replay_provider = prepare_cyclic_replay_provider(replay_destination,
        (replay_first, replay_second))
    @test Plant.execute_acquisition_provider!(replay_provider,
        path_result(path), Xoshiro(0x6111)) === replay_destination
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(copied_destination.observation))
    @test all(==(T(9)), Array(copied_destination.observation))
    @test all(==(T(10)), Array(replay_destination.observation))
    Plant.execute_acquisition_provider!(replay_provider,
        path_result(path), Xoshiro(0x6112))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(replay_destination.observation))
    @test all(==(T(11)), Array(replay_destination.observation))

    @test_throws PlantPreparationError Plant.require_path_result(
        path; backend=CPUBackend())
    @test_throws PlantPreparationError Plant.require_path_result(
        path; device=AdaptiveOpticsSim.Backends.HostComputeDevice())
    return nothing
end

function optional_device_path_batch_allocation_bytes(fixture, owner)
    run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(199_000_000),
    )
    claim = Plant.begin_optical_path_batch!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(200_000_000),
    )
    materialization_bytes = @allocated materialize_device_path_batch!(
        owner,
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    Plant.seal_optical_path_batch_materialization!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    execution_bytes = @allocated execute_device_path_batch!(
        owner,
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    Plant.complete_optical_path_batch!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    return (; materialization_bytes, execution_bytes)
end

function run_optional_device_path_batch_checks(
    ::Type{B},
    BackendArray,
) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    lifecycle = device_batch_test_fixture(
        backend=selector,
        T=Float32,
        include_unequal_rate=false,
        include_lgs=false,
    )
    @test device_path_batch_owner_count(lifecycle.prepared) == 1
    owner = device_path_batch_owner(lifecycle.prepared, 1)
    device = device_path_batch_compute_device(owner)
    @test device isa AdaptiveOpticsSim.Backends.AcceleratorComputeDevice
    @test typeof(device_path_batch_backend(owner)) === typeof(selector)
    @test device_path_batch_group_count(owner) == 2
    @test all(
        index -> path_execution_group_device_batch_owner_ordinal(
            lifecycle.prepared,
            device_path_batch_group_ordinal(owner, index),
        ) == 1,
        1:device_path_batch_group_count(owner),
    )

    implementation = owner.implementation
    retained_context = implementation.context
    retained_fft_plan = implementation.optical_batch.workspace.fft_plan
    @test AdaptiveOpticsSim.Backends._prepared_device_execution_compute_device(
        implementation.context,
    ) == device
    @test atmosphere_direction_output(
        implementation.atmosphere_batch,
    ) isa BackendArray
    @test implementation.optical_batch.workspace.field_stack isa BackendArray
    @test implementation.optical_batch.workspace.output_stack isa BackendArray
    @test all(input -> input.opd isa BackendArray,
        implementation.path_inputs)
    @test all(result -> result.values isa BackendArray,
        implementation.path_results)
    @test effective_command(
        lifecycle.prepared,
        lifecycle.state,
        :device_batch_dm_command,
    ) isa BackendArray
    @test surface_opd(
        only(lifecycle.state.controllable_optics).active,
    ) isa BackendArray
    @test Plant.sampled_aberration_opd(
        Plant.prepared_sampled_aberration(
            lifecycle.plant,
            :device_batch_static,
        ),
    ) isa BackendArray
    @test Plant.sampled_aberration_opd(
        Plant.prepared_sampled_aberration(
            lifecycle.plant,
            :device_batch_alpha_ncpa,
        ),
    ) isa BackendArray
    @test all(
        storage -> compute_device(storage) == device,
        (
            atmosphere_direction_output(implementation.atmosphere_batch),
            implementation.optical_batch.workspace.field_stack,
            implementation.optical_batch.workspace.output_stack,
        ),
    )

    allocation_bytes =
        optional_device_path_batch_allocation_bytes(lifecycle, owner)
    # Backend runtimes retain the data plane on the device but may allocate
    # bounded host launch descriptors for layer, FFT, copy, and detector
    # submissions.
    @test allocation_bytes.materialization_bytes <= 256 * 1024
    @test allocation_bytes.execution_bytes <= 256 * 1024
    @test owner.implementation.context === retained_context
    @test owner.implementation.optical_batch.workspace.fft_plan ===
        retained_fft_plan
    @test all(
        result -> all(isfinite, Array(result.values)),
        implementation.path_results,
    )

    same_epoch_oracle = device_batch_test_fixture(
        backend=selector,
        T=Float32,
        selection=Val(:none),
    )
    same_epoch_batched = prepare_device_batch_test_event_loop(
        same_epoch_oracle.plant,
        same_epoch_oracle.event_definition,
        Val(:public),
    )
    same_epoch_batched_state = PlantEventLoopState(same_epoch_batched)
    same_epoch_batched_workspace =
        PlantEventLoopWorkspace(same_epoch_batched)
    @test step_plant_events!(
        same_epoch_oracle.prepared,
        same_epoch_oracle.state,
        same_epoch_oracle.workspace,
    ) == PlantTimestamp(0)
    expected_inputs = [
        Array(prepared_path(same_epoch_oracle.plant, id).input.opd)
        for id in same_epoch_oracle.path_ids
    ]
    expected_results = [
        Array(prepared_path(same_epoch_oracle.plant, id).result.values)
        for id in same_epoch_oracle.path_ids
    ]
    @test step_plant_events!(
        same_epoch_batched,
        same_epoch_batched_state,
        same_epoch_batched_workspace,
    ) == PlantTimestamp(0)
    for (index, id) in pairs(same_epoch_oracle.path_ids)
        path = prepared_path(same_epoch_oracle.plant, id)
        @test isapprox(
            Array(path.input.opd),
            expected_inputs[index];
            rtol=2f-5,
            atol=2f-6,
        )
        @test isapprox(
            Array(path.result.values),
            expected_results[index];
            rtol=2f-4,
            atol=2f-4,
        )
    end

    independent = device_batch_test_fixture(
        backend=selector,
        T=Float32,
        selection=Val(:none),
    )
    gpu = device_batch_test_fixture(backend=selector, T=Float32)
    @test command_admission_status(
        submit_device_batch_test_command!(independent),
    ) == CommandAdmittedPending
    @test command_admission_status(
        submit_device_batch_test_command!(gpu),
    ) == CommandAdmittedPending
    horizon = PlantTimestamp(450_000_000)
    @test run_plant_events_until!(
        gpu.prepared,
        gpu.state,
        gpu.workspace,
        horizon,
    ) == run_plant_events_until!(
        independent.prepared,
        independent.state,
        independent.workspace,
        horizon,
    )
    @test scheduler_timestamp(gpu.state.scheduler) ==
        scheduler_timestamp(independent.state.scheduler)
    @test gpu.state.scheduler.revision ==
        independent.state.scheduler.revision
    @test gpu.state.scheduler.cursors ==
        independent.state.scheduler.cursors
    @test gpu.state.path_sampled == independent.state.path_sampled
    @test gpu.state.product_sequences ==
        independent.state.product_sequences
    @test gpu.state.product_ready_timestamps ==
        independent.state.product_ready_timestamps
    @test epoch_time(current_epoch(gpu.prepared.atmosphere)) ==
        epoch_time(current_epoch(independent.prepared.atmosphere))
    @test epoch_sequence(current_epoch(gpu.prepared.atmosphere)) ==
        epoch_sequence(current_epoch(independent.prepared.atmosphere))
    @test length(gpu.prepared.atmosphere_rng.streams) ==
        length(independent.prepared.atmosphere_rng.streams)
    for index in eachindex(gpu.prepared.atmosphere_rng.streams)
        @test copy(gpu.prepared.atmosphere_rng.streams[index].state) ==
            copy(independent.prepared.atmosphere_rng.streams[index].state)
    end
    compare_device_batch_test_command_state(independent, gpu)
    for id in independent.path_ids
        independent_path = prepared_path(independent.plant, id)
        gpu_path = prepared_path(gpu.plant, id)
        @test isapprox(
            sum(Array(gpu_path.result.values)),
            sum(Array(independent_path.result.values));
            rtol=3f-4,
            atol=3f-4,
        )
    end
    for id in independent.acquisition_ids
        @test acquisition_product_sequence(
            gpu.prepared,
            gpu.state,
            id,
        ) == acquisition_product_sequence(
            independent.prepared,
            independent.state,
            id,
        )
        @test acquisition_product_ready_timestamp(
            gpu.prepared,
            gpu.state,
            id,
        ) == acquisition_product_ready_timestamp(
            independent.prepared,
            independent.state,
            id,
        )
        @test isapprox(
            sum(Array(acquisition_observation(
                prepared_acquisition(gpu.plant, id),
            ))),
            sum(Array(acquisition_observation(
                prepared_acquisition(independent.plant, id),
            )));
            rtol=3f-4,
            atol=3f-4,
        )
    end
    return nothing
end

@inline function optional_wfs_plan_device_resident(
    plan::WavefrontSensors.PreparedShackHartmannOpticalFormation,
    device,
    BackendArray,
)
    propagation = plan.front_end.propagation
    arrays = (
        propagation.field,
        propagation.phasor,
        propagation.fft_stack,
        propagation.intensity_stack,
        propagation.sampled_spot_cube,
        propagation.spot_cube_accum,
        propagation.fft_asterism_stack,
    )
    return all(
        array -> array isa BackendArray &&
            compute_device(array) == device,
        arrays,
    )
end

@inline function optional_wfs_plan_device_resident(
    plan::Union{
        WavefrontSensors.PreparedPyramidOpticalFormation,
        WavefrontSensors.PreparedBioEdgeOpticalFormation,
    },
    device,
    BackendArray,
)
    front_end = plan.front_end
    propagation = front_end.propagation
    arrays = (
        front_end.modulation.phases,
        propagation.field,
        propagation.focal_field,
        propagation.pupil_field,
        propagation.phasor,
        propagation.intensity,
        propagation.temp,
        propagation.scratch,
        propagation.asterism_stack,
    )
    return all(
        array -> array isa BackendArray &&
            compute_device(array) == device,
        arrays,
    )
end

function optional_wfs_plan_device_resident(
    plan::Union{
        WavefrontSensors.PreparedShackHartmannOpticalBundleFormation,
        WavefrontSensors.PreparedPyramidOpticalBundleFormation,
        WavefrontSensors.PreparedBioEdgeOpticalBundleFormation,
    },
    device,
    BackendArray,
)
    return all(
        component -> optional_wfs_plan_device_resident(
            component,
            device,
            BackendArray,
        ),
        plan.plans,
    )
end

@inline function optional_wfs_result_device_resident(
    result::IntensityMap,
    device,
    BackendArray,
)
    return result.values isa BackendArray &&
        compute_device(result.values) == device
end

function optional_wfs_result_device_resident(
    result::OpticalProductBundle,
    device,
    BackendArray,
)
    return all(
        index -> optional_wfs_result_device_resident(
            result[index],
            device,
            BackendArray,
        ),
        1:length(result),
    )
end

@inline optional_wfs_products_approx(
    device::AbstractArray,
    oracle::AbstractArray,
) = isapprox(device, oracle; rtol=4f-4, atol=8f-5)

function optional_wfs_products_approx(
    device::Tuple,
    oracle::Tuple,
)
    length(device) == length(oracle) || return false
    return all(
        optional_wfs_products_approx(device[index], oracle[index])
        for index in eachindex(device)
    )
end

function optional_wfs_device_path_batch_allocation_bytes(
    fixture,
    ::PreparedDevicePathBatchOwner,
)
    run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(199_000_000),
    )
    completed_timestamp_bytes = @allocated run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(200_000_000),
    )
    return (; completed_timestamp_bytes)
end

function run_optional_wfs_device_model_matrix_checks(
    ::Type{B},
    BackendArray,
) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    for (family, direction, spectral) in device_model_matrix_wfs_rows()
        oracle = device_model_matrix_wfs_fixture(
            family;
            backend=CPUBackend(),
            selection=Val(:none),
            direction,
            spectral,
            T=Float32,
            r0=2f6,
        )
        device_fixture = device_model_matrix_wfs_fixture(
            family;
            backend=selector,
            selection=Val(:all),
            direction,
            spectral,
            T=Float32,
            r0=2f6,
        )
        @test device_path_batch_owner_count(device_fixture.prepared) == 1
        owner = device_path_batch_owner(device_fixture.prepared, 1)
        implementation = owner.implementation
        @test implementation isa Plant._PreparedWFSDevicePathBatch
        @test device_path_batch_group_count(owner) == 2
        device = device_path_batch_compute_device(owner)
        @test device isa AdaptiveOpticsSim.Backends.AcceleratorComputeDevice
        @test typeof(device_path_batch_backend(owner)) === typeof(selector)
        @test atmosphere_direction_output(
            implementation.atmosphere_batch,
        ) isa BackendArray
        @test compute_device(
            atmosphere_direction_output(implementation.atmosphere_batch),
        ) == device
        @test all(
            input -> input.opd isa BackendArray &&
                compute_device(input.opd) == device,
            implementation.path_inputs,
        )
        @test all(
            result -> optional_wfs_result_device_resident(
                result,
                device,
                BackendArray,
            ),
            implementation.path_results,
        )
        retained_context = implementation.context
        retained_plans = ntuple(2) do index
            group_slot = Int(owner.group_slots[index])
            plan = device_fixture.prepared.path_groups[
                group_slot].path.execution.plan
            @test optional_wfs_plan_device_resident(
                plan,
                device,
                BackendArray,
            )
            return plan
        end
        command = effective_command(
            device_fixture.prepared,
            device_fixture.state,
            :device_batch_dm_command,
        )
        optic_surface = surface_opd(
            only(device_fixture.state.controllable_optics).active,
        )
        common_aberration = Plant.sampled_aberration_opd(
            Plant.prepared_sampled_aberration(
                device_fixture.plant,
                :device_batch_static,
            ),
        )
        selected_ncpa = Plant.sampled_aberration_opd(
            Plant.prepared_sampled_aberration(
                device_fixture.plant,
                :device_batch_alpha_ncpa,
            ),
        )
        @test all(
            storage -> storage isa BackendArray &&
                compute_device(storage) == device,
            (
                command,
                optic_surface,
                common_aberration,
                selected_ncpa,
            ),
        )
        @test command_admission_status(
            submit_device_batch_test_command!(oracle),
        ) == CommandAdmittedPending
        @test command_admission_status(
            submit_device_batch_test_command!(device_fixture),
        ) == CommandAdmittedPending

        @test step_plant_events!(
            device_fixture.prepared,
            device_fixture.state,
            device_fixture.workspace,
        ) == step_plant_events!(
            oracle.prepared,
            oracle.state,
            oracle.workspace,
        ) == PlantTimestamp(0)
        device_model_matrix_copy_atmosphere_screens!(
            device_fixture.prepared.atmosphere,
            oracle.prepared.atmosphere,
        )
        horizon = PlantTimestamp(200_000_000)
        @test run_plant_events_until!(
            device_fixture.prepared,
            device_fixture.state,
            device_fixture.workspace,
            horizon,
        ) == run_plant_events_until!(
            oracle.prepared,
            oracle.state,
            oracle.workspace,
            horizon,
        )
        @test device_fixture.state.scheduler.cursors ==
            oracle.state.scheduler.cursors
        @test device_fixture.state.path_sampled ==
            oracle.state.path_sampled
        @test device_fixture.state.product_sequences ==
            oracle.state.product_sequences
        @test device_fixture.state.product_ready_timestamps ==
            oracle.state.product_ready_timestamps
        @test epoch_time(current_epoch(
            device_fixture.prepared.atmosphere,
        )) == epoch_time(current_epoch(oracle.prepared.atmosphere))
        @test epoch_sequence(current_epoch(
            device_fixture.prepared.atmosphere,
        )) == epoch_sequence(current_epoch(oracle.prepared.atmosphere))
        @test implementation.context === retained_context
        current_plans = ntuple(2) do index
            group_slot = Int(owner.group_slots[index])
            device_fixture.prepared.path_groups[
                group_slot].path.execution.plan
        end
        @test current_plans === retained_plans
        compare_device_batch_test_command_state(
            oracle,
            device_fixture;
            surface_rtol=8eps(Float32),
            surface_atol=eps(Float32) * 1f-9,
        )
        for id in oracle.path_ids
            oracle_path = prepared_path(oracle.plant, id)
            device_path = prepared_path(device_fixture.plant, id)
            @test isapprox(
                Array(device_path.input.opd),
                Array(oracle_path.input.opd);
                rtol=3f-5,
                atol=4f-6,
            )
            @test optional_wfs_products_approx(
                device_model_matrix_product_host(device_path.result),
                device_model_matrix_product_host(oracle_path.result),
            )
        end

        detector_row =
            device_model_matrix_wfs_detector_row(family, direction)
        oracle_path = prepared_path(oracle.plant, :wfs_alpha)
        device_path = prepared_path(device_fixture.plant, :wfs_alpha)
        oracle_detector_input =
            device_model_matrix_detector_facing_product(oracle_path.result)
        device_detector_input =
            device_model_matrix_detector_facing_product(device_path.result)
        detector_input_before = Array(device_detector_input.values)
        oracle_detector = device_model_matrix_execute_detector(
            detector_row;
            backend=CPUBackend(),
            T=Float32,
            input_map=oracle_detector_input,
        )
        device_detector = device_model_matrix_execute_detector(
            detector_row;
            backend=selector,
            T=Float32,
            input_map=device_detector_input,
        )
        detector_device = compute_device(device_detector.output)
        @test detector_device == device
        @test device_detector.output isa BackendArray
        @test device_detector.prepared.plan.input_values ===
            device_detector_input.values
        @test optional_detector_response_device_resident(
            device_detector.detector.params.response_model,
            device,
            BackendArray,
        )
        @test optional_detector_state_device_resident(
            device_detector.detector,
            device,
            BackendArray,
        )
        @test optional_detector_prepared_storage_device_resident(
            device_detector.prepared,
            device,
            BackendArray,
        )
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(device_detector.output),
        )
        expected_detector_output =
            device_model_matrix_expected_detector_frame(
                detector_row,
                device_detector.detector,
                detector_input_before,
                Float32,
            )
        @test isapprox(
            Array(device_detector.output),
            expected_detector_output;
            rtol=8f-5,
            atol=8f-5,
        )
        @test isapprox(
            Array(device_detector.output),
            Array(oracle_detector.output);
            rtol=8f-5,
            atol=8f-5,
        )
        @test device_detector.status_trace == oracle_detector.status_trace
        @test device_detector.event_times == oracle_detector.event_times
        @test device_model_matrix_response_metadata_signature(
            device_detector.metadata_before,
        ) == device_model_matrix_response_metadata_signature(
            device_detector.metadata_after,
        )
        @test device_detector.metadata_after.sensor ==
            device_model_matrix_detector_sensor_symbol(detector_row)
        detector_allocation_bytes =
            optional_detector_device_model_matrix_allocation_bytes(
                detector_row,
                device_detector,
            )
        @test detector_allocation_bytes <= 1024 * 1024
        @test Array(device_detector_input.values) == detector_input_before

        allocation_fixture = device_model_matrix_wfs_fixture(
            family;
            backend=selector,
            selection=Val(:all),
            direction,
            spectral,
            T=Float32,
            r0=2f6,
        )
        allocation_owner =
            device_path_batch_owner(allocation_fixture.prepared, 1)
        allocation_bytes =
            optional_wfs_device_path_batch_allocation_bytes(
                allocation_fixture,
                allocation_owner,
            )
        # This bound includes the co-scheduled ordinary detector witness path
        # as well as the two-member WFS owner and its completion barrier.
        @test allocation_bytes.completed_timestamp_bytes <= 1024 * 1024
    end
    return nothing
end

@inline optional_detector_response_device_resident(
    ::NullFrameResponse,
    ::Any,
    ::Any,
) = true

@inline function optional_detector_response_device_resident(
    response::GaussianPixelResponse,
    device,
    BackendArray,
)
    return response.kernel isa BackendArray &&
        compute_device(response.kernel) == device
end

@inline function optional_detector_response_device_resident(
    response::SampledFrameResponse,
    device,
    BackendArray,
)
    return response.kernel isa BackendArray &&
        compute_device(response.kernel) == device
end

@inline function optional_detector_response_device_resident(
    response::RectangularPixelAperture,
    device,
    BackendArray,
)
    return response.kernel_x isa BackendArray &&
        response.kernel_y isa BackendArray &&
        compute_device(response.kernel_x) == device &&
        compute_device(response.kernel_y) == device
end

function optional_detector_readout_products_device_resident(
    products::FrameReadoutProducts,
    device,
    BackendArray,
)
    for name in fieldnames(typeof(products))
        name === :read_times && continue
        value = getfield(products, name)
        if value isa AbstractArray
            value isa BackendArray || return false
            compute_device(value) == device || return false
        end
    end
    return true
end

function optional_detector_state_device_resident(
    detector::Detector,
    device,
    BackendArray,
)
    state = detector.state
    arrays = (
        state.frame,
        state.presampling_buffer,
        state.presampling_scratch,
        state.response_buffer,
        state.bin_buffer,
        state.temporal_buffer,
        state.noise_buffer,
        state.accum_buffer,
        state.latent_buffer,
    )
    all(
        array -> array isa BackendArray &&
            compute_device(array) == device,
        arrays,
    ) || return false
    state.output_buffer === nothing || (
        state.output_buffer isa BackendArray &&
        compute_device(state.output_buffer) == device
    ) || return false
    state.output_buffer_host === nothing || return false
    return optional_detector_readout_products_device_resident(
        detector.state.readout_products,
        device,
        BackendArray,
    )
end

@inline optional_detector_prepared_storage_device_resident(
    ::PreparedGlobalShutterAcquisition,
    ::Any,
    ::Any,
) = true

@inline optional_detector_prepared_storage_device_resident(
    ::PreparedRollingShutterAcquisition,
    ::Any,
    ::Any,
) = true

@inline function optional_detector_prepared_storage_device_resident(
    prepared::PreparedFrameTransferAcquisition,
    device,
    BackendArray,
)
    return prepared.storage_frame isa BackendArray &&
        compute_device(prepared.storage_frame) == device
end

function optional_detector_device_model_matrix_allocation_bytes(
    row::DeviceModelMatrixDetectorRow,
    result,
)
    warm_output = device_model_matrix_repeat_detector!(
        row,
        result,
        PlantTimestamp(2_000_000_000),
    )
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(warm_output),
    )
    return @allocated begin
        output = device_model_matrix_repeat_detector!(
            row,
            result,
            PlantTimestamp(4_000_000_000),
        )
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(output),
        )
    end
end

function run_optional_detector_device_model_matrix_checks(
    ::Type{B},
    BackendArray,
) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    for row in device_model_matrix_detector_rows()
        oracle = device_model_matrix_execute_detector(
            row;
            backend=CPUBackend(),
            T,
        )
        device_result = device_model_matrix_execute_detector(
            row;
            backend=selector,
            T,
        )
        detector = device_result.detector
        output = device_result.output
        device = compute_device(output)

        @test output isa BackendArray
        @test device isa AdaptiveOpticsSim.Backends.AcceleratorComputeDevice
        @test device_result.map.values isa BackendArray
        @test compute_device(device_result.map.values) == device
        @test device_result.prepared.plan.input_values ===
            device_result.map.values
        @test device_result.prepared.plan.detector_frame ===
            detector.state.frame
        @test optional_detector_response_device_resident(
            detector.params.response_model,
            device,
            BackendArray,
        )
        @test optional_detector_state_device_resident(
            detector,
            device,
            BackendArray,
        )
        @test optional_detector_prepared_storage_device_resident(
            device_result.prepared,
            device,
            BackendArray,
        )

        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(output),
        )
        expected = device_model_matrix_expected_detector_frame(
            row,
            detector,
            T,
        )
        @test isapprox(
            Array(output),
            expected;
            rtol=8f-6,
            atol=8f-7,
        )
        @test isapprox(
            Array(output),
            Array(oracle.output);
            rtol=8f-6,
            atol=8f-7,
        )
        @test device_result.status_trace == oracle.status_trace
        @test device_result.event_times == oracle.event_times
        @test detector_acquisition_sequence(device_result.state) ==
            detector_acquisition_sequence(oracle.state)
        @test device_model_matrix_response_metadata_signature(
            device_result.metadata_before,
        ) == device_model_matrix_response_metadata_signature(
            device_result.metadata_after,
        )
        @test device_model_matrix_response_metadata_signature(
            device_result.metadata_after,
        ) == device_model_matrix_response_metadata_signature(
            oracle.metadata_after,
        )
        @test device_result.metadata_after.sensor ==
            device_model_matrix_detector_sensor_symbol(row)
        @test device_result.metadata_after.timing_model ==
            device_model_matrix_detector_timing_symbol(row)
        @test device_result.metadata_after.sampling_mode ==
            device_model_matrix_detector_sampling_symbol(row)
        @test device_result.metadata_after.acquisition_mode ==
            device_model_matrix_detector_acquisition_symbol(row)

        frequency_x = T(0.25)
        frequency_y = T(0.125)
        @test isapprox(
            detector_mtf(detector, frequency_x, frequency_y),
            device_model_matrix_expected_detector_mtf(
                detector,
                frequency_x,
                frequency_y,
                T,
            );
            rtol=8eps(T),
            atol=8eps(T),
        )

        if row isa DeviceModelMatrixM2FrameTransferEMCCD
            @test frame_transfer_product_sequence(device_result.state) ==
                frame_transfer_product_sequence(oracle.state)
            @test acquisition_product_ready_timestamp(
                device_result.state,
            ) == acquisition_product_ready_timestamp(oracle.state)
        elseif row isa DeviceModelMatrixM5RollingCMOS
            @test rolling_opened_band_count(device_result.state) == 3
            @test rolling_closed_band_count(device_result.state) == 3
        elseif row isa DeviceModelMatrixM6UpTheRampHgCdTe
            @test detector_ramp_cube(detector) isa BackendArray
            @test detector_ramp_slope(detector) isa BackendArray
            @test detector_ramp_intercept(detector) isa BackendArray
            @test isapprox(
                Array(detector_ramp_cube(detector)),
                Array(detector_ramp_cube(oracle.detector));
                rtol=8f-6,
                atol=8f-7,
            )
            @test isapprox(
                Array(detector_ramp_slope(detector)),
                expected;
                rtol=8f-6,
                atol=8f-7,
            )
            @test isapprox(
                Array(detector_ramp_intercept(detector)),
                -expected ./ T(12);
                rtol=8f-6,
                atol=8f-7,
            )
            @test detector_ramp_times(detector) ==
                detector_ramp_times(oracle.detector)
        end

        retained_arrays = (
            detector.state.frame,
            detector.state.presampling_buffer,
            detector.state.presampling_scratch,
            detector.state.accum_buffer,
        )
        allocation_bytes =
            optional_detector_device_model_matrix_allocation_bytes(
                row,
                device_result,
            )
        @test retained_arrays === (
            detector.state.frame,
            detector.state.presampling_buffer,
            detector.state.presampling_scratch,
            detector.state.accum_buffer,
        )
        @test allocation_bytes <= 1024 * 1024
        @test optional_detector_state_device_resident(
            detector,
            device,
            BackendArray,
        )
    end
    return nothing
end

function run_optional_backend_selector_smoke(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(selector)
    T = Float32
    tel = Telescope(resolution=8, diameter=T(1), central_obstruction=T(0),
        T=T, backend=selector)
    dm = DeformableMirror(tel; n_act=2, influence_width=T(0.3), T=T, backend=selector)
    dm_dense = DeformableMirror(tel; n_act=2, influence_model=DenseInfluenceMatrix(Array(dm.state.modes)), T=T, backend=selector)
    native_model = Plant.DeformableMirrorModel(
        n_act=2,
        influence_width=T(0.3),
        T=T,
    )
    native_schema = PlantCommandSchema(
        T,
        (4,);
        id=:optional_native_dm_schema,
        version=1,
        endpoint=:optional_native_dm_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :optional_native_dm),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
    native_definition = ControllableOpticDefinition(
        :optional_native_dm,
        native_model,
        (native_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    native_implementation = Plant.prepare_controllable_optic(
        native_model,
        native_definition,
        tel,
        OptionalStaticAtmosphere(tel; T=T, backend=selector),
    )
    native_initial = BackendArray(zeros(T, 4))
    native_state = Plant.prepare_controllable_optic_state(
        native_implementation,
        native_definition,
        (CommandEndpointID(:optional_native_dm_command),),
        (native_initial,),
    )
    native_workspace = Plant.prepare_controllable_optic_workspace(
        native_implementation,
    )
    native_command = BackendArray(T[1, 0, 0, 0])
    Plant.stage_controllable_optic_command!(
        native_implementation,
        native_state,
        native_workspace,
        CommandEndpointID(:optional_native_dm_command),
        native_command,
        PlantTimestamp(1),
    )
    Plant.commit_controllable_optic_command!(
        native_implementation,
        native_state,
        native_workspace,
        CommandEndpointID(:optional_native_dm_command),
        PlantTimestamp(1),
    )
    set_command!(dm, native_command)
    update_surface!(dm)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(surface_opd(native_state.active)),
    )
    zernike_modal = ModalControllableOptic(tel, ZernikeOpticBasis([2, 3]); T=T, backend=selector)
    cartesian_modal = ModalControllableOptic(tel, CartesianTiltBasis(; scale=T(0.1)); T=T, backend=selector)
    wfs = ShackHartmannWFS(tel; n_lenslets=2, mode=Diffractive(), T=T, backend=selector)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=selector)
    sampled_response = SampledFrameResponse(
        T[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0]; T=T,
        backend=selector)
    sampled_qe = AdaptiveOpticsSim.Detectors.SampledQuantumEfficiency(
        T[0.5e-6, 0.6e-6, 0.7e-6], T[0.2, 0.8, 0.4]; T=T)
    calibration_det = Detector(noise=NoiseNone(), integration_time=T(1),
        qe=sampled_qe, response_model=sampled_response, T=T,
        backend=selector)
    pupil = PupilFunction(tel; T=T, backend=selector)
    device = compute_device(pupil.opd)
    @test device isa AdaptiveOpticsSim.Backends.AcceleratorComputeDevice
    @test typeof(AdaptiveOpticsSim.Backends.compute_device_backend(device)) ===
        typeof(selector)
    @test !isnothing(
        AdaptiveOpticsSim.Backends.compute_device_identifier(device))
    @test compute_device(pupil.amplitude) == device
    compute_device(pupil.opd)
    @test @allocated(compute_device(pupil.opd)) == 0
    @test pupil.opd isa BackendArray
    @test dm.state.coefs isa BackendArray
    @test dm.state.modes isa
        AdaptiveOpticsSim.Optics.GaussianInfluenceOperator
    @test typeof(backend(dm.state.modes)) === typeof(selector)
    @test similar(dm.state.modes, T, 2, 2) isa BackendArray
    @test AdaptiveOpticsSim.Optics.materialize_influence_matrix(dm) isa
        BackendArray
    @test dm_dense.state.coefs isa BackendArray
    @test dm_dense.state.modes isa BackendArray
    @test Array(dm_dense.state.modes) ≈ Array(dm.state.modes) atol=0 rtol=0
    @test surface_opd(native_state.active) isa BackendArray
    @test native_state.active.state.coefs isa BackendArray
    @test native_workspace.staged.state.coefs isa BackendArray
    @test Array(surface_opd(native_state.active)) ≈
        Array(surface_opd(dm)) rtol=T(2e-5) atol=T(2e-6)
    @test zernike_modal.state.coefs isa BackendArray
    @test zernike_modal.state.modes isa BackendArray
    @test cartesian_modal.state.coefs isa BackendArray
    @test cartesian_modal.state.modes isa BackendArray
    @test slopes(wfs) isa BackendArray
    @test det.state.frame isa BackendArray
    @test calibration_det.params.response_model.kernel isa BackendArray
    seed = UInt(0x51a7)
    AdaptiveOpticsSim.detector_calibration_signature(calibration_det, seed)
    optional_detector_calibration_signature_allocation_bytes(
        calibration_det, seed)
    @test optional_detector_calibration_signature_allocation_bytes(
        calibration_det, seed) == 0

    gain_map = reshape(T.(range(T(0.6), T(1.2); length=16)), 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    defect_det = Detector(noise=NoiseNone(), sensor=CMOSSensor(T=T,
            backend=selector),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(gain_map; T=T, backend=selector),
            BadPixelMask(bad_mask; T=T, backend=selector)),
        T=T, backend=selector)
    photon_rate = array_backend(fill(T(2), 4, 4))
    calibration_frame = AdaptiveOpticsSim.detector_calibration_frame!(
        defect_det, photon_rate, one(T))
    expected_frame = T(2) .* gain_map
    expected_frame[2, 3] = zero(T)
    @test Array(calibration_frame) ≈ expected_frame rtol=T(1e-6)
    @test optional_detector_calibration_signature_allocation_bytes(
        defect_det, seed) == 0

    invalid_values = array_backend(T[-1 1; 1 1])
    invalid_metadata = OpticalPlaneMetadata(FocalPlane(), invalid_values;
        coordinate_domain=AngularCoordinates(), sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    invalid_map = IntensityMap(invalid_metadata, invalid_values)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        det, invalid_map)
    return nothing
end

function run_optional_lgs_convolution_normalization(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(selector)
    T = Float32
    n = 8
    expected = reshape(collect(range(T(0.25), T(2); length=n * n)), n, n)
    expected_stack = cat(expected, reverse(expected; dims=2); dims=3)
    intensity_stack = array_backend(copy(expected_stack))
    kernel_fft = array_backend(ones(Complex{T}, n, n, 2))
    fft_stack = array_backend(zeros(Complex{T}, n, n, 2))
    fft_plan = AdaptiveOpticsSim.Backends.plan_fft_backend!(fft_stack, (1, 2))
    ifft_plan = AdaptiveOpticsSim.Backends.plan_ifft_backend!(fft_stack, (1, 2))

    WavefrontSensors.apply_lgs_convolution_stack!(
        intensity_stack, kernel_fft, fft_stack, fft_plan, ifft_plan)
    actual = Array(intensity_stack)
    @test actual ≈ expected_stack rtol=T(1e-5) atol=T(1e-6)
    @test vec(sum(actual; dims=(1, 2))) ≈
          vec(sum(expected_stack; dims=(1, 2))) rtol=T(1e-5)
    return nothing
end

function run_optional_sodium_profile_wfs(::Type{B},
    BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)
    pupil = PupilFunction(tel; T=T, backend=selector)

    for family in (:pyramid, :bioedge)
        src = LGSSource(
            na_profile=T[80000 90000 100000; 0.2 0.6 0.2],
            laser_coordinates=(T(1), T(-0.5)),
            fwhm_spot_up=T(0.8),
            photon_irradiance=one(T),
            T=T,
        )
        wfs = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=zero(T), T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=zero(T), T=T, backend=selector)

        WavefrontSensors.ensure_lgs_kernel!(wfs, pupil, src)
        propagation = wfs.front_end.propagation
        @test propagation.lgs_kernel_fft isa BackendArray
        original_tag = propagation.lgs_kernel_tag
        original_kernel = Array(propagation.lgs_kernel_fft)
        @test all(isfinite, original_kernel)
        slopes = measure!(wfs, pupil, src)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(slopes))
        @test slopes isa BackendArray
        @test all(isfinite, Array(slopes))

        src.params.na_profile[2, :] .= T[0.8, 0.1, 0.1]
        WavefrontSensors.ensure_lgs_kernel!(wfs, pupil, src)
        @test propagation.lgs_kernel_tag != original_tag
        @test !isapprox(Array(propagation.lgs_kernel_fft), original_kernel;
            rtol=T(1e-5), atol=T(1e-6))
    end
    return nothing
end

function run_optional_zernike_normalization(
    ::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    cpu_tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)
    cpu_pupil = PupilFunction(cpu_tel; T=T, backend=CPUBackend())
    gpu_pupil = PupilFunction(gpu_tel; T=T, backend=selector)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    normalization_scale = T(0.375)
    frame_host = reshape(collect(range(T(0.25), T(2); length=16)),
        4, 4)

    for normalization in (MeanValidFluxNormalization(),
            IncidenceFluxNormalization())
        cpu_wfs = ZernikeWFS(cpu_tel; pupil_samples=8, binning=2,
            normalization=normalization, T=T, backend=CPUBackend())
        gpu_wfs = ZernikeWFS(gpu_tel; pupil_samples=8, binning=2,
            normalization=normalization, T=T, backend=selector)
        fill!(cpu_wfs.estimator.state.reference_signal_2d, zero(T))
        fill!(gpu_wfs.estimator.state.reference_signal_2d, zero(T))
        frame = BackendArray(copy(frame_host))

        expected_normalization = AdaptiveOpticsSim.WavefrontSensors.zernike_normalization(
            normalization, cpu_wfs, cpu_pupil, src, frame_host,
            normalization_scale)
        actual_normalization = AdaptiveOpticsSim.WavefrontSensors.zernike_normalization(
            normalization, gpu_wfs, gpu_pupil, src, frame,
            normalization_scale)
        @test actual_normalization ≈ expected_normalization rtol=T(2e-5)

        expected = copy(AdaptiveOpticsSim.WavefrontSensors.zernike_signal!(cpu_wfs,
            cpu_pupil, frame_host, src, normalization_scale))
        actual = AdaptiveOpticsSim.WavefrontSensors.zernike_signal!(gpu_wfs, gpu_pupil,
            frame, src, normalization_scale)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(actual))
        @test gpu_wfs.estimator.state.normalization_sum isa BackendArray
        @test actual isa BackendArray
        @test Array(actual) ≈ expected rtol=T(2e-5) atol=T(2e-6)
    end

    zero_src = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=zero(T), T=T)
    zero_wfs = ZernikeWFS(gpu_tel; pupil_samples=8, binning=2,
        normalization=IncidenceFluxNormalization(), T=T,
        backend=selector)
    fill!(zero_wfs.estimator.state.reference_signal_2d, zero(T))
    zero_slopes = AdaptiveOpticsSim.WavefrontSensors.zernike_signal!(zero_wfs, gpu_pupil,
        BackendArray(copy(frame_host)), zero_src, one(T))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(zero_slopes))
    @test all(iszero, Array(zero_slopes))
    @test all(isfinite, Array(zero_slopes))
    return nothing
end

function run_optional_wfs_stage_contracts(
    ::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T, backend=selector)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(4), T=T)
    pupil = PupilFunction(tel; T=T, backend=selector)
    fill!(pupil.opd, T(2e-9))
    field = ElectricField(pupil, src; zero_padding=1, T=T)
    field_formation = prepare_pupil_field(pupil, src, field;
        center_even_grid=false)
    fill_electric_field!(field, pupil, field_formation)

    rate_values = BackendArray(zeros(T, 4, 4))
    rate_metadata = OpticalPlaneMetadata(DetectorPlane(), rate_values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    rate = IntensityMap(rate_metadata, rate_values)
    optical_model = ContractRateModel(T(3), T(1e6), rate.metadata)
    optical_plan = prepare_wfs_optical_formation(optical_model, pupil, rate)

    detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    binding = contract_detector_binding(detector, rate)
    observation = binding.observation
    acquisition_plan = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((binding,)), rate, observation)
    measurement_values = similar(observation.storage)
    measurement = WFSMeasurement(measurement_values;
        units=:detector_signal, kind=:copied_detector_frame)
    estimator_plan = prepare_wfs_estimation(ContractCopyEstimator(
        :detector_signal, :copied_detector_frame),
        observation, measurement)
    rng = Xoshiro(0x574653)

    @test @inferred(form_wfs_optical_products!(rate, pupil,
        optical_plan)) === rate
    @test @inferred(acquire_wfs_observation!(observation, rate,
        acquisition_plan, rng)) === observation
    @test @inferred(estimate_wfs_measurement!(measurement, observation,
        estimator_plan)) === measurement
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(measurement.storage))
    @test rate.values isa BackendArray
    @test observation.storage isa BackendArray
    @test measurement.storage isa BackendArray
    @test observation.metadata.device == compute_device(observation.storage)
    @test measurement.metadata.device == compute_device(measurement.storage)
    @test Array(measurement.storage) == Array(observation.storage)
    @test sum(Array(observation.storage)) ≈
        sum(Array(rate.values)) * T(0.4) * T(0.5) rtol=T(2e-6)

    direct_values = similar(rate.values)
    direct_measurement = WFSMeasurement(direct_values;
        units=:field_intensity, kind=:direct_field_measurement)
    direct_plan = prepare_wfs_estimation(ContractDirectCopyEstimator(
        :field_intensity, :direct_field_measurement), field,
        direct_measurement)
    @test wfs_measurement_path(direct_plan) isa DirectMeasurementPath
    @test @inferred(estimate_wfs_measurement!(direct_measurement, field,
        direct_plan)) === direct_measurement
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(direct_measurement.storage))
    @test direct_measurement.storage isa BackendArray
    @test Array(direct_measurement.storage) ≈ abs2.(Array(field.values))

    second_values = BackendArray(zeros(T, 4, 4))
    second_metadata = OpticalPlaneMetadata(DetectorPlane(), second_values;
        coordinate_domain=AngularCoordinates(),
        sampling=(T(0.5), T(0.5)),
        spectral=MonochromaticChannel(T(0.9e-6)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    second_rate = IntensityMap(second_metadata, second_values)
    bundle = OpticalProductBundle(rate, second_rate)
    bundle_model = ContractBundleRateModel((
        ContractRateModel(T(2), T(1e6), rate.metadata),
        ContractRateModel(T(5), T(1e6), second_rate.metadata),
    ))
    bundle_plan = prepare_wfs_optical_formation(bundle_model, pupil, bundle)
    form_wfs_optical_products!(bundle, pupil, bundle_plan)
    second_detector = Detector(noise=NoiseNone(), integration_time=T(0.7),
        qe=T(0.25), response_model=NullFrameResponse(), T=T,
        backend=selector)
    first_binding = contract_detector_binding(detector, rate)
    second_binding = contract_detector_binding(second_detector, second_rate)
    observations = (first_binding.observation, second_binding.observation)
    multi_acquisition = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((first_binding, second_binding)),
        bundle, observations)
    @test @inferred(acquire_wfs_observation!(observations, bundle,
        multi_acquisition, (Xoshiro(1), Xoshiro(2)))) === observations
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(second_binding.observation.storage))
    @test first_binding.observation.storage isa BackendArray
    @test second_binding.observation.storage isa BackendArray
    @test first_binding.observation.metadata.device ==
        second_binding.observation.metadata.device

    packed_values = BackendArray(zeros(T, 8, 4))
    regions = ((1:4, 1:4), (5:8, 1:4))
    packed = WFSObservation(packed_values; units=:detector_signal,
        layout=regions)
    packed_plan = prepare_wfs_acquisition(
        ContractPackedAcquisition(regions, T(0.2)), bundle, packed)
    @test @inferred(acquire_wfs_observation!(packed, bundle, packed_plan,
        Xoshiro(3))) === packed
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(packed.storage))
    @test packed.storage isa BackendArray
    packed_host = Array(packed.storage)
    @test packed_host[1:4, :] ≈ Array(rate.values) .* T(0.2)
    @test packed_host[5:8, :] ≈ Array(second_rate.values) .* T(0.2)

    physical_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_rate = shack_hartmann_rate_map(physical_wfs, pupil, src)
    physical_front_end = ShackHartmannOpticalFrontEnd(
        physical_wfs.front_end, src)
    @test !hasfield(typeof(physical_front_end), :sensor)
    @test physical_front_end.propagation.fft_stack isa BackendArray
    physical_optical_plan = prepare_wfs_optical_formation(
        physical_front_end, pupil, physical_rate)
    form_wfs_optical_products!(physical_rate, pupil, physical_optical_plan)
    @test physical_rate.values isa BackendArray
    @test all(isfinite, Array(physical_rate.values))
    @test sum(Array(physical_rate.values)) > zero(T)

    physical_field_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_field_rate = shack_hartmann_rate_map(physical_field_wfs, field)
    physical_field_plan = prepare_wfs_optical_formation(
        physical_field_wfs.front_end, field,
        physical_field_rate)
    form_wfs_optical_products!(physical_field_rate, field,
        physical_field_plan)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(physical_field_rate.values))
    @test physical_field_rate.values isa BackendArray
    @test isapprox(Array(physical_field_rate.values),
        Array(physical_rate.values); rtol=T(2e-5), atol=T(2e-5))

    physical_asterism = Asterism([
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(1), T=T),
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(3), T=T),
    ])
    physical_asterism_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_asterism_rate = shack_hartmann_rate_map(
        physical_asterism_wfs, pupil, physical_asterism)
    physical_asterism_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(physical_asterism_wfs.front_end,
            physical_asterism), pupil, physical_asterism_rate)
    form_wfs_optical_products!(physical_asterism_rate, pupil,
        physical_asterism_plan)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(physical_asterism_rate.values))
    @test physical_asterism_rate.values isa BackendArray
    @test isapprox(Array(physical_asterism_rate.values),
        Array(physical_rate.values); rtol=T(2e-5), atol=T(2e-5))

    physical_detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    physical_observation = WFSObservation(similar(physical_rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    physical_acquisition = prepare_wfs_acquisition(physical_detector,
        physical_rate, physical_observation)
    acquire_wfs_observation!(physical_observation, physical_rate,
        physical_acquisition, Xoshiro(0x5348))
    @test physical_observation.storage isa BackendArray
    @test sum(Array(physical_observation.storage)) ≈
        sum(Array(physical_rate.values)) * T(0.4) * T(0.5) rtol=T(2e-5)

    cpu_tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T)
    cpu_wfs = ShackHartmannWFS(cpu_tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T)
    mixed_backend_rate = shack_hartmann_rate_map(cpu_wfs, pupil, src)
    mixed_optical_error = try
        prepare_wfs_optical_formation(
            ShackHartmannOpticalFrontEnd(cpu_wfs.front_end, src), pupil,
            mixed_backend_rate)
        nothing
    catch err
        err
    end
    @test mixed_optical_error isa WFSPreparationError
    @test mixed_optical_error.reason === :backend

    cpu_observation = WFSObservation(zeros(T, size(physical_rate.values));
        units=:electron_count, layout=:lenslet_mosaic)
    gpu_measurement = WFSMeasurement(similar(slopes(physical_wfs));
        units=:pixel, kind=:centroid_slopes)
    @test_throws WFSPreparationError prepare_wfs_estimation(physical_wfs,
        cpu_observation, gpu_measurement)
    cpu_measurement = WFSMeasurement(zeros(T, length(slopes(physical_wfs)));
        units=:pixel, kind=:centroid_slopes)
    @test_throws WFSPreparationError prepare_wfs_estimation(physical_wfs,
        physical_observation, cpu_measurement)

    set_subaperture_calibration!(physical_wfs.calibration,
        zeros(T, size(physical_wfs.calibration.reference_signal_2d));
        centroid_response=one(T), wavelength=wavelength(src),
        signature=UInt(0x50485953))
    physical_measurement = WFSMeasurement(similar(slopes(physical_wfs));
        units=:pixel, kind=:centroid_slopes)
    physical_estimator = prepare_wfs_estimation(physical_wfs,
        physical_observation, physical_measurement)
    estimate_wfs_measurement!(physical_measurement, physical_observation,
        physical_estimator)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(physical_measurement.storage))
    @test physical_measurement.storage isa BackendArray
    @test all(isfinite, Array(physical_measurement.storage))

    four_pupil_cpu_tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    four_pupil_cpu = PupilFunction(four_pupil_cpu_tel; T=T)
    copyto!(four_pupil_cpu.opd, Array(pupil.opd))
    for family in (:pyramid, :bioedge)
        cpu_sensor = family === :pyramid ?
            PyramidWFS(four_pupil_cpu_tel; pupil_samples=2,
                mode=Diffractive(), modulation=0, T=T) :
            BioEdgeWFS(four_pupil_cpu_tel; pupil_samples=2,
                mode=Diffractive(), modulation=0, T=T)
        gpu_sensor = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector)
        cpu_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(cpu_sensor, src) :
            BioEdgeOpticalFrontEnd(cpu_sensor, src)
        gpu_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, src) :
            BioEdgeOpticalFrontEnd(gpu_sensor, src)
        cpu_rate = family === :pyramid ?
            pyramid_rate_map(cpu_front_end, four_pupil_cpu) :
            bioedge_rate_map(cpu_front_end, four_pupil_cpu)
        gpu_rate = family === :pyramid ?
            pyramid_rate_map(gpu_front_end, pupil) :
            bioedge_rate_map(gpu_front_end, pupil)
        cpu_plan = prepare_wfs_optical_formation(cpu_front_end,
            four_pupil_cpu, cpu_rate)
        gpu_plan = prepare_wfs_optical_formation(gpu_front_end, pupil,
            gpu_rate)
        form_wfs_optical_products!(cpu_rate, four_pupil_cpu, cpu_plan)
        form_wfs_optical_products!(gpu_rate, pupil, gpu_plan)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(gpu_rate.values))
        @test gpu_rate.values isa BackendArray
        @test isapprox(Array(gpu_rate.values), cpu_rate.values;
            rtol=T(3e-5), atol=T(3e-5))

        field_sensor = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector)
        field_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(field_sensor) :
            BioEdgeOpticalFrontEnd(field_sensor)
        field_rate = family === :pyramid ?
            pyramid_rate_map(field_front_end, field) :
            bioedge_rate_map(field_front_end, field)
        field_plan = prepare_wfs_optical_formation(field_front_end, field,
            field_rate)
        form_wfs_optical_products!(field_rate, field, field_plan)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(field_rate.values))
        @test field_rate.values isa BackendArray
        @test isapprox(Array(field_rate.values), Array(gpu_rate.values);
            rtol=T(3e-5), atol=T(3e-5))

        four_pupil_detector = Detector(noise=NoiseNone(),
            integration_time=T(0.4), qe=T(0.5),
            response_model=NullFrameResponse(), T=T, backend=selector)
        four_pupil_observation = WFSObservation(similar(gpu_rate.values);
            units=:electron_count, layout=:four_pupil_mosaic)
        four_pupil_acquisition = prepare_wfs_acquisition(
            four_pupil_detector, gpu_rate, four_pupil_observation)
        acquire_wfs_observation!(four_pupil_observation, gpu_rate,
            four_pupil_acquisition, Xoshiro(0x5042))
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(
                four_pupil_observation.storage))
        @test four_pupil_observation.storage isa BackendArray
        @test isapprox(Array(four_pupil_observation.storage),
            Array(gpu_rate.values) .* T(0.2); rtol=T(3e-5), atol=T(3e-5))

        reference = similar(gpu_sensor.estimator.state.reference_signal_2d)
        fill!(reference, zero(T))
        if family === :pyramid
            set_pyramid_calibration!(gpu_sensor, reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        else
            set_bioedge_calibration!(gpu_sensor, reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        end

        host_observation = WFSObservation(zeros(T, size(gpu_rate.values));
            units=:electron_count, layout=:four_pupil_mosaic)
        host_observation_error = try
            prepare_wfs_estimation(gpu_sensor, host_observation,
                WFSMeasurement(similar(slopes(gpu_sensor));
                    units=:dimensionless, kind=:differential_slopes))
            nothing
        catch err
            err
        end
        @test host_observation_error isa WFSPreparationError
        @test host_observation_error.reason === :backend

        host_measurement = WFSMeasurement(
            zeros(T, length(slopes(gpu_sensor))); units=:dimensionless,
            kind=:differential_slopes)
        host_measurement_error = try
            prepare_wfs_estimation(gpu_sensor, four_pupil_observation,
                host_measurement)
            nothing
        catch err
            err
        end
        @test host_measurement_error isa WFSPreparationError
        @test host_measurement_error.reason === :backend

        four_pupil_measurement = WFSMeasurement(similar(slopes(gpu_sensor));
            units=:dimensionless, kind=:differential_slopes)
        four_pupil_estimator = prepare_wfs_estimation(gpu_sensor,
            four_pupil_observation, four_pupil_measurement)
        estimate_wfs_measurement!(four_pupil_measurement,
            four_pupil_observation, four_pupil_estimator)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(
                four_pupil_measurement.storage))
        @test four_pupil_measurement.storage isa BackendArray
        @test all(isfinite, Array(four_pupil_measurement.storage))

        quantized_host = reshape(UInt16.(1:length(gpu_rate.values)),
            size(gpu_rate.values))
        quantized_storage = similar(gpu_rate.values, UInt16)
        copyto!(quantized_storage, quantized_host)
        quantized_observation = WFSObservation(quantized_storage;
            units=:adu, layout=:four_pupil_mosaic)
        quantized_measurement = WFSMeasurement(similar(slopes(gpu_sensor));
            units=:dimensionless, kind=:differential_slopes)
        quantized_estimator = prepare_wfs_estimation(gpu_sensor,
            quantized_observation, quantized_measurement)
        estimate_wfs_measurement!(quantized_measurement,
            quantized_observation, quantized_estimator)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(quantized_measurement.storage))

        cpu_reference = zeros(T,
            size(cpu_sensor.estimator.state.reference_signal_2d))
        if family === :pyramid
            set_pyramid_calibration!(cpu_sensor, cpu_reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        else
            set_bioedge_calibration!(cpu_sensor, cpu_reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        end
        cpu_quantized_observation = WFSObservation(quantized_host;
            units=:adu, layout=:four_pupil_mosaic)
        cpu_quantized_measurement = WFSMeasurement(similar(slopes(cpu_sensor));
            units=:dimensionless, kind=:differential_slopes)
        cpu_quantized_estimator = prepare_wfs_estimation(cpu_sensor,
            cpu_quantized_observation, cpu_quantized_measurement)
        estimate_wfs_measurement!(cpu_quantized_measurement,
            cpu_quantized_observation, cpu_quantized_estimator)
        @test quantized_measurement.storage isa BackendArray
        @test isapprox(Array(quantized_measurement.storage),
            cpu_quantized_measurement.storage; rtol=T(3e-5), atol=T(3e-5))

        geometric = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Geometric(), T=T,
                backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Geometric(), T=T,
                backend=selector)
        geometric_measurement = WFSMeasurement(similar(slopes(geometric));
            units=:metre, kind=:geometric_slopes)
        geometric_plan = prepare_wfs_estimation(geometric, pupil,
            geometric_measurement)
        estimate_wfs_measurement!(geometric_measurement, pupil,
            geometric_plan)
        AdaptiveOpticsSim.Backends.synchronize_backend!(
            AdaptiveOpticsSim.Backends.execution_style(
                geometric_measurement.storage))
        @test geometric_measurement.storage isa BackendArray
        @test all(isfinite, Array(geometric_measurement.storage))

        spectral_source = with_spectrum(src,
            SpectralBundle(T[0.7e-6, 0.9e-6], T[0.25, 0.75]; T=T))
        spectral_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, spectral_source) :
            BioEdgeOpticalFrontEnd(gpu_sensor, spectral_source)
        spectral_rates = family === :pyramid ?
            pyramid_rate_map(spectral_front_end, pupil) :
            bioedge_rate_map(spectral_front_end, pupil)
        spectral_plan = prepare_wfs_optical_formation(spectral_front_end,
            pupil, spectral_rates)
        form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
        @test all(product -> product.values isa BackendArray,
            spectral_rates)

        path_source = Asterism([
            Source(band=:custom, wavelength=wavelength(src),
                photon_irradiance=T(1), T=T),
            Source(band=:custom, wavelength=wavelength(src),
                photon_irradiance=T(3), T=T),
        ])
        second_pupil = PupilFunction(tel; T=T, backend=selector)
        copyto!(second_pupil.opd, pupil.opd)
        path_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, path_source) :
            BioEdgeOpticalFrontEnd(gpu_sensor, path_source)
        path_inputs = (pupil, second_pupil)
        path_rates = family === :pyramid ?
            pyramid_rate_map(path_front_end, path_inputs) :
            bioedge_rate_map(path_front_end, path_inputs)
        path_plan = prepare_wfs_optical_formation(path_front_end,
            path_inputs, path_rates)
        form_wfs_optical_products!(path_rates, path_inputs, path_plan)
        @test all(product -> product.values isa BackendArray, path_rates)
        @test isapprox(sum(Array(path_rates[1].values)) /
            sum(Array(path_rates[2].values)), T(1 / 3); rtol=T(3e-5))

        for lgs in (
            LGSSource(wavelength=wavelength(src),
                photon_irradiance=T(4), elongation_factor=T(1.8), T=T),
            LGSSource(wavelength=wavelength(src),
                photon_irradiance=T(4),
                na_profile=T[80_000 90_000 100_000; 0.2 0.6 0.2],
                laser_coordinates=(T(1), T(-0.5)),
                fwhm_spot_up=T(0.8), T=T),
        )
            cpu_lgs_sensor = family === :pyramid ?
                PyramidWFS(four_pupil_cpu_tel; pupil_samples=2,
                    mode=Diffractive(), modulation=0, T=T) :
                BioEdgeWFS(four_pupil_cpu_tel; pupil_samples=2,
                    mode=Diffractive(), modulation=0, T=T)
            gpu_lgs_sensor = family === :pyramid ?
                PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                    modulation=0, T=T, backend=selector) :
                BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                    modulation=0, T=T, backend=selector)
            cpu_lgs_front_end = family === :pyramid ?
                PyramidOpticalFrontEnd(cpu_lgs_sensor, lgs) :
                BioEdgeOpticalFrontEnd(cpu_lgs_sensor, lgs)
            gpu_lgs_front_end = family === :pyramid ?
                PyramidOpticalFrontEnd(gpu_lgs_sensor, lgs) :
                BioEdgeOpticalFrontEnd(gpu_lgs_sensor, lgs)
            cpu_lgs_rate = family === :pyramid ?
                pyramid_rate_map(cpu_lgs_front_end, four_pupil_cpu) :
                bioedge_rate_map(cpu_lgs_front_end, four_pupil_cpu)
            gpu_lgs_rate = family === :pyramid ?
                pyramid_rate_map(gpu_lgs_front_end, pupil) :
                bioedge_rate_map(gpu_lgs_front_end, pupil)
            cpu_lgs_plan = prepare_wfs_optical_formation(
                cpu_lgs_front_end, four_pupil_cpu, cpu_lgs_rate)
            gpu_lgs_plan = prepare_wfs_optical_formation(
                gpu_lgs_front_end, pupil, gpu_lgs_rate)
            form_wfs_optical_products!(cpu_lgs_rate, four_pupil_cpu,
                cpu_lgs_plan)
            form_wfs_optical_products!(gpu_lgs_rate, pupil, gpu_lgs_plan)
            AdaptiveOpticsSim.Backends.synchronize_backend!(
                AdaptiveOpticsSim.Backends.execution_style(gpu_lgs_rate.values))
            @test isapprox(Array(gpu_lgs_rate.values), cpu_lgs_rate.values;
                rtol=T(5e-5), atol=T(5e-5))
        end
    end

    cpu_measurement = WFSMeasurement(zeros(T, 4, 4);
        units=:detector_signal, kind=:cpu_copy)
    mismatch = try
        prepare_wfs_estimation(ContractCopyEstimator(
            :detector_signal, :cpu_copy), observation,
            cpu_measurement)
        nothing
    catch err
        err
    end
    @test mismatch isa WFSPreparationError
    @test mismatch.reason === :backend
    run_optional_zernike_curvature_stages(B, BackendArray)
    return nothing
end

function run_optional_zernike_curvature_stages(
    ::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    cpu_tel = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T, backend=selector)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(8), T=T)
    cpu_pupil = PupilFunction(cpu_tel; T=T)
    gpu_pupil = PupilFunction(gpu_tel; T=T, backend=selector)
    opd = reshape(T.(1:64), 8, 8) .* T(1e-10)
    copyto!(cpu_pupil.opd, opd)
    copyto!(gpu_pupil.opd, opd)

    cpu_zernike = ZernikeWFS(cpu_tel; pupil_samples=4, T=T)
    gpu_zernike = ZernikeWFS(gpu_tel; pupil_samples=4, T=T,
        backend=selector)
    cpu_zernike_rate = zernike_rate_map(
        ZernikeOpticalFrontEnd(cpu_zernike, source), cpu_pupil)
    gpu_zernike_rate = zernike_rate_map(
        ZernikeOpticalFrontEnd(gpu_zernike, source), gpu_pupil)
    cpu_zernike_plan = prepare_wfs_optical_formation(
        ZernikeOpticalFrontEnd(cpu_zernike, source), cpu_pupil,
        cpu_zernike_rate)
    gpu_zernike_plan = prepare_wfs_optical_formation(
        ZernikeOpticalFrontEnd(gpu_zernike, source), gpu_pupil,
        gpu_zernike_rate)
    form_wfs_optical_products!(cpu_zernike_rate, cpu_pupil,
        cpu_zernike_plan)
    form_wfs_optical_products!(gpu_zernike_rate, gpu_pupil,
        gpu_zernike_plan)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(gpu_zernike_rate.values))
    @test gpu_zernike_rate.values isa BackendArray
    @test isapprox(Array(gpu_zernike_rate.values),
        cpu_zernike_rate.values; rtol=T(4e-5), atol=T(4e-5))

    zernike_detector = Detector(noise=NoiseNone(),
        integration_time=T(0.25), qe=T(0.4),
        response_model=NullFrameResponse(), T=T, backend=selector)
    zernike_observation = WFSObservation(similar(gpu_zernike_rate.values);
        units=:electron_count, layout=:zernike_pupil_image)
    zernike_acquisition = prepare_wfs_acquisition(zernike_detector,
        gpu_zernike_rate, zernike_observation)
    acquire_wfs_observation!(zernike_observation, gpu_zernike_rate,
        zernike_acquisition, Xoshiro(0x5a47))
    zernike_reference = similar(gpu_zernike.estimator.state.signal_2d)
    fill!(zernike_reference, zero(T))
    set_zernike_calibration!(gpu_zernike, zernike_reference;
        wavelength_m=wavelength(source), signature=UInt(0x5a47))
    zernike_measurement = WFSMeasurement(similar(slopes(gpu_zernike));
        units=:dimensionless, kind=:normalized_pupil_signal)
    zernike_estimator = prepare_wfs_estimation(gpu_zernike,
        zernike_observation, zernike_measurement; source=source)
    estimate_wfs_measurement!(zernike_measurement, zernike_observation,
        zernike_estimator)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(zernike_measurement.storage))
    @test zernike_observation.storage isa BackendArray
    @test zernike_measurement.storage isa BackendArray
    @test all(isfinite, Array(zernike_measurement.storage))

    cpu_curvature = CurvatureWFS(cpu_tel; pupil_samples=4, T=T)
    gpu_curvature = CurvatureWFS(gpu_tel; pupil_samples=4, T=T,
        backend=selector)
    cpu_rates = curvature_rate_maps(
        CurvatureOpticalFrontEnd(cpu_curvature, source), cpu_pupil)
    gpu_rates = curvature_rate_maps(
        CurvatureOpticalFrontEnd(gpu_curvature, source), gpu_pupil)
    cpu_curvature_plan = prepare_wfs_optical_formation(
        CurvatureOpticalFrontEnd(cpu_curvature, source), cpu_pupil,
        cpu_rates)
    gpu_curvature_plan = prepare_wfs_optical_formation(
        CurvatureOpticalFrontEnd(gpu_curvature, source), gpu_pupil,
        gpu_rates)
    form_wfs_optical_products!(cpu_rates, cpu_pupil, cpu_curvature_plan)
    form_wfs_optical_products!(gpu_rates, gpu_pupil, gpu_curvature_plan)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(gpu_rates[1].values))
    @test all(rate -> rate.values isa BackendArray, gpu_rates)
    @test isapprox(Array(gpu_rates[1].values), cpu_rates[1].values;
        rtol=T(5e-5), atol=T(5e-5))
    @test isapprox(Array(gpu_rates[2].values), cpu_rates[2].values;
        rtol=T(5e-5), atol=T(5e-5))

    plus_detector = Detector(noise=NoiseNone(),
        integration_time=T(0.25), qe=T(0.4),
        response_model=NullFrameResponse(), T=T, backend=selector)
    minus_detector = Detector(noise=NoiseNone(),
        integration_time=T(0.5), qe=T(0.5),
        response_model=NullFrameResponse(), T=T, backend=selector)
    plus_observation = WFSObservation(similar(gpu_rates[1].values);
        units=:electron_count, layout=:curvature_branch_image)
    minus_observation = WFSObservation(similar(gpu_rates[2].values);
        units=:electron_count, layout=:curvature_branch_image)
    observations = (plus_observation, minus_observation)
    multiple_acquisition = prepare_wfs_acquisition(
        (plus_detector, minus_detector), gpu_rates, observations;
        source=source)
    acquire_wfs_observation!(observations, gpu_rates,
        multiple_acquisition, Xoshiro(0x4355))
    @test isapprox(Array(plus_observation.storage),
        Array(gpu_rates[1].values) .* T(0.1);
        rtol=T(3e-5), atol=T(3e-5))
    @test isapprox(Array(minus_observation.storage),
        Array(gpu_rates[2].values) .* T(0.25);
        rtol=T(3e-5), atol=T(3e-5))

    packed_detector = Detector(noise=NoiseNone(),
        integration_time=T(0.25), qe=T(0.5),
        response_model=NullFrameResponse(), T=T, backend=selector)
    packed_model = CurvaturePackedAcquisition(packed_detector)
    packed_storage = similar(gpu_rates[1].values, T, 8, 4)
    packed_observation = WFSObservation(packed_storage;
        units=:electron_count, layout=:curvature_branch_regions)
    packed_acquisition = prepare_wfs_acquisition(packed_model, gpu_rates,
        packed_observation)
    acquire_wfs_observation!(packed_observation, gpu_rates,
        packed_acquisition, Xoshiro(0x4356))
    packed_host = Array(packed_observation.storage)
    @test packed_observation.storage isa BackendArray
    @test isapprox(packed_host[1:4, :], Array(gpu_rates[1].values) .* T(0.125);
        rtol=T(3e-5), atol=T(3e-5))
    @test isapprox(packed_host[5:8, :], Array(gpu_rates[2].values) .* T(0.125);
        rtol=T(3e-5), atol=T(3e-5))

    spad = SPADArrayDetector(integration_time=T(0.25), noise=NoiseNone(),
        sensor=SPADArraySensor(pde=T(0.5), dark_count_rate=zero(T),
            fill_factor=one(T)), T=T, backend=selector)
    counting_model = CurvaturePackedAcquisition(spad;
        readout_model=CurvatureCountingReadout(), source=source)
    counting_storage = similar(gpu_rates[1].values, T, 2, 16)
    counting_observation = WFSObservation(counting_storage;
        units=:photon_count, layout=:curvature_branch_channels)
    counting_acquisition = prepare_wfs_acquisition(counting_model,
        gpu_rates, counting_observation)
    acquire_wfs_observation!(counting_observation, gpu_rates,
        counting_acquisition, Xoshiro(0x4357))
    @test counting_observation.storage isa BackendArray
    @test all(isfinite, Array(counting_observation.storage))

    curvature_reference = similar(gpu_curvature.estimator.state.signal_2d)
    fill!(curvature_reference, zero(T))
    set_curvature_calibration!(gpu_curvature, curvature_reference;
        wavelength_m=wavelength(source), signature=UInt(0x4355))
    curvature_measurement = WFSMeasurement(similar(slopes(gpu_curvature));
        units=:dimensionless, kind=:curvature_signal)
    curvature_estimator = prepare_wfs_estimation(gpu_curvature,
        observations, curvature_measurement;
        branch_rate_scales=(T(10), T(4)))
    estimate_wfs_measurement!(curvature_measurement, observations,
        curvature_estimator)
    packed_measurement = WFSMeasurement(similar(slopes(gpu_curvature));
        units=:dimensionless, kind=:curvature_signal)
    packed_estimator = prepare_wfs_estimation(gpu_curvature,
        packed_observation, packed_measurement;
        branch_rate_scales=(T(8), T(8)))
    estimate_wfs_measurement!(packed_measurement, packed_observation,
        packed_estimator)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(curvature_measurement.storage))
    @test curvature_measurement.storage isa BackendArray
    @test packed_measurement.storage isa BackendArray
    @test isapprox(Array(curvature_measurement.storage),
        Array(packed_measurement.storage); rtol=T(5e-5), atol=T(5e-5))
    return nothing
end

function run_optional_plane_product_checks(tel::Telescope,
    src::AdaptiveOpticsSim.Optics.AbstractSource,
    selector::AdaptiveOpticsSim.Backends.AbstractArrayBackend, BackendArray,
    ::Type{T}) where {T<:AbstractFloat}
    wavefront = PupilFunction(tel; T=T, backend=selector)
    field = ElectricField(wavefront, src; zero_padding=2, T=T)
    formation = prepare_pupil_field(wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test wavefront.amplitude isa BackendArray
    @test wavefront.opd isa BackendArray
    @test field.values isa BackendArray
    @test field.metadata.device == compute_device(field.values)
    field_view = @view field.values[:, :]
    for wrapper in (
        field_view,
        reshape(field_view, 1, length(field_view)),
        transpose(field_view),
        PermutedDimsArray(field_view, (2, 1)),
    )
        @test compute_device(wrapper) == field.metadata.device
    end

    wrapped_intensity_parent = similar(wavefront.opd, T, 4, 4)
    fill!(wrapped_intensity_parent, one(T))
    wrapped_intensity_view = @view wrapped_intensity_parent[:, :]
    for wrapper in (
        wrapped_intensity_view,
        reshape(wrapped_intensity_view, 1, length(wrapped_intensity_view)),
        transpose(wrapped_intensity_view),
        PermutedDimsArray(wrapped_intensity_view, (2, 1)),
    )
        @test typeof(backend(wrapper)) === typeof(selector)
        @test typeof(AdaptiveOpticsSim.Backends.array_backend_selector(
            typeof(wrapper))) === typeof(selector)
        wrapper_metadata = OpticalPlaneMetadata(FocalPlane(), wrapper;
            coordinate_domain=AngularCoordinates(),
            sampling=(one(T), one(T)),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition())
        wrapper_map = IntensityMap(wrapper_metadata, wrapper)
        @test wrapper_map.values === wrapper
        @test wrapper_map.metadata.device == field.metadata.device
    end

    prepared = prepare_direct_imaging(wavefront, src; zero_padding=2)
    formed = form_direct_image!(prepared)
    @test formed.values isa BackendArray
    @test formed.metadata.device == wavefront.metadata.device
    @test sum(Array(formed.values)) ≈
        sum(Array(pupil_photon_rate_map(tel, src))) atol=T(2e-5) rtol=T(2e-5)

    off_axis_src = Source(band=:I, magnitude=zero(T),
        coordinates=(T(0.08), T(90)), T=T)
    off_axis = prepare_direct_imaging(wavefront, off_axis_src; zero_padding=2)
    off_axis_map = form_direct_image!(off_axis)
    @test off_axis_map.values isa BackendArray
    @test off_axis_map.metadata.device == formed.metadata.device
    @test off_axis.plan.shift_samples isa NTuple{2,Int}
    @test off_axis.plan.shift_samples != (0, 0)
    @test Array(off_axis_map.values) ≈ circshift(Array(formed.values),
        off_axis.plan.shift_samples) rtol=T(2e-5) atol=T(2e-5)

    spectral_src = with_spectrum(src,
        SpectralBundle(T[0.7e-6, 0.9e-6], T[0.4, 0.6]; T=T))
    spectral = prepare_direct_imaging(wavefront, spectral_src; zero_padding=2)
    spectral_products = form_direct_image!(spectral)
    @test spectral_products isa OpticalProductBundle
    @test length(spectral_products) == 2
    @test all(product -> product.values isa BackendArray,
        spectral_products)
    @test all(product -> product.metadata.device == formed.metadata.device,
        spectral_products)
    @test spectral_products[1].metadata.sampling !=
        spectral_products[2].metadata.sampling

    extended_src = with_extended_source(src,
        GaussianDiskSourceModel(sigma_arcsec=T(0.02), n_side=5, T=T))
    extended = prepare_direct_imaging(wavefront,
        extended_source_asterism(extended_src); zero_padding=1)
    @test extended.components isa Vector
    @test extended.products isa Vector
    @test isconcretetype(eltype(extended.components))
    @test isconcretetype(eltype(extended.products))
    @test length(extended.components) == 25
    extended_map = form_direct_image!(extended)
    @test extended_map.values isa BackendArray
    @test all(product -> product.values isa BackendArray,
        extended.products)
    @test all(product -> product.metadata.device == formed.metadata.device,
        extended.products)
    @test sum(Array(extended_map.values)) ≈
        sum(Array(pupil_photon_rate_map(tel, src))) atol=T(2e-5) rtol=T(2e-5)

    host_values = zeros(T, size(formed.values))
    host_metadata = OpticalPlaneMetadata(FocalPlane(), host_values;
        coordinate_domain=AngularCoordinates(),
        sampling=formed.metadata.sampling,
        origin=formed.metadata.origin,
        spectral=formed.metadata.spectral,
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    host_output = IntensityMap(host_metadata, host_values)
    @test_throws InvalidConfiguration prepare_direct_imaging(
        wavefront, src, prepared.field, host_output)

    sum_output_values = similar(formed.values)
    first_sum_values = similar(formed.values)
    second_sum_values = similar(formed.values)
    fill!(first_sum_values, one(T))
    fill!(second_sum_values, T(2))
    sum_output = IntensityMap(formed.metadata, sum_output_values)
    first_sum_input = IntensityMap(formed.metadata, first_sum_values)
    second_sum_input = IntensityMap(formed.metadata, second_sum_values)
    sum_plan = prepare_incoherent_sum(sum_output, first_sum_input,
        second_sum_input)
    accumulate_intensity!(sum_output,
        (first_sum_input, second_sum_input), sum_plan)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(sum_output.values))
    @test sum_output.values isa BackendArray
    @test Array(sum_output.values) == fill(T(3), size(sum_output.values))

    detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    acquisition = prepare_detector_acquisition(detector, prepared.output)
    detector_frame = capture!(detector, prepared.output, acquisition;
        rng=MersenneTwister(301))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(detector_frame))
    @test detector_frame isa BackendArray
    @test sum(Array(detector_frame)) ≈
        sum(Array(prepared.output.values)) * T(0.25) atol=T(2e-5) rtol=T(2e-5)
    identical_detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    @test identical_detector.params === detector.params
    @test identical_detector.state !== detector.state
    @test_throws InvalidConfiguration capture!(identical_detector,
        prepared.output, acquisition; rng=MersenneTwister(301))

    long_detector = Detector(integration_time=T(1.0), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    long_acquisition = prepare_detector_acquisition(long_detector,
        prepared.output)
    @test acquisition.input_values === prepared.output.values
    @test long_acquisition.input_values === prepared.output.values
    short_snapshot = copy(Array(detector_frame))
    long_frame = capture!(long_detector, prepared.output,
        long_acquisition; rng=MersenneTwister(307))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(long_frame))
    @test Array(long_frame) ≈ 2 .* short_snapshot atol=T(2e-5) rtol=T(2e-5)

    incremental_detector = Detector(integration_time=T(0.5),
        noise=NoiseNone(), qe=T(0.5), response_model=NullFrameResponse(),
        T=T, backend=selector)
    incremental_acquisition = prepare_detector_acquisition(
        incremental_detector, prepared.output)
    capture!(incremental_detector, prepared.output,
        incremental_acquisition; rng=MersenneTwister(305),
        integration_duration=T(0.2))
    @test !readout_ready(incremental_detector)
    incremental_frame = capture!(incremental_detector, prepared.output,
        incremental_acquisition; rng=MersenneTwister(306),
        integration_duration=T(0.3))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(incremental_frame))
    @test incremental_frame isa BackendArray
    @test readout_ready(incremental_detector)
    @test sum(Array(incremental_frame)) ≈
        sum(Array(prepared.output.values)) * T(0.25) atol=T(2e-5) rtol=T(2e-5)
    cpu_detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=CPUBackend())
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        cpu_detector, prepared.output)
    @test_throws InvalidConfiguration capture!(cpu_detector,
        prepared.output, acquisition; rng=MersenneTwister(306))

    density_host = zeros(T, 9, 9)
    density_host[3, 5] = T(8)
    density_values = BackendArray(density_host)
    density_metadata = OpticalPlaneMetadata(FocalPlane(), density_values;
        coordinate_domain=AngularCoordinates(), sampling=(T(0.5), T(0.25)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=AdaptiveOpticsSim.Optics.PhotonRateNormalization(),
        spatial_measure=AdaptiveOpticsSim.Optics.SpatialDensityMeasure(),
        coherence=AdaptiveOpticsSim.Optics.IncoherentIntensityAddition())
    density_map = IntensityMap(density_metadata, density_values)
    response_kernel = T[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0]
    response_detector = Detector(integration_time=T(2), noise=NoiseNone(),
        qe=T(0.5), binning=3,
        response_model=SampledFrameResponse(response_kernel; T=T), T=T,
        backend=selector)
    response_acquisition = prepare_detector_acquisition(response_detector,
        density_map)
    response_frame = capture!(response_detector, density_map,
        response_acquisition; rng=MersenneTwister(302))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(response_frame))
    manual_response = zeros(T, 9, 9)
    manual_response[3, 5] = T(4.8)
    manual_response[2, 5] = T(0.8)
    manual_response[4, 5] = T(0.8)
    manual_response[3, 4] = T(0.8)
    manual_response[3, 6] = T(0.8)
    manual_binned = zeros(T, 3, 3)
    AdaptiveOpticsSim.bin2d!(manual_binned, manual_response, 3)
    @test Array(response_frame) ≈ manual_binned .* T(0.125) atol=T(1e-6) rtol=T(1e-6)

    asymmetric_kernel = T[0 0 0; 0.1 0.2 0.7; 0 0 0]
    edge_host = zeros(T, 9, 9)
    edge_host[5, end] = one(T)
    edge_values = BackendArray(edge_host)
    edge_metadata = OpticalPlaneMetadata(FocalPlane(), edge_values;
        coordinate_domain=AngularCoordinates(), sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=AdaptiveOpticsSim.Optics.PhotonRateNormalization(),
        spatial_measure=AdaptiveOpticsSim.Optics.CellIntegratedMeasure(),
        coherence=AdaptiveOpticsSim.Optics.IncoherentIntensityAddition())
    edge_map = IntensityMap(edge_metadata, edge_values)
    edge_detector = Detector(integration_time=one(T), noise=NoiseNone(),
        qe=one(T), response_model=SampledFrameResponse(asymmetric_kernel; T=T),
        T=T, backend=selector)
    edge_acquisition = prepare_detector_acquisition(edge_detector, edge_map)
    edge_frame = capture!(edge_detector, edge_map, edge_acquisition;
        rng=MersenneTwister(303))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(edge_frame))
    @test sum(Array(edge_frame)) ≈ T(0.9) atol=T(2e-6) rtol=T(2e-6)
    @test sum(Array(edge_frame)) <= sum(edge_host) + T(2e-6)

    edge_cube_host = zeros(T, 2, 9, 9)
    edge_cube_host[1, 5, end] = one(T)
    edge_cube_host[2, 5, 1] = one(T)
    edge_cube = BackendArray(edge_cube_host)
    edge_scratch = similar(edge_cube)
    edge_stack = AdaptiveOpticsSim.Detectors.capture_stack!(edge_detector, edge_cube,
        edge_scratch; rng=MersenneTwister(304))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(edge_stack))
    edge_stack_host = Array(edge_stack)
    @test sum(@view(edge_stack_host[1, :, :])) ≈ T(0.9) atol=T(2e-6) rtol=T(2e-6)
    @test sum(@view(edge_stack_host[2, :, :])) ≈ T(0.3) atol=T(2e-6) rtol=T(2e-6)

    fraunhofer = FraunhoferPropagation(field)
    propagated = propagation_output(field, fraunhofer)
    propagate_field!(propagated, field, fraunhofer)
    @test propagated.values isa BackendArray
    @test propagated.metadata.kind isa FocalPlane

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2, T=T, backend=selector)
    spatial_field = ElectricField(wavefront, src; zero_padding=2, T=T,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    spatial_formation = prepare_pupil_field(wavefront, src, spatial_field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(spatial_field, wavefront, spatial_formation)
    filtered = PupilFunction(tel; T=T, backend=selector)
    spatial_plan = prepare_spatial_filter(tel, spatial_filter, spatial_field,
        filtered)
    spatial_workspace = SpatialFilterWorkspace(spatial_filter)
    filter!(filtered, spatial_field, spatial_filter, spatial_plan,
        spatial_workspace)
    @test filtered.opd isa BackendArray
    @test filtered.amplitude isa BackendArray
    @test all(isfinite, Array(filtered.opd))
    @test all(isfinite, Array(filtered.amplitude))

    independent_path = PupilFunction(tel; T=T, backend=selector)
    independent_opd = copy(Array(independent_path.opd))
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    fill!(dm.state.coefs, T(1e-8))
    update_surface!(dm)
    apply_surface!(wavefront, dm, DMReplace())
    @test Array(wavefront.opd) ≈ Array(surface_opd(dm)) atol=zero(T) rtol=zero(T)
    @test Array(independent_path.opd) == independent_opd
    return nothing
end

function run_optional_scalable_reconstructor_checks(::Type{T}, selector,
    BackendArray) where {T<:AbstractFloat}
    rng = MersenneTwister(20260713)
    n_slopes = 12
    n_commands = 8
    D_host = randn(rng, T, n_slopes, n_commands)
    slopes_host = randn(rng, T, n_slopes)
    D = BackendArray(D_host)
    input = BackendArray(slopes_host)
    interaction = InteractionMatrix(D, T(0.1))
    dense = ModalReconstructor(interaction; gain=T(0.5))
    factorized = FactorizedReconstructor(interaction; gain=T(0.5))
    dense_out = BackendArray(zeros(T, n_commands))
    factorized_out = BackendArray(zeros(T, n_commands))
    reconstruct!(dense_out, dense, input)
    reconstruct!(factorized_out, factorized, input)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(factorized_out))
    @test Array(factorized_out) ≈ Array(dense_out) atol=T(2e-4) rtol=T(2e-4)
    @test AdaptiveOpticsSim.factorized_rank(factorized) == n_commands
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    compact = FactorizedReconstructor(interaction; gain=T(0.5), max_rank=3)
    @test AdaptiveOpticsSim.factorized_rank(compact) == 3
    @test sum(length,
        AdaptiveOpticsSim.runtime_reconstructor_storage(compact)) <
        sum(length,
            AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    controller = DiscreteIntegratorController(n_commands;
        gain=T(0.7), tau=T(0.01), T=T, backend=selector)
    controlled = ControlledReconstructor(factorized, controller; dt=T(1e-3))
    reconstruct!(factorized_out, controlled, input)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(factorized_out))
    @test all(isfinite, Array(factorized_out))
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(controlled))
    @test AdaptiveOpticsSim.reset_controller!(controlled) === controlled
    @test_throws InvalidConfiguration ControlledReconstructor(
        factorized,
        DiscreteIntegratorController(n_commands; T=T, backend=CPUBackend());
        dt=T(1e-3),
    )
    return nothing
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:sh}) where {T<:AbstractFloat}
    return ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:pyr}) where {T<:AbstractFloat}
    return PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:bio}) where {T<:AbstractFloat}
    return BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_optic(tel::Telescope, backend,
    ::Type{T}, ::Val{:tiptilt}) where {T<:AbstractFloat}
    return TipTiltMirror(
        tel;
        scale=T(0.1),
        T=T,
        backend=backend,
        label=:tiptilt,
    )
end

function _build_optional_low_order_optic(tel::Telescope, backend,
    ::Type{T}, ::Val{:steering}) where {T<:AbstractFloat}
    return ModalControllableOptic(
        tel,
        CartesianTiltBasis(; scale=T(0.1));
        labels=:steering,
        T=T,
        backend=backend,
    )
end

function _build_optional_low_order_optic(tel::Telescope, backend,
    ::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return FocusStage(
        tel;
        scale=T(0.1),
        T=T,
        backend=backend,
        label=:focus,
    )
end

function _optional_low_order_commands(::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return fill(T(1.25e-7), 1), fill(T(2.5e-7), 1), fill(T(2e-8), 16)
end

function _optional_low_order_commands(::Type{T}, ::Union{Val{:tiptilt},Val{:steering}}) where {T<:AbstractFloat}
    # These coefficients produce nanometre-scale OPD. Centimetre-scale test
    # commands make diffraction output chaotic under otherwise sub-ulp
    # Float32 backend differences and do not represent an AO operating point.
    return fill(T(1.25e-7), 2), fill(T(2.5e-7), 2), fill(T(2e-8), 16)
end

function _optional_low_order_tolerances(::Val{:focus}, ::Val{:sh}=Val(:sh))
    return (slopes_rtol=5f-3, slopes_atol=6f-3, frame_rtol=6f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Union{Val{:tiptilt},Val{:steering}}, ::Union{Val{:sh},Val{:pyr}})
    return (slopes_rtol=3f-3, slopes_atol=4f-3, frame_rtol=4f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Val{:tiptilt}, ::Val{:bio})
    return (slopes_rtol=1.5f-1, slopes_atol=4f-3, frame_rtol=6f-1, frame_atol=1f6)
end

function _build_optional_independent_optics_case(backend, ::Type{T},
    case::Val, wfs_case::Val) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=16,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    source = Source(band=:I, magnitude=0, T=T)
    pupil = PupilFunction(tel; T=T, backend=backend)
    low_order = _build_optional_low_order_optic(tel, backend, T, case)
    dm = DeformableMirror(
        tel;
        n_act=4,
        influence_width=T(0.3),
        T=T,
        backend=backend,
    )
    wfs = _build_optional_low_order_wfs(
        tel,
        backend,
        T,
        wfs_case,
    )
    detector = Detector(
        noise=NoiseNone(),
        integration_time=one(T),
        qe=one(T),
        binning=1,
        T=T,
        backend=backend,
    )
    prepare_runtime_wfs!(wfs, pupil, source)
    return (; pupil, source, low_order, dm, wfs, detector)
end

function _optional_independent_optics_snapshot!(prepared,
    low_order_command, dm_command)
    set_command!(prepared.low_order, low_order_command)
    set_command!(prepared.dm, dm_command)
    reset_opd!(prepared.pupil)
    for optic in (prepared.low_order, prepared.dm)
        update_surface!(optic)
        apply_surface!(prepared.pupil, optic, DMAdditive())
    end
    measure!(
        prepared.wfs,
        prepared.pupil,
        prepared.source,
        prepared.detector;
        rng=MersenneTwister(91),
    )
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(slopes(prepared.wfs)))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(output_frame(prepared.detector)))
    return (
        low_order_command=copy(AdaptiveOpticsSim.Optics.command_storage(
            prepared.low_order)),
        dm_command=copy(AdaptiveOpticsSim.Optics.command_storage(prepared.dm)),
        slopes=copy(slopes(prepared.wfs)),
        frame=copy(output_frame(prepared.detector)),
    )
end

function _run_optional_independent_optics_case(::Type{B}, case::Val{K},
    wfs_case::Val{W}=Val(:sh)) where {
    B<:AdaptiveOpticsSim.Backends.GPUBackendTag,K,W,
}
    T = Float32
    cpu = _build_optional_independent_optics_case(
        CPUBackend(),
        T,
        case,
        wfs_case,
    )
    gpu = _build_optional_independent_optics_case(
        backend_selector(B),
        T,
        case,
        wfs_case,
    )
    initial_low_order, updated_low_order, dm_cmd = _optional_low_order_commands(T, case)
    tol = _optional_low_order_tolerances(case, wfs_case)

    cpu_initial = _optional_independent_optics_snapshot!(
        cpu,
        initial_low_order,
        dm_cmd,
    )
    gpu_initial = _optional_independent_optics_snapshot!(
        gpu,
        initial_low_order,
        dm_cmd,
    )
    @test Array(gpu_initial.low_order_command) ==
        Array(cpu_initial.low_order_command)
    @test Array(gpu_initial.dm_command) == Array(cpu_initial.dm_command)
    @test isapprox(
        Array(gpu_initial.slopes),
        Array(cpu_initial.slopes);
        rtol=tol.slopes_rtol,
        atol=tol.slopes_atol,
    )
    @test isapprox(
        Array(gpu_initial.frame),
        Array(cpu_initial.frame);
        rtol=tol.frame_rtol,
        atol=tol.frame_atol,
    )

    cpu_updated = _optional_independent_optics_snapshot!(
        cpu,
        updated_low_order,
        dm_cmd,
    )
    gpu_updated = _optional_independent_optics_snapshot!(
        gpu,
        updated_low_order,
        dm_cmd,
    )
    @test Array(gpu_updated.dm_command) == Array(gpu_initial.dm_command)
    @test isapprox(
        Array(gpu_updated.slopes),
        Array(cpu_updated.slopes);
        rtol=tol.slopes_rtol,
        atol=tol.slopes_atol,
    )
    @test isapprox(
        Array(gpu_updated.frame),
        Array(cpu_updated.frame);
        rtol=tol.frame_rtol,
        atol=tol.frame_atol,
    )
    return nothing
end

function run_optional_independent_optics_parity(
    ::Type{B},
    BackendArray,
) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    _run_optional_independent_optics_case(B, Val(:tiptilt), Val(:sh))
    _run_optional_independent_optics_case(B, Val(:tiptilt), Val(:pyr))
    _run_optional_independent_optics_case(B, Val(:tiptilt), Val(:bio))
    _run_optional_independent_optics_case(B, Val(:steering), Val(:sh))
    _run_optional_independent_optics_case(B, Val(:focus), Val(:sh))
    return nothing
end

function run_optional_counting_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    correlation = ChannelCrosstalkModel(T(0.4))
    sensor_kwargs = (
        qe=T(1),
        dark_count_rate=T(0),
        fill_factor=T(1),
        energy_resolution=T(12),
        timing_jitter_s=T(2e-6),
        wavelength_range_m=(T(0.8e-6), T(1.4e-6)),
        correlation_model=correlation,
        T=T,
    )
    cpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=CPUBackend(),
    )
    gpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=selector,
    )
    input = zeros(T, 5, 5)
    input[3, 3] = T(10)
    cpu_output = capture!(cpu, input, MersenneTwister(19))
    gpu_output = capture!(gpu, BackendArray(input), MersenneTwister(19))
    AdaptiveOpticsSim.Backends.synchronize_backend!(AdaptiveOpticsSim.Backends.execution_style(gpu_output))
    @test gpu_output isa BackendArray
    @test isapprox(Array(gpu_output), cpu_output; rtol=1f-6, atol=1f-6)
    @test isapprox(sum(Array(gpu_output)), sum(cpu_output); rtol=1f-6, atol=1f-6)

    inside_band = Source(band=:custom, wavelength=T(1.0e-6), T=T)
    outside_band = Source(band=:custom, wavelength=T(0.55e-6), T=T)
    gpu_inside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), inside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.Backends.synchronize_backend!(AdaptiveOpticsSim.Backends.execution_style(gpu_inside))
    inside_host = Array(gpu_inside)
    gpu_outside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), outside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.Backends.synchronize_backend!(AdaptiveOpticsSim.Backends.execution_style(gpu_outside))
    outside_host = Array(gpu_outside)
    @test inside_host == fill(T(2), 2, 2)
    @test outside_host == zeros(T, 2, 2)
    return nothing
end

function run_optional_avalanche_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    input = BackendArray(fill(one(T), 128, 128))

    pc_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(10),
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(
            threshold=T(5), detection_efficiency=T(0.5)), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    pc_output = capture!(pc_detector, input; rng=MersenneTwister(2026))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(pc_output))
    pc_host = Array(pc_output)
    @test all(x -> x == zero(T) || x == one(T), pc_host)
    @test T(0.47) <= sum(pc_host) / length(pc_host) <= T(0.53)

    stochastic_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(5),
        sensor=EMCCDSensor(excess_noise_factor=T(1.4),
            multiplication_model=AdaptiveOpticsSim.Detectors.StochasticMultiplicationRegister(T(0.6)),
            T=T), response_model=NullFrameResponse(), T=T, backend=selector)
    stochastic_input = BackendArray(fill(T(50), 128, 128))
    stochastic_output = capture!(stochastic_detector, stochastic_input;
        rng=MersenneTwister(2027))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(stochastic_output))
    stochastic_host = Array(stochastic_output)
    @test all(x -> x >= zero(T), stochastic_host)
    @test isapprox(sum(stochastic_host) / length(stochastic_host), T(250);
        rtol=T(0.03))

    ramp_detector = Detector(integration_time=T(2), noise=NoiseNone(),
        qe=one(T), gain=one(T),
        sensor=HgCdTeAvalancheArraySensor(read_time=zero(T),
            sampling_mode=UpTheRampSampling(5), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    ramp_output = capture!(ramp_detector,
        BackendArray(fill(T(3), 32, 32)); rng=MersenneTwister(2028))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(ramp_output))
    @test Array(ramp_output) == fill(T(6), 32, 32)
    @test Array(detector_ramp_slope(ramp_detector)) == fill(T(3), 32, 32)
    @test maximum(abs, Array(detector_ramp_intercept(ramp_detector))) <= eps(T)
    @test size(detector_ramp_cube(ramp_detector)) == (32, 32, 5)
    @test detector_ramp_times(ramp_detector) == T[0, 0.5, 1, 1.5, 2]

    linear_apd = LinearAPDDetector(topology=APDChannelBank(4),
        integration_time=T(0.5), qe=T(0.5), avalanche_gain=T(4),
        conversion_gain=T(2), noise=NoiseNone(), T=T, backend=selector)
    linear_output = capture!(linear_apd, BackendArray(fill(T(10), 4));
        rng=MersenneTwister(2029))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(linear_output))
    @test linear_output isa BackendArray
    @test Array(linear_output) == fill(T(20), 4)
    return nothing
end

function optional_detector_event_map(
    values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.8e-6)))
    return IntensityMap(metadata, values)
end

function run_optional_detector_event_checks(::Type{B}, BackendArray) where
    {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    definition = Plant.GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    rng = MersenneTwister(2030)

    values = BackendArray(fill(T(2), 8, 8))
    map = optional_detector_event_map(values)
    detector = Detector(integration_time=one(T), qe=T(0.5),
        noise=NoiseNone(), response_model=NullFrameResponse(), T=T,
        backend=selector)
    prepared = Plant.prepare_global_shutter_acquisition(detector,
        map, definition)
    state = Plant.GlobalShutterAcquisitionState(prepared)
    start = PlantTimestamp(0)
    middle = PlantTimestamp(400_000_000)
    close = PlantTimestamp(1_000_000_000)
    Plant.begin_exposure!(prepared, state, start)
    Plant.accumulate_exposure_interval!(prepared, state, start,
        middle, rng)
    Plant.accumulate_exposure_interval!(prepared, state, middle,
        close, rng)
    Plant.close_exposure!(prepared, state, close)
    output = Plant.complete_readout!(prepared, state, close, rng)
    Plant.mark_acquisition_ready!(prepared, state, close)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(output))
    @test output isa BackendArray
    @test Array(output) == fill(one(T), 8, 8)

    ramp_values = BackendArray(fill(one(T), 4, 4))
    ramp_map = optional_detector_event_map(ramp_values)
    ramp_detector = Detector(integration_time=one(T), qe=one(T),
        noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(read_time=zero(T),
            sampling_mode=UpTheRampSampling(3), T=T),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:3, 2:3), T=T,
        backend=selector)
    ramp_prepared = Plant.prepare_global_shutter_acquisition(
        ramp_detector, ramp_map, definition)
    ramp_state = Plant.GlobalShutterAcquisitionState(
        ramp_prepared)
    ramp_middle = PlantTimestamp(500_000_000)
    Plant.begin_exposure!(ramp_prepared, ramp_state, start)
    Plant.take_nondestructive_read!(ramp_prepared, ramp_state,
        start, rng)
    Plant.accumulate_exposure_interval!(ramp_prepared, ramp_state,
        start, ramp_middle, rng)
    Plant.take_nondestructive_read!(ramp_prepared, ramp_state,
        ramp_middle, rng)
    fill!(ramp_values, T(3))
    Plant.accumulate_exposure_interval!(ramp_prepared, ramp_state,
        ramp_middle, close, rng)
    Plant.take_nondestructive_read!(ramp_prepared, ramp_state,
        close, rng)
    Plant.close_exposure!(ramp_prepared, ramp_state, close)
    ramp_output = Plant.complete_readout!(ramp_prepared,
        ramp_state, close, rng)
    Plant.mark_acquisition_ready!(ramp_prepared, ramp_state,
        close)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(ramp_output))
    @test ramp_output isa BackendArray
    @test detector_ramp_cube(ramp_detector) isa BackendArray
    @test Array(detector_ramp_cube(ramp_detector)) == cat(
        fill(zero(T), 2, 2), fill(T(0.5), 2, 2), fill(T(2), 2, 2);
        dims=3)
    @test Array(detector_ramp_slope(ramp_detector)) == fill(T(2), 2, 2)
    @test Array(ramp_output) == fill(T(2), 2, 2)

    rolling_values = BackendArray(fill(T(2), 8, 8))
    rolling_map = optional_detector_event_map(rolling_values)
    rolling_detector = Detector(integration_time=one(T), qe=one(T),
        noise=NoiseNone(), response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(T(0.1);
            row_group_size=2), T=T), T=T, backend=selector)
    rolling_prepared =
        Plant.prepare_rolling_shutter_acquisition(
            rolling_detector, rolling_map,
            Plant.RollingShutterAcquisitionDefinition(
                PlantDuration(1_000_000_000)))
    rolling_state = Plant.RollingShutterAcquisitionState(
        rolling_prepared)
    Plant.begin_exposure!(rolling_prepared, rolling_state,
        start)
    while true
        next_open = Plant.next_rolling_band_open_timestamp(
            rolling_prepared, rolling_state)
        next_close = Plant.next_rolling_band_close_timestamp(
            rolling_prepared, rolling_state)
        next_open === nothing && next_close === nothing && break
        timestamp = next_open === nothing ? next_close :
            next_close === nothing ? next_open : min(next_open, next_close)
        if Plant.integrated_through_timestamp(rolling_state) <
                timestamp &&
                Plant.rolling_opened_band_count(rolling_state) >
                    Plant.rolling_closed_band_count(rolling_state)
            Plant.accumulate_rolling_exposure_interval!(
                rolling_prepared, rolling_state,
                Plant.integrated_through_timestamp(
                    rolling_state), timestamp, rng)
        end
        next_close == timestamp &&
            Plant.close_next_rolling_band!(rolling_prepared,
                rolling_state, timestamp)
        next_open == timestamp &&
            Plant.open_next_rolling_band!(rolling_prepared,
                rolling_state, timestamp)
    end
    rolling_output = Plant.complete_readout!(rolling_prepared,
        rolling_state,
        Plant.readout_complete_timestamp(rolling_state), rng)
    Plant.mark_acquisition_ready!(rolling_prepared,
        rolling_state,
        Plant.acquisition_readiness_timestamp(rolling_state))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(rolling_output))
    @test rolling_output isa BackendArray
    @test Array(rolling_output) == fill(T(2), 8, 8)

    transfer_values = BackendArray(fill(T(3), 8, 8))
    transfer_map = optional_detector_event_map(transfer_values)
    transfer_detector = Detector(integration_time=one(T), qe=one(T),
        gain=one(T), noise=NoiseNone(), response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=T(0.1)), T=T), T=T, backend=selector)
    transfer_prepared =
        Plant.prepare_frame_transfer_acquisition(
            transfer_detector, transfer_map,
            Plant.FrameTransferAcquisitionDefinition(
                PlantDuration(1_000_000_000);
                readout_duration=PlantDuration(200_000_000)))
    transfer_state = Plant.FrameTransferAcquisitionState(
        transfer_prepared)
    Plant.begin_exposure!(transfer_prepared, transfer_state,
        start)
    transfer_close = Plant.exposure_close_timestamp(
        transfer_state)
    Plant.accumulate_exposure_interval!(transfer_prepared,
        transfer_state, start, transfer_close, rng)
    Plant.close_exposure!(transfer_prepared, transfer_state,
        transfer_close)
    transfer_complete = Plant.frame_transfer_complete_timestamp(
        transfer_state)
    Plant.complete_frame_transfer!(transfer_prepared,
        transfer_state, transfer_complete)
    transfer_readout = Plant.readout_complete_timestamp(
        transfer_state)
    transfer_output = Plant.complete_readout!(transfer_prepared,
        transfer_state, transfer_readout, rng)
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(transfer_output))
    @test transfer_prepared.storage_frame isa BackendArray
    @test transfer_output isa BackendArray
    @test Array(transfer_output) == fill(T(3), 8, 8)
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag})
    @eval import CUDA
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag})
    @eval import AMDGPU
    return nothing
end

function backend_functional(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :CUDA), :functional))
end

function backend_functional(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :AMDGPU), :functional))
end

run_optional_backend_plan_checks(::Type{<:AdaptiveOpticsSim.Backends.GPUBackendTag}, tel, backend) = nothing

function run_optional_lift_fallback_check(array_backend, ::Type{T}) where {T<:AbstractFloat}
    H_host = T[1 0; 0 2; 1 1]
    residual_host = T[2, -1, 0.5]
    H = array_backend(H_host)
    residual = array_backend(residual_host)
    rhs = array_backend(zeros(T, 2))
    damping = LiFTLevenbergMarquardt(lambda0=T(0.1),
        growth=T(10), condition_rtol=T(1e-3))
    normal = transpose(H_host) * H_host
    λ = AdaptiveOpticsSim.WavefrontSensors.damping_lambda(damping, normal)
    expected = (normal + λ * I) \ (transpose(H_host) * residual_host)
    diag = AdaptiveOpticsSim.WavefrontSensors.LiFTDiagnostics(
        T(NaN), T(NaN), T(NaN), T(NaN), zero(T), false, false)
    AdaptiveOpticsSim.WavefrontSensors.solve_lift_fallback!(
        diag, rhs, H, residual, damping)
    @test Array(rhs) ≈ expected rtol=T(1e-4) atol=T(1e-5)
    @test diag.regularization == λ
    @test diag.used_fallback
    return nothing
end

function run_optional_lift_pipeline_checks(::Type{B}, array_backend,
    ::Type{T}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag,T<:AbstractFloat}
    selector = backend_selector(B)
    lift_src = Source(band=:I, magnitude=zero(T), T=T)
    cpu_tel = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    lift_tel = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)
    zernike = ZernikeBasis(cpu_tel, 5; T=T)
    compute_zernike!(zernike, cpu_tel)
    basis_host = copy(@view zernike.modes[:, :, 2:4])
    diversity_host = T(30e-9) .* @view(zernike.modes[:, :, 5])
    truth_coefficients = T[8e-9, -4e-9, 2e-9]
    truth_opd_host = copy(diversity_host)
    @inbounds for mode_id in eachindex(truth_coefficients)
        @views @. truth_opd_host += truth_coefficients[mode_id] *
            basis_host[:, :, mode_id]
    end
    lift_basis = array_backend(basis_host)
    lift_diversity = array_backend(copy(diversity_host))
    truth_opd = array_backend(truth_opd_host)
    lift_object_kernel_host = T[0 1 0; 1 4 1; 0 1 0]
    lift_object_kernel = array_backend(lift_object_kernel_host)

    cpu_forward = prepare_lift_forward_model(cpu_tel, lift_src,
        basis_host; diversity_opd=copy(diversity_host), focal_resolution=8,
        object_kernel=lift_object_kernel_host)
    lift_forward = prepare_lift_forward_model(lift_tel, lift_src,
        lift_basis; diversity_opd=lift_diversity, focal_resolution=8,
        object_kernel=lift_object_kernel)
    cpu_rate = copy(intensity_values(evaluate_lift_forward!(cpu_forward,
        truth_opd_host)))
    lift_rate = copy(intensity_values(evaluate_lift_forward!(lift_forward,
        truth_opd)))
    @test lift_rate isa array_backend
    @test Array(lift_rate) ≈ cpu_rate rtol=T(2e-4) atol=T(1e-3)

    exposure = T(0.002)
    quantum_efficiency = T(0.8)
    rate_domain = LiFTPhotonRate(
        noise_equivalent_exposure_s=exposure,
        quantum_efficiency=quantum_efficiency)
    count_domain = LiFTExpectedCounts(exposure;
        quantum_efficiency=quantum_efficiency)
    photon_rate_per_unit = T(sum(cpu_rate))
    normalized_domain = LiFTNormalizedIntensity(photon_rate_per_unit;
        noise_equivalent_exposure_s=exposure,
        quantum_efficiency=quantum_efficiency)
    count_values = similar(lift_rate)
    normalized_values = similar(lift_rate)
    @. count_values = lift_rate * exposure * quantum_efficiency
    @. normalized_values = lift_rate / photon_rate_per_unit
    rate_observation = LiFTObservation(lift_forward, lift_rate;
        domain=rate_domain)
    count_observation = LiFTObservation(lift_forward, count_values;
        domain=count_domain)
    normalized_observation = LiFTObservation(lift_forward,
        normalized_values; domain=normalized_domain)

    for numerical in (false, true)
        rate_estimator = LiFT(lift_forward; iterations=2,
            mode_ids=(1, 2), numerical,
            solve_mode=LiFTSolveNormalEquations())
        count_estimator = LiFT(lift_forward; iterations=2,
            mode_ids=(1, 2), numerical,
            solve_mode=LiFTSolveNormalEquations())
        normalized_estimator = LiFT(lift_forward; iterations=2,
            mode_ids=(1, 2), numerical,
            solve_mode=LiFTSolveNormalEquations())
        rate_coefficients = WavefrontSensors.reconstruct(
            rate_estimator, rate_observation;
            optimize_norm=:none, check_convergence=false)
        count_coefficients = WavefrontSensors.reconstruct(count_estimator,
            count_observation; optimize_norm=:none,
            check_convergence=false)
        normalized_coefficients = WavefrontSensors.reconstruct(
            normalized_estimator,
            normalized_observation; optimize_norm=:none,
            check_convergence=false)
        @test rate_coefficients isa array_backend
        @test all(isfinite, Array(rate_coefficients))
        @test isapprox(Array(count_coefficients), Array(rate_coefficients);
            rtol=T(5e-4), atol=T(1e-11))
        @test isapprox(Array(normalized_coefficients),
            Array(rate_coefficients); rtol=T(5e-4), atol=T(1e-11))
    end

    cpu_analytic = LiFT(cpu_forward; iterations=2, mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    gpu_analytic = LiFT(lift_forward; iterations=2, mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    cpu_H = AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix(
        cpu_analytic,
        zeros(T, 3))
    gpu_H = AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix(
        gpu_analytic,
        array_backend(zeros(T, 3)))
    @test Array(gpu_H) ≈ cpu_H rtol=T(5e-4) atol=T(1e3)

    predicted_counts = similar(lift_rate)
    predict_lift_observation!(predicted_counts, lift_forward, truth_opd,
        count_domain)
    @test isapprox(Array(predicted_counts), Array(count_values);
        rtol=T(2e-5), atol=T(1e-3))
    @test_throws InvalidConfiguration LiFTObservation(lift_forward,
        Array(lift_rate))
    @test_throws InvalidConfiguration prepare_lift_forward_model(lift_tel,
        lift_src, basis_host; diversity_opd=lift_diversity,
        focal_resolution=8)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct(
        gpu_analytic,
        rate_observation; R_n=ones(T, 8, 8))
    host_lift_coefficients = zeros(T, 3)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        host_lift_coefficients,
        gpu_analytic, rate_observation)
    @test host_lift_coefficients == zeros(T, 3)
    device_lift_coefficients = similar(lift_basis, T, 3)
    fill!(device_lift_coefficients, zero(T))
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        device_lift_coefficients,
        gpu_analytic, rate_observation; coeffs0=zeros(T, 3))
    @test Array(device_lift_coefficients) == zeros(T, 3)

    convolution_source_host = reshape(T.(1:64), 8, 8)
    convolution_source = array_backend(convolution_source_host)
    dense_convolution = similar(convolution_source)
    dense_expected = similar(convolution_source_host)
    AdaptiveOpticsSim.WavefrontSensors.conv2d_same!(dense_convolution,
        convolution_source, lift_object_kernel)
    AdaptiveOpticsSim.WavefrontSensors.conv2d_same!(dense_expected,
        convolution_source_host, lift_object_kernel_host)
    @test Array(dense_convolution) ≈ dense_expected rtol=T(1e-5) atol=T(1e-5)

    row_kernel_host = T[1, 2, 1]
    col_kernel_host = T[1, 0.5, 1]
    row_kernel = array_backend(row_kernel_host)
    col_kernel = array_backend(col_kernel_host)
    separable_convolution = similar(convolution_source)
    separable_scratch = similar(convolution_source)
    separable_expected = similar(convolution_source_host)
    separable_expected_scratch = similar(convolution_source_host)
    AdaptiveOpticsSim.WavefrontSensors.conv2d_same_separable!(
        separable_convolution,
        separable_scratch, convolution_source, row_kernel, col_kernel)
    AdaptiveOpticsSim.WavefrontSensors.conv2d_same_separable!(
        separable_expected,
        separable_expected_scratch, convolution_source_host,
        row_kernel_host, col_kernel_host)
    @test isapprox(Array(separable_convolution), separable_expected;
        rtol=T(1e-5), atol=T(1e-5))
    run_optional_lift_fallback_check(array_backend, T)
    return nothing
end

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.Backends.AMDGPUBackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(backend)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    sh_large = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    det_capture = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, bits=12, full_well=T(100),
        sensor=CMOSSensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    pupil = PupilFunction(tel; T=T, backend=backend)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, pupil, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, pupil, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.WavefrontSensors.grouped_accumulation_plan(AdaptiveOpticsSim.Backends.execution_style(pyr.front_end.propagation.intensity), pyr) isa AdaptiveOpticsSim.WavefrontSensors.GroupedStaged2DPlan
    @test AdaptiveOpticsSim.WavefrontSensors.grouped_accumulation_plan(AdaptiveOpticsSim.Backends.execution_style(bio.front_end.propagation.intensity), bio) isa AdaptiveOpticsSim.WavefrontSensors.GroupedStaged2DPlan
    @test typeof(WavefrontSensors.sh_sensing_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(slopes(sh)), sh)) ===
        WavefrontSensors.ShackHartmannWFSRocmHostStatsPlan
    @test typeof(WavefrontSensors.sh_sensing_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(slopes(sh_large)), sh_large)) ===
        WavefrontSensors.ShackHartmannWFSRocmHostStatsPlan
    WavefrontSensors.prepare_sampling!(sh, pupil, src)
    sh_sub = div(tel.params.resolution, WavefrontSensors.n_lenslets(sh))
    sh_pad = size(sh.front_end.propagation.field, 1)
    sh_offset = div(sh_pad - sh_sub, 2)
    safe_intensity = WavefrontSensors.compute_intensity_safe!(
        AdaptiveOpticsSim.Backends.execution_style(sh.front_end.propagation.intensity),
        sh, pupil, src, 1, 1, sh_sub, sh_sub, sh_offset, sh_offset,
        sh_sub)
    @test safe_intensity === sh.front_end.propagation.intensity
    @test all(isfinite, Array(safe_intensity))
    @test AdaptiveOpticsSim.Detectors.detector_execution_plan(typeof(AdaptiveOpticsSim.Backends.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.Detectors.DetectorHostMirrorPlan
    capture_psf = array_backend{T}(undef, 4, 4)
    fill!(capture_psf, T(10))
    captured = capture!(det_capture, capture_psf; rng=MersenneTwister(2))
    @test captured isa array_backend
    @test maximum(Array(captured)) <= Float64(exp2(T(12)) - one(T))
    cpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=CPUBackend())
    gpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=backend)
    poisson_input = fill(T(10), 4, 4)
    cpu_poisson = capture!(cpu_poisson_det, copy(poisson_input); rng=MersenneTwister(23))
    gpu_poisson = capture!(gpu_poisson_det, array_backend(copy(poisson_input));
        rng=MersenneTwister(23))
    @test Array(gpu_poisson) == cpu_poisson
    poisson_method = which(
        AdaptiveOpticsSim.Detectors._poisson_noise_frame!,
        (
            AdaptiveOpticsSim.Detectors.DetectorHostMirrorPlan,
            typeof(gpu_poisson_det),
            typeof(MersenneTwister(1)),
            typeof(gpu_poisson_det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(poisson_method.file))
    @test AdaptiveOpticsSim.Backends.reduction_execution_plan(pyr.front_end.propagation.intensity) isa AdaptiveOpticsSim.Backends.HostMirrorReductionPlan
    @test AdaptiveOpticsSim.Backends.backend_sum_value(capture_psf) == T(160)

    phase_freqs = T[-0.2, -0.1, 0.1, 0.2]
    cpu_phase_psd = zeros(T, 4, 4)
    gpu_phase_freqs = array_backend(phase_freqs)
    gpu_phase_psd = array_backend(zeros(T, 4, 4))
    phase_args = (T(0.02), T(4pi^2), T(0.01), T(-11 / 6), zero(T), 4)
    AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(AdaptiveOpticsSim.Backends.ScalarCPUStyle(),
        cpu_phase_psd, phase_freqs, phase_args...)
    phase_style = AdaptiveOpticsSim.Backends.execution_style(gpu_phase_psd)
    AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(phase_style, gpu_phase_psd,
        gpu_phase_freqs, phase_args...)
    @test isapprox(Array(gpu_phase_psd), cpu_phase_psd; rtol=1f-6, atol=1f-7)
    phase_method = which(
        AdaptiveOpticsSim.Atmospheres._fill_phase_psd!,
        (typeof(phase_style), typeof(gpu_phase_psd), typeof(gpu_phase_freqs),
            T, T, T, T, T, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(phase_method.file))

    host_opd = reshape(T.(1:256), 16, 16) .* T(1e-9)
    host_valid = trues(4, 4)
    cpu_geometric_slopes = zeros(T, 32)
    gpu_geometric_slopes = array_backend(zeros(T, 32))
    gpu_opd = array_backend(host_opd)
    gpu_valid = array_backend(host_valid)
    AdaptiveOpticsSim.geometric_slopes!(cpu_geometric_slopes, host_opd, host_valid)
    AdaptiveOpticsSim.geometric_slopes!(gpu_geometric_slopes, gpu_opd, gpu_valid)
    @test isapprox(Array(gpu_geometric_slopes), cpu_geometric_slopes;
        rtol=1f-6, atol=1f-8)
    geometric_method = which(
        AdaptiveOpticsSim._geometric_slopes!,
        (typeof(AdaptiveOpticsSim.Backends.execution_style(gpu_geometric_slopes)),
            typeof(gpu_geometric_slopes), typeof(gpu_opd), typeof(gpu_valid),
            Int, Int, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(geometric_method.file))
    randn_method = which(
        AdaptiveOpticsSim.Detectors._randn_frame_noise!,
        (
            AdaptiveOpticsSim.Detectors.DetectorHostMirrorPlan,
            typeof(det),
            typeof(MersenneTwister(1)),
            typeof(det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(randn_method.file))
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.LayeredFresnelFieldAsyncPlan
    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.Detectors.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.Detectors.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.Detectors.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.Detectors.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.Detectors.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.Detectors.detector_read_times(cpu_windowed_det)

    gsc_mask = array_backend(fill(one(T), 8, 8))
    gsc_basis = array_backend(reshape(T.(1:192), 8, 8, 3) .* T(1e-3))
    gsc_frame = array_backend(reshape(T.(1:64), 8, 8))
    gsc = GainSensingCamera(gsc_mask, gsc_basis; T=T)
    calibrate!(gsc, gsc_frame)
    optical_gains = compute_optical_gains!(gsc, gsc_frame)
    @test optical_gains isa array_backend
    @test all(isfinite, Array(optical_gains))

    return nothing
end

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.Backends.CUDABackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(backend)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    pupil = PupilFunction(tel; T=T, backend=backend)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, pupil, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, pupil, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.WavefrontSensors.grouped_accumulation_plan(AdaptiveOpticsSim.Backends.execution_style(pyr.front_end.propagation.intensity), pyr) isa AdaptiveOpticsSim.WavefrontSensors.GroupedStackReducePlan
    @test AdaptiveOpticsSim.WavefrontSensors.grouped_accumulation_plan(AdaptiveOpticsSim.Backends.execution_style(bio.front_end.propagation.intensity), bio) isa AdaptiveOpticsSim.WavefrontSensors.GroupedStackReducePlan
    @test WavefrontSensors.sh_sensing_execution_plan(AdaptiveOpticsSim.Backends.execution_style(slopes(sh)), sh) isa WavefrontSensors.ShackHartmannWFSBatchedPlan
    @test AdaptiveOpticsSim.Detectors.detector_execution_plan(typeof(AdaptiveOpticsSim.Backends.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.Detectors.DetectorDirectPlan
    @test AdaptiveOpticsSim.Backends.reduction_execution_plan(pyr.front_end.propagation.intensity) isa AdaptiveOpticsSim.Backends.DirectReductionPlan
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.Backends.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.LayeredFresnelFieldAsyncPlan
    cpu_tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_sh = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
    gpu_sh = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)
    cpu_pupil = PupilFunction(cpu_tel; T=T, backend=CPUBackend())
    gpu_pupil = PupilFunction(gpu_tel; T=T, backend=backend)
    measure!(cpu_sh, cpu_pupil, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh, gpu_pupil, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_export = Array(WavefrontSensors.sh_exported_spot_cube(cpu_sh))
    gpu_export = Array(WavefrontSensors.sh_exported_spot_cube(gpu_sh))
    cpu_frame = Array(AdaptiveOpticsSim.WavefrontSensors.wfs_output_frame(cpu_sh, cpu_det))
    gpu_frame = Array(AdaptiveOpticsSim.WavefrontSensors.wfs_output_frame(gpu_sh, gpu_det))
    @test size(gpu_export) == size(cpu_export)
    @test isapprox(gpu_export, cpu_export; rtol=1f-5, atol=1f-4)
    @test size(gpu_frame) == size(cpu_frame)
    @test isapprox(gpu_frame, cpu_frame; rtol=1f-5, atol=1f-4)
    cpu_sh_stats = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    measure!(cpu_sh_stats, cpu_pupil, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh_stats, gpu_pupil, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_peak = WavefrontSensors.sh_safe_peak_value(cpu_sh_stats.acquisition.spot_cube)
    cpu_cutoff = WavefrontSensors.centroid_threshold(cpu_sh_stats) * cpu_peak
    WavefrontSensors.sh_signal_from_spots!(cpu_sh_stats, cpu_cutoff)
    gpu_peak = WavefrontSensors.sh_safe_peak_value(gpu_sh_stats.acquisition.spot_cube)
    gpu_cutoff = WavefrontSensors.centroid_threshold(gpu_sh_stats) * gpu_peak
    WavefrontSensors.sh_signal_from_spots_device_stats!(
        AdaptiveOpticsSim.Backends.execution_style(slopes(gpu_sh_stats)),
        gpu_sh_stats,
        gpu_cutoff,
    )
    @test isapprox(Array(slopes(gpu_sh_stats)), slopes(cpu_sh_stats); rtol=1f-5, atol=1f-4)

    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.Detectors.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.Detectors.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.Detectors.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.Detectors.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.Detectors.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.Detectors.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.Detectors.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.Detectors.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.Detectors.detector_read_times(cpu_windowed_det)

    return nothing
end

function optional_direct_imaging_batch_allocation_bytes(prepared)
    form_direct_image!(prepared)
    validation_bytes = @allocated(
        AdaptiveOpticsSim.Optics.validate_direct_imaging_batch(prepared),
    )
    completed_render_bytes = @allocated form_direct_image!(prepared)
    return (; validation_bytes, completed_render_bytes)
end

function optional_storage_range_contains(parent, child)
    parent_start = UInt(pointer(parent))
    child_start = UInt(pointer(child))
    parent_stop =
        parent_start + UInt(sizeof(eltype(parent)) * length(parent))
    child_stop =
        child_start + UInt(sizeof(eltype(child)) * length(child))
    return parent_start <= child_start && child_stop <= parent_stop
end

function run_optional_direct_imaging_batch_checks(
    tel::Telescope,
    selector::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    BackendArray,
    ::Type{T},
) where {T<:AbstractFloat}
    pupil = PupilFunction(tel; T=T, backend=selector)
    host_opd = reshape(
        collect(range(T(-20e-9), T(20e-9); length=length(pupil.opd))),
        size(pupil.opd),
    )
    copyto!(pupil.opd, host_opd)

    wavelengths = T[650e-9, 800e-9, 950e-9]
    source = Source(
        band=:custom,
        wavelength=wavelengths[2],
        photon_irradiance=T(6),
        coordinates=(T(0.08), T(90)),
        T=T,
    )
    sources = with_spectrum(
        source,
        SpectralBundle(wavelengths, T[1, 2, 3]; T=T),
    )
    prepared = prepare_direct_imaging_batch(
        pupil,
        sources;
        zero_padding=2,
    )
    products = form_direct_image!(prepared)
    serial = prepare_direct_imaging(
        pupil,
        sources;
        zero_padding=2,
    )
    serial_products = form_direct_image!(serial)
    device = compute_device(pupil.opd)

    @test products === direct_imaging_output(prepared)
    @test length(products) == length(wavelengths)
    @test prepared.workspace.field_stack isa BackendArray
    @test prepared.workspace.output_stack isa BackendArray
    @test prepared.workspace.shift_axis1 isa BackendArray
    @test prepared.workspace.shift_axis2 isa BackendArray
    @test prepared.workspace.fft_plan ===
        prepared.workspace_bindings.fft_plan
    @test prepared.fields === prepared.workspace_bindings.fields
    @test products === prepared.workspace_bindings.output
    @test all(
        index -> optional_storage_range_contains(
            prepared.workspace.field_stack,
            prepared.fields[index].values,
        ),
        eachindex(wavelengths),
    )
    @test all(
        index -> optional_storage_range_contains(
            prepared.workspace.output_stack,
            products[index].values,
        ),
        eachindex(wavelengths),
    )
    @test all(
        array -> compute_device(array) == device,
        (
            prepared.workspace.field_stack,
            prepared.workspace.output_stack,
            prepared.workspace.shift_axis1,
            prepared.workspace.shift_axis2,
        ),
    )
    @test all(
        index -> compute_device(products[index].values) == device,
        eachindex(wavelengths),
    )
    @test [
        products[index].metadata.spectral.wavelength_m
        for index in eachindex(wavelengths)
    ] == wavelengths

    for index in eachindex(wavelengths)
        @test isapprox(
            Array(products[index].values),
            Array(serial_products[index].values);
            rtol=T(5e-5),
            atol=T(5e-6),
        )
        @test products[index].metadata == serial_products[index].metadata
    end

    selected = products[2]
    short_detector = Detector(
        integration_time=T(0.25),
        noise=NoiseNone(),
        qe=T(0.5),
        response_model=NullFrameResponse(),
        T=T,
        backend=selector,
    )
    long_detector = Detector(
        integration_time=one(T),
        noise=NoiseNone(),
        qe=T(0.5),
        response_model=NullFrameResponse(),
        T=T,
        backend=selector,
    )
    short_acquisition =
        prepare_detector_acquisition(short_detector, selected)
    long_acquisition =
        prepare_detector_acquisition(long_detector, selected)
    @test short_acquisition.input_values === selected.values
    @test long_acquisition.input_values === selected.values
    short_frame = capture!(
        short_detector,
        selected,
        short_acquisition;
        rng=MersenneTwister(611),
    )
    long_frame = capture!(
        long_detector,
        selected,
        long_acquisition;
        rng=MersenneTwister(612),
    )
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(long_frame),
    )
    @test short_frame isa BackendArray
    @test long_frame isa BackendArray
    @test isapprox(
        Array(long_frame),
        T(4) .* Array(short_frame);
        rtol=T(5e-5),
        atol=T(5e-6),
    )

    allocation_bytes =
        optional_direct_imaging_batch_allocation_bytes(prepared)
    @test allocation_bytes.validation_bytes == 0
    # The numerical data plane remains device-resident. The owning GPU runtime
    # may still allocate bounded host-side launch descriptors for the field,
    # FFT, and intensity submissions.
    @test allocation_bytes.completed_render_bytes <= 128 * 1024
    return nothing
end

function optional_atmosphere_direction_batch_allocation_bytes(
    prepared,
    atm,
    epoch,
)
    render_atmosphere_directions!(prepared, atm, epoch)
    validation_bytes = @allocated(
        AdaptiveOpticsSim.Atmospheres.validate_atmosphere_direction_batch(
            prepared,
            atm,
            epoch,
        ),
    )
    completed_render_bytes =
        @allocated render_atmosphere_directions!(prepared, atm, epoch)
    return (; validation_bytes, completed_render_bytes)
end

function run_optional_atmosphere_direction_batch_checks(
    atm::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    tel::Telescope,
    sources::Asterism,
    epoch::AtmosphereEpoch,
    BackendArray,
)
    T = AdaptiveOpticsSim.Atmospheres.atmosphere_numeric_type(atm)
    n = tel.params.resolution
    output = BackendArray{T}(undef, n, n, length(sources))
    prepared = prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        output,
    )
    rendered = render_atmosphere_directions!(prepared, atm, epoch)
    @test rendered === output
    @test atmosphere_direction_output(prepared) === output
    @test prepared.workspace.shift_x isa BackendArray
    @test prepared.workspace.shift_y isa BackendArray
    @test prepared.workspace.footprint_scale isa BackendArray
    @test prepared.workspace.pupil isa BackendArray
    @test compute_device(output) == compute_device(prepared.workspace.shift_x)
    @test compute_device(output) == compute_device(prepared.workspace.shift_y)
    @test compute_device(output) ==
        compute_device(prepared.workspace.footprint_scale)
    @test compute_device(output) == compute_device(prepared.workspace.pupil)
    @test prepared.params.layer_count < prepared.params.layer_count *
        prepared.params.direction_count

    host_output = Array(output)
    @inbounds for direction in eachindex(prepared.params.sources)
        pupil = PupilFunction(tel; T=T, backend=backend(tel))
        renderer = prepare_atmosphere_renderer(
            atm,
            tel,
            prepared.params.sources[direction];
            T=T,
        )
        render_atmosphere!(pupil, renderer, atm, epoch)
        @test isapprox(
            @view(host_output[:, :, direction]),
            Array(pupil.opd);
            rtol=T(2e-5),
            atol=T(2e-6),
        )
    end

    allocation_bytes =
        optional_atmosphere_direction_batch_allocation_bytes(
            prepared,
            atm,
            epoch,
        )
    @test allocation_bytes.validation_bytes == 0
    # KernelAbstractions and the owning GPU runtime allocate small host-side
    # launch descriptors even after kernel compilation. Keep that separately
    # bounded from the zero-allocation prepared-contract validation and from
    # device-resident data-plane storage.
    @test allocation_bytes.completed_render_bytes <= 64 * 1024

    host_destination = fill(T(37), n, n, length(sources))
    @test_throws InvalidConfiguration prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        host_destination,
    )
    @test all(==(T(37)), host_destination)
    return nothing
end

function run_optional_backend_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    if pkg_path === nothing
        @info "Skipping $(backend_label(B)) smoke: $(pkg).jl is not available in this environment"
        @test true
        return nothing
    end

    import_backend_package!(B)
    if !backend_functional(B)
        @info "Skipping $(backend_label(B)) smoke: backend runtime/device is not functional on this host"
        @test true
        return nothing
    end

    AdaptiveOpticsSim.Backends.disable_scalar_backend!(B)
    backend = AdaptiveOpticsSim.Backends.gpu_backend_array_type(B)
    @test backend !== nothing
    run_optional_backend_selector_smoke(B, backend)
    run_optional_lgs_convolution_normalization(B)
    run_optional_sodium_profile_wfs(B, backend)
    run_optional_zernike_normalization(B, backend)
    run_optional_wfs_stage_contracts(B, backend)
    run_optional_prepared_plant_checks(B, backend)
    run_optional_device_path_batch_checks(B, backend)
    run_optional_wfs_device_model_matrix_checks(B, backend)
    run_optional_detector_device_model_matrix_checks(B, backend)
    run_optional_detector_event_checks(B, backend)
    run_optional_command_application_checks(B, backend)
    run_optional_controller_routing_checks(B, backend)
    run_optional_cycle_averaged_modulation_checks(B, backend)

    if get(ENV, backend_full_smoke_env(B), "0") == "1"
        include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
        run_smoke_matrix = Base.invokelatest(
            getfield, Main, :run_gpu_smoke_matrix)
        Base.invokelatest(run_smoke_matrix, B)
        @test true
        return nothing
    end

    T = Float32
    rng = MersenneTwister(7)
    BackendArray = backend
    selector = backend_selector(B)

    run_optional_counting_detector_parity(B, BackendArray)
    run_optional_avalanche_detector_parity(B, BackendArray)

    tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=selector)
    src = Source(band=:I, magnitude=0.0, T=T)
    run_optional_plane_product_checks(tel, src, selector, BackendArray, T)
    run_optional_direct_imaging_batch_checks(
        tel,
        selector,
        BackendArray,
        T,
    )

    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        T=T,
        backend=selector,
    )
    epoch = advance_by!(atm, T(1e-3); rng=rng)
    renderer = prepare_atmosphere_renderer(atm, tel, src)
    atmosphere_output = PupilFunction(tel; T=T)
    render_atmosphere!(atmosphere_output, renderer, atm, epoch)
    @test atm.layers[1].generator.state.opd isa BackendArray
    @test atmosphere_output.opd isa BackendArray
    direction_sources = Asterism(AdaptiveOpticsSim.Optics.AbstractSource[
        src,
        Source(
            band=:I,
            magnitude=zero(T),
            coordinates=(T(6), T(35)),
            T=T,
        ),
        LGSSource(
            magnitude=zero(T),
            coordinates=(T(-9), T(70)),
            altitude=T(90_000),
            T=T,
        ),
        Source(
            band=:K,
            magnitude=one(T),
            coordinates=(T(5), T(120)),
            T=T,
        ),
    ])
    run_optional_atmosphere_direction_batch_checks(
        atm,
        tel,
        direction_sources,
        epoch,
        BackendArray,
    )
    pupil = PupilFunction(tel; T=T, backend=selector)
    copyto!(pupil.opd, atmosphere_output.opd)

    inf_atm = InfiniteMultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        screen_resolution=33,
        stencil_size=35,
        T=T,
        backend=selector,
    )
    infinite_epoch = advance_by!(inf_atm, T(1e-3); rng=rng)
    infinite_renderer = prepare_atmosphere_renderer(inf_atm, tel, src)
    render_atmosphere!(atmosphere_output, infinite_renderer, inf_atm,
        infinite_epoch)
    @test inf_atm.layers[1].screen.state.screen isa BackendArray
    @test atmosphere_output.opd isa BackendArray
    run_optional_atmosphere_direction_batch_checks(
        inf_atm,
        tel,
        direction_sources,
        infinite_epoch,
        BackendArray,
    )

    prop = AtmosphericFieldPropagation(atm, pupil, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=2,
        T=T)
    field = propagate_atmosphere_field!(prop, atm, epoch)
    @test field.values isa BackendArray
    intensity = atmospheric_intensity!(prop, atm, epoch)
    @test intensity isa BackendArray

    bundle = SpectralBundle(fill(wavelength(src), 2), T[0.4, 0.6]; T=T)
    poly = with_spectrum(src, bundle)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=selector)
    slopes = measure!(sh, pupil, poly)
    @test slopes isa BackendArray

    spectral_optical_sh = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), T=T, backend=selector)
    WavefrontSensors.sampled_spots_peak!(spectral_optical_sh, pupil, poly)
    spectral_optical_spots = Array(spectral_optical_sh.acquisition.spot_cube)
    spectral_qe = AdaptiveOpticsSim.Detectors.SampledQuantumEfficiency(
        T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.2, 0.8])
    spectral_exposure = T(2.5)
    spectral_detector = Detector(noise=NoiseNone(), qe=spectral_qe,
        integration_time=spectral_exposure, binning=1,
        response_model=NullFrameResponse(), T=T, backend=selector)
    spectral_detector_sh = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), T=T, backend=selector)
    WavefrontSensors.sampled_spots_peak!(spectral_detector_sh, pupil, poly,
        spectral_detector, MersenneTwister(149))
    expected_spectral_scale = spectral_exposure *
        T(AdaptiveOpticsSim.Detectors.qe_at(spectral_qe, wavelength(src)))
    @test Array(spectral_detector_sh.acquisition.spot_cube) ≈
        spectral_optical_spots .* expected_spectral_scale rtol=5e-5

    distinct = with_spectrum(src, SpectralBundle(
        T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6];
        T=T))
    @test_throws InvalidConfiguration measure!(sh, pupil, distinct)

    science_src = Source(band=:K, magnitude=1.0, coordinates=(4.0, 90.0), T=T)
    split_dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    split_det = Detector(noise=NoiseNone(), integration_time=one(T), qe=one(T),
        binning=1, T=T, backend=selector)
    science_pupil = PupilFunction(tel; T=T, backend=selector)
    science_renderer = prepare_atmosphere_renderer(
        atm,
        tel,
        science_src,
    )
    render_atmosphere!(
        science_pupil,
        science_renderer,
        atm,
        epoch,
    )
    fill!(AdaptiveOpticsSim.Optics.command_storage(split_dm), T(1e-8))
    update_surface!(split_dm)
    apply_surface!(science_pupil, split_dm, DMAdditive())
    science_imaging = prepare_direct_imaging(
        science_pupil,
        science_src;
        zero_padding=1,
    )
    science_rate = form_direct_image!(science_imaging)
    split_acquisition = prepare_detector_acquisition(split_det, science_rate)
    split_frame = capture!(
        split_det,
        science_rate,
        split_acquisition;
        rng=MersenneTwister(17),
    )
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(split_frame))
    @test split_frame isa BackendArray
    @test all(isfinite, Array(split_frame))

    shared_detector_a = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_detector_b = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_acquisition_a = prepare_detector_acquisition(
        shared_detector_a, science_rate)
    shared_acquisition_b = prepare_detector_acquisition(
        shared_detector_b, science_rate)
    shared_frame_a = capture!(
        shared_detector_a,
        science_rate,
        shared_acquisition_a;
        rng=MersenneTwister(18),
    )
    shared_frame_b = capture!(
        shared_detector_b,
        science_rate,
        shared_acquisition_b;
        rng=MersenneTwister(18),
    )
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(shared_frame_a))
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(shared_frame_b))
    @test shared_frame_a isa BackendArray
    @test shared_frame_b isa BackendArray
    @test Array(shared_frame_a) ≈ Array(shared_frame_b) atol=0 rtol=0

    run_optional_lift_pipeline_checks(B, BackendArray, T)
    run_optional_backend_plan_checks(B, tel, selector)
    run_optional_independent_optics_parity(B, BackendArray)
    run_optional_scalable_reconstructor_checks(T, selector, BackendArray)

    curv = CurvatureWFS(tel; pupil_samples=4, T=T, backend=selector)
    curv_slopes = measure!(curv, pupil, src, atm)
    @test curv_slopes isa BackendArray

    return nothing
end
