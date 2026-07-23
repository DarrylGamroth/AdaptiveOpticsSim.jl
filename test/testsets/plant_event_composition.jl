const AOSEvents = AdaptiveOpticsSim

function event_test_intensity_map(values::AbstractMatrix{T}) where
    {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.75e-6)))
    return IntensityMap(metadata, values)
end

function run_rolling_event_frame!(prepared, state, rng,
    start::PlantTimestamp)
    begin_exposure!(prepared, state, start)
    while true
        next_open = next_rolling_band_open_timestamp(prepared, state)
        next_close = next_rolling_band_close_timestamp(prepared, state)
        next_open === nothing && next_close === nothing && break
        timestamp = next_open === nothing ? next_close :
            next_close === nothing ? next_open : min(next_open, next_close)
        if integrated_through_timestamp(state) < timestamp &&
                rolling_opened_band_count(state) >
                    rolling_closed_band_count(state)
            accumulate_rolling_exposure_interval!(prepared, state,
                integrated_through_timestamp(state), timestamp, rng)
        end
        next_close == timestamp &&
            close_next_rolling_band!(prepared, state, timestamp)
        next_open == timestamp &&
            open_next_rolling_band!(prepared, state, timestamp)
    end
    output = complete_readout!(prepared, state,
        readout_complete_timestamp(state), rng)
    mark_acquisition_ready!(prepared, state,
        acquisition_readiness_timestamp(state))
    return output
end

function run_frame_transfer_event!(prepared, state, rng,
    start::PlantTimestamp)
    close = start + prepared.definition.exposure_duration
    transfer = close + prepared.transfer_duration
    readout = transfer + prepared.definition.readout_duration
    begin_exposure!(prepared, state, start)
    accumulate_exposure_interval!(prepared, state, start, close, rng)
    close_exposure!(prepared, state, close)
    complete_frame_transfer!(prepared, state, transfer)
    return complete_readout!(prepared, state, readout, rng)
end

function event_test_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

struct EventCompositionPathModel{R}
    zero_padding::Int
    revision::R
end

struct EventCompositionPathExecution{E,C}
    imaging::E
    executions::C
end

abstract type EventCompositionSensorKind end
struct EventGlobalCMOS <: EventCompositionSensorKind end
struct EventRollingCMOS <: EventCompositionSensorKind end
struct EventCCD <: EventCompositionSensorKind end
struct EventFrameTransferEMCCD <: EventCompositionSensorKind end
struct EventHgCdTeRamp <: EventCompositionSensorKind end

struct EventCompositionAcquisitionModel{
    T<:AbstractFloat,K<:EventCompositionSensorKind}
    exposure::T
    kind::K
end

Plant.plant_model_definition_style(
    ::Type{<:EventCompositionPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:EventCompositionAcquisitionModel}) = ColdPlantModelDefinition()

function Plant.validate_path_execution_binding(
    execution::EventCompositionPathExecution, input, result)
    return Plant.validate_path_execution_binding(
        execution.imaging, input, result)
end

function Plant.execute_path!(result, input,
    execution::EventCompositionPathExecution)
    Plant.validate_path_execution_binding(execution, input,
        result)
    execution.executions[] += 1
    return Plant.execute_path!(result, input,
        execution.imaging)
end

function Plant.prepare_path_executor(
    model::EventCompositionPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    execution = EventCompositionPathExecution(imaging, Ref(0))
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        execution;
        materialization=prepare_pupil_opd_materialization(atmosphere,
            telescope, source, pupil),
        optical_model=(kind=:event_composition_direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=model.revision,
    )
end

@inline event_composition_sensor(::EventGlobalCMOS, ::Type{T}) where {T} =
    CMOSSensor(timing_model=GlobalShutter(), T=T)
@inline event_composition_sensor(::EventRollingCMOS, ::Type{T}) where {T} =
    CMOSSensor(timing_model=RollingShutter(T(0.005)), T=T)
@inline event_composition_sensor(::EventCCD, ::Type{T}) where {T} =
    CCDSensor(T=T)
@inline event_composition_sensor(::EventFrameTransferEMCCD,
    ::Type{T}) where {T} = EMCCDSensor(
    acquisition_mode=FrameTransferAcquisition(transfer_time=T(0.02), T=T),
    T=T)
@inline event_composition_sensor(::EventHgCdTeRamp, ::Type{T}) where {T} =
    HgCdTeAvalancheArraySensor(sampling_mode=UpTheRampSampling(3),
        read_time=zero(T), T=T)

function Plant.prepare_acquisition_provider(
    model::EventCompositionAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), gain=one(T),
        response_model=NullFrameResponse(),
        sensor=event_composition_sensor(model.kind, T), T=T,
        backend=path.key.backend)
    execution = FrameAcquisitionExecution(detector, path.result)
    metadata = (kind=:event_composition_frame,
        units=:detected_electrons, geometry=path.result.metadata,
        detector=detector_export_metadata(detector))
    products = AcquisitionProducts(execution.observation; metadata)
    return prepare_full_optical_provider(execution, products)
end

function event_composition_fixture(; reverse_order::Bool=false,
    unbound_trigger_consumer::Bool=false,
    faulted_trigger_fanout::Bool=false)
    T = Float64
    telescope = Telescope(resolution=4, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[1], wind_speed=T[7], wind_direction=T[35],
        altitude=T[0], layer_ids=(:ground,), T=T)
    science_source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(80), coordinates=(T(0), T(0)), T=T)
    ngs_source = Source(band=:custom, wavelength=T(0.7e-6),
        photon_irradiance=T(65), coordinates=(T(2), T(35)), T=T)
    lgs_source = LGSSource(wavelength=T(589e-9),
        photon_irradiance=T(70), coordinates=(T(-3), T(80)),
        altitude=T(90_000), T=T)
    paths = (
        OpticalPathDefinition(:science, science_source,
            EventCompositionPathModel(1, UInt(1))),
        OpticalPathDefinition(:ngs, ngs_source,
            EventCompositionPathModel(1, UInt(2))),
        OpticalPathDefinition(:lgs, lgs_source,
            EventCompositionPathModel(1, UInt(3))),
    )
    acquisitions = (
        AcquisitionDefinition(:science_cmos, :science,
            EventCompositionAcquisitionModel(T(0.2), EventGlobalCMOS())),
        AcquisitionDefinition(:science_ccd, :science,
            EventCompositionAcquisitionModel(T(0.3), EventCCD())),
        AcquisitionDefinition(:science_rolling, :science,
            EventCompositionAcquisitionModel(T(0.15), EventRollingCMOS())),
        AcquisitionDefinition(:ngs_saphira, :ngs,
            EventCompositionAcquisitionModel(T(0.2), EventHgCdTeRamp())),
        AcquisitionDefinition(:lgs_emccd, :lgs,
            EventCompositionAcquisitionModel(T(0.2),
                EventFrameTransferEMCCD())),
    )
    plant = prepare_plant(PlantDefinition(; telescope, atmosphere,
        paths, acquisitions);
        run_seed=0x7900)

    trigger_faults = faulted_trigger_fanout ? TriggerFaultTrace(
        TriggerFaultTraceEntry(2, :common_trigger_drop;
            action=DropTriggerEdge),
        TriggerFaultTraceEntry(3, :common_trigger_phase_error;
            phase_step=PlantTimeOffset(20_000_000),
            jitter=PlantTimeOffset(10_000_000),
            timestamp_label_offset=PlantTimeOffset(7_000_000)),
        TriggerFaultTraceEntry(4, :common_trigger_duplicate;
            action=DuplicateTriggerEdge,
            duplicate_delay=PlantDuration(250_000_000)),
    ) : TriggerFaultTrace()
    trigger_source = TriggerSourceDefinition(:ngs_camera_trigger,
        PeriodicSchedule(
            period_ns=faulted_trigger_fanout ? 600_000_000 : 400_000_000,
            phase_ns=50_000_000); faults=trigger_faults)
    trigger_consumer = TriggerConsumerDefinition(:ngs_camera,
        TriggerSourceID(:ngs_camera_trigger))
    trigger_consumers = faulted_trigger_fanout ? (
        trigger_consumer,
        TriggerConsumerDefinition(:science_camera,
            TriggerSourceID(:ngs_camera_trigger)),
    ) : unbound_trigger_consumer ? (
        trigger_consumer,
        TriggerConsumerDefinition(:unbound_camera,
            TriggerSourceID(:ngs_camera_trigger)),
    ) : (trigger_consumer,)
    trigger_topology = prepare_trigger_topology((trigger_source,), (),
        trigger_consumers;
        in_flight_capacity=(unbound_trigger_consumer ||
            faulted_trigger_fanout) ? 8 : 1)
    samples = (
        OpticalSampleDefinition(:science,
            PeriodicSchedule(period_ns=100_000_000, phase_ns=0)),
        OpticalSampleDefinition(:ngs,
            PeriodicSchedule(period_ns=125_000_000, phase_ns=50_000_000)),
        OpticalSampleDefinition(:lgs,
            PeriodicSchedule(period_ns=200_000_000, phase_ns=25_000_000)),
    )
    events = (
        DetectorEventDefinition(:science_cmos,
            GlobalShutterAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(20_000_000),
                readiness_delay=PlantDuration(10_000_000)),
            faulted_trigger_fanout ?
                TriggeredAcquisitionStart(:science_camera) :
                PeriodicAcquisitionStart(PeriodicSchedule(
                    period_ns=500_000_000, phase_ns=0))),
        DetectorEventDefinition(:science_ccd,
            GlobalShutterAcquisitionDefinition(PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=600_000_000, phase_ns=50_000_000))),
        DetectorEventDefinition(:science_rolling,
            RollingShutterAcquisitionDefinition(PlantDuration(150_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=500_000_000, phase_ns=250_000_000))),
        DetectorEventDefinition(:ngs_saphira,
            GlobalShutterAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(10_000_000)),
            TriggeredAcquisitionStart(:ngs_camera)),
        DetectorEventDefinition(:lgs_emccd,
            FrameTransferAcquisitionDefinition(PlantDuration(200_000_000);
                readout_duration=PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=350_000_000, phase_ns=25_000_000))),
    )
    ordered_samples = reverse_order ? reverse(samples) : samples
    ordered_events = reverse_order ? reverse(events) : events
    definition = PlantEventLoopDefinition(ordered_samples, ordered_events;
        trigger_topology)
    prepared = prepare_plant_event_loop(plant, definition)
    return plant, prepared, PlantEventLoopState(prepared),
        PlantEventLoopWorkspace(prepared)
end

@inline function run_event_composition_window!(prepared, state, workspace,
    stop::PlantTimestamp)
    return run_plant_events_until!(prepared, state, workspace, stop)
end

function event_composition_storage_signature(prepared, state, workspace)
    return (
        Base.summarysize(prepared.scheduler.definitions),
        Base.summarysize(prepared.actions),
        Base.summarysize(prepared.paths),
        Base.summarysize(prepared.acquisitions),
        Base.summarysize(state.scheduler.cursors),
        Base.summarysize(state.acquisitions),
        Base.summarysize(state.path_sampled),
        Base.summarysize(state.product_sequences),
        Base.summarysize(state.product_ready_timestamps),
        Base.summarysize(workspace.scheduler.due_slots),
        Base.summarysize(workspace.due_paths),
    )
end

@testset "Rolling-shutter row-band event lifecycle" begin
    values = reshape(collect(1.0:24.0), 6, 4)
    map = event_test_intensity_map(values)
    exposure = PlantDuration(1_000_000_000)
    line = PlantDuration(200_000_000)
    definition = RollingShutterAcquisitionDefinition(exposure;
        readiness_delay=PlantDuration(17))

    detector = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2)))
    mtf_before = detector_mtf(detector, 0.15, -0.2)
    prepared = prepare_rolling_shutter_acquisition(detector, map,
        definition)
    state = RollingShutterAcquisitionState(prepared)

    tiny_line_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(1.0e-12;
            row_group_size=2)))
    tiny_line_error = event_test_error() do
        prepare_rolling_shutter_acquisition(tiny_line_detector, map,
            definition)
    end
    @test tiny_line_error isa DetectorAcquisitionError
    @test tiny_line_error.reason == :unrepresentable_line_duration

    @test rolling_band_count(prepared) == 3
    @test rolling_band_rows(prepared, 1) == 1:2
    @test rolling_band_rows(prepared, 2) == 3:4
    @test rolling_band_rows(prepared, 3) == 5:6
    begin_exposure!(prepared, state, PlantTimestamp(100))
    @test rolling_opened_band_count(state) == 1
    @test rolling_closed_band_count(state) == 0
    @test next_rolling_band_open_timestamp(prepared, state) ==
        PlantTimestamp(100) + line
    @test next_rolling_band_close_timestamp(prepared, state) ==
        PlantTimestamp(100) + exposure
    @test readout_complete_timestamp(state) ==
        PlantTimestamp(100) + exposure + 3line

    first_open = next_rolling_band_open_timestamp(prepared, state)
    accumulate_rolling_exposure_interval!(prepared, state,
        PlantTimestamp(100), first_open, Xoshiro(700))
    open_next_rolling_band!(prepared, state, first_open)
    @test rolling_opened_band_count(state) == 2

    # Restart on a fresh owner for the complete constant-rate oracle.
    oracle_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2)))
    oracle_prepared = prepare_rolling_shutter_acquisition(oracle_detector,
        map, definition)
    oracle_state = RollingShutterAcquisitionState(oracle_prepared)
    output = run_rolling_event_frame!(oracle_prepared, oracle_state,
        Xoshiro(701), PlantTimestamp(0))
    @test output == values
    @test rolling_opened_band_count(oracle_state) == 3
    @test rolling_closed_band_count(oracle_state) == 3
    @test detector_acquisition_status(oracle_state) ==
        DetectorAcquisitionReady
    @test detector_mtf(oracle_detector, 0.15, -0.2) == mtf_before

    reset_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.2;
            row_group_size=2, exposure_mode=GlobalResetExposure())))
    reset_prepared = prepare_rolling_shutter_acquisition(reset_detector,
        event_test_intensity_map(ones(6, 2)), definition)
    reset_state = RollingShutterAcquisitionState(reset_prepared)
    reset_output = run_rolling_event_frame!(reset_prepared, reset_state,
        Xoshiro(702), PlantTimestamp(0))
    @test reset_output == [
        1.0 1.0
        1.0 1.0
        1.2 1.2
        1.2 1.2
        1.4 1.4
        1.4 1.4
    ]
    @test rolling_band_open_timestamp(reset_prepared, reset_state, 3) ==
        PlantTimestamp(0)
    @test rolling_band_close_timestamp(reset_prepared, reset_state, 3) ==
        PlantTimestamp(1_400_000_000)

    if coverage_instrumented()
        @test_skip "rolling event allocation gate disabled under coverage instrumentation"
    else
        next_start = acquisition_readiness_timestamp(oracle_state) +
            PlantDuration(1)
        run_rolling_event_frame!(oracle_prepared, oracle_state,
            Xoshiro(703), next_start)
        following_start = acquisition_readiness_timestamp(oracle_state) +
            PlantDuration(1)
        allocation_rng = Xoshiro(704)
        @test @allocated(run_rolling_event_frame!(oracle_prepared,
            oracle_state, allocation_rng, following_start)) == 0
    end
end

@testset "Frame-transfer storage and overlapping acquisition" begin
    values = ones(3, 3)
    map = event_test_intensity_map(values)
    detector = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
        gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=0.1)))
    prepared = prepare_frame_transfer_acquisition(detector, map,
        FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(800_000_000)))
    state = FrameTransferAcquisitionState(prepared)
    rng = Xoshiro(710)
    mtf_before = detector_mtf(detector, 0.2, 0.1)

    tiny_transfer_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=1.0e-12)))
    tiny_transfer_error = event_test_error() do
        prepare_frame_transfer_acquisition(tiny_transfer_detector, map,
            FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000)))
    end
    @test tiny_transfer_error isa DetectorAcquisitionError
    @test tiny_transfer_error.reason == :unrepresentable_transfer_duration

    @test frame_transfer_storage_capacity(prepared) == 1
    @test !Base.mightalias(prepared.storage_frame, detector.state.frame)
    @test !Base.mightalias(prepared.storage_frame,
        detector.state.accum_buffer)
    @test !Base.mightalias(prepared.storage_frame, output_frame(detector))

    first_start = PlantTimestamp(0)
    first_close = PlantTimestamp(1_000_000_000)
    first_transfer = PlantTimestamp(1_100_000_000)
    first_readout = PlantTimestamp(1_900_000_000)
    begin_exposure!(prepared, state, first_start)
    accumulate_exposure_interval!(prepared, state, first_start,
        first_close, rng)
    close_exposure!(prepared, state, first_close)
    complete_frame_transfer!(prepared, state, first_transfer)
    @test frame_transfer_image_ready(state)
    @test frame_transfer_readout_pending(state)
    @test frame_transfer_storage_sequence(state) == 1

    regression = event_test_error() do
        begin_exposure!(prepared, state, first_transfer - PlantDuration(1))
    end
    @test regression isa DetectorAcquisitionError
    @test regression.reason == :time_regression

    slow_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition(
            transfer_time=0.1)))
    slow_prepared = prepare_frame_transfer_acquisition(slow_detector, map,
        FrameTransferAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(1_500_000_000)))
    slow_state = FrameTransferAcquisitionState(slow_prepared)
    begin_exposure!(slow_prepared, slow_state, first_start)
    accumulate_exposure_interval!(slow_prepared, slow_state, first_start,
        first_close, Xoshiro(711))
    close_exposure!(slow_prepared, slow_state, first_close)
    complete_frame_transfer!(slow_prepared, slow_state, first_transfer)
    capacity_error = event_test_error() do
        begin_exposure!(slow_prepared, slow_state, first_transfer)
    end
    @test capacity_error isa DetectorAcquisitionError
    @test capacity_error.reason == :storage_capacity
    @test detector_acquisition_sequence(slow_state) == 1
    @test frame_transfer_storage_sequence(slow_state) == 1

    fill!(values, 2.0)
    second_start = first_transfer
    second_close = PlantTimestamp(2_100_000_000)
    second_transfer = PlantTimestamp(2_200_000_000)
    begin_exposure!(prepared, state, second_start)
    accumulate_exposure_interval!(prepared, state, second_start,
        first_readout, rng)
    first_output = copy(complete_readout!(prepared, state, first_readout,
        rng))
    @test first_output == ones(3, 3)
    @test frame_transfer_image_sequence(state) == 2
    @test frame_transfer_product_sequence(state) == 1

    accumulate_exposure_interval!(prepared, state, first_readout,
        second_close, rng)
    close_exposure!(prepared, state, second_close)
    complete_frame_transfer!(prepared, state, second_transfer)
    second_output = copy(complete_readout!(prepared, state,
        PlantTimestamp(3_000_000_000), rng))
    @test second_output == fill(2.0, 3, 3)
    @test frame_transfer_product_sequence(state) == 2
    @test frame_transfer_storage_empty(state)
    @test detector_mtf(detector, 0.2, 0.1) == mtf_before

    warm_start = PlantTimestamp(3_100_000_000)
    run_frame_transfer_event!(prepared, state, rng, warm_start)
    allocation_start = PlantTimestamp(5_100_000_000)
    if coverage_instrumented()
        @test_skip "frame-transfer allocation gate disabled under coverage instrumentation"
    else
        @test @allocated(run_frame_transfer_event!(prepared, state, rng,
            allocation_start)) == 0
    end
end

@testset "Multi-rate plant event composition" begin
    plant, prepared, state, workspace = event_composition_fixture()
    @test plant_event_path_count(prepared) == 3
    @test plant_event_acquisition_count(prepared) == 5
    @test plant_event_generator_count(prepared) == 1 + 3 + 5 * 5
    @test next_plant_event_timestamp(prepared, state, workspace) ==
        PlantTimestamp(0)

    unbound_error = event_test_error() do
        event_composition_fixture(unbound_trigger_consumer=true)
    end
    @test unbound_error isa PlantScheduleError
    @test unbound_error.component == :plant_event_loop
    @test unbound_error.reason == :unbound_trigger_consumer

    _, atomic_prepared, atomic_state, atomic_workspace =
        event_composition_fixture()
    run_plant_events_until!(atomic_prepared, atomic_state,
        atomic_workspace, PlantTimestamp(299_999_999))
    ccd_slot = findfirst(acquisition ->
        acquisition.id == AcquisitionID(:science_ccd),
        atomic_prepared.acquisitions)
    rolling_slot = findfirst(acquisition ->
        acquisition.id == AcquisitionID(:science_rolling),
        atomic_prepared.acquisitions)
    ccd_state = atomic_state.acquisitions[ccd_slot]
    ccd_progress = integrated_through_timestamp(ccd_state)
    rolling_lifecycle = atomic_prepared.acquisitions[
        rolling_slot].lifecycle
    rolling_lifecycle.detector.state.readout_ready = true
    preflight_error = event_test_error() do
        step_plant_events!(atomic_prepared, atomic_state, atomic_workspace)
    end
    @test preflight_error isa DetectorAcquisitionError
    @test preflight_error.reason == :detector_state_changed
    @test integrated_through_timestamp(ccd_state) == ccd_progress

    horizon = PlantTimestamp(1_500_000_000)
    timestamp_count = run_plant_events_until!(prepared, state, workspace,
        horizon)
    @test timestamp_count > 0
    @test acquisition_product_sequence(prepared, state,
        :science_cmos) == 3
    @test acquisition_product_sequence(prepared, state,
        :science_ccd) == 2
    @test acquisition_product_sequence(prepared, state,
        :science_rolling) == 3
    @test acquisition_product_sequence(prepared, state,
        :ngs_saphira) == 4
    @test acquisition_product_sequence(prepared, state,
        :lgs_emccd) == 3
    @test acquisition_product_ready_timestamp(prepared, state,
        :science_cmos) == PlantTimestamp(1_220_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :science_ccd) == PlantTimestamp(950_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :ngs_saphira) == PlantTimestamp(1_460_000_000)
    @test acquisition_product_ready_timestamp(prepared, state,
        :lgs_emccd) == PlantTimestamp(1_245_000_000)

    @test prepared_path(plant, :science).execution.executions[] == 16
    @test prepared_path(plant, :ngs).execution.executions[] == 12
    @test prepared_path(plant, :lgs).execution.executions[] == 8
    science_cmos = prepared_acquisition(plant, :science_cmos)
    science_ccd = prepared_acquisition(plant, :science_ccd)
    @test science_cmos.path_result === science_ccd.path_result
    @test all(isfinite, acquisition_observation(science_cmos))
    @test all(isfinite, acquisition_observation(science_ccd))

    reordered_plant, reordered, reordered_state, reordered_workspace =
        event_composition_fixture(reverse_order=true)
    @test run_plant_events_until!(reordered, reordered_state,
        reordered_workspace, horizon) == timestamp_count
    for id in (:science_cmos, :science_ccd, :science_rolling,
            :ngs_saphira, :lgs_emccd)
        @test acquisition_product_sequence(reordered, reordered_state, id) ==
            acquisition_product_sequence(prepared, state, id)
        @test acquisition_product_ready_timestamp(reordered,
            reordered_state, id) == acquisition_product_ready_timestamp(
            prepared, state, id)
        @test acquisition_observation(prepared_acquisition(reordered_plant,
            id)) == acquisition_observation(prepared_acquisition(plant, id))
    end


    storage_before = event_composition_storage_signature(prepared, state,
        workspace)
    run_event_composition_window!(prepared, state, workspace,
        PlantTimestamp(2_000_000_000))
    if coverage_instrumented()
        @test_skip "plant event-loop allocation gate disabled under coverage instrumentation"
    else
        processed = Ref(0)
        allocated = @allocated processed[] =
            run_event_composition_window!(prepared, state, workspace,
                PlantTimestamp(2_500_000_000))
        # The heterogeneous orchestration barrier deliberately avoids
        # specializing on the complete instrument tuple. Keep its bounded
        # dynamic-dispatch cost below 2 KiB per processed timestamp while the
        # detector and optical kernels remain allocation-free after warmup.
        @test allocated <= 2_048 * processed[]
    end
    long_run_timestamps = run_event_composition_window!(prepared, state,
        workspace, PlantTimestamp(25_000_000_000))
    @test 0 < long_run_timestamps < 2_000
    @test event_composition_storage_signature(prepared, state, workspace) ==
        storage_before
    @test plant_event_generator_count(prepared) == 29

    faulted_plant, faulted, faulted_state, faulted_workspace =
        event_composition_fixture(faulted_trigger_fanout=true)
    faulted_timestamp_count = run_plant_events_until!(faulted,
        faulted_state, faulted_workspace, PlantTimestamp(3_000_000_000))
    @test faulted_timestamp_count > 0
    @test acquisition_product_sequence(faulted, faulted_state,
        :science_cmos) == 5
    @test acquisition_product_sequence(faulted, faulted_state,
        :ngs_saphira) == 5
    @test acquisition_product_ready_timestamp(faulted, faulted_state,
        :science_cmos) == PlantTimestamp(2_690_000_000)
    @test acquisition_product_ready_timestamp(faulted, faulted_state,
        :ngs_saphira) == PlantTimestamp(2_680_000_000)
    @test all(isfinite, acquisition_observation(prepared_acquisition(
        faulted_plant, :science_cmos)))
    @test all(isfinite, acquisition_observation(prepared_acquisition(
        faulted_plant, :ngs_saphira)))
end
