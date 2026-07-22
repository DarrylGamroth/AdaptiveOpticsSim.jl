function plant_detector_intensity_map(
    values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.55e-6)))
    return IntensityMap(metadata, values)
end

function caught_detector_acquisition_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function run_prepared_detector_event_cycle!(prepared, state, rng,
    start::PlantTimestamp)
    close = start + prepared.definition.exposure_duration
    complete = close + prepared.definition.readout_duration
    ready = complete + prepared.definition.readiness_delay
    begin_exposure!(prepared, state, start)
    accumulate_exposure_interval!(prepared, state, start, close, rng)
    close_exposure!(prepared, state, close)
    complete_readout!(prepared, state, complete, rng)
    mark_acquisition_ready!(prepared, state, ready)
    return nothing
end

function run_prepared_ramp_event_cycle!(prepared, state, rng,
    start::PlantTimestamp)
    close = start + prepared.definition.exposure_duration
    begin_exposure!(prepared, state, start)
    take_nondestructive_read!(prepared, state, start, rng)
    first_stop = start + nondestructive_read_offset(prepared, 2)
    accumulate_exposure_interval!(prepared, state, start, first_stop, rng)
    take_nondestructive_read!(prepared, state, first_stop, rng)
    accumulate_exposure_interval!(prepared, state, first_stop, close, rng)
    take_nondestructive_read!(prepared, state, close, rng)
    close_exposure!(prepared, state, close)
    complete = close + prepared.definition.readout_duration
    complete_readout!(prepared, state, complete, rng)
    mark_acquisition_ready!(prepared, state,
        complete + prepared.definition.readiness_delay)
    return nothing
end

@testset "Prepared detector event definitions" begin
    zero_error = caught_detector_acquisition_error() do
        GlobalShutterAcquisitionDefinition(PlantDuration(0))
    end
    @test zero_error isa DetectorAcquisitionError
    @test zero_error.reason == :invalid_definition

    values = ones(4, 4)
    map = plant_detector_intensity_map(values)
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(20),
        readiness_delay=PlantDuration(30))
    detector = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0)
    plan = prepare_detector_acquisition(detector, map)
    prepared = prepare_global_shutter_acquisition(detector, map, plan,
        definition)
    state = GlobalShutterAcquisitionState(prepared)

    @test detector_acquisition_status(state) == DetectorAcquisitionReady
    @test detector_acquisition_sequence(state) == 0
    @test nondestructive_read_count(prepared) == 0
    @test_throws DetectorAcquisitionError nondestructive_read_offset(prepared, 1)
    @test isnothing(next_nondestructive_read_timestamp(prepared, state))

    mismatch = Detector(integration_time=0.5, noise=NoiseNone())
    mismatch_error = caught_detector_acquisition_error() do
        prepare_global_shutter_acquisition(mismatch, map, definition)
    end
    @test mismatch_error isa DetectorAcquisitionError
    @test mismatch_error.reason == :exposure_duration

    rolling = Detector(integration_time=1.0, noise=NoiseNone(),
        sensor=CMOSSensor(timing_model=RollingShutter(1e-6)))
    @test size(rolling.state.frame) == (1, 1)
    rolling_error = caught_detector_acquisition_error() do
        prepare_global_shutter_acquisition(rolling, map, definition)
    end
    @test rolling_error isa DetectorAcquisitionError
    @test rolling_error.reason == :unsupported_timing
    @test size(rolling.state.frame) == (1, 1)

    quantized_sensor = HgCdTeAvalancheArraySensor(
        sampling_mode=UpTheRampSampling(4), read_time=0.0)
    quantized_detector = Detector(integration_time=1e-8,
        noise=NoiseNone(), sensor=quantized_sensor)
    quantized_map = plant_detector_intensity_map(ones(2, 2))
    quantized_definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(10))
    quantized = prepare_global_shutter_acquisition(quantized_detector,
        quantized_map, quantized_definition)
    @test nondestructive_read_count(quantized) == 4
    @test Tuple(plant_nanoseconds(nondestructive_read_offset(quantized, i))
        for i in 1:4) == (0, 3, 6, 10)
    @test detector_ramp_times(quantized_detector) == [0.0, 3e-9, 6e-9, 1e-8]
    capture!(quantized_detector, ones(2, 2); rng=Xoshiro(99))
    @test detector_ramp_times(quantized_detector)[2] > 3e-9
    quantized_state = GlobalShutterAcquisitionState(quantized)
    begin_exposure!(quantized, quantized_state, PlantTimestamp(0))
    @test detector_ramp_times(quantized_detector) ==
        [0.0, 3e-9, 6e-9, 1e-8]

    spacing_detector = Detector(integration_time=1e-9, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=UpTheRampSampling(3), read_time=0.0))
    spacing_error = caught_detector_acquisition_error() do
        prepare_global_shutter_acquisition(spacing_detector,
            plant_detector_intensity_map(ones(2, 2)),
            GlobalShutterAcquisitionDefinition(PlantDuration(1)))
    end
    @test spacing_error isa DetectorAcquisitionError
    @test spacing_error.reason == :unrepresentable_read_schedule

    read_time_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=UpTheRampSampling(3), read_time=0.1))
    readout_error = caught_detector_acquisition_error() do
        prepare_global_shutter_acquisition(read_time_detector,
            plant_detector_intensity_map(ones(2, 2)),
            GlobalShutterAcquisitionDefinition(PlantDuration(1_000_000_000);
                readout_duration=PlantDuration(99_999_999)))
    end
    @test readout_error isa DetectorAcquisitionError
    @test readout_error.reason == :readout_duration
    exact_readout = prepare_global_shutter_acquisition(read_time_detector,
        plant_detector_intensity_map(ones(2, 2)),
        GlobalShutterAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(100_000_000)))
    @test nondestructive_read_count(exact_readout) == 3
end

@testset "Global-shutter event lifecycle and frame oracle" begin
    kernel = [0.0 0.1 0.0; 0.2 0.4 0.1; 0.0 0.2 0.0]
    coupling = [0.0 0.02 0.0; 0.01 0.94 0.01; 0.0 0.02 0.0]
    values = reshape(collect(1.0:25.0), 5, 5)
    map = plant_detector_intensity_map(values)
    detector_keywords = (
        integration_time=1.0,
        qe=0.6,
        noise=NoiseNone(),
        response_model=SampledFrameResponse(kernel),
        charge_coupling_model=InterpixelCapacitance(coupling),
    )
    oracle_detector = Detector(; detector_keywords...)
    event_detector = Detector(; detector_keywords...)
    oracle_output = copy(capture!(oracle_detector, map.values;
        rng=Xoshiro(400)))
    mtf_before = detector_mtf(event_detector, 0.2, -0.15)

    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(10),
        readiness_delay=PlantDuration(20))
    prepared = prepare_global_shutter_acquisition(event_detector, map,
        definition)
    state = GlobalShutterAcquisitionState(prepared)
    rng = Xoshiro(400)
    start = PlantTimestamp(1_000)
    middle = PlantTimestamp(400_001_000)
    close = PlantTimestamp(1_000_001_000)

    @test @inferred(begin_exposure!(prepared, state, start)) === nothing
    @test detector_acquisition_status(state) == DetectorExposureActive
    @test detector_acquisition_sequence(state) == 1
    @test exposure_start_timestamp(state) == start
    @test exposure_close_timestamp(state) == close
    @test readout_complete_timestamp(state) == close + PlantDuration(10)
    @test acquisition_readiness_timestamp(state) == close + PlantDuration(30)
    @test !readout_ready(event_detector)

    retrigger_rng = Xoshiro(77)
    reference_rng = Xoshiro(77)
    accumulator_before = copy(event_detector.state.accum_buffer)
    retrigger_error = caught_detector_acquisition_error() do
        begin_exposure!(prepared, state, start + PlantDuration(1))
    end
    @test retrigger_error isa DetectorAcquisitionError
    @test retrigger_error.reason == :retrigger
    @test event_detector.state.accum_buffer == accumulator_before
    @test rand(retrigger_rng) == rand(reference_rng)

    after_close_error = caught_detector_acquisition_error() do
        accumulate_exposure_interval!(prepared, state, start,
            close + PlantDuration(1), retrigger_rng)
    end
    @test after_close_error isa DetectorAcquisitionError
    @test after_close_error.reason == :interval_after_close
    @test integrated_through_timestamp(state) == start
    @test event_detector.state.accum_buffer == accumulator_before

    early_close_error = caught_detector_acquisition_error() do
        close_exposure!(prepared, state, close)
    end
    @test early_close_error isa DetectorAcquisitionError
    @test early_close_error.reason == :incomplete_integration
    @test detector_acquisition_status(state) == DetectorExposureActive

    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    @test integrated_through_timestamp(state) == close
    accumulator_at_close = copy(event_detector.state.accum_buffer)

    boundary_error = caught_detector_acquisition_error() do
        accumulate_exposure_interval!(prepared, state, close,
            close + PlantDuration(1), rng)
    end
    @test boundary_error isa DetectorAcquisitionError
    @test boundary_error.reason == :interval_after_close
    @test event_detector.state.accum_buffer == accumulator_at_close

    close_exposure!(prepared, state, close)
    @test detector_acquisition_status(state) == DetectorReadoutPending
    wrong_readout_error = caught_detector_acquisition_error() do
        complete_readout!(prepared, state, close, rng)
    end
    @test wrong_readout_error isa DetectorAcquisitionError
    @test wrong_readout_error.reason == :readout_not_due
    @test detector_acquisition_status(state) == DetectorReadoutPending
    @test event_detector.state.accum_buffer == accumulator_at_close

    output = complete_readout!(prepared, state, close + PlantDuration(10), rng)
    @test output ≈ oracle_output atol=1e-12 rtol=1e-12
    @test detector_acquisition_status(state) == DetectorReadoutComplete
    @test !readout_ready(event_detector)
    @test detector_mtf(event_detector, 0.2, -0.15) == mtf_before
    @test all(isfinite, output)
    @test sum(output) <= sum(values) * 0.6

    wrong_ready_error = caught_detector_acquisition_error() do
        mark_acquisition_ready!(prepared, state, close + PlantDuration(29))
    end
    @test wrong_ready_error isa DetectorAcquisitionError
    @test wrong_ready_error.reason == :readiness_not_due
    @test !readout_ready(event_detector)
    @test mark_acquisition_ready!(prepared, state,
        close + PlantDuration(30)) === output_frame(event_detector)
    @test readout_ready(event_detector)
    @test detector_acquisition_status(state) == DetectorAcquisitionReady

    sequence_before_regression = detector_acquisition_sequence(state)
    time_regression = caught_detector_acquisition_error() do
        begin_exposure!(prepared, state, start)
    end
    @test time_regression isa DetectorAcquisitionError
    @test time_regression.reason == :time_regression
    @test detector_acquisition_sequence(state) == sequence_before_regression
    @test detector_acquisition_status(state) == DetectorAcquisitionReady
    @test readout_ready(event_detector)

    other_detector = Detector(integration_time=1.0, noise=NoiseNone())
    other_map = plant_detector_intensity_map(ones(2, 2))
    other_prepared = prepare_global_shutter_acquisition(other_detector,
        other_map, GlobalShutterAcquisitionDefinition(
            PlantDuration(1_000_000_000)))
    foreign_error = caught_detector_acquisition_error() do
        begin_exposure!(other_prepared, state, PlantTimestamp(0))
    end
    @test foreign_error isa DetectorAcquisitionError
    @test foreign_error.reason == :foreign_state
end

@testset "Scheduled evolving-charge up-the-ramp reads" begin
    values = fill(1.0, 4, 4)
    map = plant_detector_intensity_map(values)
    detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=UpTheRampSampling(3), read_time=0.0),
        readout_window=FrameWindow(2:3, 2:3))
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(100),
        readiness_delay=PlantDuration(200))
    prepared = prepare_global_shutter_acquisition(detector, map, definition)
    state = GlobalShutterAcquisitionState(prepared)
    rng = Xoshiro(501)
    start = PlantTimestamp(5_000)
    middle = start + PlantDuration(500_000_000)
    close = start + PlantDuration(1_000_000_000)

    @test isnothing(next_nondestructive_read_timestamp(prepared, state))
    begin_exposure!(prepared, state, start)
    @test next_nondestructive_read_timestamp(prepared, state) == start
    missed_initial = caught_detector_acquisition_error() do
        accumulate_exposure_interval!(prepared, state, start,
            start + PlantDuration(1), rng)
    end
    @test missed_initial isa DetectorAcquisitionError
    @test missed_initial.reason == :missed_nondestructive_read
    @test integrated_through_timestamp(state) == start
    @test all(iszero, detector.state.accum_buffer)

    invalid_read_rng = Xoshiro(502)
    reference_read_rng = Xoshiro(502)
    ramp_cube_before = copy(detector_ramp_cube(detector))
    wrong_read_time = caught_detector_acquisition_error() do
        take_nondestructive_read!(prepared, state,
            start + PlantDuration(1), invalid_read_rng)
    end
    @test wrong_read_time isa DetectorAcquisitionError
    @test wrong_read_time.reason == :nondestructive_read_not_due
    @test detector_ramp_cube(detector) == ramp_cube_before
    @test rand(invalid_read_rng) == rand(reference_read_rng)

    take_nondestructive_read!(prepared, state, start, rng)
    @test next_nondestructive_read_timestamp(prepared, state) == middle
    @test all(iszero, detector_ramp_cube(detector)[:, :, 1])
    incomplete_read = caught_detector_acquisition_error() do
        take_nondestructive_read!(prepared, state, middle, rng)
    end
    @test incomplete_read isa DetectorAcquisitionError
    @test incomplete_read.reason == :incomplete_integration
    @test next_nondestructive_read_timestamp(prepared, state) == middle
    crossed_read = caught_detector_acquisition_error() do
        accumulate_exposure_interval!(prepared, state, start,
            middle + PlantDuration(1), rng)
    end
    @test crossed_read isa DetectorAcquisitionError
    @test crossed_read.reason == :missed_nondestructive_read

    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    take_nondestructive_read!(prepared, state, middle, rng)
    @test all(==(0.5), detector_ramp_cube(detector)[:, :, 2])
    @test all(==(0.5), detector.state.accum_buffer)

    fill!(values, 3.0)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    @test all(==(2.0), detector.state.accum_buffer)
    fill!(values, 99.0)
    close_exposure!(prepared, state, close)
    @test detector_acquisition_status(state) == DetectorReadoutPending

    missing_final = caught_detector_acquisition_error() do
        complete_readout!(prepared, state, close + PlantDuration(100), rng)
    end
    @test missing_final isa DetectorAcquisitionError
    @test missing_final.reason == :missing_nondestructive_read
    @test all(==(2.0), detector.state.accum_buffer)

    take_nondestructive_read!(prepared, state, close, rng)
    @test isnothing(next_nondestructive_read_timestamp(prepared, state))
    @test all(==(2.0), detector_ramp_cube(detector)[:, :, 3])
    output = complete_readout!(prepared, state,
        close + PlantDuration(100), rng)
    @test all(==(2.0), detector_ramp_slope(detector))
    @test all(isapprox.(
        detector_ramp_intercept(detector), -1 / 6; atol=1e-15, rtol=0))
    @test all(==(2.0), output)
    @test detector_ramp_times(detector) == [0.0, 0.5, 1.0]
    @test !readout_ready(detector)
    mark_acquisition_ready!(prepared, state, close + PlantDuration(300))
    @test readout_ready(detector)
end

@testset "Ramp snapshots retain prior avalanche realization" begin
    values = fill(4.0, 2, 2)
    map = plant_detector_intensity_map(values)
    detector = Detector(integration_time=1.0, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=1.0,
            excess_noise_factor=1.4,
            sampling_mode=UpTheRampSampling(3),
            read_time=0.0,
        ))
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    prepared = prepare_global_shutter_acquisition(detector, map, definition)
    state = GlobalShutterAcquisitionState(prepared)
    rng = Xoshiro(503)
    start = PlantTimestamp(0)
    middle = PlantTimestamp(500_000_000)
    close = PlantTimestamp(1_000_000_000)

    begin_exposure!(prepared, state, start)
    take_nondestructive_read!(prepared, state, start, rng)
    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    take_nondestructive_read!(prepared, state, middle, rng)
    middle_read = copy(@view detector_ramp_cube(detector)[:, :, 2])
    fill!(values, 0.0)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    take_nondestructive_read!(prepared, state, close, rng)

    @test detector_ramp_cube(detector)[:, :, 3] == middle_read
    close_exposure!(prepared, state, close)
    complete_readout!(prepared, state, close, rng)
    mark_acquisition_ready!(prepared, state, close)
end

@testset "Detector event inference and warmed allocation" begin
    values = fill(2.0, 8, 8)
    map = plant_detector_intensity_map(values)
    detector = Detector(integration_time=1.0, qe=0.5,
        noise=NoiseNone())
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    prepared = prepare_global_shutter_acquisition(detector, map, definition)
    state = GlobalShutterAcquisitionState(prepared)
    rng = Xoshiro(601)
    run_prepared_detector_event_cycle!(prepared, state, rng,
        PlantTimestamp(0))
    @test @inferred(run_prepared_detector_event_cycle!(prepared, state, rng,
        PlantTimestamp(1_000_000_000))) === nothing
    if !coverage_instrumented()
        @test @allocated(run_prepared_detector_event_cycle!(prepared, state,
            rng, PlantTimestamp(2_000_000_000))) == 0
    end

    ramp_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=UpTheRampSampling(3), read_time=0.0))
    ramp_map = plant_detector_intensity_map(fill(1.0, 8, 8))
    ramp_prepared = prepare_global_shutter_acquisition(ramp_detector,
        ramp_map, definition)
    ramp_state = GlobalShutterAcquisitionState(ramp_prepared)
    ramp_rng = Xoshiro(602)
    run_prepared_ramp_event_cycle!(ramp_prepared, ramp_state, ramp_rng,
        PlantTimestamp(0))
    @test @inferred(run_prepared_ramp_event_cycle!(ramp_prepared, ramp_state,
        ramp_rng, PlantTimestamp(1_000_000_000))) === nothing
    if !coverage_instrumented()
        @test @allocated(run_prepared_ramp_event_cycle!(ramp_prepared,
            ramp_state, ramp_rng, PlantTimestamp(2_000_000_000))) == 0
    end
end
