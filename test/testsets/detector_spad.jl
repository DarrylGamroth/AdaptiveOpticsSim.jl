function test_spad_poisson_moments(samples, expected_mean;
    sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_mean / sample_count)
    variance_limit = sigma_limit * sqrt(
        (expected_mean + 2 * expected_mean^2) / (sample_count - 1))
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_mean) <= variance_limit
end

struct IncompleteSPADContractDetector <: AbstractCountingDetector end

@testset "SPAD counting dispatch contracts" begin
    incomplete = IncompleteSPADContractDetector()
    @test_throws InvalidConfiguration counting_mean_response_model(incomplete)
    @test counting_detection_efficiency(incomplete, Float32) === 1.0f0

    null_response = NullCountingMeanResponse()
    afterpulse = FirstOrderAfterpulseMeanResponse(1 // 4)
    redistribution = NearestNeighborCountRedistribution(2 // 5)
    composite = CompositeCountingMeanResponse(afterpulse, redistribution)

    @test afterpulse.mean_afterpulses_per_detection == 0.25
    @test redistribution.redistribution_fraction == 0.4
    @test_throws InvalidConfiguration CompositeCountingMeanResponse(
        (null_response, 1))

    @test !_supports_first_order_afterpulse_mean_response(null_response)
    @test _supports_first_order_afterpulse_mean_response(afterpulse)
    @test _supports_first_order_afterpulse_mean_response(composite)
    @test !_supports_nearest_neighbor_count_redistribution(null_response)
    @test _supports_nearest_neighbor_count_redistribution(redistribution)
    @test _supports_nearest_neighbor_count_redistribution(composite)

    @test _afterpulse_stage_count(null_response) == 0
    @test _afterpulse_stage_count(afterpulse) == 1
    @test _afterpulse_stage_count(composite) == 1
    @test _redistribution_stage_count(null_response) == 0
    @test _redistribution_stage_count(redistribution) == 1
    @test _redistribution_stage_count(composite) == 1

    @test counting_mean_response_symbol(null_response) == :none
    @test counting_mean_response_symbol(afterpulse) ==
        :first_order_afterpulse_mean_response
    @test counting_mean_response_symbol(redistribution) ==
        :nearest_neighbor_count_redistribution
    @test counting_mean_response_symbol(composite) == :composite

    detector = SPADArrayDetector((1, 1))
    @test counting_layout(detector) == :pixel_counts
end

@testset "SPAD ownership and deterministic radiometry" begin
    sensor = SPADArraySensor(active_area_detection_efficiency=0.5,
        dark_count_rate=0.0, fill_factor=0.8)
    detector = SPADArrayDetector((2, 8); integration_time=1.0,
        noise=NoiseNone(), sensor=sensor)
    output = capture!(detector, fill(10.0, 2, 8); rng=Xoshiro(9101))
    @test output == fill(4.0, 2, 8)
    @test output_frame(detector) === output
    @test channel_output(detector) === output
    @test detector.params.dimensions == (2, 8)
    @test isconcretetype(typeof(sensor))
    @test !ismutabletype(typeof(sensor))
    @test !ismutabletype(typeof(detector.params))
    @test detector.state isa Detectors.CountingDetectorState
    @test detector.workspace isa Detectors.CountingDetectorWorkspace
    @test detector.products isa Detectors.CountingDetectorProducts
    @test isconcretetype(typeof(detector.state))
    @test isconcretetype(typeof(detector.workspace))
    @test isconcretetype(typeof(detector.products))
    @test Detectors.counting_array(detector) === detector.products.counts
    @test !Base.mightalias(detector.workspace.noise_buffer,
        detector.products.counts)

    metadata = detector_export_metadata(detector)
    @test metadata isa CountingDetectorExportMetadata
    @test metadata.sensor == :spad_array
    @test metadata.detection_efficiency == 0.5
    @test metadata.fill_factor == 0.8
    @test metadata.gain == 1.0
    @test metadata.readout.layout == :pixel_counts
    @test metadata.readout.output_size == (2, 8)
    @test metadata.readout.n_channels == 16
    @test !supports_channel_gain_map(detector)
    @test supports_counting_noise(detector)
    @test !supports_dead_time(detector)

    gated = SPADArrayDetector((2, 8); integration_time=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25), sensor=sensor)
    @test counting_exposure_time(gated) == 0.5
    @test capture!(gated, fill(10.0, 2, 8), Xoshiro(9102)) ==
        fill(2.0, 2, 8)
    @test supports_counting_gating(gated)

    gated_dark = SPADArrayDetector((1, 1); integration_time=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25),
        sensor=SPADArraySensor(active_area_detection_efficiency=0.0,
            dark_count_rate=4.0))
    @test capture!(gated_dark, zeros(1, 1), Xoshiro(9103)) ==
        fill(2.0, 1, 1)

    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    temperature = 250.0
    thermal = SPADArrayDetector((1, 1); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=0.0,
            dark_count_rate=10.0),
        thermal_model=FixedTemperature(temperature_K=temperature,
            dark_count_law=arrhenius))
    expected_dark_rate = 10.0 * exp(6000.0 * (inv(300.0) - inv(temperature)))
    @test effective_dark_count_rate(thermal) ≈ expected_dark_rate
    @test capture!(thermal, zeros(1, 1), Xoshiro(9104)) ≈
        fill(expected_dark_rate, 1, 1)
    @test detector_export_metadata(thermal).dark_count_law == :arrhenius
end

@testset "SPAD dead-time mean laws" begin
    dead_time = 0.1
    dimensionless_rates = (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)
    for (model, expected_response) in (
        (NoDeadTime(), x -> x / dead_time),
        (NonParalyzableDeadTime(dead_time), x -> (x / dead_time) / (1 + x)),
        (ParalyzableDeadTime(dead_time), x -> (x / dead_time) * exp(-x)),
    )
        detector = SPADArrayDetector((1, 1); integration_time=1.0,
            noise=NoiseNone(), sensor=SPADArraySensor(
                active_area_detection_efficiency=1.0,
                dead_time_model=model))
        for x in dimensionless_rates
            incident_rate = x / dead_time
            observed = only(capture!(detector, fill(incident_rate, 1, 1),
                Xoshiro(9110)))
            @test observed ≈ expected_response(x) rtol=1e-13 atol=1e-13
        end
    end

    nonparalyzable = SPADArrayDetector((1, 1); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            dead_time_model=NonParalyzableDeadTime(dead_time)))
    @test only(capture!(nonparalyzable, fill(1e8, 1, 1), Xoshiro(9111))) ≈
        inv(dead_time) rtol=1e-6

    paralyzable = SPADArrayDetector((1, 1); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            dead_time_model=ParalyzableDeadTime(dead_time)))
    maximum_response = only(capture!(paralyzable,
        fill(inv(dead_time), 1, 1), Xoshiro(9112)))
    @test maximum_response ≈ exp(-1) / dead_time
    @test only(capture!(paralyzable, fill(100 / dead_time, 1, 1),
        Xoshiro(9113))) < 1e-38
end

@testset "SPAD deterministic mean-response models" begin
    afterpulse = SPADArrayDetector((2, 8); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=FirstOrderAfterpulseMeanResponse(1.5)))
    @test capture!(afterpulse, fill(4.0, 2, 8), Xoshiro(9120)) ==
        fill(10.0, 2, 8)
    @test supports_first_order_afterpulse_mean_response(afterpulse)
    afterpulse_metadata = detector_export_metadata(afterpulse)
    @test afterpulse_metadata.mean_response_model ==
        :first_order_afterpulse_mean_response
    @test afterpulse_metadata.mean_afterpulses_per_detection == 1.5

    redistribution = SPADArrayDetector((3, 3); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=NearestNeighborCountRedistribution(0.4)))
    center_input = zeros(3, 3)
    center_input[2, 2] = 10.0
    center_output = copy(capture!(redistribution, center_input,
        Xoshiro(9121)))
    @test center_output == [0.0 1.0 0.0; 1.0 6.0 1.0; 0.0 1.0 0.0]
    @test sum(center_output) == sum(center_input)

    corner_input = zeros(3, 3)
    corner_input[1, 1] = 10.0
    corner_output = copy(capture!(redistribution, corner_input,
        Xoshiro(9122)))
    @test corner_output[1, 1] == 6.0
    @test corner_output[2, 1] == 2.0
    @test corner_output[1, 2] == 2.0
    @test sum(corner_output) == sum(corner_input)
    @test supports_nearest_neighbor_count_redistribution(redistribution)
    redistribution_metadata = detector_export_metadata(redistribution)
    @test redistribution_metadata.mean_response_model ==
        :nearest_neighbor_count_redistribution
    @test redistribution_metadata.redistribution_fraction == 0.4

    composite = SPADArrayDetector((3, 3); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=CompositeCountingMeanResponse(
                FirstOrderAfterpulseMeanResponse(0.25),
                NearestNeighborCountRedistribution(0.4))))
    composite_output = capture!(composite, center_input, Xoshiro(9123))
    @test composite_output == 1.25 .* center_output
    @test sum(composite_output) == 1.25 * sum(center_input)
    composite_metadata = detector_export_metadata(composite)
    @test composite_metadata.mean_response_model == :composite
    @test composite_metadata.mean_afterpulses_per_detection == 0.25
    @test composite_metadata.redistribution_fraction == 0.4

    duplicate_afterpulse = CompositeCountingMeanResponse(
        FirstOrderAfterpulseMeanResponse(0.2),
        CompositeCountingMeanResponse(
            NearestNeighborCountRedistribution(0.1),
            FirstOrderAfterpulseMeanResponse(0.3)))
    @test_throws InvalidConfiguration SPADArrayDetector((3, 3);
        sensor=SPADArraySensor(mean_response_model=duplicate_afterpulse))
    duplicate_redistribution = CompositeCountingMeanResponse(
        NearestNeighborCountRedistribution(0.1),
        NearestNeighborCountRedistribution(0.2))
    @test_throws InvalidConfiguration SPADArrayDetector((3, 3);
        sensor=SPADArraySensor(mean_response_model=duplicate_redistribution))

    singleton = SPADArrayDetector((1, 1); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=NearestNeighborCountRedistribution(0.4)))
    @test capture!(singleton, fill(10.0, 1, 1), Xoshiro(9124)) ==
        fill(10.0, 1, 1)
end

@testset "SPAD ordered deterministic pipeline" begin
    detector = SPADArrayDetector((1, 1); integration_time=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25),
        sensor=SPADArraySensor(
            active_area_detection_efficiency=0.5,
            dark_count_rate=1.0,
            dead_time_model=NonParalyzableDeadTime(0.1),
            mean_response_model=FirstOrderAfterpulseMeanResponse(0.2)))
    @test only(capture!(detector, fill(8.0, 1, 1), Xoshiro(9125))) ≈
        2.0 rtol=1e-14 atol=1e-14
end

@testset "SPAD Poisson surrogate moments" begin
    dimensions = (128, 128)
    cases = (
        (SPADArrayDetector(dimensions; noise=NoisePhoton(),
            sensor=SPADArraySensor(active_area_detection_efficiency=1.0)),
            fill(20.0, dimensions), 20.0),
        (SPADArrayDetector(dimensions; integration_time=2.0,
            noise=NoisePhoton(), sensor=SPADArraySensor(
                active_area_detection_efficiency=0.0,
                dark_count_rate=6.0)), zeros(dimensions), 12.0),
        (SPADArrayDetector(dimensions; noise=NoisePhoton(),
            sensor=SPADArraySensor(active_area_detection_efficiency=0.5,
                dark_count_rate=3.0)), fill(10.0, dimensions), 8.0),
        (SPADArrayDetector(dimensions; noise=NoisePhoton(),
            sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
                dead_time_model=NonParalyzableDeadTime(0.05))),
            fill(100.0, dimensions), 100 / 6),
    )
    for (index, (detector, input, expected_mean)) in enumerate(cases)
        samples = vec(copy(capture!(detector, input,
            Xoshiro(9130 + index))))
        test_spad_poisson_moments(samples, expected_mean)
    end
end

@testset "SPAD fixed topology, failures, replay, and hot path" begin
    detector = SPADArrayDetector((2, 3); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0),
        output_type=UInt16)
    input = fill(2.6, 2, 3)
    @test capture!(detector, input, Xoshiro(9140)) == fill(UInt16(3), 2, 3)
    @test capture!(detector, fill(1e9, 2, 3), Xoshiro(9140)) ==
        fill(typemax(UInt16), 2, 3)

    counts = detector.products.counts
    before = copy(counts)
    @test_throws DimensionMismatchError capture!(detector,
        fill(1.0, 3, 2), Xoshiro(9141))
    @test detector.products.counts === counts
    @test detector.products.counts == before

    for invalid_input in (
        fill(-1.0, 2, 3), fill(Inf, 2, 3), fill(NaN, 2, 3))
        @test_throws InvalidConfiguration capture!(detector,
            invalid_input, Xoshiro(9142))
    end
    for bad in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration SPADArrayDetector((1, 1);
            integration_time=bad)
    end
    for bad_dimensions in ((0, 1), (1, 0), (-1, 1))
        @test_throws InvalidConfiguration SPADArrayDetector(bad_dimensions)
    end
    for bad in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration SPADArraySensor(
            active_area_detection_efficiency=bad)
    end
    for bad in (-1.0, Inf, NaN)
        @test_throws InvalidConfiguration SPADArraySensor(
            dark_count_rate=bad)
    end
    for bad in (0.0, -0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration SPADArraySensor(fill_factor=bad)
    end
    for bad in (-1.0, Inf, NaN)
        @test_throws InvalidConfiguration SPADArraySensor(
            dead_time_model=NonParalyzableDeadTime(bad))
        @test_throws InvalidConfiguration SPADArraySensor(
            dead_time_model=ParalyzableDeadTime(bad))
        @test_throws InvalidConfiguration SPADArraySensor(
            mean_response_model=FirstOrderAfterpulseMeanResponse(bad))
    end
    for bad in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration SPADArraySensor(
            mean_response_model=NearestNeighborCountRedistribution(bad))
    end
    for bad in (0.0, -0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration SPADArrayDetector((1, 1);
            gate_model=DutyCycleGate(bad))
    end
    @test_throws InvalidConfiguration SPADArrayDetector((1, 1);
        noise=NoiseReadout(1.0))

    function replay_detector()
        return SPADArrayDetector((16, 16); integration_time=0.5,
            noise=NoisePhoton(), gate_model=DutyCycleGate(0.75),
            sensor=SPADArraySensor(active_area_detection_efficiency=0.7,
                dark_count_rate=0.25,
                dead_time_model=NonParalyzableDeadTime(1e-3),
                mean_response_model=FirstOrderAfterpulseMeanResponse(0.1)))
    end
    replay_input = fill(40.0, 16, 16)
    @test capture!(replay_detector(), replay_input, Xoshiro(9150)) ==
        capture!(replay_detector(), replay_input, Xoshiro(9150))

    allocation_detector = SPADArrayDetector((8, 8); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=0.5,
            mean_response_model=NearestNeighborCountRedistribution(0.1)))
    allocation_input = fill(2.0, 8, 8)
    rng = Xoshiro(9151)
    capture!(allocation_detector, allocation_input, rng)
    @test @inferred(capture!(allocation_detector,
        allocation_input, rng)) isa Matrix{Float64}
    @test_detector_allocation @allocated(capture!(allocation_detector,
        allocation_input, rng)) == 0
end
