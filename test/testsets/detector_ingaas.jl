function test_ingaas_poisson_moments(samples, expected_mean;
    sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_mean / sample_count)
    variance_limit = sigma_limit * expected_mean *
        sqrt(2 / (sample_count - 1))
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_mean) <= variance_limit
end

@testset "InGaAs construction and explicit response" begin
    sensor = InGaAsSensor(glow_rate=3.0)
    @test isconcretetype(typeof(sensor))
    @test !ismutabletype(typeof(sensor))
    @test sensor.glow_rate == 3.0
    @test sensor.persistence_model isa NullPersistence
    @test supports_sensor_glow(sensor)
    @test supports_detector_defect_maps(sensor)
    @test supports_detector_persistence(sensor)
    @test supports_detector_nonlinearity(sensor)

    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=-1.0)
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=Inf)
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=NaN)
    for model in (
        ExponentialPersistence(-0.1, 0.5),
        ExponentialPersistence(1.1, 0.5),
        ExponentialPersistence(0.5, -0.1),
        ExponentialPersistence(0.5, 1.1),
        ExponentialPersistence(NaN, 0.5),
        ExponentialPersistence(0.5, NaN),
    )
        @test_throws InvalidConfiguration InGaAsSensor(
            persistence_model=model)
    end

    impulse = zeros(7, 7)
    impulse[4, 4] = 1.0
    default_detector = Detector(noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor())
    default_output = copy(capture!(
        default_detector, impulse, Xoshiro(9301)))
    default_metadata = detector_export_metadata(default_detector)
    @test default_output == impulse
    @test default_metadata.frame_response == :none
    @test !supports_detector_mtf(default_detector)
    @test detector_mtf(default_detector, 0.2, 0.3) == 1.0

    explicit_detector = Detector(noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor(),
        response_model=GaussianPixelResponse(response_width_px=0.4))
    explicit_output = copy(capture!(
        explicit_detector, impulse, Xoshiro(9302)))
    @test detector_export_metadata(explicit_detector).frame_response ==
        :gaussian
    @test supports_detector_mtf(explicit_detector)
    @test detector_mtf(explicit_detector, 0.2, 0.3) < 1.0
    @test explicit_output != impulse
end

@testset "InGaAs deterministic frame pipeline" begin
    detector = Detector(
        integration_time=2.0,
        noise=NoiseNone(),
        qe=0.5,
        gain=2.0,
        full_well=5.0,
        sensor=InGaAsSensor(),
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
    )
    output = capture!(detector, fill(20.0, 4, 4), Xoshiro(9310))
    @test output == fill(10.0, 4, 4)
    @test detector_export_metadata(detector).nonlinearity_model ==
        :saturating
end


@testset "InGaAs glow and dark-current moments" begin
    detector = Detector(
        integration_time=2.0,
        noise=NoiseNone(),
        qe=1.0,
        dark_current=1.5,
        sensor=InGaAsSensor(glow_rate=2.5),
        response_model=NullFrameResponse(),
    )
    samples = vec(copy(capture!(
        detector, zeros(128, 128), Xoshiro(9320))))
    test_ingaas_poisson_moments(samples, 8.0)

    short = Detector(integration_time=0.25, noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor(glow_rate=20.0),
        response_model=NullFrameResponse())
    long = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor(glow_rate=20.0),
        response_model=NullFrameResponse())
    prepare_detector_buffers!(short, (2, 2))
    prepare_detector_buffers!(long, (2, 2))
    fill!(short.state.frame, 0.0)
    fill!(long.state.frame, 0.0)
    apply_sensor_statistics!(short.params.sensor, short,
        Xoshiro(9321), 0.25)
    apply_sensor_statistics!(long.params.sensor, long,
        Xoshiro(9321), 0.25)
    @test short.state.frame == long.state.frame
end

@testset "InGaAs charge-domain exponential persistence" begin
    function persistence_detector(gain)
        return Detector(
            integration_time=1.0,
            noise=NoiseNone(),
            qe=1.0,
            gain=gain,
            response_model=NullFrameResponse(),
            sensor=InGaAsSensor(persistence_model=
                ExponentialPersistence(0.25, 0.5)),
        )
    end

    unit_gain = persistence_detector(1.0)
    four_gain = persistence_detector(4.0)
    illumination = fill(8.0, 4, 4)
    dark = zeros(4, 4)
    expected_charge = (8.0, 2.0, 1.5)
    expected_latent = (2.0, 1.5, 1.125)
    for (frame_index, input) in enumerate((illumination, dark, dark))
        unit_output = copy(capture!(
            unit_gain, input, Xoshiro(9330 + frame_index)))
        four_output = copy(capture!(
            four_gain, input, Xoshiro(9330 + frame_index)))
        @test unit_output == fill(expected_charge[frame_index], 4, 4)
        @test four_output == 4 .* unit_output
        @test unit_gain.state.latent_buffer ==
            fill(expected_latent[frame_index], 4, 4)
        @test four_gain.state.latent_buffer ==
            unit_gain.state.latent_buffer
    end
    @test detector_export_metadata(unit_gain).persistence_model ==
        :exponential

    split = persistence_detector(1.0)
    whole = persistence_detector(1.0)
    capture!(split, illumination, Xoshiro(9340))
    capture!(whole, illumination, Xoshiro(9340))
    capture!(split, dark; rng=Xoshiro(9341), integration_duration=0.5)
    split_frame = copy(capture!(split, dark;
        rng=Xoshiro(9342), integration_duration=0.5))
    whole_frame = copy(capture!(whole, dark, Xoshiro(9341)))
    @test split_frame == whole_frame
    @test split.state.latent_buffer == whole.state.latent_buffer

    cube = fill(3.0, 2, 4, 4)
    original = copy(cube)
    @test_throws InvalidConfiguration capture_stack!(
        split, cube, similar(cube), Xoshiro(9343))
    @test cube == original
end

@testset "InGaAs replay, inference, and prepared allocation" begin
    function replay_detector()
        return Detector(
            integration_time=0.5,
            noise=NoisePhotonReadout(0.25),
            qe=0.75,
            dark_current=0.5,
            sensor=InGaAsSensor(
                glow_rate=0.25,
                persistence_model=ExponentialPersistence(0.1, 0.4)),
            response_model=NullFrameResponse(),
        )
    end
    input = fill(4.0, 16, 16)
    first = replay_detector()
    second = replay_detector()
    @test capture!(first, input, Xoshiro(9350)) ==
        capture!(second, input, Xoshiro(9350))
    @test first.state.latent_buffer == second.state.latent_buffer

    allocation_detector = Detector(
        integration_time=0.5,
        noise=NoiseNone(),
        qe=1.0,
        sensor=InGaAsSensor(
            persistence_model=ExponentialPersistence(0.1, 0.4)),
        response_model=NullFrameResponse(),
    )
    rate_map = detector_test_intensity_map(input)
    plan = prepare_detector_acquisition(allocation_detector, rate_map)
    rng = Xoshiro(9351)
    capture!(allocation_detector, rate_map, plan, rng)
    @test @inferred(capture!(
        allocation_detector, rate_map, plan, rng)) isa Matrix{Float64}
    @test_detector_allocation @allocated(capture!(
        allocation_detector, rate_map, plan, rng)) == 0
end
