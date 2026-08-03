function skipper_noise_samples(n_samples::Int, seed::Int)
    detector = Detector(
        integration_time=1.0,
        noise=NoiseReadout(4.0),
        qe=1.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(n_samples),
            sample_duration=2e-6),
        response_model=NullFrameResponse(),
    )
    return vec(copy(capture!(
        detector, zeros(128, 128), Xoshiro(seed))))
end

function test_skipper_gaussian_moments(
    samples, expected_variance; sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit =
        sigma_limit * sqrt(expected_variance / sample_count)
    variance_limit = sigma_limit * expected_variance *
        sqrt(2 / (sample_count - 1))
    @test abs(mean(samples)) <= mean_limit
    @test abs(var(samples) - expected_variance) <= variance_limit
end

@testset "Skipper CCD bounded independent-read contract" begin
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=SkipperSampling(0))
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=SkipperSampling(2), sample_duration=-1.0)
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=SkipperSampling(2), sample_duration=Inf)

    detector = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        gain=2.0,
        full_well=10.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(8),
            sample_duration=2e-6),
        response_model=NullFrameResponse(),
    )
    input = fill(20.0, 8, 8)
    output = capture!(detector, input, Xoshiro(9101))
    @test output == fill(20.0, 8, 8)
    @test detector.workspace.readout.baseline_frame == fill(10.0, 8, 8)
    @test detector.workspace.readout.sample_sum == fill(80.0, 8, 8)

    products = readout_products(detector)
    @test products isa SkipperReadoutProducts
    @test products.sample_count == 8
    @test detector.params.sensor.sample_duration == 2e-6
    @test products.mean_frame == output
    @test detector_signal_frame(detector) == output
    @test detector_combined_frame(detector) == output
    @test detector_read_cube(detector) === nothing
    @test fieldnames(typeof(products)) == (:mean_frame, :sample_count)
    @test size(products.mean_frame) == size(input)

    metadata = detector_export_metadata(detector)
    @test metadata.sampling_mode == :skipper
    @test metadata.sampling_reads == 8
    @test metadata.sampling_reference_reads == 0
    @test metadata.sampling_signal_reads == 8
    @test metadata.sampling_read_time == 2e-6
    @test metadata.sampling_wallclock_time == 1.0 + 16e-6
    @test metadata.provides_signal_frame
    @test !metadata.provides_read_cube
    @test supports_nondestructive_reads(detector.params.sensor)
    @test supports_multi_read_readout_products(detector.params.sensor)

    many_samples = Detector(
        noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(64),
            sample_duration=2e-6),
        response_model=NullFrameResponse())
    capture!(many_samples, zeros(8, 8), Xoshiro(9102))
    many_products = readout_products(many_samples)
    @test typeof(many_products) === typeof(products)
    @test size(many_products.mean_frame) == size(products.mean_frame)
    @test detector_read_cube(many_samples) === nothing

    cube = fill(3.0, 2, 8, 8)
    original = copy(cube)
    @test_throws InvalidConfiguration capture_stack!(
        detector, cube, similar(cube), Xoshiro(9103))
    @test cube == original
end

@testset "Skipper CCD independent Gaussian read-noise qualification" begin
    single_sample_variance = 16.0
    for (case_index, n_samples) in enumerate((1, 4, 16, 64))
        samples = skipper_noise_samples(n_samples, 9110 + case_index)
        test_skipper_gaussian_moments(
            samples, single_sample_variance / n_samples)
    end
end

@testset "Skipper CCD prepared execution" begin
    input = fill(2.0, 16, 16)
    detector = Detector(
        integration_time=0.5,
        noise=NoiseReadout(0.25),
        qe=1.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(16),
            sample_duration=1e-6),
        response_model=NullFrameResponse(),
    )
    rng = Xoshiro(9120)
    @test @inferred(capture!(detector, input, rng)) isa Matrix{Float64}
    @test_detector_allocation @allocated(
        capture!(detector, input, rng)) == 0

    transition = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        response_model=NullFrameResponse(),
        sensor=CCDSensor(sampling_mode=SkipperSampling(2)))
    capture!(transition, ones(2, 2), Xoshiro(9121))
    @test transition.workspace.readout.baseline_frame == ones(2, 2)
    capture!(transition, zeros(2, 2);
        rng=Xoshiro(9122), integration_duration=0.5)
    incremental_frame = copy(capture!(
        transition, zeros(2, 2);
        rng=Xoshiro(9123), integration_duration=0.5))
    @test all(iszero, incremental_frame)
    @test readout_ready(transition)
end
