@testset "InGaAs detector" begin
    zero_psf = zeros(4, 4)
    det_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=InGaAsSensor(glow_rate=3.0))
    frame_ingaas = capture!(det_ingaas, zero_psf; rng=MersenneTwister(13))
    @test sum(frame_ingaas) > 0
    @test supports_sensor_glow(det_ingaas.params.sensor)
    @test detector_export_metadata(det_ingaas).frame_response == :gaussian
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=-1.0)
    det_ingaas_persist = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        sensor=InGaAsSensor(persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(det_ingaas_persist, fill(4.0, 4, 4); rng=MersenneTwister(121))
    persisted = capture!(det_ingaas_persist, zeros(4, 4); rng=MersenneTwister(122))
    @test sum(persisted) ≈ 32.0

    persistence_sensor = InGaAsSensor(
        persistence_model=ExponentialPersistence(0.5, 0.0))
    persistence_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=1, response_model=NullFrameResponse(),
        sensor=persistence_sensor)
    persistence_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=1, response_model=NullFrameResponse(),
        sensor=persistence_sensor)
    persistence_seed = fill(4.0, 4, 4)
    persistence_dark = zeros(4, 4)
    capture!(persistence_whole, persistence_seed; rng=MersenneTwister(123))
    capture!(persistence_split, persistence_seed; rng=MersenneTwister(123))
    persistence_whole_frame = copy(capture!(persistence_whole,
        persistence_dark; rng=MersenneTwister(124)))
    persistence_partial_frame = copy(capture!(persistence_split,
        persistence_dark; rng=MersenneTwister(124), integration_duration=0.5))
    persistence_split_frame = copy(capture!(persistence_split,
        persistence_dark; rng=MersenneTwister(124), integration_duration=0.5))
    @test persistence_partial_frame == persistence_whole_frame
    @test persistence_split_frame == persistence_whole_frame
    @test persistence_split.state.latent_buffer ==
        persistence_whole.state.latent_buffer

    persistence_prepared_whole = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_prepared_split = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_rate_map = detector_test_intensity_map(persistence_dark)
    persistence_whole_plan = prepare_detector_acquisition(
        persistence_prepared_whole, persistence_rate_map)
    persistence_split_plan = prepare_detector_acquisition(
        persistence_prepared_split, persistence_rate_map)
    capture!(persistence_prepared_whole, persistence_seed;
        rng=MersenneTwister(125))
    capture!(persistence_prepared_split, persistence_seed;
        rng=MersenneTwister(125))
    persistence_prepared_whole_frame = copy(capture!(
        persistence_prepared_whole, persistence_rate_map,
        persistence_whole_plan; rng=MersenneTwister(126)))
    capture!(persistence_prepared_split, persistence_rate_map,
        persistence_split_plan; rng=MersenneTwister(126), integration_duration=0.5)
    persistence_prepared_split_frame = copy(capture!(
        persistence_prepared_split, persistence_rate_map,
        persistence_split_plan; rng=MersenneTwister(126), integration_duration=0.5))
    @test persistence_prepared_split_frame ==
        persistence_prepared_whole_frame
    @test persistence_prepared_split.state.latent_buffer ==
        persistence_prepared_whole.state.latent_buffer

    persistence_allocation_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_allocation_plan = prepare_detector_acquisition(
        persistence_allocation_detector, persistence_rate_map)
    capture!(persistence_allocation_detector, persistence_seed;
        rng=Xoshiro(127))
    @test_detector_allocation prepared_first_incremental_capture_allocations(
        persistence_allocation_detector, persistence_rate_map,
        persistence_allocation_plan, Xoshiro(128), 0.5) == 0
    persist_meta = detector_export_metadata(det_ingaas_persist)
    @test persist_meta.persistence_model == :exponential
    det_ingaas_nonlinear = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=InGaAsSensor())
    nonlinear_frame = capture!(det_ingaas_nonlinear, fill(10.0, 2, 2); rng=MersenneTwister(123))
    @test nonlinear_frame == fill(5.0, 2, 2)
    nonlinear_meta = detector_export_metadata(det_ingaas_nonlinear)
    @test nonlinear_meta.nonlinearity_model == :saturating
    @test supports_detector_persistence(det_ingaas_persist.params.sensor)
    @test supports_detector_nonlinearity(det_ingaas_nonlinear.params.sensor)
    @test_throws InvalidConfiguration InGaAsSensor(persistence_model=ExponentialPersistence(1.1, 0.0))
    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        nonlinearity_model=SaturatingFrameNonlinearity(-0.1), sensor=InGaAsSensor())

end
