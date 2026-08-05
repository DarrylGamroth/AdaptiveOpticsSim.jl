function pe05a_rate_map(values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(T(0.25), T(0.25)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.75e-6)))
    return IntensityMap(metadata, values)
end

function pe05a_captured_error(f)
    try
        f()
        return nothing
    catch error
        return error
    end
end

function assert_wfs_acquisition_owner_conformance!(
    prepared, observation, optical_product, rng)
    plan = @inferred WavefrontSensors.wfs_acquisition_plan(prepared)
    @test plan isa WavefrontSensors.AbstractWFSAcquisitionPlan
    @test isconcretetype(typeof(plan))
    @test isconcretetype(typeof(prepared))
    @test @inferred(WavefrontSensors.validate_wfs_acquisition_binding(
        observation, optical_product, prepared)) === nothing
    @test @inferred(acquire_wfs_observation!(
        observation, optical_product, prepared, rng)) === observation
    return plan
end

@testset "PE-05A generic WFS acquisition ownership" begin
    T = Float64
    coverage_enabled = coverage_instrumented()
    rate = pe05a_rate_map(reshape(T.(1:4), 2, 2))

    detector = Detector(exposure_duration=T(0.25), qe=one(T),
        noise=NoiseNone(), response_model=NullFrameResponse(), T=T)
    observation = WFSObservation(zeros(T, 2, 2);
        units=:electron_count, layout=:detector_frame)
    prepared = prepare_wfs_acquisition(detector, rate, observation)
    rng = Xoshiro(0x5045303541)
    plan = assert_wfs_acquisition_owner_conformance!(
        prepared, observation, rate, rng)
    @test plan isa WavefrontSensors.WFSDetectorAcquisitionPlan
    @test !ismutabletype(typeof(plan))
    @test plan.detector_plan === Detectors.detector_acquisition_plan(
        prepared.acquisition)
    @test plan.observation_metadata === observation.metadata
    @test all(field_type -> !(field_type <: AbstractArray),
        fieldtypes(typeof(plan)))
    @test validate_wfs_target(prepared, HostComputeDevice()) === prepared

    state_before = Detectors.detector_acquisition_state(
        prepared.acquisition)
    workspace_before = Detectors.detector_acquisition_workspace(
        prepared.acquisition)
    products_before = Detectors.detector_acquisition_products(
        prepared.acquisition)
    wrong_shape = WFSObservation(zeros(T, 3, 3);
        units=:electron_count, layout=:detector_frame)
    shape_error = pe05a_captured_error() do
        prepare_wfs_acquisition(detector, rate, wrong_shape)
    end
    @test shape_error isa WFSPreparationError
    @test shape_error.reason === :shape
    @test Detectors.detector_acquisition_state(prepared.acquisition) ===
        state_before
    @test Detectors.detector_acquisition_workspace(prepared.acquisition) ===
        workspace_before
    @test Detectors.detector_acquisition_products(prepared.acquisition) ===
        products_before
    @test WavefrontSensors.validate_wfs_acquisition_binding(
        observation, rate, prepared) === nothing

    input_alias = WFSObservation(rate.values;
        units=:electron_count, layout=:detector_frame)
    input_alias_error = pe05a_captured_error() do
        prepare_wfs_acquisition(detector, rate, input_alias)
    end
    @test input_alias_error isa WFSPreparationError
    @test input_alias_error.reason === :ownership
    detector_alias = WFSObservation(output_frame(detector);
        units=:electron_count, layout=:detector_frame)
    detector_alias_error = pe05a_captured_error() do
        prepare_wfs_acquisition(detector, rate, detector_alias)
    end
    @test detector_alias_error isa WFSPreparationError
    @test detector_alias_error.reason === :ownership

    observation_before_foreign = copy(observation.storage)
    rng_before_foreign = copy(rng)
    foreign_rate = pe05a_rate_map(copy(rate.values))
    foreign_error = pe05a_captured_error() do
        acquire_wfs_observation!(observation, foreign_rate, prepared, rng)
    end
    @test foreign_error isa WFSPreparationError
    @test foreign_error.reason === :prepared_binding
    @test observation.storage == observation_before_foreign
    @test rand(rng) == rand(rng_before_foreign)

    spad = SPADArrayDetector((2, 2); exposure_duration=T(0.5),
        noise=NoiseNone(), sensor=SPADArraySensor(
            active_area_detection_efficiency=one(T),
            dark_count_rate=zero(T), fill_factor=one(T), T=T), T=T)
    counting_observation = WFSObservation(zeros(T, 2, 2);
        units=:photon_count, layout=:counting_channels)
    counting = prepare_wfs_acquisition(spad, rate, counting_observation)
    counting_plan = assert_wfs_acquisition_owner_conformance!(
        counting, counting_observation, rate, rng)
    @test counting_plan isa WavefrontSensors.WFSCountingAcquisitionPlan
    @test !ismutabletype(typeof(counting_plan))
    @test counting_plan.input_metadata === rate.metadata
    @test counting_plan.observation_metadata === counting_observation.metadata
    @test all(field_type -> !(field_type <: AbstractArray),
        fieldtypes(typeof(counting_plan)))
    @test validate_wfs_target(counting, HostComputeDevice()) === counting

    counts_before = counting_array(spad)
    noise_before = counting_noise_buffer(spad)
    host_before = counting_host_buffer(spad)
    wrong_counting_shape = WFSObservation(zeros(T, 3, 3);
        units=:photon_count, layout=:counting_channels)
    counting_shape_error = pe05a_captured_error() do
        prepare_wfs_acquisition(spad, rate, wrong_counting_shape)
    end
    @test counting_shape_error isa WFSPreparationError
    @test counting_shape_error.reason === :shape
    @test counting_array(spad) === counts_before
    @test counting_noise_buffer(spad) === noise_before
    @test counting_host_buffer(spad) === host_before
    @test WavefrontSensors.validate_wfs_acquisition_binding(
        counting_observation, rate, counting) === nothing

    counting_input_alias = WFSObservation(rate.values;
        units=:photon_count, layout=:counting_channels)
    counting_input_alias_error = pe05a_captured_error() do
        prepare_wfs_acquisition(spad, rate, counting_input_alias)
    end
    @test counting_input_alias_error isa WFSPreparationError
    @test counting_input_alias_error.reason === :ownership
    counting_detector_alias = WFSObservation(counting_array(spad);
        units=:photon_count, layout=:counting_channels)
    counting_detector_alias_error = pe05a_captured_error() do
        prepare_wfs_acquisition(spad, rate, counting_detector_alias)
    end
    @test counting_detector_alias_error isa WFSPreparationError
    @test counting_detector_alias_error.reason === :ownership

    fanout_detectors = (
        Detector(exposure_duration=T(0.25), qe=one(T), noise=NoiseNone(),
            response_model=NullFrameResponse(), T=T),
        Detector(exposure_duration=T(0.5), qe=one(T), noise=NoiseNone(),
            response_model=NullFrameResponse(), T=T),
    )
    fanout_rates = (
        pe05a_rate_map(copy(rate.values)),
        pe05a_rate_map(copy(rate.values)),
    )
    fanout_observations = (
        WFSObservation(zeros(T, 2, 2);
            units=:electron_count, layout=:detector_frame),
        WFSObservation(zeros(T, 2, 2);
            units=:electron_count, layout=:detector_frame),
    )
    fanout = prepare_wfs_acquisition(
        fanout_detectors, fanout_rates, fanout_observations)
    fanout_plan = assert_wfs_acquisition_owner_conformance!(
        fanout, fanout_observations, fanout_rates, rng)
    @test fanout_plan isa
        WavefrontSensors.WFSMultipleDetectorAcquisitionPlan
    @test fanout_plan.component_plans === map(
        WavefrontSensors.wfs_acquisition_plan, fanout.acquisitions)
    @test all(field_type -> !(field_type <: AbstractArray),
        fieldtypes(typeof(fanout_plan)))
    @test validate_wfs_target(fanout, HostComputeDevice()) === fanout

    atomic_observations = (
        WFSObservation(zeros(T, 2, 2);
            units=:electron_count, layout=:detector_frame),
        WFSObservation(zeros(T, 3, 3);
            units=:electron_count, layout=:detector_frame),
    )
    atomic_error = pe05a_captured_error() do
        prepare_wfs_acquisition(
            fanout_detectors, fanout_rates, atomic_observations)
    end
    @test atomic_error isa WFSPreparationError
    @test atomic_error.reason === :shape
    @test WavefrontSensors.validate_wfs_acquisition_binding(
        fanout_observations, fanout_rates, fanout) === nothing

    duplicate_detector_error = pe05a_captured_error() do
        prepare_wfs_acquisition(
            (first(fanout_detectors), first(fanout_detectors)),
            fanout_rates, atomic_observations)
    end
    @test duplicate_detector_error isa WFSPreparationError
    @test duplicate_detector_error.reason === :ownership
    shared_observation = WFSObservation(zeros(T, 2, 2);
        units=:electron_count, layout=:detector_frame)
    shared_observation_error = pe05a_captured_error() do
        prepare_wfs_acquisition(fanout_detectors, fanout_rates,
            (shared_observation, shared_observation))
    end
    @test shared_observation_error isa WFSPreparationError
    @test shared_observation_error.reason === :ownership

    acquire_wfs_observation!(
        fanout_observations, fanout_rates, fanout, rng)
    if coverage_enabled
        @test_skip "PE-05A fan-out allocation assertion is disabled under coverage instrumentation"
    else
        @test @allocated(acquire_wfs_observation!(
            fanout_observations, fanout_rates, fanout, rng)) == 0
    end

    prepare_wfs_acquisition(last(fanout_detectors), last(fanout_rates),
        last(fanout_observations))
    fanout_before = map(value -> copy(value.storage), fanout_observations)
    rng_before_stale = copy(rng)
    stale_error = pe05a_captured_error() do
        acquire_wfs_observation!(fanout_observations, fanout_rates,
            fanout, rng)
    end
    @test stale_error isa WFSPreparationError
    @test stale_error.reason === :prepared_binding
    @test first(fanout_observations).storage == first(fanout_before)
    @test last(fanout_observations).storage == last(fanout_before)
    @test rand(rng) == rand(rng_before_stale)

    acquire_wfs_observation!(observation, rate, prepared, rng)
    acquire_wfs_observation!(counting_observation, rate, counting, rng)
    if coverage_enabled
        @test_skip "PE-05A allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(acquire_wfs_observation!(
            observation, rate, prepared, rng)) == 0
        @test @allocated(acquire_wfs_observation!(
            counting_observation, rate, counting, rng)) == 0
    end
end
