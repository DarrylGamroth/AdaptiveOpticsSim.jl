function pe04_direct_memory_storage(::Type{T}) where {T}
    T isa Union && return any(pe04_direct_memory_storage, Base.uniontypes(T))
    return T isa DataType && T <: Core.GenericMemory
end

function pe04_assert_owner_storage_contract(owner)
    @test isconcretetype(typeof(owner))
    for field_type in fieldtypes(typeof(owner))
        @test field_type !== Any
        @test !isabstracttype(field_type)
        @test !pe04_direct_memory_storage(field_type)
    end
    return nothing
end

function pe04_state_values(state)
    return (
        accum_buffer=copy(state.accum_buffer),
        latent_buffer=copy(state.latent_buffer),
        integrated_time=state.integrated_time,
        readout_ready=state.readout_ready,
    )
end

function pe04_state_values_match(state, values)
    return state.accum_buffer == values.accum_buffer &&
        state.latent_buffer == values.latent_buffer &&
        state.integrated_time == values.integrated_time &&
        state.readout_ready == values.readout_ready
end

function pe04_rejected_capture_preserves_state_and_rng(prepared, rng)
    before = pe04_state_values(detector_acquisition_state(prepared))
    expected_rng = copy(rng)
    error = try
        capture!(prepared, rng)
        nothing
    catch caught
        caught
    end
    return error isa InvalidConfiguration &&
        pe04_state_values_match(detector_acquisition_state(prepared), before) &&
        rand(rng) == rand(expected_rng)
end

function pe04_captured_error(f)
    try
        f()
        return nothing
    catch error
        return error
    end
end

function pe04_owner_fields(owner)
    return ntuple(index -> getfield(owner, index), fieldcount(typeof(owner)))
end

function pe04_frame_detector(;
    sensor=CCDSensor(),
    noise=NoiseNone(),
    readout_window=nothing,
)
    return Detector(
        exposure_duration=1.0,
        noise=noise,
        qe=1.0,
        sensor=sensor,
        readout_window=readout_window,
        response_model=NullFrameResponse(),
    )
end

@testset "PE-04 prepared detector owner boundaries" begin
    input = detector_test_intensity_map(reshape(Float64.(1:16), 4, 4))
    detector = pe04_frame_detector()
    prepared = prepare_detector_acquisition(detector, input)
    plan = detector_acquisition_plan(prepared)
    state = detector_acquisition_state(prepared)
    workspace = detector_acquisition_workspace(prepared)
    products = detector_acquisition_products(prepared)

    @test Detectors.detector_acquisition_detector(prepared) === detector
    @test Detectors.detector_acquisition_input(prepared) === input
    @test Detectors.detector_acquisition_plan(prepared) === plan
    @test Detectors.detector_acquisition_state(prepared) === state
    @test Detectors.detector_acquisition_workspace(prepared) === workspace
    @test Detectors.detector_acquisition_products(prepared) === products
    @test state === detector.state
    @test workspace === detector.workspace
    @test products === detector.products
    @test state !== workspace
    @test state !== products
    @test workspace !== products
    @test fieldnames(typeof(plan)) == (
        :detector_params,
        :input_metadata,
        :input_shape,
        :frame_shape,
        :output_shape,
        :rate_scale,
        :quantum_efficiency,
    )
    @test !ismutabletype(typeof(plan))
    @test all(name -> !(getfield(plan, name) isa AbstractArray),
        fieldnames(typeof(plan)))
    @test isdisjoint(Set(fieldnames(typeof(plan))), Set((
        :detector,
        :input,
        :state,
        :workspace,
        :products,
    )))

    for owner in (
        plan,
        prepared,
        state,
        workspace,
        products,
        prepared.state_binding,
        prepared.workspace_binding,
        prepared.product_binding,
    )
        pe04_assert_owner_storage_contract(owner)
    end

    @test !Base.mightalias(state.accum_buffer,
        workspace.presampling_buffer)
    @test !Base.mightalias(state.accum_buffer, products.frame)
    @test !Base.mightalias(workspace.presampling_buffer, products.frame)
    @test !Base.mightalias(input.values, state.accum_buffer)
    @test !Base.mightalias(input.values, workspace.presampling_buffer)
    @test !Base.mightalias(input.values, products.frame)

    rng = Xoshiro(23040)
    @test @inferred(capture!(prepared, rng)) === output_frame(detector)
    @test_detector_allocation prepared_detector_capture_allocations(
        prepared, rng) == 0

    @test !applicable(capture!, detector, input, plan, Xoshiro(1))
    @test !applicable(capture!, detector, input.values, plan, Xoshiro(1))
end

@testset "PE-04 exact bindings reject replacement atomically" begin
    input = detector_test_intensity_map(fill(2.0, 3, 3))
    detector = pe04_frame_detector()
    prepared = prepare_detector_acquisition(detector, input)

    for owner_name in (:state, :workspace, :products)
        original = getfield(detector, owner_name)
        setfield!(detector, owner_name, deepcopy(original))
        @test pe04_rejected_capture_preserves_state_and_rng(
            prepared, Xoshiro(23041))
        setfield!(detector, owner_name, original)
    end

    state = detector.state
    original_accumulator = state.accum_buffer
    state.accum_buffer = similar(original_accumulator)
    fill!(state.accum_buffer, 11.0)
    @test pe04_rejected_capture_preserves_state_and_rng(
        prepared, Xoshiro(23042))
    state.accum_buffer = original_accumulator

    workspace = detector.workspace
    original_scratch = workspace.presampling_scratch
    workspace.presampling_scratch = similar(original_scratch)
    fill!(workspace.presampling_scratch, 12.0)
    @test pe04_rejected_capture_preserves_state_and_rng(
        prepared, Xoshiro(23043))
    workspace.presampling_scratch = original_scratch

    products = detector.products
    original_frame = products.frame
    products.frame = similar(original_frame)
    fill!(products.frame, 13.0)
    @test pe04_rejected_capture_preserves_state_and_rng(
        prepared, Xoshiro(23044))
    products.frame = original_frame

    foreign_input = detector_test_intensity_map(copy(input.values))
    before = pe04_state_values(detector.state)
    rng = Xoshiro(23045)
    expected_rng = copy(rng)
    @test_throws InvalidConfiguration Detectors._require_prepared_acquisition(
        prepared, foreign_input)
    @test pe04_state_values_match(detector.state, before)
    @test rand(rng) == rand(expected_rng)

    foreign_detector = pe04_frame_detector()
    @test_throws InvalidConfiguration Detectors._require_prepared_detector_binding(
        foreign_detector, prepared)
    @test capture!(prepared, Xoshiro(23046)) === output_frame(detector)
end

@testset "PE-04 preparation failure and workspace recreation" begin
    input = detector_test_intensity_map(fill(3.0, 3, 3))
    detector = pe04_frame_detector()
    prepared = prepare_detector_acquisition(detector, input)
    owners_before = (detector.state, detector.workspace, detector.products)
    state_before = pe04_state_values(detector.state)
    workspace_before = copy(detector.workspace.presampling_buffer)
    products_before = copy(detector.products.frame)

    invalid_input = detector_test_intensity_map([
        -1.0 1.0 1.0
         1.0 1.0 1.0
         1.0 1.0 1.0
    ])
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        detector, invalid_input)
    @test detector.state === owners_before[1]
    @test detector.workspace === owners_before[2]
    @test detector.products === owners_before[3]
    @test pe04_state_values_match(detector.state, state_before)
    @test isequal(detector.workspace.presampling_buffer, workspace_before)
    @test isequal(detector.products.frame, products_before)

    original_workspace_buffer = detector.workspace.presampling_buffer
    detector.workspace.presampling_buffer = detector.state.accum_buffer
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        detector, input)
    @test detector.state === owners_before[1]
    @test detector.workspace === owners_before[2]
    @test detector.products === owners_before[3]
    detector.workspace.presampling_buffer = original_workspace_buffer
    @test capture!(prepared, Xoshiro(23047)) === output_frame(detector)

    stochastic_input = detector_test_intensity_map(fill(7.0, 4, 4))
    rebuilt = pe04_frame_detector(noise=NoisePhoton())
    control = pe04_frame_detector(noise=NoisePhoton())
    stale = prepare_detector_acquisition(rebuilt, stochastic_input)
    control_prepared = prepare_detector_acquisition(
        control, detector_test_intensity_map(copy(stochastic_input.values)))
    previous_workspace = rebuilt.workspace
    rebuilt.workspace = deepcopy(previous_workspace)
    rebuilt_prepared = prepare_detector_acquisition(rebuilt, stochastic_input)
    @test detector_acquisition_workspace(rebuilt_prepared) !==
        previous_workspace
    @test_throws InvalidConfiguration capture!(stale, Xoshiro(23048))
    @test copy(capture!(rebuilt_prepared, Xoshiro(23049))) ==
        copy(capture!(control_prepared, Xoshiro(23049)))
end

@testset "PE-04 lifecycle preparation is transactional" begin
    input = detector_test_intensity_map(fill(2.0, 4, 4))

    global_detector = pe04_frame_detector()
    global_owner = prepare_detector_acquisition(global_detector, input)
    global_owners = (global_detector.state, global_detector.workspace,
        global_detector.products)
    @test_throws DetectorAcquisitionError prepare_global_shutter_acquisition(
        global_detector, input,
        GlobalShutterAcquisitionDefinition(PlantDuration(500_000_000)))
    @test (global_detector.state, global_detector.workspace,
        global_detector.products) === global_owners
    @test capture!(global_owner, Xoshiro(23070)) ===
        output_frame(global_detector)

    rolling_detector = pe04_frame_detector(sensor=CMOSSensor(
        timing_model=RollingShutter(1e-6)))
    rolling_owner = prepare_detector_acquisition(rolling_detector, input)
    rolling_owners = (rolling_detector.state, rolling_detector.workspace,
        rolling_detector.products)
    @test_throws DetectorAcquisitionError prepare_rolling_shutter_acquisition(
        rolling_detector, input,
        RollingShutterAcquisitionDefinition(PlantDuration(500_000_000)))
    @test (rolling_detector.state, rolling_detector.workspace,
        rolling_detector.products) === rolling_owners
    @test capture!(rolling_owner, Xoshiro(23071)) ===
        output_frame(rolling_detector)

    transfer_detector = pe04_frame_detector(sensor=EMCCDSensor(
        excess_noise_factor=1.0,
        acquisition_mode=FrameTransferAcquisition(transfer_duration=0.01)))
    transfer_owner = prepare_detector_acquisition(transfer_detector, input)
    transfer_owners = (transfer_detector.state, transfer_detector.workspace,
        transfer_detector.products)
    @test_throws DetectorAcquisitionError prepare_frame_transfer_acquisition(
        transfer_detector, input,
        FrameTransferAcquisitionDefinition(PlantDuration(500_000_000)))
    @test (transfer_detector.state, transfer_detector.workspace,
        transfer_detector.products) === transfer_owners
    @test capture!(transfer_owner, Xoshiro(23072)) ===
        output_frame(transfer_detector)
end


@testset "PE-04 readout workspace and product ownership" begin
    input = detector_test_intensity_map(fill(5.0, 4, 4))

    skipper = pe04_frame_detector(
        sensor=CCDSensor(sampling_mode=SkipperSampling(4)))
    skipper_prepared = prepare_detector_acquisition(skipper, input)
    capture!(skipper_prepared, Xoshiro(23050))
    @test skipper.products.readout isa SkipperReadoutProducts
    @test skipper.workspace.readout isa Detectors.SkipperReadoutWorkspace
    pe04_assert_owner_storage_contract(skipper.products.readout)
    pe04_assert_owner_storage_contract(skipper.workspace.readout)
    @test fieldnames(typeof(skipper.products.readout)) ==
        (:mean_frame, :sample_count)
    @test fieldnames(typeof(skipper.workspace.readout)) ==
        (:baseline_frame, :sample_sum)
    @test !Base.mightalias(skipper.products.readout.mean_frame,
        skipper.workspace.readout.baseline_frame)

    fowler = pe04_frame_detector(sensor=HgCdTeSensor(
        read_duration=0.1, sampling_mode=FowlerSampling(2)),
        readout_window=FrameWindow(2:3, 2:3))
    fowler_prepared = prepare_detector_acquisition(fowler, input)
    capture!(fowler_prepared, Xoshiro(23051))
    @test fowler.products.readout isa MultiReadFrameReadoutProducts
    @test fowler.workspace.readout isa
        Detectors.MultiReadFrameReadoutWorkspace
    pe04_assert_owner_storage_contract(fowler.products.readout)
    pe04_assert_owner_storage_contract(fowler.workspace.readout)
    @test :reference_average ∉ fieldnames(typeof(fowler.products.readout))
    @test :signal_average ∉ fieldnames(typeof(fowler.products.readout))
    @test detector_read_offsets_s(fowler) isa FixedSizeVector

    ramp = pe04_frame_detector(sensor=HgCdTeSensor(
        read_duration=0.1, sampling_mode=UpTheRampSampling(3)),
        readout_window=FrameWindow(2:3, 2:3))
    ramp_prepared = prepare_detector_acquisition(ramp, input)
    capture!(ramp_prepared, Xoshiro(23052))
    @test ramp.products.readout isa UpTheRampReadoutProducts
    @test ramp.workspace.readout isa Detectors.UpTheRampReadoutWorkspace
    pe04_assert_owner_storage_contract(ramp.products.readout)
    pe04_assert_owner_storage_contract(ramp.workspace.readout)
    @test :slope ∉ fieldnames(typeof(ramp.products.readout))
    @test :cube ∉ fieldnames(typeof(ramp.products.readout))
    @test :slope_frame ∉ fieldnames(typeof(ramp.workspace.readout))
    @test :read_cube ∉ fieldnames(typeof(ramp.workspace.readout))
    @test detector_read_offsets_s(ramp) isa FixedSizeVector
    @test !Base.mightalias(ramp.products.readout.slope_frame,
        ramp.workspace.readout.slope)

    lifecycle = prepare_global_shutter_acquisition(
        ramp,
        input,
        GlobalShutterAcquisitionDefinition(
            PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(100_000_000)),
    )
    @test lifecycle.read_offsets isa FixedSizeVector
    @test lifecycle.read_offset_binding isa FixedSizeVector
    @test !applicable(prepare_global_shutter_acquisition,
        lifecycle.acquisition, lifecycle.definition)
    lifecycle_state = GlobalShutterAcquisitionState(lifecycle)
    @test lifecycle_state isa GlobalShutterAcquisitionState
    @test_throws DetectorAcquisitionError GlobalShutterAcquisitionState(
        lifecycle)
    lifecycle_fields = pe04_owner_fields(lifecycle_state)
    @test !applicable(GlobalShutterAcquisitionState, lifecycle_fields...)
    forged_lifecycle_state = GlobalShutterAcquisitionState(
        AdaptiveOpticsSim.Plant._DETECTOR_LIFECYCLE_OWNER_TOKEN,
        lifecycle_fields...)
    forged_lifecycle_error = pe04_captured_error() do
        begin_exposure!(lifecycle, forged_lifecycle_state, PlantTimestamp(0))
    end
    @test forged_lifecycle_error isa DetectorAcquisitionError
    @test forged_lifecycle_error.reason == :foreign_state

    rolling_detector = pe04_frame_detector(sensor=CMOSSensor(
        timing_model=RollingShutter(1e-6)))
    rolling = prepare_rolling_shutter_acquisition(
        rolling_detector,
        input,
        RollingShutterAcquisitionDefinition(PlantDuration(1_000_000_000)),
    )
    rolling_state = RollingShutterAcquisitionState(rolling)
    @test rolling_state isa RollingShutterAcquisitionState
    @test_throws DetectorAcquisitionError RollingShutterAcquisitionState(
        rolling)
    rolling_fields = pe04_owner_fields(rolling_state)
    @test !applicable(RollingShutterAcquisitionState, rolling_fields...)
    forged_rolling_state = RollingShutterAcquisitionState(
        AdaptiveOpticsSim.Plant._DETECTOR_LIFECYCLE_OWNER_TOKEN,
        rolling_fields...)
    forged_rolling_error = pe04_captured_error() do
        begin_exposure!(rolling, forged_rolling_state, PlantTimestamp(0))
    end
    @test forged_rolling_error isa DetectorAcquisitionError
    @test forged_rolling_error.reason == :foreign_state

    transfer_detector = pe04_frame_detector(sensor=EMCCDSensor(
        excess_noise_factor=1.0,
        acquisition_mode=FrameTransferAcquisition(transfer_duration=0.01)))
    transfer = prepare_frame_transfer_acquisition(
        transfer_detector,
        input,
        FrameTransferAcquisitionDefinition(
            PlantDuration(1_000_000_000)),
    )
    transfer_state = FrameTransferAcquisitionState(transfer)
    @test_throws DetectorAcquisitionError FrameTransferAcquisitionState(
        transfer)
    transfer_fields = pe04_owner_fields(transfer_state)
    @test !applicable(FrameTransferAcquisitionState, transfer_fields...)
    forged_transfer_state = FrameTransferAcquisitionState(
        AdaptiveOpticsSim.Plant._DETECTOR_LIFECYCLE_OWNER_TOKEN,
        transfer_fields...)
    forged_transfer_error = pe04_captured_error() do
        begin_exposure!(transfer, forged_transfer_state, PlantTimestamp(0))
    end
    @test forged_transfer_error isa DetectorAcquisitionError
    @test forged_transfer_error.reason == :foreign_state
    @test !hasfield(typeof(transfer), :storage_frame)
    @test hasfield(typeof(transfer_state), :storage_frame)
    @test !Base.mightalias(transfer_state.storage_frame,
        transfer_detector.state.accum_buffer)
    @test !Base.mightalias(transfer_state.storage_frame,
        transfer_detector.workspace.presampling_buffer)
    @test !Base.mightalias(transfer_state.storage_frame,
        transfer_detector.products.frame)
end

@testset "PE-04 representative detector family storage" begin
    input = detector_test_intensity_map(fill(1.0, 2, 2))
    for sensor in (
        CCDSensor(),
        CMOSSensor(),
        EMCCDSensor(excess_noise_factor=1.0),
        HgCdTeSensor(),
    )
        detector = pe04_frame_detector(sensor=sensor)
        prepared = prepare_detector_acquisition(detector, input)
        capture!(prepared, Xoshiro(23060))
        pe04_assert_owner_storage_contract(detector.state)
        pe04_assert_owner_storage_contract(detector.workspace)
        pe04_assert_owner_storage_contract(detector.products)
    end

    spad = SPADArrayDetector((2, 2); noise=NoiseNone())
    mkid = MKIDArrayDetector(noise=NoiseNone())
    for detector in (spad, mkid)
        pe04_assert_owner_storage_contract(detector.state)
        pe04_assert_owner_storage_contract(detector.workspace)
        pe04_assert_owner_storage_contract(detector.products)
    end
    spad_input = fill(1.0, 2, 2)
    spad_rng = Xoshiro(23061)
    capture!(spad, spad_input, spad_rng)
    @test @inferred(capture!(spad, spad_input, spad_rng)) ===
        output_frame(spad)
    @test_detector_allocation @allocated(
        capture!(spad, spad_input, spad_rng)) == 0
end
