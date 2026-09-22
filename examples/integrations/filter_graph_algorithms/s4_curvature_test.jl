using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms
using JuliaFilterGraph

function replace_s4_curvature_observation(
    observation::S4CurvatureObservation;
    values=observation.values,
    layout=observation.layout,
    calibration_signature=observation.calibration_signature,
    sequence=observation.sequence,
    exposure_duration=observation.exposure_duration,
    model_timestamp_nanoseconds=observation.model_timestamp_nanoseconds,
)
    return S4CurvatureObservation(
        values,
        layout,
        calibration_signature,
        sequence,
        exposure_duration,
        model_timestamp_nanoseconds,
    )
end

function verify_s4_curvature_frame_rejection!(prepared, transform, expected)
    reset_s4_curvature_exchange!(prepared)
    observation = produce_s4_curvature_frame!(prepared)
    retained_signal = copy(prepared.plant_image.signal)
    retained_association = (
        prepared.frame_state.published_sequence,
        prepared.frame_state.published_exposure_duration,
        prepared.frame_state.published_model_timestamp_nanoseconds,
    )
    @test process_s4_curvature_image!(prepared, transform(observation)) ===
        expected
    @test prepared.plant_image.signal == retained_signal
    @test (
        prepared.frame_state.published_sequence,
        prepared.frame_state.published_exposure_duration,
        prepared.frame_state.published_model_timestamp_nanoseconds,
    ) == retained_association
    @test prepared.frame_state.blocked
    @test prepared.frame_state.outstanding
    @test_throws ArgumentError produce_s4_curvature_frame!(prepared)
    return nothing
end

@testset "AOS Curvature plant and FGA 0.5 paired-signal boundary" begin
    prepared = prepare_s4_curvature()
    @test Base.pkgversion(FilterGraphAlgorithms) == v"0.5.0"
    @test Base.pkgversion(JuliaFilterGraph) == v"0.2.3"
    @test prepared.fga_tree_claim ==
        "0ce4ceae4904446fba170ee25324d65b435c06e6"
    @test prepared.jfg_tree_claim ==
        "1893cea2d29f0cd89303939308e23ce7546bc240"

    frozen_frame = frozen_aos_curvature_frame()
    frozen_channels = frozen_aos_curvature_channels()
    expected_fga = Float32[
        10 20 2 7
        30 40 9 5
    ]
    expected_signal = Float32[
        0.6857143, 0.79491526, 0.0, 0.6604651,
    ]
    @test @inferred(transpose_aos_curvature_frame!(
        prepared.frozen_fga_frame, frozen_frame,
    )) === prepared.frozen_fga_frame
    @test prepared.frozen_fga_frame == expected_fga
    @test prepared.image.plan.branch_origins == ((0, 0), (0, 2))
    @test prepared.channels.plan.branch_rows == (0, 1)
    @test prepared.channels.plan.channel_offsets == collect(0:3)
    @test prepared.image.plan.support == Bool[
        true false
        true true
    ]
    @test prepared.image.plan.reference_signal == Float32[
        0.1 0.0
        -0.1 0.2
    ]
    @test prepared.image.plan.branch_scales == (1.25f0, 0.75f0)
    @test prepared.image.plan.calibration_signature == UInt64(0x43555256)

    @test @inferred(Union{Nothing,CurvaturePairedSignalFailure},
        process_s4_curvature_image!(prepared.frozen_fga_frame,
            prepared.image, frozen_frame)) === nothing
    @test prepared.image.signal ≈ expected_signal rtol=2f-7
    @test @inferred(Union{Nothing,CurvaturePairedSignalFailure},
        process_s4_curvature_channels!(prepared.channels,
            frozen_channels)) === nothing
    @test prepared.channels.signal ≈ expected_signal rtol=2f-7
    @test prepared.channels.signal == prepared.image.signal
    reordered_channels = frozen_channels[:, [4, 3, 2, 1]]
    @test process_s4_curvature_channels!(prepared.channels,
        reordered_channels) === nothing
    @test !isapprox(prepared.channels.signal, expected_signal; rtol=2f-7)
    @test process_s4_curvature_channels!(prepared.channels,
        frozen_channels) === nothing
    @test prepared.channels.signal ≈ expected_signal rtol=2f-7

    retained_image = copy(prepared.image.signal)
    prepared.image.signature[1] += UInt64(1)
    @test process_s4_curvature_image!(prepared.frozen_fga_frame,
        prepared.image, frozen_frame) ===
        CurvaturePairedSignalCalibrationMismatch
    @test prepared.image.signal == retained_image
    prepared.image.signature[1] = prepared.image.plan.calibration_signature
    nonfinite_frame = copy(frozen_frame)
    nonfinite_frame[2, 2] = NaN32
    @test process_s4_curvature_image!(prepared.frozen_fga_frame,
        prepared.image, nonfinite_frame) ===
        CurvaturePairedSignalNonFiniteInput
    @test prepared.image.signal == retained_image

    retained_channels = copy(prepared.channels.signal)
    prepared.channels.signature[1] += UInt64(1)
    @test process_s4_curvature_channels!(prepared.channels,
        frozen_channels) === CurvaturePairedSignalCalibrationMismatch
    @test prepared.channels.signal == retained_channels
    prepared.channels.signature[1] =
        prepared.channels.plan.calibration_signature
    nonfinite_channels = copy(frozen_channels)
    nonfinite_channels[2, 4] = NaN32
    @test process_s4_curvature_channels!(prepared.channels,
        nonfinite_channels) === CurvaturePairedSignalNonFiniteInput
    @test prepared.channels.signal == retained_channels

    frame = @inferred produce_s4_curvature_frame!(prepared)
    @test frame.values === observation_storage(prepared.frame_observation)
    @test size(frame.values) == (4, 2)
    @test eltype(frame.values) === Float32
    @test all(isfinite, frame.values)
    @test frame.layout === :curvature_branch_regions
    @test frame.calibration_signature == UInt64(0x43555256)
    @test frame.sequence == UInt64(1)
    @test frame.exposure_duration == 1.0f0
    @test frame.model_timestamp_nanoseconds == Int64(1_000_000_000)
    @test @inferred(process_s4_curvature_image!(prepared, frame)) ===
        S4CurvatureAccepted
    @test all(isfinite, prepared.plant_image.signal)
    @test prepared.frame_state.published_sequence == frame.sequence
    @test prepared.frame_state.published_exposure_duration ==
        frame.exposure_duration
    @test prepared.frame_state.published_model_timestamp_nanoseconds ==
        frame.model_timestamp_nanoseconds

    channels = @inferred produce_s4_curvature_channels!(prepared)
    @test channels.values === observation_storage(prepared.channel_observation)
    @test size(channels.values) == (2, 4)
    @test eltype(channels.values) === Float32
    @test all(isfinite, channels.values)
    @test channels.layout === :curvature_branch_channels
    @test channels.calibration_signature == UInt64(0x43555256)
    @test channels.sequence == UInt64(1)
    @test channels.exposure_duration == 1.0f0
    @test channels.model_timestamp_nanoseconds == Int64(1_000_000_000)
    @test @inferred(process_s4_curvature_channels!(prepared, channels)) ===
        S4CurvatureAccepted
    @test prepared.plant_channels.signal ≈ prepared.plant_image.signal
    @test prepared.channel_state.published_sequence == channels.sequence
    @test prepared.channel_state.published_exposure_duration ==
        channels.exposure_duration
    @test prepared.channel_state.published_model_timestamp_nanoseconds ==
        channels.model_timestamp_nanoseconds

    @test @inferred(step_s4_curvature_frame!(prepared)) == UInt64(2)
    @test @inferred(step_s4_curvature_channels!(prepared)) == UInt64(2)
    @test @allocated(step_s4_curvature_frame!(prepared)) == 0
    @test @allocated(step_s4_curvature_channels!(prepared)) == 0
    raw_frame = observation_storage(prepared.frame_observation)
    raw_channels = observation_storage(prepared.channel_observation)
    @test @allocated(process_s4_curvature_image!(prepared.plant_fga_frame,
        prepared.plant_image, raw_frame)) == 0
    @test @allocated(process_s4_curvature_channels!(
        prepared.plant_channels, raw_channels)) == 0
end

@testset "Curvature observation association rejects before FGA publication" begin
    prepared = prepare_s4_curvature()
    @test step_s4_curvature_frame!(prepared) == UInt64(1)

    verify_s4_curvature_frame_rejection!(prepared, observation ->
        replace_s4_curvature_observation(observation;
            sequence=observation.sequence + UInt64(1)),
        S4CurvatureRejectedSequence)
    verify_s4_curvature_frame_rejection!(prepared, observation ->
        replace_s4_curvature_observation(observation;
            model_timestamp_nanoseconds=
                observation.model_timestamp_nanoseconds + Int64(1)),
        S4CurvatureRejectedModelTimestamp)
    verify_s4_curvature_frame_rejection!(prepared, observation ->
        replace_s4_curvature_observation(observation;
            exposure_duration=observation.exposure_duration / 2),
        S4CurvatureRejectedExposure)
    verify_s4_curvature_frame_rejection!(prepared, observation ->
        replace_s4_curvature_observation(observation;
            layout=:curvature_branch_channels),
        S4CurvatureRejectedLayout)
    verify_s4_curvature_frame_rejection!(prepared, observation ->
        replace_s4_curvature_observation(observation;
            calibration_signature=
                observation.calibration_signature + UInt64(1)),
        S4CurvatureRejectedCalibrationSignature)
    verify_s4_curvature_frame_rejection!(prepared, observation -> begin
        values = copy(observation.values)
        values[1] = NaN32
        replace_s4_curvature_observation(observation; values)
    end, S4CurvatureRejectedEstimator)
end
