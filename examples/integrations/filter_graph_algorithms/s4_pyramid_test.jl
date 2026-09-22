using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage

function replace_s4_frame(
    frame::S4DetectorFrame;
    values=frame.values,
    schema=frame.schema,
    identity_signature=frame.identity_signature,
    calibration_signature=frame.calibration_signature,
    order_signature=frame.order_signature,
    sequence=frame.sequence,
    timestamp_nanoseconds=frame.timestamp_nanoseconds,
    terminal=frame.terminal,
    discontinuity=frame.discontinuity,
    corrupted=frame.corrupted,
)
    return S4DetectorFrame(
        values,
        schema,
        identity_signature,
        calibration_signature,
        order_signature,
        sequence,
        timestamp_nanoseconds,
        terminal,
        discontinuity,
        corrupted,
    )
end

function verify_s4_rejection!(prepared, transform, expected)
    reset_s4_pyramid!(prepared)
    frame = produce_s4_detector_frame!(prepared)
    before = copy(prepared.adopted_command)
    @test process_s4_detector_frame!(prepared, transform(frame)) === expected
    @test prepared.adopted_command == before
    @test prepared.state.blocked
    @test prepared.state.outstanding
    @test_throws ArgumentError produce_s4_detector_frame!(prepared)
    return nothing
end

@testset "AOS physical Pyramid plant and FGA S4 estimator boundary" begin
    prepared = prepare_s4_pyramid()
    identity = prepared.identity
    @test identity.detector_axes == (:x, :y)
    @test identity.estimator_frame_axes == (:row, :column)
    @test identity.pupil_origins == ((0, 0), (2, 0), (2, 2), (0, 2))
    @test identity.pupil_order ==
          (:q1_top_left, :q2_bottom_left, :q3_bottom_right, :q4_top_right)
    @test identity.support_order == ((1, 1), (2, 1), (1, 2), (2, 2))
    @test identity.component_order == (:x, :y)
    @test identity.normalization === :mean_valid_flux
    @test identity.numeric_type === Float32
    @test identity.detector_units === :electron_count
    @test identity.signal_units === :normalized_pyramid_i4q
    @test identity.command_units === :metre
    @test identity.fga_version == v"0.4.0"
    @test identity.fga_tree_claim ==
          "598a10402690250190ad1dd976b03577281e6388"
    @test !iszero(identity.signature)
    @test !iszero(identity.calibration_signature)
    @test !iszero(identity.order_signature)

    frame = @inferred produce_s4_detector_frame!(prepared)
    @test frame.values === observation_storage(prepared.observation)
    @test eltype(frame.values) === Float32
    @test size(frame.values) == (4, 4)
    @test frame.schema == S4_DETECTOR_FRAME_SCHEMA
    @test frame.identity_signature == identity.signature
    @test frame.calibration_signature == identity.calibration_signature
    @test frame.order_signature == identity.order_signature
    @test frame.sequence == UInt64(1)
    @test frame.timestamp_nanoseconds == 1_000_000
    @test frame.terminal
    @test all(isfinite, frame.values)

    oracle_signal = zeros(Float32, 4, 2)
    oracle_divisor = zeros(Float32, 1)
    reference = prepared.image.plan.reference_signal
    gain = prepared.image.plan.optical_gain
    transposed = permutedims(frame.values, (2, 1))
    direct_s4_signal!(oracle_signal, oracle_divisor, transposed, reference, gain)
    @test @inferred(process_s4_detector_frame!(prepared, frame)) === S4FrameAccepted
    @test prepared.signal ≈ oracle_signal rtol=2f-6
    @test prepared.divisor ≈ oracle_divisor rtol=2f-6
    oracle_residual = 0.25f0 * sum(@view oracle_signal[:, 1])
    @test prepared.residual ≈ Float32[oracle_residual] rtol=2f-6
    @test prepared.candidate_command ≈ Float32[-oracle_residual] rtol=2f-6
    @test prepared.adopted_command == prepared.candidate_command
    @test !all(iszero, prepared.adopted_command)

    reset_s4_pyramid!(prepared)
    @test @inferred(step_s4_pyramid!(prepared)) == UInt64(1)
    @test @inferred(step_s4_pyramid!(prepared)) == UInt64(2)
    @test @allocated(step_s4_pyramid!(prepared)) == 0
end

@testset "S4 rejects invalid boundary frames before command adoption" begin
    prepared = prepare_s4_pyramid()
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; sequence=frame.sequence + UInt64(1)),
        S4FrameRejectedSequence,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(
            frame;
            timestamp_nanoseconds=frame.timestamp_nanoseconds + 1,
        ),
        S4FrameRejectedTimestamp,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; terminal=false),
        S4FrameRejectedIncomplete,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; discontinuity=true),
        S4FrameRejectedDiscontinuous,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; corrupted=true),
        S4FrameRejectedCorrupted,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; schema="wrong-schema/1"),
        S4FrameRejectedSchema,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(
            frame;
            identity_signature=frame.identity_signature + UInt64(1),
        ),
        S4FrameRejectedIdentity,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(
            frame;
            calibration_signature=frame.calibration_signature + UInt64(1),
        ),
        S4FrameRejectedCalibrationSignature,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(
            frame;
            order_signature=frame.order_signature + UInt64(1),
        ),
        S4FrameRejectedOrder,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; values=Float64.(frame.values)),
        S4FrameRejectedNumericType,
    )
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; values=zeros(Float32, 2, 8)),
        S4FrameRejectedShape,
    )
    verify_s4_rejection!(prepared, frame -> begin
        values = copy(frame.values)
        values[1] = NaN32
        replace_s4_frame(frame; values)
    end, S4FrameRejectedNonFinite)
    # A validly typed all-zero detector frame reaches FGA and is rejected by
    # its mean-valid-flux normalization before any command is adopted.
    verify_s4_rejection!(
        prepared,
        frame -> replace_s4_frame(frame; values=zeros(Float32, 4, 4)),
        S4FrameRejectedFGA,
    )
end

@testset "S4 reset discards blocked exchange and FGA control state" begin
    prepared = prepare_s4_pyramid()
    frame = produce_s4_detector_frame!(prepared)
    @test process_s4_detector_frame!(
        prepared,
        replace_s4_frame(frame; discontinuity=true),
    ) === S4FrameRejectedDiscontinuous
    reset_s4_pyramid!(prepared)
    @test !prepared.state.blocked
    @test !prepared.state.outstanding
    @test prepared.state.sequence == 0
    @test all(iszero, prepared.adopted_command)
    @test step_s4_pyramid!(prepared) == UInt64(1)
end
