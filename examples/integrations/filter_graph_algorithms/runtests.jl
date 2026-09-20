using Test

include("s1_lockstep.jl")
using .AOSFGALockstep
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms: SampleMetadata
using JuliaFilterGraph
using LinearAlgebra

function replace_frame(
    frame::S1DetectorFrame;
    values=frame.values,
    schema=frame.schema,
    calibration_signature=frame.calibration_signature,
    sequence=frame.sequence,
    timestamp_nanoseconds=frame.timestamp_nanoseconds,
    terminal=frame.terminal,
    discontinuity=frame.discontinuity,
    corrupted=frame.corrupted,
)
    return S1DetectorFrame(
        values,
        schema,
        calibration_signature,
        sequence,
        timestamp_nanoseconds,
        terminal,
        discontinuity,
        corrupted,
    )
end

function verify_rejection!(prepared, transform, expected)
    reset_s1_lockstep!(prepared)
    frame = produce_detector_frame!(prepared)
    before = copy(prepared.adopted_command)
    rejected = transform(frame)
    @test process_detector_frame!(prepared, rejected) === expected
    @test prepared.adopted_command == before
    @test prepared.state.blocked
    @test prepared.state.outstanding
    @test_throws ArgumentError produce_detector_frame!(prepared)
    return nothing
end

@testset "AOS plant and FGA RTC S1 lockstep" begin
    prepared = prepare_s1_lockstep()

    frame = @inferred produce_detector_frame!(prepared)
    @test frame.values === observation_storage(prepared.observation)
    @test eltype(frame.values) === Float32
    @test size(frame.values) == (8, 8)
    @test frame.schema == DETECTOR_FRAME_SCHEMA
    @test frame.calibration_signature ==
        DETECTOR_FRAME_CALIBRATION_SIGNATURE
    @test frame.sequence == UInt64(1)
    @test frame.timestamp_nanoseconds == 1_000_000
    @test frame.terminal
    @test all(isfinite, frame.values)
    @test all(iszero, prepared.applied_command)

    @test @inferred(process_detector_frame!(prepared, frame)) === FrameAccepted
    @test !prepared.state.outstanding
    @test !all(iszero, prepared.adopted_command)
    first_command = copy(prepared.adopted_command)
    metadata = JuliaFilterGraph.output_metadata(prepared.graph).demanded
    @test metadata isa SampleMetadata
    @test metadata.sequence == frame.sequence
    @test metadata.has_timestamp
    @test metadata.timestamp_nanoseconds == frame.timestamp_nanoseconds

    @test @inferred(step_lockstep!(prepared)) == UInt64(2)
    @test prepared.applied_command == first_command
    @test @allocated(step_lockstep!(prepared)) == 0

    reset_s1_lockstep!(prepared)
    residual_norms = Float32[]
    for sequence in UInt64(1):UInt64(9)
        @test step_lockstep!(prepared) == sequence
        push!(residual_norms, norm(prepared.outputs.slopes))
    end
    @test all(<(0), diff(residual_norms))
    @test last(residual_norms) < first(residual_norms) * 0.02f0
    @test prepared.adopted_command[1] ≈ -3.0f-8 rtol=0.02f0
end

@testset "S1 rejects invalid frames before command adoption" begin
    prepared = prepare_s1_lockstep()

    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; sequence=frame.sequence + UInt64(1)),
        FrameRejectedSequence,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(
            frame;
            timestamp_nanoseconds=frame.timestamp_nanoseconds + 1,
        ),
        FrameRejectedTimestamp,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; terminal=false),
        FrameRejectedIncomplete,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; discontinuity=true),
        FrameRejectedDiscontinuous,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; corrupted=true),
        FrameRejectedCorrupted,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; schema="test.wrong-frame-schema/1"),
        FrameRejectedSchema,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(
            frame;
            calibration_signature=frame.calibration_signature + UInt64(1),
        ),
        FrameRejectedCalibrationSignature,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; values=Float64.(frame.values)),
        FrameRejectedNumericType,
    )
    verify_rejection!(
        prepared,
        frame -> replace_frame(frame; values=zeros(Float32, 7, 8)),
        FrameRejectedShape,
    )
    verify_rejection!(
        prepared,
        frame -> begin
            values = copy(frame.values)
            values[1] = Float32(NaN)
            replace_frame(frame; values)
        end,
        FrameRejectedNonFinite,
    )
end
