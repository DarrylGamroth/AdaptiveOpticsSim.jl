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
    identity = prepared.calibration.identity
    @test frame.calibration_signature == identity.signature
    @test !iszero(identity.signature)
    @test identity.detector_axes == (:x, :y)
    @test identity.estimator_frame_axes == (:row, :column)
    @test identity.subaperture_order == ((0, 0), (0, 4), (4, 0), (4, 4))
    @test identity.slope_pair_order == (:x, :y)
    @test identity.pdm_actuator_order == (1,)
    @test identity.detector_units === :electron_count
    @test identity.slope_units === :pixel
    @test identity.pdm_command_units === :metre
    @test identity.numeric_type === Float32
    @test !iszero(identity.plant_signature)
    @test !iszero(identity.estimator_signature)
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

@testset "S1 calibration identity binds configuration and interpretation" begin
    plant_signature = AOSFGALockstep._plant_signature()
    @test plant_signature !=
        AOSFGALockstep._plant_signature(telescope_resolution=9)
    @test plant_signature !=
        AOSFGALockstep._plant_signature(detector_units=:photon_count)

    graph_path = joinpath(@__DIR__, "s1-shwfs-f32.conf")
    reference_slopes = zeros(Float32, 4, 2)
    reconstructor_matrix = zeros(Float32, 1, 8)
    controller_to_vdm = ones(Float32, 1, 1)
    active_to_full_vdm = ones(Float32, 1, 1)
    vdm_to_pdm = ones(Float32, 1, 1)
    estimator_signature = AOSFGALockstep._estimator_signature(
        graph_path,
        reference_slopes,
        reconstructor_matrix,
        controller_to_vdm,
        active_to_full_vdm,
        vdm_to_pdm,
    )
    changed_reconstructor = copy(reconstructor_matrix)
    changed_reconstructor[1] = 1.0f0
    @test estimator_signature != AOSFGALockstep._estimator_signature(
        graph_path,
        reference_slopes,
        changed_reconstructor,
        controller_to_vdm,
        active_to_full_vdm,
        vdm_to_pdm,
    )

    identity = AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature,
    )
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature + UInt64(1),
        estimator_signature,
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        slope_pair_order=(:y, :x),
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        pdm_command_units=:volt,
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        numeric_type=Float64,
    ).signature
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
