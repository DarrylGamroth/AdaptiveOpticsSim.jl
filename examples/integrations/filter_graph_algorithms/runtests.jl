using Test

include("s1_lockstep.jl")
include("rtc_parity.jl")
using .AOSFGALockstep
using AdaptiveOpticsCalibration.Reconstructors: reconstructor
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms: SampleMetadata
using JuliaFilterGraph
using LinearAlgebra

# Bit-exact cold-calibration products captured from pre-S2 commit d0cfead.
const S1_EXPECTED_INTERACTION = reshape(
    Float32[
        1.7699274f6,
        1.7699275f6,
        -1.7699278f6,
        1.769927f6,
        1.7699278f6,
        -1.769929f6,
        -1.7699292f6,
        -1.7699292f6,
    ],
    8,
    1,
)

const S1_EXPECTED_RECONSTRUCTOR = reshape(
    Float32[
        7.0624324f-8,
        7.062432f-8,
        -7.062433f-8,
        7.0624296f-8,
        7.062433f-8,
        -7.062438f-8,
        -7.062439f-8,
        -7.062439f-8,
    ],
    1,
    8,
)

# Deterministic accepted-frame trajectory from the AOS plant / FGA RTC boundary
# at AOS commit d6a683008d9e40ca5bd90e5ccf748e5aa3384caa.  Residual norms are in
# controller reconstruction coordinates; commands are PDM metres.
const S1_EXPECTED_RESIDUAL_NORMS = Float32[
    0.15002562,
    0.09018742,
    0.05413042,
    0.032470427,
    0.019473683,
    0.011678018,
    0.007003085,
    0.0041994397,
    0.002518242,
]
const S1_EXPECTED_DEMANDED_COMMANDS = Float32[
    -1.1987396e-8,
    -1.9193582e-8,
    -2.3518728e-8,
    -2.6113192e-8,
    -2.7669184e-8,
    -2.8602285e-8,
    -2.9161848e-8,
    -2.9497393e-8,
    -2.9698606e-8,
]

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

    @test prepared.calibration.interaction_matrix == S1_EXPECTED_INTERACTION
    @test reconstructor(prepared.calibration.reconstructor_product) ==
        S1_EXPECTED_RECONSTRUCTOR

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
    @test identity.subaperture_order == Tuple(
        prepared.graph.nodes[1].prepared.plan.regions.origins,
    )
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
    for index in eachindex(S1_EXPECTED_RESIDUAL_NORMS)
        sequence = UInt64(index)
        @test step_lockstep!(prepared) == sequence
        @test norm(prepared.outputs.slopes) ≈ S1_EXPECTED_RESIDUAL_NORMS[index] rtol=1f-6
        @test prepared.outputs.demanded[1] ≈
              S1_EXPECTED_DEMANDED_COMMANDS[index] rtol=1f-6
        @test prepared.adopted_command[1] == prepared.outputs.demanded[1]
        expected_applied = if isone(index)
            0.0f0
        else
            S1_EXPECTED_DEMANDED_COMMANDS[index - 1]
        end
        @test prepared.applied_command[1] ≈ expected_applied rtol=1f-6
    end
    @test all(<(0), diff(S1_EXPECTED_RESIDUAL_NORMS))
    @test last(S1_EXPECTED_RESIDUAL_NORMS) <
          first(S1_EXPECTED_RESIDUAL_NORMS) * 0.02f0
    @test prepared.adopted_command[1] ≈ -3.0f-8 rtol=0.02f0
end

@testset "S1 calibration identity binds configuration and interpretation" begin
    plant_signature = AOSFGALockstep._plant_signature()
    @test plant_signature !=
        AOSFGALockstep._plant_signature(telescope_resolution=9)
    @test plant_signature !=
        AOSFGALockstep._plant_signature(detector_units=:photon_count)
    @test plant_signature !=
        AOSFGALockstep._plant_signature(telescope_fov_arcsec=1.0f0)
    @test plant_signature !=
        AOSFGALockstep._plant_signature(telescope_reflectivity=0.9f0)
    @test plant_signature !=
        AOSFGALockstep._plant_signature(source_band=:I)
    @test plant_signature != AOSFGALockstep._plant_signature(
        source_coordinates_arcsec_deg=(1.0f0, 0.0f0),
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        source_radiometry=AOSFGALockstep.NormalizedTestSource(),
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        actuator_coordinates=((0.0f0, 0.0f0), (0.1f0, 0.0f0)),
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_valid_threshold=0.2f0,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_cog_threshold=0.02f0,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_convolution_threshold=0.1f0,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_diffraction_padding=3,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_half_pixel_shift=true,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_pixel_scale=0.1f0,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        shack_hartmann_shannon_sampling=false,
    )
    @test plant_signature != AOSFGALockstep._plant_signature(
        observation_layout=:detector_frame,
    )
    @test AOSFGALockstep._plant_signature(
        detector_metadata=(gain=1.0f0,),
    ) != AOSFGALockstep._plant_signature(
        detector_metadata=(gain=2.0f0,),
    )

    graph_path = joinpath(@__DIR__, "s1-shwfs-f32.conf")
    reference_slopes = zeros(Float32, 4, 2)
    reconstructor_matrix = zeros(Float32, 1, 8)
    controller_to_vdm = ones(Float32, 1, 1)
    active_to_full_vdm = ones(Float32, 1, 1)
    vdm_to_pdm = ones(Float32, 1, 1)
    subaperture_order = ((0, 0), (0, 4), (4, 0), (4, 4))
    estimator_signature = AOSFGALockstep._estimator_signature(
        graph_path,
        reference_slopes,
        reconstructor_matrix,
        controller_to_vdm,
        active_to_full_vdm,
        vdm_to_pdm,
        subaperture_order,
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
        subaperture_order,
    )
    @test estimator_signature != AOSFGALockstep._estimator_signature(
        graph_path,
        reference_slopes,
        reconstructor_matrix,
        controller_to_vdm,
        active_to_full_vdm,
        vdm_to_pdm,
        reverse(subaperture_order),
    )

    identity = AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature,
        subaperture_order=subaperture_order,
    )
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature + UInt64(1),
        estimator_signature,
        subaperture_order=subaperture_order,
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        subaperture_order=subaperture_order,
        slope_pair_order=(:y, :x),
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        subaperture_order=subaperture_order,
        pdm_command_units=:volt,
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        subaperture_order=subaperture_order,
        numeric_type=Float64,
    ).signature
    @test identity.signature != AOSFGALockstep._calibration_identity(
        plant_signature,
        estimator_signature;
        subaperture_order=reverse(subaperture_order),
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
