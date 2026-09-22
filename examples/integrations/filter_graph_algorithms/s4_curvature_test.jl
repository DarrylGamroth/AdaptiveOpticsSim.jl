using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms
using JuliaFilterGraph

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
        10 20 2 4
        30 40 6 8
    ]
    expected_signal = Float32[
        0.6857143, 0.8857143, 0.0, 0.5857143,
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
    @test frame === observation_storage(prepared.frame_observation)
    @test size(frame) == (4, 2)
    @test eltype(frame) === Float32
    @test all(isfinite, frame)
    @test @inferred(Union{Nothing,CurvaturePairedSignalFailure},
        process_s4_curvature_image!(prepared.plant_fga_frame,
            prepared.plant_image, frame)) === nothing
    @test all(isfinite, prepared.plant_image.signal)

    channels = @inferred produce_s4_curvature_channels!(prepared)
    @test channels === observation_storage(prepared.channel_observation)
    @test size(channels) == (2, 4)
    @test eltype(channels) === Float32
    @test all(isfinite, channels)
    @test @inferred(Union{Nothing,CurvaturePairedSignalFailure},
        process_s4_curvature_channels!(prepared.plant_channels,
            channels)) === nothing
    @test prepared.plant_channels.signal ≈ prepared.plant_image.signal

    produce_s4_curvature_frame!(prepared)
    process_s4_curvature_image!(prepared.plant_fga_frame,
        prepared.plant_image, frame)
    produce_s4_curvature_channels!(prepared)
    process_s4_curvature_channels!(prepared.plant_channels, channels)
    @test @allocated(produce_s4_curvature_frame!(prepared)) == 0
    @test @allocated(process_s4_curvature_image!(prepared.plant_fga_frame,
        prepared.plant_image, frame)) == 0
    @test @allocated(produce_s4_curvature_channels!(prepared)) == 0
    @test @allocated(process_s4_curvature_channels!(
        prepared.plant_channels, channels)) == 0
end
