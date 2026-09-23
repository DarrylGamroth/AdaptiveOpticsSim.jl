using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms
using JuliaFilterGraph

@testset "AOS Zernike plant and FGA 0.5 pupil-signal boundary" begin
    prepared = prepare_s4_zernike()
    @test Base.pkgversion(FilterGraphAlgorithms) == v"0.5.2"
    @test Base.pkgversion(JuliaFilterGraph) == v"0.2.3"
    @test prepared.fga_tree_claim == "922a65f3ac9486abeb7baa3e29f08794133b5bc6"
    @test prepared.jfg_tree_claim == "1893cea2d29f0cd89303939308e23ce7546bc240"

    # The asymmetric fixture freezes the AOS (x,y) to FGA
    # (row=y,column=x) transfer, support order, reference image, calibration
    # signature, and both normalization policies.
    frozen_aos = frozen_aos_zernike_frame()
    expected_fga = Float32[
        2 5 11
        3 7 13
        17 19 23
    ]
    @test prepared.mean.plan.support == Bool[
        true false true
        true true false
        false true true
    ]
    @test prepared.mean.plan.reference_pupil_image == Float32[
        0.05 0.0 0.1
        -0.2 0.15 0.0
        0.0 -0.05 0.25
    ]
    @test @inferred(transpose_aos_zernike_frame!(
        prepared.frozen_fga_frame, frozen_aos,
    )) === prepared.frozen_fga_frame
    @test prepared.frozen_fga_frame == expected_fga

    @test @inferred(Union{Nothing,ZernikePupilSignalFailure},
        process_s4_zernike!(prepared.frozen_fga_frame, prepared.mean,
            frozen_aos)) === nothing
    @test prepared.mean.divisor ≈ Float32[10.833333] rtol=2f-7
    @test prepared.mean.signal ≈ Float32[
        0.13461539, 0.47692308, 0.49615386,
        1.8038461, 0.9153846, 1.8730769,
    ] rtol=2f-7

    # The calibration signature travels with the acquired frame. Expected FGA
    # rejections must not publish a partial signal or normalization divisor.
    @test prepared.mean.signature == UInt64[prepared.mean.plan.calibration_signature]
    retained_signal = copy(prepared.mean.signal)
    retained_divisor = copy(prepared.mean.divisor)
    prepared.mean.signature[1] += UInt64(1)
    @test @inferred(Union{Nothing,ZernikePupilSignalFailure},
        process_s4_zernike!(prepared.frozen_fga_frame, prepared.mean,
            frozen_aos)) === ZernikePupilSignalCalibrationMismatch
    @test prepared.mean.signal == retained_signal
    @test prepared.mean.divisor == retained_divisor
    prepared.mean.signature[1] = prepared.mean.plan.calibration_signature

    nonfinite_aos = copy(frozen_aos)
    nonfinite_aos[2, 3] = NaN32
    @test @inferred(Union{Nothing,ZernikePupilSignalFailure},
        process_s4_zernike!(prepared.frozen_fga_frame, prepared.mean,
            nonfinite_aos)) === ZernikePupilSignalNonFiniteInput
    @test prepared.mean.signal == retained_signal
    @test prepared.mean.divisor == retained_divisor

    @test @inferred(Union{Nothing,ZernikePupilSignalFailure},
        process_s4_zernike!(prepared.frozen_fga_frame, prepared.incidence,
            frozen_aos)) === nothing
    @test prepared.incidence.divisor == Float32[13.5]
    @test prepared.incidence.signal ≈ Float32[
        0.09814815, 0.42222223, 0.36851853,
        1.4574074, 0.71481484, 1.4537038,
    ] rtol=2f-7

    frame = @inferred produce_s4_zernike_frame!(prepared)
    @test frame === observation_storage(prepared.observation)
    @test size(frame) == (3, 3)
    @test eltype(frame) === Float32
    @test all(isfinite, frame)
    @test @inferred(Union{Nothing,ZernikePupilSignalFailure},
        process_s4_zernike!(prepared.plant_fga_frame,
            prepared.plant_signal, frame)) === nothing
    @test all(isfinite, prepared.plant_signal.signal)
    @test prepared.plant_signal.divisor[1] > 0

    produce_s4_zernike_frame!(prepared)
    process_s4_zernike!(prepared.plant_fga_frame,
        prepared.plant_signal, frame)
    @test @allocated(produce_s4_zernike_frame!(prepared)) == 0
    @test @allocated(process_s4_zernike!(prepared.plant_fga_frame,
        prepared.plant_signal, frame)) == 0
end
