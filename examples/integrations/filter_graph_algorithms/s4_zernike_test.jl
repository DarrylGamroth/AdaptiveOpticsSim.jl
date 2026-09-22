using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms
using JuliaFilterGraph

@testset "AOS Zernike plant and FGA 0.4 pupil-signal boundary" begin
    prepared = prepare_s4_zernike()
    @test Base.pkgversion(FilterGraphAlgorithms) == v"0.4.0"
    @test Base.pkgversion(JuliaFilterGraph) == v"0.2.2"
    @test prepared.fga_tree_claim == "598a10402690250190ad1dd976b03577281e6388"
    @test prepared.jfg_tree_claim == "edaa8b14e2a2b9266cff5674ecbe9c9f6c20c17c"

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
