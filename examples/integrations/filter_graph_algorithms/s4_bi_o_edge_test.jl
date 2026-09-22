using Test
using AdaptiveOpticsSim.WavefrontSensors: observation_storage
using FilterGraphAlgorithms
using JuliaFilterGraph

@testset "AOS Bi-O-edge plant and FGA 0.3 image-estimator boundary" begin
    prepared = prepare_s4_bi_o_edge()
    @test Base.pkgversion(FilterGraphAlgorithms) == v"0.3.0"
    @test Base.pkgversion(JuliaFilterGraph) == v"0.2.1"

    # The frozen fixture is deliberately asymmetric. It establishes the AOS
    # (x,y) -> FGA (row=y,column=x) transfer, q1/q2/q3/q4 order TL/BL/BR/TR,
    # zero-based FGA origins, non-rectangular support, nonzero reference, and
    # unequal gains independently of the physical plant.
    frozen_aos = frozen_aos_bi_o_edge_frame()
    expected_fga = Float32[
        10 30 5 7
        20 40 6 8
        1 3 9 11
        2 4 10 12
    ]
    @test prepared.mean.plan.origins == ((0, 0), (2, 0), (2, 2), (0, 2))
    @test prepared.mean.plan.support == Bool[true true; false true]
    @test prepared.mean.plan.reference_signal == Float32[0.1 -0.1; 0.3 -0.3; 0.4 -0.4]
    @test prepared.mean.plan.optical_gain == Float32[1.5 0.25; 0.5 1.25; 2.0 0.75]
    @test @inferred(transpose_aos_bi_o_edge_frame!(prepared.frozen_fga_frame, frozen_aos)) ===
          prepared.frozen_fga_frame
    @test prepared.frozen_fga_frame == expected_fga

    @test @inferred(Union{Nothing, BiOEdgeImageFailure}, process_s4_bi_o_edge!(
        prepared.frozen_fga_frame, prepared.mean, frozen_aos,
    )) === nothing
    @test prepared.mean.divisor ≈ Float32[46.666668] rtol=2f-7
    @test vec(prepared.mean.signal) ≈ Float32[
        0.010714274, 0.09642856, 0.5714285,
        0.008928573, 0.77678573, 0.68571424,
    ] rtol=2f-7

    @test @inferred(Union{Nothing, BiOEdgeImageFailure}, process_s4_bi_o_edge!(
        prepared.frozen_fga_frame, prepared.incidence, frozen_aos,
    )) === nothing
    @test prepared.incidence.divisor == Float32[50]
    @test vec(prepared.incidence.signal) ≈ Float32[
        0.0, 0.08, 0.47999996,
        0.010000001, 0.75, 0.65999997,
    ] rtol=2f-7

    frame = @inferred produce_s4_bi_o_edge_frame!(prepared)
    @test frame === observation_storage(prepared.observation)
    @test size(frame) == (8, 8)
    @test eltype(frame) === Float32
    @test all(isfinite, frame)
    @test @inferred(Union{Nothing, BiOEdgeImageFailure}, process_s4_bi_o_edge!(
        prepared.plant_fga_frame, prepared.plant_image, frame,
    )) === nothing
    @test all(isfinite, prepared.plant_image.signal)
    @test prepared.plant_image.divisor[1] > 0

    # Warmed physical plant formation, preallocated transpose, and FGA image
    # estimation retain the CPU zero-allocation contract independently.
    produce_s4_bi_o_edge_frame!(prepared)
    process_s4_bi_o_edge!(prepared.plant_fga_frame, prepared.plant_image, frame)
    @test @allocated(produce_s4_bi_o_edge_frame!(prepared)) == 0
    @test @allocated(process_s4_bi_o_edge!(
        prepared.plant_fga_frame, prepared.plant_image, frame,
    )) == 0
end
