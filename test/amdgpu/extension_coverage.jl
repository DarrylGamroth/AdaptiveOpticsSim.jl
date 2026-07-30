using AMDGPU

include(joinpath(@__DIR__, "..", "backend_extension_coverage_common.jl"))

AMDGPU.functional() ||
    error("AMDGPU extension coverage requires a functional AMDGPU device")

@testset "AMDGPU extension coverage" begin
    run_backend_extension_coverage(Backends.AMDGPUBackendTag)

    original = ComplexF32.(reshape(1:16, 4, 4))
    fft_buffer = AMDGPU.ROCArray(original)
    forward_plan = Backends.plan_fft_backend!(fft_buffer)
    inverse_plan = Backends.plan_ifft_backend!(fft_buffer)
    Backends.execute_fft_plan!(fft_buffer, forward_plan)
    Backends.execute_fft_plan!(fft_buffer, inverse_plan)
    @test Array(fft_buffer) ≈ original rtol=2f-5 atol=2f-5

    build_backend = Calibration.default_build_backend(fft_buffer)
    build_input = Float32[1 2; 3 4]
    @test Calibration.prepare_build_matrix(
        build_backend, build_input) == build_input
    @test Backends.reduction_execution_plan(
        fft_buffer) isa Backends.HostMirrorReductionPlan

    style = Backends.execution_style(fft_buffer)
    noise = AMDGPU.zeros(Float32, 8)
    Backends.randn_backend_async!(
        style, MersenneTwister(31), noise)
    Backends._randn_backend!(
        style, MersenneTwister(32), noise)
    @test all(isfinite, Array(noise))

    detector = AdaptiveOpticsSim.Detectors.Detector(
        noise=AdaptiveOpticsSim.Detectors.NoiseNone(),
        T=Float32,
        backend=Backends.AMDGPUBackend(),
    )
    frame_noise = AMDGPU.zeros(Float32, 2, 2)
    cube_noise = AMDGPU.zeros(Float32, 2, 2, 2)
    host_plan = AdaptiveOpticsSim.Detectors.DetectorHostMirrorPlan()
    AdaptiveOpticsSim.Detectors._randn_frame_noise!(
        host_plan, detector, MersenneTwister(33), frame_noise)
    AdaptiveOpticsSim.Detectors._randn_frame_noise!(
        host_plan, detector, MersenneTwister(34), cube_noise)
    @test all(isfinite, Array(frame_noise))
    @test all(isfinite, Array(cube_noise))

    frequencies = AMDGPU.ROCArray(Float32[-0.2, -0.1, 0.1, 0.2])
    phase_psd = AMDGPU.zeros(Float32, 4, 4)
    AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(
        Backends.execution_style(phase_psd),
        phase_psd,
        frequencies,
        0.02f0,
        Float32(4pi^2),
        0.01f0,
        Float32(-11 / 6),
        0.0f0,
        4,
    )
    @test all(isfinite, Array(phase_psd))

    opd = AMDGPU.ROCArray(reshape(Float32.(1:16), 4, 4))
    valid_mask = AMDGPU.ROCArray(trues(2, 2))
    edge_mask = AMDGPU.ROCArray(trues(4, 4))
    geometric_slopes = AMDGPU.zeros(Float32, 8)
    edge_slopes = AMDGPU.zeros(Float32, 8)
    AdaptiveOpticsSim.WavefrontSensors._geometric_slopes!(
        Backends.execution_style(geometric_slopes),
        geometric_slopes,
        opd,
        valid_mask,
        2,
        2,
        4,
    )
    AdaptiveOpticsSim.WavefrontSensors._edge_geometric_slopes!(
        Backends.execution_style(edge_slopes),
        edge_slopes,
        opd,
        valid_mask,
        edge_mask,
        2,
        2,
        4,
    )
    @test all(isfinite, Array(geometric_slopes))
    @test all(isfinite, Array(edge_slopes))

    left_host = Float32[1 2; 3 4]
    right_host = Float32[5 6; 7 8]
    left = AMDGPU.ROCArray(left_host)
    right = AMDGPU.ROCArray(right_host)
    @test Array(Backends.backend_matmul(left, right)) ≈
        left_host * right_host
    @test Array(Backends.backend_matmul_transpose_right(left, right)) ≈
        left_host * transpose(right_host)

    inverse_host = Float32[2 0; 0 1]
    inverse_input = AMDGPU.ROCArray(inverse_host)
    for policy in (
        Calibration.ExactPseudoInverse(),
        Calibration.TSVDInverse(rtol=1.0f-6),
        Calibration.TikhonovInverse(0.1f0),
    )
        inverse, stats = Calibration.inverse_operator(
            build_backend, inverse_input, policy)
        @test size(inverse) == reverse(size(inverse_host))
        @test all(isfinite, Array(inverse))
        @test stats.effective_rank == 2
    end

    gram_host = Float32[2 0; 0 4]
    division = AdaptiveOpticsSim.Tomography.stable_hermitian_right_division(
        build_backend,
        AMDGPU.ROCArray(left_host),
        AMDGPU.ROCArray(gram_host),
    )
    @test Array(division) ≈ left_host / gram_host rtol=1.0f-5 atol=1.0f-6
end
