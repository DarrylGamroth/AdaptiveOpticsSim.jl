function lift_hot_path_allocation_bytes!(H, coefficients, lift,
    observation)
    WavefrontSensors.lift_interaction_matrix!(H, lift, coefficients)
    interaction_bytes = @allocated WavefrontSensors.lift_interaction_matrix!(
        H, lift, coefficients)
    WavefrontSensors.reconstruct!(coefficients, lift, observation;
        optimize_norm=:none, check_convergence=false)
    reconstruction_bytes = @allocated WavefrontSensors.reconstruct!(
        coefficients, lift, observation; optimize_norm=:none,
        check_convergence=false)
    return interaction_bytes, reconstruction_bytes
end

@testset "LiFT" begin
    fixed_damping = LiFTLevenbergMarquardt(lambda0=Float32(1e-4))
    @test fixed_damping isa LiFTLevenbergMarquardt{Float64}
    @test fixed_damping.lambda0 == Float64(Float32(1e-4))
    @test LiFTLevenbergMarquardt(Float32(1e-4), 10.0, Float32(1e-3)) isa
        LiFTLevenbergMarquardt{Float64}
    @test LiFTLevenbergMarquardt(lambda0=1f-4, growth=10f0,
        condition_rtol=1f-3) isa LiFTLevenbergMarquardt{Float32}

    adaptive_damping = LiFTAdaptiveLevenbergMarquardt(
        lambda0=Float32(1e-4), min_lambda=Float32(1e-8))
    @test adaptive_damping isa LiFTAdaptiveLevenbergMarquardt{Float64}
    @test adaptive_damping.lambda0 == Float64(Float32(1e-4))
    @test adaptive_damping.min_lambda == Float64(Float32(1e-8))
    @test WavefrontSensors.effective_damping(adaptive_damping, 1f-4) isa
        LiFTLevenbergMarquardt{Float64}

    fallback_H = [1.0 0.0; 0.0 2.0; 1.0 1.0]
    fallback_residual = [2.0, -1.0, 0.5]
    fallback_policy = LiFTLevenbergMarquardt(
        lambda0=0.1, growth=10.0, condition_rtol=1e-6)
    fallback_lambda = WavefrontSensors.damping_lambda(
        fallback_policy, transpose(fallback_H) * fallback_H)
    expected_fallback = (transpose(fallback_H) * fallback_H +
        fallback_lambda * I) \ (transpose(fallback_H) * fallback_residual)
    fallback_rhs = zeros(2)
    fallback_diagnostics = WavefrontSensors.LiFTDiagnostics(
        NaN, NaN, NaN, NaN, 0.0, false, false)
    WavefrontSensors.solve_lift_fallback!(fallback_diagnostics,
        fallback_rhs, fallback_H, fallback_residual, fallback_policy)
    @test fallback_rhs ≈ expected_fallback rtol=1e-12 atol=1e-12
    @test fallback_diagnostics.regularization == fallback_lambda
    @test fallback_diagnostics.used_fallback

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(8, 8, 3)
    diversity = zeros(8, 8)
    forward = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8)
    lift = LiFT(forward; iterations=2, mode_ids=(1, 2), numerical=true)
    H = lift_interaction_matrix(lift, zeros(3))
    @test size(H) == (64, 2)
    @test maximum(abs, H) > 0
    kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    forward_analytic = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8, object_kernel=kernel)
    lift_analytic = LiFT(forward_analytic; iterations=2,
        mode_ids=(1, 2), numerical=false)
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3))
    @test size(H_analytic) == (64, 2)
    psf = reference_direct_image(tel, src; zero_padding=1)
    observation = LiFTObservation(forward, copy(psf))
    coeffs = WavefrontSensors.reconstruct(lift, observation)
    @test length(coeffs) == 2
    @test fieldnames(typeof(lift)) == (:forward, :params, :state)
    @test !hasfield(typeof(lift), :detector)
    @test !hasfield(typeof(lift), :telescope)
    @test observation.values !== intensity_values(lift_forward_output(forward))
    @test lift.params.mode_ids == (1, 2)
    @test_throws InvalidConfiguration LiFT(forward; mode_ids=())
    @test_throws InvalidConfiguration LiFT(forward; mode_ids=(1, 1))
    @test_throws DimensionMismatchError LiFT(forward; mode_ids=(1, 4))
    @test_throws MethodError LiFT(tel, src, basis, det;
        diversity_opd=diversity)
    @test_throws MethodError WavefrontSensors.reconstruct(lift, psf)
    @test WavefrontSensors.effective_solve_mode(
        AdaptiveOpticsSim.Backends.ScalarCPUStyle(), LiFTSolveAuto()) isa
        LiFTSolveQR
    diag = diagnostics(lift)
    @test diag.used_qr isa Bool
    @test isfinite(diag.residual_norm)
    @test isfinite(diag.weighted_residual_norm)
    @test isfinite(diag.update_norm)
    @test isfinite(diag.condition_ratio) || isinf(diag.condition_ratio)
    lift_normal = LiFT(forward; iterations=2, numerical=true,
        mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    coeffs_normal = WavefrontSensors.reconstruct(lift_normal, observation)
    @test length(coeffs_normal) == 2
    @test all(isfinite, coeffs_normal)
    @test !diagnostics(lift_normal).used_qr

    allocation_values = copy(intensity_values(evaluate_lift_forward!(
        forward, diversity)))
    allocation_observation = LiFTObservation(forward, allocation_values)
    allocation_lift = LiFT(forward; iterations=1, mode_ids=(1, 2),
        numerical=false, solve_mode=LiFTSolveNormalEquations())
    allocation_coefficients = zeros(3)
    allocation_H = allocation_lift.state.H_buffer
    @test (@inferred WavefrontSensors.reconstruct!(
        allocation_coefficients, allocation_lift, allocation_observation;
        optimize_norm=:none, check_convergence=false)) ===
        allocation_coefficients
    if coverage_instrumented()
        @test_skip "LiFT allocation assertions are disabled under coverage instrumentation"
    else
        interaction_bytes, reconstruction_bytes =
            lift_hot_path_allocation_bytes!(allocation_H,
                allocation_coefficients, allocation_lift,
                allocation_observation)
        @test interaction_bytes == 0
        @test reconstruction_bytes == 0
    end

    tel32 = Telescope(resolution=8, diameter=8f0,
        central_obstruction=0f0, T=Float32)
    src32 = Source(band=:I, magnitude=0f0, T=Float32)
    det32 = Detector(noise=NoiseNone(), psf_sampling=1, T=Float32)
    basis32 = rand(MersenneTwister(29), Float32, 8, 8, 3)
    forward32 = prepare_lift_forward_model(tel32, src32, basis32;
        diversity_opd=zeros(Float32, 8, 8), focal_resolution=8)
    lift_normal32 = LiFT(forward32; iterations=2,
        mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    psf32 = reference_direct_image(tel32, src32; zero_padding=1)
    observation32 = LiFTObservation(forward32, copy(psf32))
    coeffs_normal32 = WavefrontSensors.reconstruct(
        lift_normal32, observation32)
    @test all(isfinite, coeffs_normal32)
    @test all(isfinite, @view(lift_normal32.state.normal_buffer[1:2, 1:2]))

    lift_damped = LiFT(forward; iterations=2, numerical=true,
        mode_ids=(1, 2),
        damping=LiFTLevenbergMarquardt())
    coeffs_damped = WavefrontSensors.reconstruct(lift_damped, observation)
    @test length(coeffs_damped) == 2
    @test all(isfinite, coeffs_damped)
    @test diagnostics(lift_damped).regularization >= 0
    lift_adaptive = LiFT(forward; iterations=2, numerical=true,
        mode_ids=(1, 2),
        damping=LiFTAdaptiveLevenbergMarquardt())
    coeffs_adaptive = WavefrontSensors.reconstruct(lift_adaptive, observation)
    @test length(coeffs_adaptive) == 2
    @test all(isfinite, coeffs_adaptive)
    @test diagnostics(lift_adaptive).regularization >= 0

    zernike = ZernikeBasis(tel, 4)
    compute_zernike!(zernike, tel)
    adaptive_basis = copy(@view zernike.modes[:, :, 2:3])
    adaptive_diversity = 50e-9 .* @view(zernike.modes[:, :, 4])
    adaptive_truth = [10e-9, -5e-9]
    adaptive_truth_opd = adaptive_diversity .+
        adaptive_truth[1] .* @view(adaptive_basis[:, :, 1]) .+
        adaptive_truth[2] .* @view(adaptive_basis[:, :, 2])
    adaptive_forward = prepare_lift_forward_model(tel, src, adaptive_basis;
        diversity_opd=adaptive_diversity, focal_resolution=8)
    lift_adaptive_truth = LiFT(adaptive_forward; iterations=3,
        numerical=false, solve_mode=LiFTSolveNormalEquations(),
        damping=LiFTAdaptiveLevenbergMarquardt())
    adaptive_target = copy(intensity_values(evaluate_lift_forward!(
        adaptive_forward, adaptive_truth_opd)))
    adaptive_observation = LiFTObservation(adaptive_forward,
        adaptive_target)
    adaptive_estimate = WavefrontSensors.reconstruct(lift_adaptive_truth,
        adaptive_observation; check_convergence=false)
    @test adaptive_estimate ≈ adaptive_truth rtol=2e-3 atol=eps(Float64)

    det_binned = Detector(noise=NoiseNone(), psf_sampling=1, binning=2)
    frame_binned = capture!(det_binned, psf; rng=MersenneTwister(3))
    binned_mapping = LiFTFrameMapping(binning=2)
    binned_forward = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8,
        mapping=binned_mapping)
    lift_binned = LiFT(binned_forward; iterations=2,
        mode_ids=(1, 2), numerical=true)
    binned_observation = LiFTObservation(binned_forward,
        copy(frame_binned); domain=LiFTExpectedCounts(
            det_binned.params.integration_time;
            quantum_efficiency=det_binned.params.qe))
    coeffs_binned = WavefrontSensors.reconstruct(
        lift_binned, binned_observation)
    @test length(coeffs_binned) == 2
    @test all(isfinite, coeffs_binned)
    @test lift_observation_contract(binned_forward).rate_metadata.sampling ==
        2 .* lift_observation_contract(forward).rate_metadata.sampling
    binned_model_rate = copy(intensity_values(evaluate_lift_forward!(
        binned_forward, diversity)))
    unmapped_model_rate = copy(intensity_values(evaluate_lift_forward!(
        forward, diversity)))
    expected_binned_rate = similar(binned_model_rate)
    bin2d!(expected_binned_rate, unmapped_model_rate, 2)
    @test binned_model_rate ≈ expected_binned_rate rtol=1e-12 atol=1e-8
    @test sum(binned_model_rate) ≈ sum(unmapped_model_rate) rtol=1e-12

    lift_response = GaussianPixelResponse(response_width_px=0.6)
    response_mapping = LiFTFrameMapping(response=lift_response)
    response_forward = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8,
        mapping=response_mapping)
    response_model_rate = copy(intensity_values(evaluate_lift_forward!(
        response_forward, diversity)))
    expected_response_rate = copy(unmapped_model_rate)
    response_scratch = similar(expected_response_rate)
    AdaptiveOpticsSim.Detectors.apply_response!(AdaptiveOpticsSim.Backends.ScalarCPUStyle(),
        lift_response, expected_response_rate, response_scratch)
    @test response_model_rate ≈ expected_response_rate rtol=1e-12 atol=1e-8
    @test supports_detector_mtf(lift_response)
    @test 0 <= detector_mtf(lift_response, 0.25, 0.25) <= 1

    det_readout = Detector(noise=NoiseReadout(1e-3), psf_sampling=1)
    lift_readout = LiFT(forward; iterations=2,
        mode_ids=(1, 2), numerical=true)
    readout_observation = LiFTObservation(forward, copy(psf);
        readout_noise_std=det_readout.noise.sigma)
    coeffs_readout = WavefrontSensors.reconstruct(
        lift_readout, readout_observation)
    @test length(coeffs_readout) == 2
    @test all(isfinite, coeffs_readout)

    lift_weighted = LiFT(forward; iterations=1,
        mode_ids=(1, 2), numerical=false)
    initial_model = copy(intensity_values(evaluate_lift_forward!(
        forward, diversity)))
    initial_model .*= sum(psf) / sum(initial_model)
    expected_weights = sqrt.(inv.(max.(initial_model .+ 1e-6, eps(Float64))))
    WavefrontSensors.reconstruct(lift_weighted, readout_observation;
        R_n=:model,
        check_convergence=false)
    @test lift_weighted.state.weight_buffer ≈ vec(expected_weights)

    exposure_time = 0.25
    quantum_efficiency = 0.8
    rate_domain = LiFTPhotonRate(
        noise_equivalent_exposure_s=exposure_time,
        quantum_efficiency=quantum_efficiency)
    count_domain = LiFTExpectedCounts(exposure_time;
        quantum_efficiency=quantum_efficiency)
    rate_observation = LiFTObservation(adaptive_forward,
        copy(adaptive_target); domain=rate_domain)
    count_values = adaptive_target .* (exposure_time * quantum_efficiency)
    count_observation = LiFTObservation(adaptive_forward,
        count_values; domain=count_domain)
    rate_estimator = LiFT(adaptive_forward; iterations=3,
        solve_mode=LiFTSolveNormalEquations())
    count_estimator = LiFT(adaptive_forward; iterations=3,
        solve_mode=LiFTSolveNormalEquations())
    rate_estimate = WavefrontSensors.reconstruct(
        rate_estimator, rate_observation;
        optimize_norm=:none, check_convergence=false)
    count_estimate = WavefrontSensors.reconstruct(
        count_estimator, count_observation;
        optimize_norm=:none, check_convergence=false)
    @test count_estimate ≈ rate_estimate rtol=1e-11 atol=1e-18

    predicted_counts = similar(adaptive_target)
    predict_lift_observation!(predicted_counts, adaptive_forward,
        adaptive_truth_opd, count_domain)
    @test predicted_counts ≈ count_values rtol=1e-12
    normalized_domain = LiFTNormalizedIntensity(sum(adaptive_target))
    predicted_normalized = similar(adaptive_target)
    predict_lift_observation!(predicted_normalized, adaptive_forward,
        adaptive_truth_opd, normalized_domain)
    @test sum(predicted_normalized) ≈ 1.0 rtol=1e-12

    mismatched_source = Source(band=:H, magnitude=0.0)
    wavelength_mismatch = prepare_lift_forward_model(tel,
        mismatched_source, basis; diversity_opd=diversity,
        focal_resolution=8)
    wavelength_observation = LiFTObservation(wavelength_mismatch,
        copy(psf))
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct(lift,
        wavelength_observation)
    unchanged_coefficients = fill(1.0, 3)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        unchanged_coefficients,
        lift, wavelength_observation)
    @test unchanged_coefficients == fill(1.0, 3)
    preprocessing_mismatch = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8,
        mapping=LiFTFrameMapping())
    preprocessing_observation = LiFTObservation(preprocessing_mismatch,
        copy(psf))
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct(lift,
        preprocessing_observation)
    @test_throws DimensionMismatchError WavefrontSensors.reconstruct!(
        unchanged_coefficients,
        lift, observation; R_n=zeros(7, 8))
    @test unchanged_coefficients == fill(1.0, 3)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        unchanged_coefficients,
        lift, observation; R_n=zeros(Float32, 8, 8))
    @test unchanged_coefficients == fill(1.0, 3)
    mismatched_coefficients = fill(1.0f0, 3)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        mismatched_coefficients,
        lift, observation)
    @test mismatched_coefficients == fill(1.0f0, 3)
    @test_throws DimensionMismatchError WavefrontSensors.reconstruct!(
        unchanged_coefficients,
        lift, observation; coeffs0=zeros(2))
    @test unchanged_coefficients == fill(1.0, 3)
    @test_throws InvalidConfiguration WavefrontSensors.reconstruct!(
        unchanged_coefficients,
        lift, observation; coeffs0=zeros(Float32, 3))
    @test unchanged_coefficients == fill(1.0, 3)
    @test_throws DimensionMismatchError LiFTObservation(
        lift_observation_contract(forward), zeros(7, 8))
    @test_throws InvalidConfiguration LiFTExpectedCounts(0.0)

    frozen_basis = rand(MersenneTwister(0x1f7), 8, 8, 3)
    frozen_forward = prepare_lift_forward_model(tel, src, frozen_basis;
        diversity_opd=diversity, focal_resolution=8)
    frozen_numerical = LiFT(frozen_forward; iterations=2,
        mode_ids=(1, 2), numerical=true)
    frozen_analytic = LiFT(frozen_forward; iterations=2,
        mode_ids=(1, 2), numerical=false)
    frozen_H_numerical = lift_interaction_matrix(frozen_numerical,
        zeros(3))
    frozen_H_analytic = lift_interaction_matrix(frozen_analytic,
        zeros(3))
    @test norm(frozen_H_numerical) ≈ 8.049469666421768e16 rtol=1e-12
    @test norm(frozen_H_analytic) ≈ 8.049503480686035e16 rtol=1e-12

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    sep_forward = prepare_lift_forward_model(tel, src, basis;
        diversity_opd=diversity, focal_resolution=8,
        object_kernel=sep_kernel)
    lift_sep = LiFT(sep_forward; iterations=2, numerical=false)
    @test sep_forward.model.object_kernel isa
        WavefrontSensors.LiFTSeparableObjectKernel
    dense_conv = similar(psf)
    sep_conv = similar(psf)
    tmp_conv = similar(psf)
    WavefrontSensors.conv2d_same!(dense_conv, psf, sep_kernel)
    WavefrontSensors.conv2d_same_separable!(
        sep_conv,
        tmp_conv,
        psf,
        sep_forward.model.object_kernel.row,
        sep_forward.model.object_kernel.col,
    )
    @test isapprox(sep_conv, dense_conv; rtol=1e-6, atol=1e-6)
    single_pixel_convolution = zeros(1, 1)
    WavefrontSensors.conv2d_same!(single_pixel_convolution,
        fill(2.0, 1, 1), sep_kernel)
    @test single_pixel_convolution == fill(2.0, 1, 1)
end
