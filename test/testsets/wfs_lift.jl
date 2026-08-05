function lift_hot_path_allocation_bytes!(H, full_coefficients, lift)
    WavefrontSensors.reconstruct!(lift)
    WavefrontSensors.lift_interaction_matrix!(H, lift, full_coefficients)
    interaction_bytes = @allocated WavefrontSensors.lift_interaction_matrix!(
        H, lift, full_coefficients)
    reconstruction_bytes = @allocated WavefrontSensors.reconstruct!(lift)
    return interaction_bytes, reconstruction_bytes
end

function prepare_test_lift(forward, observation; mode_ids=nothing, kwargs...)
    definition = LiFT(; mode_ids, kwargs...)
    count = mode_ids === nothing ?
        size(WavefrontSensors.lift_forward_plan(forward).basis, 3) :
        length(mode_ids)
    coefficients = similar(observation.values, eltype(observation.values), count)
    fill!(coefficients, zero(eltype(coefficients)))
    return prepare_lift_estimator(definition, forward, observation,
        coefficients)
end

function replace_lift_optical_rate_workspace(workspace, optical_rate)
    return WavefrontSensors.LiFTForwardWorkspace(workspace.propagation,
        optical_rate, workspace.amplitude_buffer, workspace.field_scratch,
        workspace.focal_buffer, workspace.mode_buffer,
        workspace.conjugate_field_buffer, workspace.response_buffer,
        workspace.response_scratch, workspace.sampling_buffer,
        workspace.mapped_rate_buffer, workspace.output_work_buffer,
        workspace.convolution_buffer, workspace.convolution_scratch,
        workspace.opd_work_buffer)
end

@testset "LiFT" begin
    @test Docs.hasdoc(WavefrontSensors, :evaluate_lift_forward!)
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
    fallback_diagnostics = WavefrontSensors.LiFTDiagnosticsWorkspace(
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
    forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8)
    @test forward.output.values !== forward.workspace.optical_rate_buffer
    @test !AdaptiveOpticsSim.WavefrontSensors._wfs_storage_mightalias(
        forward.output.values, forward.workspace.optical_rate_buffer)
    psf = reference_direct_image(tel, src; zero_padding=1)
    observation = LiFTObservation(forward, copy(psf))
    lift = prepare_test_lift(forward, observation; iterations=2,
        mode_ids=(1, 2), jacobian_method=LiFTNumericalJacobian())
    H = lift_interaction_matrix(lift, zeros(3))
    @test size(H) == (64, 2)
    @test maximum(abs, H) > 0
    @test_throws DimensionMismatchError lift_interaction_matrix!(
        zeros(63, 2), lift, zeros(3))
    @test_throws DimensionMismatchError lift_interaction_matrix!(
        zeros(64, 2), lift, zeros(2))
    untouched_H = zeros(64, 2)
    @test_throws InvalidConfiguration lift_interaction_matrix!(
        untouched_H, lift, zeros(3); rate_scale=-1)
    @test all(iszero, untouched_H)
    aliased_H = zeros(64, 2)
    @test_throws InvalidConfiguration lift_interaction_matrix!(
        aliased_H, lift, @view(vec(aliased_H)[1:3]))
    kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    forward_analytic = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8, object_kernel=kernel)
    analytic_observation = LiFTObservation(forward_analytic, copy(psf))
    lift_analytic = prepare_test_lift(forward_analytic,
        analytic_observation; iterations=2, mode_ids=(1, 2),
        jacobian_method=LiFTAnalyticJacobian())
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3))
    @test size(H_analytic) == (64, 2)
    coeffs = WavefrontSensors.reconstruct(lift)
    @test length(coeffs) == 2
    @test fieldnames(typeof(lift)) == (:forward, :plan, :workspace,
        :observation, :coefficients, :initial_coefficients, :backend, :device)
    @test !hasfield(typeof(lift), :detector)
    @test !hasfield(typeof(lift), :telescope)
    @test observation.values !== intensity_values(lift_forward_output(forward))
    @test collect(lift.plan.mode_ids) == [1, 2]
    @test typeof(LiFT(mode_ids=1:2).mode_ids) ===
        typeof(LiFT(mode_ids=1:3).mode_ids)
    @test_throws InvalidConfiguration LiFT(mode_ids=())
    @test_throws InvalidConfiguration LiFT(mode_ids=(1, 1))
    @test_throws DimensionMismatchError prepare_lift_estimator(
        LiFT(mode_ids=(1, 4)), forward, observation, zeros(2))
    @test_throws MethodError LiFT(forward; iterations=2)
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
    lift_normal = prepare_test_lift(forward, observation; iterations=2,
        jacobian_method=LiFTNumericalJacobian(), mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    coeffs_normal = WavefrontSensors.reconstruct(lift_normal)
    @test length(coeffs_normal) == 2
    @test all(isfinite, coeffs_normal)
    @test !diagnostics(lift_normal).used_qr

    allocation_values = copy(intensity_values(evaluate_lift_forward!(
        forward)))
    allocation_observation = LiFTObservation(forward, allocation_values)
    allocation_lift = prepare_test_lift(forward,
        allocation_observation; iterations=1, mode_ids=(1, 2),
        jacobian_method=LiFTAnalyticJacobian(),
        solve_mode=LiFTSolveNormalEquations(),
        model_scaling=LiFTPhysicalRatePreservation(),
        check_convergence=false)
    allocation_coefficients = allocation_lift.workspace.full_coefficients_buffer
    allocation_H = allocation_lift.workspace.H_buffer
    @test (@inferred WavefrontSensors.reconstruct!(allocation_lift)) ===
        allocation_lift.coefficients
    if coverage_instrumented()
        @test_skip "LiFT allocation assertions are disabled under coverage instrumentation"
    else
        interaction_bytes, reconstruction_bytes =
            lift_hot_path_allocation_bytes!(allocation_H,
                allocation_coefficients, allocation_lift)
        @test interaction_bytes == 0
        @test reconstruction_bytes == 0
    end

    tel32 = Telescope(resolution=8, diameter=8f0,
        central_obstruction=0f0, T=Float32)
    src32 = Source(band=:I, magnitude=0f0, T=Float32)
    det32 = Detector(noise=NoiseNone(), psf_sampling=1, T=Float32)
    basis32 = rand(MersenneTwister(29), Float32, 8, 8, 3)
    diversity32 = zeros(Float32, 8, 8)
    forward32 = prepare_lift_forward_model(tel32, src32, basis32,
        diversity32; diversity_opd=diversity32, focal_resolution=8)
    psf32 = reference_direct_image(tel32, src32; zero_padding=1)
    observation32 = LiFTObservation(forward32, copy(psf32))
    lift_normal32 = prepare_test_lift(forward32, observation32;
        iterations=2, mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations())
    coeffs_normal32 = WavefrontSensors.reconstruct(lift_normal32)
    @test all(isfinite, coeffs_normal32)
    @test all(isfinite,
        @view(lift_normal32.workspace.normal_buffer[1:2, 1:2]))

    lift_damped = prepare_test_lift(forward, observation; iterations=2,
        jacobian_method=LiFTNumericalJacobian(), mode_ids=(1, 2),
        damping=LiFTLevenbergMarquardt())
    coeffs_damped = WavefrontSensors.reconstruct(lift_damped)
    @test length(coeffs_damped) == 2
    @test all(isfinite, coeffs_damped)
    @test diagnostics(lift_damped).regularization >= 0
    lift_adaptive = prepare_test_lift(forward, observation; iterations=2,
        jacobian_method=LiFTNumericalJacobian(), mode_ids=(1, 2),
        damping=LiFTAdaptiveLevenbergMarquardt())
    coeffs_adaptive = WavefrontSensors.reconstruct(lift_adaptive)
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
    adaptive_forward = prepare_lift_forward_model(tel, src, adaptive_basis,
        adaptive_truth_opd;
        diversity_opd=adaptive_diversity, focal_resolution=8)
    adaptive_target = copy(intensity_values(evaluate_lift_forward!(
        adaptive_forward)))
    adaptive_observation = LiFTObservation(adaptive_forward,
        adaptive_target)
    lift_adaptive_truth = prepare_test_lift(adaptive_forward,
        adaptive_observation; iterations=3,
        jacobian_method=LiFTAnalyticJacobian(),
        solve_mode=LiFTSolveNormalEquations(),
        damping=LiFTAdaptiveLevenbergMarquardt(),
        check_convergence=false)
    adaptive_estimate = WavefrontSensors.reconstruct(lift_adaptive_truth)
    @test adaptive_estimate ≈ adaptive_truth rtol=2e-3 atol=eps(Float64)
    initialized_product = zeros(2)
    initialized_lift = prepare_lift_estimator(
        LiFT(iterations=1, solve_mode=LiFTSolveNormalEquations(),
            model_scaling=LiFTPhysicalRatePreservation(),
            check_convergence=false),
        adaptive_forward, adaptive_observation, initialized_product;
        initial_coefficients=copy(adaptive_truth))
    WavefrontSensors.reconstruct!(initialized_lift)
    @test initialized_product ≈ adaptive_truth rtol=2e-3 atol=eps(Float64)

    det_binned = Detector(noise=NoiseNone(), psf_sampling=1, binning=2)
    frame_binned = capture!(det_binned, psf; rng=MersenneTwister(3))
    binned_mapping = LiFTFrameMapping(binning=2)
    binned_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8,
        mapping=binned_mapping)
    binned_observation = LiFTObservation(binned_forward,
        copy(frame_binned); domain=LiFTExpectedCounts(
            det_binned.params.exposure_duration;
            quantum_efficiency=det_binned.params.qe))
    lift_binned = prepare_test_lift(binned_forward, binned_observation;
        iterations=2, mode_ids=(1, 2),
        jacobian_method=LiFTNumericalJacobian())
    coeffs_binned = WavefrontSensors.reconstruct(lift_binned)
    @test length(coeffs_binned) == 2
    @test all(isfinite, coeffs_binned)
    @test lift_observation_contract(binned_forward).rate_metadata.sampling ==
        2 .* lift_observation_contract(forward).rate_metadata.sampling
    binned_model_rate = copy(intensity_values(evaluate_lift_forward!(
        binned_forward)))
    unmapped_model_rate = copy(intensity_values(evaluate_lift_forward!(
        forward)))
    expected_binned_rate = similar(binned_model_rate)
    bin2d!(expected_binned_rate, unmapped_model_rate, 2)
    @test binned_model_rate ≈ expected_binned_rate rtol=1e-12 atol=1e-8
    @test sum(binned_model_rate) ≈ sum(unmapped_model_rate) rtol=1e-12

    lift_response = GaussianPixelResponse(response_width_px=0.6)
    response_mapping = LiFTFrameMapping(response=lift_response)
    response_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8,
        mapping=response_mapping)
    @test response_forward.output.values !==
        response_forward.workspace.mapped_rate_buffer
    response_model_rate = copy(intensity_values(evaluate_lift_forward!(
        response_forward)))
    prepared_response = WavefrontSensors.lift_forward_plan(
        response_forward).mapping.response
    @test prepared_response.kernel !== lift_response.kernel
    prepared_response_snapshot = copy(prepared_response.kernel)
    fill!(lift_response.kernel, 0)
    @test prepared_response.kernel == prepared_response_snapshot
    @test intensity_values(evaluate_lift_forward!(response_forward)) ==
        response_model_rate
    expected_response_rate = copy(unmapped_model_rate)
    response_scratch = similar(expected_response_rate)
    AdaptiveOpticsSim.Detectors.apply_response!(AdaptiveOpticsSim.Backends.ScalarCPUStyle(),
        prepared_response, expected_response_rate, response_scratch)
    @test response_model_rate ≈ expected_response_rate rtol=1e-12 atol=1e-8
    @test supports_detector_mtf(prepared_response)
    @test 0 <= detector_mtf(prepared_response, 0.25, 0.25) <= 1

    det_readout = Detector(noise=NoiseReadout(1e-3), psf_sampling=1)
    readout_observation = LiFTObservation(forward, copy(psf);
        readout_noise_std=det_readout.noise.sigma)
    lift_readout = prepare_test_lift(forward, readout_observation;
        iterations=2, mode_ids=(1, 2),
        jacobian_method=LiFTNumericalJacobian())
    coeffs_readout = WavefrontSensors.reconstruct(lift_readout)
    @test length(coeffs_readout) == 2
    @test all(isfinite, coeffs_readout)

    initial_model = copy(intensity_values(evaluate_lift_forward!(
        forward)))
    initial_model .*= sum(psf) / sum(initial_model)
    expected_weights = sqrt.(inv.(max.(initial_model .+ 1e-6, eps(Float64))))
    lift_weighted = prepare_test_lift(forward, readout_observation;
        iterations=1, mode_ids=(1, 2),
        jacobian_method=LiFTAnalyticJacobian(),
        weighting=LiFTInitialModelWeighting(), check_convergence=false)
    WavefrontSensors.reconstruct!(lift_weighted)
    @test lift_weighted.workspace.weight_buffer ≈ vec(expected_weights)

    initial_model_weights = zeros(length(initial_model))
    iterative_model_weights = similar(initial_model_weights)
    WavefrontSensors.update_weights!(initial_model_weights,
        LiFTInitialModelWeighting(), initial_model,
        readout_observation.metadata, 1)
    fixed_initial_weights = copy(initial_model_weights)
    WavefrontSensors.update_weights!(initial_model_weights,
        LiFTInitialModelWeighting(), 2 .* initial_model,
        readout_observation.metadata, 2)
    @test initial_model_weights == fixed_initial_weights
    WavefrontSensors.update_weights!(iterative_model_weights,
        LiFTIterativeModelWeighting(), initial_model,
        readout_observation.metadata, 1)
    first_iterative_weights = copy(iterative_model_weights)
    WavefrontSensors.update_weights!(iterative_model_weights,
        LiFTIterativeModelWeighting(), 2 .* initial_model,
        readout_observation.metadata, 2)
    @test iterative_model_weights != first_iterative_weights

    external_variance = ones(8, 8)
    variance_lift = prepare_test_lift(forward, readout_observation;
        iterations=1, mode_ids=(1, 2),
        weighting=LiFTVarianceMapWeighting(external_variance),
        check_convergence=false)
    prepared_variance = WavefrontSensors.lift_estimation_plan(
        variance_lift).weighting.variance
    @test prepared_variance !== external_variance
    fill!(external_variance, 2)
    @test all(isone, prepared_variance)

    exposure_duration = 0.25
    quantum_efficiency = 0.8
    rate_domain = LiFTPhotonRate(
        noise_equivalent_exposure_s=exposure_duration,
        quantum_efficiency=quantum_efficiency)
    count_domain = LiFTExpectedCounts(exposure_duration;
        quantum_efficiency=quantum_efficiency)
    rate_observation = LiFTObservation(adaptive_forward,
        copy(adaptive_target); domain=rate_domain)
    count_values = adaptive_target .* (exposure_duration * quantum_efficiency)
    count_observation = LiFTObservation(adaptive_forward,
        count_values; domain=count_domain)
    rate_estimator = prepare_test_lift(adaptive_forward, rate_observation;
        iterations=3, solve_mode=LiFTSolveNormalEquations(),
        model_scaling=LiFTPhysicalRatePreservation(),
        check_convergence=false)
    count_estimator = prepare_test_lift(adaptive_forward,
        count_observation; iterations=3,
        solve_mode=LiFTSolveNormalEquations(),
        model_scaling=LiFTPhysicalRatePreservation(),
        check_convergence=false)
    rate_estimate = WavefrontSensors.reconstruct(rate_estimator)
    count_estimate = WavefrontSensors.reconstruct(count_estimator)
    @test count_estimate ≈ rate_estimate rtol=1e-11 atol=1e-18

    predicted_counts = similar(adaptive_target)
    predict_lift_observation!(predicted_counts, adaptive_forward,
        count_domain)
    @test predicted_counts ≈ count_values rtol=1e-12
    prediction_output_snapshot = copy(adaptive_forward.output.values)
    @test_throws InvalidConfiguration predict_lift_observation!(
        adaptive_forward.output.values, adaptive_forward, count_domain)
    @test adaptive_forward.output.values == prediction_output_snapshot
    prediction_input_snapshot = copy(adaptive_forward.input)
    @test_throws InvalidConfiguration predict_lift_observation!(
        adaptive_forward.input, adaptive_forward, count_domain)
    @test adaptive_forward.input == prediction_input_snapshot
    normalized_domain = LiFTNormalizedIntensity(sum(adaptive_target))
    predicted_normalized = similar(adaptive_target)
    predict_lift_observation!(predicted_normalized, adaptive_forward,
        normalized_domain)
    @test sum(predicted_normalized) ≈ 1.0 rtol=1e-12

    mismatched_source = Source(band=:H, magnitude=0.0)
    wavelength_mismatch = prepare_lift_forward_model(tel,
        mismatched_source, basis, diversity; diversity_opd=diversity,
        focal_resolution=8)
    wavelength_observation = LiFTObservation(wavelength_mismatch,
        copy(psf))
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, wavelength_observation, zeros(2))
    preprocessing_mismatch = prepare_lift_forward_model(tel, src, basis,
        diversity;
        diversity_opd=diversity, focal_resolution=8,
        mapping=LiFTFrameMapping())
    preprocessing_observation = LiFTObservation(preprocessing_mismatch,
        copy(psf))
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, preprocessing_observation, zeros(2))
    @test_throws DimensionMismatchError prepare_lift_estimator(
        LiFT(mode_ids=(1, 2),
            weighting=LiFTVarianceMapWeighting(zeros(7, 8))),
        forward, observation, zeros(2))
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2),
            weighting=LiFTVarianceMapWeighting(zeros(Float32, 8, 8))),
        forward, observation, zeros(2))
    @test_throws DimensionMismatchError prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, zeros(3))
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, zeros(Float32, 2))
    @test_throws DimensionMismatchError prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, zeros(2);
        initial_coefficients=zeros(2))
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, zeros(2);
        initial_coefficients=zeros(Float32, 3))
    aliased_coefficients = @view observation.values[1:2]
    observation_snapshot = copy(observation.values)
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, aliased_coefficients)
    @test observation.values == observation_snapshot
    aliased_initial = @view observation.values[1:3]
    @test_throws InvalidConfiguration prepare_lift_estimator(
        LiFT(mode_ids=(1, 2)), forward, observation, zeros(2);
        initial_coefficients=aliased_initial)
    @test_throws DimensionMismatchError LiFTObservation(
        lift_observation_contract(forward), zeros(7, 8))
    @test_throws InvalidConfiguration LiFTExpectedCounts(0.0)
    @test_throws MethodError evaluate_lift_forward!(forward, diversity)

    independent_lift = prepare_test_lift(forward, observation;
        mode_ids=(1, 2), check_convergence=false)
    @test independent_lift.forward.plan === lift.forward.plan
    @test independent_lift.forward.input == lift.forward.plan.diversity_opd
    @test all(isfinite, intensity_values(evaluate_lift_forward!(
        independent_lift.forward)))
    @test independent_lift.forward.workspace !== lift.forward.workspace
    @test independent_lift.workspace !== lift.workspace
    @test independent_lift.coefficients !== lift.coefficients

    malformed_forward = PreparedLiFTForward(forward.plan, forward.workspace,
        forward.output.values, forward.output, forward.backend, forward.device)
    malformed_output_snapshot = copy(forward.output.values)
    @test_throws InvalidConfiguration evaluate_lift_forward!(malformed_forward)
    @test forward.output.values == malformed_output_snapshot
    workspace_output = IntensityMap(forward.output.metadata,
        forward.workspace.optical_rate_buffer)
    workspace_output_snapshot = copy(workspace_output.values)
    malformed_workspace_output = PreparedLiFTForward(forward.plan,
        forward.workspace, forward.input, workspace_output, forward.backend,
        forward.device)
    @test_throws InvalidConfiguration evaluate_lift_forward!(
        malformed_workspace_output)
    @test workspace_output.values == workspace_output_snapshot
    malformed_workspace = replace_lift_optical_rate_workspace(
        forward.workspace, zeros(7, 8))
    malformed_workspace_forward = PreparedLiFTForward(forward.plan,
        malformed_workspace, forward.input, forward.output, forward.backend,
        forward.device)
    malformed_workspace_snapshot = copy(forward.output.values)
    @test_throws DimensionMismatchError evaluate_lift_forward!(
        malformed_workspace_forward)
    @test forward.output.values == malformed_workspace_snapshot

    concurrent_first = prepare_test_lift(forward, allocation_observation;
        iterations=1, mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations(),
        model_scaling=LiFTPhysicalRatePreservation(),
        check_convergence=false)
    concurrent_second = prepare_test_lift(forward, allocation_observation;
        iterations=1, mode_ids=(1, 2),
        solve_mode=LiFTSolveNormalEquations(),
        model_scaling=LiFTPhysicalRatePreservation(),
        check_convergence=false)
    first_task = Threads.@spawn WavefrontSensors.reconstruct!(concurrent_first)
    second_task = Threads.@spawn WavefrontSensors.reconstruct!(concurrent_second)
    fetch(first_task)
    fetch(second_task)
    @test concurrent_first.coefficients ≈ concurrent_second.coefficients
    @test concurrent_first.forward.workspace !==
        concurrent_second.forward.workspace

    @test_throws InvalidConfiguration prepare_lift_forward_model(
        tel, src, basis, diversity; diversity_opd=diversity,
        focal_resolution=8, object_kernel=zeros(3, 3))
    @test_throws InvalidConfiguration prepare_lift_forward_model(
        tel, src, basis, diversity; diversity_opd=diversity,
        focal_resolution=8, object_kernel=fill(-1.0, 3, 3))

    frozen_basis = rand(MersenneTwister(0x1f7), 8, 8, 3)
    frozen_forward = prepare_lift_forward_model(tel, src, frozen_basis,
        diversity;
        diversity_opd=diversity, focal_resolution=8)
    frozen_observation = LiFTObservation(frozen_forward, copy(psf))
    frozen_numerical = prepare_test_lift(frozen_forward,
        frozen_observation; iterations=2, mode_ids=(1, 2),
        jacobian_method=LiFTNumericalJacobian())
    frozen_analytic = prepare_test_lift(frozen_forward,
        frozen_observation; iterations=2, mode_ids=(1, 2),
        jacobian_method=LiFTAnalyticJacobian())
    frozen_H_numerical = lift_interaction_matrix(frozen_numerical,
        zeros(3))
    frozen_H_analytic = lift_interaction_matrix(frozen_analytic,
        zeros(3))
    @test norm(frozen_H_numerical) ≈ 8.049469666421768e16 rtol=1e-12
    @test norm(frozen_H_analytic) ≈ 8.049503480686035e16 rtol=1e-12

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    sep_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8,
        object_kernel=sep_kernel)
    @test sep_forward.plan.object_kernel isa
        WavefrontSensors.LiFTSeparableObjectKernel
    dense_conv = similar(psf)
    sep_conv = similar(psf)
    tmp_conv = similar(psf)
    WavefrontSensors.conv2d_same!(dense_conv, psf, sep_kernel)
    WavefrontSensors.conv2d_same_separable!(
        sep_conv,
        tmp_conv,
        psf,
        sep_forward.plan.object_kernel.row,
        sep_forward.plan.object_kernel.col,
    )
    @test isapprox(sep_conv, dense_conv; rtol=1e-6, atol=1e-6)
    single_pixel_convolution = zeros(1, 1)
    WavefrontSensors.conv2d_same!(single_pixel_convolution,
        fill(2.0, 1, 1), sep_kernel)
    @test single_pixel_convolution == fill(2.0, 1, 1)
end
