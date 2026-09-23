function replace_lift_optical_rate_workspace(workspace, optical_rate)
    return WavefrontSensors.LiFTForwardWorkspace(workspace.propagation,
        optical_rate, workspace.amplitude_buffer, workspace.field_scratch,
        workspace.focal_buffer, workspace.mode_buffer,
        workspace.conjugate_field_buffer, workspace.response_buffer,
        workspace.response_scratch, workspace.sampling_buffer,
        workspace.mapped_rate_buffer, workspace.convolution_buffer,
        workspace.convolution_scratch)
end

function prepare_aoc_lift(forward, observation, domain;
    initial_coefficients=nothing, read_noise_std=zero(eltype(observation)),
    kwargs...)
    model = LiFTForwardModel(forward)
    specification = AOCPhaseRetrieval.LiFTSpecification(model, domain;
        read_noise_std)
    plan = AdaptiveOpticsCalibration.prepare(
        AOCPhaseRetrieval.LiFT(; kwargs...), specification)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = AOCPhaseRetrieval.LiFTInputs(observation; initial_coefficients)
    return (; model, specification, plan, result, workspace, inputs)
end

function process_aoc_lift!(prepared)
    return AdaptiveOpticsCalibration.process!(prepared.result,
        prepared.workspace, prepared.plan, prepared.inputs)
end

function process_aoc_observation(forward, observation; kwargs...)
    specification = AOCPhaseRetrieval.LiFTSpecification(forward, observation)
    plan = AdaptiveOpticsCalibration.prepare(
        AOCPhaseRetrieval.LiFT(; kwargs...), specification)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = AOCPhaseRetrieval.LiFTInputs(observation.values)
    AdaptiveOpticsCalibration.process!(result, workspace, plan, inputs)
    return (; specification, plan, result, workspace, inputs)
end

@testset "LiFT physical forward and AdaptiveOpticsCalibration inverse" begin
    @test Docs.hasdoc(WavefrontSensors, :evaluate_lift_forward!)

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    basis = rand(MersenneTwister(41), 8, 8, 3)
    diversity = zeros(8, 8)
    forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8)

    @test forward.output.values !== forward.workspace.optical_rate_buffer
    @test !AdaptiveOpticsSim.WavefrontSensors._wfs_storage_mightalias(
        forward.output.values, forward.workspace.optical_rate_buffer)
    rate_output = evaluate_lift_forward!(forward)
    @test rate_output === lift_forward_output(forward)
    @test all(isfinite, intensity_values(rate_output))
    @test all(>=(0), intensity_values(rate_output))
    rate_snapshot = copy(intensity_values(rate_output))
    @test_throws InvalidConfiguration evaluate_lift_forward!(forward;
        rate_scale=-1)
    @test intensity_values(rate_output) == rate_snapshot
    @test lift_observation_contract(forward).rate_metadata.dimensions == (8, 8)
    @test lift_observation_contract(forward).rate_metadata.normalization isa
        PhotonRateNormalization
    rate_observation = LiFTObservation(forward,
        copy(intensity_values(rate_output)))
    @test rate_observation.metadata.contract === lift_observation_contract(forward)
    @test rate_observation.metadata.domain isa LiFTPhotonRate
    @test_throws DimensionMismatchError LiFTObservation(
        lift_observation_contract(forward), zeros(7, 8))
    @test_throws InvalidConfiguration LiFTExpectedCounts(0.0)

    count_domain = LiFTExpectedCounts(0.25; quantum_efficiency=0.8)
    count_prediction = similar(intensity_values(rate_output))
    predict_lift_observation!(count_prediction, forward, count_domain)
    @test count_prediction ≈ 0.2 .* intensity_values(rate_output) rtol=1e-12
    normalized_prediction = similar(count_prediction)
    predict_lift_observation!(normalized_prediction, forward,
        LiFTNormalizedIntensity(sum(intensity_values(rate_output))))
    @test sum(normalized_prediction) ≈ 1.0 rtol=1e-12
    output_snapshot = copy(forward.output.values)
    @test_throws InvalidConfiguration predict_lift_observation!(
        forward.output.values, forward, count_domain)
    @test forward.output.values == output_snapshot
    input_snapshot = copy(forward.input)
    @test_throws InvalidConfiguration predict_lift_observation!(
        forward.input, forward, count_domain)
    @test forward.input == input_snapshot

    binned_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8,
        mapping=LiFTFrameMapping(binning=2))
    binned_rate = copy(intensity_values(evaluate_lift_forward!(binned_forward)))
    unbinned_rate = copy(intensity_values(evaluate_lift_forward!(forward)))
    expected_binned_rate = similar(binned_rate)
    bin2d!(expected_binned_rate, unbinned_rate, 2)
    @test binned_rate ≈ expected_binned_rate rtol=1e-12 atol=1e-8
    @test sum(binned_rate) ≈ sum(unbinned_rate) rtol=1e-12
    @test lift_observation_contract(binned_forward).rate_metadata.sampling ==
        2 .* lift_observation_contract(forward).rate_metadata.sampling
    @test lift_observation_contract(binned_forward).preprocessing_signature !=
        lift_observation_contract(forward).preprocessing_signature

    lift_response = GaussianPixelResponse(response_width_px=0.6)
    response_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8,
        mapping=LiFTFrameMapping(response=lift_response))
    response_rate = copy(intensity_values(evaluate_lift_forward!(response_forward)))
    prepared_response = lift_forward_plan(response_forward).mapping.response
    @test prepared_response.kernel !== lift_response.kernel
    prepared_response_snapshot = copy(prepared_response.kernel)
    fill!(lift_response.kernel, 0)
    @test prepared_response.kernel == prepared_response_snapshot
    expected_response_rate = copy(unbinned_rate)
    response_scratch = similar(expected_response_rate)
    AdaptiveOpticsSim.Detectors.apply_response!(
        AdaptiveOpticsSim.Backends.ScalarCPUStyle(), prepared_response,
        expected_response_rate, response_scratch)
    @test response_rate ≈ expected_response_rate rtol=1e-12 atol=1e-8
    @test intensity_values(evaluate_lift_forward!(response_forward)) == response_rate
    @test supports_detector_mtf(prepared_response)
    @test 0 <= detector_mtf(prepared_response, 0.25, 0.25) <= 1
    @test lift_observation_contract(response_forward).preprocessing_signature !=
        lift_observation_contract(forward).preprocessing_signature

    dense_kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    dense_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8, object_kernel=dense_kernel)
    dense_rate = copy(intensity_values(evaluate_lift_forward!(dense_forward)))
    expected_dense_rate = similar(unbinned_rate)
    WavefrontSensors.conv2d_same!(expected_dense_rate, unbinned_rate,
        dense_kernel)
    @test dense_rate ≈ expected_dense_rate rtol=1e-12 atol=1e-8

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    sep_forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=8, object_kernel=sep_kernel)
    @test sep_forward.plan.object_kernel isa
        WavefrontSensors.LiFTSeparableObjectKernel
    dense_convolution = similar(unbinned_rate)
    separable_convolution = similar(unbinned_rate)
    convolution_scratch = similar(unbinned_rate)
    WavefrontSensors.conv2d_same!(dense_convolution, unbinned_rate, sep_kernel)
    WavefrontSensors.conv2d_same_separable!(separable_convolution,
        convolution_scratch, unbinned_rate, sep_forward.plan.object_kernel.row,
        sep_forward.plan.object_kernel.col)
    @test separable_convolution ≈ dense_convolution rtol=1e-6 atol=1e-6

    malformed_forward = PreparedLiFTForward(forward.plan, forward.workspace,
        forward.output.values, forward.output, forward.backend, forward.device)
    malformed_snapshot = copy(forward.output.values)
    @test_throws InvalidConfiguration evaluate_lift_forward!(malformed_forward)
    @test forward.output.values == malformed_snapshot
    plan_aliased_forward = PreparedLiFTForward(forward.plan, forward.workspace,
        forward.plan.diversity_opd, forward.output, forward.backend,
        forward.device)
    @test_throws InvalidConfiguration evaluate_lift_forward!(plan_aliased_forward)
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
    @test_throws DimensionMismatchError evaluate_lift_forward!(
        malformed_workspace_forward)

    zernike = ZernikeBasis(tel, 4)
    compute_zernike!(zernike, tel)
    adaptive_basis = copy(@view zernike.modes[:, :, 2:3])
    adaptive_diversity = 50e-9 .* @view(zernike.modes[:, :, 4])
    adaptive_truth = [10e-9, -5e-9]
    adaptive_opd = adaptive_diversity .+
        adaptive_truth[1] .* @view(adaptive_basis[:, :, 1]) .+
        adaptive_truth[2] .* @view(adaptive_basis[:, :, 2])
    adaptive_forward = prepare_lift_forward_model(tel, src, adaptive_basis,
        adaptive_opd; diversity_opd=adaptive_diversity, focal_resolution=8)
    adaptive_rate = copy(intensity_values(evaluate_lift_forward!(
        adaptive_forward)))

    aoc_model = LiFTForwardModel(adaptive_forward)
    @test aoc_model isa AOCPhaseRetrieval.AbstractLiFTForwardModel{Float64}
    @test fieldnames(typeof(aoc_model)) == (:plan,)
    @test aoc_model.plan === lift_forward_plan(adaptive_forward)
    @test isconcretetype(typeof(aoc_model))
    @test AOCPhaseRetrieval.coefficient_count(aoc_model) == 2
    @test AOCPhaseRetrieval.observation_axes(aoc_model) == axes(adaptive_rate)
    aoc_model_workspace = AOCPhaseRetrieval.allocate_model_workspace(aoc_model)
    aoc_prediction = AOCPhaseRetrieval.allocate_photon_rate(aoc_model)
    @test (@inferred AOCPhaseRetrieval.predict_photon_rate!(aoc_prediction,
        aoc_model, aoc_model_workspace, adaptive_truth)) === aoc_prediction
    @test aoc_prediction ≈ adaptive_rate rtol=2e-3 atol=eps(Float64)

    common = (; iterations=3,
        solve_mode=AOCPhaseRetrieval.LiFTSolveNormalEquations(),
        damping=AOCPhaseRetrieval.LiFTAdaptiveLevenbergMarquardt(),
        mode_indices=(1, 2),
        model_scaling=AOCPhaseRetrieval.LiFTPhysicalRatePreservation(),
        check_convergence=false)
    rate_domain = AOCPhaseRetrieval.LiFTPhotonRate()
    analytic = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
        rate_domain; common...,
        jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian())
    @test isconcretetype(typeof(analytic.workspace))
    @test isconcretetype(typeof(analytic.workspace.model_workspace))
    @test (@inferred process_aoc_lift!(analytic)) === analytic.result
    @test AOCPhaseRetrieval.lift_coefficients(analytic.result) ≈
        adaptive_truth rtol=2e-3 atol=eps(Float64)
    analytic_diagnostics = AOCPhaseRetrieval.lift_diagnostics(analytic.result)
    @test analytic_diagnostics.iterations == 3
    @test analytic_diagnostics.regularization >= 0
    @test isfinite(analytic_diagnostics.residual_norm)
    @test isfinite(analytic_diagnostics.weighted_residual_norm)

    numerical = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
        rate_domain; common...,
        jacobian_method=AOCPhaseRetrieval.LiFTNumericalJacobian(1e-9))
    @test process_aoc_lift!(numerical) === numerical.result
    @test AOCPhaseRetrieval.lift_coefficients(numerical.result) ≈
        adaptive_truth rtol=2e-3 atol=eps(Float64)
    @test AOCPhaseRetrieval.lift_coefficients(numerical.result) ≈
        AOCPhaseRetrieval.lift_coefficients(analytic.result) rtol=2e-3 atol=eps(Float64)

    selected_options = merge(common, (; mode_indices=(2, 1),
        jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian()))
    selected_order = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
        rate_domain; selected_options...)
    process_aoc_lift!(selected_order)
    @test AOCPhaseRetrieval.lift_coefficients(selected_order.result) ≈
        reverse(adaptive_truth) rtol=2e-3 atol=eps(Float64)

    initial = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
        rate_domain; iterations=1,
        solve_mode=AOCPhaseRetrieval.LiFTSolveNormalEquations(),
        model_scaling=AOCPhaseRetrieval.LiFTPhysicalRatePreservation(),
        check_convergence=true, initial_coefficients=copy(adaptive_truth))
    process_aoc_lift!(initial)
    @test AOCPhaseRetrieval.lift_coefficients(initial.result) ≈
        adaptive_truth rtol=2e-3 atol=eps(Float64)
    @test AOCPhaseRetrieval.lift_diagnostics(initial.result).converged

    exposure_duration = 0.25
    quantum_efficiency = 0.8
    readout_noise = 2e-3
    aos_rate_domain = LiFTPhotonRate(
        noise_equivalent_exposure_s=exposure_duration,
        quantum_efficiency=quantum_efficiency)
    aos_count_domain = LiFTExpectedCounts(exposure_duration;
        quantum_efficiency=quantum_efficiency)
    aos_normalized_domain = LiFTNormalizedIntensity(sum(adaptive_rate);
        noise_equivalent_exposure_s=exposure_duration,
        quantum_efficiency=quantum_efficiency)
    aos_rate_observation = LiFTObservation(adaptive_forward, copy(adaptive_rate);
        domain=aos_rate_domain, readout_noise_std=readout_noise)
    aos_count_observation = LiFTObservation(adaptive_forward,
        adaptive_rate .* (exposure_duration * quantum_efficiency);
        domain=aos_count_domain, readout_noise_std=readout_noise)
    aos_normalized_observation = LiFTObservation(adaptive_forward,
        adaptive_rate ./ sum(adaptive_rate); domain=aos_normalized_domain,
        readout_noise_std=readout_noise)
    observation_options = merge(common, (;
        jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian()))
    rate_from_observation = process_aoc_observation(adaptive_forward,
        aos_rate_observation; observation_options...)
    count_from_observation = process_aoc_observation(adaptive_forward,
        aos_count_observation; observation_options...)
    normalized_from_observation = process_aoc_observation(adaptive_forward,
        aos_normalized_observation; observation_options...)
    @test rate_from_observation.specification.observation_domain isa
        AOCPhaseRetrieval.LiFTPhotonRate
    @test count_from_observation.specification.observation_domain isa
        AOCPhaseRetrieval.LiFTExpectedCounts
    @test normalized_from_observation.specification.observation_domain isa
        AOCPhaseRetrieval.LiFTNormalizedIntensity
    @test rate_from_observation.specification.read_noise_std == readout_noise
    @test rate_from_observation.specification.observation_domain.noise_equivalent_exposure_s ==
        exposure_duration
    @test count_from_observation.specification.observation_domain.exposure_duration_s ==
        exposure_duration
    @test normalized_from_observation.specification.observation_domain.photon_rate_per_unit ==
        sum(adaptive_rate)
    for prepared in (rate_from_observation, count_from_observation,
        normalized_from_observation)
        @test AOCPhaseRetrieval.lift_coefficients(prepared.result) ≈
            AOCPhaseRetrieval.lift_coefficients(analytic.result) rtol=2e-3 atol=eps(Float64)
    end

    mismatched_source = Source(band=:H, magnitude=0.0)
    wavelength_forward = prepare_lift_forward_model(tel, mismatched_source,
        adaptive_basis, adaptive_opd; diversity_opd=adaptive_diversity,
        focal_resolution=8)
    wavelength_observation = LiFTObservation(wavelength_forward,
        copy(intensity_values(evaluate_lift_forward!(wavelength_forward))))
    @test_throws InvalidConfiguration AOCPhaseRetrieval.LiFTSpecification(
        adaptive_forward, wavelength_observation)
    preprocessing_forward = prepare_lift_forward_model(tel, src,
        adaptive_basis, adaptive_opd; diversity_opd=adaptive_diversity,
        focal_resolution=8, mapping=LiFTFrameMapping())
    preprocessing_observation = LiFTObservation(preprocessing_forward,
        copy(intensity_values(evaluate_lift_forward!(preprocessing_forward))))
    @test_throws InvalidConfiguration AOCPhaseRetrieval.LiFTSpecification(
        adaptive_forward, preprocessing_observation)
    malformed_observation = LiFTObservation(aos_rate_observation.metadata,
        zeros(7, 8))
    @test_throws DimensionMismatchError AOCPhaseRetrieval.LiFTSpecification(
        adaptive_forward, malformed_observation)

    for scaling in (
        AOCPhaseRetrieval.LiFTTotalRateMatching(),
        AOCPhaseRetrieval.LiFTPeakRateMatching(),
        AOCPhaseRetrieval.LiFTPhysicalRatePreservation(),
    )
        scaled = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
            rate_domain; iterations=2,
            jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian(),
            solve_mode=AOCPhaseRetrieval.LiFTSolveNormalEquations(),
            model_scaling=scaling, check_convergence=false)
        process_aoc_lift!(scaled)
        @test all(isfinite, AOCPhaseRetrieval.lift_coefficients(scaled.result))
    end

    for weighting in (
        AOCPhaseRetrieval.LiFTInitialModelWeighting(),
        AOCPhaseRetrieval.LiFTIterativeModelWeighting(),
        AOCPhaseRetrieval.LiFTReadNoiseWeighting(),
        AOCPhaseRetrieval.LiFTVarianceMapWeighting(ones(size(adaptive_rate))),
    )
        weighting_options = merge(common, (;
            jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian(),
            weighting=weighting))
        weighted = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
            rate_domain; weighting_options...)
        process_aoc_lift!(weighted)
        @test all(isfinite, AOCPhaseRetrieval.lift_coefficients(weighted.result))
    end
    external_variance = ones(size(adaptive_rate))
    variance_method = AOCPhaseRetrieval.LiFT(iterations=1,
        mode_indices=(1, 2),
        weighting=AOCPhaseRetrieval.LiFTVarianceMapWeighting(external_variance))
    variance_plan = AdaptiveOpticsCalibration.prepare(variance_method,
        AOCPhaseRetrieval.LiFTSpecification(aoc_model, rate_domain))
    fill!(external_variance, 2)
    @test all(isone, variance_plan.weighting.variance)

    for solve_mode in (
        AOCPhaseRetrieval.LiFTSolveAuto(),
        AOCPhaseRetrieval.LiFTSolveQR(),
        AOCPhaseRetrieval.LiFTSolveNormalEquations(),
    )
        solve_options = merge(common, (;
            jacobian_method=AOCPhaseRetrieval.LiFTAnalyticJacobian(),
            solve_mode=solve_mode))
        solved = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
            rate_domain; solve_options...)
        process_aoc_lift!(solved)
        diagnostics = AOCPhaseRetrieval.lift_diagnostics(solved.result)
        @test all(isfinite, AOCPhaseRetrieval.lift_coefficients(solved.result))
        @test diagnostics.regularization >= 0
        @test diagnostics.used_qr == !(solve_mode isa
            AOCPhaseRetrieval.LiFTSolveNormalEquations)
    end
    for damping in (
        AOCPhaseRetrieval.LiFTDampingNone(),
        AOCPhaseRetrieval.LiFTLevenbergMarquardt(),
        AOCPhaseRetrieval.LiFTAdaptiveLevenbergMarquardt(),
    )
        damping_options = merge(common, (;
            jacobian_method=AOCPhaseRetrieval.LiFTNumericalJacobian(1e-9),
            damping=damping))
        damped = prepare_aoc_lift(adaptive_forward, copy(adaptive_rate),
            rate_domain; damping_options...)
        process_aoc_lift!(damped)
        @test all(isfinite, AOCPhaseRetrieval.lift_coefficients(damped.result))
        @test AOCPhaseRetrieval.lift_diagnostics(damped.result).regularization >= 0
    end

    aliased_coefficients = @view numerical.workspace.model_workspace.opd[1:2]
    @test_throws InvalidConfiguration AOCPhaseRetrieval.predict_photon_rate!(
        aoc_prediction, aoc_model, numerical.workspace.model_workspace,
        aliased_coefficients)
    invalid_jacobian = similar(numerical.workspace.jacobian)
    @test_throws DimensionMismatchError AOCPhaseRetrieval.analytic_photon_rate_jacobian!(
        invalid_jacobian, aoc_model, numerical.workspace.model_workspace,
        adaptive_truth, [1, 3])
    @test_throws ArgumentError AdaptiveOpticsCalibration.prepare(
        AOCPhaseRetrieval.LiFT(mode_indices=(1, 1)),
        AOCPhaseRetrieval.LiFTSpecification(aoc_model, rate_domain))

    process_aoc_lift!(analytic)
    result_snapshot = copy(AOCPhaseRetrieval.lift_coefficients(analytic.result))
    diagnostics_snapshot = AOCPhaseRetrieval.lift_diagnostics(analytic.result)
    invalid_observation = copy(adaptive_rate)
    invalid_observation[1] = -1
    invalid_inputs = AOCPhaseRetrieval.LiFTInputs(invalid_observation)
    @test_throws ArgumentError AdaptiveOpticsCalibration.process!(analytic.result,
        analytic.workspace, analytic.plan, invalid_inputs)
    @test AOCPhaseRetrieval.lift_coefficients(analytic.result) == result_snapshot
    @test AOCPhaseRetrieval.lift_diagnostics(analytic.result) == diagnostics_snapshot

    if coverage_instrumented()
        @test_skip "AdaptiveOpticsCalibration LiFT allocation assertion is disabled under coverage instrumentation"
    else
        process_aoc_lift!(analytic)
        @test @allocated(process_aoc_lift!(analytic)) == 0
    end

    tel32 = Telescope(resolution=8, diameter=8f0, central_obstruction=0f0,
        T=Float32)
    src32 = Source(band=:I, magnitude=0f0, T=Float32)
    basis32 = rand(MersenneTwister(29), Float32, 8, 8, 3)
    diversity32 = zeros(Float32, 8, 8)
    forward32 = prepare_lift_forward_model(tel32, src32, basis32, diversity32;
        diversity_opd=diversity32, focal_resolution=8)
    rate32 = copy(intensity_values(evaluate_lift_forward!(forward32)))
    float32_inverse = prepare_aoc_lift(forward32, rate32,
        AOCPhaseRetrieval.LiFTPhotonRate(); iterations=2,
        jacobian_method=AOCPhaseRetrieval.LiFTNumericalJacobian(1f-9),
        solve_mode=AOCPhaseRetrieval.LiFTSolveNormalEquations(),
        mode_indices=(1, 2), check_convergence=false)
    process_aoc_lift!(float32_inverse)
    @test eltype(AOCPhaseRetrieval.lift_coefficients(
        float32_inverse.result)) === Float32
    @test all(isfinite, AOCPhaseRetrieval.lift_coefficients(
        float32_inverse.result))

    @test_throws InvalidConfiguration prepare_lift_forward_model(
        tel, src, basis, diversity; diversity_opd=diversity,
        focal_resolution=8, object_kernel=zeros(3, 3))
    @test_throws InvalidConfiguration prepare_lift_forward_model(
        tel, src, basis, diversity; diversity_opd=diversity,
        focal_resolution=8, object_kernel=fill(-1.0, 3, 3))
end
