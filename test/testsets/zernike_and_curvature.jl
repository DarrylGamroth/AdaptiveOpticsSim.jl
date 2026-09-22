@inline wfs_optical_rate_storage(wfs::ZernikeWFS) =
    wfs.acquisition.products.frame
@inline wfs_optical_rate_storage(wfs) =
    wfs.acquisition.products.frame

@testset "Zernike optical front end" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    source = Source(band=:custom, wavelength=0.75e-6, photon_irradiance=10.0)
    wfs = ZernikeWFS(tel; pupil_samples=8, diffraction_padding=2)
    front_end = ZernikeOpticalFrontEnd(wfs, source)
    rate = zernike_rate_map(front_end, pupil)
    plan = prepare_wfs_optics(front_end, pupil, rate)
    form_wfs_optical_products!(rate, pupil, plan)
    @test all(isfinite, rate.values)
    @test all(>=(0.0), rate.values)
    if coverage_instrumented()
        @test_skip "Zernike optical allocation assertion is disabled under coverage instrumentation"
    else
        @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0
    end

    attenuated_tel = Telescope(resolution=32, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    attenuated_pupil = PupilFunction(attenuated_tel)
    attenuated_wfs = ZernikeWFS(attenuated_tel; pupil_samples=8,
        diffraction_padding=2)
    attenuated_front_end = ZernikeOpticalFrontEnd(attenuated_wfs, source)
    attenuated_rate = zernike_rate_map(attenuated_front_end, attenuated_pupil)
    attenuated_plan = prepare_wfs_optics(attenuated_front_end, attenuated_pupil,
        attenuated_rate)
    form_wfs_optical_products!(attenuated_rate, attenuated_pupil, attenuated_plan)
    @test sum(attenuated_rate.values) ≈ 0.25 * sum(rate.values) rtol=1e-12
end

@testset "Zernike and Curvature source-composition support boundaries" begin
    telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    pupil = PupilFunction(telescope)
    source = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    spectral = with_spectrum(source,
        SpectralBundle([0.70e-6, 0.80e-6], [0.5, 0.5]))
    extended = with_extended_source(source,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))

    for expanded in (spectral, extended)
        zernike = ZernikeWFS(telescope; pupil_samples=2)
        front_end = ZernikeOpticalFrontEnd(zernike, expanded)
        rate = zernike_rate_map(front_end, pupil)
        @test_throws UnsupportedAlgorithm prepare_wfs_optics(
            front_end, pupil, rate)
        curvature = CurvatureWFS(telescope; pupil_samples=2)
        curvature_front_end = CurvatureOpticalFrontEnd(curvature, expanded)
        curvature_rates = curvature_rate_maps(curvature_front_end, pupil)
        @test_throws UnsupportedAlgorithm prepare_wfs_optics(
            curvature_front_end, pupil, curvature_rates)
    end
end

@testset "Curvature pupil-reflectivity throughput" begin
    transmission = 0.25
    source = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    full_telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    attenuated_telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=transmission)
    full_pupil = PupilFunction(full_telescope)
    attenuated_pupil = PupilFunction(attenuated_telescope)

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        full_sensor = CurvatureWFS(full_telescope; pupil_samples=2,
            diffraction_padding=2)
        attenuated_sensor = CurvatureWFS(attenuated_telescope;
            pupil_samples=2, diffraction_padding=2)
        curvature_intensity!(style, full_sensor, full_pupil, source)
        curvature_intensity!(style, attenuated_sensor, attenuated_pupil,
            source)
        full_rate = sum(wfs_optical_rate_storage(full_sensor))
        @test full_rate > 0
        @test sum(wfs_optical_rate_storage(attenuated_sensor)) ≈
            transmission * full_rate rtol=1e-12
    end
end

@testset "Curvature diffractive photon-rate conservation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    pupil = PupilFunction(tel)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    expected_two_branch_rate = 2 * sum(pupil_photon_rate_map(tel, src))

    for style in (ScalarCPUStyle(), KA_CPU_STYLE), padding in (1, 2, 3)
        wfs = CurvatureWFS(tel; pupil_samples=2,
            diffraction_padding=padding)
        curvature_intensity!(style, wfs, pupil, src)
        @test sum(wfs.front_end.propagation.workspace.intensity_stack) ≈
            expected_two_branch_rate atol=1e-10 rtol=1e-12
        if padding == 1
            @test sum(wfs_optical_rate_storage(wfs)) ≈
                expected_two_branch_rate atol=1e-10 rtol=1e-12
        else
            @test 0 < sum(wfs_optical_rate_storage(wfs)) <=
                expected_two_branch_rate
        end
    end
end

@testset "Curvature KernelAbstractions CPU stage parity" begin
    T = Float64
    tel = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:64), 8, 8) .* T(1e-10)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)

    scalar_sensor = CurvatureWFS(tel; pupil_samples=2, T=T)
    accelerated_sensor = CurvatureWFS(tel; pupil_samples=2, T=T)
    @test backend(accelerated_sensor) isa CPUBackend
    scalar_front_end = CurvatureOpticalFrontEnd(scalar_sensor, source)
    accelerated_front_end = CurvatureOpticalFrontEnd(
        accelerated_sensor, source)

    scalar_fields = copy(
        WavefrontSensors._form_curvature_branch_fields!(
            ScalarCPUStyle(), scalar_front_end, pupil))
    accelerated_fields = copy(
        WavefrontSensors._form_curvature_branch_fields!(
            KA_CPU_STYLE, accelerated_front_end, pupil))
    @test accelerated_fields ≈ scalar_fields rtol=T(2e-12) atol=T(2e-12)

    scalar_rates = curvature_rate_maps(scalar_front_end, pupil)
    accelerated_rates = curvature_rate_maps(accelerated_front_end, pupil)
    scalar_workspace = scalar_sensor.front_end.propagation.workspace
    accelerated_workspace = accelerated_sensor.front_end.propagation.workspace
    intensity = reshape(T.(1:length(scalar_workspace.intensity_stack)),
        size(scalar_workspace.intensity_stack))
    copyto!(scalar_workspace.intensity_stack, intensity)
    copyto!(accelerated_workspace.intensity_stack, intensity)
    WavefrontSensors._sample_curvature_rate_planes!(
        ScalarCPUStyle(), scalar_rates, scalar_front_end)
    WavefrontSensors._sample_curvature_rate_planes!(
        KA_CPU_STYLE, accelerated_rates, accelerated_front_end)
    @test accelerated_rates[1].values ≈ scalar_rates[1].values
    @test accelerated_rates[2].values ≈ scalar_rates[2].values

    padding = scalar_sensor.front_end.propagation.plan.defocus_pair.diffraction_padding
    field = ElectricField(pupil, source; zero_padding=padding, T=T)
    field_plan = prepare_pupil_field(pupil, source, field;
        center_even_grid=false)
    fill_electric_field!(field, pupil, field_plan)
    scalar_field_stack = copy(
        WavefrontSensors._form_curvature_field_input!(
            ScalarCPUStyle(), scalar_workspace, field))
    accelerated_field_stack = copy(
        WavefrontSensors._form_curvature_field_input!(
            KA_CPU_STYLE, accelerated_workspace, field))
    @test accelerated_field_stack ≈ scalar_field_stack

    accelerated_rates = curvature_rate_maps(accelerated_front_end, pupil)
    optics = prepare_wfs_optics(accelerated_front_end, pupil,
        accelerated_rates)
    @test validate_wfs_target(optics, HostComputeDevice()) === optics
end

@testset "Curvature optical front end" begin
    T = Float64
    telescope = Telescope(resolution=32, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    sensor = CurvatureWFS(telescope; pupil_samples=8,
        defocus_rms_nm=T(500), T=T)

    @test !hasfield(typeof(sensor), :estimator)
    front_end = CurvatureOpticalFrontEnd(sensor, source)
    rates = curvature_rate_maps(front_end, pupil)
    optics = prepare_wfs_optics(front_end, pupil, rates)
    @test validate_wfs_target(optics, HostComputeDevice()) === optics
    form_wfs_optical_products!(rates, pupil, optics)
    @test size(rates[1].values) == (8, 8)
    @test size(rates[2].values) == (8, 8)
    @test all(isfinite, rates[1].values)
    @test all(isfinite, rates[2].values)
    @test all(>=(zero(T)), rates[1].values)
    @test all(>=(zero(T)), rates[2].values)
    if coverage_instrumented()
        @test_skip "Curvature optical allocation assertion is disabled under coverage instrumentation"
    else
        @test @allocated(form_wfs_optical_products!(rates, pupil,
            optics)) == 0
    end

    flat_plus = copy(rates[1].values)
    flat_minus = copy(rates[2].values)
    basis = ZernikeBasis(telescope, 5)
    compute_zernike!(basis, telescope)
    focus = @view basis.modes[:, :, 5]
    @. pupil.opd = T(5e-8) * focus
    form_wfs_optical_products!(rates, pupil, optics)
    positive_plus = copy(rates[1].values)
    positive_minus = copy(rates[2].values)
    @. pupil.opd = -T(5e-8) * focus
    form_wfs_optical_products!(rates, pupil, optics)
    @test norm(rates[1].values - flat_plus) > T(1e-6)
    @test norm(rates[2].values - flat_minus) > T(1e-6)
    @test !(rates[1].values ≈ positive_plus)
    @test !(rates[2].values ≈ positive_minus)
    fill!(pupil.opd, zero(T))

    response = CurvatureBranchResponse(T=T, plus_throughput=T(1.2),
        minus_throughput=T(0.8), plus_background=T(5),
        minus_background=one(T))
    imbalanced = CurvatureWFS(telescope; pupil_samples=8,
        defocus_rms_nm=T(500), branch_response=response, T=T)
    imbalanced_front_end = CurvatureOpticalFrontEnd(imbalanced, source)
    imbalanced_rates = curvature_rate_maps(imbalanced_front_end, pupil)
    imbalanced_optics = prepare_wfs_optics(imbalanced_front_end, pupil,
        imbalanced_rates)
    form_wfs_optical_products!(imbalanced_rates, pupil, imbalanced_optics)
    @test mean(imbalanced_rates[1].values) >
        mean(imbalanced_rates[2].values)
    @test_throws InvalidConfiguration CurvatureBranchResponse(
        plus_throughput=-1.0)

    oversampled = CurvatureWFS(telescope; pupil_samples=8,
        readout_crop_resolution=16, readout_pixels_per_sample=2, T=T)
    oversampled_front_end = CurvatureOpticalFrontEnd(oversampled, source)
    oversampled_rates = curvature_rate_maps(oversampled_front_end, pupil)
    oversampled_optics = prepare_wfs_optics(oversampled_front_end, pupil,
        oversampled_rates)
    form_wfs_optical_products!(oversampled_rates, pupil,
        oversampled_optics)
    @test size.(getfield.(oversampled_rates, :values)) ==
        ((16, 16), (16, 16))
    @test_throws InvalidConfiguration CurvatureWFS(telescope;
        pupil_samples=8, readout_crop_resolution=18,
        readout_pixels_per_sample=2, T=T)
    @test_throws InvalidConfiguration CurvatureWFS(telescope;
        pupil_samples=8, readout_model=CurvatureChannelReadout(),
        readout_pixels_per_sample=2, T=T)
end
