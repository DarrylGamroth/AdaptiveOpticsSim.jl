function prepared_direct_image_allocations(prepared)
    form_direct_image!(prepared)
    return @allocated form_direct_image!(prepared)
end

@testset "Prepared direct-science formation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.8e-6,
        photon_irradiance=2.0)
    pupil = PupilFunction(tel)
    opd = reshape(collect(range(-5e-9, 5e-9; length=256)), 16, 16)
    apply_opd!(pupil, opd)
    field = ElectricField(pupil, src; zero_padding=2)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    prepared = prepare_direct_imaging(tel, pupil, src, field, output)

    formed = @inferred form_direct_image!(output, pupil, field,
        prepared.plan, prepared.workspace)
    @test formed === output
    @test output.metadata.kind isa FocalPlane
    @test output.metadata.coordinate_domain isa AngularCoordinates
    @test output.metadata.normalization isa PhotonRateNormalization
    @test output.metadata.spatial_measure isa CellIntegratedMeasure
    @test output.metadata.coherence isa IncoherentIntensityAddition
    @test output.metadata.spectral == MonochromaticChannel(wavelength(src))
    @test sum(output.values) ≈ sum(pupil_photon_rate_map(tel, src))
    @test @allocated(form_direct_image!(output, pupil, field,
        prepared.plan, prepared.workspace)) == 0

    initial_output = copy(output.values)
    apply_opd!(pupil, 20 .* opd)
    @test form_direct_image!(output, pupil, field,
        prepared.plan, prepared.workspace) === output
    @test output.values != initial_output
    @test sum(output.values) ≈ sum(pupil_photon_rate_map(tel, src))

    output_before = copy(output.values)
    replacement_values = similar(output.values)
    fill!(replacement_values, -1.0)
    replacement_output = IntensityMap(output.metadata, replacement_values)
    @test_throws InvalidConfiguration form_direct_image!(replacement_output,
        pupil, field, prepared.plan, prepared.workspace)
    @test replacement_output.values == fill(-1.0, size(replacement_values))
    @test output.values == output_before

    replacement_field_values = copy(field.values)
    replacement_field = ElectricField(field.metadata,
        replacement_field_values)
    @test_throws InvalidConfiguration form_direct_image!(output, pupil,
        replacement_field, prepared.plan, prepared.workspace)
    @test output.values == output_before

    replacement_pupil = PupilFunction(tel)
    @test_throws InvalidConfiguration form_direct_image!(output,
        replacement_pupil, field, prepared.plan, prepared.workspace)
    @test output.values == output_before

    other_stage = prepare_direct_imaging(tel, pupil, src;
        zero_padding=2)
    @test_throws InvalidConfiguration form_direct_image!(output, pupil,
        field, prepared.plan, other_stage.workspace)
    @test output.values == output_before

    formation = prepare_pupil_field(tel, pupil, src, field)
    fill_electric_field!(field, pupil, formation)
    preformed_output = IntensityMap(field, FraunhoferPropagation(field))
    preformed = prepare_direct_imaging(src, field, preformed_output)
    @test @inferred(form_direct_image!(preformed_output, field,
        preformed.plan, preformed.workspace)) === preformed_output
    @test preformed_output.values ≈ output.values
    @test @allocated(form_direct_image!(preformed_output, field,
        preformed.plan, preformed.workspace)) == 0

    normalized = Source(band=:custom, wavelength=wavelength(src),
        normalized_power=1.0)
    normalized_field = ElectricField(pupil, normalized; zero_padding=2)
    normalized_output = IntensityMap(normalized_field,
        FraunhoferPropagation(normalized_field))
    @test_throws InvalidConfiguration prepare_direct_imaging(tel, pupil,
        normalized, normalized_field, normalized_output)

    one_sample_arcsec = focal_plane_pixel_scale_arcsec(output)
    positive_x = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=photon_irradiance(src),
        coordinates=(one_sample_arcsec, 0.0))
    positive_y = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=photon_irradiance(src),
        coordinates=(one_sample_arcsec, 90.0))
    on_axis = prepare_direct_imaging(tel, pupil, src; zero_padding=2)
    x_shifted = prepare_direct_imaging(tel, pupil, positive_x;
        zero_padding=2)
    y_shifted = prepare_direct_imaging(tel, pupil, positive_y;
        zero_padding=2)
    on_axis_values = copy(intensity_values(form_direct_image!(on_axis)))
    x_values = intensity_values(form_direct_image!(x_shifted))
    y_values = intensity_values(form_direct_image!(y_shifted))
    @test x_shifted.plan.shift_samples == (1, 0)
    @test y_shifted.plan.shift_samples == (0, 1)
    @test x_values == circshift(on_axis_values, (1, 0))
    @test y_values == circshift(on_axis_values, (0, 1))

    oriented_values = copy(field.values)
    oriented_metadata = OpticalPlaneMetadata(PupilPlane(), oriented_values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling,
        origin=field.metadata.origin,
        centering=field.metadata.centering,
        orientation=PlaneAxisOrientation((:y, :x), (-1, 1)),
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=field.metadata.coherence)
    oriented_field = ElectricField(oriented_metadata, oriented_values)
    oriented_on_axis = prepare_direct_imaging(src, oriented_field)
    oriented_x = prepare_direct_imaging(positive_x, oriented_field)
    oriented_y = prepare_direct_imaging(positive_y, oriented_field)
    oriented_reference = copy(intensity_values(
        form_direct_image!(oriented_on_axis)))
    oriented_x_values = intensity_values(form_direct_image!(oriented_x))
    oriented_y_values = intensity_values(form_direct_image!(oriented_y))
    @test oriented_x.plan.shift_samples == (0, 1)
    @test oriented_y.plan.shift_samples == (-1, 0)
    @test oriented_x_values == circshift(oriented_reference, (0, 1))
    @test oriented_y_values == circshift(oriented_reference, (-1, 0))

    unsupported_orientation = OpticalPlaneMetadata(PupilPlane(),
        oriented_values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling,
        origin=field.metadata.origin,
        centering=field.metadata.centering,
        orientation=PlaneAxisOrientation((:row, :column)),
        spectral=field.metadata.spectral,
        normalization=field.metadata.normalization,
        spatial_measure=field.metadata.spatial_measure,
        coherence=field.metadata.coherence)
    unsupported_field = ElectricField(unsupported_orientation,
        oriented_values)
    @test_throws InvalidConfiguration prepare_direct_imaging(src,
        unsupported_field)

    huge_offset = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=photon_irradiance(src),
        coordinates=(1.0e17, 0.0))
    field_before = copy(field.values)
    output_before_huge_offset = copy(output.values)
    @test_throws InvalidConfiguration prepare_direct_imaging(huge_offset,
        field, output)
    @test field.values == field_before
    @test output.values == output_before_huge_offset
end

@testset "Direct-science source composition" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    wavelength_m = 0.8e-6
    on_axis = Source(band=:custom, wavelength=wavelength_m,
        photon_irradiance=2.0)
    one_pixel_arcsec = (180 * 3600 / pi) * wavelength_m /
        tel.params.diameter / 2
    off_axis = Source(band=:custom, wavelength=wavelength_m,
        photon_irradiance=1.0, coordinates=(one_pixel_arcsec, 0.0))
    pupil = PupilFunction(tel)

    asterism = prepare_direct_imaging(tel, pupil,
        Asterism([on_axis, off_axis]); zero_padding=2)
    @test asterism.components isa Vector
    @test asterism.products isa Vector
    @test isconcretetype(eltype(asterism.components))
    @test isconcretetype(eltype(asterism.products))
    combined = @inferred form_direct_image!(asterism)
    component_maps = asterism.products
    @test length(component_maps) == 2
    @test intensity_values(combined) ≈
        intensity_values(component_maps[1]) .+
        intensity_values(component_maps[2])
    @test prepared_direct_image_allocations(asterism) == 0

    replacement_products = copy(asterism.products)
    first_product = first(replacement_products)
    replacement_products[1] = IntensityMap(first_product.metadata,
        copy(first_product.values))
    combined_before_rejection = copy(combined.values)
    @test_throws InvalidConfiguration accumulate_intensity!(combined,
        replacement_products, asterism.accumulation)
    @test combined.values == combined_before_rejection

    mixed_lgs = LGSSource(wavelength=wavelength_m,
        photon_irradiance=0.5, coordinates=(one_pixel_arcsec, 90.0))
    mixed_sources = AdaptiveOpticsSim.AbstractSource[on_axis, mixed_lgs]
    mixed = prepare_direct_imaging(tel, pupil, Asterism(mixed_sources);
        zero_padding=2)
    @test mixed.components isa Vector
    @test isconcretetype(eltype(mixed.components))
    mixed_output = form_direct_image!(mixed)
    mixed_explicit = zeros(eltype(mixed_output.values),
        size(mixed_output.values))
    for product in mixed.products
        mixed_explicit .+= product.values
    end
    @test mixed_output.values ≈ mixed_explicit
    @test prepared_direct_image_allocations(mixed) == 0

    extended = with_extended_source(on_axis,
        PointCloudSourceModel([(0.0, 0.0), (one_pixel_arcsec, 0.0)],
            [0.25, 0.75]))
    extended_prepared = prepare_direct_imaging(tel, pupil,
        extended_source_asterism(extended); zero_padding=2)
    extended_output = form_direct_image!(extended_prepared)
    @test sum(intensity_values(extended_output)) ≈
        sum(pupil_photon_rate_map(tel, on_axis))
    @test_throws UnsupportedAlgorithm prepare_direct_imaging(tel, pupil,
        extended; zero_padding=2)

    gaussian = with_extended_source(on_axis,
        GaussianDiskSourceModel(sigma_arcsec=one_pixel_arcsec,
            n_side=5))
    gaussian_prepared = prepare_direct_imaging(tel, pupil,
        extended_source_asterism(gaussian); zero_padding=2)
    @test gaussian_prepared.components isa Vector
    @test gaussian_prepared.products isa Vector
    @test isconcretetype(eltype(gaussian_prepared.components))
    @test isconcretetype(eltype(gaussian_prepared.products))
    @test length(gaussian_prepared.components) == 25
    gaussian_output = form_direct_image!(gaussian_prepared)
    gaussian_explicit = zeros(eltype(gaussian_output.values),
        size(gaussian_output.values))
    for product in gaussian_prepared.products
        gaussian_explicit .+= product.values
    end
    @test gaussian_output.values ≈ gaussian_explicit
    @test sum(gaussian_output.values) ≈
        sum(pupil_photon_rate_map(tel, on_axis))
    @test prepared_direct_image_allocations(gaussian_prepared) == 0

    tiny_tel = Telescope(resolution=4, diameter=1.0,
        central_obstruction=0.0)
    tiny_source = Source(band=:custom, wavelength=wavelength_m,
        photon_irradiance=1.0)
    dense_model = GaussianDiskSourceModel(sigma_arcsec=0.1, n_side=11)
    dense_source = with_extended_source(tiny_source, dense_model)
    dense_prepared = prepare_direct_imaging(tiny_tel,
        PupilFunction(tiny_tel), extended_source_asterism(dense_source);
        zero_padding=1)
    @test dense_prepared.components isa Vector
    @test isconcretetype(eltype(dense_prepared.components))
    @test length(dense_prepared.components) == 121
    dense_output = form_direct_image!(dense_prepared)
    @test all(isfinite, dense_output.values)
    @test sum(dense_output.values) ≈
        sum(pupil_photon_rate_map(tiny_tel, tiny_source))
    @test prepared_direct_image_allocations(dense_prepared) == 0

    zero_source = Source(band=:custom, wavelength=wavelength_m,
        photon_irradiance=0.0)
    zero_off_axis = Source(band=:custom, wavelength=wavelength_m,
        photon_irradiance=0.0,
        coordinates=(one_pixel_arcsec, 0.0))
    zero_asterism = prepare_direct_imaging(tel, pupil,
        Asterism([zero_source, zero_off_axis]); zero_padding=2)
    @test all(iszero, form_direct_image!(zero_asterism).values)
    @test prepared_direct_image_allocations(zero_asterism) == 0

    spectral_source = with_spectrum(on_axis,
        SpectralBundle([0.7e-6, 0.9e-6], [0.4, 0.6]))
    spectral = prepare_direct_imaging(tel, pupil, spectral_source;
        zero_padding=2)
    @test spectral.components isa Vector
    @test isconcretetype(eltype(spectral.components))
    products = @inferred form_direct_image!(spectral)
    @test products isa OpticalProductBundle
    @test products.products isa Vector
    @test isconcretetype(eltype(products.products))
    @test length(products) == 2
    @test products[1].metadata.sampling != products[2].metadata.sampling
    sum_values = similar(products[1].values)
    sum_output = IntensityMap(products[1].metadata, sum_values)
    @test_throws InvalidConfiguration prepare_incoherent_sum(sum_output,
        products[1], products[2])
    @test sum(sum(intensity_values(product)) for product in products) ≈
        sum(pupil_photon_rate_map(tel, on_axis))
    @test prepared_direct_image_allocations(spectral) == 0

    zero_spectral = prepare_direct_imaging(tel, pupil,
        with_spectrum(zero_source,
            SpectralBundle([0.7e-6, 0.9e-6], [0.4, 0.6]));
        zero_padding=2)
    @test all(product -> all(iszero, product.values),
        form_direct_image!(zero_spectral))
    @test prepared_direct_image_allocations(zero_spectral) == 0

    float32_tel = Telescope(resolution=8, diameter=0.1f0,
        central_obstruction=0.0f0, T=Float32)
    representable_total = 2.0 * Float64(floatmax(Float32))
    high_flux = Source(band=:custom, wavelength=0.8e-6,
        photon_irradiance=representable_total)
    high_flux_spectral = with_spectrum(high_flux,
        SpectralBundle(Float32[0.7e-6, 0.9e-6], Float32[0.5, 0.5];
            T=Float32))
    high_flux_prepared = prepare_direct_imaging(float32_tel,
        PupilFunction(float32_tel), high_flux_spectral; zero_padding=1)
    high_flux_products = form_direct_image!(high_flux_prepared)
    @test all(product -> all(isfinite, intensity_values(product)),
        high_flux_products)

    unrepresentable_total = 4.0 * Float64(floatmax(Float32))
    excessive_flux = Source(band=:custom, wavelength=0.8e-6,
        photon_irradiance=unrepresentable_total)
    excessive_spectral = with_spectrum(excessive_flux,
        SpectralBundle(Float32[0.7e-6, 0.9e-6], Float32[0.5, 0.5];
            T=Float32))
    @test_throws InvalidConfiguration prepare_direct_imaging(float32_tel,
        PupilFunction(float32_tel), excessive_spectral; zero_padding=1)

    leaf_excessive = Source(band=:custom, wavelength=0.8e-6,
        photon_irradiance=floatmax(Float64))
    @test_throws InvalidConfiguration prepare_direct_imaging(float32_tel,
        PupilFunction(float32_tel), leaf_excessive; zero_padding=1)
end

@testset "Direct-science detector fan-out" begin
    tel = Telescope(resolution=8, diameter=4.0, central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.8e-6,
        photon_irradiance=1.0)
    imaging = prepare_direct_imaging(tel, PupilFunction(tel), src;
        zero_padding=1)
    rate_map = form_direct_image!(imaging)
    rate_before = copy(rate_map.values)

    short = Detector(integration_time=0.25, noise=NoiseNone(), qe=0.5,
        response_model=NullFrameResponse(), sensor=CMOSSensor())
    long = Detector(integration_time=1.5, noise=NoiseNone(), qe=0.5,
        response_model=NullFrameResponse())
    short_plan = prepare_detector_acquisition(short, rate_map)
    long_plan = prepare_detector_acquisition(long, rate_map)
    @test short_plan.input_values === rate_map.values
    @test long_plan.input_values === rate_map.values
    short_frame = copy(capture!(short, rate_map, short_plan;
        rng=MersenneTwister(1)))
    long_frame = copy(capture!(long, rate_map, long_plan;
        rng=MersenneTwister(2)))
    @test long_frame ≈ 6 .* short_frame
    @test rate_map.values == rate_before

    response_kernel = [0.0 0.125 0.0; 0.125 0.5 0.125; 0.0 0.125 0.0]
    skipper = Detector(integration_time=0.5, noise=NoiseReadout(0.05), qe=0.25,
        response_model=SampledFrameResponse(response_kernel),
        sensor=CCDSensor(sampling_mode=SkipperSampling(3)))
    skipper_plan = prepare_detector_acquisition(skipper, rate_map)
    @test skipper_plan.input_values === rate_map.values
    skipper_frame = capture!(skipper, rate_map, skipper_plan;
        rng=MersenneTwister(4))
    @test all(isfinite, skipper_frame)
    @test skipper_frame !== short_frame
    @test readout_products(skipper) isa SkipperReadoutProducts
    @test rate_map.values == rate_before

    normalized_values = fill(2.0, 2, 2)
    normalized_metadata = OpticalPlaneMetadata(FocalPlane(),
        normalized_values; coordinate_domain=AngularCoordinates(),
        sampling=(0.1, 0.1),
        spectral=MonochromaticChannel(wavelength(src)),
        normalization=DimensionlessNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    normalized_map = IntensityMap(normalized_metadata, normalized_values)
    external = Detector(integration_time=0.5, noise=NoiseNone(), qe=0.25,
        response_model=NullFrameResponse())
    @test_throws InvalidConfiguration prepare_detector_acquisition(external,
        normalized_map)
    converted = prepare_detector_acquisition(external, normalized_map;
        normalized_to_photon_rate=40.0)
    @test capture!(external, normalized_map, converted;
        rng=MersenneTwister(3)) == fill(10.0, 2, 2)
end
