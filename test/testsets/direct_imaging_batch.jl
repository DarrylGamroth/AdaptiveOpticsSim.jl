function direct_imaging_batch_reference(
    inputs,
    sources;
    zero_padding::Int,
)
    stages = map(eachindex(sources)) do index
        prepare_direct_imaging(
            inputs[index],
            sources[index];
            zero_padding=zero_padding,
        )
    end
    products = map(stage -> form_direct_image!(stage), stages)
    return stages, products
end

function direct_imaging_batch_allocation_bytes(prepared)
    form_direct_image!(prepared)
    validation_bytes = @allocated(
        AdaptiveOpticsSim.Optics.validate_direct_imaging_batch(prepared),
    )
    render_bytes = @allocated form_direct_image!(prepared)
    return (; validation_bytes, render_bytes)
end

function direct_imaging_batch_close(actual, expected)
    T = eltype(expected)
    scale = maximum(abs, expected)
    tolerance = T(32) * eps(T)
    return isapprox(
        actual,
        expected;
        rtol=tolerance,
        atol=tolerance * scale,
    )
end

mutable struct DirectImagingCountingFFTPlan{P}
    plan::P
    submissions::Int
end

function AdaptiveOpticsSim.Backends.execute_fft_plan!(
    buffer,
    plan::DirectImagingCountingFFTPlan,
)
    plan.submissions += 1
    return AdaptiveOpticsSim.Backends.execute_fft_plan!(buffer, plan.plan)
end

@testset "Prepared direct-imaging batch parity and signature" begin
    T = Float64
    n = 12
    zero_padding = 2
    tel = Telescope(
        resolution=n,
        diameter=T(6),
        central_obstruction=T(0.1),
        T=T,
    )
    first_pupil = PupilFunction(tel; T=T)
    second_pupil = PupilFunction(first_pupil)
    first_opd = reshape(
        collect(range(T(-12e-9), T(12e-9); length=n * n)),
        n,
        n,
    )
    second_opd = reverse(first_opd; dims=1)
    apply_opd!(first_pupil, first_opd)
    apply_opd!(second_pupil, second_opd)

    wavelength_m = T(800e-9)
    output_sample_arcsec =
        T(180 * 3600 / pi) * wavelength_m / T(tel.params.diameter) /
        T(zero_padding)
    on_axis = Source(
        band=:custom,
        wavelength=wavelength_m,
        photon_irradiance=T(2),
        T=T,
    )
    off_axis = Source(
        band=:custom,
        wavelength=wavelength_m,
        photon_irradiance=T(0.75),
        coordinates=(T(25) * output_sample_arcsec, T(0)),
        T=T,
    )
    inputs = [first_pupil, second_pupil]
    sources = [on_axis, off_axis]
    prepared = prepare_direct_imaging_batch(
        inputs,
        sources;
        zero_padding=zero_padding,
    )

    @test prepared isa AdaptiveOpticsSim.Optics.PreparedDirectImagingBatch
    @test direct_imaging_output(prepared) === prepared.output
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_products(prepared) ===
        prepared.output
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_count(prepared) == 2
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_inputs(prepared)[1] ===
        first_pupil
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_inputs(prepared)[2] ===
        second_pupil
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_sources(prepared)[1] ===
        on_axis
    copied_membership =
        AdaptiveOpticsSim.Optics.direct_imaging_batch_sources(prepared)._storage
    copied_membership[1] = off_axis
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_sources(prepared)[1] ===
        on_axis
    copied_sources =
        copy(AdaptiveOpticsSim.Optics.direct_imaging_batch_sources(prepared))
    copied_sources[1] = off_axis
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_sources(prepared)[1] ===
        on_axis
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_capability(
        FraunhoferPropagation,
    ) isa AdaptiveOpticsSim.Optics.StackedFraunhoferDirectImagingBatchCapability
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_capability(
        FresnelPropagation,
    ) isa AdaptiveOpticsSim.Optics.UnsupportedDirectImagingBatchCapability

    signature = AdaptiveOpticsSim.Optics.direct_imaging_batch_signature(prepared)
    @test signature.model_type === FraunhoferPropagation
    @test signature.sample_count == 2
    @test signature.fft_dimensions == (1, 2)
    @test signature.padded_dimensions ==
        (n * zero_padding, n * zero_padding)
    @test signature.numeric_type === T
    @test signature.backend isa CPUBackend
    @test signature.device == AdaptiveOpticsSim.Backends.HostComputeDevice()
    @test signature.radiometry isa PhysicalPhotonIrradianceSource
    @test signature.product_contract.kind isa FocalPlane
    @test signature.product_contract.coordinate_domain isa AngularCoordinates
    @test signature.product_contract.normalization isa PhotonRateNormalization
    @test signature.product_contract.spatial_measure isa CellIntegratedMeasure
    @test signature.product_contract.coherence isa
        IncoherentIntensityAddition
    @test signature.samples[1].wavelength_m == wavelength_m
    @test signature.samples[2].shift_samples == (25, 0)
    @test prepared.workspace.fft_plan ===
        prepared.workspace_bindings.fft_plan
    @test parent(prepared.fields[1].values) ===
        prepared.workspace.field_stack
    @test parent(prepared.output[2].values) ===
        prepared.workspace.output_stack
    @test compute_device(prepared.output[1].values) ==
        AdaptiveOpticsSim.Backends.HostComputeDevice()

    formed = @inferred form_direct_image!(prepared)
    @test formed === prepared.output
    serial_stages, references = direct_imaging_batch_reference(
        inputs,
        sources;
        zero_padding=zero_padding,
    )
    @test length(serial_stages) == 2
    @test length(formed) == length(references)
    @test direct_imaging_batch_close(
        formed[1].values,
        references[1].values,
    )
    @test direct_imaging_batch_close(
        formed[2].values,
        references[2].values,
    )
    @test formed[1].metadata == references[1].metadata
    @test formed[2].metadata == references[2].metadata

    counting_plan = DirectImagingCountingFFTPlan(
        prepared.workspace.fft_plan,
        0,
    )
    counted_workspace = AdaptiveOpticsSim.Optics.DirectImagingBatchWorkspace(
        prepared.workspace.field_stack,
        prepared.workspace.output_stack,
        prepared.workspace.shift_axis1,
        prepared.workspace.shift_axis2,
        counting_plan,
    )
    counted_bindings =
        AdaptiveOpticsSim.Optics.DirectImagingBatchWorkspaceBindings(
            counted_workspace.field_stack,
            counted_workspace.output_stack,
            counted_workspace.shift_axis1,
            counted_workspace.shift_axis2,
            counting_plan,
            prepared.fields,
            prepared.output,
        )
    counted = AdaptiveOpticsSim.Optics.PreparedDirectImagingBatch(
        prepared.signature,
        prepared.sources,
        prepared.inputs,
        prepared.fields,
        prepared.formation_plans,
        prepared.output,
        counted_workspace,
        counted_bindings,
    )
    form_direct_image!(counted)
    @test counting_plan.submissions == 1
    form_direct_image!(counted)
    @test counting_plan.submissions == 2

    allocation_bytes = direct_imaging_batch_allocation_bytes(prepared)
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test allocation_bytes.validation_bytes == 0
        @test allocation_bytes.render_bytes == 0
    end

    apply_opd!(first_pupil, T(1.5) .* first_opd)
    apply_opd!(second_pupil, T(0.5) .* second_opd)
    form_direct_image!(prepared)
    _, updated_references = direct_imaging_batch_reference(
        inputs,
        sources;
        zero_padding=zero_padding,
    )
    @test direct_imaging_batch_close(
        prepared.output[1].values,
        updated_references[1].values,
    )
    @test direct_imaging_batch_close(
        prepared.output[2].values,
        updated_references[2].values,
    )

    shared = prepare_direct_imaging_batch(
        first_pupil,
        Asterism(sources);
        zero_padding=zero_padding,
    )
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_inputs(shared)[1] ===
        first_pupil
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_inputs(shared)[2] ===
        first_pupil
    form_direct_image!(shared)
    _, shared_references = direct_imaging_batch_reference(
        [first_pupil, first_pupil],
        sources;
        zero_padding=zero_padding,
    )
    @test direct_imaging_batch_close(
        shared.output[1].values,
        shared_references[1].values,
    )
    @test direct_imaging_batch_close(
        shared.output[2].values,
        shared_references[2].values,
    )

    single = prepare_direct_imaging_batch(
        first_pupil,
        on_axis;
        zero_padding=zero_padding,
    )
    @test AdaptiveOpticsSim.Optics.direct_imaging_batch_count(single) == 1
    @test typeof(signature.samples) ===
        typeof(AdaptiveOpticsSim.Optics.direct_imaging_batch_signature(single).samples)
    form_direct_image!(single)
    single_reference = prepare_direct_imaging(
        first_pupil,
        on_axis;
        zero_padding=zero_padding,
    )
    form_direct_image!(single_reference)
    @test direct_imaging_batch_close(
        single.output[1].values,
        direct_imaging_output(single_reference).values,
    )
end

@testset "Direct-imaging spectral batch and detector fan-out" begin
    T = Float64
    tel = Telescope(
        resolution=10,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    pupil = PupilFunction(tel; T=T)
    source = Source(
        band=:custom,
        wavelength=T(800e-9),
        photon_irradiance=T(4),
        coordinates=(T(0.04), T(35)),
        T=T,
    )
    wavelengths = T[650e-9, 800e-9, 950e-9]
    weights = T[0.2, 0.3, 0.5]
    spectral_source = with_spectrum(
        source,
        SpectralBundle(wavelengths, weights; T=T),
    )
    prepared = prepare_direct_imaging_batch(
        pupil,
        spectral_source;
        zero_padding=2,
    )
    products = form_direct_image!(prepared)
    serial = prepare_direct_imaging(
        pupil,
        spectral_source;
        zero_padding=2,
    )
    serial_products = form_direct_image!(serial)

    @test products isa OpticalProductBundle
    @test length(products) == length(wavelengths)
    @test all(
        index -> direct_imaging_batch_close(
            products[index].values,
            serial_products[index].values,
        ),
        eachindex(wavelengths),
    )
    @test all(index -> products[index].metadata ==
        serial_products[index].metadata, eachindex(wavelengths))
    @test [products[index].metadata.spectral.wavelength_m
        for index in eachindex(wavelengths)] == wavelengths
    @test allunique([
        products[index].metadata.sampling
        for index in eachindex(wavelengths)
    ])
    @test isapprox(
        sum(sum(product.values) for product in products),
        sum(pupil_photon_rate_map(tel, source)),
    )

    selected = products[2]
    selected_before = copy(selected.values)
    short = Detector(
        integration_time=T(0.25),
        noise=NoiseNone(),
        qe=T(0.5),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(),
        T=T,
    )
    long = Detector(
        integration_time=T(1.5),
        noise=NoiseNone(),
        qe=T(0.5),
        response_model=NullFrameResponse(),
        sensor=CCDSensor(),
        T=T,
    )
    short_plan = prepare_detector_acquisition(short, selected)
    long_plan = prepare_detector_acquisition(long, selected)
    @test short_plan.input_values === selected.values
    @test long_plan.input_values === selected.values
    short_frame = copy(capture!(
        short,
        selected,
        short_plan;
        rng=MersenneTwister(501),
    ))
    long_frame = copy(capture!(
        long,
        selected,
        long_plan;
        rng=MersenneTwister(502),
    ))
    @test long_frame ≈ T(6) .* short_frame
    @test selected.values == selected_before
end

@testset "Direct-imaging batch rejection and nonmutation" begin
    T = Float64
    n = 8
    tel = Telescope(
        resolution=n,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    pupil = PupilFunction(tel; T=T)
    source = Source(
        band=:custom,
        wavelength=T(750e-9),
        photon_irradiance=one(T),
        T=T,
    )
    prepared = prepare_direct_imaging_batch(pupil, source)
    form_direct_image!(prepared)
    @test AdaptiveOpticsSim.Optics.validate_direct_imaging_batch(prepared) ===
        prepared

    fill!(prepared.workspace.output_stack, T(37))
    saved_shifts = prepared.workspace.shift_axis1
    prepared.workspace.shift_axis1 = zeros(Int, 2)
    @test_throws InvalidConfiguration form_direct_image!(prepared)
    @test all(==(T(37)), prepared.workspace.output_stack)
    prepared.workspace.shift_axis1 = saved_shifts

    foreign_output = ContractDeviceArray(
        fill(T(41), size(prepared.workspace.output_stack)),
        ContractComputeDevice(7),
    )
    foreign_workspace = AdaptiveOpticsSim.Optics.DirectImagingBatchWorkspace(
        prepared.workspace.field_stack,
        foreign_output,
        prepared.workspace.shift_axis1,
        prepared.workspace.shift_axis2,
        prepared.workspace.fft_plan,
    )
    foreign_bindings =
        AdaptiveOpticsSim.Optics.DirectImagingBatchWorkspaceBindings(
            prepared.workspace.field_stack,
            foreign_output,
            prepared.workspace.shift_axis1,
            prepared.workspace.shift_axis2,
            prepared.workspace.fft_plan,
            prepared.fields,
            prepared.output,
        )
    foreign_prepared = AdaptiveOpticsSim.Optics.PreparedDirectImagingBatch(
        prepared.signature,
        prepared.sources,
        prepared.inputs,
        prepared.fields,
        prepared.formation_plans,
        prepared.output,
        foreign_workspace,
        foreign_bindings,
    )
    @test_throws InvalidConfiguration form_direct_image!(foreign_prepared)
    @test all(==(T(41)), foreign_output.storage)

    incompatible_fft_signature =
        AdaptiveOpticsSim.Optics.DirectImagingBatchCompatibilitySignature(
            prepared.signature.model_type,
            prepared.signature.capability,
            prepared.signature.input_metadata,
            prepared.signature.radiometry,
            prepared.signature.backend,
            prepared.signature.device,
            prepared.signature.product_contract,
            prepared.signature.samples,
            prepared.signature.padded_dimensions,
            (2, 3),
            prepared.signature.numeric_type,
            prepared.signature.sample_count,
        )
    incompatible_fft_prepared =
        AdaptiveOpticsSim.Optics.PreparedDirectImagingBatch(
            incompatible_fft_signature,
            prepared.sources,
            prepared.inputs,
            prepared.fields,
            prepared.formation_plans,
            prepared.output,
            prepared.workspace,
            prepared.workspace_bindings,
        )
    fill!(prepared.workspace.output_stack, T(43))
    @test_throws InvalidConfiguration form_direct_image!(
        incompatible_fft_prepared,
    )
    @test all(==(T(43)), prepared.workspace.output_stack)

    incompatible_model_signature =
        AdaptiveOpticsSim.Optics.DirectImagingBatchCompatibilitySignature(
            FresnelPropagation,
            prepared.signature.capability,
            prepared.signature.input_metadata,
            prepared.signature.radiometry,
            prepared.signature.backend,
            prepared.signature.device,
            prepared.signature.product_contract,
            prepared.signature.samples,
            prepared.signature.padded_dimensions,
            prepared.signature.fft_dimensions,
            prepared.signature.numeric_type,
            prepared.signature.sample_count,
        )
    incompatible_model_prepared =
        AdaptiveOpticsSim.Optics.PreparedDirectImagingBatch(
            incompatible_model_signature,
            prepared.sources,
            prepared.inputs,
            prepared.fields,
            prepared.formation_plans,
            prepared.output,
            prepared.workspace,
            prepared.workspace_bindings,
        )
    fill!(prepared.workspace.output_stack, T(47))
    @test_throws InvalidConfiguration form_direct_image!(
        incompatible_model_prepared,
    )
    @test all(==(T(47)), prepared.workspace.output_stack)

    @test_throws DimensionMismatchError prepare_direct_imaging_batch(
        [pupil, PupilFunction(pupil)],
        [source],
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        PupilFunction[],
        Source[],
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        pupil,
        Asterism(Source[]),
    )
    different_wavelength = Source(
        band=:custom,
        wavelength=T(900e-9),
        photon_irradiance=one(T),
        T=T,
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        pupil,
        Asterism([source, different_wavelength]),
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        pupil,
        source;
        zero_padding=0,
    )

    normalized_source = Source(
        band=:custom,
        wavelength=T(750e-9),
        normalized_power=one(T),
        T=T,
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        pupil,
        normalized_source,
    )
    lgs = LGSSource(
        wavelength=T(589e-9),
        photon_irradiance=one(T),
        T=T,
    )
    @test_throws UnsupportedAlgorithm prepare_direct_imaging_batch(
        pupil,
        Asterism(AdaptiveOpticsSim.Optics.AbstractSource[lgs]),
    )
    @test_throws UnsupportedAlgorithm prepare_direct_imaging_batch(
        pupil,
        lgs,
    )
    @test_throws UnsupportedAlgorithm prepare_direct_imaging_batch(
        FresnelPropagation,
        pupil,
        source,
    )

    different_grid = PupilFunction(Telescope(
        resolution=n,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    ))
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        [pupil, different_grid],
        [source, source],
    )
    float32_pupil = PupilFunction(Telescope(
        resolution=n,
        diameter=Float32(4),
        central_obstruction=Float32(0),
        T=Float32,
    ))
    mixed_inputs = PupilFunction[pupil, float32_pupil]
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        mixed_inputs,
        AdaptiveOpticsSim.Optics.AbstractSource[source, source],
    )
    float32_source = Source(
        band=:custom,
        wavelength=Float32(750e-9),
        photon_irradiance=one(Float32),
        T=Float32,
    )
    @test_throws InvalidConfiguration prepare_direct_imaging_batch(
        [pupil, pupil],
        AdaptiveOpticsSim.Optics.AbstractSource[source, float32_source],
    )
end
