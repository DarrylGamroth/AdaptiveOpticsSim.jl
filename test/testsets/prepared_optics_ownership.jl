function pe03_optics_hot_path_allocations(direct, spatial)
    form_direct_image!(direct)
    filter!(spatial)
    return @allocated begin
        form_direct_image!(direct)
        filter!(spatial)
    end
end

@testset "PE-03 prepared optics ownership" begin
    T = Float64
    tel = Telescope(
        resolution=8,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    src = Source(
        band=:custom,
        wavelength=T(800e-9),
        photon_irradiance=T(2),
        T=T,
    )
    pupil = PupilFunction(tel; T=T)
    apply_opd!(pupil, reshape(
        collect(range(T(-4e-9), T(4e-9); length=64)),
        8,
        8,
    ))
    field = ElectricField(pupil, src; zero_padding=2)
    formation = prepare_pupil_field(pupil, src, field)
    fill_electric_field!(field, pupil, formation)

    fraunhofer = FraunhoferPropagation(field)
    fraunhofer_plan = Optics.propagation_plan(fraunhofer)
    replacement_fraunhofer_workspace =
        Optics.FraunhoferPropagationWorkspace(fraunhofer_plan, field.values)
    first_fraunhofer = similar(field.values)
    second_fraunhofer = similar(field.values)
    propagate_field!(first_fraunhofer, field, fraunhofer_plan,
        Optics.propagation_workspace(fraunhofer))
    propagate_field!(second_fraunhofer, field, fraunhofer_plan,
        replacement_fraunhofer_workspace)
    @test first_fraunhofer ≈ second_fraunhofer
    @test replacement_fraunhofer_workspace !==
        Optics.propagation_workspace(fraunhofer)

    fresnel = FresnelPropagation(field; distance_m=T(12))
    fresnel_plan = Optics.propagation_plan(fresnel)
    replacement_fresnel_workspace =
        Optics.FresnelPropagationWorkspace(fresnel_plan, field.values)
    first_fresnel = similar(field.values)
    second_fresnel = similar(field.values)
    propagate_field!(first_fresnel, field, fresnel_plan,
        Optics.propagation_workspace(fresnel))
    propagate_field!(second_fresnel, field, fresnel_plan,
        replacement_fresnel_workspace)
    @test first_fresnel ≈ second_fresnel
    @test fresnel_plan.transfer !== field.values
    @test fresnel_plan.transfer !== first_fresnel
    @test replacement_fresnel_workspace !==
        Optics.propagation_workspace(fresnel)

    output = IntensityMap(field, fraunhofer)
    direct = prepare_direct_imaging(pupil, src, field, output)
    direct_plan = Optics.direct_imaging_plan(direct)
    direct_workspace = Optics.direct_imaging_workspace(direct)
    @test all(name -> !(getfield(direct_plan, name) isa AbstractArray),
        fieldnames(typeof(direct_plan)))
    @test direct_plan.input_plan.formation === formation
    @test direct_plan.field_metadata === field.metadata
    @test direct_plan.output_metadata === output.metadata
    @test direct_workspace.propagation !== direct_plan.propagation

    replacement_direct_workspace = Optics.DirectImagingWorkspace(
        Optics.FraunhoferPropagationWorkspace(
            direct_plan.propagation, field.values),
        similar(output.values),
    )
    replacement_direct = Optics.PreparedDirectImaging(
        Optics._PREPARED_DIRECT_IMAGING_TOKEN,
        pupil,
        field,
        output,
        direct_plan,
        replacement_direct_workspace,
    )
    first_direct = copy(intensity_values(form_direct_image!(direct)))
    fill!(output.values, zero(T))
    second_direct = copy(intensity_values(form_direct_image!(replacement_direct)))
    @test first_direct ≈ second_direct

    fill!(output.values, T(-3))
    @test_throws InvalidConfiguration Optics._form_prepared_direct_image!(
        direct,
        pupil,
        field,
        output,
        replacement_direct_workspace,
    )
    @test all(==(T(-3)), output.values)

    spatial_filter = SpatialFilter(
        tel;
        shape=SquareFilter(),
        diameter=T(4),
        zero_padding=2,
        T=T,
    )
    spatial_field = ElectricField(
        pupil,
        src;
        zero_padding=2,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination(),
    )
    spatial_formation = prepare_pupil_field(
        pupil,
        src,
        spatial_field;
        center_even_grid=false,
        amplitude_scale=one(T),
    )
    fill_electric_field!(spatial_field, pupil, spatial_formation)
    spatial_output = PupilFunction(tel; T=T)
    spatial = prepare_spatial_filter(
        tel, spatial_filter, spatial_field, spatial_output)
    spatial_plan = Optics.spatial_filter_plan(spatial)
    @test all(name -> !(getfield(spatial_plan, name) isa AbstractArray),
        fieldnames(typeof(spatial_plan)))

    foreign_aperture = Telescope(
        resolution=8,
        diameter=T(4),
        central_obstruction=T(0.5),
        T=T,
    )
    foreign_spatial_output = PupilFunction(foreign_aperture; T=T)
    foreign_amplitude_before = copy(foreign_spatial_output.amplitude)
    foreign_opd_before = copy(foreign_spatial_output.opd)
    @test Optics.aperture_revision(foreign_spatial_output) ==
        Optics.aperture_revision(tel)
    @test pupil_support(foreign_spatial_output) != pupil_mask(tel)
    @test_throws InvalidConfiguration prepare_spatial_filter(
        tel, spatial_filter, spatial_field, foreign_spatial_output)
    @test foreign_spatial_output.amplitude == foreign_amplitude_before
    @test foreign_spatial_output.opd == foreign_opd_before

    first_spatial = begin
        filter!(spatial)
        (copy(spatial_output.amplitude), copy(spatial_output.opd))
    end
    replacement_spatial_workspace =
        Optics.SpatialFilterWorkspace(spatial_filter)
    replacement_spatial = Optics.PreparedSpatialFilter(
        Optics._PREPARED_SPATIAL_FILTER_TOKEN,
        spatial_filter,
        spatial_field,
        spatial_output,
        spatial_plan,
        replacement_spatial_workspace,
    )
    fill!(spatial_output.amplitude, zero(T))
    fill!(spatial_output.opd, zero(T))
    filter!(replacement_spatial)
    @test spatial_output.amplitude ≈ first_spatial[1]
    @test spatial_output.opd ≈ first_spatial[2]

    @test @inferred(form_direct_image!(direct)) === output
    @test @inferred(filter!(spatial)) === spatial_output
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test pe03_optics_hot_path_allocations(direct, spatial) == 0
    end

    off_axis = Source(
        band=:custom,
        wavelength=wavelength(src),
        photon_irradiance=one(T),
        coordinates=(T(0.01), T(0)),
        T=T,
    )
    one_source = prepare_direct_imaging_batch(pupil, src; zero_padding=2)
    two_sources = prepare_direct_imaging_batch(
        pupil,
        Asterism([src, off_axis]);
        zero_padding=2,
    )
    @test one_source.sources isa FixedSizeVector
    @test two_sources.sources isa FixedSizeVector
    @test isconcretetype(eltype(two_sources.sources))
    @test typeof(one_source.sources) === typeof(two_sources.sources)
    @test_throws MethodError push!(two_sources.sources, src)
    @test_throws MethodError resize!(two_sources.sources, 1)
end
