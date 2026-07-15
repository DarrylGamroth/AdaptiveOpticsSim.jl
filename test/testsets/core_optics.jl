@testset "GPU backend registry" begin
    @test !gpu_backend_loaded(AdaptiveOpticsSim.CUDABackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.MetalBackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.AMDGPUBackendTag)
    @test gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.MetalBackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag) === nothing
    @test gpu_backend_name(Matrix{Float64}) === nothing
    @test backend(zeros(2, 2)) isa CPUBackend
    @test backend(zeros(2)) isa CPUBackend
    @test available_gpu_backends() == ()
    @test AdaptiveOpticsSim.GPUArrayBuildBackend(AdaptiveOpticsSim.CUDABackendTag) isa AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.CUDABackendTag}
end

function explicit_psf_cycle!(output, field, wavefront, plan, workspace)
    compute_psf!(output, field, wavefront, plan, workspace)
    return output
end

function explicit_spatial_filter_cycle!(output, field, spatial_filter, plan,
    workspace)
    filter!(output, field, spatial_filter, plan, workspace)
    return output
end

@testset "API export curation" begin
    exported = names(AdaptiveOpticsSim)
    @test length(exported) <= 370
    @test Base.isexported(AdaptiveOpticsSim, :Telescope)
    @test Base.isexported(AdaptiveOpticsSim, :ShackHartmannWFS)
    @test Base.isexported(AdaptiveOpticsSim, :Detector)
    @test Base.isexported(AdaptiveOpticsSim, :MKIDArrayDetector)
    @test Base.isexported(AdaptiveOpticsSim, :CMOSReadNoiseMap)
    @test Base.isexported(AdaptiveOpticsSim, :SkipperSampling)
    @test Base.isexported(AdaptiveOpticsSim, :GlobalShutter)
    @test Base.isexported(AdaptiveOpticsSim, :SingleRead)
    @test Base.isexported(AdaptiveOpticsSim, :AveragedNonDestructiveReads)
    @test Base.isexported(AdaptiveOpticsSim, :UpTheRampSampling)
    @test Base.isexported(AdaptiveOpticsSim, :FrameTransferAcquisition)
    @test Base.isexported(AdaptiveOpticsSim, :detector_ramp_cube)
    @test Base.isexported(AdaptiveOpticsSim, :detector_ramp_times)
    @test Base.isexported(AdaptiveOpticsSim, :InterpixelCapacitance)
    @test Base.isexported(AdaptiveOpticsSim, :ControlLoopScenario)
    @test Base.isexported(AdaptiveOpticsSim, :SimulationEnsemble)
    @test Base.isexported(AdaptiveOpticsSim, :FactorizedReconstructor)
    @test Base.isexported(AdaptiveOpticsSim, :ControlledReconstructor)
    @test Base.isexported(AdaptiveOpticsSim, :DeterministicExecution)
    @test Base.isexported(AdaptiveOpticsSim, :AcceleratedKernelsExecution)
    @test Base.isexported(AdaptiveOpticsSim, :DaggerExecution)
    @test Base.isexported(AdaptiveOpticsSim, :DeformableMirror)
    @test Base.isexported(AdaptiveOpticsSim, :PyramidWFS)
    @test Base.isexported(AdaptiveOpticsSim, :BioEdgeWFS)
    @test Base.isexported(AdaptiveOpticsSim, :ZernikeWFS)
    @test Base.isexported(AdaptiveOpticsSim, :CurvatureWFS)
    @test Base.isexported(AdaptiveOpticsSim, :influence_model)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_runtime_wfs!)
    @test Base.isexported(AdaptiveOpticsSim, :wfs_source)
    @test Base.isexported(AdaptiveOpticsSim, :science_source)
    @test Base.isexported(AdaptiveOpticsSim, :CPUHILExecutionPlan)
    @test Base.isexported(AdaptiveOpticsSim, :DeviceResidentExecutionPlan)
    @test Base.isexported(AdaptiveOpticsSim, :runtime_execution_plan)
    @test Base.isexported(AdaptiveOpticsSim, :synchronize_runtime!)
    @test Base.isexported(AdaptiveOpticsSim, :OpticalWFSChannel)
    @test Base.isexported(AdaptiveOpticsSim, :SharedOpticalArm)
    @test Base.isexported(AdaptiveOpticsSim, :SharedOpticalRuntime)
    @test Base.isexported(AdaptiveOpticsSim, :optical_arms)
    @test Base.isexported(AdaptiveOpticsSim, :subaperture_layout)
    @test Base.isexported(AdaptiveOpticsSim, :OpticalPlaneMetadata)
    @test Base.isexported(AdaptiveOpticsSim, :PupilWavefront)
    @test Base.isexported(AdaptiveOpticsSim, :ElectricField)
    @test Base.isexported(AdaptiveOpticsSim, :IrradiancePlane)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_pupil_field)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_direct_psf)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_spatial_filter)
    @test !Base.isexported(AdaptiveOpticsSim, :TelescopeParams)
    @test !Base.isexported(AdaptiveOpticsSim, :TelescopeState)
    @test !Base.isexported(AdaptiveOpticsSim, :DetectorParams)
    @test !Base.isexported(AdaptiveOpticsSim, :DetectorState)
    @test !Base.isexported(AdaptiveOpticsSim, :supports_detector_mtf)
    @test !Base.isexported(AdaptiveOpticsSim, :supports_prepared_runtime)
    @test !Base.isexported(AdaptiveOpticsSim, :backend_rand)
    @test !Base.isexported(AdaptiveOpticsSim, :backend_zeros)
    @test !Base.isexported(AdaptiveOpticsSim, :BuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :CPUBuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :GPUArrayBuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :default_build_backend)
    @test !Base.isexported(AdaptiveOpticsSim, :set_fft_provider_threads!)
    @test !Base.isexported(AdaptiveOpticsSim, :GPUBackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :AbstractRuntimeExecutionPlan)
    @test !Base.isexported(AdaptiveOpticsSim, :runtime_reconstructor_storage)
    @test !Base.isexported(AdaptiveOpticsSim, :ensemble_ownership_roots)
    @test !Base.isexported(AdaptiveOpticsSim, :run_ensemble!)
    @test !Base.isexported(AdaptiveOpticsSim, :CUDABackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :MetalBackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :AMDGPUBackendTag)
    @test AdaptiveOpticsSim.CPUBuildBackend() isa AdaptiveOpticsSim.BuildBackend
end

@testset "Telescope and PSF" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilWavefront(tel)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=2)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test size(field.values) == (64, 64)
    @test field.metadata.spectral == MonochromaticChannel(wavelength(src))
    fraunhofer = FraunhoferPropagation(field)
    centered_psf = similar(field.values, Float64)
    AdaptiveOpticsSim.centered_psf_from_field!(centered_psf, field,
        fraunhofer)
    psf = compute_psf!(tel, src; zero_padding=2)
    @test size(psf) == (64, 64)
    @test maximum(psf) > 0
    @test isfinite(sum(psf))
    @test centered_psf ≈ psf
    @test tel.state.psf_workspace !== nothing
    cached_ws = tel.state.psf_workspace
    compute_psf!(tel, src; zero_padding=2)
    @test tel.state.psf_workspace === cached_ws

    tel_dim = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2,
        pupil_reflectivity=0.25)
    psf_dim = compute_psf!(tel_dim, src; zero_padding=2)
    @test sum(psf_dim) ≈ 0.25 * sum(psf)
    fmap = flux_map(tel_dim, src)
    @test size(fmap) == size(pupil_mask(tel_dim))
    @test maximum(fmap) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"

    tel_simple = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src_simple = Source(band=:I, magnitude=0.0)
    wavefront_simple = PupilWavefront(tel_simple)
    field_simple = ElectricField(wavefront_simple, src_simple;
        zero_padding=1)
    formation_simple = prepare_pupil_field(tel_simple, wavefront_simple,
        src_simple, field_simple)
    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    direct = similar(field_simple.values)
    fill_telescope_field!(direct, tel_simple, src_simple; zero_padding=1)
    @test direct == field_simple.values

    phase_map = fill(pi / 2, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_phase!(field_simple, phase_map; units=:phase)
    @test field_simple.values ≈ baseline .* cis.(phase_map)

    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    quarter_wave_opd = fill(eltype(tel_simple.state.opd)(wavelength(src_simple) / 4), tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_phase!(field_simple, quarter_wave_opd; units=:opd)
    @test field_simple.values ≈ baseline .* cispi.(2 .* quarter_wave_opd ./ wavelength(src_simple))

    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    amplitude_map = fill(0.5, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_amplitude!(field_simple, amplitude_map)
    @test field_simple.values ≈ baseline .* amplitude_map

    intensity_buffer = similar(field_simple.values, Float64)
    intensity!(intensity_buffer, field_simple)
    @test intensity_buffer ≈ abs2.(field_simple.values)

    propagated = propagation_output(field, fraunhofer)
    propagate_field!(propagated, field, fraunhofer)
    propagated_psf = similar(propagated.values, Float64)
    @. propagated_psf = abs2(propagated.values)
    @test propagated_psf ≈ centered_psf atol=1e-10 rtol=1e-10
    @test fraunhofer.params.output_sampling_rad ≈ wavelength(src) /
        (field.metadata.dimensions[1] * field.metadata.sampling[1])
    @test propagated.metadata.kind isa FocalPlane

    fresnel = FresnelPropagation(field; distance_m=25.0)
    fresnel_out = propagation_output(field, fresnel)
    propagate_field!(fresnel_out, field, fresnel)
    @test size(fresnel_out.values) == size(field.values)
    @test isfinite(sum(abs, fresnel_out.values))
    reverse_fresnel = FresnelPropagation(fresnel_out; distance_m=-25.0,
        output_kind=PupilPlane())
    reverse_output = propagation_output(fresnel_out, reverse_fresnel)
    propagate_field!(reverse_output, fresnel_out, reverse_fresnel)
    @test reverse_output.values ≈ field.values atol=1e-9 rtol=1e-9

    zero_fresnel = FresnelPropagation(field; distance_m=0.0)
    zero_out = propagation_output(field, zero_fresnel)
    propagate_field!(zero_out, field, zero_fresnel)
    @test zero_out.values ≈ field.values atol=1e-10 rtol=1e-10

    atm_tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm_src = Source(band=:I, magnitude=0.0)
    atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, atm_tel; rng=MersenneTwister(1))
    geom_prop = AtmosphericFieldPropagation(atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldSynchronousPlan
    geom_field = propagate_atmosphere_field!(geom_prop, atm, atm_tel, atm_src)
    tel_geom = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    propagate!(atm, tel_geom, atm_src)
    collapsed_wavefront = PupilWavefront(tel_geom)
    apply_opd!(collapsed_wavefront, opd_map(tel_geom))
    collapsed = ElectricField(collapsed_wavefront, atm_src; zero_padding=1)
    collapsed_plan = prepare_pupil_field(tel_geom, collapsed_wavefront,
        atm_src, collapsed)
    fill_electric_field!(collapsed, collapsed_wavefront, collapsed_plan)
    @test geom_field.values ≈ collapsed.values atol=1e-8 rtol=1e-8

    fresnel_atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(fresnel_atm, atm_tel; rng=MersenneTwister(2))
    fresnel_prop = AtmosphericFieldPropagation(fresnel_atm, atm_tel, atm_src;
        model=LayeredFresnelAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldSynchronousPlan
    fresnel_field = propagate_atmosphere_field!(fresnel_prop, fresnel_atm, atm_tel, atm_src)
    geom_single = AtmosphericFieldPropagation(fresnel_atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    geom_single_field = propagate_atmosphere_field!(geom_single, fresnel_atm, atm_tel, atm_src)
    @test fresnel_field.values ≈ geom_single_field.values atol=1e-8 rtol=1e-8

    spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src), 1.1 * wavelength(atm_src)], [0.75, 0.25]))
    spectral_prop = AtmosphericFieldPropagation(atm, atm_tel, spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    spectral_intensity = atmospheric_intensity!(spectral_prop, atm, atm_tel, spectral)
    mono_prop = AtmosphericFieldPropagation(atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    mono_intensity = atmospheric_intensity!(mono_prop, atm, atm_tel, atm_src)
    single_spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src)], [1.0]))
    single_spectral_prop = AtmosphericFieldPropagation(atm, atm_tel, single_spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    @test atmospheric_intensity!(single_spectral_prop, atm, atm_tel, single_spectral) ≈ mono_intensity atol=1e-8 rtol=1e-8
    @test size(spectral_intensity) == size(mono_intensity)
    @test sum(spectral_intensity) > 0
end

@testset "Explicit optical products and surfaces" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    telescope_opd_before = copy(opd_map(tel))
    aperture_before = copy(pupil_mask(tel))
    reflectivity_before = copy(pupil_reflectivity(tel))

    path_a = PupilWavefront(tel)
    path_b = PupilWavefront(tel)
    static_map = OPDMap(fill(2e-9, 16, 16))
    apply_surface!(path_a, static_map, DMAdditive())
    @test path_a.opd == static_map.opd
    @test all(iszero, path_b.opd)
    @test opd_map(tel) == telescope_opd_before

    ncpa = NCPA(fill(3e-9, 16, 16), nothing, nothing)
    apply_surface!(path_b, ncpa, DMReplace())
    @test path_b.opd == ncpa.opd

    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    set_command!(dm, fill(1e-8, n_actuators(dm)))
    update_surface!(dm, tel)
    apply_surface!(path_a, dm, DMReplace())
    @test path_a.opd == surface_opd(dm)

    modal = TipTiltMirror(tel; scale=1e-8)
    set_command!(modal, [0.25, -0.5])
    update_surface!(modal, tel)
    reset_opd!(path_b)
    apply_surface!(path_b, modal, DMAdditive())
    @test path_b.opd == surface_opd(modal)
    @test pupil_mask(tel) == aperture_before
    @test pupil_reflectivity(tel) == reflectivity_before
    @test opd_map(tel) == telescope_opd_before

    reset_opd!(path_a)
    path_b_before_propagation = copy(path_b.opd)
    telescope_psf_before = copy(tel.state.psf)
    field = ElectricField(path_a, src; zero_padding=2)
    prepared = prepare_direct_psf(tel, path_a, src, field)
    explicit_psf_cycle!(prepared.output, field, path_a, prepared.plan,
        prepared.workspace)
    @test @allocated(explicit_psf_cycle!(prepared.output, field, path_a,
        prepared.plan, prepared.workspace)) == 0
    @test path_b.opd == path_b_before_propagation
    @test tel.state.psf == telescope_psf_before

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2)
    spatial_formation = prepare_pupil_field(tel, path_a, src, field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(field, path_a, spatial_formation)
    spatial_output = PupilWavefront(tel)
    spatial_plan = prepare_spatial_filter(tel, spatial_filter, field,
        spatial_output)
    spatial_workspace = SpatialFilterWorkspace(spatial_filter)
    explicit_spatial_filter_cycle!(spatial_output, field, spatial_filter,
        spatial_plan, spatial_workspace)
    @test @allocated(explicit_spatial_filter_cycle!(spatial_output, field,
        spatial_filter, spatial_plan, spatial_workspace)) == 0
end

@testset "Optical-plane compatibility validation" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilWavefront(tel)
    field = ElectricField(wavefront, src; zero_padding=1)
    values = similar(field.values)

    sampling_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=(2.0, 2.0), spectral=field.metadata.spectral)
    sampling_field = ElectricField(sampling_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, sampling_field)

    kind_metadata = OpticalPlaneMetadata(FocalPlane(), values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    kind_field = ElectricField(kind_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, kind_field)

    origin_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=field.metadata.sampling, origin=(0.0, 0.0),
        spectral=field.metadata.spectral)
    origin_field = ElectricField(origin_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, origin_field)

    centering_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        centering=(SampleCentered, SampleCentered),
        spectral=field.metadata.spectral)
    centering_field = ElectricField(centering_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, centering_field)

    orientation_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        orientation=PlaneAxisOrientation((:y, :x)),
        spectral=field.metadata.spectral)
    orientation_field = ElectricField(orientation_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, orientation_field)

    spectral_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=MonochromaticChannel(1.1 * wavelength(src)))
    spectral_field = ElectricField(spectral_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, spectral_field)

    float32_values = Matrix{ComplexF32}(undef, size(values))
    numeric_metadata = OpticalPlaneMetadata(PupilPlane(), float32_values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    numeric_field = ElectricField(numeric_metadata, float32_values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, numeric_field)

    wrong_size_values = Matrix{ComplexF64}(undef, 9, 9)
    dimension_metadata = OpticalPlaneMetadata(PupilPlane(),
        wrong_size_values; sampling=field.metadata.sampling,
        spectral=field.metadata.spectral)
    dimension_field = ElectricField(dimension_metadata, wrong_size_values)
    @test_throws DimensionMismatchError prepare_pupil_field(tel, wavefront,
        src, dimension_field)

    propagation = FraunhoferPropagation(field)
    wrong_destination_values = Matrix{Float64}(undef, 9, 9)
    wrong_destination_metadata = OpticalPlaneMetadata(FocalPlane(),
        wrong_destination_values; sampling=propagation.output_metadata.sampling,
        spectral=propagation.output_metadata.spectral)
    wrong_destination = IrradiancePlane(wrong_destination_metadata,
        wrong_destination_values)
    @test_throws DimensionMismatchError prepare_direct_psf(tel, wavefront,
        src, field, wrong_destination)

    declared_device_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral,
        device=AdaptiveOpticsSim.AcceleratorPlaneDevice(
            KernelAbstractions.CPU(), 1))
    @test_throws InvalidConfiguration ElectricField(
        declared_device_metadata, values)
end

@testset "Aperture masks" begin
    bool_mask = falses(32, 32)
    build_mask!(bool_mask, CircularAperture(radius=1.0))
    @test bool_mask[16, 16]
    @test !bool_mask[1, 1]

    annulus = falses(32, 32)
    build_mask!(annulus, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test !annulus[16, 16]
    @test count(annulus) < count(bool_mask)

    weighted = fill(-1.0, 16, 16)
    build_mask!(weighted, RectangularROI(5:8, 6:10); inside=2.0, outside=-1.0)
    @test all(weighted[5:8, 6:10] .== 2.0)
    @test weighted[4, 6] == -1.0

    spider = trues(32, 32)
    apply_mask!(spider, SpiderMask(thickness=0.08, angle_rad=pi / 2))
    @test count(spider) < length(spider)

    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.25)
    expected = falses(32, 32)
    build_mask!(expected, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test pupil_mask(tel) == expected
    apply_spiders!(tel; thickness=0.4, angles=[0.0, 90.0])
    manual = copy(expected)
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=0.0))
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=pi / 2))
    @test pupil_mask(tel) == manual

    sf_tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    sf = SpatialFilter(sf_tel; shape=SquareFilter(), diameter=4, zero_padding=2)
    @test count(x -> !iszero(x), sf.mask) > 0
    foucault = SpatialFilter(sf_tel; shape=FoucaultFilter(), diameter=4, zero_padding=2)
    foucault_count = count(x -> !iszero(x), foucault.mask)
    @test 0 < foucault_count < length(foucault.mask)

    pupil = falses(8, 8)
    pupil[1:4, 1:4] .= true
    pupil[5:8, 5:8] .= true
    valid = falses(2, 2)
    build_mask!(valid, SubapertureGridMask(threshold=0.5), pupil)
    @test valid == Bool[true false; false true]
    valid2 = similar(valid)
    AdaptiveOpticsSim.set_valid_subapertures!(valid2, pupil, 0.5)
    @test valid2 == valid
end

@testset "Zernike basis" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end
