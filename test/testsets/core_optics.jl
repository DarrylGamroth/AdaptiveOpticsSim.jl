@testset "GPU backend registry" begin
    @test !gpu_backend_loaded(AdaptiveOpticsSim.CUDABackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.MetalBackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.AMDGPUBackendTag)
    @test gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.MetalBackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag) === nothing
    @test gpu_backend_name(Matrix{Float64}) === nothing
    @test available_gpu_backends() == ()
    @test AdaptiveOpticsSim.GPUArrayBuildBackend(AdaptiveOpticsSim.CUDABackendTag) isa AdaptiveOpticsSim.GPUArrayBuildBackend{AdaptiveOpticsSim.CUDABackendTag}
end

@testset "API export curation" begin
    exported = names(AdaptiveOpticsSim)
    @test length(exported) <= 300
    @test Base.isexported(AdaptiveOpticsSim, :Telescope)
    @test Base.isexported(AdaptiveOpticsSim, :ShackHartmannWFS)
    @test Base.isexported(AdaptiveOpticsSim, :Detector)
    @test Base.isexported(AdaptiveOpticsSim, :ControlLoopScenario)
    @test Base.isexported(AdaptiveOpticsSim, :DeformableMirror)
    @test Base.isexported(AdaptiveOpticsSim, :PyramidWFS)
    @test Base.isexported(AdaptiveOpticsSim, :BioEdgeWFS)
    @test Base.isexported(AdaptiveOpticsSim, :ZernikeWFS)
    @test Base.isexported(AdaptiveOpticsSim, :CurvatureWFS)
    @test Base.isexported(AdaptiveOpticsSim, :influence_model)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_runtime_wfs!)
    @test Base.isexported(AdaptiveOpticsSim, :subaperture_layout)
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
    @test !Base.isexported(AdaptiveOpticsSim, :CUDABackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :MetalBackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :AMDGPUBackendTag)
    @test AdaptiveOpticsSim.CPUBuildBackend() isa AdaptiveOpticsSim.BuildBackend
end

@testset "Telescope and PSF" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    field = ElectricField(tel, src; zero_padding=2)
    @test size(field.state.field) == (64, 64)
    @test field.params.wavelength == wavelength(src)
    centered_psf = AdaptiveOpticsSim.centered_psf_from_field!(similar(field.state.intensity), field)
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
    @test size(fmap) == size(tel_dim.state.pupil)
    @test maximum(fmap) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"

    tel_simple = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src_simple = Source(band=:I, magnitude=0.0)
    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    direct = similar(field_simple.state.field)
    fill_telescope_field!(direct, tel_simple, src_simple; zero_padding=1)
    @test direct == field_simple.state.field

    phase_map = fill(pi / 2, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_phase!(field_simple, phase_map; units=:phase)
    @test field_simple.state.field ≈ baseline .* cis.(phase_map)

    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    opd_map = fill(eltype(tel_simple.state.opd)(wavelength(src_simple) / 4), tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_phase!(field_simple, opd_map; units=:opd)
    @test field_simple.state.field ≈ baseline .* cispi.(2 .* opd_map ./ wavelength(src_simple))

    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    amplitude_map = fill(0.5, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_amplitude!(field_simple, amplitude_map)
    @test field_simple.state.field ≈ baseline .* amplitude_map

    intensity_buffer = similar(field_simple.state.intensity)
    intensity!(intensity_buffer, field_simple)
    @test intensity_buffer ≈ abs2.(field_simple.state.field)
    @test intensity!(field_simple) === field_simple.state.intensity

    fraunhofer = FraunhoferPropagation(field)
    propagated = similar(field.state.field)
    propagate_field!(propagated, field, fraunhofer)
    propagated_psf = similar(field.state.intensity)
    @. propagated_psf = abs2(propagated)
    @test propagated_psf ≈ centered_psf atol=1e-10 rtol=1e-10
    @test fraunhofer.params.output_sampling_rad ≈ wavelength(src) / (field.params.padded_resolution * field.params.sampling_m)

    fresnel = FresnelPropagation(field; distance_m=25.0)
    fresnel_out = similar(field.state.field)
    propagate_field!(fresnel_out, field, fresnel)
    @test size(fresnel_out) == size(field.state.field)
    @test isfinite(sum(abs, fresnel_out))
    reverse_fresnel = FresnelPropagation(field; distance_m=-25.0)
    reverse_input = ElectricField(tel, src; zero_padding=2)
    copyto!(reverse_input.state.field, fresnel_out)
    propagate_field!(reverse_input, reverse_fresnel)
    @test reverse_input.state.field ≈ field.state.field atol=1e-9 rtol=1e-9

    zero_fresnel = FresnelPropagation(field; distance_m=0.0)
    zero_out = similar(field.state.field)
    propagate_field!(zero_out, field, zero_fresnel)
    @test zero_out ≈ field.state.field atol=1e-10 rtol=1e-10

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
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.state.field),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldSynchronousPlan
    geom_field = propagate_atmosphere_field!(geom_prop, atm, atm_tel, atm_src)
    tel_geom = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    propagate!(atm, tel_geom, atm_src)
    collapsed = ElectricField(tel_geom, atm_src; zero_padding=1)
    @test geom_field.state.field ≈ collapsed.state.field atol=1e-8 rtol=1e-8

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
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.state.field),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldSynchronousPlan
    fresnel_field = propagate_atmosphere_field!(fresnel_prop, fresnel_atm, atm_tel, atm_src)
    geom_single = AtmosphericFieldPropagation(fresnel_atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    geom_single_field = propagate_atmosphere_field!(geom_single, fresnel_atm, atm_tel, atm_src)
    @test fresnel_field.state.field ≈ geom_single_field.state.field atol=1e-8 rtol=1e-8

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
    @test tel.state.pupil == expected
    apply_spiders!(tel; thickness=0.4, angles=[0.0, 90.0])
    manual = copy(expected)
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=0.0))
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=pi / 2))
    @test tel.state.pupil == manual

    sf_tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    sf = SpatialFilter(sf_tel; shape=SquareFilter(), diameter=4, zero_padding=2)
    @test count(x -> !iszero(x), sf.state.mask) > 0
    foucault = SpatialFilter(sf_tel; shape=FoucaultFilter(), diameter=4, zero_padding=2)
    foucault_count = count(x -> !iszero(x), foucault.state.mask)
    @test 0 < foucault_count < length(foucault.state.mask)

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
