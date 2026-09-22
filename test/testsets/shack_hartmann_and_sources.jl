@testset "Shack-Hartmann physical subaperture geometry" begin
    telescope = Telescope(resolution=352, diameter=1.22)
    geometric = ShackHartmannWFS(telescope; n_lenslets=16, T=Float32)
    flux = ShackHartmannWFS(telescope;
        n_lenslets=16,
        valid_subaperture_policy=
            FluxThresholdValidSubapertures(light_ratio=0.5f0, T=Float32),
        T=Float32)

    geometric_mask = copy(geometric.front_end.layout.valid_mask_host)
    flux_mask = copy(flux.front_end.layout.valid_mask_host)
    @test count(geometric_mask) == 216
    @test count(flux_mask) == 208
    @test flux_mask != geometric_mask
    @test isempty(findall(flux_mask .& .!geometric_mask))

    sensor = ShackHartmannWFS(
        Telescope(resolution=4, diameter=2.0);
        n_lenslets=2,
        n_pix_subap=2,
        T=Float32,
    )
    layout = sensor.front_end.layout
    revision = subaperture_layout_revision(layout)
    @test_throws DimensionMismatchError set_valid_subapertures!(sensor,
        fill(true, 1, 1))
    parent = fill(false, 4, 4)
    parent[2:3, 2:3] .= Bool[true false; false true]
    @test @inferred(set_valid_subapertures!(sensor, @view(parent[2:3, 2:3]))) === sensor
    @test subaperture_layout_revision(layout) == revision + UInt(1)
    @test n_valid_subapertures(layout) == 2
    @test valid_subaperture_indices(layout) == CartesianIndex{2}[
        CartesianIndex(1, 1), CartesianIndex(2, 2)]
end

@testset "Shack-Hartmann physical optics and acquisition" begin
    T = Float64
    telescope = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    optics = shack_hartmann_optics(sensor, source)
    rate = shack_hartmann_rate_map(optics, pupil)
    prepared = @inferred prepare_wfs_optics(optics, pupil, rate)
    @test prepared isa WavefrontSensors.PreparedShackHartmannOptics
    @test @inferred(form_wfs_optical_products!(rate, pupil, prepared)) === rate
    @test all(isfinite, rate.values)
    @test sum(rate.values) > zero(T)
    @test rate.metadata.normalization isa PhotonRateNormalization
    @test rate.metadata.spatial_measure isa CellIntegratedMeasure
    @test rate.metadata.coherence isa IncoherentIntensityAddition

    detector = Detector(noise=NoiseNone(), exposure_duration=T(0.25),
        qe=one(T), response_model=NullFrameResponse(), T=T)
    observation = WFSObservation(similar(rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    acquisition = @inferred prepare_wfs_acquisition(detector, rate,
        observation)
    @test @inferred(acquire_wfs_observation!(observation, rate,
        acquisition, Xoshiro(0x309))) === observation
    @test observation.storage ≈ rate.values .* T(0.25) atol=0 rtol=0
    @test_throws WFSPreparationError prepare_wfs_estimation(sensor,
        observation,
        WFSMeasurement(zeros(T, 2); units=:pixel, kind=:centroid_slopes))

    if !coverage_instrumented()
        @test @allocated(form_wfs_optical_products!(rate, pupil,
            prepared)) == 0
        rng = Xoshiro(0x309)
        acquire_wfs_observation!(observation, rate, acquisition, rng)
        @test @allocated(acquire_wfs_observation!(observation, rate,
            acquisition, rng)) == 0
    end
end

@testset "Shack-Hartmann physical source composition" begin
    T = Float64
    telescope = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(6), T=T)

    asterism = Asterism([
        source,
        Source(band=:custom, wavelength=wavelength(source),
            photon_irradiance=T(4), coordinates=(T(0.1), T(-0.05)), T=T),
    ])
    asterism_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    asterism_rate = shack_hartmann_rate_map(asterism_sensor, pupil,
        asterism)
    asterism_plan = prepare_wfs_optics(
        shack_hartmann_optics(asterism_sensor, asterism), pupil,
        asterism_rate)
    form_wfs_optical_products!(asterism_rate, pupil, asterism_plan)
    @test all(isfinite, asterism_rate.values)
    @test supports_stacked_sources(asterism_sensor, asterism)
    @test supports_grouped_execution(asterism_sensor, asterism)

    lgs = LGSSource(wavelength=wavelength(source), photon_irradiance=T(6),
        elongation_factor=T(1.8), T=T)
    lgs_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    lgs_rate = shack_hartmann_rate_map(lgs_sensor, pupil, lgs)
    lgs_plan = prepare_wfs_optics(shack_hartmann_optics(lgs_sensor, lgs),
        pupil, lgs_rate)
    form_wfs_optical_products!(lgs_rate, pupil, lgs_plan)
    @test all(isfinite, lgs_rate.values)
    @test sum(lgs_rate.values) > zero(T)

    spectral = with_spectrum(source, SpectralBundle(
        T[0.9 * wavelength(source), 1.1 * wavelength(source)],
        T[0.4, 0.6]; T=T))
    spectral_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    spectral_rates = shack_hartmann_rate_map(spectral_sensor, pupil,
        spectral)
    @test spectral_rates isa OpticalProductBundle
    @test length(spectral_rates) == 2
    spectral_plan = prepare_wfs_optics(
        shack_hartmann_optics(spectral_sensor, spectral), pupil,
        spectral_rates)
    form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
    @test all(all(isfinite, product.values) for product in spectral_rates)

    extended = with_extended_source(source,
        PointCloudSourceModel([(0.0, 0.0), (0.2, 0.0)], [0.5, 0.5]))
    extended_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    extended_rate = shack_hartmann_rate_map(extended_sensor, pupil,
        extended)
    extended_plan = prepare_wfs_optics(
        shack_hartmann_optics(extended_sensor, extended), pupil,
        extended_rate)
    form_wfs_optical_products!(extended_rate, pupil, extended_plan)
    @test all(isfinite, extended_rate.values)
    @test sum(extended_rate.values) > zero(T)
end

@testset "Shack-Hartmann geometric wavefront truth" begin
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    layout_sensor = ShackHartmannWFS(telescope; n_lenslets=2,
        n_pix_subap=4, T=T)
    pupil.opd .= reshape(T.(1:64), 8, 8) .* T(1e-9)
    truth = zeros(T, 8)
    layout = layout_sensor.front_end.layout
    @test @inferred(geometric_wavefront_slopes!(truth, pupil.opd,
        layout.valid_mask, pupil.metadata.sampling)) === truth
    @test all(isfinite, truth)
    @test any(!iszero, truth)
    if !coverage_instrumented()
        @test @allocated(geometric_wavefront_slopes!(truth, pupil.opd,
            layout.valid_mask, pupil.metadata.sampling)) == 0
    end
    @test !applicable(measure!, layout_sensor, pupil, Source(band=:I,
        magnitude=0.0))
    @test !applicable(slopes, layout_sensor)
end
