struct TestExpandedSourceWrapper{S<:AdaptiveOpticsSim.AbstractSource} <:
    AdaptiveOpticsSim.AbstractSource
    source::S
end

AdaptiveOpticsSim.source_composition_style(::TestExpandedSourceWrapper) =
    AdaptiveOpticsSim.ExpandedSourceComposition()

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

function explicit_direct_image_cycle!(output, field, wavefront, plan,
    workspace)
    form_direct_image!(output, wavefront, field, plan, workspace)
    return output
end

function explicit_spatial_filter_cycle!(output, field, spatial_filter, plan,
    workspace)
    filter!(output, field, spatial_filter, plan, workspace)
    return output
end

@testset "API export curation" begin
    exported = names(AdaptiveOpticsSim)
    # Pre-HIL stage contracts intentionally add a small, documented public
    # product/protocol seam. Keep headroom bounded so unrelated internals do
    # not drift into the ordinary user namespace.
    # Pre-HIL 9 adds eleven documented component/stage names and no aliases.
    @test length(exported) <= 467
    @test Base.isexported(AdaptiveOpticsSim, :Telescope)
    @test Base.isexported(AdaptiveOpticsSim, :ShackHartmannWFS)
    @test Base.isexported(AdaptiveOpticsSim, :MicrolensArray)
    @test Base.isexported(AdaptiveOpticsSim, :microlens_array)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_microlens_propagation)
    @test !Base.isexported(AdaptiveOpticsSim,
        :PreparedMicrolensPropagation)
    @test Base.isexported(AdaptiveOpticsSim,
        :ShackHartmannDirectFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :ShackHartmannOpticalFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :shack_hartmann_rate_map)
    @test Base.isexported(AdaptiveOpticsSim, :PyramidOpticalFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :BioEdgeOpticalFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :ZernikePhaseSpot)
    @test Base.isexported(AdaptiveOpticsSim, :ZernikeOpticalFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :CurvatureDefocusPair)
    @test Base.isexported(AdaptiveOpticsSim, :CurvatureOpticalFrontEnd)
    @test Base.isexported(AdaptiveOpticsSim, :CurvaturePackedAcquisition)
    @test Base.isexported(AdaptiveOpticsSim, :CircularModulation)
    @test Base.isexported(AdaptiveOpticsSim, :SampledModulation)
    @test !Base.isexported(AdaptiveOpticsSim,
        :PreparedFocalPlaneModulation)
    @test Base.isexported(AdaptiveOpticsSim,
        :set_subaperture_calibration!)
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
    @test Base.isexported(AdaptiveOpticsSim, :photon_irradiance)
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
    @test Base.isexported(AdaptiveOpticsSim, :MetricCoordinates)
    @test Base.isexported(AdaptiveOpticsSim, :AngularCoordinates)
    @test Base.isexported(AdaptiveOpticsSim, :AchromaticSpectralCoordinate)
    @test Base.isexported(AdaptiveOpticsSim, :MonochromaticChannel)
    @test Base.isexported(AdaptiveOpticsSim, :IntegratedSpectralChannel)
    @test Base.isexported(AdaptiveOpticsSim, :PupilFunction)
    @test Base.isexported(AdaptiveOpticsSim, :ElectricField)
    @test Base.isexported(AdaptiveOpticsSim, :IntensityMap)
    @test Base.isexported(AdaptiveOpticsSim, :PhotonRateNormalization)
    @test Base.isexported(AdaptiveOpticsSim, :DimensionlessNormalization)
    @test Base.isexported(AdaptiveOpticsSim, :PointSampledMeasure)
    @test Base.isexported(AdaptiveOpticsSim, :SpatialDensityMeasure)
    @test Base.isexported(AdaptiveOpticsSim, :CellIntegratedMeasure)
    @test Base.isexported(AdaptiveOpticsSim, :CoherentFieldCombination)
    @test Base.isexported(AdaptiveOpticsSim, :IncoherentIntensityAddition)
    @test Base.isexported(AdaptiveOpticsSim, :NonCombinableProduct)
    @test Base.isexported(AdaptiveOpticsSim, :PhysicalPhotonIrradianceSource)
    @test Base.isexported(AdaptiveOpticsSim, :NormalizedTestSource)
    @test Base.isexported(AdaptiveOpticsSim, :source_radiometry)
    @test Base.isexported(AdaptiveOpticsSim, :OpticalProductBundle)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_incoherent_sum)
    @test Base.isexported(AdaptiveOpticsSim, :accumulate_intensity!)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_pupil_field)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_direct_imaging)
    @test Base.isexported(AdaptiveOpticsSim, :form_direct_image!)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_spatial_filter)
    @test Base.isexported(AdaptiveOpticsSim, :AtmosphereEpoch)
    @test Base.isexported(AdaptiveOpticsSim, :advance_by!)
    @test Base.isexported(AdaptiveOpticsSim, :advance_to!)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_atmosphere_renderer)
    @test Base.isexported(AdaptiveOpticsSim, :render_atmosphere!)
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

@testset "Optical radiometry and combination contracts" begin
    for invalid_wavelength in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration MonochromaticChannel(
            invalid_wavelength)
    end
    @test_throws InvalidConfiguration IntegratedSpectralChannel(Symbol(""))

    metadata_storage = zeros(2, 2)
    for invalid_sampling in (
        (0.0, 1.0),
        (-1.0, 1.0),
        (Inf, 1.0),
        (NaN, 1.0),
    )
        @test_throws InvalidConfiguration OpticalPlaneMetadata(
            FocalPlane(), metadata_storage;
            coordinate_domain=AngularCoordinates(),
            sampling=invalid_sampling)
    end
    for invalid_origin in ((Inf, 0.0), (NaN, 0.0))
        @test_throws InvalidConfiguration OpticalPlaneMetadata(
            FocalPlane(), metadata_storage;
            coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
            origin=invalid_origin)
    end

    host_view = @view metadata_storage[:, :]
    host_wrappers = (
        host_view,
        reshape(host_view, 1, 4),
        transpose(host_view),
        PermutedDimsArray(host_view, (2, 1)),
        reinterpret(Float32, host_view),
    )
    host_device = plane_device(metadata_storage)
    for wrapper in host_wrappers
        @test plane_device(wrapper) == host_device
        @test backend(wrapper) isa CPUBackend
        @test AdaptiveOpticsSim.array_backend_selector(typeof(wrapper)) isa
            CPUBackend
        wrapper_metadata = OpticalPlaneMetadata(FocalPlane(), wrapper;
            coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition())
        wrapper_map = IntensityMap(wrapper_metadata, wrapper)
        @test wrapper_map.values === wrapper
    end

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0,
        pupil_reflectivity=0.25)
    wavefront = PupilFunction(tel)
    @test wavefront.metadata.spectral isa AchromaticSpectralCoordinate

    physical_source = Source(band=:custom, wavelength=1.0e-6,
        photon_irradiance=3.0)
    physical_field = ElectricField(wavefront, physical_source;
        zero_padding=2)
    physical_formation = prepare_pupil_field(tel, wavefront,
        physical_source, physical_field)
    fill_electric_field!(physical_field, wavefront, physical_formation)
    expected_rate = photon_irradiance(physical_source) *
        prod(wavefront.metadata.sampling) * sum(abs2, wavefront.amplitude)
    @test source_radiometry(physical_source) isa
        PhysicalPhotonIrradianceSource
    @test physical_field.metadata.normalization isa PhotonRateNormalization
    @test physical_field.metadata.spatial_measure isa CellIntegratedMeasure
    @test physical_field.metadata.coherence isa CoherentFieldCombination
    @test sum(abs2, physical_field.values) ≈ expected_rate

    physical_propagation = FraunhoferPropagation(physical_field)
    physical_intensity = IntensityMap(physical_field, physical_propagation)
    fraunhofer_intensity_from_field!(physical_intensity, physical_field,
        physical_propagation)
    @test physical_intensity.metadata.normalization isa
        PhotonRateNormalization
    @test physical_intensity.metadata.spatial_measure isa
        CellIntegratedMeasure
    @test physical_intensity.metadata.coherence isa
        IncoherentIntensityAddition
    @test sum(physical_intensity.values) ≈ expected_rate

    normalized_source = Source(band=:custom, wavelength=1.0e-6,
        normalized_power=2.5)
    normalized_field = ElectricField(wavefront, normalized_source;
        zero_padding=2)
    normalized_formation = prepare_pupil_field(tel, wavefront,
        normalized_source, normalized_field)
    fill_electric_field!(normalized_field, wavefront, normalized_formation)
    @test source_radiometry(normalized_source) isa NormalizedTestSource
    @test_throws InvalidConfiguration photon_irradiance(normalized_source)
    @test normalized_field.metadata.normalization isa
        DimensionlessNormalization
    @test normalized_field.metadata.spatial_measure isa CellIntegratedMeasure
    @test sum(abs2, normalized_field.values) ≈ 2.5

    source_bundle = SpectralBundle(
        [wavelength(physical_source), 1.1 * wavelength(physical_source)],
        [0.5, 0.5])
    spectral_leaf = with_spectrum(physical_source, source_bundle)
    directional = Asterism([physical_source])
    extended = with_extended_source(physical_source,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))
    third_party_expansion = TestExpandedSourceWrapper(physical_source)
    @test_throws UnsupportedAlgorithm with_spectrum(spectral_leaf,
        source_bundle)
    @test_throws UnsupportedAlgorithm with_spectrum(directional,
        source_bundle)
    @test_throws UnsupportedAlgorithm with_spectrum(extended,
        source_bundle)
    @test_throws UnsupportedAlgorithm SpectralSource(physical_source,
        source_bundle)
    @test_throws UnsupportedAlgorithm SpectralSource(spectral_leaf,
        source_bundle)
    @test_throws UnsupportedAlgorithm Asterism([third_party_expansion])
    @test_throws UnsupportedAlgorithm Asterism(
        AdaptiveOpticsSim.AbstractSource[
            physical_source,
            third_party_expansion,
        ])

    normalized_propagation = FraunhoferPropagation(normalized_field)
    normalized_intensity = IntensityMap(normalized_field,
        normalized_propagation)
    fraunhofer_intensity_from_field!(normalized_intensity, normalized_field,
        normalized_propagation)
    @test normalized_intensity.metadata.normalization isa
        DimensionlessNormalization
    @test sum(normalized_intensity.values) ≈ 2.5

    normalized_direct = IntensityMap(normalized_field,
        normalized_propagation)
    @test_throws InvalidConfiguration prepare_direct_imaging(tel, wavefront,
        normalized_source, normalized_field, normalized_direct)

    @test_throws InvalidConfiguration pupil_photon_rate_map(tel,
        normalized_source)
    @test_throws InvalidConfiguration prepare_direct_imaging(tel, wavefront,
        normalized_source; zero_padding=2)

    float32_overflow = 2 * Float64(floatmax(Float32))
    for invalid_wavelength in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=invalid_wavelength, normalized_power=1.0,
            T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            wavelength=invalid_wavelength, normalized_power=1.0,
            T=Float32)
    end
    for invalid_value in (-1.0, Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, photon_irradiance=invalid_value,
            T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, normalized_power=invalid_value,
            T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            photon_irradiance=invalid_value, T=Float32)
    end
    @test source_radiometric_value(Source(band=:custom,
        wavelength=1.0e-6, normalized_power=0.0, T=Float32)) === 0.0f0

    for invalid_finite_value in (Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, magnitude=invalid_finite_value,
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6,
            coordinates=(invalid_finite_value, 0.0),
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6,
            coordinates=(1.0, invalid_finite_value),
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            magnitude=invalid_finite_value, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            coordinates=(invalid_finite_value, 0.0), T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            laser_coordinates=(invalid_finite_value, 0.0), T=Float32)
    end
    for invalid_altitude in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration LGSSource(
            altitude=invalid_altitude, T=Float32)
    end
    for invalid_nonnegative in (-1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration LGSSource(
            elongation_factor=invalid_nonnegative, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            fwhm_spot_up=invalid_nonnegative, T=Float32)
    end
    @test LGSSource(elongation_factor=0.0,
        fwhm_spot_up=0.0).params.elongation_factor == 0.0

    valid_profile = [80_000.0 90_000.0 100_000.0; 0.25 0.5 0.25]
    profile_source = LGSSource(na_profile=valid_profile, T=Float32)
    @test profile_source.params.altitude ≈ 90_000.0f0
    @test profile_source.params.na_profile == Float32.(valid_profile)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=zeros(2, 0), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 90_000.0; 0.0 0.0], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[0.0 90_000.0; 0.5 0.5], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 Inf; 0.5 0.5], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 float32_overflow; 0.5 0.5], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 90_000.0; -0.1 1.1], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 90_000.0; NaN 1.0], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 90_000.0; Inf 1.0], T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        na_profile=[80_000.0 90_000.0;
            floatmax(Float32) floatmax(Float32)], T=Float32)

    @test_throws InvalidConfiguration GeometricAtmosphericPropagation(
        chromatic_reference_wavelength=Inf, T=Float32)
    @test_throws InvalidConfiguration LayeredFresnelAtmosphericPropagation(
        chromatic_reference_wavelength=-1.0, T=Float32)
    for invalid_wavelength in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration SpectralBundle(
            [invalid_wavelength], [1.0]; T=Float32)
    end
    for invalid_weight in (-1.0, Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration SpectralBundle(
            [1.0e-6], [invalid_weight]; T=Float32)
    end
    @test_throws InvalidConfiguration SpectralBundle(
        SpectralSample{Float32}[
            SpectralSample{Float32}(1.0f-6, floatmax(Float32)),
            SpectralSample{Float32}(1.1f-6, floatmax(Float32)),
        ])
    parameterized_samples = SpectralSample{Float64}[
        SpectralSample(1.0e-6, 0.25),
        SpectralSample(1.1e-6, 0.75),
    ]
    parameterized_bundle = SpectralBundle{Float64,
        typeof(parameterized_samples)}(parameterized_samples)
    @test parameterized_bundle.samples == parameterized_samples
    @test_throws InvalidConfiguration SpectralBundle{Float64,
        Vector{SpectralSample{Float64}}}([
            SpectralSample(NaN, -1.0),
        ])
    @test_throws InvalidConfiguration SpectralBundle{Float64,
        Vector{SpectralSample{Float64}}}([
            SpectralSample(1.0e-6, 2.0),
        ])

    density_field = ElectricField(wavefront, physical_source;
        normalization=DimensionlessNormalization(),
        spatial_measure=SpatialDensityMeasure(),
        coherence=CoherentFieldCombination())
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        physical_source, density_field; amplitude_scale=1.0)

    noncombinable_field = ElectricField(wavefront, physical_source;
        normalization=DimensionlessNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=NonCombinableProduct())
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        physical_source, noncombinable_field; amplitude_scale=1.0)
    noncombinable_propagation = FraunhoferPropagation(noncombinable_field)
    @test_throws InvalidConfiguration IntensityMap(noncombinable_field,
        noncombinable_propagation)

    unspecified_values = zeros(ComplexF64, 8, 8)
    unspecified_metadata = OpticalPlaneMetadata(PupilPlane(),
        unspecified_values; coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=CoherentFieldCombination())
    unspecified_field = ElectricField(unspecified_metadata,
        unspecified_values)
    @test_throws InvalidConfiguration FraunhoferPropagation(
        unspecified_field)
    @test_throws InvalidConfiguration FresnelPropagation(unspecified_field;
        distance_m=1.0)

    function intensity_map(values;
        sampling=(1.0, 1.0), wavelength_m=1.0e-6,
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
        metadata = OpticalPlaneMetadata(FocalPlane(), values;
            coordinate_domain=AngularCoordinates(), sampling=sampling,
            spectral=MonochromaticChannel(wavelength_m),
            normalization=normalization, spatial_measure=spatial_measure,
            coherence=coherence)
        return IntensityMap(metadata, values)
    end

    output = intensity_map(zeros(4, 4))
    first_input = intensity_map(fill(1.0, 4, 4))
    second_input = intensity_map(fill(2.0, 4, 4))
    sum_plan = prepare_incoherent_sum(output, first_input, second_input)
    @test !applicable(PreparedIncoherentSum, output.metadata,
        (first_input.metadata, second_input.metadata))
    @test !applicable(
        PreparedIncoherentSum{typeof(output.metadata),
            typeof((first_input.metadata, second_input.metadata))},
        output.metadata, (first_input.metadata, second_input.metadata))
    @test @inferred(accumulate_intensity!(output,
        (first_input, second_input), sum_plan)) === output
    @test output.values == fill(3.0, 4, 4)
    @test @allocated(accumulate_intensity!(output,
        (first_input, second_input), sum_plan)) == 0

    integrated_metadata = OpticalPlaneMetadata(FocalPlane(), zeros(4, 4);
        coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
        spectral=IntegratedSpectralChannel(:science_passband),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    integrated_output = IntensityMap(integrated_metadata, zeros(4, 4))
    integrated_first = IntensityMap(integrated_metadata, fill(1.0, 4, 4))
    integrated_second = IntensityMap(integrated_metadata, fill(2.0, 4, 4))
    integrated_plan = prepare_incoherent_sum(integrated_output,
        integrated_first, integrated_second)
    accumulate_intensity!(integrated_output,
        (integrated_first, integrated_second), integrated_plan)
    @test integrated_output.values == fill(3.0, 4, 4)

    achromatic_metadata = OpticalPlaneMetadata(FocalPlane(), zeros(4, 4);
        coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
        spectral=AchromaticSpectralCoordinate(),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    achromatic_output = IntensityMap(achromatic_metadata, zeros(4, 4))
    achromatic_first = IntensityMap(achromatic_metadata, fill(1.0, 4, 4))
    achromatic_second = IntensityMap(achromatic_metadata, fill(2.0, 4, 4))
    achromatic_plan = prepare_incoherent_sum(achromatic_output,
        achromatic_first, achromatic_second)
    accumulate_intensity!(achromatic_output,
        (achromatic_first, achromatic_second), achromatic_plan)
    @test achromatic_output.values == fill(3.0, 4, 4)
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        achromatic_output, first_input)

    other_integrated_metadata = OpticalPlaneMetadata(FocalPlane(),
        fill(1.0, 4, 4); coordinate_domain=AngularCoordinates(),
        sampling=(1.0, 1.0),
        spectral=IntegratedSpectralChannel(:other_passband),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    other_integrated = IntensityMap(other_integrated_metadata,
        fill(1.0, 4, 4))
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        integrated_output, other_integrated)

    unspecified_sum_values = zeros(4, 4)
    unspecified_sum_metadata = OpticalPlaneMetadata(FocalPlane(),
        unspecified_sum_values; coordinate_domain=AngularCoordinates(),
        sampling=(1.0, 1.0), normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    unspecified_sum = IntensityMap(unspecified_sum_metadata,
        unspecified_sum_values)
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        unspecified_sum, unspecified_sum)
    identity_alias = IntensityMap{
        typeof(first_input.metadata),typeof(output.values),CPUBackend,
    }(first_input.metadata, output.values)
    view_values = @view output.values[:, :]
    view_alias = IntensityMap{
        typeof(first_input.metadata),typeof(view_values),CPUBackend,
    }(first_input.metadata, view_values)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output,
        identity_alias)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output,
        view_alias)
    fill!(output.values, 7.0)
    @test_throws InvalidConfiguration accumulate_intensity!(output,
        (view_alias, second_input), sum_plan)
    @test output.values == fill(7.0, 4, 4)
    unplanned_input = intensity_map(fill(1.0, 4, 4);
        sampling=(2.0, 2.0))
    @test_throws InvalidConfiguration accumulate_intensity!(output,
        (unplanned_input, second_input), sum_plan)
    @test output.values == fill(7.0, 4, 4)
    @test !applicable(accumulate_intensity!, output.values, physical_field)
    @test_throws DimensionMismatchError accumulate_intensity!(output,
        (first_input,), sum_plan)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output)

    incompatible_normalization = intensity_map(fill(1.0, 4, 4);
        normalization=DimensionlessNormalization())
    incompatible_measure = intensity_map(fill(1.0, 4, 4);
        spatial_measure=SpatialDensityMeasure())
    incompatible_policy = intensity_map(fill(1.0, 4, 4);
        coherence=NonCombinableProduct())
    incompatible_sampling = intensity_map(fill(1.0, 4, 4);
        sampling=(2.0, 2.0))
    incompatible_spectral = intensity_map(fill(1.0, 4, 4);
        wavelength_m=1.1e-6)
    for incompatible in (
        incompatible_normalization,
        incompatible_measure,
        incompatible_policy,
        incompatible_sampling,
        incompatible_spectral,
    )
        @test_throws InvalidConfiguration prepare_incoherent_sum(output,
            incompatible)
    end

    bundle = OpticalProductBundle(first_input, incompatible_spectral)
    @test length(bundle) == 2
    @test bundle[1] === first_input
    @test bundle[2] === incompatible_spectral

    mixed_wavelengths = Asterism([
        Source(band=:custom, wavelength=1.0e-6, photon_irradiance=1.0),
        Source(band=:custom, wavelength=1.1e-6, photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration prepare_direct_imaging(tel,
        PupilFunction(tel), mixed_wavelengths)
end

@testset "Telescope and direct imaging" begin
    for invalid_reflectivity in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration Telescope(resolution=8,
            diameter=8.0, pupil_reflectivity=invalid_reflectivity)
    end
    invalid_reflectivity_map = ones(8, 8)
    invalid_reflectivity_map[1, 1] = NaN
    @test_throws InvalidConfiguration Telescope(resolution=8,
        diameter=8.0, pupil_reflectivity=invalid_reflectivity_map)

    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilFunction(tel)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=2)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test size(field.values) == (64, 64)
    @test field.metadata.coordinate_domain isa MetricCoordinates
    @test field.metadata.spectral == MonochromaticChannel(wavelength(src))
    fraunhofer = FraunhoferPropagation(field)
    @test fraunhofer.output_metadata.coordinate_domain isa AngularCoordinates
    centered_rate = similar(field.values, Float64)
    fraunhofer_intensity_from_field!(centered_rate, field, fraunhofer)
    direct = prepare_direct_imaging(tel, wavefront, src;
        zero_padding=2)
    rate_map = form_direct_image!(direct)
    rate = intensity_values(rate_map)
    @test size(rate) == (64, 64)
    @test maximum(rate) > 0
    @test isfinite(sum(rate))
    @test centered_rate ≈ rate
    @test !hasproperty(tel.state, :psf)
    @test !hasproperty(tel.state, :psf_stack)
    @test !hasproperty(tel.state, :psf_workspace)

    tel_dim = Telescope(resolution=32, diameter=8.0, central_obstruction=0.2,
        pupil_reflectivity=0.25)
    dim_pupil = PupilFunction(tel_dim)
    dim_direct = prepare_direct_imaging(tel_dim, dim_pupil, src;
        zero_padding=2)
    rate_dim = intensity_values(form_direct_image!(dim_direct))
    @test sum(rate_dim) ≈ 0.25 * sum(rate)
    photon_rate = pupil_photon_rate_map(tel_dim, src)
    @test size(photon_rate) == size(pupil_mask(tel_dim))
    @test maximum(photon_rate) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"
    reflectivity_before_invalid_set = copy(pupil_reflectivity(tel_dim))
    @test_throws InvalidConfiguration set_pupil_reflectivity!(tel_dim, -0.1)
    @test pupil_reflectivity(tel_dim) == reflectivity_before_invalid_set

    tel_simple = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src_simple = Source(band=:I, magnitude=0.0)
    wavefront_simple = PupilFunction(tel_simple)
    field_simple = ElectricField(wavefront_simple, src_simple;
        zero_padding=1)
    formation_simple = prepare_pupil_field(tel_simple, wavefront_simple,
        src_simple, field_simple)
    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
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
    propagated_intensity = similar(propagated.values, Float64)
    @. propagated_intensity = abs2(propagated.values)
    @test propagated_intensity ≈ centered_rate atol=1e-10 rtol=1e-10
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

    atm_tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    atm_src = Source(band=:I, magnitude=0.0)
    atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance_by!(atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(1))
    geom_prop = AtmosphericFieldPropagation(atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldSynchronousPlan
    geom_field = propagate_atmosphere_field!(geom_prop, atm, atm_tel, atm_src)
    tel_geom = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    propagate!(atm, tel_geom, atm_src)
    collapsed_wavefront = PupilFunction(tel_geom)
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
    advance_by!(fresnel_atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(2))
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

    normalized_atm_src = Source(band=:custom,
        wavelength=wavelength(atm_src), normalized_power=2.0)
    normalized_spectral = with_spectrum(normalized_atm_src,
        SpectralBundle([wavelength(atm_src), 1.1 * wavelength(atm_src)],
            [0.75, 0.25]))
    normalized_spectral_prop = AtmosphericFieldPropagation(atm, atm_tel,
        normalized_spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test all(slice.field.metadata.normalization isa
        DimensionlessNormalization for slice in
        normalized_spectral_prop.state.slices)
    normalized_spectral_intensity = atmospheric_intensity!(
        normalized_spectral_prop, atm, atm_tel, normalized_spectral)
    @test sum(normalized_spectral_intensity) ≈ 2.0 rtol=1e-10
end

@testset "Explicit optical products and surfaces" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    telescope_opd_before = copy(opd_map(tel))
    aperture_before = copy(pupil_mask(tel))
    reflectivity_before = copy(pupil_reflectivity(tel))

    path_a = PupilFunction(tel)
    path_b = PupilFunction(tel)
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
    field = ElectricField(path_a, src; zero_padding=2)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    fill!(output.values, 0)
    prepared = prepare_direct_imaging(tel, path_a, src, field, output)
    explicit_direct_image_cycle!(output, field, path_a, prepared.plan,
        prepared.workspace)
    @test @allocated(explicit_direct_image_cycle!(output, field, path_a,
        prepared.plan, prepared.workspace)) == 0
    @test path_b.opd == path_b_before_propagation
    @test !hasproperty(tel.state, :psf)

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2)
    spatial_field = ElectricField(path_a, src; zero_padding=2,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    spatial_formation = prepare_pupil_field(tel, path_a, src, spatial_field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(spatial_field, path_a, spatial_formation)
    spatial_output = PupilFunction(tel)
    spatial_plan = prepare_spatial_filter(tel, spatial_filter, spatial_field,
        spatial_output)
    spatial_workspace = SpatialFilterWorkspace(spatial_filter)
    explicit_spatial_filter_cycle!(spatial_output, spatial_field, spatial_filter,
        spatial_plan, spatial_workspace)
    @test @allocated(explicit_spatial_filter_cycle!(spatial_output, spatial_field,
        spatial_filter, spatial_plan, spatial_workspace)) == 0
end

@testset "Optical-plane compatibility validation" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilFunction(tel)
    field = ElectricField(wavefront, src; zero_padding=1)
    values = similar(field.values)

    sampling_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=(2.0, 2.0), spectral=field.metadata.spectral)
    sampling_field = ElectricField(sampling_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, sampling_field)

    coordinate_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    coordinate_field = ElectricField(coordinate_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, coordinate_field)
    @test_throws InvalidConfiguration FraunhoferPropagation(coordinate_field)
    spatial_filter = SpatialFilter(tel; zero_padding=1)
    @test_throws InvalidConfiguration prepare_spatial_filter(tel,
        spatial_filter, coordinate_field, wavefront)

    kind_metadata = OpticalPlaneMetadata(FocalPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    kind_field = ElectricField(kind_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, kind_field)

    origin_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=(0.0, 0.0),
        spectral=field.metadata.spectral)
    origin_field = ElectricField(origin_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, origin_field)

    centering_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        centering=(SampleCentered, SampleCentered),
        spectral=field.metadata.spectral)
    centering_field = ElectricField(centering_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, centering_field)

    orientation_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        orientation=PlaneAxisOrientation((:y, :x)),
        spectral=field.metadata.spectral)
    orientation_field = ElectricField(orientation_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, orientation_field)

    spectral_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=MonochromaticChannel(1.1 * wavelength(src)))
    spectral_field = ElectricField(spectral_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, spectral_field)

    float32_values = Matrix{ComplexF32}(undef, size(values))
    numeric_metadata = OpticalPlaneMetadata(PupilPlane(), float32_values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    numeric_field = ElectricField(numeric_metadata, float32_values)
    @test_throws InvalidConfiguration prepare_pupil_field(tel, wavefront,
        src, numeric_field)

    wrong_size_values = Matrix{ComplexF64}(undef, 9, 9)
    dimension_metadata = OpticalPlaneMetadata(PupilPlane(),
        wrong_size_values; coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling,
        spectral=field.metadata.spectral)
    dimension_field = ElectricField(dimension_metadata, wrong_size_values)
    @test_throws DimensionMismatchError prepare_pupil_field(tel, wavefront,
        src, dimension_field)

    propagation = FraunhoferPropagation(field)
    wrong_destination_values = Matrix{Float64}(undef, 9, 9)
    wrong_destination_metadata = OpticalPlaneMetadata(FocalPlane(),
        wrong_destination_values;
        coordinate_domain=propagation.output_metadata.coordinate_domain,
        sampling=propagation.output_metadata.sampling,
        spectral=propagation.output_metadata.spectral)
    wrong_destination = IntensityMap(wrong_destination_metadata,
        wrong_destination_values)
    @test_throws DimensionMismatchError prepare_direct_imaging(tel, wavefront,
        src, field, wrong_destination)

    declared_device_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
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

    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.25)
    expected = falses(32, 32)
    build_mask!(expected, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test pupil_mask(tel) == expected
    apply_spiders!(tel; thickness=0.4, angles=[0.0, 90.0])
    manual = copy(expected)
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=0.0))
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=pi / 2))
    @test pupil_mask(tel) == manual

    sf_tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
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
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end
