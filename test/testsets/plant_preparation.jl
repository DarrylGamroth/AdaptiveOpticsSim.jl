struct PlantPreparationTestAtmosphere <: AdaptiveOpticsSim.AbstractAtmosphere end

struct DirectSciencePathModel{R}
    zero_padding::Int
    revision::R
end

struct ShackHartmannPlantPathModel{R}
    n_lenslets::Int
    n_pix_subap::Int
    revision::R
end

abstract type PlantPathRequirement end
struct MatchingPlantPath <: PlantPathRequirement end
struct MismatchedSourceGeometry <: PlantPathRequirement end
struct MismatchedRadiometry <: PlantPathRequirement end
struct MismatchedRevision <: PlantPathRequirement end
struct MismatchedBackend <: PlantPathRequirement end
struct MismatchedDevice <: PlantPathRequirement end

struct FramePlantAcquisitionModel{T<:AbstractFloat,R<:PlantPathRequirement}
    exposure::T
    requirement::R
end

struct ContractWFSPlantAcquisitionModel{T<:AbstractFloat}
    exposure::T
end

struct UnsupportedPreparedPathModel end
struct UnsupportedPreparedAcquisitionModel end
struct InvalidPreparedPathModel end
struct InvalidPreparedAcquisitionModel end

struct PlantBindingOnlyExecution{I,R}
    input::I
    result::R
end

struct UnsupportedPlantSource <: AbstractSource end

function AdaptiveOpticsSim.validate_path_execution_binding(
    execution::PlantBindingOnlyExecution, input, result)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "test path execution does not match its prepared products"))
    return nothing
end

function AdaptiveOpticsSim.execute_path!(result, input,
    execution::PlantBindingOnlyExecution)
    AdaptiveOpticsSim.validate_path_execution_binding(execution, input,
        result)
    return result
end

for model in (
    DirectSciencePathModel,
    ShackHartmannPlantPathModel,
    FramePlantAcquisitionModel,
    ContractWFSPlantAcquisitionModel,
    UnsupportedPreparedPathModel,
    UnsupportedPreparedAcquisitionModel,
    InvalidPreparedPathModel,
    InvalidPreparedAcquisitionModel,
)
    @eval AdaptiveOpticsSim.plant_model_definition_style(
        ::Type{<:$model},
    ) = ColdPlantModelDefinition()
end

function AdaptiveOpticsSim.prepare_path_executor(
    ::InvalidPreparedPathModel,
    ::OpticalPathDefinition,
    ::AbstractSource,
    ::Telescope,
    ::PlantPreparationTestAtmosphere,
)
    return nothing
end

function AdaptiveOpticsSim.prepare_acquisition_owner(
    ::InvalidPreparedAcquisitionModel,
    ::AcquisitionDefinition,
    ::PreparedPathExecutor,
)
    return nothing
end

function AdaptiveOpticsSim.prepare_path_executor(
    model::DirectSciencePathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    ::PlantPreparationTestAtmosphere,
)
    pupil = PupilFunction(telescope)
    imaging = prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        optical_model=(kind=:direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=model.revision,
    )
end

function AdaptiveOpticsSim.prepare_path_executor(
    model::ShackHartmannPlantPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    ::PlantPreparationTestAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T)
    sensor = ShackHartmannWFS(telescope;
        n_lenslets=model.n_lenslets,
        n_pix_subap=model.n_pix_subap,
        mode=Diffractive(),
        T=T,
        backend=backend(telescope))
    front_end = ShackHartmannOpticalFrontEnd(sensor.front_end, source)
    output = shack_hartmann_rate_map(front_end, pupil)
    plan = prepare_wfs_optical_formation(front_end, pupil, output)
    execution = WFSOpticalPathExecution(plan)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        pupil,
        output,
        execution;
        optical_model=(kind=:shack_hartmann,
            n_lenslets=model.n_lenslets,
            n_pix_subap=model.n_pix_subap),
        propagation_model=:microlens_fraunhofer,
        model_revisions=(definition=model.revision,
            layout=AdaptiveOpticsSim.subaperture_layout_revision(
                front_end.layout)),
    )
end

@inline _require_test_path(::MatchingPlantPath,
    path::PreparedPathExecutor) = require_path_result(path)

@inline _require_test_path(::MismatchedSourceGeometry,
    path::PreparedPathExecutor) = require_path_result(path;
        source_geometry=(kind=:different_direction,))

@inline _require_test_path(::MismatchedRadiometry,
    path::PreparedPathExecutor) = require_path_result(path;
        radiometry=(policy=:different, value=-1))

@inline _require_test_path(::MismatchedRevision,
    path::PreparedPathExecutor) = require_path_result(path;
        revisions=(telescope=path.key.revisions.telescope,
            model=:stale))

@inline _require_test_path(::MismatchedBackend,
    path::PreparedPathExecutor) = require_path_result(path;
        backend=CUDABackend())

@inline _require_test_path(::MismatchedDevice,
    path::PreparedPathExecutor) = require_path_result(path;
        device=ContractPlaneDevice(404))

function AdaptiveOpticsSim.prepare_acquisition_owner(
    model::FramePlantAcquisitionModel,
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    _require_test_path(model.requirement, path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), response_model=NullFrameResponse(),
        T=T, backend=path.key.backend)
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(execution.observation)
    return PreparedAcquisitionOwner(definition, path, execution, products)
end

function AdaptiveOpticsSim.prepare_acquisition_owner(
    model::ContractWFSPlantAcquisitionModel,
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), response_model=NullFrameResponse(),
        T=T, backend=path.key.backend)
    observation = WFSObservation(similar(path.result.values);
        units=:electron_count, layout=:contract_rate_frame)
    acquisition = prepare_wfs_acquisition(detector, path.result,
        observation)
    measurement = WFSMeasurement(similar(path.result.values, T, 1);
        units=:electron_count, kind=:channel_totals)
    estimator = prepare_wfs_estimation(
        ContractSumEstimator(:electron_count, :channel_totals),
        observation,
        measurement,
    )
    execution = WFSAcquisitionExecution(acquisition, estimator,
        observation, measurement)
    products = AcquisitionProducts(observation, measurement)
    return PreparedAcquisitionOwner(definition, path, execution, products)
end

function captured_plant_preparation_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function plant_preparation_intensity_map(values::AbstractArray;
    kind=FocalPlane(), coordinate_domain=AngularCoordinates(),
    spectral=MonochromaticChannel(convert(eltype(values), 1.0e-6)),
    normalization=PhotonRateNormalization(),
    spatial_measure=CellIntegratedMeasure(),
    coherence=IncoherentIntensityAddition())
    sampling = (one(eltype(values)), one(eltype(values)))
    metadata = OpticalPlaneMetadata(kind, values;
        coordinate_domain=coordinate_domain,
        sampling=sampling,
        spectral=spectral,
        normalization=normalization,
        spatial_measure=spatial_measure,
        coherence=coherence)
    return IntensityMap(metadata, values)
end

function assert_plant_preparation_error(f, component::Symbol,
    reason::Symbol)
    error = captured_plant_preparation_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function prepared_path_execution_allocations(path::PreparedPathExecutor)
    execute_path!(path)
    return @allocated execute_path!(path)
end

function prepared_acquisition_execution_allocations(
    owner::PreparedAcquisitionOwner, rng)
    execute_acquisition!(owner, rng)
    return @allocated execute_acquisition!(owner, rng)
end

@testset "Prepared plant paths and acquisition owners" begin
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = PlantPreparationTestAtmosphere()
    science_source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(3), T=T)
    wfs_source = Source(band=:custom, wavelength=T(0.65e-6),
        photon_irradiance=T(5), T=T)

    science_definition = OpticalPathDefinition(:science, science_source,
        DirectSciencePathModel(2, UInt(3)))
    wfs_definition = OpticalPathDefinition(:wfs, wfs_source,
        ShackHartmannPlantPathModel(2, 4, UInt(7)))
    fast_definition = AcquisitionDefinition(:fast_science, :science,
        FramePlantAcquisitionModel(T(0.25), MatchingPlantPath()))
    slow_definition = AcquisitionDefinition(:slow_science, :science,
        FramePlantAcquisitionModel(T(0.75), MatchingPlantPath()))
    wfs_acquisition_definition = AcquisitionDefinition(:wfs_frame, :wfs,
        ContractWFSPlantAcquisitionModel(T(0.5)))
    definition = PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        paths=(science=science_definition, wfs=wfs_definition),
        acquisitions=(fast_science=fast_definition,
            slow_science=slow_definition,
            wfs_frame=wfs_acquisition_definition),
    )

    plant = prepare_plant(definition)
    @test plant isa PreparedPlant
    @test prepared_paths(plant) isa Tuple
    @test prepared_acquisitions(plant) isa Tuple
    @test length(prepared_paths(plant)) == 2
    @test length(prepared_acquisitions(plant)) == 3
    @test isconcretetype(typeof(prepared_paths(plant)))
    @test isconcretetype(typeof(prepared_acquisitions(plant)))

    science_path = prepared_path(plant, :science)
    wfs_path = prepared_path(plant, OpticalPathID(:wfs))
    fast = prepared_acquisition(plant, :fast_science)
    slow = prepared_acquisition(plant, AcquisitionID(:slow_science))
    wfs = prepared_acquisition(plant, :wfs_frame)

    @test path_input(science_path) === science_path.input
    @test path_result(science_path) === science_path.result
    @test path_result_key(science_path) === science_path.key
    @test science_path.key.sampling_contract isa
        AdaptiveOpticsSim.InstantaneousOpticalSample
    @test science_path.key.backend isa CPUBackend
    @test science_path.key.device == AdaptiveOpticsSim.HostPlaneDevice()
    equivalent_key = AdaptiveOpticsSim.PathResultKey(
        science_path.key.source_geometry,
        science_path.key.spectral_sampling,
        science_path.key.radiometry,
        science_path.key.optical_model,
        science_path.key.sampling_contract,
        science_path.key.propagation_model,
        science_path.key.output_plane,
        science_path.key.revisions,
        science_path.key.backend,
        science_path.key.device,
    )
    @test equivalent_key == science_path.key
    @test isequal(equivalent_key, science_path.key)
    @test hash(equivalent_key) == hash(science_path.key)

    sodium_profile = T[80_000 90_000 100_000; 0.2 0.6 0.2]
    lgs_source = LGSSource(wavelength=T(589e-9),
        photon_irradiance=T(2), na_profile=sodium_profile,
        laser_coordinates=(T(1), T(-0.5)), elongation_factor=T(1.2),
        fwhm_spot_up=T(0.8), T=T)
    lgs_geometry = AdaptiveOpticsSim.path_source_geometry_key(lgs_source)
    @test lgs_geometry.kind === LGSSource
    @test lgs_geometry.sodium_profile == sodium_profile
    @test lgs_geometry.sodium_profile !== lgs_source.params.na_profile
    @test only(AdaptiveOpticsSim.path_source_spectral_key(
        lgs_source)).wavelength_m == wavelength(lgs_source)
    @test AdaptiveOpticsSim.path_source_radiometry_key(
        lgs_source).value == source_radiometric_value(lgs_source)
    lgs_key = AdaptiveOpticsSim.PathResultKey(
        lgs_geometry,
        AdaptiveOpticsSim.path_source_spectral_key(lgs_source),
        AdaptiveOpticsSim.path_source_radiometry_key(lgs_source),
        science_path.key.optical_model,
        science_path.key.sampling_contract,
        science_path.key.propagation_model,
        science_path.key.output_plane,
        science_path.key.revisions,
        science_path.key.backend,
        science_path.key.device,
    )
    lgs_hash = hash(lgs_key)
    keyed_results = Dict(lgs_key => :prepared)
    exposed_profile = lgs_key.source_geometry.sodium_profile
    exposed_profile[1, 1] = -one(T)
    lgs_geometry.sodium_profile[1, 1] = -T(2)
    @test hash(lgs_key) == lgs_hash
    @test keyed_results[lgs_key] === :prepared
    @test lgs_key.source_geometry.sodium_profile == sodium_profile

    spectral_source = with_spectrum(science_source, SpectralBundle(
        T[0.75e-6, 0.85e-6], T[0.4, 0.6]))
    companion_source = Source(band=:custom, wavelength=wavelength(science_source),
        photon_irradiance=T(2), coordinates=(T(0.2), T(30)), T=T)
    asterism_source = Asterism([science_source, companion_source])
    extended_source = with_extended_source(science_source,
        PointCloudSourceModel([(T(0), T(0)), (T(0.1), T(-0.05))],
            T[0.25, 0.75]))

    @test AdaptiveOpticsSim.path_source_geometry_key(spectral_source) ==
        AdaptiveOpticsSim.path_source_geometry_key(science_source)
    @test length(AdaptiveOpticsSim.path_source_spectral_key(
        spectral_source)) == 2
    @test AdaptiveOpticsSim.path_source_radiometry_key(spectral_source) ==
        AdaptiveOpticsSim.path_source_radiometry_key(science_source)
    @test length(AdaptiveOpticsSim.path_source_geometry_key(
        asterism_source)) == 2
    @test length(AdaptiveOpticsSim.path_source_spectral_key(
        asterism_source)) == 2
    @test length(AdaptiveOpticsSim.path_source_radiometry_key(
        asterism_source)) == 2
    @test length(AdaptiveOpticsSim.path_source_geometry_key(
        extended_source)) == 2
    @test AdaptiveOpticsSim.path_source_spectral_key(extended_source) ==
        AdaptiveOpticsSim.path_source_spectral_key(science_source)
    @test length(AdaptiveOpticsSim.path_source_radiometry_key(
        extended_source)) == 2

    unsupported_source = UnsupportedPlantSource()
    for (key_function, reason) in (
        (AdaptiveOpticsSim.path_source_geometry_key, :source_geometry),
        (AdaptiveOpticsSim.path_source_spectral_key, :spectral_sampling),
        (AdaptiveOpticsSim.path_source_radiometry_key, :radiometry),
    )
        assert_plant_preparation_error(
            () -> key_function(unsupported_source), :path, reason)
    end

    for (id, composite_source, revision) in (
        (:spectral_science, spectral_source, UInt(11)),
        (:asterism_science, asterism_source, UInt(12)),
    )
        composite_definition = OpticalPathDefinition(id, composite_source,
            DirectSciencePathModel(2, revision))
        composite_path = prepare_path_executor(composite_definition,
            telescope, atmosphere)
        @test execute_path!(composite_path) === composite_path.result
        assert_plant_preparation_error(
            () -> execute_path!(composite_path.result,
                PupilFunction(telescope), composite_path.execution),
            :path,
            :prepared_binding,
        )
    end

    foreign_values = ContractDeviceArray(
        zeros(Complex{T}, size(science_path.input.opd)),
        ContractPlaneDevice(17),
    )
    foreign_metadata = OpticalPlaneMetadata(PupilPlane(), foreign_values;
        coordinate_domain=MetricCoordinates(),
        sampling=science_path.input.metadata.sampling,
        origin=science_path.input.metadata.origin,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    foreign_field = ElectricField(foreign_metadata, foreign_values)
    assert_plant_preparation_error(
        () -> PreparedPathExecutor(
            science_definition,
            science_path.source,
            telescope,
            (science_path.input, foreign_field),
            science_path.result,
            science_path.execution;
            optical_model=science_path.key.optical_model,
            sampling_contract=science_path.key.sampling_contract,
            propagation_model=science_path.key.propagation_model,
            model_revisions=science_path.key.revisions.model,
        ),
        :path,
        :device,
    )

    prepared_test_path = (input, result;
            execution=PlantBindingOnlyExecution(input, result)) ->
        PreparedPathExecutor(
            science_definition,
            science_path.source,
            telescope,
            input,
            result,
            execution;
            optical_model=science_path.key.optical_model,
            sampling_contract=science_path.key.sampling_contract,
            propagation_model=science_path.key.propagation_model,
            model_revisions=science_path.key.revisions.model,
        )

    invalid_outputs = (
        (plant_preparation_intensity_map(zeros(T, 2, 2);
            kind=PupilPlane(), coordinate_domain=MetricCoordinates()),
            :output_plane),
        (plant_preparation_intensity_map(zeros(T, 2, 2);
            normalization=DimensionlessNormalization()), :radiometry),
        (plant_preparation_intensity_map(zeros(T, 2, 2);
            spatial_measure=PointSampledMeasure()), :radiometry),
        (plant_preparation_intensity_map(zeros(T, 2, 2);
            coherence=CoherentFieldCombination()), :radiometry),
        (plant_preparation_intensity_map(zeros(T, 2, 2);
            spectral=AchromaticSpectralCoordinate()), :spectral_sampling),
    )
    for (invalid_output, reason) in invalid_outputs
        assert_plant_preparation_error(
            () -> prepared_test_path(science_path.input, invalid_output),
            :path,
            reason,
        )
    end

    density_integrated_result = plant_preparation_intensity_map(
        zeros(T, 2, 2);
        spectral=IntegratedSpectralChannel(:test_passband),
        spatial_measure=SpatialDensityMeasure())
    for reusable_result in (
        (science_path.result,),
        OpticalProductBundle(science_path.result),
        density_integrated_result,
    )
        reusable_path = prepared_test_path(science_path.input,
            reusable_result)
        @test reusable_path.result === reusable_result
    end


    source_products = AdaptiveOpticsSim.AbstractOpticalProduct[
        science_path.result]
    vector_bundle = OpticalProductBundle(source_products)
    source_products[1] = density_integrated_result
    @test only(vector_bundle.products) === science_path.result
    @test vector_bundle.products isa AbstractVector
    @test !(vector_bundle.products isa Vector)
    @test_throws MethodError push!(vector_bundle.products,
        density_integrated_result)
    @test_throws Base.CanonicalIndexError vector_bundle.products[1] =
        density_integrated_result

    assert_plant_preparation_error(
        () -> prepared_test_path(science_path.input, [science_path.result]),
        :path,
        :output_plane,
    )

    assert_plant_preparation_error(
        () -> prepared_test_path(science_path.input, IntensityMap[]),
        :path,
        :output_plane,
    )
    assert_plant_preparation_error(
        () -> prepared_test_path(science_path.input, ()),
        :path,
        :output_plane,
    )
    assert_plant_preparation_error(
        () -> prepared_test_path(science_path.input, nothing),
        :path,
        :output_plane,
    )
    assert_plant_preparation_error(
        () -> prepared_test_path((), science_path.result),
        :path,
        :input_plane,
    )
    assert_plant_preparation_error(
        () -> prepared_test_path(nothing, science_path.result),
        :path,
        :input_plane,
    )

    pupil_field = ElectricField(science_path.input, science_path.source)
    tuple_input = (science_path.input, pupil_field)
    @test prepared_test_path(tuple_input, science_path.result).input ===
        tuple_input

    unbound_result = IntensityMap(science_path.result.metadata,
        copy(science_path.result.values))
    assert_plant_preparation_error(
        () -> execute_path!(unbound_result, science_path.input,
            science_path.execution),
        :path,
        :prepared_binding,
    )
    assert_plant_preparation_error(
        () -> execute_path!(science_path.result, science_path.input, nothing),
        :path,
        :unsupported_execution,
    )

    unbound_wfs_result = IntensityMap(wfs_path.result.metadata,
        copy(wfs_path.result.values))
    @test_throws WFSPreparationError PreparedPathExecutor(
        wfs_definition,
        wfs_path.source,
        telescope,
        wfs_path.input,
        unbound_wfs_result,
        wfs_path.execution;
        optical_model=wfs_path.key.optical_model,
        sampling_contract=wfs_path.key.sampling_contract,
        propagation_model=wfs_path.key.propagation_model,
        model_revisions=wfs_path.key.revisions.model,
    )

    @test !applicable(PreparedPathExecutor, science_path.definition,
        science_path.source, science_path.telescope, science_path.input,
        science_path.result, science_path.execution, science_path.key)

    @test fast.path_result === science_path.result
    @test slow.path_result === science_path.result
    @test fast.path_key === science_path.key
    @test slow.path_key === science_path.key
    @test fast.execution.detector.state !== slow.execution.detector.state
    @test acquisition_observation(fast) !== acquisition_observation(slow)
    @test !Base.mightalias(acquisition_observation(fast),
        acquisition_observation(slow))
    @test acquisition_measurement(fast) === nothing
    @test acquisition_products(fast) === fast.products

    formed_science = @inferred execute_path!(science_path)
    @test formed_science === science_path.result
    science_before = copy(science_path.result.values)
    fast_products = @inferred execute_acquisition!(fast, Xoshiro(11))
    slow_products = @inferred execute_acquisition!(slow, Xoshiro(12))
    @test fast_products === fast.products
    @test slow_products === slow.products
    @test acquisition_observation(slow) ≈
        T(3) .* acquisition_observation(fast) atol=0 rtol=0
    @test science_path.result.values == science_before

    formed_wfs = @inferred execute_path!(wfs_path)
    @test formed_wfs === wfs_path.result
    wfs_rate_before = copy(wfs_path.result.values)
    wfs_products = @inferred execute_acquisition!(wfs, Xoshiro(13))
    @test wfs_products === wfs.products
    @test acquisition_measurement(wfs).storage[1] ≈
        T(0.5) * sum(wfs_path.result.values) atol=T(1e-12) rtol=T(1e-12)
    @test wfs_path.result.values == wfs_rate_before

    detector = fast.execution.detector
    detector_frame = output_frame(detector)
    explicit_observation = similar(detector_frame)
    fill!(explicit_observation, zero(T))
    explicit_frame_execution = FrameAcquisitionExecution(detector,
        science_path.result, explicit_observation)
    explicit_frame_products = AcquisitionProducts(explicit_observation)
    @test execute_acquisition!(explicit_frame_products, science_path.result,
        explicit_frame_execution, Xoshiro(16)) === explicit_frame_products

    raw_frame_plan = prepare_detector_acquisition(detector,
        science_path.result)
    @test !applicable(FrameAcquisitionExecution, detector, raw_frame_plan,
        detector_frame)
    @test_throws MethodError FrameAcquisitionExecution(detector,
        raw_frame_plan, detector_frame)

    invalid_observations = (
        (zeros(T, 1, 1), :shape),
        (zeros(Float32, size(detector_frame)), :numeric_type),
        (ContractDeviceArray(zeros(T, size(detector_frame)),
            ContractPlaneDevice(18)), :device),
        (detector_frame, :ownership),
    )
    for (invalid_observation, reason) in invalid_observations
        assert_plant_preparation_error(
            () -> FrameAcquisitionExecution(detector, science_path.result,
                invalid_observation),
            :acquisition,
            reason,
        )
    end

    assert_plant_preparation_error(
        () -> execute_acquisition!(
            AcquisitionProducts(similar(explicit_observation)),
            science_path.result,
            explicit_frame_execution,
            Xoshiro(17),
        ),
        :acquisition,
        :prepared_binding,
    )

    observation_only_execution = WFSAcquisitionExecution(
        wfs.execution.acquisition, wfs.execution.observation)
    observation_only_products = AcquisitionProducts(
        wfs.execution.observation)
    @test execute_acquisition!(observation_only_products, wfs.path_result,
        observation_only_execution, Xoshiro(18)) === observation_only_products

    mismatched_wfs_observation = deepcopy(wfs.execution.observation)
    mismatched_wfs_execution = WFSAcquisitionExecution(
        wfs.execution.acquisition, mismatched_wfs_observation)
    @test_throws WFSPreparationError PreparedAcquisitionOwner(
        wfs_acquisition_definition, wfs_path, mismatched_wfs_execution,
        AcquisitionProducts(mismatched_wfs_observation))

    mismatched_wfs_measurement = deepcopy(wfs.execution.measurement)
    mismatched_estimator_execution = WFSAcquisitionExecution(
        wfs.execution.acquisition, wfs.execution.estimator,
        wfs.execution.observation, mismatched_wfs_measurement)
    @test_throws WFSPreparationError PreparedAcquisitionOwner(
        wfs_acquisition_definition, wfs_path, mismatched_estimator_execution,
        AcquisitionProducts(wfs.execution.observation,
            mismatched_wfs_measurement))
    assert_plant_preparation_error(
        () -> execute_acquisition!(
            AcquisitionProducts(deepcopy(wfs.execution.observation)),
            wfs.path_result,
            observation_only_execution,
            Xoshiro(19),
        ),
        :acquisition,
        :prepared_binding,
    )
    assert_plant_preparation_error(
        () -> execute_acquisition!(
            AcquisitionProducts(wfs.products.observation,
                deepcopy(wfs.products.measurement)),
            wfs.path_result,
            wfs.execution,
            Xoshiro(20),
        ),
        :acquisition,
        :prepared_binding,
    )
    assert_plant_preparation_error(
        () -> execute_acquisition!(fast.products, fast.path_result, nothing,
            Xoshiro(21)),
        :acquisition,
        :unsupported_execution,
    )

    @test !applicable(PreparedAcquisitionOwner, fast.definition,
        fast.path_key, fast.path_result, fast.execution, fast.products)
    @test !applicable(PreparedPlant, plant.definition, plant.paths,
        plant.acquisitions)

    if coverage_instrumented()
        @test_skip "prepared-plant allocation assertions are disabled under coverage instrumentation"
    else
        @test prepared_path_execution_allocations(science_path) == 0
        @test prepared_path_execution_allocations(wfs_path) == 0
        @test prepared_acquisition_execution_allocations(
            fast, Xoshiro(14)) == 0
        @test prepared_acquisition_execution_allocations(
            wfs, Xoshiro(15)) == 0
    end

    result_before_rejection = copy(science_path.result.values)
    mismatches = (
        (MismatchedSourceGeometry(), :source_geometry),
        (MismatchedRadiometry(), :radiometry),
        (MismatchedRevision(), :revision),
        (MismatchedBackend(), :backend),
        (MismatchedDevice(), :device),
    )
    for (requirement, reason) in mismatches
        incompatible = AcquisitionDefinition(
            Symbol(:incompatible_, reason),
            :science,
            FramePlantAcquisitionModel(T(0.5), requirement),
        )
        assert_plant_preparation_error(
            () -> prepare_acquisition_owner(incompatible, science_path),
            :acquisition,
            reason,
        )
        @test science_path.result.values == result_before_rejection
    end

    additional_mismatches = (
        (() -> require_path_result(science_path;
            spectral_sampling=((wavelength_m=T(1), weight=T(1)),)),
            :spectral_sampling),
        (() -> require_path_result(science_path;
            optical_model=:different_model), :optical_model),
        (() -> require_path_result(science_path;
            sampling_contract=:different_sampling), :sampling_contract),
        (() -> require_path_result(science_path;
            propagation_model=:different_propagation), :propagation_model),
        (() -> require_path_result(science_path;
            output_plane=(kind=:different_plane,)), :output_plane),
    )
    for (requirement, reason) in additional_mismatches
        assert_plant_preparation_error(requirement, :acquisition, reason)
    end

    wrong_path_acquisition = AcquisitionDefinition(:wrong_path_binding,
        :wfs, FramePlantAcquisitionModel(T(0.5), MatchingPlantPath()))
    assert_plant_preparation_error(
        () -> PreparedAcquisitionOwner(wrong_path_acquisition, science_path,
            fast.execution, fast.products),
        :acquisition,
        :unknown_path,
    )

    assert_plant_preparation_error(
        () -> prepared_path(plant, :missing),
        :path,
        :unknown_id,
    )
    assert_plant_preparation_error(
        () -> prepared_acquisition(plant, :missing),
        :acquisition,
        :unknown_id,
    )

    unsupported_path = OpticalPathDefinition(:unsupported, science_source,
        UnsupportedPreparedPathModel())
    assert_plant_preparation_error(
        () -> prepare_path_executor(unsupported_path, telescope, atmosphere),
        :path,
        :unsupported_model,
    )
    unsupported_acquisition = AcquisitionDefinition(:unsupported,
        :science, UnsupportedPreparedAcquisitionModel())
    assert_plant_preparation_error(
        () -> prepare_acquisition_owner(unsupported_acquisition,
            science_path),
        :acquisition,
        :unsupported_model,
    )

    invalid_path = OpticalPathDefinition(:invalid_preparation, science_source,
        InvalidPreparedPathModel())
    assert_plant_preparation_error(
        () -> prepare_path_executor(invalid_path, telescope, atmosphere),
        :path,
        :invalid_preparation,
    )
    invalid_acquisition = AcquisitionDefinition(:invalid_preparation,
        :science, InvalidPreparedAcquisitionModel())
    assert_plant_preparation_error(
        () -> prepare_acquisition_owner(invalid_acquisition, science_path),
        :acquisition,
        :invalid_preparation,
    )

    stale_result = copy(science_path.result.values)
    set_pupil_reflectivity!(telescope, T(0.9))
    assert_plant_preparation_error(
        () -> execute_path!(science_path),
        :path,
        :revision,
    )
    @test science_path.result.values == stale_result
end
