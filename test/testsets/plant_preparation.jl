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

for model in (
    DirectSciencePathModel,
    ShackHartmannPlantPathModel,
    FramePlantAcquisitionModel,
    ContractWFSPlantAcquisitionModel,
    UnsupportedPreparedPathModel,
    UnsupportedPreparedAcquisitionModel,
)
    @eval AdaptiveOpticsSim.plant_model_definition_style(
        ::Type{<:$model},
    ) = ColdPlantModelDefinition()
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

    stale_result = copy(science_path.result.values)
    set_pupil_reflectivity!(telescope, T(0.9))
    assert_plant_preparation_error(
        () -> execute_path!(science_path),
        :path,
        :revision,
    )
    @test science_path.result.values == stale_result
end
