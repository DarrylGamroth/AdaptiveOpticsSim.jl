# Focused preparation-only target-partition qualification.  This file is
# intentionally self-contained until the suite registry adopts it: run it with
# `include("test/runtests_head.jl"); include("test/testsets/plant_target_partitions.jl")`.

struct TargetPartitionTestBackend <: AbstractArrayBackend end

const TARGET_PARTITION_TEST_ACCELERATOR = AcceleratorComputeDevice(
    TargetPartitionTestBackend(), UInt32(1))

struct TargetPartitionTestArray{T,N} <: DenseArray{T,N}
    storage::Array{T,N}
end

Base.size(array::TargetPartitionTestArray) = size(array.storage)
Base.axes(array::TargetPartitionTestArray) = axes(array.storage)
Base.IndexStyle(::Type{<:TargetPartitionTestArray}) = IndexLinear()
Base.getindex(array::TargetPartitionTestArray, index::Int) =
    @inbounds array.storage[index]
Base.getindex(array::TargetPartitionTestArray, row::Int, column::Int) =
    @inbounds array.storage[row, column]
Base.setindex!(array::TargetPartitionTestArray, value, index::Int) =
    (@inbounds array.storage[index] = value; array)
Base.setindex!(array::TargetPartitionTestArray, value, row::Int,
    column::Int) = (@inbounds array.storage[row, column] = value; array)
Base.similar(::TargetPartitionTestArray, ::Type{T},
    dimensions::Dims{N}) where {T,N} =
    TargetPartitionTestArray(Array{T}(undef, dimensions))
Base.similar(array::TargetPartitionTestArray,
    dimensions::Dims{N}) where {N} = similar(array, eltype(array), dimensions)
Base.copyto!(destination::TargetPartitionTestArray,
    source::TargetPartitionTestArray) =
    (copyto!(destination.storage, source.storage); destination)

Backends.array_backend_selector(::Type{<:TargetPartitionTestArray}) =
    TargetPartitionTestBackend()
Backends.backend(::TargetPartitionTestArray) = TargetPartitionTestBackend()
Backends.compute_device(::TargetPartitionTestArray) =
    TARGET_PARTITION_TEST_ACCELERATOR
Backends.allocate_array(::TargetPartitionTestBackend, ::Type{T},
    dimensions::Vararg{Int,N}) where {T,N} =
    TargetPartitionTestArray(Array{T}(undef, dimensions...))
Backends.execution_style(::TargetPartitionTestArray) = ScalarCPUStyle()
Backends.compute_device_availability(
    ::AcceleratorComputeDevice{TargetPartitionTestBackend}) =
    ComputeDeviceAvailable()
Backends._with_compute_device(f::F,
    ::AcceleratorComputeDevice{TargetPartitionTestBackend}) where {F} = f()

struct TargetPartitionTestContext <:
    Backends._AbstractPreparedDeviceExecutionContext
    target::typeof(TARGET_PARTITION_TEST_ACCELERATOR)
end

Backends._prepare_device_execution_context(
    target::AcceleratorComputeDevice{TargetPartitionTestBackend}) =
    TargetPartitionTestContext(target)
Backends._prepared_device_execution_compute_device(
    context::TargetPartitionTestContext) = context.target
Backends._with_prepared_device_execution_context(f::F,
    ::TargetPartitionTestContext) where {F} = f()
Backends._synchronize_prepared_device_execution_context!(
    ::TargetPartitionTestContext) = nothing

Plant.structural_array_bytes(array::TargetPartitionTestArray,
    target::AbstractComputeDevice) = begin
    compute_device(array) == target || throw(StructuralResourceError(
        :array_storage, :wrong_device,
        "target-partition test array occupies $(compute_device(array)); expected $target"))
    UInt64(length(array)) * UInt64(sizeof(eltype(array)))
end

struct TargetPartitionTestTelescopeDefinition <: AbstractTelescopeDefinition
    revision::UInt
end

struct TargetPartitionTestTelescope{A,B<:AbstractArrayBackend} <:
    AbstractTelescope
    reflectivity::A
    revision::UInt
    selector::B
end

Optics.backend(telescope::TargetPartitionTestTelescope) = telescope.selector
Optics.pupil_reflectivity(telescope::TargetPartitionTestTelescope) =
    telescope.reflectivity
Optics.aperture_revision(telescope::TargetPartitionTestTelescope) =
    telescope.revision

function Optics._prepare_telescope(
    definition::TargetPartitionTestTelescopeDefinition,
    target::AbstractComputeDevice,
)
    selector = compute_device_backend(target)
    reflectivity = allocate_array(selector, Float64, 2, 2)
    fill!(reflectivity, 1.0)
    return TargetPartitionTestTelescope(
        reflectivity, definition.revision, selector)
end

function Optics.validate_telescope_target(
    telescope::TargetPartitionTestTelescope,
    target::AbstractComputeDevice,
)
    compute_device(telescope.reflectivity) == target || throw(
        InvalidConfiguration("target-partition test telescope is on the wrong device"))
    return telescope
end

function Plant.structural_resource_fact(
    telescope::TargetPartitionTestTelescope,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    compute_device(telescope.reflectivity) == target ||
        return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(
        id,
        target,
        structural_array_bytes(telescope.reflectivity, target),
        0,
    )
end

struct TargetPartitionTestAtmosphereDefinition <:
    AbstractTimedAtmosphereDefinition end

mutable struct TargetPartitionTestAtmosphereState
    timeline::AtmosphereTimelineState{Float64}
end

struct TargetPartitionTestAtmosphere{A} <: AbstractTimedAtmosphere
    identity::AtmosphereIdentity
    state::TargetPartitionTestAtmosphereState
    state_storage::A
end

function Atmospheres._prepare_timed_atmosphere(
    ::TargetPartitionTestAtmosphereDefinition,
    ::TargetPartitionTestTelescope,
    backend::AbstractArrayBackend,
)
    storage = allocate_array(backend, Float64, 2, 2)
    fill!(storage, 1.0)
    return TargetPartitionTestAtmosphere(
        AtmosphereIdentity(),
        TargetPartitionTestAtmosphereState(
            new_atmosphere_timeline(Float64)),
        storage,
    )
end

function Atmospheres.prepare_timed_atmosphere(
    definition::TargetPartitionTestAtmosphereDefinition,
    telescope::TargetPartitionTestTelescope,
    target::AbstractComputeDevice,
)
    selector = compute_device_backend(target)
    atmosphere = Atmospheres._prepare_timed_atmosphere(
        definition, telescope, selector)
    validate_timed_atmosphere_target(atmosphere, target)
    return atmosphere
end

function Atmospheres.validate_timed_atmosphere_target(
    atmosphere::TargetPartitionTestAtmosphere,
    target::AbstractComputeDevice,
)
    compute_device(atmosphere.state_storage) == target || throw(
        InvalidConfiguration("target-partition test atmosphere is on the wrong device"))
    return atmosphere
end

Plant.partition_atmosphere_layer_ids(
    ::TargetPartitionTestAtmosphereDefinition) = ()

function Plant.structural_resource_fact(
    atmosphere::TargetPartitionTestAtmosphere,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    compute_device(atmosphere.state_storage) == target ||
        return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(
        id,
        target,
        structural_array_bytes(atmosphere.state_storage, target),
        0,
    )
end

struct TargetPartitionTestPathModel
    label::Symbol
end

struct TargetPartitionTestAcquisitionModel
    label::Symbol
end

struct TargetPartitionTestExecution{I,R,W}
    input::I
    result::R
    workspace::W
end

Plant.plant_model_definition_style(::Type{TargetPartitionTestPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{TargetPartitionTestAcquisitionModel}) = ColdPlantModelDefinition()

function Plant.validate_path_execution_binding(
    execution::TargetPartitionTestExecution,
    input::IntensityMap,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "target-partition test execution lost its exact products"))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::TargetPartitionTestExecution,
    target::AbstractComputeDevice,
)
    all(storage -> compute_device(storage) == target,
        (execution.input.values, execution.result.values, execution.workspace)) ||
        throw(PlantPreparationError(:path, :wrong_device,
            "target-partition test execution storage is not co-located"))
    return execution
end

function Plant.structural_resource_fact(
    execution::TargetPartitionTestExecution,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    validate_path_execution_target(execution, target)
    resident = structural_array_bytes(execution.input.values, target) +
        structural_array_bytes(execution.result.values, target)
    workspace = structural_array_bytes(execution.workspace, target)
    return KnownStructuralResourceFact(id, target, resident, workspace)
end

function target_partition_test_intensity_map(
    backend::AbstractArrayBackend,
    label::Symbol,
)
    values = allocate_array(backend, Float64, 2, 2)
    fill!(values, Float64(length(String(label))))
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        values;
        coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        spectral=MonochromaticChannel(500e-9),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    return IntensityMap(metadata, values)
end

function Plant.prepare_target_local_path_resources(
    model::TargetPartitionTestPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::TargetPartitionTestTelescope,
    context,
)
    selector = backend(telescope)
    input = target_partition_test_intensity_map(selector, model.label)
    result = target_partition_test_intensity_map(selector, model.label)
    workspace = allocate_array(selector, Float64, 3)
    execution = TargetPartitionTestExecution(input, result, workspace)
    return PreparedTargetLocalPathResources(
        definition,
        source,
        telescope,
        input,
        result,
        execution;
        context,
        optical_model=(kind=:target_partition_test, label=model.label),
        propagation_model=:test_identity,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_target_local_acquisition_provider(
    model::TargetPartitionTestAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedTargetLocalPathResources,
)
    selector = path_result_key(path).backend
    product = allocate_array(selector, Float64, 5)
    fill!(product, Float64(length(String(model.label))))
    products = AcquisitionProducts(
        product;
        metadata=(kind=:target_partition_test_product, label=model.label),
    )
    return prepare_unchanged_synthetic_provider(products)
end

function target_partition_test_definition(; include_gamma::Bool=false)
    telescope = TargetPartitionTestTelescopeDefinition(UInt(1))
    source = Source(band=:I, magnitude=0.0)
    alpha = OpticalPathDefinition(
        :alpha, source, TargetPartitionTestPathModel(:alpha))
    beta = OpticalPathDefinition(
        :beta, source, TargetPartitionTestPathModel(:beta))
    gamma = OpticalPathDefinition(
        :gamma, source, TargetPartitionTestPathModel(:gamma))
    alpha_camera = AcquisitionDefinition(
        :alpha_camera, :alpha, TargetPartitionTestAcquisitionModel(:alpha))
    beta_camera = AcquisitionDefinition(
        :beta_camera, :beta, TargetPartitionTestAcquisitionModel(:beta))
    gamma_camera = AcquisitionDefinition(
        :gamma_camera, :gamma, TargetPartitionTestAcquisitionModel(:gamma))
    return PlantDefinition(
        ;
        telescope,
        atmosphere=TargetPartitionTestAtmosphereDefinition(),
        paths=include_gamma ? (alpha, beta, gamma) : (alpha, beta),
        acquisitions=include_gamma ?
            (alpha_camera, beta_camera, gamma_camera) :
            (alpha_camera, beta_camera),
    )
end

function target_partition_test_error(f, component::Symbol, reason::Symbol)
    error = try
        f()
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function target_partition_test_path_ids(partition)
    return Tuple(path_id(path.definition) for path in prepared_paths(partition))
end

function target_partition_test_acquisition_ids(partition)
    return Tuple(acquisition_id(acquisition.definition) for acquisition in
        prepared_acquisitions(partition))
end

function target_partition_test_byte_oracle(partition; authority::Bool)
    telescope_bytes = UInt64(4 * sizeof(Float64))
    atmosphere_bytes = authority ? UInt64(4 * sizeof(Float64)) : UInt64(0)
    path_resident = UInt64(2 * 2 * 2 * sizeof(Float64))
    path_workspace = UInt64(3 * sizeof(Float64))
    acquisition_resident = UInt64(5 * sizeof(Float64))
    count = UInt64(length(prepared_paths(partition)))
    return (
        resident=telescope_bytes + atmosphere_bytes +
            count * (path_resident + acquisition_resident),
        workspace=count * path_workspace,
    )
end

@testset "Preparation-only target-local partitions" begin
    definition = target_partition_test_definition()
    host = HostComputeDevice()
    accelerator = TARGET_PARTITION_TEST_ACCELERATOR
    cpu_assignment = resolve_plant_partition_assignment(
        definition, host, :beta => host, :alpha => host)
    accelerator_assignment = resolve_plant_partition_assignment(
        definition, accelerator, :alpha => accelerator, :beta => accelerator)
    mixed_assignment = resolve_plant_partition_assignment(
        definition, host, :beta => accelerator, :alpha => host)

    # This internal-token constructor is deliberately the only constructor;
    # callers cannot fabricate an all-fields assignment value.
    @test_throws MethodError ResolvedPlantPartitionAssignment(
        definition, nothing, nothing, (), UInt8(1), ())

    cpu = prepare_plant_partitions(definition, cpu_assignment; run_seed=71)
    accelerator_only = prepare_plant_partitions(
        definition, accelerator_assignment; run_seed=71)
    mixed = prepare_plant_partitions(definition, mixed_assignment; run_seed=71)

    @test length(prepared_partitions(cpu)) == 1
    @test length(prepared_partitions(accelerator_only)) == 1
    @test length(prepared_partitions(mixed)) == 2
    @test compute_device(mixed.authority) == host
    @test atmosphere_authority_target(mixed.authority) == host
    @test atmosphere_authority_identity(mixed.authority) ===
        atmosphere_identity(prepared_atmosphere(mixed.authority))
    @test all(partition -> atmosphere_authority_binding(partition) ===
        atmosphere_authority_binding(mixed.authority), prepared_partitions(mixed))
    @test length(unique(objectid.([
        atmosphere_authority_identity(atmosphere_authority_binding(partition))
        for partition in prepared_partitions(mixed)]))) == 1

    cpu_partition = only(prepared_partitions(cpu))
    accelerator_partition = only(prepared_partitions(accelerator_only))
    mixed_host = prepared_partition(mixed, host)
    mixed_accelerator = prepared_partition(mixed, accelerator)
    @test target_partition_test_path_ids(cpu_partition) ==
        (OpticalPathID(:alpha), OpticalPathID(:beta))
    @test target_partition_test_acquisition_ids(cpu_partition) ==
        (AcquisitionID(:alpha_camera), AcquisitionID(:beta_camera))
    @test target_partition_test_path_ids(mixed_host) == (OpticalPathID(:alpha),)
    @test target_partition_test_acquisition_ids(mixed_host) ==
        (AcquisitionID(:alpha_camera),)
    @test target_partition_test_path_ids(mixed_accelerator) ==
        (OpticalPathID(:beta),)
    @test target_partition_test_acquisition_ids(mixed_accelerator) ==
        (AcquisitionID(:beta_camera),)

    for partition in prepared_partitions(mixed)
        @test !hasproperty(partition, :atmosphere)
        @test !hasproperty(partition, :timeline)
        for path in prepared_paths(partition)
            @test compute_device(path_input(path).values) == compute_device(partition)
            @test compute_device(path_result(path).values) == compute_device(partition)
            @test compute_device(path.execution.workspace) == compute_device(partition)
            @test !hasproperty(path, :atmosphere)
            @test !hasproperty(path, :timeline)
        end
        for acquisition in prepared_acquisitions(partition)
            @test compute_device(acquisition_products(
                acquisition_provider(acquisition)).observation) ==
                compute_device(partition)
        end
    end
    @test compute_device(accelerator_partition) == accelerator
    @test all(path -> path_input(path).values isa TargetPartitionTestArray,
        prepared_paths(accelerator_partition))

    @test rng_replay_metadata(cpu) == rng_replay_metadata(accelerator_only)
    @test rng_replay_metadata(cpu) == rng_replay_metadata(mixed)
    @test rng_replay_metadata(mixed).run_seed == UInt64(71)
    @test length(rng_replay_metadata(mixed).owners) == 5

    for (partition, authority) in ((cpu_partition, true),
            (accelerator_partition, true), (mixed_host, true),
            (mixed_accelerator, false))
        report = partition_resource_report(partition)
        oracle = target_partition_test_byte_oracle(partition; authority)
        @test structural_resident_bytes(report) == oracle.resident
        @test structural_workspace_bytes(report) == oracle.workspace
        @test all(structural_resource_known, structural_resource_facts(report))
        @test !any(fact -> structural_resource_owner_id(fact).category ==
            :renderer, structural_resource_facts(report))
    end

    stale_definition = target_partition_test_definition()
    foreign_definition = target_partition_test_definition(include_gamma=true)
    target_partition_test_error(() -> prepare_plant_partitions(
        stale_definition, cpu_assignment; run_seed=71),
        :partition_assignment, :stale_assignment)
    target_partition_test_error(() -> prepare_plant_partitions(
        foreign_definition, cpu_assignment; run_seed=71),
        :partition_assignment, :foreign_assignment)
    target_partition_test_error(() -> prepared_partition(mixed,
        AcceleratorComputeDevice(TargetPartitionTestBackend(), UInt32(2))),
        :partition, :unknown_target)
end
