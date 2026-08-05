# Focused preparation-only target-partition qualification.

struct TargetPartitionTestBackend <: AbstractArrayBackend end

const TARGET_PARTITION_TEST_ACCELERATOR = AcceleratorComputeDevice(
    TargetPartitionTestBackend(), UInt32(1))
const TARGET_PARTITION_TEST_AVAILABLE = Ref(true)
const TARGET_PARTITION_TEST_MISPLACE_ALLOCATIONS = Ref(false)

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
Backends.array_backend_type(::TargetPartitionTestBackend) =
    TargetPartitionTestArray
Backends.backend(::TargetPartitionTestArray) = TargetPartitionTestBackend()
Backends.compute_device(::TargetPartitionTestArray) =
    TARGET_PARTITION_TEST_ACCELERATOR
function Backends.allocate_array(
    ::TargetPartitionTestBackend,
    ::Type{T},
    dimensions::Vararg{Int,N},
) where {T,N}
    TARGET_PARTITION_TEST_MISPLACE_ALLOCATIONS[] &&
        return Array{T}(undef, dimensions...)
    return TargetPartitionTestArray(Array{T}(undef, dimensions...))
end
Backends.execution_style(::TargetPartitionTestArray) = ScalarCPUStyle()
Backends.compute_device_availability(
    ::AcceleratorComputeDevice{TargetPartitionTestBackend}) =
    TARGET_PARTITION_TEST_AVAILABLE[] ? ComputeDeviceAvailable() :
    ComputeDeviceUnavailable(:test_device_offline)
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

struct TargetPartitionTestOpticModel end

struct TargetPartitionTestPreparedOptic{D<:AbstractComputeDevice}
    endpoint::CommandEndpointID
    target::D
end

mutable struct TargetPartitionTestOpticState{A<:AbstractArray}
    active::A
end

mutable struct TargetPartitionTestOpticWorkspace{A<:AbstractArray}
    staging::A
end

struct TargetPartitionTestExecution{I,R,W}
    input::I
    result::R
    workspace::W
end

struct TargetPartitionTestPupilExecution{I,R,W}
    input::I
    result::R
    workspace::W
end

Plant.plant_model_definition_style(::Type{TargetPartitionTestPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{TargetPartitionTestAcquisitionModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(::Type{TargetPartitionTestOpticModel}) =
    ColdPlantModelDefinition()

function Plant.prepare_target_local_controllable_optic(
    ::TargetPartitionTestOpticModel,
    definition::ControllableOpticDefinition,
    ::TargetPartitionTestTelescope,
    ::TargetPartitionTestAtmosphereDefinition,
    target::AbstractComputeDevice,
)
    return TargetPartitionTestPreparedOptic(
        only(command_endpoint_ids(definition)), target)
end

function Plant.prepare_controllable_optic_state(
    prepared::TargetPartitionTestPreparedOptic,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    only(endpoint_ids) == prepared.endpoint || throw(
        PlantPreparationError(:controllable_optic, :endpoint_mismatch,
            "target-partition test optic endpoint changed"))
    return TargetPartitionTestOpticState(only(initial_commands))
end


function Plant.prepare_controllable_optic_workspace(
    prepared::TargetPartitionTestPreparedOptic,
)
    staging = allocate_device_array(prepared.target, Float64, 1)
    fill!(staging, 0.0)
    return TargetPartitionTestOpticWorkspace(staging)
end


function Plant.stage_controllable_optic_command!(
    prepared::TargetPartitionTestPreparedOptic,
    ::TargetPartitionTestOpticState,
    workspace::TargetPartitionTestOpticWorkspace,
    endpoint::CommandEndpointID,
    command::AbstractVector{Float64},
    ::PlantTimestamp,
)
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "target-partition test optic received another endpoint"))
    copyto!(workspace.staging, command)
    return nothing
end


function Plant.commit_controllable_optic_command!(
    ::TargetPartitionTestPreparedOptic,
    state::TargetPartitionTestOpticState{A},
    workspace::TargetPartitionTestOpticWorkspace{A},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {A}
    previous = state.active
    state.active = workspace.staging
    workspace.staging = previous
    return nothing
end


function Plant.validate_controllable_optic_target(
    prepared::TargetPartitionTestPreparedOptic,
    target::AbstractComputeDevice,
)
    prepared.target == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "target-partition test optic occupies another target"))
    return prepared
end


function Plant.validate_controllable_optic_state_target(
    ::TargetPartitionTestPreparedOptic,
    state::TargetPartitionTestOpticState,
    target::AbstractComputeDevice,
)
    compute_device(state.active) == target || throw(PlantPreparationError(
        :controllable_optic, :wrong_device,
        "target-partition test optic state occupies another target"))
    return state
end


function Plant.validate_controllable_optic_workspace_target(
    ::TargetPartitionTestPreparedOptic,
    workspace::TargetPartitionTestOpticWorkspace,
    target::AbstractComputeDevice,
)
    compute_device(workspace.staging) == target || throw(
        PlantPreparationError(:controllable_optic, :wrong_device,
            "target-partition test optic workspace occupies another target"))
    return workspace
end


function Plant.structural_resource_fact(
    ::TargetPartitionTestPreparedOptic,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(id, target, 0, 0)
end


function Plant.structural_resource_fact(
    state::TargetPartitionTestOpticState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(
        id, target, structural_array_bytes(state.active, target), 0)
end


function Plant.structural_resource_fact(
    workspace::TargetPartitionTestOpticWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(
        id, target, 0, structural_array_bytes(workspace.staging, target))
end

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

function Plant.validate_path_execution_binding(
    execution::TargetPartitionTestPupilExecution,
    input::PupilFunction,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "target-partition pupil execution lost its exact products"))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::TargetPartitionTestPupilExecution,
    target::AbstractComputeDevice,
)
    all(storage -> compute_device(storage) == target,
        (execution.input.support,
            execution.input.amplitude,
            execution.input.opd,
            execution.result.values,
            execution.workspace)) || throw(PlantPreparationError(
        :path, :wrong_device,
        "target-partition pupil execution storage is not co-located"))
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

function target_partition_test_pupil_path(
    partition::PreparedTargetPartition,
    path::PreparedTargetLocalPathResources,
)
    telescope = prepared_telescope(partition)
    selector = backend(telescope)
    support = allocate_array(selector, Bool, 2, 2)
    amplitude = allocate_array(selector, Float64, 2, 2)
    opd = allocate_array(selector, Float64, 2, 2)
    fill!(support, true)
    fill!(amplitude, 1.0)
    fill!(opd, 0.0)
    metadata = OpticalPlaneMetadata(
        PupilPlane(),
        opd;
        coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination(),
    )
    input = PupilFunction{
        typeof(metadata),
        typeof(support),
        typeof(amplitude),
        typeof(opd),
        typeof(selector),
    }(metadata, support, amplitude, opd, aperture_revision(telescope))
    result = target_partition_test_intensity_map(
        selector, path_id(path.definition).name)
    workspace = allocate_array(selector, Float64, 3)
    execution = TargetPartitionTestPupilExecution(
        input, result, workspace)
    return PreparedTargetLocalPathResources(
        path.definition,
        path.source,
        telescope,
        input,
        result,
        execution;
        context=getfield(partition, :context),
        optical_model=(kind=:target_partition_pupil_coupling_test,
            path=path_id(path.definition)),
        propagation_model=:test_identity,
        model_revisions=UInt(1),
    )
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

function target_partition_test_sampled_aberration()
    opd = fill(2e-9, 2, 2)
    metadata = OpticalPlaneMetadata(
        PupilPlane(),
        opd;
        coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
    )
    return SampledAberrationDefinition(
        :partition_static_opd,
        OPDMap(opd),
        metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMAdditive(),
    )
end

function target_partition_test_controllable_optic(;
    visibility::AbstractPathVisibility=AllPathVisibility(),
    silence_policy::CommandSilencePolicy=CommandSilencePolicy(),
)
    schema = PlantCommandSchema(
        Float64,
        (1,);
        id=:partition_test_schema,
        version=1,
        endpoint=:partition_test_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :partition_test_actuator),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy,
    )
    return ControllableOpticDefinition(
        :partition_test_optic,
        TargetPartitionTestOpticModel(),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility,
    )
end

target_partition_test_command_endpoints() = (
    CommandEndpointConfiguration(
        :partition_test_command, [0.0]; capacity=2),
)

function target_partition_test_definition(; include_gamma::Bool=false,
    include_sampled_aberration::Bool=false,
    optic_visibility::AbstractPathVisibility=AllPathVisibility(),
    command_silence_policy::CommandSilencePolicy=CommandSilencePolicy(),
)
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
    sampled_aberrations = include_sampled_aberration ?
        (target_partition_test_sampled_aberration(),) : ()
    return PlantDefinition(
        ;
        telescope,
        atmosphere=TargetPartitionTestAtmosphereDefinition(),
        controllable_optics=(
            target_partition_test_controllable_optic(
                visibility=optic_visibility,
                silence_policy=command_silence_policy,
            ),),
        paths=include_gamma ? (alpha, beta, gamma) : (alpha, beta),
        acquisitions=include_gamma ?
            (alpha_camera, beta_camera, gamma_camera) :
            (alpha_camera, beta_camera),
        sampled_aberrations,
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
    return Tuple(path_id(path) for path in prepared_paths(partition))
end

function target_partition_test_acquisition_ids(partition)
    return Tuple(acquisition_id(acquisition.definition) for acquisition in
        prepared_acquisitions(partition))
end

function target_partition_test_byte_oracle(partition; authority::Bool,
    sampled_aberration_count::Integer=0)
    telescope_bytes = UInt64(4 * sizeof(Float64))
    atmosphere_bytes = authority ? UInt64(4 * sizeof(Float64)) : UInt64(0)
    sampled_aberration_bytes = UInt64(sampled_aberration_count) *
        UInt64(4 * sizeof(Float64))
    path_resident = UInt64(2 * 2 * 2 * sizeof(Float64))
    path_workspace = UInt64(3 * sizeof(Float64))
    acquisition_resident = UInt64(5 * sizeof(Float64))
    count = UInt64(length(prepared_paths(partition)))
    optic_count = UInt64(length(
        target_local_controllable_optic_owners(partition)))
    endpoint_count = UInt64(sum(
        owner -> length(target_local_command_endpoints(
            prepared_target_local_controllable_optic(owner))),
        target_local_controllable_optic_owners(partition);
        init=0,
    ))
    command_bytes = UInt64(sizeof(Float64))
    return (
        resident=telescope_bytes + atmosphere_bytes +
            sampled_aberration_bytes +
            count * (path_resident + acquisition_resident) +
            (optic_count + endpoint_count) * command_bytes,
        workspace=count * path_workspace +
            (optic_count + endpoint_count) * command_bytes,
    )
end

function target_partition_test_fresh_token_guards(
    assignment::ResolvedPlantPartitionAssignment,
    prepared::PreparedPlantPartitions,
)
    partition = first(prepared_partitions(prepared))
    path = first(prepared_paths(partition))
    acquisition = first(prepared_acquisitions(partition))
    authority = prepared_atmosphere_authority(prepared)
    command_authority = prepared_command_authority(prepared)
    binding = atmosphere_authority_binding(prepared)

    @test_throws ArgumentError PreparedTargetLocalPathResources(
        Plant._PreparedTargetLocalPathResourcesToken(),
        getfield(path, :definition),
        getfield(path, :source),
        getfield(path, :telescope),
        getfield(path, :context),
        getfield(path, :input),
        getfield(path, :result),
        getfield(path, :execution),
        getfield(path, :key),
    )
    @test_throws ArgumentError PreparedTargetLocalAcquisitionResources(
        Plant._PreparedTargetLocalAcquisitionResourcesToken(),
        getfield(acquisition, :definition),
        getfield(acquisition, :path_key),
        getfield(acquisition, :path_result),
        getfield(acquisition, :context),
        getfield(acquisition, :provider),
    )
    @test_throws ArgumentError AtmosphereAuthorityBinding(
        Plant._AtmosphereAuthorityBindingToken(),
        getfield(binding, :target),
        getfield(binding, :identity),
    )
    @test_throws ArgumentError PreparedAtmosphereAuthority(
        Plant._PreparedAtmosphereAuthorityToken(),
        getfield(authority, :definition_identity),
        getfield(authority, :target),
        getfield(authority, :context),
        getfield(authority, :telescope),
        getfield(authority, :atmosphere),
        getfield(authority, :rngs),
        getfield(authority, :binding),
    )
    @test_throws ArgumentError PreparedCommandAuthority(
        Plant._PreparedCommandAuthorityToken(),
        getfield(command_authority, :binding),
        getfield(command_authority, :identity),
        getfield(command_authority, :target),
        getfield(command_authority, :context),
        getfield(command_authority, :optic_definitions),
        getfield(command_authority, :endpoints),
    )
    @test_throws ArgumentError PreparedTargetPartition(
        Plant._PreparedTargetPartitionToken(),
        getfield(partition, :target),
        getfield(partition, :context),
        getfield(partition, :telescope),
        getfield(partition, :authority_binding),
        getfield(partition, :sampled_aberrations),
        getfield(partition, :paths),
        getfield(partition, :acquisitions),
        getfield(partition, :rngs),
        getfield(partition, :controllable_optics),
        getfield(partition, :resource_report),
        getfield(partition, :controllable_optic_ids),
        getfield(partition, :command_endpoint_ids),
    )
    @test_throws ArgumentError PreparedPlantPartitions(
        Plant._PreparedPlantPartitionsToken(),
        getfield(prepared, :definition),
        assignment,
        authority,
        prepared_command_authority(prepared),
        getfield(prepared, :partitions),
        getfield(prepared, :run_seed),
        getfield(prepared, :rng_derivation_version),
    )
    return nothing
end

@testset "Preparation-only target-local partitions" begin
    definition = target_partition_test_definition()
    host = HostComputeDevice()
    accelerator = TARGET_PARTITION_TEST_ACCELERATOR
    cpu_assignment = @inferred resolve_plant_partition_assignment(
        definition, host, :beta => host, :alpha => host)
    accelerator_assignment = @inferred resolve_plant_partition_assignment(
        definition, accelerator, :alpha => accelerator, :beta => accelerator)
    mixed_assignment = @inferred resolve_plant_partition_assignment(
        definition, host, :beta => accelerator, :alpha => host)
    mixed_accelerator_authority_assignment =
        @inferred resolve_plant_partition_assignment(
            definition, accelerator, :beta => accelerator, :alpha => host)

    # This internal-token constructor is deliberately the only constructor;
    # callers cannot fabricate an all-fields assignment value.
    @test_throws MethodError ResolvedPlantPartitionAssignment(
        definition, nothing, nothing, (), UInt8(1), ())

    command_endpoints = target_partition_test_command_endpoints()
    cpu = prepare_plant_partitions(
        definition, cpu_assignment;
        run_seed=71, command_authority_target=host, command_endpoints)
    accelerator_only = prepare_plant_partitions(
        definition, accelerator_assignment;
        run_seed=71, command_authority_target=accelerator, command_endpoints)
    mixed = prepare_plant_partitions(
        definition, mixed_assignment;
        run_seed=71, command_authority_target=host, command_endpoints)
    mixed_accelerator_authority = prepare_plant_partitions(
        mixed_accelerator_authority_assignment;
        run_seed=71,
        command_authority_target=accelerator,
        command_endpoints,
    )
    accelerator_paths_host_command_authority = prepare_plant_partitions(
        definition,
        accelerator_assignment;
        run_seed=71,
        command_authority_target=host,
        command_endpoints,
    )
    target_partition_test_fresh_token_guards(mixed_assignment, mixed)

    @test length(prepared_partitions(cpu)) == 1
    @test length(prepared_partitions(accelerator_only)) == 1
    @test length(prepared_partitions(mixed)) == 2
    @test compute_device(prepared_atmosphere_authority(mixed)) == host
    @test compute_device(prepared_atmosphere_authority(
        mixed_accelerator_authority)) == accelerator
    @test atmosphere_authority_target(mixed) == host
    @test atmosphere_authority_target(mixed_accelerator_authority) ==
        accelerator
    @test command_authority_target(cpu) == host
    @test command_authority_target(accelerator_only) == accelerator
    @test command_authority_target(mixed) == host
    @test command_authority_target(mixed_accelerator_authority) == accelerator
    @test atmosphere_authority_target(
        accelerator_paths_host_command_authority) == accelerator
    @test command_authority_target(
        accelerator_paths_host_command_authority) == host
    @test all(partition -> compute_device(partition) == accelerator,
        prepared_partitions(accelerator_paths_host_command_authority))
    @test prepared_command_authority(mixed) isa PreparedCommandAuthority
    @test command_authority_identity(prepared_command_authority(mixed)) ===
        command_authority_identity(mixed)

    authority_plan = prepared_command_authority(mixed)
    authority_state = @inferred CommandAuthorityState(authority_plan)
    authority_workspace = @inferred CommandAuthorityWorkspace(authority_plan)
    endpoint_id = CommandEndpointID(:partition_test_command)
    authority_endpoint =
        prepared_command_authority_endpoint(authority_plan, endpoint_id)
    @test authority_endpoint isa PreparedCommandEndpoint
    @test command_endpoint_id(authority_endpoint) == endpoint_id
    @test command_authority_endpoint_state(
        authority_plan, authority_state, endpoint_id) isa CommandEndpointState
    application_state = command_authority_application_state(
        authority_plan, authority_state, endpoint_id)
    @test effective_command(application_state) == [0.0]
    @test command_authority_disposition_workspace(
        authority_plan, authority_workspace, endpoint_id) isa
        CommandDispositionWorkspace
    @test !command_authority_failed(authority_state)

    accelerator_authority_plan =
        prepared_command_authority(accelerator_only)
    accelerator_authority_state =
        @inferred CommandAuthorityState(accelerator_authority_plan)
    accelerator_effective = effective_command(
        command_authority_application_state(
            accelerator_authority_plan,
            accelerator_authority_state,
            endpoint_id,
        ))
    @test accelerator_effective isa TargetPartitionTestArray
    @test compute_device(accelerator_effective) == accelerator

    caller_initial = initial_effective_command(only(command_endpoints))
    caller_initial[1] = 99.0
    isolated_state = CommandAuthorityState(authority_plan)
    @test effective_command(command_authority_application_state(
        authority_plan, isolated_state, endpoint_id)) == [0.0]
    caller_initial[1] = 0.0

    safe_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(5),
    )
    wrapped_definition = target_partition_test_definition(
        command_silence_policy=safe_policy)
    wrapped_assignment = resolve_plant_partition_assignment(
        wrapped_definition, host, :alpha => host, :beta => host)
    initial_parent = [-0.5, 9.0]
    safe_parent = [0.25, 9.0]
    wrapped_endpoints = (
        CommandEndpointConfiguration(
            :partition_test_command,
            @view(initial_parent[1:1]);
            capacity=2,
            safe_command=@view(safe_parent[1:1]),
        ),)
    wrapped = prepare_plant_partitions(
        wrapped_definition,
        wrapped_assignment;
        run_seed=72,
        command_authority_target=host,
        command_endpoints=wrapped_endpoints,
    )
    initial_parent[1] = 8.0
    safe_parent[1] = 7.0
    wrapped_authority = prepared_command_authority(wrapped)
    wrapped_state = CommandAuthorityState(wrapped_authority)
    wrapped_application = command_authority_application_state(
        wrapped_authority, wrapped_state, endpoint_id)
    @test effective_command(wrapped_application) == [-0.5]
    wrapped_replica = only(target_local_controllable_optic_owners(
        only(prepared_partitions(wrapped))))
    @test effective_command(target_local_command_endpoint_state(
        wrapped_replica, endpoint_id)) == [-0.5]
    wrapped_endpoint = prepared_command_authority_endpoint(
        wrapped_authority, endpoint_id)
    wrapped_endpoint_state = command_authority_endpoint_state(
        wrapped_authority, wrapped_state, endpoint_id)
    wrapped_dispositions = CommandDispositionWorkspace(wrapped_endpoint)
    apply_command_silence_transition!(
        wrapped_dispositions,
        wrapped_endpoint,
        wrapped_endpoint_state,
        wrapped_application,
        PlantTimestamp(5),
    )
    @test effective_command(wrapped_application) == [0.25]

    foreign_authority = prepared_command_authority(cpu)
    target_partition_test_error(() -> command_authority_endpoint_state(
        foreign_authority, authority_state, endpoint_id),
        :command_authority, :foreign_state)
    target_partition_test_error(() ->
        command_authority_disposition_workspace(
            foreign_authority, authority_workspace, endpoint_id),
        :command_authority, :foreign_workspace)
    @test atmosphere_authority_identity(mixed) === atmosphere_identity(
        prepared_atmosphere(prepared_atmosphere_authority(mixed)))
    @test all(partition -> atmosphere_authority_binding(partition) ===
        atmosphere_authority_binding(mixed), prepared_partitions(mixed))
    @test all(partition -> atmosphere_authority_binding(partition) ===
        atmosphere_authority_binding(mixed_accelerator_authority),
        prepared_partitions(mixed_accelerator_authority))
    @test length(unique(objectid.([
        atmosphere_authority_identity(atmosphere_authority_binding(partition))
        for partition in prepared_partitions(mixed)]))) == 1

    cpu_partition = only(prepared_partitions(cpu))
    accelerator_partition = only(prepared_partitions(accelerator_only))
    mixed_host = prepared_partition(mixed, host)
    mixed_accelerator = prepared_partition(mixed, accelerator)
    mixed_accelerator_authority_host = prepared_partition(
        mixed_accelerator_authority, host)
    mixed_accelerator_authority_accelerator = prepared_partition(
        mixed_accelerator_authority, accelerator)
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
    @test target_partition_test_path_ids(
        mixed_accelerator_authority_host) ==
        target_partition_test_path_ids(mixed_host)
    @test target_partition_test_path_ids(
        mixed_accelerator_authority_accelerator) ==
        target_partition_test_path_ids(mixed_accelerator)
    @test Tuple(partition_controllable_optic_ids(cpu_partition)) ==
        (ControllableOpticID(:partition_test_optic),)
    @test Tuple(partition_command_endpoint_ids(cpu_partition)) ==
        (CommandEndpointID(:partition_test_command),)
    @test partition_controllable_optic_ids(accelerator_partition) ==
        partition_controllable_optic_ids(cpu_partition)
    @test partition_command_endpoint_ids(accelerator_partition) ==
        partition_command_endpoint_ids(cpu_partition)
    @test partition_controllable_optic_ids(mixed_host) ==
        partition_controllable_optic_ids(cpu_partition)
    @test partition_command_endpoint_ids(mixed_host) ==
        partition_command_endpoint_ids(cpu_partition)
    @test partition_controllable_optic_ids(mixed_accelerator) ==
        partition_controllable_optic_ids(cpu_partition)
    @test partition_command_endpoint_ids(mixed_accelerator) ==
        partition_command_endpoint_ids(cpu_partition)
    @test all(partition -> length(
            target_local_controllable_optic_owners(partition)) == 1,
        prepared_partitions(mixed))
    @test length(prepared_paths(cpu_partition)) == 2
    @test length(target_local_controllable_optic_owners(cpu_partition)) == 1
    @test command_authority_identity(mixed) !==
        command_authority_identity(cpu)
    @test all(partition -> command_authority_identity(
            only(target_local_command_endpoints(
                prepared_target_local_controllable_optic(only(
                    target_local_controllable_optic_owners(partition)))))) ===
            command_authority_identity(mixed),
        prepared_partitions(mixed))

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

    selected_definition = target_partition_test_definition(
        optic_visibility=SelectedPathVisibility(:alpha))
    selected_assignment = resolve_plant_partition_assignment(
        selected_definition, host, :alpha => host, :beta => accelerator)
    selected = prepare_plant_partitions(
        selected_definition, selected_assignment;
        run_seed=71,
        command_authority_target=host,
        command_endpoints,
    )
    selected_host = prepared_partition(selected, host)
    selected_accelerator = prepared_partition(selected, accelerator)
    @test length(target_local_controllable_optic_owners(selected_host)) == 1
    @test isempty(target_local_controllable_optic_owners(
        selected_accelerator))
    @test Tuple(partition_controllable_optic_ids(selected_host)) ==
        (ControllableOpticID(:partition_test_optic),)
    @test isempty(partition_controllable_optic_ids(selected_accelerator))
    @test Tuple(partition_command_endpoint_ids(selected_host)) ==
        (CommandEndpointID(:partition_test_command),)
    @test isempty(partition_command_endpoint_ids(selected_accelerator))
    for (partition, authority, optic_count) in (
            (selected_host, true, 1),
            (selected_accelerator, false, 0))
        report = partition_resource_report(partition)
        oracle = target_partition_test_byte_oracle(partition; authority)
        @test structural_resident_bytes(report) == oracle.resident
        @test structural_workspace_bytes(report) == oracle.workspace
        @test count(fact -> structural_resource_owner_id(fact).category ===
                :controllable_optic_replica,
            structural_resource_facts(report)) == optic_count
        @test count(fact -> structural_resource_owner_id(fact).category ===
                :effective_command_replica,
            structural_resource_facts(report)) == optic_count
    end

    sampled_definition = target_partition_test_definition(
        include_sampled_aberration=true)
    sampled_assignment = resolve_plant_partition_assignment(
        sampled_definition, host, :alpha => host, :beta => accelerator)
    sampled = prepare_plant_partitions(
        sampled_definition, sampled_assignment;
        run_seed=71,
        command_authority_target=host,
        command_endpoints,
    )
    caller_opd = surface_opd(sampled_aberration_surface(
        only(sampled_aberration_definitions(sampled_definition))))
    for partition in prepared_partitions(sampled)
        aberration = only(prepared_sampled_aberrations(partition))
        prepared_opd = sampled_aberration_opd(aberration)
        @test compute_device(prepared_opd) == compute_device(partition)
        @test prepared_opd !== caller_opd
        @test Array(prepared_opd) == caller_opd

        pupil_path = target_partition_test_pupil_path(
            partition, first(prepared_paths(partition)))
        bindings = Plant._prepare_sampled_aberration_path_bindings(
            prepared_sampled_aberrations(partition), [pupil_path])
        @test prepared_sampled_aberration_binding_count(bindings) == 1
        sampled_coupling = prepared_sampled_aberration_path_coupling(
            bindings, 1)
        @test sampled_coupling isa
            PreparedIdentityPupilFootprintCoupling
        @test getfield(sampled_coupling, :destination) ===
            path_input(pupil_path)
        apply_sampled_pupil_surface!(
            path_input(pupil_path),
            prepared_opd,
            sampled_coupling,
            sampled_aberration_application(aberration.definition),
        )
        @test all(path_input(pupil_path).opd .== 2e-9)

        optic = prepared_target_local_controllable_optic(only(
            target_local_controllable_optic_owners(partition)))
        @test controllable_optic_placement(optic) ===
            controllable_optic_placement(optic.definition)
        @test controllable_optic_visibility(optic) ===
            controllable_optic_visibility(optic.definition)
        optic_coupling = prepare_controllable_optic_path_coupling(
            optic, pupil_path)
        @test optic_coupling isa PreparedDirectPupilSurfaceCoupling
        @test getfield(optic_coupling, :destination) ===
            path_input(pupil_path)

        report = partition_resource_report(partition)
        oracle = target_partition_test_byte_oracle(
            partition;
            authority=compute_device(partition) == host,
            sampled_aberration_count=1,
        )
        @test structural_resident_bytes(report) == oracle.resident
        @test structural_workspace_bytes(report) == oracle.workspace
        @test count(fact ->
                structural_resource_owner_id(fact).category ===
                    :sampled_aberration,
            structural_resource_facts(report)) == 1
    end

    sampled_host = prepared_partition(sampled, host)
    sampled_accelerator = prepared_partition(sampled, accelerator)
    host_optic = prepared_target_local_controllable_optic(only(
        target_local_controllable_optic_owners(sampled_host)))
    accelerator_pupil_path = target_partition_test_pupil_path(
        sampled_accelerator,
        first(prepared_paths(sampled_accelerator)),
    )
    target_partition_test_error(
        () -> prepare_controllable_optic_path_coupling(
            host_optic, accelerator_pupil_path),
        :controllable_optic,
        :wrong_device,
    )

    @test rng_replay_metadata(cpu) == rng_replay_metadata(accelerator_only)
    @test rng_replay_metadata(cpu) == rng_replay_metadata(mixed)
    @test rng_replay_metadata(cpu) ==
        rng_replay_metadata(mixed_accelerator_authority)
    @test rng_replay_metadata(mixed).run_seed == UInt64(71)
    @test length(rng_replay_metadata(mixed).owners) == 5

    for (partition, authority) in ((cpu_partition, true),
            (accelerator_partition, true), (mixed_host, true),
            (mixed_accelerator, false),
            (mixed_accelerator_authority_host, false),
            (mixed_accelerator_authority_accelerator, true))
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
    grown_assignment = @inferred resolve_plant_partition_assignment(
        foreign_definition,
        host,
        :alpha => host,
        :beta => host,
        :gamma => host,
    )
    grown = prepare_plant_partitions(
        grown_assignment;
        run_seed=71, command_authority_target=host, command_endpoints)
    @test typeof(grown) === typeof(cpu)
    @test typeof(only(prepared_partitions(grown))) === typeof(cpu_partition)
    @test typeof(partition_controllable_optic_ids(
        only(prepared_partitions(grown)))) ===
        typeof(partition_controllable_optic_ids(cpu_partition))
    @test typeof(partition_command_endpoint_ids(
        only(prepared_partitions(grown)))) ===
        typeof(partition_command_endpoint_ids(cpu_partition))
    target_partition_test_error(() -> prepare_plant_partitions(
        stale_definition, cpu_assignment;
        run_seed=71, command_authority_target=host),
        :partition_assignment, :stale_assignment)
    target_partition_test_error(() -> prepare_plant_partitions(
        foreign_definition, cpu_assignment;
        run_seed=71, command_authority_target=host),
        :partition_assignment, :foreign_assignment)
    second_accelerator = AcceleratorComputeDevice(
        TargetPartitionTestBackend(), UInt32(2))
    target_partition_test_error(() -> prepare_plant_partitions(
        definition,
        mixed_assignment;
        run_seed=71,
        command_authority_target=second_accelerator,
        command_endpoints,
    ), :command_authority, :multiple_accelerators)
    target_partition_test_error(() -> prepare_plant_partitions(
        definition,
        cpu_assignment;
        run_seed=71,
        command_authority_target=:host,
        command_endpoints,
    ), :command_authority, :invalid_target)
    device_initial = (
        CommandEndpointConfiguration(
            :partition_test_command,
            TargetPartitionTestArray([0.0]);
            capacity=2,
        ),)
    target_partition_test_error(() -> prepare_plant_partitions(
        definition,
        cpu_assignment;
        run_seed=71,
        command_authority_target=host,
        command_endpoints=device_initial,
    ), :command_replica, :configuration_residency)
    safe_definition = target_partition_test_definition(
        command_silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(5),
        ))
    safe_assignment = resolve_plant_partition_assignment(
        safe_definition, host, :alpha => host, :beta => host)
    device_safe = (
        CommandEndpointConfiguration(
            :partition_test_command,
            [0.0];
            capacity=2,
            safe_command=TargetPartitionTestArray([0.0]),
        ),)
    target_partition_test_error(() -> prepare_plant_partitions(
        safe_definition,
        safe_assignment;
        run_seed=71,
        command_authority_target=host,
        command_endpoints=device_safe,
    ), :command_replica, :configuration_residency)

    authority_only_accelerator = prepare_plant_partitions(
        definition,
        cpu_assignment;
        run_seed=71,
        command_authority_target=accelerator,
        command_endpoints,
    )
    try
        TARGET_PARTITION_TEST_MISPLACE_ALLOCATIONS[] = true
        @test_throws ComputeDeviceError CommandAuthorityState(
            prepared_command_authority(authority_only_accelerator))
        @test_throws ComputeDeviceError prepare_plant_partitions(
            definition,
            cpu_assignment;
            run_seed=71,
            command_authority_target=accelerator,
            command_endpoints,
        )
    finally
        TARGET_PARTITION_TEST_MISPLACE_ALLOCATIONS[] = false
    end
    try
        TARGET_PARTITION_TEST_AVAILABLE[] = false
        target_partition_test_error(() -> prepare_plant_partitions(
            mixed_assignment;
            run_seed=71, command_authority_target=host),
            :partition_assignment, :unavailable_target)
        target_partition_test_error(() -> prepare_plant_partitions(
            definition,
            cpu_assignment;
            run_seed=71,
            command_authority_target=accelerator,
            command_endpoints,
        ), :command_authority, :unavailable_target)
    finally
        TARGET_PARTITION_TEST_AVAILABLE[] = true
    end
    target_partition_test_error(() -> prepared_partition(mixed,
        AcceleratorComputeDevice(TargetPartitionTestBackend(), UInt32(2))),
        :partition, :unknown_target)
end
