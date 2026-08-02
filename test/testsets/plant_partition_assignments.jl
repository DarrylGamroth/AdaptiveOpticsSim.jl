struct PartitionAssignmentTestPathModel
    label::Symbol
end

struct PartitionAssignmentTestAcquisitionModel
    label::Symbol
end

struct PartitionAssignmentFakeBackend <: AbstractArrayBackend end

struct PartitionAssignmentDuplicateLayerAtmosphereDefinition <:
    AbstractTimedAtmosphereDefinition end

Plant.plant_model_definition_style(::Type{PartitionAssignmentTestPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{PartitionAssignmentTestAcquisitionModel}) = ColdPlantModelDefinition()

Backends.compute_device_availability(
    ::AcceleratorComputeDevice{PartitionAssignmentFakeBackend}) =
    ComputeDeviceAvailable()

Plant.partition_atmosphere_layer_ids(
    ::PartitionAssignmentDuplicateLayerAtmosphereDefinition) =
    (AtmosphereLayerID(:duplicate), AtmosphereLayerID(:duplicate))

for name in (
    :ResolvedPlantPartitionAssignment,
    :assigned_path_ids,
    :atmosphere_authority_target,
    :partition_target,
    :partition_targets,
    :resolve_plant_partition_assignment,
)
    isdefined(@__MODULE__, name) && continue
    @eval const $(name) = getfield(Plant, $(QuoteNode(name)))
end

function partition_assignment_test_definition(; reverse_paths::Bool=false,
    path_suffix::Symbol=:base, include_gamma::Bool=false,
    include_paths::Bool=true,
    atmosphere_definition=KolmogorovAtmosphereDefinition(r0=0.2, L0=25.0))
    telescope = TelescopeDefinition(
        resolution=8,
        diameter=8.0,
        central_obstruction=0.0,
        revision=1,
    )
    source = Source(band=:I, magnitude=0.0)
    alpha = OpticalPathDefinition(:alpha, source,
        PartitionAssignmentTestPathModel(path_suffix))
    beta = OpticalPathDefinition(:beta, source,
        PartitionAssignmentTestPathModel(path_suffix))
    alpha_camera = AcquisitionDefinition(:alpha_camera, :alpha,
        PartitionAssignmentTestAcquisitionModel(:alpha))
    beta_camera = AcquisitionDefinition(:beta_camera, :beta,
        PartitionAssignmentTestAcquisitionModel(:beta))
    gamma = OpticalPathDefinition(:gamma, source,
        PartitionAssignmentTestPathModel(path_suffix))
    gamma_camera = AcquisitionDefinition(:gamma_camera, :gamma,
        PartitionAssignmentTestAcquisitionModel(:gamma))
    base_paths = include_paths ?
        (include_gamma ? (alpha, beta, gamma) : (alpha, beta)) : ()
    base_acquisitions = include_paths ?
        (include_gamma ? (alpha_camera, beta_camera, gamma_camera) :
         (alpha_camera, beta_camera)) : ()
    paths = reverse_paths ? reverse(base_paths) : base_paths
    acquisitions = reverse_paths ? reverse(base_acquisitions) :
        base_acquisitions
    return PlantDefinition(;
        telescope,
        atmosphere=atmosphere_definition,
        paths,
        acquisitions,
    )
end

function partition_assignment_test_error(f, reason::Symbol)
    error = try
        f()
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === :partition_assignment
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

@testset "Caller-resolved plant partition assignments" begin
    definition = partition_assignment_test_definition()
    host = HostComputeDevice()
    assignment = resolve_plant_partition_assignment(
        definition, host, :beta => host, OpticalPathID(:alpha) => host)

    @test assignment isa ResolvedPlantPartitionAssignment
    @test plant_definition(assignment) === definition
    @test atmosphere_authority_target(assignment) == host
    @test Tuple(partition_targets(assignment)) == (host,)
    @test assigned_path_ids(assignment, host) ==
        (OpticalPathID(:alpha), OpticalPathID(:beta))
    @test partition_target(assignment, :alpha) == host
    @test partition_target(assignment, OpticalPathID(:beta)) == host
    @test Plant._require_current_partition_assignment_definition(
        assignment, definition) === assignment
    @test_throws ArgumentError ResolvedPlantPartitionAssignment(
        Plant._ResolvedPlantPartitionAssignmentToken(),
        getfield(assignment, :definition),
        getfield(assignment, :definition_identity),
        getfield(assignment, :topology),
        getfield(assignment, :targets),
        getfield(assignment, :atmosphere_authority_target_ordinal),
        getfield(assignment, :paths),
    )

    reordered_definition = partition_assignment_test_definition(
        reverse_paths=true)
    reordered = resolve_plant_partition_assignment(
        reordered_definition, host, :alpha => host, :beta => host)
    @test assigned_path_ids(reordered, host) == assigned_path_ids(assignment,
        host)
    @test partition_targets(reordered) == partition_targets(assignment)

    authority_only_definition = partition_assignment_test_definition(
        include_paths=false)
    authority_only = resolve_plant_partition_assignment(
        authority_only_definition, host)
    @test Tuple(partition_targets(authority_only)) == (host,)
    @test assigned_path_ids(authority_only, host) == ()
    @test atmosphere_authority_target(authority_only) == host

    layered_definition = partition_assignment_test_definition(
        atmosphere_definition=MultiLayerAtmosphereDefinition(
            r0=0.2,
            L0=25.0,
            fractional_cn2=(0.4, 0.6),
            wind_speed=(5.0, 10.0),
            wind_direction=(0.0, 90.0),
            altitude=(5000.0, 0.0),
            layer_ids=(:high, :ground),
        ),
    )
    layered = resolve_plant_partition_assignment(
        layered_definition, host, :alpha => host, :beta => host)
    @test Tuple(getfield(layered, :topology).atmosphere_layers) ==
        (AtmosphereLayerID(:ground), AtmosphereLayerID(:high))
    duplicate_layer_definition = partition_assignment_test_definition(
        atmosphere_definition=
            PartitionAssignmentDuplicateLayerAtmosphereDefinition(),
    )
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(duplicate_layer_definition, host,
            :alpha => host, :beta => host),
        :duplicate_atmosphere_layer_id)

    fake_accelerator = AcceleratorComputeDevice(
        PartitionAssignmentFakeBackend(), UInt32(1))
    mixed = resolve_plant_partition_assignment(
        definition, host, :beta => fake_accelerator, :alpha => host)
    @test Tuple(partition_targets(mixed)) == (host, fake_accelerator)
    @test assigned_path_ids(mixed, host) == (OpticalPathID(:alpha),)
    @test assigned_path_ids(mixed, fake_accelerator) == (OpticalPathID(:beta),)
    @test partition_target(mixed, :beta) == fake_accelerator

    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, host,
            :alpha => host), :missing_path)
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, host,
            :alpha => host, :alpha => host, :beta => host), :duplicate_path)
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, host,
            :alpha => host, :missing => host), :unknown_path)
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, host,
            :alpha => host, :beta => :not_a_target), :invalid_target)
    partition_assignment_test_error(() -> partition_target(assignment,
        :missing), :unknown_path)
    partition_assignment_test_error(() -> assigned_path_ids(assignment,
        AcceleratorComputeDevice(CUDABackend(), 0)), :unknown_target)
    unavailable = AcceleratorComputeDevice(CUDABackend(), 0)
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, unavailable,
            :alpha => unavailable, :beta => unavailable), :unavailable_target)
    second_accelerator = AcceleratorComputeDevice(
        PartitionAssignmentFakeBackend(), UInt32(2))
    partition_assignment_test_error(() ->
        resolve_plant_partition_assignment(definition, host,
            :alpha => fake_accelerator, :beta => second_accelerator),
        :multiple_accelerators)

    same_topology_new_definition = partition_assignment_test_definition()
    partition_assignment_test_error(() ->
        Plant._require_current_partition_assignment_definition(assignment,
            same_topology_new_definition), :stale_assignment)
    foreign_definition = partition_assignment_test_definition(
        include_gamma=true)
    grown_assignment = @inferred resolve_plant_partition_assignment(
        foreign_definition,
        host,
        :alpha => host,
        :beta => host,
        :gamma => host,
    )
    @test typeof(grown_assignment) === typeof(assignment)
    @test typeof(getfield(grown_assignment, :topology)) ===
        typeof(getfield(assignment, :topology))
    @test typeof(getfield(grown_assignment, :paths)) ===
        typeof(getfield(assignment, :paths))
    partition_assignment_test_error(() ->
        Plant._require_current_partition_assignment_definition(assignment,
            foreign_definition), :foreign_assignment)
end
