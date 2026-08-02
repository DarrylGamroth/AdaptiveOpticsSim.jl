# Focused native pupil-OPD publication fixture. The fake accelerator permits
# scalar access only while its cold telescope is being constructed; every
# prepared publication, transfer, apply, and path-execution boundary runs with
# scalar indexing disabled.

const PATH_INPUT_TEST_SELECTED_DEVICE = Ref{UInt32}(0)

function Backends._with_compute_device(
    f::F,
    target::AcceleratorComputeDevice{HandoffTestBackend},
) where {F}
    previous_device = PATH_INPUT_TEST_SELECTED_DEVICE[]
    PATH_INPUT_TEST_SELECTED_DEVICE[] =
        UInt32(compute_device_identifier(target))
    try
        return f()
    finally
        PATH_INPUT_TEST_SELECTED_DEVICE[] = previous_device
    end
end

function Backends.allocate_array(
    ::HandoffTestBackend,
    ::Type{T},
    dimensions::Vararg{Int,N},
) where {T,N}
    ordinal = PATH_INPUT_TEST_SELECTED_DEVICE[]
    iszero(ordinal) && error(
        "fake path-input accelerator allocation requires exact device selection")
    target = AcceleratorComputeDevice(HandoffTestBackend(), ordinal)
    return HandoffTestArray(Array{T}(undef, dimensions...), target)
end

Backends.array_backend_type(::HandoffTestBackend) = HandoffTestArray

function Base.fill!(array::HandoffTestArray, value)
    fill!(array.storage, value)
    return array
end

function Base.copy(array::HandoffTestArray)
    return HandoffTestArray(copy(array.storage), array.device)
end

function Base.copyto!(
    destination::HandoffTestArray,
    source::HandoffTestArray,
)
    destination.device == source.device || throw(ArgumentError(
        "fake device-to-device copy requires one exact device"))
    copyto!(destination.storage, source.storage)
    return destination
end

function Base.copyto!(destination::HandoffTestArray, source::Array)
    HANDOFF_TEST_COLD_SCALAR_INDEXING[] || error(
        "unprepared host-to-fake-device copy is forbidden")
    copyto!(destination.storage, source)
    return destination
end

function Base.copyto!(destination::Array, source::HandoffTestArray)
    HANDOFF_TEST_COLD_SCALAR_INDEXING[] || error(
        "unprepared fake-device-to-host copy is forbidden")
    copyto!(destination, source.storage)
    return destination
end

function Base.copyto!(
    destination::HandoffTestArray,
    broadcasted::Base.Broadcast.Broadcasted,
)
    HANDOFF_TEST_COLD_SCALAR_INDEXING[] || error(
        "fake path-input accelerator broadcast is cold-preparation only")
    copyto!(destination.storage, broadcasted)
    return destination
end

function Plant.structural_array_bytes(
    array::HandoffTestArray,
    target::AbstractComputeDevice,
)
    compute_device(array) == target || throw(StructuralResourceError(
        :array_storage,
        :wrong_device,
        "path-input test array occupies $(compute_device(array)); expected $target",
    ))
    return UInt64(length(array)) * UInt64(sizeof(eltype(array)))
end

struct PathInputPublicationTestPathModel end
struct PathInputPublicationUnsupportedTestPathModel end
struct PathInputPublicationWrongRendererSourceTestPathModel end
struct PathInputPublicationWrongTelescopeTestPathModel end

Plant.plant_model_definition_style(
    ::Type{PathInputPublicationTestPathModel},
) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{PathInputPublicationUnsupportedTestPathModel},
) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{PathInputPublicationWrongRendererSourceTestPathModel},
) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{PathInputPublicationWrongTelescopeTestPathModel},
) = ColdPlantModelDefinition()

function Plant.prepare_pupil_opd_publication_materialization(
    ::PathInputPublicationTestPathModel,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    return prepare_pupil_opd_materialization(
        atmosphere, telescope, source, pupil)
end

function Plant.prepare_pupil_opd_publication_materialization(
    ::PathInputPublicationWrongRendererSourceTestPathModel,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    return path_input_publication_wrong_renderer_source_materialization(
        atmosphere, telescope, source, pupil)
end

function path_input_publication_wrong_renderer_source_materialization(
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    wrong_source = Source(
        band=:custom,
        wavelength=800e-9,
        photon_irradiance=1.0,
        coordinates=(18.0, 75.0),
    )
    renderer = prepare_atmosphere_renderer(
        atmosphere, telescope, wrong_source; T=eltype(pupil.opd))
    return PreparedPupilOPDMaterialization(
        Plant._PREPARED_PUPIL_OPD_MATERIALIZATION_TOKEN,
        renderer,
        pupil,
        source,
        telescope,
    )
end

function Plant.prepare_pupil_opd_publication_materialization(
    ::PathInputPublicationWrongTelescopeTestPathModel,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    return path_input_publication_wrong_telescope_materialization(
        atmosphere, telescope, source, pupil)
end

function path_input_publication_wrong_telescope_materialization(
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    other = Telescope(
        resolution=telescope.params.resolution,
        diameter=telescope.params.diameter,
        central_obstruction=0.25,
        fov_arcsec=telescope.params.fov_arcsec,
        T=eltype(pupil.opd),
        backend=backend(telescope),
    )
    AdaptiveOpticsSim.Optics.advance_aperture_revision!(other)
    return prepare_pupil_opd_materialization(
        atmosphere, other, source, pupil)
end

struct PathInputPublicationTestExecution{I,R}
    input::I
    result::R
end

function path_input_publication_test_result(pupil::PupilFunction)
    values = similar(pupil.opd)
    fill!(values, zero(eltype(values)))
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        values;
        coordinate_domain=MetricCoordinates(),
        sampling=pupil.metadata.sampling,
        origin=pupil.metadata.origin,
        orientation=pupil.metadata.orientation,
        spectral=MonochromaticChannel(eltype(values)(800e-9)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    return IntensityMap(metadata, values)
end

function path_input_publication_test_path_parts(telescope::Telescope)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    result = path_input_publication_test_result(pupil)
    return pupil, result, PathInputPublicationTestExecution(pupil, result)
end

function with_path_input_publication_cold_scalar_indexing(f::F) where {F}
    previous = HANDOFF_TEST_COLD_SCALAR_INDEXING[]
    HANDOFF_TEST_COLD_SCALAR_INDEXING[] = true
    try
        return f()
    finally
        HANDOFF_TEST_COLD_SCALAR_INDEXING[] = previous
    end
end

function Plant.validate_path_execution_binding(
    execution::PathInputPublicationTestExecution,
    input::PupilFunction,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(
            :path,
            :prepared_binding,
            "path-input publication test execution lost its exact products",
        ))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::PathInputPublicationTestExecution,
    target::AbstractComputeDevice,
)
    compute_device(execution.input.opd) == target &&
        compute_device(execution.result.values) == target || throw(
        PlantPreparationError(
            :path,
            :wrong_device,
            "path-input publication test execution occupies another target",
        ))
    return execution
end

function Plant.execute_path!(
    result::IntensityMap,
    input::PupilFunction,
    execution::PathInputPublicationTestExecution,
)
    validate_path_execution_binding(execution, input, result)
    copyto!(result.values, input.opd)
    return result
end

function Plant.structural_resource_fact(
    execution::PathInputPublicationTestExecution,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    validate_path_execution_target(execution, target)
    resident = structural_array_bytes(execution.input.support, target) +
        structural_array_bytes(execution.input.amplitude, target) +
        structural_array_bytes(execution.input.opd, target) +
        structural_array_bytes(execution.result.values, target)
    return KnownStructuralResourceFact(id, target, resident, UInt64(0))
end

function Plant.prepare_path_executor(
    ::PathInputPublicationTestPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    context,
)
    pupil, result, execution =
        path_input_publication_test_path_parts(telescope)
    materialization = prepare_pupil_opd_materialization(
        atmosphere, telescope, source, pupil)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        result,
        execution;
        context,
        materialization,
        optical_model=:path_input_publication_test,
        propagation_model=:copy_pupil_opd,
        model_revisions=UInt(1),
    )
end

function path_input_publication_test_target_local_resources(
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    context,
    optical_model::Symbol,
)
    pupil, result, execution =
        with_path_input_publication_cold_scalar_indexing() do
            path_input_publication_test_path_parts(telescope)
        end
    return PreparedTargetLocalPathResources(
        definition,
        source,
        telescope,
        pupil,
        result,
        execution;
        context,
        optical_model,
        propagation_model=:copy_pupil_opd,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_target_local_path_resources(
    ::PathInputPublicationTestPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    context,
)
    return path_input_publication_test_target_local_resources(
        definition, source, telescope, context,
        :path_input_publication_test)
end

function Plant.prepare_target_local_path_resources(
    ::PathInputPublicationUnsupportedTestPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    context,
)
    return path_input_publication_test_target_local_resources(
        definition, source, telescope, context,
        :unsupported_path_input_publication_test)
end

function Plant.prepare_target_local_path_resources(
    ::Union{
        PathInputPublicationWrongRendererSourceTestPathModel,
        PathInputPublicationWrongTelescopeTestPathModel,
    },
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    context,
)
    return path_input_publication_test_target_local_resources(
        definition, source, telescope, context,
        :invalid_path_input_publication_test)
end

function path_input_publication_test_definition(;
    alpha_model=PathInputPublicationTestPathModel(),
    beta_model=PathInputPublicationTestPathModel(),
)
    telescope = TelescopeDefinition(
        resolution=8,
        diameter=4.0,
        central_obstruction=0.0,
        revision=1,
    )
    atmosphere = MultiLayerAtmosphereDefinition(
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.65, 0.35],
        wind_speed=[7.0, 11.0],
        wind_direction=[20.0, 125.0],
        altitude=[0.0, 5_000.0],
        layer_ids=(:ground, :high),
    )
    alpha = OpticalPathDefinition(
        :alpha,
        Source(
            band=:custom,
            wavelength=800e-9,
            photon_irradiance=60.0,
            coordinates=(0.0, 0.0),
        ),
        alpha_model,
    )
    beta = OpticalPathDefinition(
        :beta,
        Source(
            band=:custom,
            wavelength=800e-9,
            photon_irradiance=45.0,
            coordinates=(4.0, 35.0),
        ),
        beta_model,
    )
    return PlantDefinition(; telescope, atmosphere, paths=(alpha, beta))
end

function path_input_publication_test_partitions(;
    alpha_target::AbstractComputeDevice=HostComputeDevice(),
    beta_target::AbstractComputeDevice=HANDOFF_TEST_ACCELERATOR,
    authority_target::AbstractComputeDevice=HostComputeDevice(),
    run_seed::Integer=0x218,
    alpha_model=PathInputPublicationTestPathModel(),
    beta_model=PathInputPublicationTestPathModel(),
)
    definition = path_input_publication_test_definition(;
        alpha_model, beta_model)
    assignment = resolve_plant_partition_assignment(
        definition,
        authority_target,
        :alpha => alpha_target,
        :beta => beta_target,
    )
    return with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed,
            command_authority_target=authority_target,
        )
    end
end

function path_input_publication_test_path(
    partitions::PreparedPlantPartitions,
    id::OpticalPathID,
)
    target = partition_target(getfield(partitions, :assignment), id)
    partition = prepared_partition(partitions, target)
    for path in prepared_paths(partition)
        path_id(path.definition) == id && return path
    end
    error("missing prepared path-input publication test path $id")
end

@inline path_input_publication_test_path(
    partitions::PreparedPlantPartitions,
    name::Symbol,
) = path_input_publication_test_path(partitions, OpticalPathID(name))

function path_input_publication_test_epoch!(
    partitions::PreparedPlantPartitions,
    timestamp::PlantTimestamp,
)
    authority = prepared_atmosphere_authority(partitions)
    atmosphere = prepared_atmosphere(authority)
    rng = Plant._prepared_atmosphere_rng(
        atmosphere, getfield(authority, :rngs))
    return advance_to!(
        atmosphere,
        plant_time_seconds(timestamp, Float64),
        rng,
    )
end

function path_input_publication_test_oracle(
    partitions::PreparedPlantPartitions,
    id::OpticalPathID,
)
    authority = prepared_atmosphere_authority(partitions)
    definition = path_definition(plant_definition(partitions), id)
    return prepare_path_executor(
        definition,
        prepared_telescope(authority),
        prepared_atmosphere(authority),
        getfield(authority, :context),
    )
end

@inline path_input_publication_test_oracle(
    partitions::PreparedPlantPartitions,
    name::Symbol,
) = path_input_publication_test_oracle(partitions, OpticalPathID(name))

function execute_path_input_publication_test_path!(
    path::PreparedTargetLocalPathResources,
)
    return Plant._with_completed_prepared_device_execution_context(
        getfield(path, :context)) do
        execute_path!(
            path_result(path),
            path_input(path),
            getfield(path, :execution),
        )
    end
end
