function captured_pupil_opd_publication_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_pupil_opd_publication_error(f, reason::Symbol)
    error = captured_pupil_opd_publication_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === :pupil_opd_publication
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function assert_path_materialization_binding_error(f)
    error = captured_pupil_opd_publication_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === :path
        @test error.reason === :prepared_binding
        @test !isempty(error.msg)
    end
    return error
end

function forged_pupil_opd_publication(
    publication::MaterializedPupilOPDPublication;
    route_identity=getfield(publication, :route_identity),
    atmosphere_identity=getfield(publication, :atmosphere_identity),
    epoch=getfield(publication, :epoch),
    timestamp=getfield(publication, :timestamp),
    path=getfield(publication, :path),
    route_generation=getfield(publication, :route_generation),
    handoff_slot=getfield(publication, :handoff_slot),
    handoff_generation=getfield(publication, :handoff_generation),
)
    return MaterializedPupilOPDPublication(
        Plant._MATERIALIZED_PUPIL_OPD_PUBLICATION_TOKEN,
        route_identity,
        atmosphere_identity,
        epoch,
        timestamp,
        path,
        route_generation,
        handoff_slot,
        handoff_generation,
    )
end

function pupil_opd_publication_allocation_sample!(
    direct::Plant.PreparedDirectPupilOPDPublicationRoute,
    remote::Plant.PreparedRemotePupilOPDPublicationRoute,
    epoch::AtmosphereEpoch,
    timestamp::PlantTimestamp,
)
    direct_output = prepare_pupil_opd_publication_output(direct)
    remote_output = prepare_pupil_opd_publication_output(remote)
    direct_materialize = @allocated begin
        @assert materialize_pupil_opd_publication!(
            direct_output, direct, epoch, timestamp) ==
            PupilOPDPublicationSucceeded
    end
    remote_materialize = @allocated begin
        @assert materialize_pupil_opd_publication!(
            remote_output, remote, epoch, timestamp) ==
            PupilOPDPublicationSucceeded
    end
    direct_publication = direct_output[]
    remote_publication = remote_output[]
    submit = @allocated begin
        @assert submit_pupil_opd_publication!(
            remote, remote_publication) == PupilOPDPublicationSucceeded
    end
    complete = @allocated begin
        @assert try_complete_pupil_opd_publication!(
            remote, remote_publication) == PupilOPDPublicationSucceeded
    end
    apply = @allocated begin
        @assert apply_pupil_opd_publication!(
            remote, remote_publication) == PupilOPDPublicationSucceeded
    end
    reclaim_remote = @allocated begin
        @assert reclaim_pupil_opd_publication!(
            remote, remote_publication) == PupilOPDPublicationSucceeded
    end
    reclaim_direct = @allocated begin
        @assert reclaim_pupil_opd_publication!(
            direct, direct_publication) == PupilOPDPublicationSucceeded
    end
    return (
        direct_materialize,
        remote_materialize,
        submit,
        complete,
        apply,
        reclaim_remote,
        reclaim_direct,
    )
end

@testset "Authority-owned pupil-OPD publication routes" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions()
    direct = prepare_pupil_opd_publication_route(partitions, :alpha)
    remote = prepare_pupil_opd_publication_route(
        partitions, OpticalPathID(:beta))
    direct_path = path_input_publication_test_path(partitions, :alpha)
    remote_path = path_input_publication_test_path(partitions, :beta)

    @test direct isa Plant.PreparedDirectPupilOPDPublicationRoute
    @test remote isa Plant.PreparedRemotePupilOPDPublicationRoute
    @test isconcretetype(typeof(direct))
    @test isconcretetype(typeof(remote))
    @test !ismutabletype(typeof(direct))
    @test !ismutabletype(typeof(remote))
    @test ismutabletype(typeof(getfield(direct, :state)))
    @test ismutabletype(typeof(getfield(remote, :state)))
    @test pupil_opd_publication_handoff_capacity(direct) == 0
    @test pupil_opd_publication_handoff_capacity(remote) == 1
    @test pupil_opd_publication_authority_target(direct) ==
        HostComputeDevice()
    @test pupil_opd_publication_authority_target(remote) ==
        HostComputeDevice()
    @test pupil_opd_publication_route_target(direct) == HostComputeDevice()
    @test pupil_opd_publication_route_target(remote) ==
        HANDOFF_TEST_ACCELERATOR
    @test pupil_opd_publication_path_id(direct) == OpticalPathID(:alpha)
    @test pupil_opd_publication_path_id(remote) == OpticalPathID(:beta)
    @test pupil_opd_publication_route_identity(direct) !==
        pupil_opd_publication_route_identity(remote)
    @test fieldcount(
        typeof(pupil_opd_publication_route_identity(remote))) == 0
    @test_throws MethodError MaterializedPupilOPDPublication(
        nothing,
        nothing,
        nothing,
        PlantTimestamp(0),
        OpticalPathID(:beta),
        UInt64(1),
        UInt32(1),
        UInt64(1),
    )
    @test :PreparedPupilOPDPublicationRoute in names(Plant)
    @test :PreparedDirectPupilOPDPublicationRoute ∉ names(Plant)
    @test :PreparedRemotePupilOPDPublicationRoute ∉ names(Plant)
    @test :PreparedPupilOPDPublicationRoute ∉ names(AdaptiveOpticsSim)

    direct_output = @inferred prepare_pupil_opd_publication_output(direct)
    remote_output = @inferred prepare_pupil_opd_publication_output(remote)
    @test isconcretetype(eltype(typeof(direct_output)))
    @test isconcretetype(eltype(typeof(remote_output)))
    @test Base.allocatedinline(eltype(typeof(remote_output)))

    alpha_oracle = path_input_publication_test_oracle(partitions, :alpha)
    beta_oracle = path_input_publication_test_oracle(partitions, :beta)
    replay_before = rng_replay_metadata(partitions)
    timestamp = PlantTimestamp(0)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)
    materialize_path_input!(alpha_oracle, epoch)
    materialize_path_input!(beta_oracle, epoch)

    remote_before = copy(path_input(remote_path).opd.storage)
    remote_support_before = copy(path_input(remote_path).support.storage)
    remote_amplitude_before = copy(path_input(remote_path).amplitude.storage)
    @test @inferred(materialize_pupil_opd_publication!(
        direct_output, direct, epoch, timestamp)) ==
        PupilOPDPublicationSucceeded
    @test @inferred(materialize_pupil_opd_publication!(
        remote_output, remote, epoch, timestamp)) ==
        PupilOPDPublicationSucceeded
    direct_publication = direct_output[]
    remote_publication = remote_output[]

    @test pupil_opd_publication_route_identity(direct_publication) ===
        pupil_opd_publication_route_identity(direct)
    @test pupil_opd_publication_route_identity(remote_publication) ===
        pupil_opd_publication_route_identity(remote)
    @test pupil_opd_publication_atmosphere_identity(remote_publication) ===
        atmosphere_authority_identity(partitions)
    @test pupil_opd_publication_epoch(remote_publication) == epoch
    @test pupil_opd_publication_timestamp(remote_publication) == timestamp
    @test pupil_opd_publication_path_id(remote_publication) ==
        OpticalPathID(:beta)
    @test :source ∉ fieldnames(typeof(remote_publication))
    @test path_input(direct_path).opd == path_input(alpha_oracle).opd
    @test getfield(remote, :source_pupil).opd == path_input(beta_oracle).opd
    @test path_input(direct_path).opd != path_input(beta_oracle).opd
    @test path_input(remote_path).opd.storage == remote_before
    @test materialize_pupil_opd_publication!(
        remote_output, remote, epoch, timestamp) ==
        PupilOPDPublicationRouteBusy
    @test apply_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationNotCompleted
    @test submit_pupil_opd_publication!(remote, :wrong_product) ==
        PupilOPDPublicationRejected

    HANDOFF_TEST_PENDING_OBSERVATIONS[] = UInt8(1)
    @test @inferred(submit_pupil_opd_publication!(
        remote, remote_publication)) == PupilOPDPublicationSucceeded
    @test path_input(remote_path).opd.storage == remote_before
    @test @inferred(try_complete_pupil_opd_publication!(
        remote, remote_publication)) ==
        PupilOPDPublicationTransferPending
    @test apply_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationNotCompleted
    @test @inferred(try_complete_pupil_opd_publication!(
        remote, remote_publication)) == PupilOPDPublicationSucceeded
    @test path_input(remote_path).opd.storage == remote_before
    @test @inferred(apply_pupil_opd_publication!(
        remote, remote_publication)) == PupilOPDPublicationSucceeded
    @test apply_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationAlreadyApplied
    @test path_input(remote_path).opd.storage == path_input(beta_oracle).opd
    @test path_input(remote_path).support.storage == remote_support_before
    @test path_input(remote_path).amplitude.storage == remote_amplitude_before
    @test HANDOFF_TEST_COLD_SCALAR_INDEXING[] == false

    execute_path!(alpha_oracle)
    execute_path!(beta_oracle)
    execute_path_input_publication_test_path!(direct_path)
    execute_path_input_publication_test_path!(remote_path)
    @test path_result(direct_path).values == path_result(alpha_oracle).values
    @test path_result(remote_path).values.storage ==
        path_result(beta_oracle).values
    @test rng_replay_metadata(partitions) == replay_before

    @test @inferred(submit_pupil_opd_publication!(
        direct, direct_publication)) == PupilOPDPublicationSucceeded
    @test @inferred(try_complete_pupil_opd_publication!(
        direct, direct_publication)) == PupilOPDPublicationSucceeded
    @test @inferred(apply_pupil_opd_publication!(
        direct, direct_publication)) == PupilOPDPublicationSucceeded
    @test @inferred(reclaim_pupil_opd_publication!(
        remote, remote_publication)) == PupilOPDPublicationSucceeded
    @test @inferred(reclaim_pupil_opd_publication!(
        direct, direct_publication)) == PupilOPDPublicationSucceeded
end

@testset "Pupil-OPD publication validation and retained inputs" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions()
    direct = prepare_pupil_opd_publication_route(partitions, :alpha)
    remote = prepare_pupil_opd_publication_route(partitions, :beta)
    direct_output = prepare_pupil_opd_publication_output(direct)
    remote_output = prepare_pupil_opd_publication_output(remote)
    remote_path = path_input_publication_test_path(partitions, :beta)
    timestamp = PlantTimestamp(0)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)

    unchanged_direct = copy(path_input(
        path_input_publication_test_path(partitions, :alpha)).opd)
    assert_pupil_opd_publication_error(:timestamp) do
        materialize_pupil_opd_publication!(
            direct_output, direct, epoch, PlantTimestamp(1))
    end
    @test path_input(
        path_input_publication_test_path(partitions, :alpha)).opd ==
        unchanged_direct

    @test materialize_pupil_opd_publication!(
        direct_output, direct, epoch, timestamp) ==
        PupilOPDPublicationSucceeded
    @test materialize_pupil_opd_publication!(
        remote_output, remote, epoch, timestamp) ==
        PupilOPDPublicationSucceeded
    remote_publication = remote_output[]
    source_snapshot = copy(getfield(remote, :source_pupil).opd)

    other_identity = AtmosphereIdentity()
    other_epoch = AtmosphereEpoch(
        other_identity, epoch_time(epoch), epoch_sequence(epoch))
    wrong_publications = (
        forged_pupil_opd_publication(
            remote_publication;
            route_identity=pupil_opd_publication_route_identity(direct),
        ),
        forged_pupil_opd_publication(
            remote_publication;
            atmosphere_identity=other_identity,
            epoch=other_epoch,
        ),
        forged_pupil_opd_publication(
            remote_publication;
            epoch=AtmosphereEpoch(
                atmosphere_authority_identity(partitions),
                epoch_time(epoch),
                epoch_sequence(epoch) + one(UInt64),
            ),
        ),
        forged_pupil_opd_publication(
            remote_publication; timestamp=PlantTimestamp(1)),
        forged_pupil_opd_publication(
            remote_publication; path=OpticalPathID(:alpha)),
        forged_pupil_opd_publication(
            remote_publication;
            route_generation=getfield(
                remote_publication, :route_generation) + one(UInt64),
        ),
        forged_pupil_opd_publication(
            remote_publication; handoff_slot=UInt32(0)),
        forged_pupil_opd_publication(
            remote_publication;
            handoff_generation=getfield(
                remote_publication, :handoff_generation) + one(UInt64),
        ),
    )
    for wrong in wrong_publications
        @test submit_pupil_opd_publication!(remote, wrong) ==
            PupilOPDPublicationRejected
    end
    @test path_input(remote_path).opd.storage != source_snapshot

    # Deliberately advance before submission to prove that the route retains an
    # OPD product rather than treating an epoch token as retained storage.
    # The mixed-composition barrier added by #219 will forbid this ordering.
    next_timestamp = PlantTimestamp(1_000_000)
    next_epoch = path_input_publication_test_epoch!(
        partitions, next_timestamp)
    @test submit_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationSucceeded
    @test try_complete_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationSucceeded
    @test apply_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationSucceeded
    @test path_input(remote_path).opd.storage == source_snapshot
    @test reclaim_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationSucceeded
    @test reclaim_pupil_opd_publication!(direct, direct_output[]) ==
        PupilOPDPublicationSucceeded

    @test_throws AtmosphereEpochError begin
        materialize_pupil_opd_publication!(
            remote_output, remote, epoch, timestamp)
    end
    @test materialize_pupil_opd_publication!(
        remote_output, remote, next_epoch, next_timestamp) ==
        PupilOPDPublicationSucceeded
    @test submit_pupil_opd_publication!(remote, remote_publication) ==
        PupilOPDPublicationRejected
end

@testset "Pupil-OPD publication requires explicit model opt-in" begin
    partitions = path_input_publication_test_partitions(
        beta_model=PathInputPublicationUnsupportedTestPathModel())
    error = assert_pupil_opd_publication_error(:unsupported_model) do
        prepare_pupil_opd_publication_route(partitions, :beta)
    end
    @test occursin("does not opt in", error.msg)
    @test HANDOFF_TEST_COLD_SCALAR_INDEXING[] == false
end

@testset "Pupil-OPD publication rejects disguised optical bindings" begin
    direct_wrong_source_partitions = path_input_publication_test_partitions(
        alpha_model=PathInputPublicationWrongRendererSourceTestPathModel())
    source_error = assert_path_materialization_binding_error() do
        prepare_pupil_opd_publication_route(
            direct_wrong_source_partitions, :alpha)
    end
    @test occursin("source geometry", source_error.msg)

    remote_wrong_source_partitions = path_input_publication_test_partitions(
        beta_model=PathInputPublicationWrongRendererSourceTestPathModel())
    remote_source_error = assert_path_materialization_binding_error() do
        prepare_pupil_opd_publication_route(
            remote_wrong_source_partitions, :beta)
    end
    @test occursin("source geometry", remote_source_error.msg)

    direct_wrong_telescope_partitions = path_input_publication_test_partitions(
        alpha_model=PathInputPublicationWrongTelescopeTestPathModel())
    telescope_error = assert_path_materialization_binding_error() do
        prepare_pupil_opd_publication_route(
            direct_wrong_telescope_partitions, :alpha)
    end
    @test occursin("exact telescope", telescope_error.msg)

    remote_wrong_telescope_partitions = path_input_publication_test_partitions(
        beta_model=PathInputPublicationWrongTelescopeTestPathModel())
    remote_telescope_error = assert_path_materialization_binding_error() do
        prepare_pupil_opd_publication_route(
            remote_wrong_telescope_partitions, :beta)
    end
    @test occursin("exact telescope", remote_telescope_error.msg)
    @test HANDOFF_TEST_COLD_SCALAR_INDEXING[] == false
end

@testset "Prepared path executors reject disguised optical bindings" begin
    partitions = path_input_publication_test_partitions()
    path = path_input_publication_test_oracle(partitions, :alpha)
    wrong_source_materialization =
        path_input_publication_wrong_renderer_source_materialization(
            getfield(path, :atmosphere),
            getfield(path, :telescope),
            getfield(path, :source),
            path_input(path),
        )
    constructor_error = assert_path_materialization_binding_error() do
        PreparedPathExecutor(
            getfield(path, :definition),
            getfield(path, :source),
            getfield(path, :telescope),
            getfield(path, :atmosphere),
            path_input(path),
            path_result(path),
            getfield(path, :execution);
            context=getfield(path, :context),
            materialization=wrong_source_materialization,
            optical_model=:invalid_path_input_publication_test,
            propagation_model=:copy_pupil_opd,
            model_revisions=UInt(1),
        )
    end
    @test occursin("source geometry", constructor_error.msg)

    wrong_telescope_materialization =
        path_input_publication_wrong_telescope_materialization(
            getfield(path, :atmosphere),
            getfield(path, :telescope),
            getfield(path, :source),
            path_input(path),
        )
    telescope_error = assert_path_materialization_binding_error() do
        PreparedPathExecutor(
            getfield(path, :definition),
            getfield(path, :source),
            getfield(path, :telescope),
            getfield(path, :atmosphere),
            path_input(path),
            path_result(path),
            getfield(path, :execution);
            context=getfield(path, :context),
            materialization=wrong_telescope_materialization,
            optical_model=:invalid_path_input_publication_test,
            propagation_model=:copy_pupil_opd,
            model_revisions=UInt(1),
        )
    end
    @test occursin("exact telescope", telescope_error.msg)

    forged = PreparedPathExecutor(
        Plant._PREPARED_PATH_EXECUTOR_TOKEN,
        getfield(path, :definition),
        getfield(path, :source),
        getfield(path, :telescope),
        getfield(path, :atmosphere),
        getfield(path, :context),
        path_input(path),
        path_result(path),
        wrong_source_materialization,
        getfield(path, :execution),
        path_result_key(path),
    )
    timestamp = PlantTimestamp(0)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)
    before = copy(path_input(path).opd)
    current_error = assert_path_materialization_binding_error() do
        materialize_path_input!(forged, epoch)
    end
    @test occursin("source geometry", current_error.msg)
    @test path_input(path).opd == before
    @test HANDOFF_TEST_COLD_SCALAR_INDEXING[] == false
end

@testset "Remote pupil-OPD apply revalidates target bindings" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions(run_seed=0x222)
    remote = prepare_pupil_opd_publication_route(partitions, :beta)
    output = prepare_pupil_opd_publication_output(remote)
    timestamp = PlantTimestamp(0)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)
    @test materialize_pupil_opd_publication!(
        output, remote, epoch, timestamp) == PupilOPDPublicationSucceeded
    publication = output[]
    @test submit_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationSucceeded
    @test try_complete_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationSucceeded

    path = path_input_publication_test_path(partitions, :beta)
    before = copy(path_input(path).opd.storage)
    AdaptiveOpticsSim.Optics.advance_aperture_revision!(
        getfield(path, :telescope))
    assert_pupil_opd_publication_error(:revision) do
        apply_pupil_opd_publication!(remote, publication)
    end
    @test path_input(path).opd.storage == before
    @test apply_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationUncertain
    @test reclaim_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationUncertain
end

@testset "Pupil-OPD publication transfer failures" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions()
    remote = prepare_pupil_opd_publication_route(partitions, :beta)
    output = prepare_pupil_opd_publication_output(remote)
    timestamp = PlantTimestamp(0)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)
    @test materialize_pupil_opd_publication!(
        output, remote, epoch, timestamp) == PupilOPDPublicationSucceeded
    publication = output[]

    HANDOFF_TEST_FAIL_SUBMISSION[] = true
    @test submit_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationTransferFailed
    HANDOFF_TEST_FAIL_SUBMISSION[] = false
    @test try_complete_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationNotSubmitted
    @test reclaim_pupil_opd_publication!(remote, publication) ==
        PupilOPDPublicationSucceeded

    next_timestamp = PlantTimestamp(1_000_000)
    next_epoch = path_input_publication_test_epoch!(
        partitions, next_timestamp)
    @test materialize_pupil_opd_publication!(
        output, remote, next_epoch, next_timestamp) ==
        PupilOPDPublicationSucceeded
    next_publication = output[]
    @test submit_pupil_opd_publication!(remote, next_publication) ==
        PupilOPDPublicationSucceeded
    HANDOFF_TEST_FAIL_COMPLETION[] = true
    @test try_complete_pupil_opd_publication!(remote, next_publication) ==
        PupilOPDPublicationTransferFailed
    HANDOFF_TEST_FAIL_COMPLETION[] = false
    @test apply_pupil_opd_publication!(remote, next_publication) ==
        PupilOPDPublicationNotCompleted
    @test reclaim_pupil_opd_publication!(remote, next_publication) ==
        PupilOPDPublicationSucceeded

    uncertain_partitions = path_input_publication_test_partitions(
        run_seed=0x219)
    uncertain = prepare_pupil_opd_publication_route(
        uncertain_partitions, :beta)
    uncertain_output = prepare_pupil_opd_publication_output(uncertain)
    uncertain_epoch = path_input_publication_test_epoch!(
        uncertain_partitions, timestamp)
    @test materialize_pupil_opd_publication!(
        uncertain_output, uncertain, uncertain_epoch, timestamp) ==
        PupilOPDPublicationSucceeded
    HANDOFF_TEST_THROW_SUBMISSION[] = true
    @test_throws ErrorException submit_pupil_opd_publication!(
        uncertain, uncertain_output[])
    HANDOFF_TEST_THROW_SUBMISSION[] = false
    @test try_complete_pupil_opd_publication!(
        uncertain, uncertain_output[]) == PupilOPDPublicationUncertain
    @test apply_pupil_opd_publication!(
        uncertain, uncertain_output[]) == PupilOPDPublicationUncertain
    @test reclaim_pupil_opd_publication!(
        uncertain, uncertain_output[]) == PupilOPDPublicationUncertain
end

@testset "Pupil-OPD publication inference and allocation" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions(run_seed=0x220)
    direct = prepare_pupil_opd_publication_route(partitions, :alpha)
    remote = prepare_pupil_opd_publication_route(partitions, :beta)

    warm_timestamp = PlantTimestamp(0)
    warm_epoch = path_input_publication_test_epoch!(
        partitions, warm_timestamp)
    pupil_opd_publication_allocation_sample!(
        direct, remote, warm_epoch, warm_timestamp)

    timestamp = PlantTimestamp(1_000_000)
    epoch = path_input_publication_test_epoch!(partitions, timestamp)
    allocations = pupil_opd_publication_allocation_sample!(
        direct, remote, epoch, timestamp)
    if !coverage_instrumented()
        @test allocations == (0, 0, 0, 0, 0, 0, 0)
    end
end

@testset "Pupil-OPD publication routes are explicit opt-in resources" begin
    definition = PlantDefinition(;
        telescope=TelescopeDefinition(
            resolution=4,
            diameter=4.0,
            central_obstruction=0.0,
            revision=1,
        ),
        atmosphere=KolmogorovAtmosphereDefinition(r0=0.2, L0=25.0),
        paths=(OpticalPathDefinition(
            :reduced,
            Source(band=:I, magnitude=0.0),
            PathInputPublicationTestPathModel(),
        ),),
    )
    assignment = resolve_plant_partition_assignment(
        definition, HostComputeDevice(), :reduced => HostComputeDevice())
    partitions = prepare_plant_partitions(
        definition,
        assignment;
        run_seed=0x221,
        command_authority_target=HostComputeDevice(),
    )
    path = only(prepared_paths(only(prepared_partitions(partitions))))
    @test !hasproperty(path, :materialization)
    @test !hasproperty(path, :handoff)
    @test !hasproperty(partitions, :path_input_routes)
    atmosphere = prepared_atmosphere(
        prepared_atmosphere_authority(partitions))
    @test !getfield(atmosphere, :state).timeline.initialized
end
