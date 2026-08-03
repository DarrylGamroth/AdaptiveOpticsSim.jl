@inline acquisition_conformance_snapshot(
    state::Plant.GlobalShutterAcquisitionState,
) = (
    state.status,
    state.sequence,
    state.exposure_start,
    state.exposure_close,
    state.integrated_through,
    state.readout_complete,
    state.readiness,
    state.next_read_index,
)

@inline acquisition_conformance_snapshot(
    state::Plant.RollingShutterAcquisitionState,
) = (
    state.status,
    state.sequence,
    state.frame_start,
    state.integrated_through,
    state.readout_complete,
    state.readiness,
    state.opened_bands,
    state.closed_bands,
)

@inline acquisition_conformance_snapshot(
    state::Plant.FrameTransferAcquisitionState,
) = (
    state.image_status,
    state.storage_status,
    state.sequence,
    state.image_sequence,
    state.storage_sequence,
    state.product_sequence,
    state.exposure_start,
    state.exposure_close,
    state.integrated_through,
    state.transfer_complete,
    state.storage_readout_complete,
    state.product_ready,
)

@inline acquisition_conformance_snapshot(
    state::Plant.DirectMeasurementAcquisitionState,
) = (
    state.status,
    state.sequence,
    state.exposure_start,
    state.exposure_close,
    state.integrated_through,
    state.readout_complete,
    state.readiness,
)

function assert_prepared_acquisition_lifecycle_conformance(
    prepared::P,
    state::S,
    foreign_state::S,
    target::D,
    expected_error::Type{E};
    detector_backed::Bool,
) where {
    P<:Plant.AbstractPreparedAcquisitionLifecycle,
    S<:Plant.AbstractAcquisitionLifecycleState,
    D<:Backends.AbstractComputeDevice,
    E<:Exception,
}
    @test (prepared isa Plant.AbstractPreparedDetectorAcquisitionLifecycle) ==
        detector_backed
    @test @inferred(Plant._require_exact_acquisition_lifecycle_target(
        prepared, target)) === prepared
    @test applicable(Plant.begin_exposure!, prepared, state,
        zero(Plant.PlantTimestamp))
    @test applicable(Plant.complete_readout!, prepared, state,
        zero(Plant.PlantTimestamp), Xoshiro(0))

    owner = Plant.StructuralResourceOwnerID(:acquisition, :pe01_interface)
    fact = Plant.structural_resource_fact(prepared, state, owner, target)
    @test fact isa Plant.AbstractStructuralResourceFact
    @test Backends.compute_device(fact) == target

    before = acquisition_conformance_snapshot(foreign_state)
    caught = try
        Plant.begin_exposure!(prepared, foreign_state,
            zero(Plant.PlantTimestamp))
        nothing
    catch error
        error
    end
    @test caught isa E
    if caught isa E
        @test caught.reason === :foreign_state
    end
    @test acquisition_conformance_snapshot(foreign_state) == before
    return fact
end
