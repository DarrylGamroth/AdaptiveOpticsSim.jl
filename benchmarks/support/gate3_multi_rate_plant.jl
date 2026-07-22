module Gate3MultiRatePlantBenchmark

using AdaptiveOpticsSim
using SHA

include(joinpath(@__DIR__, "gate2_serial_plant.jl"))

const AOS = AdaptiveOpticsSim
const Gate2 = Gate2SerialPlantBenchmark

abstract type MultiRateSensorKind end
struct GlobalCMOSKind <: MultiRateSensorKind end
struct RollingCMOSKind <: MultiRateSensorKind end
struct CCDKind <: MultiRateSensorKind end
struct FrameTransferEMCCDKind <: MultiRateSensorKind end
struct HgCdTeRampKind <: MultiRateSensorKind end

struct MultiRateAcquisitionModel{
    T<:AbstractFloat,K<:MultiRateSensorKind}
    exposure_s::T
    quantum_efficiency::T
    kind::K
    rolling_line_s::T
    rolling_row_group_size::Int
    transfer_s::T
    nondestructive_reads::Int
end

AOS.plant_model_definition_style(::Type{<:MultiRateAcquisitionModel}) =
    AOS.ColdPlantModelDefinition()

@inline function multi_rate_sensor(model::MultiRateAcquisitionModel{
    T,<:GlobalCMOSKind}) where {T}
    return AOS.CMOSSensor(timing_model=AOS.GlobalShutter(), T=T)
end

@inline function multi_rate_sensor(model::MultiRateAcquisitionModel{
    T,<:RollingCMOSKind}) where {T}
    return AOS.CMOSSensor(timing_model=AOS.RollingShutter(
        model.rolling_line_s;
        row_group_size=model.rolling_row_group_size), T=T)
end

@inline function multi_rate_sensor(model::MultiRateAcquisitionModel{
    T,<:CCDKind}) where {T}
    return AOS.CCDSensor(T=T)
end

@inline function multi_rate_sensor(model::MultiRateAcquisitionModel{
    T,<:FrameTransferEMCCDKind}) where {T}
    return AOS.EMCCDSensor(acquisition_mode=AOS.FrameTransferAcquisition(
        transfer_time=model.transfer_s, T=T), T=T)
end

@inline function multi_rate_sensor(model::MultiRateAcquisitionModel{
    T,<:HgCdTeRampKind}) where {T}
    return AOS.HgCdTeAvalancheArraySensor(
        sampling_mode=AOS.UpTheRampSampling(model.nondestructive_reads),
        read_time=zero(T), T=T)
end

function AOS.prepare_acquisition_provider(
    model::MultiRateAcquisitionModel,
    ::AOS.AcquisitionDefinition,
    path::AOS.PreparedPathExecutor,
)
    AOS.require_path_result(path)
    T = eltype(AOS._first_path_result(path.result).values)
    detector = AOS.Detector(
        integration_time=T(model.exposure_s),
        noise=AOS.NoiseNone(),
        qe=T(model.quantum_efficiency),
        gain=one(T),
        response_model=AOS.NullFrameResponse(),
        sensor=multi_rate_sensor(model),
        T=T,
        backend=path.key.backend,
    )
    execution = AOS.FrameAcquisitionExecution(detector, path.result)
    metadata = (
        kind=:gate3_multi_rate_frame,
        units=:detected_electrons,
        geometry=path.result.metadata,
        detector=AOS.detector_export_metadata(detector),
        semantics=:complete_acquisition,
    )
    products = AOS.AcquisitionProducts(execution.observation; metadata)
    return AOS.prepare_full_optical_provider(execution, products)
end

mutable struct MultiRatePlantOperation{P,S,W}
    prepared::P
    state::S
    workspace::W
    processed_timestamps::UInt64
end

@inline function (operation::MultiRatePlantOperation)()
    timestamp = AOS.step_plant_events!(operation.prepared, operation.state,
        operation.workspace)
    timestamp === nothing && error("Gate 3 multi-rate plant exhausted events")
    operation.processed_timestamps += UInt64(1)
    return timestamp
end

function run_timestamp_window!(operation::MultiRatePlantOperation,
    timestamps::Int)
    timestamps >= 0 || error("timestamp window must be nonnegative")
    last_timestamp = zero(AOS.PlantTimestamp)
    @inbounds for _ in 1:timestamps
        last_timestamp = operation()
    end
    return last_timestamp
end

@inline ns_seconds(value::Integer) = Float64(value) / 1.0e9

function multi_rate_plant_definition(raw::AbstractDict;
    reverse_declarations::Bool=false)
    T = Float64
    telescope = AOS.Telescope(
        resolution=Int(raw["resolution"]),
        diameter=T(raw["diameter_m"]),
        central_obstruction=T(raw["central_obstruction"]),
        T=T,
    )
    atmosphere = AOS.MultiLayerAtmosphere(telescope;
        r0=T(raw["r0_m"]),
        L0=T(raw["outer_scale_m"]),
        fractional_cn2=T.(raw["fractional_cn2"]),
        wind_speed=T.(raw["wind_speed_m_per_s"]),
        wind_direction=T.(raw["wind_direction_deg"]),
        altitude=T.(raw["layer_altitude_m"]),
        layer_ids=Tuple(Symbol.(raw["layer_ids"])),
        T=T,
    )
    science_source = AOS.Source(
        band=:custom,
        wavelength=T(raw["science_wavelength_m"]),
        photon_irradiance=T(raw["science_photon_irradiance"]),
        coordinates=(T(raw["science_radius_arcsec"]),
            T(raw["science_azimuth_deg"])),
        T=T,
    )
    ngs_source = AOS.Source(
        band=:custom,
        wavelength=T(raw["ngs_wavelength_m"]),
        photon_irradiance=T(raw["ngs_photon_irradiance"]),
        coordinates=(T(raw["ngs_radius_arcsec"]),
            T(raw["ngs_azimuth_deg"])),
        T=T,
    )
    lgs_source = AOS.LGSSource(
        wavelength=T(raw["lgs_wavelength_m"]),
        photon_irradiance=T(raw["lgs_photon_irradiance"]),
        coordinates=(T(raw["lgs_radius_arcsec"]),
            T(raw["lgs_azimuth_deg"])),
        altitude=T(raw["lgs_altitude_m"]),
        T=T,
    )
    paths = (
        AOS.OpticalPathDefinition(:science, science_source,
            Gate2.DirectSciencePathModel(
                Int(raw["science_zero_padding"]), UInt(1))),
        AOS.OpticalPathDefinition(:ngs_shack_hartmann, ngs_source,
            Gate2.ShackHartmannPathModel(Int(raw["sh_n_lenslets"]),
                Int(raw["sh_n_pix_subap"]), UInt(2))),
        AOS.OpticalPathDefinition(:lgs_pyramid, lgs_source,
            Gate2.PyramidPathModel(Int(raw["pyramid_pupil_samples"]),
                T(raw["pyramid_modulation_lambda_over_d"]),
                Int(raw["pyramid_modulation_points"]), UInt(3))),
    )
    qe = T(raw["quantum_efficiency"])
    zero_s = zero(T)
    acquisitions = (
        AOS.AcquisitionDefinition(:science_cmos, :science,
            MultiRateAcquisitionModel(
                ns_seconds(raw["science_cmos_exposure_ns"]), qe,
                GlobalCMOSKind(), zero_s, 1, zero_s, 1)),
        AOS.AcquisitionDefinition(:science_ccd, :science,
            MultiRateAcquisitionModel(
                ns_seconds(raw["science_ccd_exposure_ns"]), qe,
                CCDKind(), zero_s, 1, zero_s, 1)),
        AOS.AcquisitionDefinition(:science_rolling, :science,
            MultiRateAcquisitionModel(
                ns_seconds(raw["science_rolling_exposure_ns"]), qe,
                RollingCMOSKind(),
                ns_seconds(raw["science_rolling_line_ns"]),
                Int(raw["science_rolling_row_group_size"]), zero_s, 1)),
        AOS.AcquisitionDefinition(:ngs_saphira, :ngs_shack_hartmann,
            MultiRateAcquisitionModel(
                ns_seconds(raw["ngs_saphira_exposure_ns"]), qe,
                HgCdTeRampKind(), zero_s, 1, zero_s,
                Int(raw["ngs_saphira_reads"]))),
        AOS.AcquisitionDefinition(:lgs_emccd, :lgs_pyramid,
            MultiRateAcquisitionModel(
                ns_seconds(raw["lgs_emccd_exposure_ns"]), qe,
                FrameTransferEMCCDKind(), zero_s, 1,
                ns_seconds(raw["lgs_emccd_transfer_ns"]), 1)),
    )
    ordered_paths = reverse_declarations ? reverse(paths) : paths
    ordered_acquisitions = reverse_declarations ? reverse(acquisitions) :
        acquisitions
    return AOS.PlantDefinition(; telescope, atmosphere,
        paths=ordered_paths, acquisitions=ordered_acquisitions)
end

@inline scaled_period(raw::AbstractDict, name::AbstractString,
    multiplier::Int) = Int64(raw[name]) * multiplier

function trigger_topology(raw::AbstractDict; faulted::Bool=false,
    period_multiplier::Int=1)
    source_faults = if faulted
        AOS.TriggerFaultTrace(
            AOS.TriggerFaultTraceEntry(raw["fault_drop_sequence"],
                :common_trigger_drop; action=AOS.DropTriggerEdge),
            AOS.TriggerFaultTraceEntry(raw["fault_shift_sequence"],
                :common_trigger_phase_error;
                phase_step=AOS.PlantTimeOffset(
                    raw["fault_phase_step_ns"]),
                jitter=AOS.PlantTimeOffset(raw["fault_jitter_ns"]),
                timestamp_label_offset=AOS.PlantTimeOffset(
                    raw["fault_label_offset_ns"])),
            AOS.TriggerFaultTraceEntry(raw["fault_duplicate_sequence"],
                :common_trigger_duplicate;
                action=AOS.DuplicateTriggerEdge,
                duplicate_delay=AOS.PlantDuration(
                    raw["fault_duplicate_delay_ns"])),
        )
    else
        AOS.TriggerFaultTrace()
    end
    source = AOS.TriggerSourceDefinition(:camera_trigger,
        AOS.PeriodicSchedule(
            period_ns=scaled_period(raw, "trigger_period_ns",
                period_multiplier),
            phase_ns=Int64(raw["trigger_phase_ns"]));
        faults=source_faults)
    source_id = AOS.TriggerSourceID(:camera_trigger)
    science_link = AOS.TriggerLinkDefinition(:science_camera_link,
        source_id; propagation_delay=AOS.PlantDuration(
            raw["science_trigger_delay_ns"]))
    ngs_link = AOS.TriggerLinkDefinition(:ngs_camera_link, source_id;
        propagation_delay=AOS.PlantDuration(raw["ngs_trigger_delay_ns"]))
    consumers = (
        AOS.TriggerConsumerDefinition(:science_camera,
            AOS.TriggerLinkID(:science_camera_link)),
        AOS.TriggerConsumerDefinition(:ngs_camera,
            AOS.TriggerLinkID(:ngs_camera_link)),
    )
    return AOS.prepare_trigger_topology((source,),
        (science_link, ngs_link), consumers;
        in_flight_capacity=Int(raw["trigger_in_flight_capacity"]))
end

function multi_rate_event_definition(raw::AbstractDict;
    reverse_declarations::Bool=false, faulted::Bool=false,
    period_multiplier::Int=1)
    schedule(period, phase) = AOS.PeriodicSchedule(
        period_ns=Int64(period) * period_multiplier,
        phase_ns=Int64(phase))
    samples = (
        AOS.OpticalSampleDefinition(:science,
            schedule(raw["science_sample_period_ns"],
                raw["science_sample_phase_ns"])),
        AOS.OpticalSampleDefinition(:ngs_shack_hartmann,
            schedule(raw["ngs_sample_period_ns"],
                raw["ngs_sample_phase_ns"])),
        AOS.OpticalSampleDefinition(:lgs_pyramid,
            schedule(raw["lgs_sample_period_ns"],
                raw["lgs_sample_phase_ns"])),
    )
    events = (
        AOS.DetectorEventDefinition(:science_cmos,
            AOS.GlobalShutterAcquisitionDefinition(AOS.PlantDuration(
                raw["science_cmos_exposure_ns"]);
                readout_duration=AOS.PlantDuration(
                    raw["science_cmos_readout_ns"]),
                readiness_delay=AOS.PlantDuration(
                    raw["science_cmos_readiness_ns"])),
            AOS.TriggeredAcquisitionStart(:science_camera)),
        AOS.DetectorEventDefinition(:science_ccd,
            AOS.GlobalShutterAcquisitionDefinition(AOS.PlantDuration(
                raw["science_ccd_exposure_ns"])),
            AOS.PeriodicAcquisitionStart(schedule(
                raw["science_ccd_period_ns"],
                raw["science_ccd_phase_ns"]))),
        AOS.DetectorEventDefinition(:science_rolling,
            AOS.RollingShutterAcquisitionDefinition(AOS.PlantDuration(
                raw["science_rolling_exposure_ns"])),
            AOS.PeriodicAcquisitionStart(schedule(
                raw["science_rolling_period_ns"],
                raw["science_rolling_phase_ns"]))),
        AOS.DetectorEventDefinition(:ngs_saphira,
            AOS.GlobalShutterAcquisitionDefinition(AOS.PlantDuration(
                raw["ngs_saphira_exposure_ns"]);
                readout_duration=AOS.PlantDuration(
                    raw["ngs_saphira_readout_ns"])),
            AOS.TriggeredAcquisitionStart(:ngs_camera)),
        AOS.DetectorEventDefinition(:lgs_emccd,
            AOS.FrameTransferAcquisitionDefinition(AOS.PlantDuration(
                raw["lgs_emccd_exposure_ns"]);
                readout_duration=AOS.PlantDuration(
                    raw["lgs_emccd_readout_ns"])),
            AOS.PeriodicAcquisitionStart(schedule(
                raw["lgs_emccd_period_ns"],
                raw["lgs_emccd_phase_ns"]))),
    )
    ordered_samples = reverse_declarations ? reverse(samples) : samples
    ordered_events = reverse_declarations ? reverse(events) : events
    return AOS.PlantEventLoopDefinition(ordered_samples, ordered_events;
        trigger_topology=trigger_topology(raw; faulted,
            period_multiplier))
end

function prepare_multi_rate_operation(raw::AbstractDict;
    reverse_declarations::Bool=false, faulted::Bool=false,
    period_multiplier::Int=1)
    definition = multi_rate_plant_definition(raw; reverse_declarations)
    plant = AOS.prepare_plant(definition;
        run_seed=UInt64(raw["run_seed"]),
        rng_derivation_version=Int(raw["rng_derivation_version"]))
    event_definition = multi_rate_event_definition(raw;
        reverse_declarations, faulted, period_multiplier)
    prepared = AOS.prepare_plant_event_loop(plant, event_definition)
    state = AOS.PlantEventLoopState(prepared)
    workspace = AOS.PlantEventLoopWorkspace(prepared)
    return plant, MultiRatePlantOperation(prepared, state, workspace,
        UInt64(0))
end

const ACQUISITION_IDS = (
    :science_cmos,
    :science_ccd,
    :science_rolling,
    :ngs_saphira,
    :lgs_emccd,
)

function observation_array(plant, id::Symbol)
    observation = AOS.acquisition_observation(
        AOS.prepared_acquisition(plant, id))
    values = observation isa AOS.WFSObservation ? observation.storage :
        observation
    return Array(values)
end

function observation_sha256(plant, id::Symbol)
    values = observation_array(plant, id)
    return bytes2hex(SHA.sha256(reinterpret(UInt8, vec(values))))
end

function product_snapshot(plant, operation::MultiRatePlantOperation)
    return map(ACQUISITION_IDS) do id
        (
            id=id,
            sequence=AOS.acquisition_product_sequence(operation.prepared,
                operation.state, id),
            ready_timestamp=AOS.acquisition_product_ready_timestamp(
                operation.prepared, operation.state, id),
            observation_sha256=observation_sha256(plant, id),
        )
    end
end

function product_sequence_vector(operation::MultiRatePlantOperation)
    return UInt64[AOS.acquisition_product_sequence(operation.prepared,
        operation.state, id) for id in ACQUISITION_IDS]
end

function storage_signature(operation::MultiRatePlantOperation)
    prepared = operation.prepared
    state = operation.state
    workspace = operation.workspace
    return (
        Base.summarysize(prepared.scheduler.definitions),
        Base.summarysize(prepared.actions),
        Base.summarysize(prepared.paths),
        Base.summarysize(prepared.acquisitions),
        Base.summarysize(prepared.trigger_topology),
        Base.summarysize(state.scheduler.cursors),
        Base.summarysize(state.acquisitions),
        Base.summarysize(state.path_sampled),
        Base.summarysize(state.product_sequences),
        Base.summarysize(state.product_ready_timestamps),
        Base.summarysize(state.trigger),
        Base.summarysize(workspace.scheduler.due_slots),
        Base.summarysize(workspace.due_paths),
        Base.summarysize(workspace.trigger),
        Base.summarysize(workspace.delivery),
    )
end

@inline timestamp_ns(timestamp::AOS.PlantTimestamp) =
    AOS.plant_nanoseconds(timestamp)
@inline timestamp_ns(::Nothing) = nothing

function snapshot_dict(snapshot)
    result = Dict{String,Any}()
    for entry in snapshot
        result[String(entry.id)] = Dict{String,Any}(
            "sequence" => Int64(entry.sequence),
            "product_ready_timestamp_ns" => timestamp_ns(
                entry.ready_timestamp),
            "observation_sha256" => entry.observation_sha256,
        )
    end
    return result
end

function validate_replay_and_reordering(raw::AbstractDict)
    canonical_plant, canonical = prepare_multi_rate_operation(raw)
    reordered_plant, reordered = prepare_multi_rate_operation(raw;
        reverse_declarations=true)
    horizon = AOS.PlantTimestamp(raw["replay_horizon_ns"])
    canonical_count = AOS.run_plant_events_until!(canonical.prepared,
        canonical.state, canonical.workspace, horizon)
    reordered_count = AOS.run_plant_events_until!(reordered.prepared,
        reordered.state, reordered.workspace, horizon)
    canonical_snapshot = product_snapshot(canonical_plant, canonical)
    reordered_snapshot = product_snapshot(reordered_plant, reordered)
    canonical_count == reordered_count || error(
        "Gate 3 timestamp count changed after declaration reordering")
    canonical_snapshot == reordered_snapshot || error(
        "Gate 3 products changed after declaration reordering")

    expected_sequences = (UInt64(3), UInt64(2), UInt64(3), UInt64(3),
        UInt64(3))
    Tuple(entry.sequence for entry in canonical_snapshot) ==
        expected_sequences || error(
        "Gate 3 representative product sequence changed")
    expected_ready_ns = (
        1_490_000_000,
        950_000_000,
        1_480_000_000,
        1_460_000_000,
        1_245_000_000,
    )
    Tuple(timestamp_ns(entry.ready_timestamp) for entry in
        canonical_snapshot) == expected_ready_ns || error(
        "Gate 3 representative product readiness changed")

    science_owners = map((:science_cmos, :science_ccd,
        :science_rolling)) do id
        AOS.prepared_acquisition(canonical_plant, id)
    end
    all(owner -> owner.path_result === first(science_owners).path_result,
        science_owners) || error(
        "Gate 3 science acquisitions do not reuse one path result")
    all(entry -> begin
            values = observation_array(canonical_plant, entry.id)
            all(isfinite, values) && sum(values) > 0
        end, canonical_snapshot) || error(
        "Gate 3 representative plant produced an invalid product")

    epoch = AOS.current_epoch(canonical.prepared.atmosphere)
    reordered_epoch = AOS.current_epoch(reordered.prepared.atmosphere)
    AOS.epoch_time(epoch) == AOS.epoch_time(reordered_epoch) || error(
        "Gate 3 replay atmosphere times differ")
    AOS.epoch_sequence(epoch) == AOS.epoch_sequence(reordered_epoch) ||
        error("Gate 3 replay atmosphere sequences differ")

    return Dict{String,Any}(
        "declaration_reordering_replay" => true,
        "processed_timestamp_count" => canonical_count,
        "path_count" => AOS.plant_event_path_count(canonical.prepared),
        "acquisition_count" => AOS.plant_event_acquisition_count(
            canonical.prepared),
        "generator_count" => AOS.plant_event_generator_count(
            canonical.prepared),
        "science_path_reused" => true,
        "atmosphere_model_time_s" => AOS.epoch_time(epoch),
        "atmosphere_epoch_sequence" => Int64(AOS.epoch_sequence(epoch)),
        "products" => snapshot_dict(canonical_snapshot),
    )
end

function validate_faulted_trigger_fanout(raw::AbstractDict)
    plant, operation = prepare_multi_rate_operation(raw; faulted=true)
    horizon = AOS.PlantTimestamp(raw["fault_horizon_ns"])
    timestamp_count = AOS.run_plant_events_until!(operation.prepared,
        operation.state, operation.workspace, horizon)
    snapshot = product_snapshot(plant, operation)
    sequence_by_id = Dict(entry.id => entry.sequence for entry in snapshot)
    sequence_by_id[:science_cmos] == UInt64(5) || error(
        "faulted common trigger did not produce five science CMOS products")
    sequence_by_id[:ngs_saphira] == UInt64(5) || error(
        "faulted common trigger did not produce five NGS products")
    ready_by_id = Dict(entry.id => timestamp_ns(entry.ready_timestamp)
        for entry in snapshot)
    ready_by_id[:science_cmos] == 2_710_000_000 || error(
        "faulted science CMOS readiness changed")
    ready_by_id[:ngs_saphira] == 2_680_000_000 || error(
        "faulted NGS readiness changed")
    return Dict{String,Any}(
        "integrated_fault_trace" => true,
        "common_source_fanout" => true,
        "drop_sequence" => Int(raw["fault_drop_sequence"]),
        "phase_step_and_jitter_sequence" =>
            Int(raw["fault_shift_sequence"]),
        "duplicate_sequence" => Int(raw["fault_duplicate_sequence"]),
        "processed_timestamp_count" => timestamp_count,
        "products" => snapshot_dict(snapshot),
    )
end

function validate_fixed_storage(raw::AbstractDict, horizon_ns::Integer,
    max_timestamps::Integer)
    _, operation = prepare_multi_rate_operation(raw)
    AOS.run_plant_events_until!(operation.prepared, operation.state,
        operation.workspace, AOS.PlantTimestamp(2_000_000_000))
    before = storage_signature(operation)
    timestamp_count = AOS.run_plant_events_until!(operation.prepared,
        operation.state, operation.workspace,
        AOS.PlantTimestamp(horizon_ns))
    after = storage_signature(operation)
    before == after || error(
        "Gate 3 prepared storage changed with simulated run length")
    timestamp_count > 0 || error("Gate 3 long run processed no timestamps")
    timestamp_count <= max_timestamps || error(
        "Gate 3 long run exceeded its logical-timestamp bound")
    return Dict{String,Any}(
        "fixed_storage" => true,
        "storage_component_bytes" => collect(before),
        "additional_processed_timestamps" => timestamp_count,
        "maximum_processed_timestamps" => Int(max_timestamps),
        "horizon_ns" => Int64(horizon_ns),
    )
end

function validate_direct_jump(raw::AbstractDict, multiplier::Int)
    multiplier > 1 || error("direct-jump multiplier must exceed one")
    _, operation = prepare_multi_rate_operation(raw;
        period_multiplier=multiplier)
    processed = AOS.run_plant_events_until!(operation.prepared,
        operation.state, operation.workspace,
        AOS.PlantTimestamp(1_000_000_000))
    next_timestamp = AOS.next_plant_event_timestamp(operation.prepared,
        operation.state, operation.workspace)
    expected = AOS.PlantTimestamp(scaled_period(raw,
        "science_sample_period_ns", multiplier))
    next_timestamp == expected || error(
        "Gate 3 long-period schedule did not jump to the next event")
    return Dict{String,Any}(
        "direct_logical_timestamp_jump" => true,
        "period_multiplier" => multiplier,
        "initial_window_timestamp_count" => processed,
        "next_timestamp_ns" => timestamp_ns(next_timestamp),
        "empty_base_ticks_visited" => 0,
    )
end

end # module
