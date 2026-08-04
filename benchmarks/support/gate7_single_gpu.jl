module Gate7SingleGPUBenchmark

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Plant
using AdaptiveOpticsSim.WavefrontSensors

const AOS = AdaptiveOpticsSim
const AOSPlant = AdaptiveOpticsSim.Plant

struct Gate7ShackHartmannPathModel
    n_lenslets::Int
    n_pix_subap::Int
end

struct Gate7WitnessPathModel
    zero_padding::Int
end

struct Gate7WitnessAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

AOSPlant.plant_model_definition_style(
    ::Type{Gate7ShackHartmannPathModel},
) = AOSPlant.ColdPlantModelDefinition()

AOSPlant.plant_model_definition_style(
    ::Type{Gate7WitnessPathModel},
) = AOSPlant.ColdPlantModelDefinition()

AOSPlant.plant_model_definition_style(
    ::Type{<:Gate7WitnessAcquisitionModel},
) = AOSPlant.ColdPlantModelDefinition()

function AOSPlant.prepare_path_executor(
    model::Gate7ShackHartmannPathModel,
    definition::OpticalPathDefinition,
    source::AOS.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AOS.Atmospheres.AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(
        telescope; T, backend=AOS.Backends.backend(telescope))
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=model.n_lenslets,
        n_pix_subap=model.n_pix_subap,
        mode=Diffractive(),
        T=T,
        backend=AOS.Backends.backend(telescope),
    )
    front_end = ShackHartmannOpticalFrontEnd(sensor.front_end, source)
    output = shack_hartmann_rate_map(front_end, pupil)
    plan = prepare_wfs_optics(front_end, pupil, output)
    execution = AOSPlant.WFSOpticalPathExecution(plan)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        output,
        execution;
        context=context,
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere,
            telescope,
            source,
            pupil,
        ),
        optical_model=(
            kind=:gate7_benchmark_shack_hartmann,
            n_lenslets=model.n_lenslets,
            n_pix_subap=model.n_pix_subap,
        ),
        propagation_model=:diffractive_shack_hartmann,
        model_revisions=UInt(1),
    )
end

function AOSPlant.prepare_path_executor(
    model::Gate7WitnessPathModel,
    definition::OpticalPathDefinition,
    source::AOS.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AOS.Atmospheres.AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(
        telescope; T, backend=AOS.Backends.backend(telescope))
    execution = prepare_direct_imaging(
        pupil,
        source;
        zero_padding=model.zero_padding,
    )
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(execution),
        execution;
        context=context,
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere,
            telescope,
            source,
            pupil,
        ),
        optical_model=(
            kind=:gate7_benchmark_witness,
            zero_padding=model.zero_padding,
        ),
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function AOSPlant.prepare_acquisition_provider(
    model::Gate7WitnessAcquisitionModel,
    ::AcquisitionDefinition,
    path::AOSPlant.PreparedPathExecutor,
)
    AOSPlant.require_path_result(path)
    result = AOSPlant.path_result(path)
    T = eltype(result.values)
    detector = Detector(
        integration_time=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=AOSPlant.path_result_key(path).backend,
    )
    execution = AOSPlant.FrameAcquisitionExecution(detector, result)
    products = AOSPlant.AcquisitionProducts(
        execution.observation;
        metadata=(
            kind=:gate7_benchmark_witness,
            units=:detected_electrons,
            geometry=result.metadata,
            detector=detector_export_metadata(detector),
        ),
    )
    return prepare_full_optical_provider(execution, products)
end

@inline function _prepare_gate7_event_loop(
    plant,
    definition,
    ::Val{:independent},
)
    return AOSPlant._prepare_plant_event_loop(
        plant,
        definition,
        Val(:none),
    )
end

@inline function _prepare_gate7_event_loop(
    plant,
    definition,
    ::Val{:batched},
)
    return prepare_plant_event_loop(plant, definition)
end

mutable struct Gate7Operation{P,E,S,W}
    plant::P
    prepared::E
    state::S
    workspace::W
    completed_steps::Int
end

@inline gate7_path_ids() = (:wfs_alpha, :wfs_beta)

function prepare_gate7_operation(
    selector::AOS.Backends.AbstractArrayBackend,
    workload::AbstractDict,
    mode::Val,
)
    T = Float32
    resolution = Int(workload["pupil_resolution"])
    telescope = AOS.Optics.TelescopeDefinition(
        resolution=resolution,
        diameter=T(workload["telescope_diameter_m"]),
        central_obstruction=zero(T),
        revision=UInt(1),
        T=T,
    )
    atmosphere = AOS.Atmospheres.MultiLayerAtmosphereDefinition(;
        r0=T(workload["r0_m"]),
        L0=T(workload["outer_scale_m"]),
        fractional_cn2=T.(workload["fractional_cn2"]),
        wind_speed=T.(workload["wind_speed_m_s"]),
        wind_direction=T.(workload["wind_direction_deg"]),
        altitude=T.(workload["layer_altitude_m"]),
        layer_ids=(:ground, :high),
        T=T,
    )
    wavelength_m = T(workload["wavelength_m"])
    source_alpha = Source(
        band=:custom,
        wavelength=wavelength_m,
        photon_irradiance=T(workload["alpha_photon_irradiance"]),
        coordinates=(
            T(workload["alpha_offset_arcsec"]),
            T(workload["alpha_azimuth_deg"]),
        ),
        T=T,
    )
    source_beta = Source(
        band=:custom,
        wavelength=wavelength_m,
        photon_irradiance=T(workload["beta_photon_irradiance"]),
        coordinates=(
            T(workload["beta_offset_arcsec"]),
            T(workload["beta_azimuth_deg"]),
        ),
        T=T,
    )
    witness_source = Source(
        band=:custom,
        wavelength=wavelength_m,
        photon_irradiance=one(T),
        coordinates=(zero(T), zero(T)),
        T=T,
    )
    wfs_model = Gate7ShackHartmannPathModel(
        Int(workload["lenslets"]),
        Int(workload["samples_per_subaperture"]),
    )
    paths = (
        OpticalPathDefinition(:wfs_alpha, source_alpha, wfs_model),
        OpticalPathDefinition(:wfs_beta, source_beta, wfs_model),
        OpticalPathDefinition(
            :declaration_witness,
            witness_source,
            Gate7WitnessPathModel(1),
        ),
    )
    acquisitions = (
        AcquisitionDefinition(
            :declaration_witness_camera,
            :declaration_witness,
            Gate7WitnessAcquisitionModel(T(0.0005)),
        ),
    )
    definition = PlantDefinition(
        ;
        telescope,
        atmosphere,
        paths,
        acquisitions,
    )
    target = AOS.Backends.compute_device(
        AOS.Backends.allocate_array(selector, UInt8, 1),
    )
    plant = prepare_plant(
        definition,
        target;
        run_seed=UInt64(workload["run_seed"]),
    )
    sample_period_ns = Int64(workload["sample_period_ns"])
    witness_phase_ns = Int64(workload["witness_phase_ns"])
    witness_period_ns = Int64(workload["witness_period_ns"])
    samples = (
        OpticalSampleDefinition(
            :wfs_alpha,
            PeriodicSchedule(period_ns=sample_period_ns),
        ),
        OpticalSampleDefinition(
            :wfs_beta,
            PeriodicSchedule(period_ns=sample_period_ns),
        ),
        OpticalSampleDefinition(
            :declaration_witness,
            PeriodicSchedule(
                period_ns=witness_period_ns,
                phase_ns=witness_phase_ns,
            ),
        ),
    )
    acquisition_events = (
        AcquisitionEventDefinition(
            :declaration_witness_camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(500_000);
                readout_duration=PlantDuration(100_000),
                readiness_delay=PlantDuration(100_000),
            ),
            PeriodicAcquisitionStart(
                PeriodicSchedule(
                    period_ns=witness_period_ns,
                    phase_ns=witness_phase_ns + 100_000,
                ),
            ),
        ),
    )
    event_definition = PlantEventLoopDefinition(samples, acquisition_events)
    prepared = _prepare_gate7_event_loop(plant, event_definition, mode)
    return Gate7Operation(
        plant,
        prepared,
        AOSPlant.PlantEventLoopState(prepared),
        AOSPlant.PlantEventLoopWorkspace(prepared),
        0,
    )
end

@inline function _gate7_path(operation::Gate7Operation, id::Symbol)
    return AOSPlant.prepared_path(operation.plant, id)
end

function synchronize_gate7_outputs!(operation::Gate7Operation)
    @inbounds for id in gate7_path_ids()
        values = _gate7_path(operation, id).result.values
        AOS.Backends.synchronize_backend!(
            AOS.Backends.execution_style(values))
    end
    return operation
end

function step_gate7_device_ready!(operation::Gate7Operation)
    timestamp = step_plant_events!(
        operation.prepared,
        operation.state,
        operation.workspace,
    )
    timestamp === nothing &&
        error("Gate 7 benchmark exhausted its periodic event stream")
    synchronize_gate7_outputs!(operation)
    operation.completed_steps += 1
    return timestamp
end

function gate7_host_buffers(operation::Gate7Operation)
    return map(gate7_path_ids()) do id
        values = _gate7_path(operation, id).result.values
        Matrix{eltype(values)}(undef, size(values))
    end
end

function copy_gate7_outputs!(
    host_buffers::Tuple,
    operation::Gate7Operation,
)
    @inbounds for index in eachindex(host_buffers)
        values =
            _gate7_path(operation, gate7_path_ids()[index]).result.values
        copyto!(host_buffers[index], values)
    end
    return host_buffers
end

struct Gate7DeviceReadyBoundary{O}
    operation::O
end

struct Gate7HostReadyBoundary{O,H}
    operation::O
    host_buffers::H
end

struct Gate7TransferOnlyBoundary{O,H}
    operation::O
    host_buffers::H
end

@inline function run_gate7_boundary!(boundary::Gate7DeviceReadyBoundary)
    return step_gate7_device_ready!(boundary.operation)
end

@inline function run_gate7_boundary!(boundary::Gate7HostReadyBoundary)
    timestamp = step_gate7_device_ready!(boundary.operation)
    copy_gate7_outputs!(boundary.host_buffers, boundary.operation)
    return timestamp
end

@inline function run_gate7_boundary!(boundary::Gate7TransferOnlyBoundary)
    copy_gate7_outputs!(boundary.host_buffers, boundary.operation)
    return nothing
end

function prepare_gate7_boundary(
    id::AbstractString,
    selector::AOS.Backends.AbstractArrayBackend,
    workload::AbstractDict,
)
    if id == "independent_device_ready"
        operation =
            prepare_gate7_operation(selector, workload, Val(:independent))
        return Gate7DeviceReadyBoundary(operation)
    elseif id == "batched_device_ready"
        operation =
            prepare_gate7_operation(selector, workload, Val(:batched))
        return Gate7DeviceReadyBoundary(operation)
    elseif id == "batched_host_ready"
        operation =
            prepare_gate7_operation(selector, workload, Val(:batched))
        return Gate7HostReadyBoundary(
            operation,
            gate7_host_buffers(operation),
        )
    elseif id == "transfer_only"
        operation =
            prepare_gate7_operation(selector, workload, Val(:batched))
        step_gate7_device_ready!(operation)
        return Gate7TransferOnlyBoundary(
            operation,
            gate7_host_buffers(operation),
        )
    end
    error("unknown Gate 7 benchmark boundary '$id'")
end

function gate7_residency(operation::Gate7Operation)
    path_records = Dict{String,Any}()
    expected_device = nothing
    @inbounds for id in gate7_path_ids()
        path = _gate7_path(operation, id)
        input_device = AOS.Backends.compute_device(path.input.opd)
        output_device = AOS.Backends.compute_device(path.result.values)
        input_device == output_device ||
            error("Gate 7 benchmark path $id spans compute devices")
        if expected_device === nothing
            expected_device = input_device
        else
            input_device == expected_device ||
                error("Gate 7 benchmark paths occupy different devices")
        end
        path_records[String(id)] = Dict{String,Any}(
            "input_type" => string(typeof(path.input.opd)),
            "output_type" => string(typeof(path.result.values)),
            "input_device" => string(input_device),
            "output_device" => string(output_device),
            "output_shape" => collect(size(path.result.values)),
        )
    end
    layer_records = Vector{Dict{String,Any}}(
        undef,
        length(operation.prepared.atmosphere.layers),
    )
    @inbounds for index in eachindex(operation.prepared.atmosphere.layers)
        screen = operation.prepared.atmosphere.layers[
            index].generator.state.opd
        AOS.Backends.compute_device(screen) == expected_device ||
            error("Gate 7 atmosphere layer $index occupies another device")
        layer_records[index] = Dict{String,Any}(
            "type" => string(typeof(screen)),
            "device" => string(AOS.Backends.compute_device(screen)),
            "shape" => collect(size(screen)),
        )
    end
    owner_count = AOSPlant.device_path_batch_owner_count(operation.prepared)
    owner_record = if iszero(owner_count)
        Dict{String,Any}(
            "count" => 0,
            "group_count" => 0,
            "device" => "not applicable",
            "implementation_type" => "not applicable",
        )
    else
        owner_count == 1 ||
            error("Gate 7 benchmark expected exactly one device owner")
        owner = AOSPlant.device_path_batch_owner(operation.prepared, 1)
        AOSPlant.device_path_batch_compute_device(owner) == expected_device ||
            error("Gate 7 benchmark owner occupies another device")
        AOSPlant.device_path_batch_group_count(owner) == 2 ||
            error("Gate 7 benchmark owner must retain exactly two WFS paths")
        Dict{String,Any}(
            "count" => 1,
            "group_count" => AOSPlant.device_path_batch_group_count(owner),
            "device" => string(
                AOSPlant.device_path_batch_compute_device(owner),
            ),
            "implementation_type" => string(typeof(owner.implementation)),
            "context_type" => string(typeof(owner.implementation.context)),
        )
    end
    return Dict{String,Any}(
        "compute_device" => string(expected_device),
        "paths" => path_records,
        "atmosphere_layers" => layer_records,
        "owner" => owner_record,
    )
end

function gate7_rng_states(operation::Gate7Operation)
    return map(operation.prepared.atmosphere_rng.streams) do stream
        copy(stream.state)
    end
end

function gate7_submission_proxy(
    operation::Gate7Operation,
    ::Val{:independent},
)
    owner_count =
        AOSPlant.device_path_batch_owner_count(operation.prepared)
    iszero(owner_count) ||
        error("Gate 7 independent proxy unexpectedly has a device owner")
    path_count = length(gate7_path_ids())
    return Dict{String,Any}(
        "top_level_path_submissions" => path_count,
        "device_owner_submissions" => 0,
        "atmosphere_direction_render_calls" => path_count,
        "wfs_optics_calls" => path_count,
    )
end

function gate7_submission_proxy(
    operation::Gate7Operation,
    ::Val{:batched},
)
    owner_count =
        AOSPlant.device_path_batch_owner_count(operation.prepared)
    owner_count == 1 ||
        error("Gate 7 batched proxy requires exactly one device owner")
    owner = AOSPlant.device_path_batch_owner(operation.prepared, 1)
    group_count = AOSPlant.device_path_batch_group_count(owner)
    group_count == length(gate7_path_ids()) ||
        error("Gate 7 batched proxy owner has an unexpected group count")
    return Dict{String,Any}(
        "top_level_path_submissions" => 0,
        "device_owner_submissions" => owner_count,
        "atmosphere_direction_render_calls" => owner_count,
        "wfs_optics_calls" => group_count,
    )
end

function gate7_submission_proxy_evidence(
    independent::Gate7Operation,
    batched::Gate7Operation,
    expected::AbstractDict,
)
    observed = Dict{String,Any}(
        "independent" =>
            gate7_submission_proxy(independent, Val(:independent)),
        "batched" => gate7_submission_proxy(batched, Val(:batched)),
    )
    for mode in ("independent", "batched")
        observed[mode] == expected[mode] || error(
            "Gate 7 $mode submission proxy differs from its contract",
        )
    end
    independent_proxy = observed["independent"]
    batched_proxy = observed["batched"]
    batched_proxy["atmosphere_direction_render_calls"] <
        independent_proxy["atmosphere_direction_render_calls"] ||
        error("Gate 7 batch did not reduce atmosphere render calls")
    batched_proxy["wfs_optics_calls"] ==
        independent_proxy["wfs_optics_calls"] ||
        error("Gate 7 batch changed the modeled WFS optics count")
    return Dict{String,Any}(
        "passed" => true,
        "derivation" =>
            "prepared Plant path/owner topology; not a vendor launch trace",
        "observed" => observed,
    )
end

function gate7_correctness_evidence(
    independent::Gate7Operation,
    batched::Gate7Operation;
    rtol::Real,
    atol::Real,
)
    independent.completed_steps == batched.completed_steps ||
        error("Gate 7 parity contexts completed different step counts")
    independent.state.scheduler.cursors ==
        batched.state.scheduler.cursors ||
        error("Gate 7 parity scheduler cursors differ")
    independent.state.scheduler.revision ==
        batched.state.scheduler.revision ||
        error("Gate 7 parity scheduler revisions differ")
    independent.state.path_sampled == batched.state.path_sampled ||
        error("Gate 7 parity sampled-path state differs")
    independent_epoch = current_epoch(independent.prepared.atmosphere)
    batched_epoch = current_epoch(batched.prepared.atmosphere)
    epoch_time(independent_epoch) == epoch_time(batched_epoch) ||
        error("Gate 7 parity atmosphere times differ")
    epoch_sequence(independent_epoch) == epoch_sequence(batched_epoch) ||
        error("Gate 7 parity atmosphere sequences differ")
    gate7_rng_states(independent) == gate7_rng_states(batched) ||
        error("Gate 7 parity atmosphere RNG states differ")
    AOSPlant.device_path_batch_owner_count(independent.prepared) == 0 ||
        error("Gate 7 independent baseline unexpectedly has a device owner")
    AOSPlant.device_path_batch_owner_count(batched.prepared) == 1 ||
        error("Gate 7 batched context has no device owner")

    path_records = Dict{String,Any}()
    @inbounds for id in gate7_path_ids()
        independent_path = _gate7_path(independent, id)
        batched_path = _gate7_path(batched, id)
        independent_opd = Array(independent_path.input.opd)
        batched_opd = Array(batched_path.input.opd)
        isapprox(batched_opd, independent_opd; rtol, atol) ||
            error("Gate 7 parity pupil OPD differs for $id")
        independent_output = Array(independent_path.result.values)
        batched_output = Array(batched_path.result.values)
        isapprox(batched_output, independent_output; rtol, atol) ||
            error("Gate 7 parity WFS product differs for $id")
        difference =
            Float64.(batched_output) .- Float64.(independent_output)
        reference = Float64.(independent_output)
        reference_norm = sqrt(sum(abs2, reference))
        path_records[String(id)] = Dict{String,Any}(
            "output_shape" => collect(size(reference)),
            "cpu_visible_checksum" => sum(reference),
            "batched_checksum" => sum(Float64, batched_output),
            "maximum_absolute_error" => maximum(abs, difference),
            "relative_l2_error" => iszero(reference_norm) ? 0.0 :
                sqrt(sum(abs2, difference)) / reference_norm,
        )
    end
    return Dict{String,Any}(
        "passed" => true,
        "rtol" => Float64(rtol),
        "atol" => Float64(atol),
        "completed_steps" => independent.completed_steps,
        "scheduler_revision" => independent.state.scheduler.revision,
        "atmosphere_epoch_sequence" => epoch_sequence(independent_epoch),
        "paths" => path_records,
    )
end

end # module
