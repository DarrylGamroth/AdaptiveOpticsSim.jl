module Gate4CommandPlantBenchmark

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant
using SHA

const AOS = AdaptiveOpticsSim
const AOSPlant = AdaptiveOpticsSim.Plant

struct CommandPlaneOpticModel{T<:AbstractFloat}
    axis::Symbol
    gain::T
end

struct CommandSciencePathModel end

struct CommandFrameAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

struct PreparedCommandPlaneOptic{T<:AbstractFloat,A<:AbstractMatrix{T}}
    endpoint::AOSPlant.CommandEndpointID
    pattern::A
end

mutable struct CommandPlaneOpticState{T<:AbstractFloat}
    visible::T
end

mutable struct CommandPlaneOpticWorkspace{T<:AbstractFloat}
    staged::T
end

AOSPlant.plant_model_definition_style(
    ::Type{<:CommandPlaneOpticModel}) = AOSPlant.ColdPlantModelDefinition()
AOSPlant.plant_model_definition_style(
    ::Type{CommandSciencePathModel}) = AOSPlant.ColdPlantModelDefinition()
AOSPlant.plant_model_definition_style(
    ::Type{<:CommandFrameAcquisitionModel}) =
    AOSPlant.ColdPlantModelDefinition()

function AOSPlant.prepare_controllable_optic(
    model::CommandPlaneOpticModel{T},
    definition::AOSPlant.ControllableOpticDefinition,
    telescope::AOS.Telescope,
    ::AOS.AbstractAtmosphere,
) where {T}
    endpoint = only(AOSPlant.command_endpoint_ids(definition))
    resolution = size(AOS.pupil_reflectivity(telescope), 1)
    pattern = Matrix{T}(undef, resolution, resolution)
    center = T(resolution + 1) / T(2)
    scale = model.gain / T(resolution)
    if model.axis === :x
        @inbounds for column in axes(pattern, 2), row in axes(pattern, 1)
            pattern[row, column] = scale * (T(column) - center)
        end
    elseif model.axis === :y
        @inbounds for column in axes(pattern, 2), row in axes(pattern, 1)
            pattern[row, column] = scale * (T(row) - center)
        end
    else
        throw(AOSPlant.PlantPreparationError(:controllable_optic,
            :invalid_axis, "command benchmark optic axis must be :x or :y"))
    end
    return PreparedCommandPlaneOptic(endpoint, pattern)
end

function AOSPlant.prepare_controllable_optic_state(
    prepared::PreparedCommandPlaneOptic{T},
    ::AOSPlant.ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
) where {T}
    only(endpoint_ids) == prepared.endpoint ||
        throw(AOSPlant.PlantPreparationError(:controllable_optic,
            :prepared_binding, "command benchmark endpoint binding changed"))
    return CommandPlaneOpticState(T(only(initial_commands)))
end

AOSPlant.prepare_controllable_optic_workspace(
    ::PreparedCommandPlaneOptic{T}) where {T} =
    CommandPlaneOpticWorkspace(zero(T))

function AOSPlant.stage_controllable_optic_command!(
    prepared::PreparedCommandPlaneOptic{T},
    ::CommandPlaneOpticState{T},
    workspace::CommandPlaneOpticWorkspace{T},
    endpoint::AOSPlant.CommandEndpointID,
    command,
    ::AOSPlant.PlantTimestamp,
) where {T}
    endpoint == prepared.endpoint ||
        throw(AOSPlant.PlantCommandError(:physical_application,
            :endpoint_mismatch, "command benchmark received another endpoint"))
    workspace.staged = T(command)
    return nothing
end

function AOSPlant.commit_controllable_optic_command!(
    ::PreparedCommandPlaneOptic{T},
    state::CommandPlaneOpticState{T},
    workspace::CommandPlaneOpticWorkspace{T},
    ::AOSPlant.CommandEndpointID,
    ::AOSPlant.PlantTimestamp,
) where {T}
    state.visible = workspace.staged
    return nothing
end

function AOSPlant.apply_controllable_optic_surface!(
    input::AOS.PupilFunction,
    prepared::PreparedCommandPlaneOptic,
    state::CommandPlaneOpticState,
)
    @. input.opd += state.visible * prepared.pattern
    return input
end

function AOSPlant.prepare_path_executor(
    ::CommandSciencePathModel,
    definition::AOSPlant.OpticalPathDefinition,
    source::AOS.AbstractSource,
    telescope::AOS.Telescope,
    atmosphere::AOS.AbstractTimedAtmosphere,
)
    T = eltype(AOS.pupil_reflectivity(telescope))
    pupil = AOS.PupilFunction(telescope; T, backend=AOS.backend(telescope))
    imaging = AOS.prepare_direct_imaging(pupil, source; zero_padding=1)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        AOS.direct_imaging_output(imaging),
        imaging;
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere, telescope, source, pupil),
        optical_model=:gate4_command_direct_imaging,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function AOSPlant.prepare_acquisition_provider(
    model::CommandFrameAcquisitionModel,
    ::AOSPlant.AcquisitionDefinition,
    path::AOSPlant.PreparedPathExecutor,
)
    AOSPlant.require_path_result(path)
    T = eltype(path.result.values)
    detector = AOS.Detector(
        integration_time=T(model.exposure_s),
        noise=AOS.NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=AOS.NullFrameResponse(),
        sensor=AOS.CMOSSensor(timing_model=AOS.GlobalShutter(), T=T),
        T=T,
        backend=path.key.backend,
    )
    execution = AOSPlant.FrameAcquisitionExecution(detector, path.result)
    products = AOSPlant.AcquisitionProducts(execution.observation;
        metadata=(
            kind=:gate4_command_frame,
            units=:detected_electrons,
            geometry=path.result.metadata,
            semantics=:complete_acquisition,
        ))
    return AOSPlant.prepare_full_optical_provider(execution, products)
end

function command_schema(endpoint::Symbol)
    return AOSPlant.PlantCommandSchema(
        Float64,
        ();
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=AOSPlant.CommandBasis(:modal, endpoint),
        basis_revision=1,
        semantics=AOSPlant.AbsoluteCommand,
        bounds=AOSPlant.UniformCommandBounds(-1.0, 1.0),
        value_policy=AOSPlant.CommandValuePolicy(
            range_stage=AOSPlant.EnforceOnApplication,
            out_of_range=AOSPlant.RejectInvalidCommand),
        sequence_policy=AOSPlant.CommandSequencePolicy(),
        effective_time_policy=AOSPlant.CommandEffectiveTimePolicy(
            supersession=AOSPlant.PreservePendingCommands),
        silence_policy=AOSPlant.CommandSilencePolicy(),
    )
end

@inline ns_seconds(value::Integer) = Float64(value) / 1.0e9

function command_plant_definition(raw::AbstractDict;
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
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction=T[0],
        altitude=T[0],
        layer_ids=(:ground,),
        T=T,
    )
    source = AOS.Source(
        band=:custom,
        wavelength=T(raw["science_wavelength_m"]),
        photon_irradiance=T(raw["science_photon_irradiance"]),
        T=T,
    )
    path = AOSPlant.OpticalPathDefinition(
        :science, source, CommandSciencePathModel())
    acquisition = AOSPlant.AcquisitionDefinition(
        :camera, :science,
        CommandFrameAcquisitionModel(
            ns_seconds(raw["exposure_duration_ns"])))
    first_schema = command_schema(:a_low_order)
    second_schema = command_schema(:b_high_order)
    first = AOSPlant.ControllableOpticDefinition(
        :low_order,
        CommandPlaneOpticModel(:x, T(raw["low_order_opd_gain_m"])),
        (first_schema,))
    second = AOSPlant.ControllableOpticDefinition(
        :high_order,
        CommandPlaneOpticModel(:y, T(raw["high_order_opd_gain_m"])),
        (second_schema,))
    optics = reverse_declarations ? (second, first) : (first, second)
    definition = AOSPlant.PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        controllable_optics=optics,
        paths=(path,),
        acquisitions=(acquisition,),
    )
    configurations = (
        AOSPlant.CommandEndpointConfiguration(
            :a_low_order, 0.0;
            capacity=Int(raw["command_capacity"])),
        AOSPlant.CommandEndpointConfiguration(
            :b_high_order, 0.0;
            capacity=Int(raw["command_capacity"])),
    )
    reverse_declarations && (configurations = reverse(configurations))
    return definition, configurations, first_schema, second_schema
end

function command_event_definition(raw::AbstractDict)
    optical_schedule = AOSPlant.PeriodicSchedule(
        period_ns=Int64(raw["optical_sample_period_ns"]),
        phase_ns=0)
    acquisition_schedule = AOSPlant.PeriodicSchedule(
        period_ns=Int64(raw["acquisition_period_ns"]),
        phase_ns=0)
    return AOSPlant.PlantEventLoopDefinition(
        (AOSPlant.OpticalSampleDefinition(
            :science, optical_schedule),),
        (AOSPlant.AcquisitionEventDefinition(
            :camera,
            AOSPlant.GlobalShutterAcquisitionDefinition(
                AOSPlant.PlantDuration(raw["exposure_duration_ns"])),
            AOSPlant.PeriodicAcquisitionStart(acquisition_schedule)),),
    )
end

function prepare_command_plant(raw::AbstractDict;
    reverse_declarations::Bool=false)
    definition, configurations, first_schema, second_schema =
        command_plant_definition(raw; reverse_declarations)
    plant = AOSPlant.prepare_plant(definition;
        run_seed=UInt64(raw["run_seed"]),
        command_endpoints=configurations)
    prepared = AOSPlant.prepare_plant_event_loop(
        plant, command_event_definition(raw))
    return plant, prepared, AOSPlant.PlantEventLoopState(prepared),
        AOSPlant.PlantEventLoopWorkspace(prepared), first_schema,
        second_schema
end

mutable struct CommandPlantOperation{P,S,W,A,B}
    prepared::P
    state::S
    workspace::W
    first_schema::A
    second_schema::B
    first_stride::UInt64
    second_stride::UInt64
    cycles::UInt64
    first_sequence::UInt64
    second_sequence::UInt64
    presented::UInt64
    terminal::UInt64
    applied::UInt64
    maximum_pending::UInt64
    maximum_dispositions::UInt64
end

function prepare_command_operation(raw::AbstractDict;
    reverse_declarations::Bool=false)
    plant, prepared, state, workspace, first_schema, second_schema =
        prepare_command_plant(raw; reverse_declarations)
    first_stride = UInt64(raw["low_order_command_stride"])
    second_stride = UInt64(raw["high_order_command_stride"])
    iszero(first_stride) &&
        error("low-order command stride must be positive")
    iszero(second_stride) &&
        error("high-order command stride must be positive")
    operation = CommandPlantOperation(
        prepared,
        state,
        workspace,
        first_schema,
        second_schema,
        first_stride,
        second_stride,
        UInt64(0),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        UInt64(0),
    )
    return plant, operation
end

@inline function _is_due(cycle::UInt64, stride::UInt64)
    return rem(cycle - UInt64(1), stride) == UInt64(0)
end

@inline function _command_value(sequence::UInt64, magnitude::Float64)
    return isodd(sequence) ? magnitude : -magnitude
end

@inline function _require_admitted(status, label::AbstractString)
    status in (AOSPlant.CommandAdmittedPending,
        AOSPlant.CommandAdmittedReady) ||
        error("Gate 4 $label command was not admitted")
    return nothing
end

@inline function _record_pending_depth!(operation::CommandPlantOperation)
    pending = UInt64(0)
    @inbounds for state in operation.state.command_endpoints
        pending += UInt64(AOSPlant.pending_command_count(state))
    end
    operation.maximum_pending = max(operation.maximum_pending, pending)
    return pending
end

@inline function _consume_terminal_dispositions!(
    operation::CommandPlantOperation)
    count = AOSPlant.command_disposition_count(operation.workspace)
    operation.maximum_dispositions = max(
        operation.maximum_dispositions, UInt64(count))
    @inbounds for index in 1:count
        disposition = AOSPlant.command_disposition(
            operation.workspace, index)
        AOSPlant.command_terminal_kind(disposition) ==
            AOSPlant.AppliedCommand || error(
                "Gate 4 benchmark received a non-applied terminal disposition")
        operation.applied += UInt64(1)
    end
    operation.terminal += UInt64(count)
    AOSPlant.clear_command_dispositions!(operation.workspace)
    return count
end

@inline function _admit_due_commands!(
    operation::CommandPlantOperation,
    cycle::UInt64,
    effective_timestamp::AOSPlant.PlantTimestamp,
    presentation_timestamp::AOSPlant.PlantTimestamp,
)
    if _is_due(cycle, operation.first_stride)
        operation.first_sequence += UInt64(1)
        admission = AOSPlant.admit_plant_command!(
            operation.prepared,
            operation.state,
            operation.workspace,
            AOSPlant.PlantCommand(
                operation.first_schema,
                operation.first_sequence,
                effective_timestamp,
                _command_value(operation.first_sequence, 0.25)),
            presentation_timestamp,
        )
        _require_admitted(AOSPlant.command_admission_status(admission),
            "low-order")
        operation.presented += UInt64(1)
    end
    if _is_due(cycle, operation.second_stride)
        operation.second_sequence += UInt64(1)
        admission = AOSPlant.admit_plant_command!(
            operation.prepared,
            operation.state,
            operation.workspace,
            AOSPlant.PlantCommand(
                operation.second_schema,
                operation.second_sequence,
                effective_timestamp,
                _command_value(operation.second_sequence, 0.15)),
            presentation_timestamp,
        )
        _require_admitted(AOSPlant.command_admission_status(admission),
            "high-order")
        operation.presented += UInt64(1)
    end
    return nothing
end

@inline function (operation::CommandPlantOperation)()
    cycle = operation.cycles + UInt64(1)
    next_timestamp = AOSPlant.next_plant_event_timestamp(
        operation.prepared, operation.state, operation.workspace)
    next_timestamp === nothing &&
        error("Gate 4 command benchmark exhausted plant events")
    current = AOSPlant.scheduler_timestamp(operation.state.scheduler)
    presentation_timestamp = next_timestamp == current ? current :
        current + AOSPlant.PlantDuration(1)
    _admit_due_commands!(operation, cycle, next_timestamp,
        presentation_timestamp)
    _record_pending_depth!(operation)
    processed = AOSPlant.step_plant_events!(
        operation.prepared, operation.state, operation.workspace)
    processed == next_timestamp ||
        error("Gate 4 command benchmark processed another timestamp")
    _consume_terminal_dispositions!(operation)
    operation.cycles = cycle
    return processed
end

function run_command_cycles!(operation::CommandPlantOperation, cycles::Int)
    cycles >= 0 || error("command-cycle count must be nonnegative")
    last_timestamp = zero(AOSPlant.PlantTimestamp)
    @inbounds for _ in 1:cycles
        last_timestamp = operation()
    end
    return last_timestamp
end

function storage_signature(operation::CommandPlantOperation)
    prepared = operation.prepared
    state = operation.state
    workspace = operation.workspace
    return (
        Base.summarysize(prepared.scheduler.definitions),
        Base.summarysize(prepared.actions),
        Base.summarysize(prepared.optics),
        Base.summarysize(prepared.command_endpoints),
        Base.summarysize(prepared.paths),
        Base.summarysize(prepared.acquisitions),
        Base.summarysize(state.scheduler.cursors),
        Base.summarysize(state.command_endpoints),
        Base.summarysize(state.command_applications),
        Base.summarysize(state.controllable_optics),
        Base.summarysize(state.acquisitions),
        Base.summarysize(workspace.scheduler.due_slots),
        Base.summarysize(workspace.command_endpoints),
        Base.summarysize(workspace.controllable_optics),
        Base.summarysize(workspace.command_dispositions),
        Base.summarysize(workspace.transaction_endpoint_slots),
        Base.summarysize(workspace.transaction_admissions),
        Base.summarysize(workspace.transaction_claims),
        Base.summarysize(workspace.transaction_staged),
        Base.summarysize(workspace.due_paths),
    )
end

function endpoint_capacities(operation::CommandPlantOperation)
    return Int[
        AOSPlant.command_endpoint_capacity(binding.binding.endpoint)
        for binding in operation.prepared.command_endpoints
    ]
end

function observation_array(plant)
    return Array(AOSPlant.acquisition_observation(
        AOSPlant.prepared_acquisition(plant, :camera)))
end

function path_array(plant)
    return Array(AOSPlant.path_result(
        AOSPlant.prepared_path(plant, :science)).values)
end

function array_sha256(array)
    return bytes2hex(SHA.sha256(reinterpret(UInt8, vec(array))))
end

function validate_independent_exposure(raw::AbstractDict;
    reverse_declarations::Bool=false)
    plant, prepared, state, workspace, first_schema, second_schema =
        prepare_command_plant(raw; reverse_declarations)
    sample_period_ns = Int64(raw["optical_sample_period_ns"])
    exposure_duration_ns = Int64(raw["exposure_duration_ns"])
    exposure_duration_ns == 3 * sample_period_ns ||
        error("independent exposure oracle requires three optical intervals")
    first_timestamp = AOSPlant.PlantTimestamp(sample_period_ns)
    second_timestamp = AOSPlant.PlantTimestamp(2 * sample_period_ns)
    close_timestamp = AOSPlant.PlantTimestamp(exposure_duration_ns)
    first_admission = AOSPlant.admit_plant_command!(
        prepared, state, workspace,
        AOSPlant.PlantCommand(first_schema, 1,
            first_timestamp, 0.35),
        AOSPlant.PlantTimestamp(0))
    second_admission = AOSPlant.admit_plant_command!(
        prepared, state, workspace,
        AOSPlant.PlantCommand(second_schema, 1,
            second_timestamp, -0.2),
        AOSPlant.PlantTimestamp(0))
    AOSPlant.command_admission_status(first_admission) ==
        AOSPlant.CommandAdmittedPending ||
        error("independent low-order command was not admitted")
    AOSPlant.command_admission_status(second_admission) ==
        AOSPlant.CommandAdmittedPending ||
        error("independent high-order command was not admitted")

    AOSPlant.step_plant_events!(prepared, state, workspace) ==
        AOSPlant.PlantTimestamp(0) ||
        error("independent exposure missed its initial timestamp")
    baseline = copy(AOSPlant.path_result(
        AOSPlant.prepared_path(plant, :science)).values)
    AOSPlant.step_plant_events!(prepared, state, workspace) ==
        first_timestamp ||
        error("independent exposure missed the first command timestamp")
    first_rate = copy(AOSPlant.path_result(
        AOSPlant.prepared_path(plant, :science)).values)
    AOSPlant.step_plant_events!(prepared, state, workspace) ==
        second_timestamp ||
        error("independent exposure missed the second command timestamp")
    second_rate = copy(AOSPlant.path_result(
        AOSPlant.prepared_path(plant, :science)).values)
    AOSPlant.step_plant_events!(prepared, state, workspace) ==
        close_timestamp ||
        error("independent exposure missed its close timestamp")
    observation = observation_array(plant)
    interval_seconds = ns_seconds(sample_period_ns)
    expected = interval_seconds .* (
        baseline .+ first_rate .+ second_rate)
    isapprox(observation, expected; atol=1e-10, rtol=1e-10) ||
        error("command-during-exposure integration changed")
    AOSPlant.command_disposition_count(workspace) == 2 ||
        error("independent commands did not produce two dispositions")
    @inbounds for index in 1:2
        AOSPlant.command_terminal_kind(
            AOSPlant.command_disposition(workspace, index)) ==
            AOSPlant.AppliedCommand ||
            error("independent command was not applied")
    end
    return (
        baseline=baseline,
        first_rate=first_rate,
        second_rate=second_rate,
        observation=observation,
    )
end

function validate_independent_exposure_and_reordering(raw::AbstractDict)
    canonical = validate_independent_exposure(raw)
    reordered = validate_independent_exposure(raw;
        reverse_declarations=true)
    canonical == reordered ||
        error("command exposure changed after declaration reordering")
    return Dict{String,Any}(
        "independent_command_timing" => true,
        "command_during_exposure" => true,
        "declaration_reordering_replay" => true,
        "observation_sha256" => array_sha256(canonical.observation),
    )
end

function validate_atomic_application(raw::AbstractDict)
    _, prepared, state, workspace, first_schema, second_schema =
        prepare_command_plant(raw)
    timestamp = AOSPlant.PlantTimestamp(0)
    transaction = AOSPlant.PlantCommandTransaction(
        AOSPlant.PlantCommand(first_schema, 1, timestamp, 0.3),
        AOSPlant.PlantCommand(second_schema, 1, timestamp, -0.2),
    )
    admission = AOSPlant.admit_plant_command_transaction!(
        prepared, state, workspace, transaction, timestamp)
    _require_admitted(AOSPlant.command_admission_status(admission),
        "atomic transaction")
    AOSPlant.command_transaction_member_count(admission) == 2 ||
        error("atomic command transaction member count changed")
    AOSPlant.step_plant_events!(prepared, state, workspace) == timestamp ||
        error("atomic command transaction missed its effective timestamp")
    AOSPlant.effective_command(prepared, state, :a_low_order) == 0.3 ||
        error("atomic low-order command was not visible")
    AOSPlant.effective_command(prepared, state, :b_high_order) == -0.2 ||
        error("atomic high-order command was not visible")
    AOSPlant.command_disposition_count(workspace) == 2 ||
        error("atomic command transaction did not produce two dispositions")
    @inbounds for index in 1:2
        disposition = AOSPlant.command_disposition(workspace, index)
        AOSPlant.command_terminal_kind(disposition) ==
            AOSPlant.AppliedCommand ||
            error("atomic command transaction member was not applied")
        AOSPlant.command_terminal_timestamp(disposition) == timestamp ||
            error("atomic command terminal timestamp changed")
    end
    return Dict{String,Any}(
        "atomic_multi_optic_application" => true,
        "transaction_members" => 2,
        "terminal_dispositions" => 2,
        "effective_timestamp_ns" => 0,
    )
end

function validate_operation_replay(raw::AbstractDict, cycles::Int)
    canonical_plant, canonical = prepare_command_operation(raw)
    reordered_plant, reordered = prepare_command_operation(raw;
        reverse_declarations=true)
    run_command_cycles!(canonical, cycles)
    run_command_cycles!(reordered, cycles)
    canonical.presented == canonical.terminal == canonical.applied ||
        error("canonical command replay lost terminal accounting")
    reordered.presented == reordered.terminal == reordered.applied ||
        error("reordered command replay lost terminal accounting")
    canonical.presented == reordered.presented ||
        error("command presentation count changed after reordering")
    AOSPlant.effective_command(canonical.prepared, canonical.state,
        :a_low_order) ==
        AOSPlant.effective_command(reordered.prepared, reordered.state,
            :a_low_order) ||
        error("low-order command replay changed after reordering")
    AOSPlant.effective_command(canonical.prepared, canonical.state,
        :b_high_order) ==
        AOSPlant.effective_command(reordered.prepared, reordered.state,
            :b_high_order) ||
        error("high-order command replay changed after reordering")
    path_array(canonical_plant) == path_array(reordered_plant) ||
        error("command-responsive optical path changed after reordering")
    observation_array(canonical_plant) ==
        observation_array(reordered_plant) ||
        error("command-responsive acquisition changed after reordering")
    return Dict{String,Any}(
        "cycles" => cycles,
        "presented_commands" => Int64(canonical.presented),
        "terminal_dispositions" => Int64(canonical.terminal),
        "applied_commands" => Int64(canonical.applied),
        "declaration_reordering_replay" => true,
        "path_sha256" => array_sha256(path_array(canonical_plant)),
        "observation_sha256" =>
            array_sha256(observation_array(canonical_plant)),
    )
end

function validate_fixed_storage(raw::AbstractDict, warmup_cycles::Int,
    long_run_cycles::Int)
    _, operation = prepare_command_operation(raw)
    run_command_cycles!(operation, warmup_cycles)
    before = storage_signature(operation)
    presented_before = operation.presented
    terminal_before = operation.terminal
    run_command_cycles!(operation, long_run_cycles)
    after = storage_signature(operation)
    before == after ||
        error("Gate 4 prepared storage changed with command-cycle count")
    presented_delta = operation.presented - presented_before
    terminal_delta = operation.terminal - terminal_before
    presented_delta == terminal_delta ||
        error("Gate 4 long run lost a terminal command disposition")
    operation.presented == operation.terminal == operation.applied ||
        error("Gate 4 cumulative command accounting is not exact")
    all(state -> AOSPlant.pending_command_count(state) == 0,
        operation.state.command_endpoints) ||
        error("Gate 4 long run retained pending commands")
    capacities = endpoint_capacities(operation)
    operation.maximum_pending <= UInt64(length(capacities)) ||
        error("Gate 4 pending depth exceeded one command per endpoint")
    operation.maximum_dispositions <= UInt64(length(capacities)) ||
        error("Gate 4 terminal ledger exceeded the endpoint count")
    return Dict{String,Any}(
        "fixed_storage" => true,
        "storage_component_bytes" => collect(before),
        "endpoint_count" =>
            AOSPlant.plant_event_command_endpoint_count(operation.prepared),
        "endpoint_capacities" => capacities,
        "maximum_pending_commands" => Int64(operation.maximum_pending),
        "maximum_terminal_dispositions" =>
            Int64(operation.maximum_dispositions),
        "long_run_cycles" => long_run_cycles,
        "presented_commands" => Int64(presented_delta),
        "terminal_dispositions" => Int64(terminal_delta),
        "exact_terminal_accounting" => true,
    )
end

end # module
