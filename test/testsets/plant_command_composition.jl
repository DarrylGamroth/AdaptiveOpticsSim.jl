struct CommandCompositionOpticModel{T<:AbstractFloat}
    axis::Symbol
    gain::T
end

struct CommandCompositionPathModel end

struct CommandCompositionAcquisitionModel{T<:AbstractFloat}
    exposure::T
end

struct ArrayInitialCommandOpticModel end

struct PreparedCommandCompositionOptic{T<:AbstractFloat,
    A<:AbstractMatrix{T}}
    endpoint::CommandEndpointID
    pattern::A
end

struct PreparedArrayInitialCommandOptic{T<:AbstractFloat}
    endpoint::CommandEndpointID
    command_count::Int
end

mutable struct CommandCompositionOpticState{T<:AbstractFloat}
    visible::T
end

mutable struct ArrayInitialCommandOpticState{A<:AbstractVector}
    visible::A
end

mutable struct CommandCompositionOpticWorkspace{T<:AbstractFloat}
    staged::T
end

mutable struct ArrayInitialCommandOpticWorkspace{T<:AbstractFloat}
    staged::Vector{T}
end

Plant.plant_model_definition_style(
    ::Type{<:CommandCompositionOpticModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{ArrayInitialCommandOpticModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{CommandCompositionPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:CommandCompositionAcquisitionModel}) =
    ColdPlantModelDefinition()

function Plant.prepare_controllable_optic(
    model::CommandCompositionOpticModel{T},
    definition::ControllableOpticDefinition,
    telescope::Telescope,
    ::AdaptiveOpticsSim.AbstractAtmosphere,
) where {T}
    endpoint = only(command_endpoint_ids(definition))
    n = size(pupil_reflectivity(telescope), 1)
    pattern = Matrix{T}(undef, n, n)
    center = T(n + 1) / T(2)
    scale = model.gain / T(n)
    if model.axis === :x
        @inbounds for column in 1:n, row in 1:n
            pattern[row, column] = scale * (T(column) - center)
        end
    elseif model.axis === :y
        @inbounds for column in 1:n, row in 1:n
            pattern[row, column] = scale * (T(row) - center)
        end
    else
        throw(PlantPreparationError(:controllable_optic, :invalid_axis,
            "test controllable optic axis must be :x or :y"))
    end
    return PreparedCommandCompositionOptic(endpoint, pattern)
end

function Plant.prepare_controllable_optic(
    ::ArrayInitialCommandOpticModel,
    definition::ControllableOpticDefinition,
    ::Telescope,
    ::AdaptiveOpticsSim.AbstractAtmosphere,
)
    schema = only(command_schemas(definition))
    T = command_numeric_type(schema)
    dimensions = command_dimensions(schema)
    length(dimensions) == 1 || throw(PlantPreparationError(
        :controllable_optic, :invalid_dimensions,
        "array-initial test optic requires one-dimensional commands"))
    return PreparedArrayInitialCommandOptic{T}(
        command_endpoint_id(schema), only(dimensions))
end

function Plant.prepare_controllable_optic_state(
    prepared::PreparedCommandCompositionOptic{T},
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
) where {T}
    only(endpoint_ids) == prepared.endpoint ||
        throw(PlantPreparationError(:controllable_optic,
            :prepared_binding, "test optic endpoint binding changed"))
    return CommandCompositionOpticState(T(only(initial_commands)))
end

function Plant.prepare_controllable_optic_state(
    prepared::PreparedArrayInitialCommandOptic{T},
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
) where {T}
    only(endpoint_ids) == prepared.endpoint ||
        throw(PlantPreparationError(:controllable_optic,
            :prepared_binding, "array-initial test endpoint binding changed"))
    initial = only(initial_commands)
    length(initial) == prepared.command_count ||
        throw(PlantPreparationError(:controllable_optic,
            :prepared_binding, "array-initial test command size changed"))
    return ArrayInitialCommandOpticState(initial)
end

function Plant.prepare_controllable_optic_workspace(
    ::PreparedCommandCompositionOptic{T}) where {T}
    return CommandCompositionOpticWorkspace(zero(T))
end

function Plant.prepare_controllable_optic_workspace(
    prepared::PreparedArrayInitialCommandOptic{T}) where {T}
    return ArrayInitialCommandOpticWorkspace(
        Vector{T}(undef, prepared.command_count))
end

function Plant.stage_controllable_optic_command!(
    prepared::PreparedCommandCompositionOptic{T},
    ::CommandCompositionOpticState{T},
    workspace::CommandCompositionOpticWorkspace{T},
    endpoint::CommandEndpointID,
    command,
    ::PlantTimestamp,
) where {T}
    endpoint == prepared.endpoint ||
        throw(PlantCommandError(:physical_application,
            :endpoint_mismatch, "test optic received another endpoint"))
    command == T(0.9) &&
        throw(PlantCommandError(:physical_application,
            :test_physical_application_failure,
            "test optic rejected its staged physical response"))
    workspace.staged = T(command)
    return nothing
end

function Plant.stage_controllable_optic_command!(
    prepared::PreparedArrayInitialCommandOptic,
    ::ArrayInitialCommandOpticState,
    workspace::ArrayInitialCommandOpticWorkspace,
    endpoint::CommandEndpointID,
    command::AbstractVector,
    ::PlantTimestamp)
    endpoint == prepared.endpoint ||
        throw(PlantCommandError(:physical_application,
            :endpoint_mismatch, "array-initial test received another endpoint"))
    copyto!(workspace.staged, command)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::PreparedCommandCompositionOptic{T},
    state::CommandCompositionOpticState{T},
    workspace::CommandCompositionOpticWorkspace{T},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {T}
    state.visible = workspace.staged
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::PreparedArrayInitialCommandOptic,
    state::ArrayInitialCommandOpticState,
    workspace::ArrayInitialCommandOpticWorkspace,
    ::CommandEndpointID,
    ::PlantTimestamp)
    copyto!(state.visible, workspace.staged)
    return nothing
end

function Plant.apply_controllable_optic_surface!(
    input::PupilFunction,
    prepared::PreparedCommandCompositionOptic,
    state::CommandCompositionOpticState)
    @. input.opd += state.visible * prepared.pattern
    return input
end

function Plant.apply_controllable_optic_surface!(
    input::PupilFunction,
    ::PreparedArrayInitialCommandOptic,
    state::ArrayInitialCommandOpticState)
    offset = sum(state.visible)
    @. input.opd += offset
    return input
end

function Plant.prepare_path_executor(
    ::CommandCompositionPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source; zero_padding=1)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        materialization=prepare_pupil_opd_materialization(atmosphere,
            telescope, source, pupil),
        optical_model=:command_composition_direct_imaging,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::CommandCompositionAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T), T=T,
        backend=path.key.backend)
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(execution.observation;
        metadata=(kind=:command_composition_frame,
            units=:detected_electrons,
            geometry=path.result.metadata))
    return prepare_full_optical_provider(execution, products)
end

function command_composition_schema(endpoint::Symbol;
    silence_policy=CommandSilencePolicy(),
    effective_time_policy=CommandEffectiveTimePolicy(
        supersession=PreservePendingCommands))
    return PlantCommandSchema(
        Float64,
        ();
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:modal, endpoint),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(-1.0, 1.0),
        value_policy=CommandValuePolicy(
            range_stage=EnforceOnApplication,
            out_of_range=RejectInvalidCommand),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy,
        silence_policy,
    )
end

function command_composition_fixture(; reverse_order::Bool=false,
    first_silence_policy=CommandSilencePolicy(),
    first_effective_time_policy=CommandEffectiveTimePolicy(
        supersession=PreservePendingCommands),
    first_safe_command=nothing)
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[1], wind_speed=T[0], wind_direction=T[0],
        altitude=T[0], layer_ids=(:ground,), T=T)
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(100), T=T)
    path = OpticalPathDefinition(:science, source,
        CommandCompositionPathModel())
    acquisition = AcquisitionDefinition(:camera, :science,
        CommandCompositionAcquisitionModel(T(0.3)))
    first_schema = command_composition_schema(:a_woofer;
        silence_policy=first_silence_policy,
        effective_time_policy=first_effective_time_policy)
    second_schema = command_composition_schema(:b_tweeter)
    first = ControllableOpticDefinition(:woofer,
        CommandCompositionOpticModel(:x, T(8e-7)),
        (first_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility())
    second = ControllableOpticDefinition(:tweeter,
        CommandCompositionOpticModel(:y, T(6e-7)),
        (second_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility())
    optics = reverse_order ? (second, first) : (first, second)
    configurations = (
        CommandEndpointConfiguration(:a_woofer, 0.0; capacity=4,
            safe_command=first_safe_command),
        CommandEndpointConfiguration(:b_tweeter, 0.0; capacity=4),
    )
    reverse_order && (configurations = reverse(configurations))
    definition = PlantDefinition(; telescope, atmosphere,
        controllable_optics=optics, paths=(path,),
        acquisitions=(acquisition,))
    plant = prepare_plant(definition; run_seed=0x7a00,
        command_endpoints=configurations)
    loop_definition = PlantEventLoopDefinition(
        (OpticalSampleDefinition(:science,
            PeriodicSchedule(period_ns=100_000_000, phase_ns=0)),),
        (AcquisitionEventDefinition(:camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=1_000_000_000, phase_ns=0))),),
    )
    prepared = prepare_plant_event_loop(plant, loop_definition)
    return plant, prepared, PlantEventLoopState(prepared),
        PlantEventLoopWorkspace(prepared), first_schema, second_schema
end

function array_initial_command_fixture()
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[1], wind_speed=T[0], wind_direction=T[0],
        altitude=T[0], layer_ids=(:ground,), T=T)
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(100), T=T)
    path = OpticalPathDefinition(:science, source,
        CommandCompositionPathModel())
    acquisition = AcquisitionDefinition(:camera, :science,
        CommandCompositionAcquisitionModel(T(0.3)))
    schema = PlantCommandSchema(
        T,
        (2,);
        id=:array_initial_schema,
        version=1,
        endpoint=:array_initial,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:modal, :array_initial),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(T(-1), T(1)),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=PreservePendingCommands),
        silence_policy=CommandSilencePolicy(),
    )
    optic = ControllableOpticDefinition(:array_initial_optic,
        ArrayInitialCommandOpticModel(), (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility())
    definition = PlantDefinition(; telescope, atmosphere,
        controllable_optics=(optic,), paths=(path,),
        acquisitions=(acquisition,))
    initial = T[0.1, -0.2]
    plant = prepare_plant(definition; run_seed=0x7a01,
        command_endpoints=(
            CommandEndpointConfiguration(:array_initial, initial;
                capacity=2),
        ))
    loop_definition = PlantEventLoopDefinition(
        (OpticalSampleDefinition(:science,
            PeriodicSchedule(period_ns=100_000_000, phase_ns=0)),),
        (AcquisitionEventDefinition(:camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(300_000_000)),
            PeriodicAcquisitionStart(PeriodicSchedule(
                period_ns=1_000_000_000, phase_ns=0))),),
    )
    return initial, plant, prepare_plant_event_loop(plant, loop_definition)
end

function captured_command_composition_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function command_composition_submit!(prepared, state, workspace, schema,
    sequence, effective_ns, value)
    return admit_plant_command!(prepared, state, workspace,
        PlantCommand(schema, sequence, PlantTimestamp(effective_ns), value),
        PlantTimestamp(0))
end

function command_composition_run_independent(; reverse_order=false)
    plant, prepared, state, workspace, first_schema, second_schema =
        command_composition_fixture(; reverse_order)
    first_admission = command_composition_submit!(prepared, state, workspace,
        first_schema, 1, 100_000_000, 0.35)
    second_admission = command_composition_submit!(prepared, state, workspace,
        second_schema, 1, 200_000_000, -0.2)
    @test command_admission_status(first_admission) ==
        CommandAdmittedPending
    @test command_admission_status(second_admission) ==
        CommandAdmittedPending

    @test step_plant_events!(prepared, state, workspace) ==
        PlantTimestamp(0)
    baseline = copy(path_result(prepared_path(plant, :science)).values)
    @test step_plant_events!(prepared, state, workspace) ==
        PlantTimestamp(100_000_000)
    first_rate = copy(path_result(prepared_path(plant, :science)).values)
    @test effective_command(prepared, state, :a_woofer) == 0.35
    @test effective_command(prepared, state, :b_tweeter) == 0.0
    @test step_plant_events!(prepared, state, workspace) ==
        PlantTimestamp(200_000_000)
    second_rate = copy(path_result(prepared_path(plant, :science)).values)
    @test effective_command(prepared, state, :a_woofer) == 0.35
    @test effective_command(prepared, state, :b_tweeter) == -0.2
    @test step_plant_events!(prepared, state, workspace) ==
        PlantTimestamp(300_000_000)
    observation = copy(acquisition_observation(
        prepared_acquisition(plant, :camera)))
    return (; plant, prepared, state, workspace, baseline, first_rate,
        second_rate, observation)
end

@inline function command_composition_step_allocations(
    prepared, state, workspace)
    return @allocated step_plant_events!(prepared, state, workspace)
end

@testset "Prepared command composition API and independent timing" begin
    @test Base.isexported(Plant, :CommandEndpointConfiguration)
    @test !Base.isexported(AdaptiveOpticsSim,
        :CommandEndpointConfiguration)
    for name in (
        :PreparedControllableOptic,
        :prepare_controllable_optic,
        :prepare_controllable_optic_state,
        :prepare_controllable_optic_workspace,
        :stage_controllable_optic_command!,
        :commit_controllable_optic_command!,
        :apply_controllable_optic_surface!,
        :PlantCommandTransaction,
        :PlantCommandTransactionAdmission,
        :admit_plant_command_transaction!,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end

    result = command_composition_run_independent()
    @test plant_event_controllable_optic_count(result.prepared) == 2
    @test plant_event_command_endpoint_count(result.prepared) == 2
    @test length(Plant.prepared_controllable_optics(result.plant)) == 2
    @test length(Plant.prepared_command_endpoints(result.plant)) == 2
    @test Plant.prepared_controllable_optic(result.plant, :woofer) ===
        Plant.prepared_controllable_optics(result.plant)[2]
    @test Plant.prepared_command_endpoint(result.plant, :a_woofer) ===
        Plant.prepared_command_endpoints(result.plant)[1]
    @test command_disposition_count(result.workspace) == 2
    @test all(index -> command_terminal_kind(
            command_disposition(result.workspace, index)) ==
            AppliedCommand, 1:2)
    @test result.baseline != result.first_rate
    @test result.first_rate != result.second_rate
    expected = 0.1 .* (result.baseline .+ result.first_rate .+
        result.second_rate)
    @test result.observation ≈ expected atol=1e-10 rtol=1e-10

    reordered = command_composition_run_independent(reverse_order=true)
    @test reordered.baseline == result.baseline
    @test reordered.first_rate == result.first_rate
    @test reordered.second_rate == result.second_rate
    @test reordered.observation == result.observation
end

@testset "Prepared configuration and physical-state ownership" begin
    plant, _, _, _, _, _ = command_composition_fixture()
    definition = plant.definition

    duplicate_error = captured_command_composition_error() do
        prepare_plant(definition; run_seed=0x7a02,
            command_endpoints=(
                CommandEndpointConfiguration(:a_woofer, 0.0; capacity=2),
                CommandEndpointConfiguration(:a_woofer, 0.0; capacity=2),
            ))
    end
    @test duplicate_error isa PlantPreparationError
    @test duplicate_error.reason == :duplicate_configuration

    missing_error = captured_command_composition_error() do
        prepare_plant(definition; run_seed=0x7a03,
            command_endpoints=(
                CommandEndpointConfiguration(:a_woofer, 0.0; capacity=2),
                CommandEndpointConfiguration(:extra, 0.0; capacity=2),
            ))
    end
    @test missing_error isa PlantPreparationError
    @test missing_error.reason == :missing_configuration

    named_error = captured_command_composition_error() do
        prepare_plant(definition; run_seed=0x7a04,
            command_endpoints=(
                wrong=CommandEndpointConfiguration(
                    :a_woofer, 0.0; capacity=2),
                b_tweeter=CommandEndpointConfiguration(
                    :b_tweeter, 0.0; capacity=2),
            ))
    end
    @test named_error isa PlantPreparationError
    @test named_error.reason == :identity_mismatch

    _, timed, timed_state, timed_workspace, timed_first, timed_second =
        command_composition_fixture()
    overtake_admission_error = captured_command_composition_error() do
        admit_plant_command!(timed, timed_state, timed_workspace,
            PlantCommand(timed_first, 1, PlantTimestamp(100_000_000), 0.1),
            PlantTimestamp(1))
    end
    @test overtake_admission_error isa PlantScheduleError
    @test overtake_admission_error.reason ==
        :command_admission_overtakes_event
    @test command_disposition_count(timed_workspace) == 0
    overtake_transaction_error = captured_command_composition_error() do
        Plant.admit_plant_command_transaction!(timed, timed_state,
            timed_workspace, Plant.PlantCommandTransaction(
                PlantCommand(timed_first, 1,
                    PlantTimestamp(100_000_000), 0.1),
                PlantCommand(timed_second, 1,
                    PlantTimestamp(100_000_000), -0.1),
            ), PlantTimestamp(1))
    end
    @test overtake_transaction_error isa PlantScheduleError
    @test overtake_transaction_error.reason ==
        :command_admission_overtakes_event
    @test command_disposition_count(timed_workspace) == 0
    @test next_plant_event_timestamp(timed, timed_state, timed_workspace) ==
        PlantTimestamp(0)
    @test step_plant_events!(timed, timed_state, timed_workspace) ==
        PlantTimestamp(0)

    elapsed_admission_error = captured_command_composition_error() do
        admit_plant_command!(timed, timed_state, timed_workspace,
            PlantCommand(timed_first, 1, PlantTimestamp(100_000_000), 0.1),
            PlantTimestamp(0))
    end
    @test elapsed_admission_error isa PlantScheduleError
    @test elapsed_admission_error.reason == :command_admission_time_elapsed
    @test command_disposition_count(timed_workspace) == 0
    between_event_admission = admit_plant_command!(
        timed, timed_state, timed_workspace,
        PlantCommand(timed_first, 1, PlantTimestamp(100_000_000), 0.1),
        PlantTimestamp(1))
    @test command_admission_status(between_event_admission) ==
        CommandAdmittedPending
    @test command_scheduled_timestamp(
        command_order_key(between_event_admission)) ==
        PlantTimestamp(100_000_000)
    @test step_plant_events!(timed, timed_state, timed_workspace) ==
        PlantTimestamp(100_000_000)
    @test command_disposition_count(timed_workspace) == 1
    clear_command_dispositions!(timed_workspace)
    regressed_admission_error = captured_command_composition_error() do
        admit_plant_command!(timed, timed_state, timed_workspace,
            PlantCommand(timed_first, 2, PlantTimestamp(200_000_000), 0.2),
            PlantTimestamp(50_000_000))
    end
    @test regressed_admission_error isa PlantScheduleError
    @test regressed_admission_error.reason ==
        :command_admission_time_regression
    @test command_disposition_count(timed_workspace) == 0

    initial, _, prepared = array_initial_command_fixture()
    initial .= 0.8
    first_state = PlantEventLoopState(prepared)
    second_state = PlantEventLoopState(prepared)
    first_visible = first_state.controllable_optics[1].visible
    second_visible = second_state.controllable_optics[1].visible
    @test first_visible == [0.1, -0.2]
    @test second_visible == [0.1, -0.2]
    @test first_visible !== second_visible
    first_visible[1] = 0.7
    @test second_visible == [0.1, -0.2]
end

@testset "Equal-time right continuity and explicit atomic transactions" begin
    plant, prepared, state, workspace, first_schema, second_schema =
        command_composition_fixture()
    command_composition_submit!(prepared, state, workspace, first_schema,
        1, 0, 0.4)
    @test step_plant_events!(prepared, state, workspace) ==
        PlantTimestamp(0)
    @test effective_command(prepared, state, :a_woofer) == 0.4
    @test command_terminal_timestamp(command_disposition(workspace, 1)) ==
        PlantTimestamp(0)
    @test path_result(prepared_path(plant, :science)).values !=
        command_composition_run_independent().baseline

    _, atomic, atomic_state, atomic_workspace, atomic_first, atomic_second =
        command_composition_fixture()
    transaction = Plant.PlantCommandTransaction(
        PlantCommand(atomic_second, 1, PlantTimestamp(100_000_000), -0.3),
        PlantCommand(atomic_first, 1, PlantTimestamp(100_000_000), 0.25),
    )
    admission = Plant.admit_plant_command_transaction!(atomic, atomic_state,
        atomic_workspace, transaction, PlantTimestamp(0))
    @test command_admission_status(admission) == CommandAdmittedPending
    @test Plant.command_transaction_id(admission) !== nothing
    @test Plant.command_transaction_member_count(admission) == 2
    collision = admit_plant_command!(atomic, atomic_state, atomic_workspace,
        PlantCommand(atomic_first, 2, PlantTimestamp(100_000_000), 0.5),
        PlantTimestamp(0))
    @test command_admission_status(collision) ==
        CommandTerminatedOnAdmission
    @test command_disposition_count(atomic_workspace) == 1
    @test command_terminal_kind(
        command_disposition(atomic_workspace, 1)) == RejectedCommand
    @test command_disposition_reason(
        command_disposition(atomic_workspace, 1)).name ==
        :equal_time_endpoint_conflict
    clear_command_dispositions!(atomic_workspace)
    @test step_plant_events!(atomic, atomic_state, atomic_workspace) ==
        PlantTimestamp(0)
    @test step_plant_events!(atomic, atomic_state, atomic_workspace) ==
        PlantTimestamp(100_000_000)
    @test effective_command(atomic, atomic_state, :a_woofer) == 0.25
    @test effective_command(atomic, atomic_state, :b_tweeter) == -0.3
    @test command_disposition_count(atomic_workspace) == 2
    @test all(index -> command_terminal_kind(
            command_disposition(atomic_workspace, index)) ==
            AppliedCommand, 1:2)

    _, rejected, rejected_state, rejected_workspace, rejected_first,
        rejected_second = command_composition_fixture()
    rejected_transaction = Plant.PlantCommandTransaction(
        PlantCommand(rejected_first, 1, PlantTimestamp(100_000_000), 0.2),
        PlantCommand(rejected_second, 1, PlantTimestamp(100_000_000), 2.0),
    )
    Plant.admit_plant_command_transaction!(rejected, rejected_state,
        rejected_workspace, rejected_transaction, PlantTimestamp(0))
    step_plant_events!(rejected, rejected_state, rejected_workspace)
    step_plant_events!(rejected, rejected_state, rejected_workspace)
    @test effective_command(rejected, rejected_state, :a_woofer) == 0.0
    @test effective_command(rejected, rejected_state, :b_tweeter) == 0.0
    @test command_disposition_count(rejected_workspace) == 2
    @test all(index -> command_terminal_kind(
            command_disposition(rejected_workspace, index)) ==
            RejectedCommand, 1:2)
    reasons = Set(command_disposition_reason(
        command_disposition(rejected_workspace, index)).name for index in 1:2)
    @test :atomic_transaction_aborted in reasons
    @test :out_of_range_rejected in reasons

    _, conflict, conflict_state, conflict_workspace, conflict_first,
        conflict_second = command_composition_fixture()
    command_composition_submit!(conflict, conflict_state, conflict_workspace,
        conflict_first, 1, 100_000_000, 0.1)
    conflicting_transaction = Plant.PlantCommandTransaction(
        PlantCommand(conflict_first, 2, PlantTimestamp(100_000_000), 0.2),
        PlantCommand(conflict_second, 1, PlantTimestamp(100_000_000), -0.2),
    )
    conflict_admission = Plant.admit_plant_command_transaction!(
        conflict, conflict_state, conflict_workspace,
        conflicting_transaction, PlantTimestamp(0))
    @test command_admission_status(conflict_admission) ==
        CommandTerminatedOnAdmission
    @test Plant.command_transaction_id(conflict_admission) === nothing
    @test command_scheduled_timestamp(conflict_admission) === nothing
    @test command_disposition_count(conflict_workspace) == 2
    @test all(index -> command_terminal_kind(
            command_disposition(conflict_workspace, index)) ==
            RejectedCommand, 1:2)
    @test all(index -> command_disposition_reason(
            command_disposition(conflict_workspace, index)).name ==
            :atomic_transaction_aborted_equal_time_endpoint_conflict, 1:2)
    clear_command_dispositions!(conflict_workspace)
    step_plant_events!(conflict, conflict_state, conflict_workspace)
    step_plant_events!(conflict, conflict_state, conflict_workspace)
    @test effective_command(conflict, conflict_state, :a_woofer) == 0.1
    @test effective_command(conflict, conflict_state, :b_tweeter) == 0.0

    fail_late = CommandEffectiveTimePolicy(
        AllowFutureCommand,
        FailOnLateCommand,
        PreservePendingCommands,
    )
    _, failed, failed_state, failed_workspace, failed_first, failed_second =
        command_composition_fixture(
            first_effective_time_policy=fail_late)
    failed_transaction = Plant.PlantCommandTransaction(
        PlantCommand(failed_first, 1, PlantTimestamp(0), 0.2),
        PlantCommand(failed_second, 1, PlantTimestamp(0), -0.2),
    )
    @test step_plant_events!(failed, failed_state, failed_workspace) ==
        PlantTimestamp(0)
    failure = try
        Plant.admit_plant_command_transaction!(failed, failed_state,
            failed_workspace, failed_transaction, PlantTimestamp(1))
        nothing
    catch error
        error
    end
    @test failure isa PlantCommandError
    @test failure.reason == :late_command
    @test command_disposition_count(failed_workspace) == 2
    @test all(index -> command_terminal_kind(
            command_disposition(failed_workspace, index)) ==
            FailedCommand, 1:2)
    @test all(index -> command_disposition_reason(
            command_disposition(failed_workspace, index)).name ==
            :atomic_transaction_aborted_late_command, 1:2)
    @test effective_command(failed, failed_state, :a_woofer) == 0.0
    @test effective_command(failed, failed_state, :b_tweeter) == 0.0

    _, late, late_state, late_workspace, late_first, _ =
        command_composition_fixture(first_effective_time_policy=fail_late)
    @test step_plant_events!(late, late_state, late_workspace) ==
        PlantTimestamp(0)
    late_failure = try
        admit_plant_command!(late, late_state, late_workspace,
            PlantCommand(late_first, 1, PlantTimestamp(0), 0.2),
            PlantTimestamp(1))
        nothing
    catch error
        error
    end
    @test late_failure isa PlantCommandError
    @test late_failure.reason == :late_command
    @test command_disposition_count(late_workspace) == 1
    @test command_terminal_kind(command_disposition(late_workspace, 1)) ==
        FailedCommand
    @test command_disposition_reason(
        command_disposition(late_workspace, 1)).name == :late_command

    _, routed_failure, routed_failure_state, routed_failure_workspace,
        routed_failure_first, _ = command_composition_fixture()
    command_composition_submit!(routed_failure, routed_failure_state,
        routed_failure_workspace, routed_failure_first, 1, 100_000_000, 0.9)
    step_plant_events!(routed_failure, routed_failure_state,
        routed_failure_workspace)
    physical_failure = try
        step_plant_events!(routed_failure, routed_failure_state,
            routed_failure_workspace)
        nothing
    catch error
        error
    end
    @test physical_failure isa PlantCommandError
    @test physical_failure.reason == :test_physical_application_failure
    @test command_disposition_count(routed_failure_workspace) == 1
    @test command_terminal_kind(
        command_disposition(routed_failure_workspace, 1)) == FailedCommand
    @test command_disposition_reason(
        command_disposition(routed_failure_workspace, 1)).name ==
        :test_physical_application_failure
    @test effective_command(
        routed_failure, routed_failure_state, :a_woofer) == 0.0

    safe_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(150_000_000),
    )
    safe_plant, safe_loop, safe_state, safe_workspace, _, _ =
        command_composition_fixture(
            first_silence_policy=safe_policy,
            first_safe_command=-0.45,
        )
    @test step_plant_events!(safe_loop, safe_state, safe_workspace) ==
        PlantTimestamp(0)
    safe_baseline =
        copy(path_result(prepared_path(safe_plant, :science)).values)
    @test step_plant_events!(safe_loop, safe_state, safe_workspace) ==
        PlantTimestamp(100_000_000)
    @test path_result(prepared_path(safe_plant, :science)).values ==
        safe_baseline
    @test step_plant_events!(safe_loop, safe_state, safe_workspace) ==
        PlantTimestamp(150_000_000)
    @test effective_command(safe_loop, safe_state, :a_woofer) == -0.45
    @test command_disposition_count(safe_workspace) == 0
    @test step_plant_events!(safe_loop, safe_state, safe_workspace) ==
        PlantTimestamp(200_000_000)
    @test path_result(prepared_path(safe_plant, :science)).values !=
        safe_baseline
end

@testset "Composed command path allocation budget" begin
    # Warm the same concrete command/optic/event types before measuring a
    # second independently prepared run.
    command_composition_run_independent()
    plant, prepared, state, workspace, first_schema, _ =
        command_composition_fixture()
    command_composition_submit!(prepared, state, workspace, first_schema,
        1, 100_000_000, 0.15)
    step_plant_events!(prepared, state, workspace)
    if coverage_instrumented()
        @test_skip "command-composition allocation gate disabled under " *
            "coverage instrumentation"
    else
        @test next_plant_event_timestamp(prepared, state, workspace) ==
            PlantTimestamp(100_000_000)
        allocated = command_composition_step_allocations(prepared, state,
            workspace)
        @test allocated <= 4_096
    end
    @test all(isfinite, path_result(prepared_path(plant, :science)).values)
end
