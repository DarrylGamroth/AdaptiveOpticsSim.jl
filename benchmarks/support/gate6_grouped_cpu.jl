module Gate6GroupedCPUBenchmark

import AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant

const AOSPlant = AdaptiveOpticsSim.Plant

@enum WorkerPhase::UInt8 begin
    MaterializationPhase = 0x01
    ExecutionPhase = 0x02
end

struct WorkerCommand
    phase::WorkerPhase
    claim::AOSPlant.OpticalPathBatchClaim
    measure_allocations::Bool
end

struct WorkerCompletion
    phase::WorkerPhase
    exception::Any
    allocated_bytes::Int
    task_id::UInt
end

"""
Benchmark-only fixed-owner coordinator for the qualified core batch lifecycle.

The channels make the experimental ownership boundary explicit and keep one
long-lived task per prepared path group. They are not a proposed production HIL
queue; AdaptiveOpticsHIL owns the eventual bounded SPSC implementation.
"""
mutable struct FixedOwnerBatchExecutor{P,S,W} <:
    AOSPlant.AbstractOpticalPathBatchExecutor
    prepared::P
    state::S
    workspace::W
    commands::Memory{Channel{Union{Nothing,WorkerCommand}}}
    completions::Memory{Channel{WorkerCompletion}}
    tasks::Memory{Task}
    due_ordinals::Memory{Int}
    task_ids::Memory{UInt}
    batch_count::Int
    materialization_count::Int
    execution_count::Int
    measure_allocations::Bool
    measured_call_count::Int
    measured_allocation_bytes::Int
    maximum_call_allocation_bytes::Int
    closed::Bool
end

function execute_worker_command!(
    ordinal::Int,
    prepared,
    state,
    workspace,
    command::WorkerCommand,
)
    if command.phase == MaterializationPhase
        AOSPlant.materialize_path_execution_group!(
            prepared, state, workspace, command.claim, ordinal)
    else
        AOSPlant.execute_path_execution_group!(
            prepared, state, workspace, command.claim, ordinal)
    end
    return nothing
end

function worker_loop!(
    ordinal::Int,
    prepared,
    state,
    workspace,
    commands::Channel{Union{Nothing,WorkerCommand}},
    completions::Channel{WorkerCompletion},
)
    task_id = objectid(current_task())
    while true
        command = take!(commands)
        command === nothing && break
        allocated_bytes = 0
        exception = try
            if command.measure_allocations
                allocated_bytes = @allocated execute_worker_command!(
                    ordinal, prepared, state, workspace, command)
            else
                execute_worker_command!(
                    ordinal, prepared, state, workspace, command)
            end
            nothing
        catch caught
            caught
        end
        put!(completions, WorkerCompletion(
            command.phase, exception, allocated_bytes, task_id))
    end
    return nothing
end

function FixedOwnerBatchExecutor(
    prepared::AOSPlant.PreparedPlantEventLoop,
    state::AOSPlant.PlantEventLoopState,
    workspace::AOSPlant.PlantEventLoopWorkspace,
)
    count = AOSPlant.path_execution_group_count(prepared)
    commands =
        Memory{Channel{Union{Nothing,WorkerCommand}}}(undef, count)
    completions = Memory{Channel{WorkerCompletion}}(undef, count)
    tasks = Memory{Task}(undef, count)
    @inbounds for ordinal in 1:count
        command = Channel{Union{Nothing,WorkerCommand}}(1)
        completion = Channel{WorkerCompletion}(1)
        commands[ordinal] = command
        completions[ordinal] = completion
        tasks[ordinal] = Threads.@spawn worker_loop!(
            ordinal, prepared, state, workspace, command, completion)
    end
    task_ids = Memory{UInt}(undef, count)
    fill!(task_ids, UInt(0))
    return FixedOwnerBatchExecutor(
        prepared,
        state,
        workspace,
        commands,
        completions,
        tasks,
        Memory{Int}(undef, count),
        task_ids,
        0,
        0,
        0,
        false,
        0,
        0,
        0,
        false,
    )
end

function require_executor_owners(
    executor::FixedOwnerBatchExecutor,
    prepared,
    state,
    workspace,
)
    prepared === executor.prepared ||
        error("grouped benchmark received a foreign prepared plant")
    state === executor.state ||
        error("grouped benchmark received a foreign state owner")
    workspace === executor.workspace ||
        error("grouped benchmark received a foreign workspace owner")
    executor.closed &&
        error("grouped benchmark executor is closed")
    return nothing
end

function dispatch_phase!(
    executor::FixedOwnerBatchExecutor,
    claim::AOSPlant.OpticalPathBatchClaim,
    due_count::Int,
    phase::WorkerPhase,
    reverse_order::Bool,
)
    @inbounds for index in 1:due_count
        due_index = reverse_order ? due_count - index + 1 : index
        ordinal = executor.due_ordinals[due_index]
        put!(executor.commands[ordinal], WorkerCommand(
            phase, claim, executor.measure_allocations))
    end
    @inbounds for index in 1:due_count
        ordinal = executor.due_ordinals[index]
        completion = take!(executor.completions[ordinal])
        completion.phase == phase ||
            error("grouped benchmark worker completed the wrong phase")
        known_task_id = executor.task_ids[ordinal]
        if iszero(known_task_id)
            executor.task_ids[ordinal] = completion.task_id
        else
            completion.task_id == known_task_id ||
                error("prepared path group changed task owner")
        end
        completion.exception === nothing ||
            throw(completion.exception)
        if executor.measure_allocations
            executor.measured_call_count += 1
            executor.measured_allocation_bytes +=
                completion.allocated_bytes
            executor.maximum_call_allocation_bytes = max(
                executor.maximum_call_allocation_bytes,
                completion.allocated_bytes,
            )
        end
    end
    if phase == MaterializationPhase
        executor.materialization_count += due_count
    else
        executor.execution_count += due_count
    end
    return nothing
end

function AOSPlant.execute_optical_path_batch!(
    executor::FixedOwnerBatchExecutor,
    prepared::AOSPlant.PreparedPlantEventLoop,
    state::AOSPlant.PlantEventLoopState,
    workspace::AOSPlant.PlantEventLoopWorkspace,
    timestamp::AOSPlant.PlantTimestamp,
)
    require_executor_owners(executor, prepared, state, workspace)
    claim = AOSPlant.begin_optical_path_batch!(
        prepared, state, workspace, timestamp)
    due_count = AOSPlant.optical_path_batch_due_group_count(
        prepared, state, workspace, claim)
    @inbounds for index in 1:due_count
        executor.due_ordinals[index] =
            AOSPlant.optical_path_batch_due_group_ordinal(
                prepared, state, workspace, claim, index)
    end
    reverse_materialization = iseven(executor.batch_count)
    dispatch_phase!(executor, claim, due_count, MaterializationPhase,
        reverse_materialization)
    AOSPlant.seal_optical_path_batch_materialization!(
        prepared, state, workspace, claim)
    dispatch_phase!(executor, claim, due_count, ExecutionPhase,
        !reverse_materialization)
    timestamp = AOSPlant.complete_optical_path_batch!(
        prepared, state, workspace, claim)
    executor.batch_count += 1
    return timestamp
end

function reset_allocation_measurement!(
    executor::FixedOwnerBatchExecutor;
    enabled::Bool=true,
)
    executor.measure_allocations = enabled
    executor.measured_call_count = 0
    executor.measured_allocation_bytes = 0
    executor.maximum_call_allocation_bytes = 0
    return executor
end

function close_executor!(executor::FixedOwnerBatchExecutor)
    executor.closed && return executor
    @inbounds for command in executor.commands
        put!(command, nothing)
    end
    @inbounds for task in executor.tasks
        wait(task)
    end
    executor.closed = true
    return executor
end

mutable struct GroupedPlantOperation{O,E}
    serial_storage::O
    executor::E
end

function GroupedPlantOperation(operation)
    executor = FixedOwnerBatchExecutor(
        operation.prepared, operation.state, operation.workspace)
    return GroupedPlantOperation(operation, executor)
end

@inline function (operation::GroupedPlantOperation)()
    storage = operation.serial_storage
    timestamp = AOSPlant.step_plant_events!(
        storage.prepared,
        storage.state,
        storage.workspace,
        operation.executor,
    )
    timestamp === nothing &&
        error("Gate 6 grouped plant exhausted events")
    storage.processed_timestamps += UInt64(1)
    return timestamp
end

function run_timestamp_window!(
    operation::GroupedPlantOperation,
    timestamps::Int,
)
    timestamps >= 0 ||
        error("timestamp window must be nonnegative")
    last_timestamp = zero(AOSPlant.PlantTimestamp)
    @inbounds for _ in 1:timestamps
        last_timestamp = operation()
    end
    return last_timestamp
end

end # module
