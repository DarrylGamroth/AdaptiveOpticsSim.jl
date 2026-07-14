abstract type AbstractExecutionPolicy end

"""Execute coarse work directly on the calling task."""
struct SequentialExecution <: AbstractExecutionPolicy end

"""Execute independent coarse work on Julia threads."""
struct ThreadedExecution <: AbstractExecutionPolicy end

"""Delegate independent backend work to backend-owned streams."""
struct BackendStreamExecution <: AbstractExecutionPolicy end

"""
    DeterministicExecution()

Execute ensemble members in a fixed serial order. Construction of an ensemble
with this policy requires a single-threaded Julia process and sets FFT and BLAS
thread counts to one.
"""
struct DeterministicExecution <: AbstractExecutionPolicy end

"""
    AcceleratedKernelsExecution(; max_tasks=Threads.nthreads(), min_members_per_task=1)

Opt-in coarse CPU task partitioning through AcceleratedKernels.jl. This policy
requires the `AcceleratedKernels` weak dependency to be loaded.
"""
struct AcceleratedKernelsExecution <: AbstractExecutionPolicy
    max_tasks::Int
    min_members_per_task::Int

    function AcceleratedKernelsExecution(max_tasks::Integer,
        min_members_per_task::Integer)
        max_tasks > 0 ||
            throw(InvalidConfiguration("max_tasks must be greater than zero"))
        min_members_per_task > 0 ||
            throw(InvalidConfiguration("min_members_per_task must be greater than zero"))
        return new(Int(max_tasks), Int(min_members_per_task))
    end
end

AcceleratedKernelsExecution(; max_tasks::Integer=Threads.nthreads(),
    min_members_per_task::Integer=1) =
    AcceleratedKernelsExecution(max_tasks, min_members_per_task)

"""
    DaggerExecution(; scope=nothing)

Opt-in Dagger.jl scheduling for coarse task graphs, multi-process runs, and
data-local ensemble execution. It is not a real-time inner-loop policy. The
optional `scope` is a Dagger processor scope or scope constraint.

Each `run_ensemble!` call fetches the updated members at its completion. A
distributed callable should therefore encompass a trajectory or sweep unit,
not one latency-sensitive simulation step.
"""
struct DaggerExecution{S} <: AbstractExecutionPolicy
    scope::S
end

DaggerExecution(; scope=nothing) = DaggerExecution{typeof(scope)}(scope)

struct ParallelConfig
    threads::Int
    fft_threads::Int
    blas_threads::Int
end

function with_parallel_config(cfg::ParallelConfig, f::Function)
    if Threads.nthreads() != cfg.threads
        @warn "start Julia with JULIA_NUM_THREADS=$(cfg.threads)" current=Threads.nthreads()
    end
    set_fft_provider_threads!(cfg.fft_threads)
    BLAS.set_num_threads(cfg.blas_threads)
    return f()
end
