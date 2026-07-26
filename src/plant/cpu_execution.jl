const _CPU_EXECUTION_COMPONENT = :cpu_execution

@noinline function _cpu_execution_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(
        _CPU_EXECUTION_COMPONENT,
        reason,
        String(message),
    ))
end

abstract type _AbstractCPUExecutionMode end
struct _DeterministicCPUExecutionMode <: _AbstractCPUExecutionMode end
struct _GroupedCPUExecutionMode <: _AbstractCPUExecutionMode end

struct _CPUExecutionBudgetToken end
const _CPU_EXECUTION_BUDGET_TOKEN = _CPUExecutionBudgetToken()

"""
Immutable execution-context declaration for prepared CPU path groups.

`cpu_context_count` is the complete CPU-context reservation admitted to this
execution domain. `julia_thread_count` is the required Julia default-pool
size, `outer_owner_count` is the maximum number of simultaneously executing
path-group owners, and `group_julia_thread_count` is the total Julia threads
one owner may use for group-local numerical work such as an explicitly
threaded KA/AcceleratedKernels policy. The FFT/BLAS counts are the
provider-reported total threads that one owner may use in the corresponding
numerical call.
"""
struct CPUExecutionBudget{M<:_AbstractCPUExecutionMode}
    cpu_context_count::Int
    julia_thread_count::Int
    outer_owner_count::Int
    group_julia_thread_count::Int
    fft_thread_count::Int
    blas_thread_count::Int
    mode::M

    function CPUExecutionBudget(
        cpu_context_count::Int,
        julia_thread_count::Int,
        outer_owner_count::Int,
        group_julia_thread_count::Int,
        fft_thread_count::Int,
        blas_thread_count::Int,
        mode::M,
        ::_CPUExecutionBudgetToken,
    ) where {M<:_AbstractCPUExecutionMode}
        return new{M}(
            cpu_context_count,
            julia_thread_count,
            outer_owner_count,
            group_julia_thread_count,
            fft_thread_count,
            blas_thread_count,
            mode,
        )
    end
end

"""
Observed execution environment used to validate a `CPUExecutionBudget`.

Construction records values only; it does not change Julia, FFT-provider, or
BLAS process-global settings. The caller supplies its admitted CPU-context
capacity and observed FFT-provider thread count because neither can be inferred
portably from the host or from a backend-neutral optical plan.
"""
struct CPUExecutionEnvironment
    available_cpu_context_count::Int
    julia_thread_count::Int
    fft_thread_count::Int
    blas_thread_count::Int
end

@inline function _cpu_execution_count(
    value::Integer,
    reason::Symbol,
    label::AbstractString,
)
    value > 0 || _cpu_execution_error(
        reason,
        "$label must be greater than zero",
    )
    value <= typemax(Int) || _cpu_execution_error(
        reason,
        "$label exceeds the supported Int range",
    )
    return Int(value)
end

@inline _cpu_execution_count(
    ::Bool,
    reason::Symbol,
    label::AbstractString,
) = _cpu_execution_error(
    reason,
    "$label must be an integer count, not Bool",
)

@inline function _require_cpu_execution_fit(
    outer_owner_count::Int,
    inner_thread_count::Int,
    cpu_context_count::Int,
    reason::Symbol,
    label::AbstractString,
)
    outer_owner_count <= div(cpu_context_count, inner_thread_count) ||
        _cpu_execution_error(
            reason,
            "$label oversubscribes the admitted CPU execution contexts",
        )
    return nothing
end

function _grouped_cpu_execution_budget(
    cpu_context_count::Int,
    julia_thread_count::Int,
    outer_owner_count::Int,
    group_julia_thread_count::Int,
    fft_thread_count::Int,
    blas_thread_count::Int,
)
    julia_thread_count <= cpu_context_count || _cpu_execution_error(
        :julia_thread_oversubscription,
        "Julia default-pool threads exceed the admitted CPU execution contexts",
    )
    outer_owner_count <= julia_thread_count || _cpu_execution_error(
        :outer_owner_oversubscription,
        "outer path-group owners exceed the declared Julia default-pool threads",
    )
    outer_owner_count <= div(
        julia_thread_count, group_julia_thread_count) ||
        _cpu_execution_error(
            :group_julia_thread_oversubscription,
            "concurrent group-local Julia threads exceed the declared Julia default pool",
        )
    _require_cpu_execution_fit(
        outer_owner_count,
        fft_thread_count,
        cpu_context_count,
        :fft_thread_oversubscription,
        "concurrent FFT threads",
    )
    _require_cpu_execution_fit(
        outer_owner_count,
        blas_thread_count,
        cpu_context_count,
        :blas_thread_oversubscription,
        "concurrent BLAS threads",
    )
    return CPUExecutionBudget(
        cpu_context_count,
        julia_thread_count,
        outer_owner_count,
        group_julia_thread_count,
        fft_thread_count,
        blas_thread_count,
        _GroupedCPUExecutionMode(),
        _CPU_EXECUTION_BUDGET_TOKEN,
    )
end

"""
    deterministic_cpu_execution_budget()

Declare the canonical deterministic CPU policy: one admitted CPU context, one
Julia thread, one serial path-group owner, one FFT-provider thread, and one
BLAS thread.
"""
function deterministic_cpu_execution_budget()
    return CPUExecutionBudget(
        1,
        1,
        1,
        1,
        1,
        1,
        _DeterministicCPUExecutionMode(),
        _CPU_EXECUTION_BUDGET_TOKEN,
    )
end

"""
    grouped_cpu_execution_budget(;
        cpu_context_count,
        julia_thread_count,
        outer_owner_count,
        group_julia_thread_count=1,
        fft_thread_count=1,
        blas_thread_count=1)

Declare a bounded coarse CPU path-group execution policy. Construction rejects
outer ownership or nested group-local Julia, FFT, and BLAS thread counts whose
worst-case simultaneous use exceeds the relevant Julia-pool or CPU-context
capacity.
"""
function grouped_cpu_execution_budget(;
    cpu_context_count::Integer,
    julia_thread_count::Integer,
    outer_owner_count::Integer,
    group_julia_thread_count::Integer=1,
    fft_thread_count::Integer=1,
    blas_thread_count::Integer=1,
)
    contexts = _cpu_execution_count(
        cpu_context_count,
        :invalid_cpu_context_count,
        "CPU execution-context count",
    )
    julia_threads = _cpu_execution_count(
        julia_thread_count,
        :invalid_julia_thread_count,
        "Julia default-pool thread count",
    )
    owners = _cpu_execution_count(
        outer_owner_count,
        :invalid_outer_owner_count,
        "outer path-group owner count",
    )
    group_julia_threads = _cpu_execution_count(
        group_julia_thread_count,
        :invalid_group_julia_thread_count,
        "group-local Julia thread count",
    )
    fft_threads = _cpu_execution_count(
        fft_thread_count,
        :invalid_fft_thread_count,
        "FFT-provider thread count",
    )
    blas_threads = _cpu_execution_count(
        blas_thread_count,
        :invalid_blas_thread_count,
        "BLAS thread count",
    )
    return _grouped_cpu_execution_budget(
        contexts,
        julia_threads,
        owners,
        group_julia_threads,
        fft_threads,
        blas_threads,
    )
end

"""
    CPUExecutionEnvironment(;
        available_cpu_context_count,
        fft_thread_count,
        julia_thread_count=Threads.nthreads(),
        blas_thread_count=BLAS.get_num_threads())

Record the observed resource boundary without changing any process-global
thread setting.
"""
function CPUExecutionEnvironment(;
    available_cpu_context_count::Integer,
    fft_thread_count::Integer,
    julia_thread_count::Integer=Threads.nthreads(),
    blas_thread_count::Integer=BLAS.get_num_threads(),
)
    return CPUExecutionEnvironment(
        _cpu_execution_count(
            available_cpu_context_count,
            :invalid_available_cpu_context_count,
            "available CPU execution-context count",
        ),
        _cpu_execution_count(
            julia_thread_count,
            :invalid_observed_julia_thread_count,
            "observed Julia default-pool thread count",
        ),
        _cpu_execution_count(
            fft_thread_count,
            :invalid_observed_fft_thread_count,
            "observed FFT-provider thread count",
        ),
        _cpu_execution_count(
            blas_thread_count,
            :invalid_observed_blas_thread_count,
            "observed BLAS thread count",
        ),
    )
end

@inline is_deterministic_cpu_execution(
    ::CPUExecutionBudget{_DeterministicCPUExecutionMode},
) = true
@inline is_deterministic_cpu_execution(
    ::CPUExecutionBudget{_GroupedCPUExecutionMode},
) = false

"""
    validate_cpu_execution_budget(budget, environment)

Validate the declared CPU budget against an explicitly observed environment.
This operation is mutation-free and returns `budget` on success.
"""
function validate_cpu_execution_budget(
    budget::CPUExecutionBudget,
    environment::CPUExecutionEnvironment,
)
    budget.cpu_context_count <=
        environment.available_cpu_context_count ||
        _cpu_execution_error(
            :cpu_context_capacity,
            "declared CPU execution contexts exceed the admitted environment capacity",
        )
    budget.julia_thread_count == environment.julia_thread_count ||
        _cpu_execution_error(
            :julia_thread_mismatch,
            "observed Julia default-pool thread count does not match the declared CPU budget",
        )
    budget.fft_thread_count == environment.fft_thread_count ||
        _cpu_execution_error(
            :fft_thread_mismatch,
            "observed FFT-provider thread count does not match the declared CPU budget",
        )
    budget.blas_thread_count == environment.blas_thread_count ||
        _cpu_execution_error(
            :blas_thread_mismatch,
            "observed BLAS thread count does not match the declared CPU budget",
        )
    return budget
end
