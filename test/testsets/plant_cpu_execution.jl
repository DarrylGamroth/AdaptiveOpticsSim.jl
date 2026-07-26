function cpu_execution_test_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_cpu_execution_error(f, reason::Symbol)
    error = cpu_execution_test_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component == :cpu_execution
        @test error.reason == reason
        @test !isempty(error.msg)
    end
    return nothing
end

@testset "Prepared CPU execution-budget contract" begin
    deterministic = @inferred deterministic_cpu_execution_budget()
    @test deterministic isa CPUExecutionBudget
    @test isbitstype(typeof(deterministic))
    @test is_deterministic_cpu_execution(deterministic)
    @test deterministic.cpu_context_count == 1
    @test deterministic.julia_thread_count == 1
    @test deterministic.outer_owner_count == 1
    @test deterministic.group_julia_thread_count == 1
    @test deterministic.fft_thread_count == 1
    @test deterministic.blas_thread_count == 1

    grouped = @inferred grouped_cpu_execution_budget(
        cpu_context_count=8,
        julia_thread_count=4,
        outer_owner_count=3,
        group_julia_thread_count=1,
        fft_thread_count=2,
        blas_thread_count=2,
    )
    @test grouped isa CPUExecutionBudget
    @test isbitstype(typeof(grouped))
    @test !is_deterministic_cpu_execution(grouped)
    @test grouped.cpu_context_count == 8
    @test grouped.julia_thread_count == 4
    @test grouped.outer_owner_count == 3
    @test grouped.group_julia_thread_count == 1
    @test grouped.fft_thread_count == 2
    @test grouped.blas_thread_count == 2

    inner_parallel = grouped_cpu_execution_budget(
        cpu_context_count=8,
        julia_thread_count=1,
        outer_owner_count=1,
        fft_thread_count=8,
        blas_thread_count=8,
    )
    @test inner_parallel.outer_owner_count == 1
    @test inner_parallel.fft_thread_count == 8
    @test inner_parallel.blas_thread_count == 8

    julia_inner_parallel = grouped_cpu_execution_budget(
        cpu_context_count=8,
        julia_thread_count=8,
        outer_owner_count=1,
        group_julia_thread_count=8,
    )
    @test julia_inner_parallel.outer_owner_count == 1
    @test julia_inner_parallel.group_julia_thread_count == 8

    invalid_declarations = (
        (
            :invalid_cpu_context_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=0,
                julia_thread_count=1,
                outer_owner_count=1,
            ),
        ),
        (
            :invalid_julia_thread_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=1,
                julia_thread_count=false,
                outer_owner_count=1,
            ),
        ),
        (
            :invalid_outer_owner_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=1,
                julia_thread_count=1,
                outer_owner_count=big(typemax(Int)) + 1,
            ),
        ),
        (
            :invalid_fft_thread_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=1,
                julia_thread_count=1,
                outer_owner_count=1,
                fft_thread_count=0,
            ),
        ),
        (
            :invalid_group_julia_thread_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=1,
                julia_thread_count=1,
                outer_owner_count=1,
                group_julia_thread_count=0,
            ),
        ),
        (
            :invalid_blas_thread_count,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=1,
                julia_thread_count=1,
                outer_owner_count=1,
                blas_thread_count=-1,
            ),
        ),
        (
            :julia_thread_oversubscription,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=3,
                julia_thread_count=4,
                outer_owner_count=3,
            ),
        ),
        (
            :outer_owner_oversubscription,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=4,
                julia_thread_count=2,
                outer_owner_count=3,
            ),
        ),
        (
            :fft_thread_oversubscription,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=4,
                julia_thread_count=4,
                outer_owner_count=3,
                fft_thread_count=2,
            ),
        ),
        (
            :group_julia_thread_oversubscription,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=8,
                julia_thread_count=4,
                outer_owner_count=3,
                group_julia_thread_count=2,
            ),
        ),
        (
            :blas_thread_oversubscription,
            () -> grouped_cpu_execution_budget(
                cpu_context_count=4,
                julia_thread_count=4,
                outer_owner_count=3,
                blas_thread_count=2,
            ),
        ),
    )
    for (reason, operation) in invalid_declarations
        assert_cpu_execution_error(operation, reason)
    end

    for name in (
        :CPUExecutionBudget,
        :CPUExecutionEnvironment,
        :deterministic_cpu_execution_budget,
        :grouped_cpu_execution_budget,
        :is_deterministic_cpu_execution,
        :validate_cpu_execution_budget,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
end

@testset "CPU execution environment validation is mutation-free" begin
    julia_threads_before = Threads.nthreads()
    blas_threads_before = BLAS.get_num_threads()
    environment = @inferred CPUExecutionEnvironment(
        available_cpu_context_count=8,
        julia_thread_count=4,
        fft_thread_count=2,
        blas_thread_count=2,
    )
    @test isbitstype(typeof(environment))
    @test environment.available_cpu_context_count == 8
    @test environment.julia_thread_count == 4
    @test environment.fft_thread_count == 2
    @test environment.blas_thread_count == 2

    grouped = grouped_cpu_execution_budget(
        cpu_context_count=8,
        julia_thread_count=4,
        outer_owner_count=3,
        fft_thread_count=2,
        blas_thread_count=2,
    )
    @test @inferred(validate_cpu_execution_budget(
        grouped, environment)) === grouped

    mismatches = (
        (
            :cpu_context_capacity,
            CPUExecutionEnvironment(
                available_cpu_context_count=7,
                julia_thread_count=4,
                fft_thread_count=2,
                blas_thread_count=2,
            ),
        ),
        (
            :julia_thread_mismatch,
            CPUExecutionEnvironment(
                available_cpu_context_count=8,
                julia_thread_count=3,
                fft_thread_count=2,
                blas_thread_count=2,
            ),
        ),
        (
            :fft_thread_mismatch,
            CPUExecutionEnvironment(
                available_cpu_context_count=8,
                julia_thread_count=4,
                fft_thread_count=1,
                blas_thread_count=2,
            ),
        ),
        (
            :blas_thread_mismatch,
            CPUExecutionEnvironment(
                available_cpu_context_count=8,
                julia_thread_count=4,
                fft_thread_count=2,
                blas_thread_count=1,
            ),
        ),
    )
    for (reason, observed) in mismatches
        assert_cpu_execution_error(
            () -> validate_cpu_execution_budget(grouped, observed),
            reason,
        )
    end

    actual = CPUExecutionEnvironment(
        available_cpu_context_count=max(1, Threads.nthreads()),
        fft_thread_count=1,
    )
    if Threads.nthreads() == 1 && BLAS.get_num_threads() == 1
        deterministic = deterministic_cpu_execution_budget()
        @test validate_cpu_execution_budget(
            deterministic, actual) === deterministic
    else
        assert_cpu_execution_error(
            () -> validate_cpu_execution_budget(
                deterministic_cpu_execution_budget(), actual),
            Threads.nthreads() == 1 ?
                :blas_thread_mismatch : :julia_thread_mismatch,
        )
    end
    @test Threads.nthreads() == julia_threads_before
    @test BLAS.get_num_threads() == blas_threads_before

    invalid_environment = (
        (
            :invalid_available_cpu_context_count,
            () -> CPUExecutionEnvironment(
                available_cpu_context_count=false,
                fft_thread_count=1,
            ),
        ),
        (
            :invalid_observed_julia_thread_count,
            () -> CPUExecutionEnvironment(
                available_cpu_context_count=1,
                julia_thread_count=0,
                fft_thread_count=1,
            ),
        ),
        (
            :invalid_observed_fft_thread_count,
            () -> CPUExecutionEnvironment(
                available_cpu_context_count=1,
                fft_thread_count=0,
            ),
        ),
        (
            :invalid_observed_blas_thread_count,
            () -> CPUExecutionEnvironment(
                available_cpu_context_count=1,
                fft_thread_count=1,
                blas_thread_count=0,
            ),
        ),
    )
    for (reason, operation) in invalid_environment
        assert_cpu_execution_error(operation, reason)
    end
end
