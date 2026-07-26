include(joinpath(
    @__DIR__,
    "..",
    "..",
    "benchmarks",
    "support",
    "gate5_optical_placement.jl",
))

const Gate5Closure = Gate5OpticalPlacementBenchmark

function gate5_quick_workload()
    contract = TOML.parsefile(joinpath(
        @__DIR__,
        "..",
        "..",
        "benchmarks",
        "contracts",
        "gate5_optical_placement.toml",
    ))
    workload = deepcopy(contract["workload"])
    workload["path_count"] = 4
    workload["resolution"] = 8
    return contract, workload
end

@testset "Gate 5 integrated optical-placement closure" begin
    contract, workload = gate5_quick_workload()

    oracle = Gate5Closure.validate_integrated_oracle(workload)
    @test oracle["passed"]
    @test oracle["path_count"] == 4
    @test oracle["wfs_path_count"] == 2
    @test oracle["science_path_count"] == 2
    @test oracle["maximum_opd_error_m"] <= oracle["tolerance"]
    @test oracle["declaration_order_opd_error_m"] <= oracle["tolerance"]
    @test oracle["declaration_order_output_error"] <=
        oracle["output_order_tolerance"]
    @test length(oracle["path_hashes"]) == 4

    finite_support = Gate5Closure.validate_finite_support(workload)
    @test finite_support["passed"]
    @test finite_support["maximum_error_m"] <=
        finite_support["tolerance"]
    @test finite_support["supported_samples"] == 9
    @test finite_support["zeroed_samples"] == 16

    fixed_storage =
        Gate5Closure.validate_fixed_storage(workload, 4, 12)
    @test fixed_storage["storage_signature_stable"]
    @test fixed_storage["deterministic_product_hashes_stable"]
    @test fixed_storage["final_cycle"] == 16

    topology = only(Gate5Closure.topology_characterization(
        workload,
        (4,),
        8,
    ))
    @test topology["counts_passed"]
    @test topology["path_count"] == 4
    @test topology["controllable_optic_count"] == 5
    @test topology["sampled_aberration_count"] == 3
    @test topology["controllable_binding_count"] ==
        topology["expected_controllable_binding_count"] == 14
    @test topology["sampled_aberration_binding_count"] ==
        topology["expected_sampled_aberration_binding_count"] == 6
    @test topology["first_cycle_completed"]

    invalid_altitude = deepcopy(workload)
    invalid_altitude["common_dm_altitudes_m"][1] = 1.0
    @test_throws ErrorException Gate5Closure.gate5_plant_definition(
        invalid_altitude,
    )

    operation, expected, _ =
        Gate5Closure.prepare_gate5_operation(workload)
    Gate5Closure.run_gate5_cycles!(operation, 20)
    @test Gate5Closure.numerical_oracle(
        operation,
        expected,
    ).maximum_error <= Float64(workload["numerical_atol"])
    allocation_cycles = 20
    if coverage_instrumented()
        @test_skip "Gate 5 allocation gate disabled under coverage instrumentation"
    else
        GC.gc()
        allocated = @allocated Gate5Closure.run_gate5_cycles!(
            operation,
            allocation_cycles,
        )
        @test allocated <=
            Int(contract["max_alloc_bytes_per_cycle"]) * allocation_cycles
    end
end
