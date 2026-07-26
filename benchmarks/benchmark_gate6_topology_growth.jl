using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant
using Dates
using HdrHistogram
using InteractiveUtils
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "gate5_optical_placement.jl"))
include(joinpath(@__DIR__, "support", "hdr_histogram_artifact.jl"))

const AOSPlant = AdaptiveOpticsSim.Plant
const Gate5Topology = Gate5OpticalPlacementBenchmark
const HdrArtifact = HdrHistogramArtifact
const GATE6_TOPOLOGY_CONTRACT_PATH = get(
    ENV,
    "AOS_GATE6_TOPOLOGY_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate6_topology_growth.toml"),
)
const GATE6_TOPOLOGY_OUTPUT_PATH =
    get(ENV, "AOS_GATE6_TOPOLOGY_OUTPUT", "")

function configure_gate6_topology_probe!()
    Threads.nthreads() == 1 ||
        error("Gate 6 topology-growth evidence requires one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

@inline gate6_topology_sha(value::AbstractString) =
    bytes2hex(SHA.sha256(codeunits(value)))

function gate6_topology_type_fingerprint(value)
    display = string(typeof(value))
    return Dict{String,Any}(
        "sha256" => gate6_topology_sha(display),
        "display_bytes" => ncodeunits(display),
    )
end

function gate6_topology_native_code_proxy(function_value, argument_types)
    io = IOBuffer()
    code_native(
        io,
        function_value,
        argument_types;
        debuginfo=:none,
        binary=false,
    )
    code = String(take!(io))
    return Dict{String,Any}(
        "text_bytes" => ncodeunits(code),
        "sha256" => gate6_topology_sha(code),
    )
end

function gate6_topology_method_instance_count(
    function_value, argument_types)
    method = which(function_value, argument_types)
    return count(
        specialization -> specialization isa Core.MethodInstance,
        Base.specializations(method),
    )
end

@inline function gate6_topology_control_statement(statement)
    return statement isa Core.ReturnNode ||
        statement isa Core.GotoNode ||
        statement isa Core.GotoIfNot ||
        (statement isa Expr && statement.head in (
            :gc_preserve_begin,
            :gc_preserve_end,
            :boundscheck,
            :inbounds,
        ))
end

function gate6_topology_any_value_count(
    function_value, argument_types)
    typed = only(Base.code_typed(
        function_value, argument_types; optimize=true))
    code = typed.first
    count_any = count(==(Any), code.slottypes)
    for index in eachindex(code.code)
        code.ssavaluetypes[index] === Any || continue
        gate6_topology_control_statement(code.code[index]) && continue
        count_any += 1
    end
    return count_any
end

function gate6_topology_all_acquisition_ids(plant)
    return Symbol[
        AOSPlant.acquisition_id(definition).name
        for definition in AOSPlant.acquisition_definitions(
            plant.definition)
    ]
end

function gate6_topology_histogram_summary(
    histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64,
    wall_end_ns::UInt64,
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    encoded = HdrArtifact.verified_sparse_histogram(
        histogram,
        lowest_ns,
        highest_ns,
        significant_figures,
        samples,
    )
    wall_ns = Int64(wall_end_ns - wall_start_ns)
    return Dict{String,Any}(
        "samples" => samples,
        "completed_operations" => samples,
        "failed_operations" => 0,
        "monotonic_start_ns" => Int64(wall_start_ns),
        "monotonic_end_ns" => Int64(wall_end_ns),
        "wall_ns" => wall_ns,
        "throughput_hz" => 1.0e9 * samples / wall_ns,
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" =>
            HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" =>
            HdrHistogram.value_at_percentile(histogram, 90.0),
        "p99_ns" =>
            HdrHistogram.value_at_percentile(histogram, 99.0),
        "p99_9_ns" =>
            HdrHistogram.value_at_percentile(histogram, 99.9),
        "max_ns" => max(histogram),
        encoded...,
    )
end

function measure_gate6_topology_service!(
    operation::Gate5Topology.Gate5OpticalOperation,
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns, highest_ns, significant_figures)
    GC.gc()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        operation()
        HdrHistogram.record_value!(
            histogram, Int64(time_ns() - start))
    end
    wall_end = time_ns()
    return gate6_topology_histogram_summary(
        histogram,
        wall_start,
        wall_end,
        samples,
        lowest_ns,
        highest_ns,
        significant_figures,
    )
end

function gate6_topology_registry_fingerprints(
    operation::Gate5Topology.Gate5OpticalOperation,
    selection,
)
    plant = operation.plant
    return Dict{String,Any}(
        "plant_definition" =>
            gate6_topology_type_fingerprint(plant.definition),
        "prepared_plant" =>
            gate6_topology_type_fingerprint(plant),
        "prepared_rngs" =>
            gate6_topology_type_fingerprint(plant.rngs),
        "controllable_optic_registry" =>
            gate6_topology_type_fingerprint(plant.controllable_optics),
        "command_endpoint_registry" =>
            gate6_topology_type_fingerprint(plant.command_endpoints),
        "path_registry" =>
            gate6_topology_type_fingerprint(plant.paths),
        "acquisition_registry" =>
            gate6_topology_type_fingerprint(plant.acquisitions),
        "selection" =>
            gate6_topology_type_fingerprint(selection),
        "selection_path_registry" =>
            gate6_topology_type_fingerprint(selection.paths),
        "selection_acquisition_registry" =>
            gate6_topology_type_fingerprint(selection.acquisitions),
        "event_loop" =>
            gate6_topology_type_fingerprint(operation.event_loop),
        "event_state" =>
            gate6_topology_type_fingerprint(operation.state),
        "event_workspace" =>
            gate6_topology_type_fingerprint(operation.workspace),
    )
end

function gate6_topology_code_proxies(
    operation::Gate5Topology.Gate5OpticalOperation,
    selection,
)
    epoch = current_epoch(
        AOSPlant.plant_atmosphere(operation.plant.definition))
    event_types = Tuple{
        typeof(operation.event_loop),
        typeof(operation.state),
        typeof(operation.workspace),
    }
    selection_types = Tuple{typeof(selection),typeof(epoch)}

    path = first(AOSPlant.prepared_paths(operation.plant))
    acquisition =
        first(AOSPlant.prepared_acquisitions(operation.plant))
    acquisition_execution =
        acquisition.provider.implementation.execution
    acquisition_products = acquisition.provider.products
    acquisition_rng = AOSPlant.rng_stream_state(
        first(operation.plant.rngs.acquisitions), Val(:detector))
    path_kernel_types = Tuple{
        typeof(path.result),
        typeof(path.input),
        typeof(path.execution),
    }
    acquisition_kernel_types = Tuple{
        typeof(acquisition_products),
        typeof(acquisition.path_result),
        typeof(acquisition_execution),
        typeof(acquisition_rng),
    }

    event_native = gate6_topology_native_code_proxy(
        AOSPlant.step_plant_events!, event_types)
    selection_native = gate6_topology_native_code_proxy(
        AOSPlant.execute_acquisition_selection!, selection_types)
    return Dict{String,Any}(
        "event_entry" => Dict{String,Any}(
            "native_code" => event_native,
            "method_instance_count" =>
                gate6_topology_method_instance_count(
                    AOSPlant.step_plant_events!, event_types),
        ),
        "selection_entry" => Dict{String,Any}(
            "native_code" => selection_native,
            "method_instance_count" =>
                gate6_topology_method_instance_count(
                    AOSPlant.execute_acquisition_selection!,
                    selection_types),
        ),
        "stage_kernels" => Dict{String,Any}(
            "path_any_slots" => gate6_topology_any_value_count(
                AOSPlant.execute_path!, path_kernel_types),
            "acquisition_any_slots" =>
                gate6_topology_any_value_count(
                    AOSPlant.execute_acquisition!,
                    acquisition_kernel_types),
        ),
    )
end

function gate6_topology_counts(operation, profile)
    plant = operation.plant
    counts = Dict{String,Any}(
        "path_count" => length(plant.paths),
        "acquisition_count" => length(plant.acquisitions),
        "controllable_optic_count" =>
            length(plant.controllable_optics),
        "command_endpoint_count" =>
            length(plant.command_endpoints),
        "sampled_aberration_count" =>
            length(plant.sampled_aberrations),
    )
    passed = all(
        counts[key] == Int(profile[key])
        for key in keys(counts)
    )
    counts["passed"] = passed
    return counts
end

function run_gate6_topology_child(
    contract_path::AbstractString,
    label::AbstractString,
    path_count::Int,
    first_use_mode::AbstractString,
    repetition::Int,
    samples::Int,
    warmup_cycles::Int,
    allocation_cycles::Int,
    output_path::AbstractString,
)
    configure_gate6_topology_probe!()
    contract = TOML.parsefile(contract_path)
    gate5_contract = TOML.parsefile(joinpath(
        dirname(contract_path),
        String(contract["gate5_workload_contract"]),
    ))
    workload = deepcopy(gate5_contract["workload"])
    profile = only(filter(contract["topologies"]) do candidate
        String(candidate["label"]) == label
    end)

    GC.gc()
    timed = @timed begin
        operation, expected, _ =
            Gate5Topology.prepare_gate5_operation(
                workload;
                path_count,
                resolution=Int(contract["topology_resolution"]),
            )
        selection = AOSPlant.prepare_acquisition_selection(
            operation.plant,
            gate6_topology_all_acquisition_ids(operation.plant),
        )
        (; operation, expected, selection)
    end
    operation = timed.value.operation
    expected = timed.value.expected
    selection = timed.value.selection

    first_start = time_ns()
    if first_use_mode == "event"
        operation()
    elseif first_use_mode == "selection"
        AOSPlant.execute_acquisition_selection_at!(selection, 0.0)
    else
        error("unknown Gate 6 topology first-use mode $first_use_mode")
    end
    first_invocation_ns = Int64(time_ns() - first_start)

    allocation = Dict{String,Any}()
    service = Dict{String,Any}()
    numerical_oracle = Dict{String,Any}()
    if first_use_mode == "event"
        Gate5Topology.run_gate5_cycles!(operation, warmup_cycles)
        Gate5Topology.run_gate5_cycles!(operation, allocation_cycles)
        GC.gc()
        allocated = Int64(@allocated Gate5Topology.run_gate5_cycles!(
            operation, allocation_cycles))
        allocation = Dict{String,Any}(
            "cycles" => allocation_cycles,
            "window_bytes" => allocated,
            "bytes_per_cycle" => allocated / allocation_cycles,
        )
        oracle = Gate5Topology.numerical_oracle(operation, expected)
        numerical_oracle = Dict{String,Any}(
            "maximum_error_m" => oracle.maximum_error,
            "wfs_maximum_error_m" => oracle.wfs_maximum_error,
            "science_maximum_error_m" => oracle.science_maximum_error,
            "passed" => oracle.maximum_error <=
                Float64(workload["numerical_atol"]),
        )
        service = measure_gate6_topology_service!(
            operation,
            samples,
            Int64(contract["histogram_lowest_ns"]),
            Int64(contract["histogram_highest_ns"]),
            Int(contract["histogram_significant_figures"]),
        )
    end

    result = Dict{String,Any}(
        "label" => label,
        "path_count" => path_count,
        "first_use_mode" => first_use_mode,
        "repetition" => repetition,
        "preparation" => Dict{String,Any}(
            "elapsed_ns" => round(Int64, timed.time * 1.0e9),
            "allocated_bytes" => Int64(timed.bytes),
            "gc_ns" => round(Int64, timed.gctime * 1.0e9),
        ),
        "first_invocation_ns" => first_invocation_ns,
        "prepared_summary_bytes" => Base.summarysize((
            operation.plant,
            operation.event_loop,
            operation.state,
            operation.workspace,
            selection,
        )),
        "counts" => gate6_topology_counts(operation, profile),
        "type_fingerprints" =>
            gate6_topology_registry_fingerprints(operation, selection),
        "code_proxies" =>
            gate6_topology_code_proxies(operation, selection),
        "allocation" => allocation,
        "service" => service,
        "numerical_oracle" => numerical_oracle,
    )
    open(output_path, "w") do io
        TOML.print(io, result; sorted=true)
    end
    return result
end

gate6_topology_median_integer(values) =
    round(Int64, median(collect(values)))

function gate6_topology_profile_summary(probes, profile)
    label = String(profile["label"])
    selected = filter(probe -> probe["label"] == label, probes)
    event_probes = filter(
        probe -> probe["first_use_mode"] == "event", selected)
    selection_probes = filter(
        probe -> probe["first_use_mode"] == "selection", selected)
    storage_values =
        Int64[probe["prepared_summary_bytes"] for probe in selected]
    preparation_times = Int64[
        probe["preparation"]["elapsed_ns"] for probe in selected]
    preparation_allocations = Int64[
        probe["preparation"]["allocated_bytes"] for probe in selected]
    event_first = Int64[
        probe["first_invocation_ns"] for probe in event_probes]
    selection_first = Int64[
        probe["first_invocation_ns"] for probe in selection_probes]
    event_allocations = Float64[
        probe["allocation"]["bytes_per_cycle"]
        for probe in event_probes
    ]
    service_p50 = Int64[
        probe["service"]["p50_ns"] for probe in event_probes]
    service_p99 = Int64[
        probe["service"]["p99_ns"] for probe in event_probes]
    return Dict{String,Any}(
        "label" => label,
        "path_count" => Int(profile["path_count"]),
        "probe_count" => length(selected),
        "median_preparation_elapsed_ns" =>
            gate6_topology_median_integer(preparation_times),
        "maximum_preparation_elapsed_ns" =>
            maximum(preparation_times),
        "median_preparation_allocated_bytes" =>
            gate6_topology_median_integer(preparation_allocations),
        "maximum_preparation_allocated_bytes" =>
            maximum(preparation_allocations),
        "median_prepared_summary_bytes" =>
            gate6_topology_median_integer(storage_values),
        "maximum_prepared_summary_bytes" => maximum(storage_values),
        "median_first_event_invocation_ns" =>
            gate6_topology_median_integer(event_first),
        "maximum_first_event_invocation_ns" => maximum(event_first),
        "median_first_selection_invocation_ns" =>
            gate6_topology_median_integer(selection_first),
        "maximum_first_selection_invocation_ns" =>
            maximum(selection_first),
        "maximum_warmed_allocation_bytes_per_cycle" =>
            maximum(event_allocations),
        "maximum_warmed_allocation_bytes_per_path_per_cycle" =>
            maximum(event_allocations) / Int(profile["path_count"]),
        "median_service_p50_ns" =>
            gate6_topology_median_integer(service_p50),
        "median_service_p99_ns" =>
            gate6_topology_median_integer(service_p99),
    )
end

function gate6_topology_distinct_fingerprints(probes)
    first_fingerprints = first(probes)["type_fingerprints"]
    return Dict{String,Any}(
        component => length(unique(
            probe["type_fingerprints"][component]["sha256"]
            for probe in probes
        ))
        for component in keys(first_fingerprints)
    )
end

function gate6_topology_code_growth(probes)
    result = Dict{String,Any}()
    for entry in ("event_entry", "selection_entry")
        sizes = Int64[
            probe["code_proxies"][entry]["native_code"]["text_bytes"]
            for probe in probes
        ]
        instances = Int[
            probe["code_proxies"][entry]["method_instance_count"]
            for probe in probes
        ]
        result[entry] = Dict{String,Any}(
            "minimum_native_code_text_bytes" => minimum(sizes),
            "maximum_native_code_text_bytes" => maximum(sizes),
            "native_code_size_ratio" =>
                maximum(sizes) / minimum(sizes),
            "minimum_method_instance_count" => minimum(instances),
            "maximum_method_instance_count" => maximum(instances),
            "method_instance_count_growth" =>
                maximum(instances) - minimum(instances),
            "distinct_native_code_hashes" => length(unique(
                probe["code_proxies"][entry]["native_code"]["sha256"]
                for probe in probes
            )),
        )
    end
    return result
end

function gate6_topology_command_output(command)
    try
        return readchomp(command)
    catch
        return "unknown"
    end
end

function gate6_topology_environment()
    cpu = first(Sys.cpu_info())
    status = gate6_topology_command_output(`git status --porcelain=v1`)
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "git_commit" =>
            gate6_topology_command_output(`git rev-parse HEAD`),
        "git_dirty" => !isempty(status) && status != "unknown",
        "git_status_porcelain" => status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" =>
            string(Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" =>
            string(Base.pkgversion(HdrHistogram)),
        "active_project" => something(
            Base.active_project(), "unknown"),
        "architecture" => string(Sys.ARCH),
        "cpu_target" => string(Sys.CPU_NAME),
        "cpu_model" => cpu.model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "fft_threads" => 1,
        "command" => "julia --threads=1 --project=benchmarks " *
            "benchmarks/benchmark_gate6_topology_growth.jl",
    )
end

function write_gate6_topology_artifact(
    path::AbstractString, artifact)
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return path
end

function run_gate6_topology_growth_benchmark()
    configure_gate6_topology_probe!()
    contract = TOML.parsefile(GATE6_TOPOLOGY_CONTRACT_PATH)
    profiles = contract["topologies"]
    repetitions =
        Int(contract["fresh_process_repetitions_per_first_use_mode"])
    samples = Int(contract["service_samples_per_run"])
    warmup_cycles = Int(contract["warmup_cycles"])
    allocation_cycles = Int(contract["allocation_cycles"])
    quick = get(ENV, "AOS_GATE6_TOPOLOGY_QUICK", "0") == "1"
    if quick
        isempty(GATE6_TOPOLOGY_OUTPUT_PATH) || error(
            "refusing durable topology-growth output in quick mode")
        repetitions = 1
        samples = min(samples, 100)
        warmup_cycles = min(warmup_cycles, 10)
        allocation_cycles = min(allocation_cycles, 10)
    end
    if !quick
        repetitions == Int(contract["service_runs"]) ||
            error("fresh event probes must match configured service runs")
    end

    probes = Dict{String,Any}[]
    mktempdir() do directory
        for profile in profiles
            label = String(profile["label"])
            path_count = Int(profile["path_count"])
            for first_use_mode in ("event", "selection")
                for repetition in 1:repetitions
                    output_path = joinpath(
                        directory,
                        "$(label)-$(first_use_mode)-$(repetition).toml",
                    )
                    command = `$(Base.julia_cmd()) --threads=1 --startup-file=no --project=$(@__DIR__) $(@__FILE__) --child $(abspath(GATE6_TOPOLOGY_CONTRACT_PATH)) $label $path_count $first_use_mode $repetition $samples $warmup_cycles $allocation_cycles $output_path`
                    run(command)
                    push!(probes, TOML.parsefile(output_path))
                end
            end
        end
    end

    summaries = Dict{String,Any}[
        gate6_topology_profile_summary(probes, profile)
        for profile in profiles
    ]
    smallest = first(summaries)
    largest = last(summaries)
    preparation_time_ratio =
        largest["median_preparation_elapsed_ns"] /
        smallest["median_preparation_elapsed_ns"]
    preparation_allocation_ratio =
        largest["median_preparation_allocated_bytes"] /
        smallest["median_preparation_allocated_bytes"]
    storage_ratio =
        largest["median_prepared_summary_bytes"] /
        smallest["median_prepared_summary_bytes"]
    distinct_fingerprints =
        gate6_topology_distinct_fingerprints(probes)
    code_growth = gate6_topology_code_growth(probes)
    maximum_any_slots = maximum(
        max(
            probe["code_proxies"]["stage_kernels"]["path_any_slots"],
            probe["code_proxies"]["stage_kernels"][
                "acquisition_any_slots"],
        )
        for probe in probes
    )

    gates = Dict{String,Any}(
        "exact_topology_counts" =>
            all(probe["counts"]["passed"] for probe in probes),
        "numerical_oracle" => all(
            probe["first_use_mode"] != "event" ||
                probe["numerical_oracle"]["passed"]
            for probe in probes
        ),
        "preparation_elapsed" => all(
            probe["preparation"]["elapsed_ns"] <=
                Int64(contract["max_preparation_elapsed_ns"])
            for probe in probes
        ),
        "preparation_allocated" => all(
            probe["preparation"]["allocated_bytes"] <=
                Int64(contract["max_preparation_allocated_bytes"])
            for probe in probes
        ),
        "first_event_invocation" => all(
            summary["maximum_first_event_invocation_ns"] <=
                Int64(contract["max_first_event_invocation_ns"])
            for summary in summaries
        ),
        "first_selection_invocation" => all(
            summary["maximum_first_selection_invocation_ns"] <=
                Int64(contract["max_first_selection_invocation_ns"])
            for summary in summaries
        ),
        "prepared_summary" => all(
            probe["prepared_summary_bytes"] <=
                Int64(contract["max_prepared_summary_bytes"])
            for probe in probes
        ),
        "preparation_time_growth" =>
            preparation_time_ratio <=
                Float64(contract[
                    "max_largest_to_small_preparation_time_ratio"]),
        "preparation_allocation_growth" =>
            preparation_allocation_ratio <=
                Float64(contract[
                    "max_largest_to_small_preparation_allocation_ratio"]),
        "prepared_storage_growth" =>
            storage_ratio <= Float64(contract[
                "max_largest_to_small_prepared_storage_ratio"]),
        "warmed_allocation" => all(
            summary[
                "maximum_warmed_allocation_bytes_per_path_per_cycle"] <=
                Int64(contract[
                    "max_warmed_alloc_bytes_per_path_per_cycle"])
            for summary in summaries
        ),
        "registry_type_cardinality" => all(
            count <= Int(contract[
                "max_distinct_registry_type_fingerprints"])
            for count in values(distinct_fingerprints)
        ),
        "native_code_growth" => all(
            code_growth[entry]["native_code_size_ratio"] <=
                Float64(contract[
                    "max_entry_native_code_size_ratio"])
            for entry in keys(code_growth)
        ),
        "method_instance_growth" => all(
            code_growth[entry]["method_instance_count_growth"] <=
                Int(contract[
                    "max_entry_method_instance_count_growth"])
            for entry in keys(code_growth)
        ),
        "stage_kernel_inference" =>
            maximum_any_slots <=
                Int(contract["max_stage_kernel_any_slots"]),
    )
    all_gates_passed = all(values(gates))
    environment = gate6_topology_environment()
    if !isempty(GATE6_TOPOLOGY_OUTPUT_PATH) &&
        environment["git_dirty"]
        error("refusing durable Gate 6 topology evidence from a dirty tree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(
            abspath(GATE6_TOPOLOGY_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" =>
            environment["git_commit"],
        "configured_fresh_process_repetitions_per_first_use_mode" =>
            repetitions,
        "configured_service_samples_per_run" => samples,
        "configured_service_runs" => repetitions,
        "quick_mode" => quick,
        "contract" => contract["contract"],
        "scope_exclusions" => contract["scope_exclusions"],
        "topologies" => profiles,
        "histogram" => Dict{String,Any}(
            "lowest_ns" =>
                Int64(contract["histogram_lowest_ns"]),
            "highest_ns" =>
                Int64(contract["histogram_highest_ns"]),
            "significant_figures" =>
                Int(contract["histogram_significant_figures"]),
            "p99_claim_supported" => false,
            "interpretation" =>
                "diagnostic closed-loop service-cost distributions",
        ),
        "environment" => environment,
        "probes" => probes,
        "profile_summaries" => summaries,
        "growth" => Dict{String,Any}(
            "largest_to_small_preparation_time_ratio" =>
                preparation_time_ratio,
            "largest_to_small_preparation_allocation_ratio" =>
                preparation_allocation_ratio,
            "largest_to_small_prepared_storage_ratio" =>
                storage_ratio,
        ),
        "distinct_registry_type_fingerprints" =>
            distinct_fingerprints,
        "code_growth" => code_growth,
        "maximum_stage_kernel_any_slots" => maximum_any_slots,
        "gates" => gates,
        "all_gates_passed" => all_gates_passed,
    )
    for summary in summaries
        println(
            summary["label"],
            ": prep=",
            summary["median_preparation_elapsed_ns"],
            " ns storage=",
            summary["median_prepared_summary_bytes"],
            " bytes p50=",
            summary["median_service_p50_ns"],
            " ns",
        )
    end
    println("all_gates_passed=", all_gates_passed)
    all_gates_passed ||
        error("Gate 6 topology-growth acceptance gate failed")
    write_gate6_topology_artifact(
        GATE6_TOPOLOGY_OUTPUT_PATH, artifact)
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isempty(ARGS) && first(ARGS) == "--child"
        length(ARGS) == 10 ||
            error("Gate 6 topology child expects nine arguments")
        run_gate6_topology_child(
            ARGS[2],
            ARGS[3],
            parse(Int, ARGS[4]),
            ARGS[5],
            parse(Int, ARGS[6]),
            parse(Int, ARGS[7]),
            parse(Int, ARGS[8]),
            parse(Int, ARGS[9]),
            ARGS[10],
        )
    else
        run_gate6_topology_growth_benchmark()
    end
end
