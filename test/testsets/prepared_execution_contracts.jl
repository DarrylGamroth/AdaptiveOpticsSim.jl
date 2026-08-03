const PE_CONTRACT_ROOT = joinpath(dirname(@__DIR__), "contracts")
const PE_REPOSITORY_ROOT = dirname(dirname(@__DIR__))

@testset "PE-00 source scanner semantics" begin
    @test pe_normalize_repository_path(raw"src\plant\preparation.jl") ==
        "src/plant/preparation.jl"
    @test pe_normalize_repository_path("src/plant/preparation.jl") ==
        "src/plant/preparation.jl"

    function scan_snippet(source)
        expression = Meta.parseall(source; filename="synthetic.jl")
        memory_sites = PEDebtSite[]
        any_sites = PEDebtSite[]
        _pe_scan_expression!(memory_sites, any_sites, expression,
            "synthetic.jl", 1, "<top-level>", "<expression>")
        return (memory=memory_sites, stored_any=any_sites)
    end

    permitted = scan_snippet("""
        wildcard(::Any) = nothing
        bounded(x::Tuple{P,Vararg{Any,N}}) where {P,N} = x
        unconstrained(x::T) where {T<:Any} = x
        invoke(wildcard, Tuple{Any}, 1)
    """)
    @test isempty(permitted.memory)
    @test isempty(permitted.stored_any)

    prohibited = scan_snippet("""
        struct BadOwner
            untyped
            scalar::Any
            values::Vector{Any}
            fixed::Memory{Any}
            qualified::Base.Memory{Int}
            erased_fixed::FixedSizeArray{Any}
            const qualified_const::Base.Memory{Int}
            const erased_const::Vector{Any} = Any[]
        end
        function prepare_bad(values::Vector{Any})::Any
            scratch = Any[]
            fixed = Memory{Any}(undef, 1)
            qualified = Base.Memory{Int}(undef, 1)
            erased_fixed = FixedSizeArray{Any}(undef, 1)
            mapped = Dict{String,Any}()
            return values, scratch, fixed, qualified, erased_fixed, mapped
        end
        function return_bad()::Vector{Any}
            implicit_array = []
            implicit_dict = Dict()
            implicit_id_dict = IdDict()
            implicit_weak_dict = WeakKeyDict()
            implicit_set = Set()
            implicit_vector = Vector()
            implicit_sized_vector = Vector(undef, 1)
            implicit_matrix = Matrix(undef, 1, 1)
            return implicit_array, implicit_dict, implicit_id_dict,
                implicit_weak_dict, implicit_set, implicit_vector,
                implicit_sized_vector, implicit_matrix
        end
    """)
    @test count(site -> site.category == "field_stored_any",
        prohibited.stored_any) == 6
    @test count(site -> site.category == "signature_erased_storage",
        prohibited.stored_any) == 1
    @test count(site -> site.category == "prepared_return_any",
        prohibited.stored_any) == 1
    @test count(site -> site.category == "any_array_literal",
        prohibited.stored_any) == 2
    @test count(site -> site.category == "erased_storage_constructor",
        prohibited.stored_any) == 3
    @test count(site -> site.category ==
            "implicit_any_storage_constructor",
        prohibited.stored_any) == 7
    @test count(site -> site.category == "implicit_any_array_literal",
        prohibited.stored_any) == 1
    @test count(site -> site.category == "return_erased_storage",
        prohibited.stored_any) == 1
    @test count(site -> site.category == "field_type",
        prohibited.memory) == 3
    @test count(site -> site.category == "constructor",
        prohibited.memory) == 2
end

@testset "PE-00 prepared execution contract" begin
    contract_path = joinpath(
        PE_CONTRACT_ROOT,
        "prepared_execution_ownership.toml",
    )
    @test isfile(contract_path)
    contract = TOML.parsefile(contract_path)
    @test contract["schema_version"] == 1
    @test contract["name"] == "prepared_execution_ownership"
    @test contract["status"] == "characterized_baseline"

    decisions = contract["interface_decisions"]
    @test decisions["work_issue"] == 229
    @test decisions["delivery_state"] == "implemented"
    decision_columns = String.(decisions["columns"])
    @test decision_columns == [
        "owner",
        "name",
        "decision",
        "contract_kind",
        "follow_on_gate",
        "rationale",
    ]
    allowed_decisions = Set(String.(decisions["allowed_decisions"]))
    allowed_contract_kinds =
        Set(String.(decisions["allowed_contract_kinds"]))
    decision_rows = Dict{
        Tuple{String,String},
        NamedTuple{
            (:decision, :contract_kind, :follow_on_gate, :rationale),
            Tuple{String,String,String,String},
        },
    }()
    for row in decisions["entries"]
        @test length(row) == length(decision_columns)
        values = Dict(decision_columns .=> row)
        key = (String(values["owner"]), String(values["name"]))
        @test !haskey(decision_rows, key)
        @test values["decision"] in allowed_decisions
        @test values["contract_kind"] in allowed_contract_kinds
        @test !isempty(values["follow_on_gate"])
        @test !isempty(values["rationale"])
        decision_rows[key] = (
            decision=String(values["decision"]),
            contract_kind=String(values["contract_kind"]),
            follow_on_gate=String(values["follow_on_gate"]),
            rationale=String(values["rationale"]),
        )
    end
    @test length(decision_rows) == 17

    retained_interfaces = Set((
        ("Backends", "_AbstractPreparedDeviceExecutionContext"),
        ("Calibration", "AbstractCalibrationCommandPlan"),
        ("Plant", "AbstractPreparedAcquisitionLifecycle"),
        ("Plant", "AbstractPreparedDetectorAcquisitionLifecycle"),
        ("Plant", "PreparedPupilOPDPublicationRoute"),
        ("Plant", "_PreparedEffectiveCommandPublicationRoute"),
        ("WavefrontSensors", "AbstractPreparedFourPupilLGS"),
    ))
    removed_interfaces = Set((
        ("Optics", "AbstractDirectImagingInputPlan"),
        ("Plant", "_PreparedReducedOrderEventResponse"),
        ("WavefrontSensors",
            "AbstractPreparedCurvatureObservationMapping"),
    ))
    strategy_interfaces = Set((
        ("Atmospheres", "AbstractAtmosphericFieldExecutionStrategy"),
        ("Backends", "AbstractReductionExecutionStrategy"),
        ("Detectors", "AbstractDetectorExecutionStrategy"),
        ("WavefrontSensors", "AbstractGroupedAccumulationStrategy"),
        ("WavefrontSensors", "AbstractShackHartmannWFSSensingStrategy"),
    ))
    deferred_interfaces = Set((
        ("Plant", "_AbstractPreparedMixedSerialEventAcquisition"),
        ("Plant", "_AbstractPreparedMixedSerialPath"),
    ))
    decision_keys(value) = Set(key for (key, row) in decision_rows
        if row.decision == value)
    @test decision_keys("retain") == retained_interfaces
    @test decision_keys("remove_now") == removed_interfaces
    @test decision_keys("reclassify_strategy") == strategy_interfaces
    @test decision_keys("defer_removal") == deferred_interfaces
    @test all(decision_rows[key].follow_on_gate == "PE-01"
        for key in union(retained_interfaces, removed_interfaces))
    @test all(decision_rows[key].follow_on_gate == "PE-02"
        for key in strategy_interfaces)
    @test all(decision_rows[key].follow_on_gate == "PE-06"
        for key in deferred_interfaces)

    suite_names = Set(spec.name for spec in TEST_SUITE_SPECS)
    interface_contracts = Dict{
        Tuple{String,String},
        NamedTuple{
            (:implementations, :required_methods),
            Tuple{Vector{String},Vector{String}},
        },
    }()
    conformance_source_parts = String[]
    for (root, _, files) in walkdir(joinpath(PE_REPOSITORY_ROOT, "test"))
        for file in files
            endswith(file, ".jl") || continue
            push!(conformance_source_parts,
                read(joinpath(root, file), String))
        end
    end
    conformance_source = join(conformance_source_parts, '\n')
    for entry in contract["interface_contracts"]
        key = (String(entry["owner"]), String(entry["name"]))
        @test !haskey(interface_contracts, key)
        @test entry["visibility"] in
            ("internal", "internal_extension", "qualified_public")
        @test !isempty(entry["implementations"])
        @test all(!isempty, entry["implementations"])
        @test length(unique(entry["implementations"])) ==
            length(entry["implementations"])
        @test !isempty(entry["required_methods"])
        @test all(!isempty, entry["required_methods"])
        @test length(unique(entry["required_methods"])) ==
            length(entry["required_methods"])
        @test entry["valid_defaults"] isa Vector
        @test !isempty(entry["ownership"])
        @test !isempty(entry["aliasing"])
        @test !isempty(entry["failure_behavior"])
        @test !isempty(entry["reentrancy"])
        @test !isempty(entry["backend_traits"])
        @test !isempty(entry["conformance_helpers"])
        @test all(helper -> occursin(helper, conformance_source),
            entry["conformance_helpers"])
        @test !isempty(entry["conformance_suites"])
        @test all(in(suite_names), entry["conformance_suites"])
        interface_contracts[key] = (
            implementations=String.(entry["implementations"]),
            required_methods=String.(entry["required_methods"]),
        )
    end
    @test Set(keys(interface_contracts)) == retained_interfaces

    inventory = contract["type_inventory"]
    @test Set(String.(inventory["vocabulary_tokens"])) ==
        PE_VOCABULARY
    columns = String.(inventory["columns"])
    @test columns == [
        "name",
        "declaration",
        "field_count",
        "current_role",
        "target_roles",
        "mutable",
        "persistent_state",
        "product_visibility",
        "exact_binding",
        "backend_relevance",
        "cpu_hot_path",
        "migration_gate",
    ]
    @test !isempty(inventory["target_roles_semantics"])

    expected = Dict{Tuple{String,String},Tuple{String,Bool,Int}}()
    mixed_target_count = 0
    allowed_current = Set(String.(inventory["allowed_current_roles"]))
    allowed_target = Set(String.(inventory["allowed_target_roles"]))
    allowed_visibility = Set(String.(inventory["allowed_product_visibility"]))
    allowed_backend = Set(String.(inventory["allowed_backend_relevance"]))
    allowed_gates = Set(String.(inventory["allowed_migration_gates"]))
    @test issorted(String(file["path"]) for file in inventory["files"])
    for file in inventory["files"]
        path = String(file["path"])
        @test isfile(joinpath(PE_REPOSITORY_ROOT, path))
        @test !isempty(file["owner"])
        @test issorted(String(row[1]) for row in file["entries"])
        for row in file["entries"]
            @test length(row) == length(columns)
            values = Dict(columns .=> row)
            key = (path, String(values["name"]))
            @test !haskey(expected, key)
            @test values["declaration"] in
                ("struct", "mutable_struct", "abstract_type",
                    "primitive_type", "type_alias", "enum")
            @test values["field_count"] isa Integer
            @test values["current_role"] in allowed_current
            @test values["target_roles"] isa Vector
            @test !isempty(values["target_roles"])
            @test all(in(allowed_target), values["target_roles"])
            @test length(unique(values["target_roles"])) ==
                length(values["target_roles"])
            mixed_target_count += length(values["target_roles"]) > 1
            @test values["mutable"] isa Bool
            @test values["persistent_state"] isa Bool
            @test values["product_visibility"] in allowed_visibility
            @test values["exact_binding"] isa Bool
            @test values["backend_relevance"] in allowed_backend
            @test values["cpu_hot_path"] isa Bool
            @test values["migration_gate"] in allowed_gates
            values["persistent_state"] &&
                @test "state" in values["target_roles"]
            values["product_visibility"] == "caller_visible" &&
                @test "product" in values["target_roles"]
            values["target_roles"] == ["binding_token"] && begin
                @test values["field_count"] == 0
                @test values["exact_binding"]
            end
            if "strategy" in values["target_roles"]
                if values["declaration"] == "abstract_type"
                    @test values["field_count"] == -1
                else
                    @test values["field_count"] == 0
                end
            end
            values["current_role"] == "plan" &&
                    values["field_count"] == 0 &&
                @test !isdisjoint(values["target_roles"],
                    ("strategy", "binding_token"))
            expected[key] = (
                String(values["declaration"]),
                values["mutable"],
                values["field_count"],
            )
        end
    end
    @test mixed_target_count == inventory["baseline_mixed_type_count"]
    @test length(expected) == inventory["baseline_type_count"]

    inventory_by_name = Dict{
        String,
        NamedTuple{
            (:path, :declaration, :current_role, :target_roles,
                :migration_gate),
            Tuple{String,String,String,Vector{String},String},
        },
    }()
    for file in inventory["files"]
        path = String(file["path"])
        for row in file["entries"]
            values = Dict(columns .=> row)
            name = String(values["name"])
            @test !haskey(inventory_by_name, name)
            inventory_by_name[name] = (
                path=path,
                declaration=String(values["declaration"]),
                current_role=String(values["current_role"]),
                target_roles=String.(values["target_roles"]),
                migration_gate=String(values["migration_gate"]),
            )
        end
    end

    owner_modules = Dict(
        "Atmospheres" => Atmospheres,
        "Backends" => Backends,
        "Calibration" => Calibration,
        "Detectors" => Detectors,
        "Optics" => Optics,
        "Plant" => Plant,
        "WavefrontSensors" => WavefrontSensors,
    )
    for key in retained_interfaces
        owner, name = key
        @test haskey(inventory_by_name, name)
        @test inventory_by_name[name].declaration == "abstract_type"
        root = getfield(owner_modules[owner], Symbol(name))
        @test isabstracttype(root)
        contract_entry = interface_contracts[key]
        @test all(implementation ->
                haskey(inventory_by_name, implementation) ||
                isdefined(owner_modules[owner], Symbol(implementation)),
            contract_entry.implementations)
        @test all(method_name -> isdefined(owner_modules[owner],
                Symbol(method_name)), contract_entry.required_methods)
    end
    for (owner, name) in removed_interfaces
        @test !haskey(inventory_by_name, name)
        @test !isdefined(owner_modules[owner], Symbol(name))
    end
    for (_, name) in strategy_interfaces
        @test inventory_by_name[name].current_role == "strategy"
        @test inventory_by_name[name].target_roles == ["strategy"]
        @test inventory_by_name[name].migration_gate == "PE-02"
        @test endswith(name, "Strategy")
    end
    for (_, name) in deferred_interfaces
        @test inventory_by_name[name].declaration == "abstract_type"
        @test inventory_by_name[name].migration_gate == "PE-06"
    end

    observed = Dict(
        (declaration.path, declaration.name) =>
            (declaration.declaration, declaration.mutable,
                declaration.field_count)
        for declaration in pe_type_declarations(PE_REPOSITORY_ROOT)
    )
    @test expected == observed

    debt = pe_debt_sites(PE_REPOSITORY_ROOT)
    for (section_name, sites) in (
        "direct_memory" => debt.memory,
        "stored_any" => debt.stored_any,
    )
        section = contract[section_name]
        @test section["policy"] == "ratchet"
        site_columns = String.(section["columns"])
        baseline = Dict{Tuple{String,String,String,String},Int}()
        for file in section["files"]
            path = String(file["path"])
            @test isfile(joinpath(PE_REPOSITORY_ROOT, path))
            @test !isempty(file["owner"])
            for row in file["entries"]
                @test length(row) == length(site_columns)
                values = Dict(site_columns .=> row)
                key = (
                    path,
                    String(values["category"]),
                    String(values["scope"]),
                    String(values["binding"]),
                )
                @test !haskey(baseline, key)
                @test values["baseline_line"] > 0
                @test values["count"] > 0
                @test values["hot_path"] isa Bool
                @test !isempty(values["replacement"])
                section_name == "direct_memory" &&
                    @test values["backend_relevance"] ==
                        "host_ownership_storage"
                baseline[key] = values["count"]
            end
        end
        observed_counts = pe_debt_counts(sites)
        observed_keys = Set(keys(observed_counts))
        baseline_keys = Set(keys(baseline))
        @test isempty(setdiff(observed_keys, baseline_keys))
        @test all(observed_counts[key] <= baseline[key]
            for key in intersect(observed_keys, baseline_keys))
        @test sum(values(observed_counts)) <= section["baseline_total"]
        @test section["baseline_total"] == sum(values(baseline))
    end

    evidence = contract["evidence_index"]
    @test Set(String.(evidence["classes"])) == Set((
        "numerical",
        "inference",
        "allocation",
        "storage",
        "package_load",
        "cpu",
        "cuda",
        "amdgpu",
    ))
    for entry in evidence["entries"]
        @test entry["class"] in evidence["classes"]
        @test isfile(joinpath(PE_REPOSITORY_ROOT, entry["path"]))
        @test !isempty(entry["claim_limit"])
    end

    hil = contract["adaptive_optics_hil"]
    @test length(hil["source_revision"]) == 40
    hil_inventory_path = joinpath(
        PE_REPOSITORY_ROOT, hil["import_inventory"])
    @test isfile(hil_inventory_path)
    hil_inventory = TOML.parsefile(hil_inventory_path)
    @test hil_inventory["schema_version"] == 1
    @test hil_inventory["status"] == "frozen_source_snapshot"
    @test hil_inventory["repository"] == hil["repository"]
    @test hil_inventory["source_revision"] == hil["source_revision"]
    @test hil_inventory["adaptive_optics_sim_pin"] ==
        hil["adaptive_optics_sim_pin"]
    @test all(hash -> length(hash) == 64, (
        hil_inventory["project_files"]["package_sha256"],
        hil_inventory["project_files"]["benchmarks_sha256"],
    ))
    @test Set((
        "PreparedPlantEventLoop",
        "PlantEventLoopState",
        "PlantEventLoopWorkspace",
    )) ⊆ Set(String.(hil["affected_owner_types"]))
    affected = Set(String.(hil["affected_owner_types"]))
    @test affected == Set(String.(
        hil_inventory["ownership"]["affected_symbols"]))
    @test !isempty(hil["construction_assumptions"])
    @test hil["construction_assumptions"] ==
        hil_inventory["ownership"]["construction_assumptions"]
    referenced = Set{String}()
    for entry in hil_inventory["reference_files"]
        @test entry["scope"] in hil_inventory["scope"]
        @test startswith(entry["path"], entry["scope"] * "/")
        @test length(entry["sha256"]) == 64
        @test !isempty(entry["symbols"])
        union!(referenced, String.(entry["symbols"]))
    end
    @test referenced == affected
    hil_memory = hil_inventory["direct_memory"]
    @test hil_memory["migration_gate"] == "PE-08"
    memory_files = Dict(String(row[1]) => Int(row[2])
        for row in hil_memory["files"])
    @test length(memory_files) == length(hil_memory["files"])
    memory_count(prefix) = sum((count for (path, count) in memory_files
        if startswith(path, prefix)); init=0)
    @test hil_memory["source_total"] == memory_count("src/")
    @test hil_memory["test_total"] == memory_count("test/")
    @test hil_memory["benchmarks_total"] == memory_count("benchmarks/")

    fixed = contract["fixed_size_arrays"]
    @test fixed["core_dependency"] === false
    @test fixed["characterized_version"] == "1.3.0"
    @test isfile(joinpath(PE_REPOSITORY_ROOT, fixed["project"]))
    @test isfile(joinpath(PE_REPOSITORY_ROOT, fixed["contract"]))
    @test isfile(joinpath(PE_REPOSITORY_ROOT, fixed["artifact"]))
    fixed_contract = TOML.parsefile(joinpath(
        PE_REPOSITORY_ROOT, fixed["contract"]))
    fixed_artifact = TOML.parsefile(joinpath(
        PE_REPOSITORY_ROOT, fixed["artifact"]))
    @test fixed_artifact["correctness_passed"]
    @test !fixed_artifact["repository_dirty"]
    @test length(fixed_artifact["repository_head"]) == 40
    @test fixed_artifact["julia_version"] ==
        fixed_contract["julia_version"]
    @test fixed_artifact["fresh_process_runs"] ==
        fixed_contract["fresh_process_runs"]
    @test fixed_artifact["fixed_size_arrays_version"] ==
        fixed["characterized_version"]
    @test fixed_artifact["julia_threads"] == 1
    @test fixed_artifact["logical_cpu_threads"] > 0
    for key in (
        "kernel",
        "architecture",
        "cpu_target",
        "cpu_model",
        "julia_cpu_target_env",
    )
        @test !isempty(fixed_artifact[key])
    end
    @test all(iszero, values(
        fixed_artifact["warmed_allocated_bytes"]))
    @test all(identity, values(fixed_artifact["inference"]))
    @test fixed_artifact["shape_contract"][
        "runtime_vector_lengths_share_type"]
    @test fixed_artifact["shape_contract"][
        "runtime_matrix_shapes_share_type"]
    @test fixed_artifact["shape_contract"]["resize_rejected"]
    @test fixed_artifact["shape_contract"]["push_rejected"]
    @test fixed_artifact["element_type_contract"][
        "any_is_constructible"]
    @test fixed_artifact["element_type_contract"][
        "concrete_is_accepted"]
    @test fixed_artifact["element_type_contract"]["any_is_rejected"]
    @test !fixed_artifact["storage_comparison"][
        "dependency_backing_inspected"]
end
