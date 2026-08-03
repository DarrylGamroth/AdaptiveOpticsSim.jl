const PE_CONTRACT_ROOT = joinpath(dirname(@__DIR__), "contracts")
const PE_REPOSITORY_ROOT = dirname(dirname(@__DIR__))

@testset "PE-00 source scanner semantics" begin
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
            "strategy" in values["target_roles"] &&
                @test values["field_count"] == 0
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
