using SHA

const NAMESPACE_CONTRACT_ROOT =
    joinpath(dirname(@__DIR__), "contracts")
const REPOSITORY_ROOT = dirname(dirname(@__DIR__))

function declared_surface(path::AbstractString)
    exported = String[]
    qualified_public = String[]
    for line in eachline(path)
        stripped = strip(line)
        startswith(stripped, "export ") ||
            startswith(stripped, "public ") || continue
        kind, bindings = split(stripped, ' '; limit=2)
        destination = kind == "export" ? exported : qualified_public
        for binding in strip.(split(bindings, ','))
            binding in destination || push!(destination, binding)
        end
    end
    return (exports=exported, public=qualified_public)
end

file_sha256(path::AbstractString) = bytes2hex(SHA.sha256(codeunits(
    replace(read(path, String), "\r\n" => "\n"),
)))

function runtime_surface(mod::Module)
    bindings = names(mod; all=true, imported=false)
    exported = sort!(String[
        String(binding) for binding in bindings
        if Base.isexported(mod, binding)
    ])
    qualified_public = sort!(String[
        String(binding) for binding in bindings
        if Base.ispublic(mod, binding) && !Base.isexported(mod, binding)
    ])
    return (exports=exported, public=qualified_public)
end

binding_parentmodule(value::Enum) = parentmodule(typeof(value))
binding_parentmodule(value) = parentmodule(value)

function assert_acyclic_owner_graph(graph)
    dependencies = Dict(
        String(node["name"]) => String.(node["dependencies"])
        for node in graph
    )
    temporary = Set{String}()
    permanent = Set{String}()
    function visit(owner)
        owner in permanent && return
        @test owner ∉ temporary
        push!(temporary, owner)
        for dependency in dependencies[owner]
            @test haskey(dependencies, dependency)
            visit(dependency)
        end
        delete!(temporary, owner)
        push!(permanent, owner)
        return
    end
    foreach(visit, keys(dependencies))
    @test length(permanent) == length(dependencies)
end

function relative_import_surfaces(path::AbstractString)
    imports = Dict{String,Vector{String}}()
    provider = nothing
    for line in eachline(path)
        stripped = strip(line)
        if startswith(stripped, "import ..")
            declaration = split(stripped, ':'; limit=2)
            length(declaration) == 2 || continue
            provider = String(chop(
                declaration[1];
                head=length("import .."),
                tail=0,
            ))
            imports[provider] = String[]
            append!(imports[provider], strip.(filter(
                !isempty,
                split(declaration[2], ','),
            )))
        elseif !isnothing(provider) && isempty(stripped)
            provider = nothing
        elseif !isnothing(provider)
            append!(imports[provider], strip.(filter(
                !isempty,
                split(stripped, ','),
            )))
        end
    end
    foreach(unique!, values(imports))
    return imports
end

function binding_occurs(text::AbstractString, binding::AbstractString)
    escaped = replace(binding, r"([\\.^$|?*+(){}\\[\\]])" => s"\\\1")
    return occursin(
        Regex("(?<![A-Za-z0-9_!])" * escaped * "(?![A-Za-z0-9_!])"),
        text,
    )
end

function adaptive_optics_sim_extension_hooks(path::AbstractString)
    text = read(path, String)
    patterns = (
        "Root" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?AdaptiveOpticsSim\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Backends" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Backends\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Detectors" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Detectors\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Atmospheres" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Atmospheres\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "WavefrontSensors" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?WavefrontSensors\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Calibration" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Calibration\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Control" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Control\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Tomography" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Tomography\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
        "Ensembles" =>
            r"(?m)^(?:@inline\s+)?(?:function\s+)?(?:AdaptiveOpticsSim\.)?Ensembles\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()",
    )
    owners = Dict{String,String}()
    for (owner, pattern) in patterns
        for matched in eachmatch(pattern, text)
            hook = String(matched.captures[1])
            @test get(owners, hook, owner) == owner
            owners[hook] = owner
        end
    end
    return owners
end

function adaptive_optics_sim_root_references(path::AbstractString)
    text = read(path, String)
    aliases = ["AdaptiveOpticsSim"]
    append!(aliases, String(matched.captures[1]) for matched in eachmatch(
        r"(?m)^\s*const\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*AdaptiveOpticsSim\s*(?:#.*)?$",
        text,
    ))
    references = Pair{String,String}[]
    for alias in unique(aliases)
        pattern = Regex(
            "(?<![A-Za-z0-9_])" * alias *
            raw"\.([A-Za-z_][A-Za-z0-9_!]*)",
        )
        append!(references,
            alias => String(matched.captures[1])
            for matched in eachmatch(pattern, text))
    end
    return unique(references)
end

function recursive_key_values(value, sought::Set{String}, found=Pair{String,Any}[])
    if value isa AbstractDict
        for (key, child) in value
            String(key) in sought && push!(found, String(key) => child)
            recursive_key_values(child, sought, found)
        end
    elseif value isa AbstractVector
        foreach(child -> recursive_key_values(child, sought, found), value)
    end
    return found
end

function seed_assignments(value, path::AbstractString="")
    assignments = String[]
    if value isa AbstractDict
        for (key, child) in value
            child_path = isempty(path) ? String(key) : "$path.$key"
            if occursin("seed", lowercase(String(key))) &&
                    !(child isa AbstractDict || child isa AbstractVector)
                push!(assignments, "$child_path=$child")
            end
            append!(assignments, seed_assignments(child, child_path))
        end
    elseif value isa AbstractVector
        for (index, child) in enumerate(value)
            append!(assignments,
                seed_assignments(child, "$path[$(index - 1)]"))
        end
    end
    return assignments
end

@testset "NS-00B namespace authority" begin
    authority_path =
        joinpath(NAMESPACE_CONTRACT_ROOT, "namespace_authority.toml")
    authority = TOML.parsefile(authority_path)
    migration_state = TOML.parsefile(joinpath(
        NAMESPACE_CONTRACT_ROOT,
        "namespace_migration_state.toml",
    ))
    @test authority["schema_version"] == 1
    @test authority["name"] == "namespace_authority"
    @test authority["status"] == "authoritative"
    @test authority["implementation_phase"] == "pre_namespace_refactor"
    @test authority["semantics"]["compatibility"] ==
        "no root forwarding aliases, property forwarding, state views, or compatibility adapters"
    @test migration_state["schema_version"] == 1
    @test migration_state["name"] == "namespace_migration_state"
    @test migration_state["status"] == "complete"
    @test migration_state["authority"] == basename(authority_path)

    current = authority["current_root"]
    current_path = joinpath(REPOSITORY_ROOT, current["source"])
    @test isfile(current_path)
    @test length(current["sha256"]) == 64
    @test length(current["exports"]) == length(unique(current["exports"]))
    @test length(current["public"]) == length(unique(current["public"]))
    baseline_divergence = migration_state["baseline_divergence"]
    @test baseline_divergence["root_surface"] == current["source"]
    @test file_sha256(current_path) != current["sha256"]

    canonical_owners = Set((
        "Root",
        "Backends",
        "Optics",
        "Atmospheres",
        "Detectors",
        "WavefrontSensors",
        "Plant",
        "Calibration",
        "Control",
        "Tomography",
        "Ensembles",
    ))
    allowlists = authority["api_allowlists"]
    @test Set(String(entry["owner"]) for entry in allowlists) ==
        canonical_owners
    implemented_owners =
        Set(String.(migration_state["implemented_owners"]))
    partial_owners =
        Set(String.(get(migration_state, "partial_owners", String[])))
    @test implemented_owners == setdiff(canonical_owners, Set(("Root",)))
    @test isempty(partial_owners)
    @test isdisjoint(implemented_owners, partial_owners)
    @test "Root" ∉ implemented_owners
    @test "Root" ∉ partial_owners
    @test Set(keys(migration_state["owner_sources"])) ==
        implemented_owners
    @test Set(keys(get(
        migration_state, "partial_owner_sources", Dict{String,Any}()))) ==
        partial_owners
    @test Set(keys(get(
        migration_state, "partial_owner_bindings", Dict{String,Any}()))) ==
        partial_owners

    partial_bindings = Dict{String,Dict{String,Set{String}}}()
    for owner in partial_owners
        stage = migration_state["partial_owner_bindings"][owner]
        @test Set(keys(stage)) == Set(("exports", "public", "internal"))
        bindings = Dict(
            visibility => Set(String.(stage[visibility]))
            for visibility in ("exports", "public", "internal")
        )
        @test all(visibility -> length(bindings[visibility]) ==
            length(stage[visibility]), keys(bindings))
        @test isdisjoint(bindings["exports"], bindings["public"])
        @test isdisjoint(bindings["exports"], bindings["internal"])
        @test isdisjoint(bindings["public"], bindings["internal"])
        partial_bindings[owner] = bindings
    end

    function binding_is_migrated(owner::AbstractString,
        binding::AbstractString)
        owner in implemented_owners && return true
        owner in partial_owners || return false
        return any(
            binding in partial_bindings[owner][visibility]
            for visibility in ("exports", "public", "internal")
        )
    end

    owner_by_binding = Dict{String,Tuple{String,String}}()
    for entry in allowlists
        owner = String(entry["owner"])
        for visibility in ("exports", "public")
            names = String.(entry[visibility])
            @test length(names) == length(unique(names))
            for binding in names
                @test !haskey(owner_by_binding, binding)
                owner_by_binding[binding] = (owner, visibility)
            end
        end
    end

    migration = authority["migration"]
    introduced_root = Set(String.(migration["new_root_exports"]))
    introduced_public = Set(String.(migration["new_domain_public"]))
    export_to_public = Set(String.(migration["export_to_public"]))
    removed_exports = Set(String.(migration["removed_exports"]))
    removed_public = Set(String.(migration["removed_public"]))
    @test removed_exports ⊆ Set(String.(current["exports"]))
    @test removed_public ⊆ Set(String.(current["public"]))
    @test isdisjoint(removed_exports, introduced_root)
    @test isdisjoint(removed_public, introduced_public)
    @test isdisjoint(introduced_public, export_to_public)
    @test all(binding -> owner_by_binding[binding] == ("Root", "exports"),
        introduced_root)
    @test all(introduced_public) do qualified
        owner, binding = split(qualified, '.'; limit=2)
        owner_by_binding[binding] == (owner, "public")
    end
    @test all(export_to_public) do qualified
        owner, binding = split(qualified, '.'; limit=2)
        binding in current["exports"] &&
            owner_by_binding[binding] == (owner, "public")
    end
    export_to_public_bindings = Set(
        last(split(qualified, '.'; limit=2))
        for qualified in export_to_public
    )

    existing_target = Set{String}()
    for (binding, (owner, visibility)) in owner_by_binding
        owner == "Plant" && continue
        binding in introduced_root && continue
        "$owner.$binding" in union(introduced_public,
            export_to_public) && continue
        push!(existing_target, binding)
    end
    @test existing_target == Set([
        (binding for binding in current["exports"]
            if binding ∉ export_to_public_bindings &&
                binding ∉ removed_exports)...;
        (binding for binding in current["public"]
            if binding ∉ removed_public)...;
    ])
    @test all(binding -> binding in removed_exports ||
            binding in export_to_public_bindings ||
            owner_by_binding[binding][2] == "exports",
        current["exports"])
    @test all(binding -> binding in removed_public ||
            owner_by_binding[binding][2] == "public",
        current["public"])

    implemented_exports = Set{String}()
    implemented_public = Set{String}()
    for owner in implemented_owners
        allowlist = only(entry for entry in allowlists
            if entry["owner"] == owner)
        union!(implemented_exports, String.(allowlist["exports"]))
        union!(implemented_public, String.(allowlist["public"]))
    end
    for owner in partial_owners
        allowlist = only(entry for entry in allowlists
            if entry["owner"] == owner)
        @test partial_bindings[owner]["exports"] ⊆
            Set(String.(allowlist["exports"]))
        @test partial_bindings[owner]["public"] ⊆
            Set(String.(allowlist["public"]))
        union!(implemented_exports, partial_bindings[owner]["exports"])
        union!(implemented_public, partial_bindings[owner]["public"])
    end
    expected_root_exports = Set(
        binding for binding in current["exports"]
        if binding ∉ implemented_exports &&
            binding ∉ export_to_public_bindings &&
            binding ∉ removed_exports
    )
    union!(expected_root_exports,
        intersect(introduced_root,
            union(implemented_owners, partial_owners)))
    expected_root_public = Set(
        binding for binding in current["public"]
        if binding ∉ implemented_public && binding ∉ removed_public
    )

    declared_root = declared_surface(current_path)
    @test Set(declared_root.exports) == expected_root_exports
    @test Set(declared_root.public) == expected_root_public
    root_runtime = runtime_surface(AdaptiveOpticsSim)
    @test Set(root_runtime.exports) ==
        union(expected_root_exports, Set(("AdaptiveOpticsSim",)))
    @test Set(root_runtime.public) == expected_root_public
    domain_modules = Set(
        getfield(AdaptiveOpticsSim, Symbol(owner))
        for owner in setdiff(canonical_owners, Set(("Root",)))
    )
    root_domain_aliases = Symbol[]
    for binding in names(AdaptiveOpticsSim; all=true, imported=true)
        isdefined(AdaptiveOpticsSim, binding) || continue
        value = getfield(AdaptiveOpticsSim, binding)
        value_owner = try
            binding_parentmodule(value)
        catch
            AdaptiveOpticsSim
        end
        value_owner in domain_modules &&
            push!(root_domain_aliases, binding)
    end
    @test isempty(root_domain_aliases)

    benchmark_paths = String[]
    for (directory, _, files) in walkdir(
        joinpath(REPOSITORY_ROOT, "benchmarks"))
        append!(benchmark_paths,
            joinpath(directory, file)
            for file in files if endswith(file, ".jl"))
    end
    for path in benchmark_paths
        for (alias, binding) in adaptive_optics_sim_root_references(path)
            @test isdefined(AdaptiveOpticsSim, Symbol(binding))
        end
    end

    for owner in implemented_owners
        allowlist = only(entry for entry in allowlists
            if entry["owner"] == owner)
        source = joinpath(
            REPOSITORY_ROOT,
            migration_state["owner_sources"][owner],
        )
        declared_owner = declared_surface(source)
        @test declared_owner.exports == allowlist["exports"]
        @test declared_owner.public == allowlist["public"]

        owner_module = getfield(AdaptiveOpticsSim, Symbol(owner))
        owner_runtime = runtime_surface(owner_module)
        @test Set(owner_runtime.exports) ==
            Set([allowlist["exports"]; owner])
        @test Set(owner_runtime.public) == Set(allowlist["public"])
        for visibility in ("exports", "public")
            for binding in allowlist[visibility]
                value = getfield(owner_module, Symbol(binding))
                @test binding_parentmodule(value) === owner_module
                @test binding ∉ root_runtime.exports
                @test binding ∉ root_runtime.public
            end
        end
    end
    for owner in partial_owners
        source = joinpath(
            REPOSITORY_ROOT,
            migration_state["partial_owner_sources"][owner],
        )
        declared_owner = declared_surface(source)
        @test Set(declared_owner.exports) ==
            partial_bindings[owner]["exports"]
        @test Set(declared_owner.public) ==
            partial_bindings[owner]["public"]

        owner_module = getfield(AdaptiveOpticsSim, Symbol(owner))
        owner_runtime = runtime_surface(owner_module)
        @test Set(owner_runtime.exports) ==
            union(partial_bindings[owner]["exports"], Set((owner,)))
        @test Set(owner_runtime.public) ==
            partial_bindings[owner]["public"]
        for visibility in ("exports", "public", "internal")
            for binding in partial_bindings[owner][visibility]
                value = getfield(owner_module, Symbol(binding))
                @test binding_parentmodule(value) === owner_module
                @test binding ∉ root_runtime.exports
                @test binding ∉ root_runtime.public
            end
        end
    end

    graph = authority["owner_graph"]
    @test Set(String(node["name"]) for node in graph) ==
        setdiff(canonical_owners, Set(("Root",)))
    assert_acyclic_owner_graph(graph)

    physical = authority["physical_ownership"]
    @test physical["ncpa_model_owner"] == "Optics"
    @test physical["ncpa_visibility_owner"] == "Plant"
    @test physical["kl_zernike_m2c_synthesis_owner"] == "Calibration"
    @test owner_by_binding["NCPA"] == ("Optics", "exports")
    @test owner_by_binding["KLBasis"] == ("Calibration", "exports")
    @test owner_by_binding["ZernikeModalBasis"] ==
        ("Calibration", "exports")

    terminology = authority["terminology"]
    glossary = read(joinpath(REPOSITORY_ROOT, terminology["glossary"]),
        String)
    @test all(name -> haskey(owner_by_binding, name),
        terminology["canonical_api_names"])
    @test all(term -> occursin(term, glossary),
        terminology["canonical_glossary_terms"])

    plant_contract = authority["plant_surface"]
    plant_path = joinpath(REPOSITORY_ROOT, plant_contract["source"])
    @test file_sha256(plant_path) == plant_contract["sha256"]
    plant_declared = declared_surface(plant_path)
    plant_allowlist = only(entry for entry in allowlists
        if entry["owner"] == "Plant")
    @test plant_declared.exports == plant_allowlist["exports"]
    @test plant_declared.public == plant_allowlist["public"]
    plant_runtime = runtime_surface(Plant)
    @test Set(plant_runtime.exports) ==
        Set([plant_declared.exports; "Plant"])
    @test Set(plant_runtime.public) == Set(plant_declared.public)

    import_contracts = authority["cross_module_imports"]
    imported_bindings = String[]
    expected_imports = Dict{String,Set{String}}()
    plant_body = join(
        read(path, String) for path in readdir(
            joinpath(REPOSITORY_ROOT, "src", "plant");
            join=true,
        ) if endswith(path, ".jl") && basename(path) != "plant.jl"
    )
    for entry in import_contracts
        @test entry["consumer"] == "Plant"
        @test entry["provider"] in canonical_owners
        @test entry["action"] in ("retain", "remove_unused_import")
        append!(imported_bindings, String.(entry["bindings"]))
        if entry["action"] == "retain"
            @test all(binding -> binding_occurs(plant_body, binding),
                entry["bindings"])
            for binding in String.(entry["bindings"])
                provider = binding_is_migrated(
                    String(entry["provider"]), binding) ?
                    String(entry["provider"]) : "AdaptiveOpticsSim"
                push!(get!(expected_imports, provider, Set{String}()),
                    binding)
            end
        else
            @test all(binding -> !binding_occurs(plant_body, binding),
                entry["bindings"])
        end
    end
    @test length(imported_bindings) == length(unique(imported_bindings))
    actual_imports = relative_import_surfaces(joinpath(
        REPOSITORY_ROOT,
        "src",
        "plant",
        "plant.jl",
    ))
    @test Set(keys(actual_imports)) == Set(keys(expected_imports))
    for (provider, bindings) in expected_imports
        @test Set(actual_imports[provider]) == bindings
    end

    extension_files = authority["extension_files"]
    observed_hooks = Set{String}()
    for entry in extension_files
        path = joinpath(REPOSITORY_ROOT, entry["path"])
        hooks = adaptive_optics_sim_extension_hooks(path)
        @test Set(keys(hooks)) ==
            Set(String.(entry["adaptive_optics_sim_hooks"]))
        union!(observed_hooks, keys(hooks))
    end
    hook_owners = Dict{String,String}()
    for entry in authority["extension_hook_owners"]
        for hook in entry["hooks"]
            @test !haskey(hook_owners, hook)
            hook_owners[hook] = entry["owner"]
        end
    end
    @test Set(keys(hook_owners)) == observed_hooks
    @test all(owner -> owner in canonical_owners, values(hook_owners))
    migrated_extension_sources = Set{String}()
    for entry in extension_files
        path = joinpath(REPOSITORY_ROOT, entry["path"])
        hooks = adaptive_optics_sim_extension_hooks(path)
        for (hook, actual_owner) in hooks
            canonical_owner = hook_owners[hook]
            expected_owner = binding_is_migrated(canonical_owner, hook) ?
                canonical_owner : "Root"
            @test actual_owner == expected_owner
            actual_owner == "Root" ||
                push!(migrated_extension_sources, entry["path"])
        end
    end
    @test migrated_extension_sources ==
        Set(String.(baseline_divergence["extension_sources"]))
end

@testset "NS-00B AdaptiveOpticsHIL import inventory" begin
    authority = TOML.parsefile(joinpath(
        NAMESPACE_CONTRACT_ROOT,
        "namespace_authority.toml",
    ))
    inventory = TOML.parsefile(joinpath(
        NAMESPACE_CONTRACT_ROOT,
        "adaptiveopticshil_imports.toml",
    ))
    @test inventory["schema_version"] == 1
    @test inventory["status"] == "frozen_source_snapshot"
    @test length(inventory["source_revision"]) == 40
    @test inventory["scope"] == ["src", "test", "benchmarks"]

    owner_by_binding = Dict(
        String(binding) => String(entry["owner"])
        for entry in authority["api_allowlists"]
        for visibility in ("exports", "public")
        for binding in entry[visibility]
    )
    plant_allowlist = only(entry for entry in authority["api_allowlists"]
        if entry["owner"] == "Plant")
    plant_names =
        Set([plant_allowlist["exports"]; plant_allowlist["public"]])
    for entry in inventory["import_files"]
        @test entry["scope"] in inventory["scope"]
        @test startswith(entry["path"], entry["scope"] * "/")
        for access in entry["access"]
            @test access in (
                "root_using_wildcard",
                "root_import_module",
                "root_using_explicit",
                "root_import_explicit",
                "plant_using_wildcard",
                "plant_import_module",
                "plant_using_explicit",
                "plant_import_explicit",
            )
        end
        root_names = [
            entry["root_explicit"];
            entry["root_qualified"];
        ]
        @test all(name -> haskey(owner_by_binding, name), root_names)
        plant_references = [
            entry["plant_explicit"];
            entry["plant_qualified"];
        ]
        @test all(in(plant_names), plant_references)
    end
    @test any(
        entry -> "root_using_wildcard" in entry["access"],
        inventory["import_files"],
    )
    @test any(
        entry -> "plant_using_wildcard" in entry["access"],
        inventory["import_files"],
    )
end

@testset "NS-00B correctness and performance baselines" begin
    characterization = TOML.parsefile(joinpath(
        NAMESPACE_CONTRACT_ROOT,
        "namespace_characterization.toml",
    ))
    migration_state = TOML.parsefile(joinpath(
        NAMESPACE_CONTRACT_ROOT,
        "namespace_migration_state.toml",
    ))
    @test characterization["schema_version"] == 1
    @test characterization["status"] == "frozen_baseline_index"
    migrated_extension_tests = Set(String.(
        migration_state["baseline_divergence"]["extension_tests"],
    ))

    for fixture in characterization["numerical_fixtures"]
        path = joinpath(REPOSITORY_ROOT, fixture["manifest"])
        @test file_sha256(path) == fixture["sha256"]
        manifest = TOML.parsefile(path)
        @test length(manifest["cases"]) == fixture["case_count"]
        @test all(values(manifest["cases"])) do case
            all(field -> haskey(case, field),
                fixture["tolerance_fields"])
        end
        @test Set(seed_assignments(manifest)) ==
            Set(String.(fixture["seed_assignments"]))
    end

    for baseline in characterization["allocation_baselines"]
        path = joinpath(REPOSITORY_ROOT, baseline["contract"])
        @test file_sha256(path) == baseline["sha256"]
        contract = TOML.parsefile(path)
        allocation_values = recursive_key_values(
            contract,
            Set(String.(baseline["allocation_budget_keys"])),
        )
        @test Set(first.(allocation_values)) ==
            Set(String.(baseline["allocation_budget_keys"]))
        @test all(pair -> pair.second isa Real && pair.second >= 0,
            allocation_values)
        seed_values = recursive_key_values(
            contract,
            Set(String.(baseline["seed_keys"])),
        )
        @test Set(first.(seed_values)) ==
            Set(String.(baseline["seed_keys"]))
        @test !isempty(baseline["timer_key"])
    end

    for baseline in characterization["extension_baselines"]
        path = joinpath(REPOSITORY_ROOT, baseline["path"])
        @test isfile(path)
        if baseline["path"] in migrated_extension_tests
            @test file_sha256(path) != baseline["sha256"]
        else
            @test file_sha256(path) == baseline["sha256"]
        end
    end
    @test migrated_extension_tests ⊆ Set(String(
        baseline["path"]) for baseline in
        characterization["extension_baselines"])

    latency_index = characterization["package_latency"]
    latency_contract = TOML.parsefile(joinpath(
        REPOSITORY_ROOT,
        latency_index["contract"],
    ))
    @test latency_contract["fresh_process_runs"] >= 5
    @test latency_contract["contract"]["timer"] ==
        "Base.time_ns around each boundary; process launch and result serialization are excluded"
    @test latency_contract["workload"]["seed_policy"] ==
        "no random operation occurs"
    @test isfile(joinpath(REPOSITORY_ROOT, latency_index["harness"]))

    artifact = TOML.parsefile(joinpath(
        REPOSITORY_ROOT,
        latency_index["artifact"],
    ))
    @test artifact["correctness_passed"]
    @test !artifact["repository_dirty"]
    @test !artifact["package_source_dirty"]
    @test artifact["fresh_process_runs"] ==
        latency_contract["fresh_process_runs"]
    for key in (
        "load_samples_ns",
        "first_call_samples_ns",
        "warmed_call_samples_ns",
    )
        @test length(artifact[key]) == artifact["fresh_process_runs"]
        @test all(>(0), artifact[key])
    end

    latency_manifest = TOML.parsefile(joinpath(
        REPOSITORY_ROOT,
        latency_index["manifest"],
    ))
    @test latency_manifest["closure"]["status"] == "characterized"
    latency_entry = only(latency_manifest["artifacts"])
    @test latency_entry["correctness_passed"]
    @test latency_entry["sha256"] ==
        file_sha256(joinpath(REPOSITORY_ROOT, latency_index["artifact"]))
    @test latency_entry["authority_revision"] ==
        artifact["repository_head"]
end
