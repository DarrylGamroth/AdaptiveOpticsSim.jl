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

file_sha256(path::AbstractString) =
    bytes2hex(SHA.sha256(read(path)))

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

function plant_import_surface(path::AbstractString)
    text = read(path, String)
    block = split(split(
        text,
        "import ..AdaptiveOpticsSim:";
        limit=2,
    )[2], "\n\ninclude("; limit=2)[1]
    return unique!(strip.(filter(
        !isempty,
        split(replace(block, '\n' => ' '), ','),
    )))
end

function binding_occurs(text::AbstractString, binding::AbstractString)
    escaped = replace(binding, r"([\\.^$|?*+(){}\\[\\]])" => s"\\\1")
    return occursin(
        Regex("(?<![A-Za-z0-9_!])" * escaped * "(?![A-Za-z0-9_!])"),
        text,
    )
end

function adaptive_optics_sim_extension_hooks(path::AbstractString)
    pattern = r"(?m)^(?:@inline\s+)?(?:function\s+)?AdaptiveOpticsSim\.([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\{|\()"
    return sort!(unique!(String[
        match.captures[1] for match in eachmatch(pattern, read(path, String))
    ]))
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
    @test authority["schema_version"] == 1
    @test authority["name"] == "namespace_authority"
    @test authority["status"] == "authoritative"
    @test authority["implementation_phase"] == "pre_namespace_refactor"
    @test authority["semantics"]["compatibility"] ==
        "no root forwarding aliases, property forwarding, state views, or compatibility adapters"

    current = authority["current_root"]
    current_path = joinpath(REPOSITORY_ROOT, current["source"])
    @test file_sha256(current_path) == current["sha256"]
    declared = declared_surface(current_path)
    @test declared.exports == current["exports"]
    @test declared.public == current["public"]

    runtime = runtime_surface(AdaptiveOpticsSim)
    @test Set(runtime.exports) ==
        Set([current["exports"]; "AdaptiveOpticsSim"])
    @test Set(runtime.public) == Set(current["public"])

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
    @test isempty(migration["removed_exports"])
    @test isempty(migration["removed_public"])
    @test all(binding -> owner_by_binding[binding] == ("Root", "exports"),
        introduced_root)
    @test all(introduced_public) do qualified
        owner, binding = split(qualified, '.'; limit=2)
        owner_by_binding[binding] == (owner, "public")
    end

    existing_target = Set{String}()
    for (binding, (owner, visibility)) in owner_by_binding
        owner == "Plant" && continue
        binding in introduced_root && continue
        "$owner.$binding" in introduced_public && continue
        push!(existing_target, binding)
    end
    @test existing_target ==
        Set([current["exports"]; current["public"]])
    @test all(binding -> owner_by_binding[binding][2] == "exports",
        current["exports"])
    @test all(binding -> owner_by_binding[binding][2] == "public",
        current["public"])

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
        else
            @test all(binding -> !binding_occurs(plant_body, binding),
                entry["bindings"])
        end
    end
    @test length(imported_bindings) == length(unique(imported_bindings))
    @test Set(imported_bindings) ==
        Set(plant_import_surface(joinpath(
            REPOSITORY_ROOT,
            "src",
            "plant",
            "plant.jl",
        )))

    extension_files = authority["extension_files"]
    observed_hooks = Set{String}()
    for entry in extension_files
        path = joinpath(REPOSITORY_ROOT, entry["path"])
        hooks = adaptive_optics_sim_extension_hooks(path)
        @test hooks == sort!(String.(entry["adaptive_optics_sim_hooks"]))
        union!(observed_hooks, hooks)
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
    @test characterization["schema_version"] == 1
    @test characterization["status"] == "frozen_baseline_index"

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
        @test file_sha256(path) == baseline["sha256"]
    end

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
end
