# Version 2 recipes describe one composable optical subsystem. Version 1 CAL2
# preparation remains owned by zemax_preparation.jl.
function component_model(path)
    filesize(path) <= 1024^2 || throw(ArgumentError("ZEMAX component recipe exceeds 1 MiB"))
    document = read_zemax_yaml(path)
    model_fields(document, ("schema_version", "catalogs", "conversion", "branches", "assumptions"), "ZEMAX component recipe")
    document["schema_version"] === 2 || throw(ArgumentError("unsupported ZEMAX component recipe schema_version"))
    conversion = document["conversion"]
    model_fields(conversion, ("first_surface", "initial_direction", "powered_mirrors"), "component conversion")
    first_surface = model_integer(conversion["first_surface"], "first_surface")
    initial_direction = conversion["initial_direction"]
    initial_direction isa Int && initial_direction in (-1, 1) ||
        throw(ArgumentError("initial_direction must be -1 or 1"))
    mirrors = conversion["powered_mirrors"]
    mirrors isa AbstractVector && length(mirrors) <= 256 ||
        throw(ArgumentError("powered_mirrors must contain at most 256 surface IDs"))
    powered_mirrors = Int[model_integer(id, "powered mirror surface"; minimum=first_surface) for id in mirrors]
    length(unique(powered_mirrors)) == length(powered_mirrors) ||
        throw(ArgumentError("duplicate powered mirror surface"))
    selections = document["branches"]
    selections isa AbstractDict && 1 <= length(selections) <= 16 ||
        throw(ArgumentError("component branches must contain 1 to 16 entries"))
    branches = Dict{String,@NamedTuple{configuration_id::Int,last_surface::Int}}()
    for (name, selection) in selections
        model_string(name, "component branch name")
        occursin(r"^[a-z][a-z0-9_-]*$", name) || throw(ArgumentError("component branch must be a lower-case name"))
        model_fields(selection, ("configuration_id", "last_surface"), "component branch $name")
        configuration_id = model_integer(selection["configuration_id"], "configuration_id"; minimum=1)
        last_surface = model_integer(selection["last_surface"], "last_surface"; minimum=first_surface)
        branches[name] = (; configuration_id, last_surface)
    end
    references = document["catalogs"]
    references isa AbstractVector && length(references) <= 64 ||
        throw(ArgumentError("catalogs must contain at most 64 references"))
    catalogs = Dict{String,String}[]
    for reference in references
        model_fields(reference, ("path", "sha256"), "catalog reference")
        push!(catalogs, Dict("path" => model_path(reference["path"], "catalog path"),
            "sha256" => model_string(reference["sha256"], "catalog sha256")))
    end
    length(unique(c["path"] for c in catalogs)) == length(catalogs) ||
        throw(ArgumentError("duplicate component catalog path"))
    authored_assumptions = document["assumptions"]
    authored_assumptions isa AbstractVector && length(authored_assumptions) <= 256 ||
        throw(ArgumentError("assumptions must contain at most 256 entries"))
    assumptions = String[model_string(s, "component assumption") for s in authored_assumptions]
    return (; first_surface, initial_direction, powered_mirrors, branches, catalogs, assumptions)
end

function component_catalogs(root, references)
    glasses = Dict{String,GlassDocument}()
    for reference in references
        path = physical_reference(root, reference)
        current = nothing
        for line in eachline(path)
            t = String.(split(strip(line)))
            isempty(t) && continue
            if t[1] == "NM"
                length(t) >= 4 || throw(ArgumentError("invalid AGF NM record in $path"))
                current = GlassDocument("catalog" => reference["path"], "name" => t[2],
                    "formula" => parse(Int, t[3]), "raw_records" => String[strip(line)])
                haskey(glasses, t[2]) || (glasses[t[2]] = current)
            elseif current !== nothing && t[1] in ("CD", "TD", "LD")
                haskey(current, t[1]) && throw(ArgumentError("duplicate AGF $(t[1]) record in $path"))
                current[t[1]] = parse.(Float64, t[2:end])
                push!(current["raw_records"], strip(line))
            end
        end
    end
    return glasses
end

# YAML sequences can decode as Vector{Any}, including mixed integer/float
# coefficients. Restore the same bounded storage used by the AGF parser.
function frozen_glass(document)
    return GlassDocument(key => (key == "raw_records" ? String[item for item in value] :
        key in ("CD", "TD", "LD") ? Float64[item for item in value] : value)
        for (key, value) in document)
end

function component_settings(settings, diameter_m)
    model_fields(settings, ("resolution", "pupil_resolution", "wavelength_um", "temperature_c"), "component settings")
    resolution = model_integer(settings["resolution"], "resolution"; minimum=1)
    iseven(resolution) || throw(ArgumentError("component resolution must be even"))
    pupil_resolution = model_integer(settings["pupil_resolution"], "pupil_resolution"; minimum=1)
    pupil_resolution <= resolution || throw(ArgumentError("pupil_resolution exceeds resolution"))
    for key in ("wavelength_um", "temperature_c")
        value = settings[key]
        value isa Real && !(value isa Bool) && isfinite(value) && isfinite(Float64(value)) ||
            throw(ArgumentError("$key must be finite"))
    end
    Float64(settings["wavelength_um"]) > 0 || throw(ArgumentError("wavelength_um must be positive"))
    diameter_m isa Real && !(diameter_m isa Bool) && isfinite(diameter_m) &&
        isfinite(Float64(diameter_m)) && Float64(diameter_m) > 0 ||
        throw(ArgumentError("component diameter_m must be positive and finite"))
    return Dict{String,Union{Int,Float64}}("resolution" => resolution,
        "pupil_resolution" => pupil_resolution, "wavelength_um" => Float64(settings["wavelength_um"]),
        "temperature_c" => Float64(settings["temperature_c"]))
end

function component_prescription(root, instrument)
    reference = get(instrument, "optical_prescription", nothing)
    reference isa AbstractDict && all(key -> key in ("format", "path", "sha256", "catalog_archive"), keys(reference)) &&
        all(key -> haskey(reference, key), ("format", "path", "sha256")) ||
        throw(ArgumentError("invalid optical_prescription reference"))
    reference["format"] == "zemax" || throw(ArgumentError("optical_prescription format must be zemax"))
    physical_reference(root, reference)
    if haskey(reference, "catalog_archive")
        model_fields(reference["catalog_archive"], ("path", "sha256"), "catalog_archive")
        physical_reference(root, reference["catalog_archive"])
    end
    return reference
end

# Authored graphs stay symbolic. RTC expands settings, profiles, responses and
# component imports after loading these exact raw copies.
function component_preparation_inputs(root, definition_path, instrument_path, definition, instrument;
        catalog_dir=nothing, frozen=false)
    catalog_dir === nothing || throw(ArgumentError("component recipes declare explicit hashed catalogs; catalog_dir is only for standalone preparation"))
    reference = component_prescription(root, instrument)
    source_path = physical_reference(root, reference)
    components = get(definition, "components", nothing)
    components isa AbstractDict && 1 <= length(components) <= 256 ||
        throw(ArgumentError("component preparation requires 1 to 256 authored components"))
    sources = Set([relpath(definition_path, root), relpath(instrument_path, root), relpath(source_path, root)])
    haskey(reference, "catalog_archive") && push!(sources, reference["catalog_archive"]["path"])
    ids = Int[]
    catalogs = Dict{String,Dict{String,GlassDocument}}()
    glass_path = optical_path(root, ZEMAX_GLASSES_PATH; generated=true)
    resolved_path = optical_path(root, ZEMAX_RESOLVED_PATH; generated=true)
    execution_path = optical_path(root, ZEMAX_EXECUTION_PATH; generated=true)
    frozen_catalog = frozen ? read_zemax_yaml(glass_path) : nothing
    if frozen
        model_fields(frozen_catalog, ("schema_version", "recipes"), "component glass snapshot")
        frozen_catalog["schema_version"] === 2 || throw(ArgumentError("unsupported component glass snapshot; regenerate the optical preparation"))
    end
    for name in sort!(collect(keys(components)))
        component = components[name]
        component isa AbstractDict || throw(ArgumentError("component $name must be a mapping"))
        haskey(component, "preparation") || continue
        preparation = component["preparation"]
        model_fields(preparation, ("kind", "recipe", "branch"), "component preparation")
        preparation["kind"] == "zemax" || throw(ArgumentError("unsupported component preparation kind"))
        recipe = model_path(preparation["recipe"], "component recipe")
        recipe_path = optical_path(root, recipe; generated=false)
        branch = model_string(preparation["branch"], "component branch")
        model = component_model(recipe_path)
        haskey(model.branches, branch) || throw(ArgumentError("unknown component branch $branch"))
        push!(ids, model.branches[branch].configuration_id)
        push!(sources, relpath(recipe_path, root))
        for catalog in model.catalogs
            push!(sources, relpath(physical_reference(root, catalog), root))
        end
        if !haskey(catalogs, recipe)
            catalogs[recipe] = frozen ? Dict{String,GlassDocument}(name => frozen_glass(glass)
                for (name, glass) in frozen_catalog["recipes"][recipe]) : component_catalogs(root, model.catalogs)
        end
    end
    isempty(ids) && throw(ArgumentError("authored graphs require at least one ZEMAX component preparation"))
    frozen && Set(keys(frozen_catalog["recipes"])) != Set(keys(catalogs)) &&
        throw(ArgumentError("component glass snapshot recipes changed; regenerate the optical preparation"))
    selected_ids = sort!(unique(ids))
    resolved = frozen ? read_zemax_yaml(resolved_path) : Zemax.inspect(source_path; configurations=selected_ids)
    get(resolved, "schema_version", nothing) === 1 &&
        get(resolved, "selected_configuration_ids", nothing) == selected_ids &&
        resolved["source"]["basename"] == basename(source_path) &&
        resolved["source"]["sha256"] == reference["sha256"] ||
        throw(ArgumentError("component resolved snapshot does not match the source recipe; regenerate the optical preparation"))
    # Verify selected branches exist before any generated files are replaced.
    for id in selected_ids
        count(branch -> branch["configuration_id"] == id, resolved["configurations"]) == 1 ||
            throw(ArgumentError("missing selected component configuration; regenerate the optical preparation"))
    end
    declarations = definition["graphs"]
    declarations isa AbstractDict && 1 <= length(declarations) <= 64 ||
        throw(ArgumentError("authored graphs must contain 1 to 64 modes"))
    graph_paths = Dict{String,String}()
    graphs = map(Tuple(sort!(collect(keys(declarations))))) do name
        model_string(name, "graph mode")
        occursin(r"^[a-z][a-z0-9_-]*$", name) && name != "simulation" ||
            throw(ArgumentError("graph modes must be lower-case filename stems distinct from simulation"))
        authored = declarations[name]
        graph = if authored isa AbstractString
            path = optical_path(root, model_path(authored, "authored graph"); generated=false)
            push!(sources, relpath(path, root))
            read_zemax_yaml(path)
        else
            authored
        end
        graph isa AbstractDict || throw(ArgumentError("authored graph must be a mapping or root-relative YAML path"))
        relative = "generated/$name.yaml"
        graph_paths[name] = relative
        (relative, ordered(graph), Dict{String,Matrix{Float64}}())
    end
    # Explicit file bindings are also immutable preparation inputs.
    for binding in values(get(definition, "bindings", Dict()))
        get(binding, "kind", nothing) == "fits" || continue
        push!(sources, relpath(physical_reference(root, binding), root))
    end
    execution = deepcopy(definition)
    execution["graphs"] = graph_paths
    execution["optical_model"] = Dict("format" => "zemax", "source" => "../simulation/definition.yaml")
    execution["instrument"] = relpath(instrument_path, dirname(execution_path))
    push!(sources, ZEMAX_RESOLVED_PATH, ZEMAX_GLASSES_PATH)
    sources = sort!(collect(sources))
    outputs = vcat([first(graph) for graph in graphs],
        [ZEMAX_RESOLVED_PATH, ZEMAX_GLASSES_PATH, ZEMAX_PROVENANCE_PATH, ZEMAX_EXECUTION_PATH])
    validate_optical_outputs(root, sources, outputs; regenerated=[ZEMAX_RESOLVED_PATH, ZEMAX_GLASSES_PATH])
    catalog = ordered(Dict("schema_version" => 2, "recipes" => catalogs))
    return (; root, catalog, sources, base_sources=sources, variant_paths=Set{String}(),
        graphs, glass_path, resolved_path, resolved=ordered(resolved),
        execution_path, execution=ordered(execution), write_catalog=true)
end

function Zemax.prepare_component(root::AbstractString, recipe::AbstractString, branch::AbstractString;
        settings, diameter_m, frozen=false)
    root = realpath(root)
    frozen && Zemax.validate_provenance(root)
    selected_settings = component_settings(settings, diameter_m)
    instrument_path = optical_path(root, "instrument.yaml"; generated=false)
    instrument = read_zemax_yaml(instrument_path)
    instrument isa AbstractDict || throw(ArgumentError("instrument must be a mapping"))
    reference = component_prescription(root, instrument)
    source_path = physical_reference(root, reference)
    recipe_path = optical_path(root, recipe; generated=false)
    model = component_model(recipe_path)
    haskey(model.branches, branch) || throw(ArgumentError("unknown component branch $branch"))
    selection = model.branches[branch]
    resolved = frozen ? read_zemax_yaml(optical_path(root, ZEMAX_RESOLVED_PATH; generated=true)) :
        Zemax.inspect(source_path; configurations=[selection.configuration_id])
    glasses = frozen ? read_zemax_yaml(optical_path(root, ZEMAX_GLASSES_PATH; generated=true))["recipes"][recipe] :
        component_catalogs(root, model.catalogs)
    selected = only(filter(b -> b["configuration_id"] == selection.configuration_id, resolved["configurations"]))
    prepared = Zemax.component_train(selected, model.first_surface,
        selection.last_surface, model.initial_direction, model.powered_mirrors,
        selected_settings, glasses, Float64(diameter_m))
    sources = Dict{String,String}(relpath(path, root) => hashfile(path)
        for path in (instrument_path, source_path, recipe_path))
    for catalog in model.catalogs
        sources[catalog["path"]] = catalog["sha256"]
    end
    if haskey(reference, "catalog_archive")
        archive = reference["catalog_archive"]
        model_fields(archive, ("path", "sha256"), "catalog_archive")
        archive_path = physical_reference(root, archive)
        sources[relpath(archive_path, root)] = hashfile(archive_path)
    end
    assumptions = vcat(model.assumptions, [
        "Centered scalar paraxial conversion uses saved evaluated ZMX values, not a Zemax ray trace or solve reevaluation.",
        "Coordinate-break translations along the axis are retained; transverse decenters, tilts, aperture offsets and conic/aspheric aberrations are omitted.",
        "Air and vacuum use index 1; glass propagation uses catalog index at the authored wavelength and temperature.",
        "Only CLAP circular clear apertures and manual DIAM mode 1 are applied; automatic footprints, SQAP and other aperture shapes are omitted.",
        "Powered mirrors use unfolded reduced-angle power -2 × incoming direction × index × curvature; every mirror reverses the axial direction.",
        "Propagation stops at the selected last surface before its outgoing thickness; external atmosphere, controllable mirrors, sensors and detectors remain authored graph components.",
    ])
    return (; prepared.config, prepared.diagnostics, sources, assumptions)
end
