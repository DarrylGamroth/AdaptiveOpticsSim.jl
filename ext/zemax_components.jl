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

function Zemax.prepare_component(root::AbstractString, recipe::AbstractString, branch::AbstractString;
        settings, diameter_m)
    root = realpath(root)
    selected_settings = component_settings(settings, diameter_m)
    instrument_path = optical_path(root, "instrument.yaml"; generated=false)
    instrument = read_zemax_yaml(instrument_path)
    instrument isa AbstractDict || throw(ArgumentError("instrument must be a mapping"))
    reference = get(instrument, "optical_prescription", nothing)
    reference isa AbstractDict && all(key -> key in ("format", "path", "sha256", "catalog_archive"), keys(reference)) &&
        all(key -> haskey(reference, key), ("format", "path", "sha256")) ||
        throw(ArgumentError("invalid optical_prescription reference"))
    reference["format"] == "zemax" || throw(ArgumentError("optical_prescription format must be zemax"))
    source_path = physical_reference(root, reference)
    recipe_path = optical_path(root, recipe; generated=false)
    model = component_model(recipe_path)
    haskey(model.branches, branch) || throw(ArgumentError("unknown component branch $branch"))
    selection = model.branches[branch]
    resolved = Zemax.inspect(source_path; configurations=[selection.configuration_id])
    glasses = component_catalogs(root, model.catalogs)
    prepared = Zemax.component_train(only(resolved["configurations"]), model.first_surface,
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
