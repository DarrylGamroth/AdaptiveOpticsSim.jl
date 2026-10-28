"""Bounded sequential Zemax inspection and paraxial optical conversion."""
module Zemax

using SHA
using LinearAlgebra

export inspect, write_inspection, load_model, prepare, execution_graphs, refractive_index, train, validate_provenance
export prepare_component

# These cold document builders retain the string-key mapping API without erasing
# stored values. Each document has a finite set of scalar and child shapes.
const SurfaceDocument = Dict{String,Union{Int,Float64,String,Vector{String},Dict{String,Float64}}}
const FieldDocument = Dict{String,Union{Int,Float64,String}}
const BranchDocument = Dict{String,Union{Int,String,Vector{Int},Vector{String},Vector{SurfaceDocument},Vector{FieldDocument},Dict{String,Float64}}}
const SourceDocument = Dict{String,Union{Int,String}}
const InspectionDocument = Dict{String,Union{Int,String,Vector{Int},SourceDocument,Vector{BranchDocument},Vector{String}}}
const CircleDocument = Dict{String,Union{Int,Float64}}
const StepDocument = Dict{String,Union{Int,String,Float64,Vector{Int},Vector{Float64},Vector{CircleDocument}}}
const ChannelDocument = Dict{String,Union{Float64,Vector{StepDocument}}}
const GridDocument = Dict{String,Union{Vector{Int},Vector{Float64}}}
const TrainConfiguration = Dict{String,Union{Int,Float64,Vector{StepDocument},Vector{ChannelDocument},GridDocument}}
const PortDocument = Dict{String,String}
const TrainNode = Dict{String,Union{String,TrainConfiguration,PortDocument}}
const TrainDocument = Dict{String,Union{Int,String,Vector{TrainNode},Vector{PortDocument}}}

# Bounded sequential-ZMX inspection: evaluated file values, never a Zemax solver.
number(s) = uppercase(s) == "INFINITY" ? Inf : parse(Float64, s)
tokens(s) = String[m.match for m in eachmatch(r"\"[^\"]*\"|\S+", s)]

function read_zmx(path)
    bytes = read(path)
    length(bytes) >= 2 && bytes[1:2] == [0xff, 0xfe] ||
        error("Expected UTF-16LE ZMX with BOM")
    iseven(length(bytes)) || error("Truncated UTF-16LE input")
    units = [UInt16(bytes[i]) | UInt16(bytes[i+1]) << 8 for i in 3:2:length(bytes)]
    return transcode(String, units), bytes
end

function surface(id, records)
    result = SurfaceDocument("id" => id, "comment" => "", "type" => "STANDARD",
        "curvature_inv_m" => 0.0, "thickness_m" => 0.0, "glass_after" => "",
        "semi_diameter_m" => 0.0, "parameters_native" => Dict{String,Float64}(),
        "raw_records" => records)
    for line in records
        t = tokens(line)
        isempty(t) && continue
        key = t[1]
        key == "COMM" && (result["comment"] = String(strip(line[5:end])))
        key == "TYPE" && (result["type"] = t[2])
        key == "CURV" && (result["curvature_inv_m"] = number(t[2]) * 1000)
        key == "DISZ" && (result["thickness_m"] = number(t[2]) / 1000)
        key == "GLAS" && (result["glass_after"] = t[2])
        key == "DIAM" && (result["semi_diameter_m"] = number(t[2]) / 1000)
        key == "PARM" && (result["parameters_native"][t[2]] = number(t[3]))
        # Preserve SQAP exactly: no aperture-shape or pickup evaluation is implied.
        key == "SQAP" && (result["square_aperture_native_record"] = line)
    end
    return result
end

# ZMX MCE status codes and pickup semantics follow SETMCOPERAND and
# Solve Types (multiple configuration editor) in the OpticStudio help:
# https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/SETMCOPERAND.html
# Row numbers include MOFF comments; each configuration has one cell per row.
function mce_cells(rows, count)
    cells = Dict{Tuple{Int,Int},Vector{String}}()
    ordinals = zeros(Int, count)
    for line in rows
        t = tokens(line)
        length(t) >= 5 || error("Malformed MCE record: $line")
        config = parse(Int, t[3])
        1 <= config <= count || error("Unknown MCE configuration: $line")
        ordinals[config] += 1
        cells[(ordinals[config], config)] = t
    end
    if any(t -> parse(Int, t[5]) in (2, 3), values(cells))
        all(==(first(ordinals)), ordinals) || error("MCE pickups require complete configuration rows")
        for row in 1:first(ordinals)
            definition = cells[(row, 1)][1:2]
            all(config -> cells[(row, config)][1:2] == definition, 1:count) ||
                error("Inconsistent MCE operand definition at row $row")
        end
    end
    return cells
end

function planar_air_thickness(target, surfaces, rows, system_records)
    # Thermal THIC in air normally depends on edge thickness and neighboring
    # radii. Restrict the zero-temperature-change identity to two planar air
    # surfaces, with one shared explicit environment and no local overrides.
    # See OpticStudio: Thickness (THIC) Values; Defining Multiple Temperature
    # and Pressure Values. Nonzero thermal analysis is outside this inspector.
    # https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Thickness_THIC_Values.html
    env = filter(line -> startswith(line, "ENVD "), system_records)
    length(env) == 1 || return false
    t = tokens(only(env))
    length(t) == 4 && number(t[4]) == 0 || return false
    all(isfinite, number.(t[2:3])) || return false
    any(line -> first(tokens(line)) in ("TEMP", "PRES", "TCEX", "CRVT", "PRAM"), rows) && return false
    # The edge-thickness proof requires these two actual neighboring surfaces.
    # Conservatively reject any configuration that can ignore either endpoint;
    # do not infer adjacency through an ignored-surface sequence.
    for line in rows
        r = tokens(line)
        r[1] in ("IGNR", "IGNM") || continue
        first_id = parse(Int, r[2])
        r[1] == "IGNM" && length(r) < 6 && return false
        last_id = r[1] == "IGNM" ? parse(Int, r[6]) : first_id
        if any(id -> first_id <= id <= last_id, (target, target + 1))
            parse(Int, r[5]) == 0 && number(r[4]) == 0 || return false
        end
    end
    for id in (target, target + 1)
        haskey(surfaces, id) || return false
        s = surfaces[id]
        s["type"] == "STANDARD" && s["curvature_inv_m"] == 0 &&
            s["glass_after"] in ("", "VACUUM") || return false
        any(line -> first(tokens(line)) in ("ENVD", "GLAS", "TCE", "PZUP", "PPAR"), s["raw_records"]) && return false
        any(line -> begin
            r = tokens(line)
            r[1] == "CURV" && length(r) >= 3 && parse(Int, r[3]) != 0
        end, s["raw_records"]) && return false
        any(line -> begin
            r = tokens(line)
            r[1] == "GLSS" && parse(Int, r[2]) == id && strip(r[4], '"') ∉ ("", "VACUUM")
        end, rows) && return false
    end
    return true
end

function mce_value(row, config, cells, surfaces, rows, system_records, visiting=Set{Tuple{Int,Int}}())
    key = (row, config)
    key ∈ visiting && error("Cyclic MCE pickup at row $row configuration $config")
    haskey(cells, key) || error("Unknown MCE pickup row $row configuration $config")
    t = cells[key]
    status = parse(Int, t[5])
    status == 0 && return number(t[4])
    status in (2, 3) || error("Selected MCE pickup/solve unsupported: $(join(t, ' '))")
    t[1] == "THIC" && length(t) == 12 && t[6:7] == ["0", "0"] && t[12] == "0" ||
        error("Only bounded THIC MCE pickups are supported: $(join(t, ' '))")
    target = parse(Int, t[2])
    status == 3 && !planar_air_thickness(target, surfaces, rows, system_records) &&
        error("Thermal THIC pickup requires shared fixed environment and planar air surfaces")
    from_config, from_row = parse.(Int, t[8:9])
    1 <= from_config <= config && 1 <= from_row <= row || error("Invalid MCE pickup source")
    haskey(cells, (from_row, from_config)) || error("Unknown MCE pickup source")
    cells[(from_row, from_config)][1] == "THIC" || error("MCE pickup source must be THIC")
    scale, offset = number.(t[10:11])
    all(isfinite, (scale, offset)) || error("Nonfinite MCE pickup scale/offset")
    push!(visiting, key)
    value = mce_value(from_row, from_config, cells, surfaces, rows, system_records, visiting) * scale + offset
    delete!(visiting, key)
    return value
end

"""
    inspect(path; configurations=nothing)

Extract MM sequential surface snapshots and explicit multi-configuration overrides.
Status-zero values and bounded THIC linear pickups are evaluated. Thermal THIC
pickups require a shared fixed environment and adjacent planar air surfaces;
other solves are rejected. Surface solves remain saved evaluated values.
Omitted configurations select the saved configuration. IGNM intervals retain
union semantics; IGNR sets an individual surface's status in editor order.
System aperture and field-type overrides remain native metadata: SATP 3 denotes
float by stop size, APER is the stop semi-diameter in mm, and FLTP 1 denotes
object height in mm. Surface semi-diameters do not define illuminated pupils.
Parameters remain in native type-dependent units (mm/degrees for coordinate breaks).
"""
function inspect(path; configurations=nothing)
    source, bytes = read_zmx(path)
    lines = split(source, '\n')
    surfaces = Dict{Int,SurfaceDocument}()
    rows = String[]
    names = Dict{Int,String}()
    system_records = String[]
    notes = String[]
    count = saved = 0
    id = nothing
    records = String[]
    in_mce = false
    unit_mm = false
    for raw in lines
        line = strip(raw)
        isempty(line) && continue
        t = tokens(line)
        if t[1] == "UNIT"
            t[2] == "MM" || error("Only MM ZMX units are supported")
            unit_mm = true
            push!(system_records, line)
        elseif t[1] == "SURF"
            id === nothing || (surfaces[id] = surface(id, records))
            id = parse(Int, t[2]); records = String[]
        elseif t[1] == "MNUM"
            id === nothing || (surfaces[id] = surface(id, records))
            id = nothing
            count, saved = parse.(Int, t[2:3]); in_mce = true
        elseif in_mce
            push!(rows, line)
            if t[1] == "MOFF"
                config = parse(Int, t[3])
                name = strip(t[4], '"')
                isempty(name) || (names[config] = name)
            end
        elseif id !== nothing && startswith(raw, " ")
            push!(records, line)
        elseif id !== nothing
            # Surface definitions end when unindented system/editor records resume.
            surfaces[id] = surface(id, records); id = nothing
            push!(system_records, line)
            t[1] == "NOTE" && push!(notes, line)
        else
            push!(system_records, line)
            t[1] == "NOTE" && push!(notes, line)
        end
    end
    id === nothing || (surfaces[id] = surface(id, records))
    unit_mm || error("Missing UNIT MM")
    count > 0 && 1 <= saved <= count || error("Missing or invalid MNUM")
    isempty(surfaces) && error("No surfaces")
    configurations = configurations === nothing ? [saved] : configurations
    length(unique(configurations)) == length(configurations) || error("Repeated configuration")
    all(c -> 1 <= c <= count, configurations) || error("Unknown configuration id")
    cells = mce_cells(rows, count)
    branches = BranchDocument[]
    for config in configurations
        resolved = deepcopy(surfaces)
        ignored = Set{Int}()
        fields = FieldDocument[]
        selected_rows = String[]
        system_overrides = Dict{String,Float64}()
        row = 0
        for line in rows
            t = tokens(line)
            length(t) >= 5 || error("Malformed MCE record: $line")
            row_config = parse(Int, t[3])
            1 <= row_config <= count || error("Unknown MCE configuration: $line")
            row_config == config || continue
            row += 1
            op = t[1]
            push!(selected_rows, line)
            if op == "MOFF"
                parse(Int, t[5]) == 0 || error("MOFF solve unsupported")
                continue
            end
            op in ("IGNM", "IGNR", "THIC", "CRVT", "PRAM", "SDIA", "XFIE", "YFIE", "GLSS", "SATP", "STPS", "APER", "FLTP") ||
                error("Unsupported selected MCE operand: $op")
            length(t) >= 6 || error("Malformed MCE record: $line")
            target = parse(Int, t[2])
            if op == "GLSS"
                parse(Int, t[5]) == 0 || error("GLSS pickup/solve unsupported")
                haskey(resolved, target) || error("Unknown surface $target")
                startswith(t[4], '"') && endswith(t[4], '"') || error("GLSS requires quoted glass")
                resolved[target]["glass_after"] = String(strip(t[4], '"'))
                continue
            end
            value = mce_value(row, config, cells, surfaces, rows, system_records)
            if op in ("SATP", "STPS", "APER", "FLTP")
                target == 0 || error("System MCE operand requires ignored target 0")
                op == "SATP" && value != 3 && error("Only SATP 3 (float by stop size) is supported")
                op == "FLTP" && value != 1 && error("Only FLTP 1 (object height) is supported")
                op == "APER" && !(isfinite(value) && value > 0) && error("Invalid APER")
                op == "STPS" && !(isinteger(value) && haskey(resolved, Int(value)) &&
                    minimum(keys(resolved)) < value < maximum(keys(resolved))) && error("Invalid STPS")
                system_overrides[op] = value
                continue
            end
            if op in ("XFIE", "YFIE")
                push!(fields, FieldDocument("operand" => op, "field_id" => target, "value_native" => value))
                continue
            end
            haskey(resolved, target) || error("Unknown surface $target")
            if op == "IGNR"
                value in (0, 1) || error("IGNR must be zero or one")
                value == 1 ? push!(ignored, target) : delete!(ignored, target)
            elseif op == "IGNM"
                last = parse(Int, t[6])
                target <= last && all(haskey(resolved, i) for i in target:last) ||
                    error("Invalid IGNM surface range")
                value in (0, 1) || error("IGNM must be zero or one")
                value == 1 && union!(ignored, target:last)
            elseif op == "PRAM"
                resolved[target]["parameters_native"][t[6]] = value
            else
                key, factor = op == "THIC" ? ("thickness_m", 0.001) :
                    op == "CRVT" ? ("curvature_inv_m", 1000.0) : ("semi_diameter_m", 0.001)
                resolved[target][key] = value * factor
            end
        end
        haskey(system_overrides, "APER") && get(system_overrides, "SATP", NaN) != 3 &&
            error("APER requires explicit SATP 3 (float by stop size)")
        branch = BranchDocument("configuration_id" => config,
            "name" => get(names, config, ""), "ignored_surface_ids" => sort!(collect(ignored)),
            "active_surfaces" => [resolved[i] for i in sort!(collect(keys(resolved))) if i ∉ ignored],
            "field_overrides_native" => fields, "raw_mce_records" => selected_rows)
        isempty(system_overrides) || (branch["system_overrides_native"] = system_overrides)
        push!(branches, branch)
    end
    result = InspectionDocument("schema_version" => 1, "source" => SourceDocument("basename" => basename(path),
        "sha256" => bytes2hex(sha256(bytes)), "encoding" => "UTF-16LE with BOM",
        "units" => "MM", "surface_count" => length(surfaces)),
        "configuration_count" => count, "saved_configuration_id" => saved,
        "selected_configuration_ids" => Int[configurations...],
        "scope" => "Sequential ZMX inspection snapshot; includes housings and coordinate breaks; not a physical-optics model or ray trace.",
        "value_semantics" => "Saved evaluated surface values plus explicit status-zero MCE overrides. Surface solves/pickups and SQAP records are preserved, not reevaluated. Parameters and fields retain native ZMX units.",
        "configurations" => branches)
    result["raw_system_records"] = system_records
    result["notes"] = notes
    if any(branch -> haskey(branch, "system_overrides_native"), branches)
        result["system_override_semantics"] = "SATP 3: float by stop size; STPS: stop surface id; APER: stop semi-diameter in mm; FLTP 1: object-height fields in mm. These records do not identify an illuminated pupil diameter."
    end
    if any(branch -> any(line -> parse(Int, tokens(line)[5]) in (2, 3), branch["raw_mce_records"]), branches)
        result["value_semantics"] = "Saved evaluated surface values plus explicit MCE overrides and bounded THIC pickups (scale × source + offset). Thermal THIC is evaluated only for a shared fixed environment and adjacent planar air surfaces. Surface solves/pickups remain saved values. System aperture and field-type metadata retain native ZMX units; surface semi-diameters are not illuminated pupil diameters."
    end
    return result
end

"""
    refractive_index(glass, wavelength_um, temperature_c)

Evaluate AGF Sellmeier 1 coefficients with the catalog TD thermal correction.
Wavelength is in micrometres and temperature in degrees Celsius. Air and vacuum
index conventions are the caller's responsibility.
"""
function refractive_index(glass, wavelength_um, temperature_c)
    glass["formula"] == 2 || error("Only Sellmeier 1 (AGF formula 2) is supported")
    c = glass["CD"]; lambda2 = wavelength_um^2
    haskey(glass, "LD") && !(glass["LD"][1] <= wavelength_um <= glass["LD"][2]) &&
        error("Wavelength outside catalog range for $(glass["name"])")
    n = sqrt(1 + sum(c[i] * lambda2 / (lambda2 - c[i+1]) for i in 1:2:5))
    if haskey(glass, "TD")
        d0, d1, d2, e0, e1, ltk, tref = glass["TD"]
        dt = temperature_c - tref
        n += (n^2 - 1) / (2n) * (d0*dt + d1*dt^2 + d2*dt^3 +
            (e0*dt + e1*dt^2) / (lambda2 - ltk^2))
    end
    return n
end

step(op, name; kwargs...) = StepDocument("op" => op, "name" => name,
    (string(k) => v for (k,v) in kwargs)...)

function ideal_focal_mask(masks)
    radius = masks["focal_plane"]["radius_m"]
    radius isa Real && !(radius isa Bool) && isfinite(radius) && radius > 0 ||
        throw(ArgumentError("ideal focal mask radius_m must be positive and finite"))
    return step("circular_obscuration", "ideal_focal_mask"; radius_m=radius)
end

function ideal_lyot_mask(masks, mode)
    circles = masks["lyot"]["circles"]
    circles isa AbstractVector && 1 <= length(circles) <= 16 ||
        throw(ArgumentError("ideal Lyot mask requires 1:16 circles"))
    for circle in circles
        circle isa AbstractDict && Set(keys(circle)) == Set(("radius_m","x_m","y_m")) ||
            throw(ArgumentError("invalid ideal Lyot circle fields"))
        all(value -> value isa Real && !(value isa Bool) && isfinite(value), values(circle)) &&
            circle["radius_m"] > 0 || throw(ArgumentError("invalid ideal Lyot circle geometry"))
    end
    return step("circular_mask", "ideal_lyot_mask"; mode,
        circles=CircleDocument[CircleDocument(circle) for circle in circles])
end

"""
    train(branch, name, selection, model, settings, masks, glasses, diameter; lyot_mode, instrument_id="")

Convert an explicitly selected sequential snapshot into a centered scalar
paraxial optical graph document and diagnostic ABCD matrices. Coordinate breaks
retain signed axial translation; transverse offsets and tilts are omitted.
The authored conversion rules select idealized input-lens and OAP surfaces and
masks. `selection = model["branches"][name]` supplies the configuration, detector
and Lyot surface; `lyot_mode` explicitly selects `transmit` or `reflect`.
`prepare` maps the physical Lyot stop port and validates the complete selection.
Its `mask_steps` keyword supplies sampled-mask steps and their ordinary cold
construction bindings when physical FITS maps replace the ideal masks.
Other surfaces must be STANDARD, COORDBRK, or zero-curvature TILTSURF;
TILTSURF tilts are omitted by the centered approximation. Curved mirrors outside the
selected OAP are rejected. Curvature is in inverse metres. Coordinate-break
parameters retain native units during inspection; tilts and offsets are omitted
by this centered conversion.
"""
function train(branch, name, selection, model, settings, masks, glasses, diameter;
        lyot_mode, instrument_id="", mask_steps=Dict{String,StepDocument}())
    lyot_mode in ("transmit", "reflect") || throw(ArgumentError("invalid lyot_mode"))
    rule = model["conversion"]
    surfaces = filter(s -> rule["first_surface"] <= s["id"] <= selection["detector_surface"], branch["active_surfaces"])
    isempty(surfaces) && error("Empty selected train")
    last(surfaces)["id"] == selection["detector_surface"] || error("Detector is not active")
    steps = StepDocument[]; n = 1.0; direction = 1.0; matrix = Matrix{Float64}(I, 2, 2)
    diagnostics = Dict{String,Matrix{Float64}}()
    for s in surfaces
        id = s["id"]; label = "surface_$(id)"; glass = s["glass_after"]
        surface_type = s["type"]
        (surface_type in ("STANDARD", "COORDBRK") ||
            (surface_type == "TILTSURF" && s["curvature_inv_m"] == 0) ||
            id == rule["input_lens_surface"] || id == rule["oap_surface"]) ||
            error("Unsupported optical surface type $surface_type at surface $id")
        (glass != "MIRROR" || s["curvature_inv_m"] == 0 || id == rule["oap_surface"]) ||
            error("Unmodeled curved mirror at surface $id")
        if id == rule["apodizer_surface"]
            push!(steps, haskey(mask_steps,"apodizer") ? mask_steps["apodizer"] :
                step("multiply", "ideal_unity_apodizer"; amplitude=1.0))
        elseif id == rule["focal_mask_surface"]
            push!(steps, haskey(mask_steps,"focal_plane") ? mask_steps["focal_plane"] :
                ideal_focal_mask(masks))
        elseif id == selection["lyot_surface"]
            push!(steps, haskey(mask_steps,"lyot") ? mask_steps["lyot"] : ideal_lyot_mask(masks,lyot_mode))
        end
        phi = 0.0
        if id == rule["input_lens_surface"]
            f = s["parameters_native"][rule["input_focal_parameter"]] / 1000
            phi = 1/f
        elseif glass == "MIRROR"
            id == rule["oap_surface"] && (phi = 2s["curvature_inv_m"])
            direction *= -1
        elseif s["type"] != "COORDBRK"
            after = glass in ("", "VACUUM") ? 1.0 :
                refractive_index(glasses[glass], settings["wavelength_um"], settings["temperature_c"])
            phi = (after - n) * s["curvature_inv_m"] * direction
            n = after
        end
        if phi != 0
            push!(steps, step("lens", "$(label)_power"; focal_length_m=1/phi))
            matrix = [1.0 0.0; -phi 1.0] * matrix
        end
        if s["type"] != "COORDBRK" && id != rule["oap_surface"]
            clap = filter(r -> startswith(r, "CLAP "), s["raw_records"])
            diam = filter(r -> startswith(r, "DIAM "), s["raw_records"])
            radius = if !isempty(clap)
                parse(Float64, split(only(clap))[3]) / 1000
            elseif !isempty(diam) && split(only(diam))[3] == "1"
                s["semi_diameter_m"]
            else
                0.0
            end
            radius > 0 && push!(steps, step("circular_aperture", "$(label)_clear_aperture"; radius_m=radius))
        end
        if id in (rule["focal_mask_surface"], selection["lyot_surface"], selection["detector_surface"])
            diagnostics[string(id)] = copy(matrix)
        end
        id == selection["detector_surface"] && break
        distance = direction * s["thickness_m"] / n
        if distance != 0
            push!(steps, step("propagate", "$(label)_propagation"; distance_m=distance))
            matrix = [1.0 distance; 0.0 1.0] * matrix
        end
    end
    config = TrainConfiguration("resolution" => settings["resolution"], "pupil_resolution" => settings["pupil_resolution"],
        "wavelength_um" => settings["wavelength_um"], "diameter_m" => diameter, "steps" => steps)
    node = TrainNode("type" => "proper_optical_train_f64", "name" => "optics", "config" => config)
    if !isempty(mask_steps)
        bindings = PortDocument()
        for mask in values(mask_steps), field in ("amplitude","phase")
            haskey(mask,field) && (bindings[mask[field]] = mask[field])
        end
        node["bindings"] = bindings
    end
    graph = TrainDocument("schema_version" => 1, "name" => isempty(instrument_id) ? name : "$(instrument_id)_$(name)",
        "nodes" => [node],
        "inputs" => [PortDocument("name" => port, "binding" => port, "destination" => "optics.$port") for port in ("pupil_opd", "pupil_amplitude")],
        "outputs" => [PortDocument("name" => "intensity", "source" => "optics.output")])
    return graph, diagnostics
end

"""
    write_inspection(output, data; yamlfmt="yamlfmt", formatter_config=nothing)

Write a generated YAML inspection document and format it with standalone
`yamlfmt`. Load YAML.jl to activate this optional writer.
"""
function write_inspection end

"""
    load_model(path)

Read one authored version-1 ZEMAX YAML conversion recipe, rejecting duplicate keys,
unknown fields and invalid conversion or branch selections. Return a mapping
with bounded concrete value types. Load YAML.jl to activate this optional loader.
The source simulation definition selects this recipe by a relative path.
"""
function load_model end

"""
    prepare(root; catalog_dir=nothing, yamlfmt="yamlfmt", formatter_config=nothing)

Generate paraxial optical graphs, frozen ZMX inspection, a shared runtime
definition at `generated/simulation.yaml`, and a provenance receipt from an
instrument's authored simulation definition and selected conversion recipe.
All derived outputs belong under `generated/`. Load YAML.jl to activate
preparation. Supply extracted AGF catalogs through `catalog_dir` to generate glass
data; otherwise reuse the frozen `generated/optics/glasses.yaml` snapshot for
standalone recipes. Graph-bearing definitions with version-2 component recipes
use their explicit hashed catalogs and preserve raw authored graph documents.
Use `yamlfmt=nothing` to write ordered YAML without an external formatter.

Standalone simulation definitions may add `preparation.variants`, a mapping of at most
16 lower-case names (at most 64 bytes) to `(branch, settings)` mappings. Each
branch must name an existing physical optical path. Settings strictly override
declared simulation defaults, while conversion, spectral/filter routing and
physical mask-port selection retain that base branch. Variants produce separate
named graphs under `generated/` and are included in execution and provenance.
They do not modify the conversion recipe or the base graphs. For example,
`variants.scc_unfringed = {branch: scc, settings: {masks: {lyot: {circles: [...]}}}}`
declares explicit static mask geometry; no reference aperture is inferred.
"""
function prepare end

"""
    execution_graphs(root; settings=Dict())

Validate frozen Zemax preparation and replan optical graphs for selected cold
settings in memory. Authored files and generated snapshots remain unchanged.
Declared static variants reapply their authored overrides to those selected
settings and remain bound to their base physical branch.
Load YAML.jl to activate this preparation boundary.
"""
function execution_graphs end

"""
    validate_provenance(root[, receipt])

Validate the exact authored and frozen snapshot source set, selected graphs,
prepared execution definition, local paths, and content hashes. This reads frozen
artifacts without regenerating files or invoking formatters. Load YAML.jl to
activate this optional document validation.
"""
function validate_provenance end

"""
    prepare_component(root, recipe, branch; settings, diameter_m, frozen=false,
                      phase_screens=(), pupil_grid=nothing)

Resolve a selected ZEMAX subsystem into native PROPER node configuration for
composition inside an ordinary simulation graph. The version-2 recipe selects
surfaces, explicit catalog inputs and approximation assumptions. Source and
catalog hashes are checked on every cold preparation; no files are written.
With `frozen=true`, validate generated provenance and use its resolved ZMX and
glass snapshots while recomputing the configuration for the selected settings
and incident diameter. Missing or stale preparation raises an error.
Optional `phase_screens` declare runtime OPD inputs by active surface ID, source
`size`, `axes`, `diameter_m`, `schema` and `outside_opd_m`. Optional `pupil_grid`
declares entrance `size`, `axes` and `diameter_m`. Both derive native physical
sampling and fractional pixel centers; no sampled arrays enter this boundary.
Load YAML.jl to activate this optional boundary.
"""
function prepare_component end

include("zemax_components.jl")

end # module Zemax
