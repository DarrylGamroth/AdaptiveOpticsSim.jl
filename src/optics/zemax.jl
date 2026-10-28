"""Bounded sequential Zemax inspection and paraxial optical conversion."""
module Zemax

using SHA
using LinearAlgebra

export inspect, write_inspection, prepare, refractive_index, train, validate_provenance

# These cold document builders retain the string-key mapping API without erasing
# stored values. Each document has a finite set of scalar and child shapes.
const SurfaceDocument = Dict{String,Union{Int,Float64,String,Vector{String},Dict{String,Float64}}}
const FieldDocument = Dict{String,Union{Int,Float64,String}}
const BranchDocument = Dict{String,Union{Int,String,Vector{Int},Vector{String},Vector{SurfaceDocument},Vector{FieldDocument}}}
const SourceDocument = Dict{String,Union{Int,String}}
const InspectionDocument = Dict{String,Union{Int,String,Vector{Int},SourceDocument,Vector{BranchDocument}}}
const CircleDocument = Dict{String,Union{Int,Float64}}
const StepDocument = Dict{String,Union{String,Float64,Vector{CircleDocument}}}
const TrainConfiguration = Dict{String,Union{Int,Float64,Vector{StepDocument}}}
const TrainNode = Dict{String,Union{String,TrainConfiguration}}
const PortDocument = Dict{String,String}
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

"""
    inspect(path; configurations=nothing)

Extract MM sequential surface snapshots and explicit multi-configuration overrides.
Only status-zero MCE values are accepted for selected configurations; pickups and
solves in surface records remain saved evaluated values with their raw records.
Omitted configurations select the saved configuration.
IGNM intervals include both ends and combine by union. Parameters remain in
native type-dependent units (for example mm and degrees for coordinate breaks).
"""
function inspect(path; configurations=nothing)
    source, bytes = read_zmx(path)
    lines = split(source, '\n')
    surfaces = Dict{Int,SurfaceDocument}()
    rows = String[]
    names = Dict{Int,String}()
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
        elseif t[1] == "SURF"
            id === nothing || (surfaces[id] = surface(id, records))
            id = parse(Int, t[2]); records = String[]
        elseif t[1] == "MNUM"
            id === nothing || (surfaces[id] = surface(id, records))
            id = nothing
            count, saved = parse.(Int, t[2:3]); in_mce = true
        elseif in_mce
            if t[1] == "MOFF"
                config = parse(Int, t[3])
                name = strip(t[4], '"')
                isempty(name) || (names[config] = name)
            else
                push!(rows, line)
            end
        elseif id !== nothing && startswith(raw, " ")
            push!(records, line)
        elseif id !== nothing
            # Surface definitions end when unindented system/editor records resume.
            surfaces[id] = surface(id, records); id = nothing
        end
    end
    id === nothing || (surfaces[id] = surface(id, records))
    unit_mm || error("Missing UNIT MM")
    count > 0 && 1 <= saved <= count || error("Missing or invalid MNUM")
    isempty(surfaces) && error("No surfaces")
    configurations = configurations === nothing ? [saved] : configurations
    length(unique(configurations)) == length(configurations) || error("Repeated configuration")
    all(c -> 1 <= c <= count, configurations) || error("Unknown configuration id")
    branches = BranchDocument[]
    for config in configurations
        resolved = deepcopy(surfaces)
        ignored = Set{Int}()
        fields = FieldDocument[]
        selected_rows = String[]
        for line in rows
            t = tokens(line)
            length(t) >= 6 || error("Malformed MCE record: $line")
            row_config = parse(Int, t[3])
            1 <= row_config <= count || error("Unknown MCE configuration: $line")
            row_config == config || continue
            op = t[1]
            op in ("IGNM", "THIC", "CRVT", "PRAM", "SDIA", "XFIE", "YFIE") ||
                error("Unsupported selected MCE operand: $op")
            parse(Int, t[5]) == 0 || error("Selected MCE pickup/solve unsupported: $line")
            push!(selected_rows, line)
            target = parse(Int, t[2]); value = number(t[4])
            if op in ("XFIE", "YFIE")
                push!(fields, FieldDocument("operand" => op, "field_id" => target, "value_native" => value))
                continue
            end
            haskey(resolved, target) || error("Unknown surface $target")
            if op == "IGNM"
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
        push!(branches, BranchDocument("configuration_id" => config,
            "name" => get(names, config, ""), "ignored_surface_ids" => sort!(collect(ignored)),
            "active_surfaces" => [resolved[i] for i in sort!(collect(keys(resolved))) if i ∉ ignored],
            "field_overrides_native" => fields, "raw_mce_records" => selected_rows))
    end
    return InspectionDocument("schema_version" => 1, "source" => SourceDocument("basename" => basename(path),
        "sha256" => bytes2hex(sha256(bytes)), "encoding" => "UTF-16LE with BOM",
        "units" => "MM", "surface_count" => length(surfaces)),
        "configuration_count" => count, "saved_configuration_id" => saved,
        "selected_configuration_ids" => Int[configurations...],
        "scope" => "Sequential ZMX inspection snapshot; includes housings and coordinate breaks; not a physical-optics model or ray trace.",
        "value_semantics" => "Saved evaluated surface values plus explicit status-zero MCE overrides. Surface solves/pickups and SQAP records are preserved, not reevaluated. Parameters and fields retain native ZMX units.",
        "configurations" => branches)
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

"""
    train(branch, name, selection, model, settings, masks, glasses, diameter; instrument_id="")

Convert an explicitly selected sequential snapshot into a centered scalar
paraxial optical graph document and diagnostic ABCD matrices. Coordinate breaks
retain signed axial translation; transverse offsets and tilts are omitted.
The authored conversion rules select idealized input-lens and OAP surfaces and
masks. Other surfaces must be STANDARD, COORDBRK, or zero-curvature TILTSURF;
TILTSURF tilts are omitted by the centered approximation. Curved mirrors outside the
selected OAP are rejected. Curvature is in inverse metres. Coordinate-break
parameters retain native units during inspection; tilts and offsets are omitted
by this centered conversion.
"""
function train(branch, name, selection, model, settings, masks, glasses, diameter; instrument_id="")
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
            push!(steps, step("multiply", "ideal_unity_apodizer"; amplitude=1.0))
        elseif id == rule["focal_mask_surface"]
            push!(steps, step("circular_obscuration", "ideal_focal_mask"; radius_m=masks["focal_plane"]["radius_m"]))
        elseif id == rule["lyot_surfaces"][name]
            push!(steps, step("circular_mask", "ideal_lyot_mask";
                mode=rule["lyot_modes"][name],
                circles=CircleDocument[CircleDocument(c) for c in masks["lyot"]["circles"]]))
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
        if id in (rule["focal_mask_surface"], rule["lyot_surfaces"][name], selection["detector_surface"])
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
    graph = TrainDocument("schema_version" => 1, "name" => isempty(instrument_id) ? name : "$(instrument_id)_$(name)",
        "nodes" => [TrainNode("type" => "proper_optical_train_f64", "name" => "optics", "config" => config)],
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
    prepare(root; catalog_dir=nothing, yamlfmt="yamlfmt", formatter_config=nothing)

Generate paraxial optical graph documents and a provenance receipt from an
instrument's authored simulation data. Load YAML.jl to activate preparation.
An extracted archive catalog directory is needed only to regenerate glass data.
"""
function prepare end

"""
    validate_provenance(root[, receipt])

Validate required authored sources, selected graphs, local paths, and exact
content hashes before preparing a generated optical model. Load YAML.jl to
activate this optional document validation.
"""
function validate_provenance end

end # module Zemax
