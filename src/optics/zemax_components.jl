# Centered paraxial component conversion shares the native PROPER step schema
# with `train`, while leaving application graph ownership with the caller.

function component_aperture(s)
    records = s["raw_records"]
    clap = filter(r -> startswith(r, "CLAP "), records)
    diam = filter(r -> startswith(r, "DIAM "), records)
    if !isempty(clap)
        t = tokens(only(clap))
        length(t) == 4 && number(t[2]) == 0 ||
            throw(ArgumentError("unsupported circular aperture at surface $(s["id"])"))
        radius = number(t[3]) / 1000
    elseif !isempty(diam) && tokens(only(diam))[3] == "1"
        radius = s["semi_diameter_m"]
    else
        return 0.0
    end
    isfinite(radius) && radius >= 0 || throw(ArgumentError("invalid manual aperture radius"))
    return radius
end

function component_train(branch, first_surface, last_surface, initial_direction,
        powered_mirrors, settings, glasses, diameter_m; phase_screens=())
    all_surfaces = branch["active_surfaces"]
    ids = Int[s["id"] for s in all_surfaces]
    issorted(ids) && length(unique(ids)) == length(ids) ||
        throw(ArgumentError("active surfaces must be unique and ordered"))
    first_surface in ids && last_surface in ids ||
        throw(ArgumentError("component boundaries must be active surfaces"))
    surfaces = filter(s -> first_surface <= s["id"] <= last_surface, all_surfaces)
    # A component begins in air/vacuum. A source cut inside glass requires a
    # different explicit input-medium contract, which this recipe does not have.
    medium = ""
    for s in all_surfaces
        s["id"] >= first_surface && break
        s["type"] == "COORDBRK" && continue
        s["glass_after"] == "MIRROR" || (medium = s["glass_after"])
    end
    medium in ("", "VACUUM") ||
        throw(ArgumentError("component input medium must be air/vacuum; selected cut is inside $medium"))
    for id in filter(id -> id <= last_surface, powered_mirrors)
        matches = filter(s -> s["id"] == id, surfaces)
        length(matches) == 1 && only(matches)["glass_after"] == "MIRROR" &&
            only(matches)["curvature_inv_m"] != 0 ||
            throw(ArgumentError("powered_mirrors entry $id must select an active curved mirror"))
    end
    all(screen -> screen.surface in ids && first_surface <= screen.surface <= last_surface,
        phase_screens) || throw(ArgumentError("phase screens must select active surfaces inside the component"))
    steps = StepDocument[]
    diagnostics = Dict{String,Matrix{Float64}}()
    matrix = Matrix{Float64}(I, 2, 2)
    n = 1.0
    direction = Float64(initial_direction)
    for s in surfaces
        id, surface_type, glass = s["id"], s["type"], s["glass_after"]
        c = s["curvature_inv_m"]
        isfinite(c) && isfinite(s["thickness_m"]) ||
            throw(ArgumentError("nonfinite component geometry at surface $id"))
        (surface_type == "STANDARD" ||
            (surface_type in ("COORDBRK", "TILTSURF") && c == 0)) ||
            throw(ArgumentError("unsupported component surface type $surface_type at surface $id"))
        surface_type == "COORDBRK" && glass != "" &&
            throw(ArgumentError("coordinate break cannot change medium at surface $id"))
        phi = 0.0
        if glass == "MIRROR"
            c == 0 || id in powered_mirrors ||
                throw(ArgumentError("curved mirror $id requires explicit powered_mirrors admission"))
            # Unfolded reduced-angle convention: a concave mirror has positive
            # power for either incoming axial direction. The outgoing translation
            # changes sign at every reflection, including planar mirrors.
            phi = -2 * direction * n * c
            direction = -direction
        elseif surface_type != "COORDBRK"
            after = if glass in ("", "VACUUM")
                1.0
            else
                haskey(glasses, glass) || throw(ArgumentError("unknown component glass $glass at surface $id"))
                refractive_index(glasses[glass], settings["wavelength_um"], settings["temperature_c"])
            end
            isfinite(after) && after > 0 || throw(ArgumentError("invalid refractive index at surface $id"))
            phi = (after - n) * c * direction
            n = after
        end
        isfinite(phi) || throw(ArgumentError("nonfinite optical power at surface $id"))
        label = "surface_$id"
        if phi != 0
            push!(steps, step("lens", "$(label)_power"; focal_length_m=1 / phi))
            matrix = [1.0 0.0; -phi 1.0] * matrix
        end
        if surface_type != "COORDBRK"
            radius = component_aperture(s)
            radius > 0 && push!(steps, step("circular_aperture", "$(label)_clear_aperture"; radius_m=radius))
        end
        diagnostics[string(id)] = copy(matrix)
        for screen in phase_screens
            screen.surface == id || continue
            push!(steps, step("add_opd", "$(label)_$(screen.input)";
                input=screen.input, size=collect(screen.size), axes=collect(screen.axes),
                sampling_m=collect(screen.sampling_m), center_pixel=collect(screen.center_pixel),
                schema=screen.schema, outside_opd_m=screen.outside_opd_m))
        end
        id == last_surface && break
        distance = direction * s["thickness_m"] / n
        isfinite(distance) || throw(ArgumentError("nonfinite propagation distance at surface $id"))
        if distance != 0
            push!(steps, step("propagate", "$(label)_propagation"; distance_m=distance))
            matrix = [1.0 distance; 0.0 1.0] * matrix
        end
    end
    config = TrainConfiguration("resolution" => settings["resolution"],
        "pupil_resolution" => settings["pupil_resolution"],
        "wavelength_um" => settings["wavelength_um"], "diameter_m" => diameter_m,
        "steps" => steps)
    return (; config, diagnostics)
end
