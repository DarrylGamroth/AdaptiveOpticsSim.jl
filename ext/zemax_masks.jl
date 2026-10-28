# FITS decoding belongs to the instrument artifact loader. Conversion emits
# ordinary cold graph bindings and records the selected physical source hashes.
function mask_number(value, label; positive=false)
    value isa Real && !(value isa Bool) && isfinite(value) && (!positive || value > 0) ||
        throw(ArgumentError("invalid mask $label"))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("mask $label is not representable"))
    return result
end

function mask_pair(value, label; positive=false)
    value isa AbstractVector && length(value) == 2 ||
        throw(ArgumentError("mask $label requires [row,column] values"))
    return Float64[mask_number(item,label;positive) for item in value]
end

function mask_artifact(root, reference)
    model_fields(reference,("path","sha256","hdu","axes","type","shape"),"FITS mask artifact")
    path = physical_reference(root,reference)
    model_integer(reference["hdu"],"FITS mask hdu";minimum=1)
    reference["type"] in ("Float32","Float64") || throw(ArgumentError("FITS masks require floating-point bindings"))
    shape = reference["shape"]
    shape isa AbstractVector && length(shape) == 2 || throw(ArgumentError("FITS mask shape must have rank two"))
    foreach(value -> model_integer(value,"FITS mask extent";minimum=1),shape)
    axes = reference["axes"]
    axes isa AbstractVector && length(axes) == 2 &&
        all(value -> value isa Int,axes) && Set(axes) == Set((1,2)) ||
        throw(ArgumentError("FITS mask axes must explicitly permute [1,2]"))
    return path
end

function sampled_mask_steps(root, instrument, settings, name)
    steps = Dict{String,Zemax.StepDocument}()
    sources = String[]
    choices = settings["masks"]
    for role in ("apodizer","focal_plane","lyot")
        haskey(choices,role) || continue
        choice = choices[role]
        choice isa AbstractDict || throw(ArgumentError("mask selection must be a mapping"))
        kind = get(choice,"kind",role == "apodizer" ? "unity" : "ideal")
        if kind != "fits"
            kind == (role == "apodizer" ? "unity" : "ideal") ||
                throw(ArgumentError("unsupported $role mask kind"))
            continue
        end
        model_fields(choice,("kind","mask"),"$role mask selection")
        id = model_string(choice["mask"],"physical mask name")
        masks = get(instrument,"optical_masks",nothing)
        masks isa AbstractDict && haskey(masks,id) || throw(ArgumentError("unknown physical optical mask: $id"))
        responses = masks[id]
        responses isa AbstractDict && !isempty(responses) &&
            all(port -> port in ("transmitted","reflected"),keys(responses)) ||
            throw(ArgumentError("invalid physical mask responses"))
        port = role == "lyot" ? instrument["optical_paths"][name]["lyot_stop_port"] : "transmitted"
        haskey(responses,port) || throw(ArgumentError("physical mask $id has no $port response"))
        response = responses[port]
        required = ("amplitude","sampling_m","center_pixel","outside_amplitude","outside_phase_rad")
        response isa AbstractDict && all(key -> key in (required...,"phase","wavelength_m"),keys(response)) &&
            all(key -> haskey(response,key),required) || throw(ArgumentError("invalid physical mask response fields"))
        haskey(response,"phase") && !haskey(response,"wavelength_m") &&
            throw(ArgumentError("phase masks require a reference wavelength_m"))
        if haskey(response,"wavelength_m")
            wavelength = mask_number(response["wavelength_m"],"wavelength_m";positive=true)
            isapprox(wavelength,settings["wavelength_um"]*1e-6;rtol=1e-12,atol=0) ||
                throw(ArgumentError("mask reference wavelength differs from the simulation wavelength"))
        end
        amplitude = response["amplitude"]
        push!(sources,relpath(mask_artifact(root,amplitude),root))
        outside = mask_number(response["outside_amplitude"],"outside_amplitude")
        0 <= outside <= 1 || throw(ArgumentError("outside_amplitude must be in [0,1]"))
        prefix = "$(role)_$(port)"
        selected = Zemax.step("sampled_mask","sampled_$role";
            amplitude="$(prefix)_amplitude",sampling_m=mask_pair(response["sampling_m"],"sampling_m";positive=true),
            center_pixel=mask_pair(response["center_pixel"],"center_pixel"),outside_amplitude=outside,
            outside_phase_rad=mask_number(response["outside_phase_rad"],"outside_phase_rad"))
        if haskey(response,"phase")
            phase = response["phase"]
            push!(sources,relpath(mask_artifact(root,phase),root))
            phase["shape"] == amplitude["shape"] || throw(ArgumentError("mask phase and amplitude shapes differ"))
            selected["phase"] = "$(prefix)_phase"
        end
        steps[role] = selected
    end
    return (;steps,sources)
end
