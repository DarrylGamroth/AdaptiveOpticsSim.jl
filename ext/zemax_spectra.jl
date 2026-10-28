# Cold source quadrature and physical-filter evaluation. Each channel gets its
# own dispersive Zemax conversion; propagation remains the ordinary PROPER node.
function selected_settings(defaults, overrides)
    overrides isa AbstractDict || throw(ArgumentError("optical settings must be a mapping"))
    result = deepcopy(defaults)
    for (key, value) in overrides
        haskey(defaults, key) || throw(ArgumentError("unknown optical setting $key"))
        result[key] = value isa AbstractDict && defaults[key] isa AbstractDict ?
            selected_settings(defaults[key], value) : deepcopy(value)
    end
    return result
end

function filter_transmission(instrument, name, wavelength_m)
    name == "none" && return 1.0
    filters = get(instrument, "optical_filters", Dict())
    haskey(filters, name) || throw(ArgumentError("unknown optical filter $name"))
    curve = filters[name]
    model_fields(curve, ("wavelength_m", "transmission"), "optical filter $name")
    wavelengths, transmission = curve["wavelength_m"], curve["transmission"]
    wavelengths isa AbstractVector && transmission isa AbstractVector &&
        2 <= length(wavelengths) == length(transmission) <= 4096 ||
        throw(ArgumentError("optical filter requires matching wavelength/transmission samples"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x > 0, wavelengths) &&
        all(>(0), diff(wavelengths)) || throw(ArgumentError("filter wavelengths must increase strictly"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && 0 <= x <= 1, transmission) ||
        throw(ArgumentError("filter intensity transmission must be in [0,1]"))
    first(wavelengths) <= wavelength_m <= last(wavelengths) || return 0.0
    wavelength_m == last(wavelengths) && return Float64(last(transmission))
    index = searchsortedlast(wavelengths, wavelength_m)
    fraction = (wavelength_m-wavelengths[index])/(wavelengths[index+1]-wavelengths[index])
    return Float64((1-fraction)*transmission[index] + fraction*transmission[index+1])
end

function spectral_channels(settings, instrument, name)
    spectra = get(settings, "spectra", nothing)
    spectra === nothing && return nothing
    spectra isa AbstractDict && Set(keys(spectra)) == Set(keys(instrument["optical_paths"])) ||
        throw(ArgumentError("spectra must cover exactly the physical optical paths"))
    selected = spectra[name]
    model_fields(selected, ("wavelength_um", "fractional_bandwidth", "samples", "filter"), "spectrum $name")
    wavelength = mask_number(selected["wavelength_um"], "spectral wavelength_um"; positive=true)
    width = mask_number(selected["fractional_bandwidth"], "fractional_bandwidth")
    0 <= width < 2 || throw(ArgumentError("fractional_bandwidth must be in [0,2)"))
    count = model_integer(selected["samples"], "spectral samples"; minimum=1)
    count <= 64 && (width == 0 ? count == 1 : count >= 2) ||
        throw(ArgumentError("monochromatic spectra require one sample; finite bands require 2:64"))
    filter = model_string(selected["filter"], "spectral filter")
    allowed = get(instrument["optical_paths"][name], "filters", String[])
    (isempty(allowed) ? filter == "none" : filter in allowed) ||
        throw(ArgumentError("filter is not installed in optical path $name"))
    # Midpoint quadrature for a flat photon-number density per vacuum wavelength.
    # Filter loss is applied after source weights normalize, never renormalized.
    channels = [(;wavelength_um=wavelength*(1-width/2+width*(i-0.5)/count),
        weight=1.0/count) for i in 1:count]
    return (; wavelength_um=wavelength, filter, channels)
end

function spectral_train(root, branch, name, selection, model, settings, instrument, glasses, mode)
    spectral = spectral_channels(settings, instrument, name)
    reference = deepcopy(settings)
    spectral === nothing || (reference["wavelength_um"] = spectral.wavelength_um)
    masks = sampled_mask_steps(root, instrument, reference, name)
    graph, diagnostics = Zemax.train(branch, name, selection, model, reference, reference["masks"],
        glasses, instrument["optical_input"]["diameter_m"]; lyot_mode=mode,
        instrument_id=instrument["id"], mask_steps=masks.steps)
    spectral === nothing && return (;graph, diagnostics, sources=masks.sources)
    if length(spectral.channels) == 1 &&
            filter_transmission(instrument,spectral.filter,spectral.wavelength_um*1e-6) == 1
        return (;graph, diagnostics, sources=masks.sources)
    end
    channels = Zemax.ChannelDocument[]
    for channel in spectral.channels
        configured = deepcopy(reference)
        configured["wavelength_um"] = channel.wavelength_um
        # Single-wavelength measured phase data must reject other wavelengths.
        selected_masks = sampled_mask_steps(root, instrument, configured, name)
        optical, _ = Zemax.train(branch, name, selection, model, configured, configured["masks"],
            glasses, instrument["optical_input"]["diameter_m"]; lyot_mode=mode,
            instrument_id=instrument["id"], mask_steps=selected_masks.steps)
        push!(channels, Zemax.ChannelDocument("wavelength_um"=>channel.wavelength_um,
            "weight"=>channel.weight, "transmission"=>filter_transmission(instrument,
                spectral.filter, channel.wavelength_um*1e-6), "steps"=>only(optical["nodes"])["config"]["steps"]))
    end
    only(graph["nodes"])["config"]["channels"] = channels
    return (; graph, diagnostics, sources=masks.sources)
end
