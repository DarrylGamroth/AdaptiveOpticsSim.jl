default_gate0_reference_root() = get(ENV, "ADAPTIVEOPTICS_GATE0_REFERENCE_ROOT",
    joinpath(@__DIR__, "reference_data_gate0"))

has_gate0_reference_bundle(root::AbstractString=default_gate0_reference_root()) =
    has_reference_bundle(root)

load_gate0_reference_bundle(root::AbstractString=default_gate0_reference_root()) =
    load_reference_bundle(root)

function gate0_telescope(case::ReferenceCase)
    tel = build_reference_telescope(case.config["telescope"])
    if haskey(case.config, "opd")
        apply_reference_opd!(tel, case.config["opd"])
    end
    return tel
end

# Frozen Gate 0 optical references predate rate-based optical products. Keep
# their former telescope duration explicit in this comparison adapter only.
gate0_legacy_optical_duration(case::ReferenceCase) =
    Float64(case.config["telescope"]["sampling_time"])

function gate0_spatial_filter_shape(name)
    normalized = lowercase(String(name))
    normalized == "circular" && return CircularFilter()
    normalized == "square" && return SquareFilter()
    normalized == "foucault" && return FoucaultFilter()
    throw(InvalidConfiguration("unknown Gate 0 spatial-filter shape '$name'"))
end

function gate0_telescope_planes(case::ReferenceCase)
    tel = gate0_telescope(case)
    return cat(
        Float64.(pupil_mask(tel)),
        Array(pupil_reflectivity(tel)),
        Array(tel.state.opd);
        dims=3,
    )
end

function gate0_electric_field(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    wavefront = PupilFunction(tel; T=Float64)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=zero_padding,
        T=Float64)
    plan = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, plan)
    duration = gate0_legacy_optical_duration(case)
    amplitude_duration_scale = sqrt(duration)
    return cat(real.(field.values) .* amplitude_duration_scale,
        imag.(field.values) .* amplitude_duration_scale,
        abs2.(field.values) .* duration; dims=3)
end

function gate0_spatial_filter(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    cfg = case.config["spatial_filter"]
    sf = SpatialFilter(tel;
        shape=gate0_spatial_filter_shape(get(cfg, "shape", "circular")),
        diameter=Float64(cfg["diameter"]),
        zero_padding=Int(get(cfg, "zero_padding", 2)),
        T=Float64,
    )
    wavefront = PupilFunction(tel; T=Float64)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src;
        zero_padding=sf.params.zero_padding, T=Float64,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    formation = prepare_pupil_field(tel, wavefront, src, field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(field, wavefront, formation)
    output = PupilFunction(tel; T=Float64)
    plan = prepare_spatial_filter(tel, sf, field, output)
    workspace = SpatialFilterWorkspace(sf)
    filter!(output, field, sf, plan, workspace)
    phase = output.opd .* (2pi / wavelength(src))
    return cat(Array(phase), Array(output.amplitude); dims=3)
end

function gate0_optic_surface(case::ReferenceCase)
    tel = gate0_telescope(case)
    cfg = case.config["controllable_optic"]
    optic = build_reference_controllable_optic(cfg, tel)
    command = Float64.(cfg["command"])
    set_command!(optic, command)
    application = lowercase(String(get(cfg, "application", "additive")))
    mode = if application == "additive"
        DMAdditive()
    elseif application == "replace"
        DMReplace()
    else
        throw(InvalidConfiguration("unknown Gate 0 optic application '$application'"))
    end
    apply!(optic, tel, mode)
    return Array(tel.state.opd)
end

function gate0_radiometric_chain(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    detectors = [build_reference_detector(cfg) for cfg in
        case.config["detectors"]]
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    wavefront = PupilFunction(tel; T=Float64)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=zero_padding,
        T=Float64)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    photon_rate = pupil_photon_rate_map(tel, src)
    size(photon_rate) == size(field.values) || throw(DimensionMismatchError(
        "Gate 0 radiometric fixture requires zero_padding=1"))
    prepared = prepare_reference_direct_imaging(tel, src;
        zero_padding=zero_padding)
    rate_map = execute_reference_direct_imaging!(prepared, tel)
    image_rate = copy(intensity_values(rate_map))
    seed = Int(get(case.config["compute"], "seed", 1))
    acquisitions = [prepare_detector_acquisition(detector, rate_map)
        for detector in detectors]
    frames = [copy(capture!(detector, rate_map, acquisitions[index];
        rng=MersenneTwister(seed + index - 1)))
        for (index, detector) in enumerate(detectors)]
    duration = gate0_legacy_optical_duration(case)
    return cat(Array(photon_rate) .* duration,
        abs2.(field.values) .* duration, Array(image_rate) .* duration,
        (frame .* duration for frame in frames)...; dims=3)
end

function gate0_spectral_psf(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    cfg = case.config["spectrum"]
    bundle = SpectralBundle(Float64.(cfg["wavelengths"]),
        Float64.(cfg["weights"]); T=Float64)
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    spectral = with_spectrum(src, bundle)
    prepared = prepare_reference_direct_imaging(tel, spectral;
        zero_padding=zero_padding)
    products = execute_reference_direct_imaging!(prepared, tel)
    n = tel.params.resolution * zero_padding
    stack = Array{Float64}(undef, n, n, length(products) + 1)
    combined = zeros(Float64, n, n)
    for (index, product) in enumerate(products)
        plane = intensity_values(product)
        @views stack[:, :, index] .= plane
        # Frozen legacy characterization only. Production preserves these
        # distinct wavelength grids in `products` rather than index-summing.
        combined .+= plane
    end
    @views stack[:, :, end] .= combined
    stack .*= gate0_legacy_optical_duration(case)
    return stack
end
function gate0_direct_science(case::ReferenceCase)
    tel = gate0_telescope(case)
    source_cfgs = case.config["sources"]
    sources = [build_reference_source(cfg) for cfg in source_cfgs]
    length(sources) == 2 || throw(InvalidConfiguration(
        "Gate 0 direct-science fixture currently requires two sources"))
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    n = tel.params.resolution * zero_padding
    stack = Array{Float64}(undef, n, n, 5)
    first_prepared = prepare_reference_direct_imaging(tel, sources[1];
        zero_padding=zero_padding)
    second_prepared = prepare_reference_direct_imaging(tel, sources[2];
        zero_padding=zero_padding)
    @views stack[:, :, 1] .= intensity_values(
        execute_reference_direct_imaging!(first_prepared, tel))
    @views stack[:, :, 2] .= intensity_values(
        execute_reference_direct_imaging!(second_prepared, tel))
    combined_prepared = prepare_reference_direct_imaging(tel,
        Asterism(sources); zero_padding=zero_padding)
    combined = intensity_values(execute_reference_direct_imaging!(
        combined_prepared, tel))
    component_products = map(direct_imaging_output,
        direct_imaging_components(combined_prepared))
    @views begin
        stack[:, :, 3] .= intensity_values(component_products[1])
        stack[:, :, 4] .= intensity_values(component_products[2])
        stack[:, :, 5] .= combined
    end
    stack .*= gate0_legacy_optical_duration(case)
    return stack
end

function gate0_atmosphere_directions(case::ReferenceCase)
    tel = gate0_telescope(case)
    atmosphere = build_reference_atmosphere(case.config["atmosphere"], tel)
    seed = Int(get(case.config["atmosphere"], "advance_seed", 1))
    steps = Int(get(case.config["atmosphere"], "advance_steps", 1))
    duration = Float64(case.config["atmosphere"]["legacy_elapsed_duration"])
    step_duration = duration / steps
    rng = MersenneTwister(seed)
    for _ in 1:steps
        advance_by!(atmosphere, step_duration; rng=rng)
    end
    T = eltype(tel.state.opd)
    shifts = zeros(T, length(atmosphere.layers))
    scales = ones(T, length(atmosphere.layers))
    epoch_device = similar(tel.state.opd)
    AdaptiveOpticsSim.accumulate_rendered_layers!(epoch_device,
        atmosphere.layers, shifts, shifts, scales)
    epoch = copy(Array(epoch_device))
    onaxis = build_reference_source(case.config["onaxis_source"])
    offaxis = build_reference_source(case.config["offaxis_source"])
    lgs = build_reference_source(case.config["lgs_source"])
    propagate!(atmosphere, tel, onaxis)
    onaxis_opd = copy(Array(tel.state.opd))
    propagate!(atmosphere, tel, offaxis)
    offaxis_opd = copy(Array(tel.state.opd))
    propagate!(atmosphere, tel, lgs)
    lgs_opd = copy(Array(tel.state.opd))
    return cat(epoch, onaxis_opd, offaxis_opd, lgs_opd; dims=3)
end

function compute_gate0_actual(case::ReferenceCase)
    case.kind === :gate0_telescope_planes && return gate0_telescope_planes(case)
    case.kind === :gate0_electric_field && return gate0_electric_field(case)
    case.kind === :gate0_spatial_filter && return gate0_spatial_filter(case)
    case.kind === :gate0_optic_surface && return gate0_optic_surface(case)
    case.kind === :gate0_radiometric_chain && return gate0_radiometric_chain(case)
    case.kind === :gate0_spectral_psf && return gate0_spectral_psf(case)
    case.kind === :gate0_direct_science && return gate0_direct_science(case)
    case.kind === :gate0_atmosphere_directions &&
        return gate0_atmosphere_directions(case)
    throw(InvalidConfiguration("unsupported Gate 0 characterization kind '$(case.kind)'"))
end

function gate0_comparison_actual(case::ReferenceCase, actual)
    case.kind === :gate0_direct_science || return actual

    # The frozen pre-HIL fixture captured the retired row/column swap in the
    # legacy off-axis shift. Keep that historical artifact immutable and adapt
    # only its comparison view; production direct imaging follows the declared
    # (:x, :y) axis order and is tested independently in direct_science.jl.
    comparison = copy(actual)
    @views begin
        comparison[:, :, 2] .= permutedims(actual[:, :, 2], (2, 1))
        comparison[:, :, 4] .= permutedims(actual[:, :, 4], (2, 1))
        comparison[:, :, 5] .= comparison[:, :, 3] .+ comparison[:, :, 4]
    end
    return comparison
end

function validate_gate0_case(case::ReferenceCase)
    expected = load_reference_array(case)
    actual = compute_gate0_actual(case)
    if size(actual) != size(expected)
        return (ok=false, actual=actual, expected=expected, maxabs=Inf)
    end
    comparison = gate0_comparison_actual(case, actual)
    delta = abs.(comparison .- expected)
    maxabs = isempty(delta) ? 0.0 : maximum(delta)
    ok = isapprox(comparison, expected; atol=case.atol, rtol=case.rtol)
    return (ok=ok, actual=actual, expected=expected, maxabs=maxabs)
end
