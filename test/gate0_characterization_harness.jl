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
        Float64.(tel.state.pupil),
        Array(tel.state.pupil_reflectivity),
        Array(tel.state.opd);
        dims=3,
    )
end

function gate0_electric_field(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    n = tel.params.resolution * zero_padding
    field = Matrix{ComplexF64}(undef, n, n)
    fill_telescope_field!(field, tel, src; zero_padding=zero_padding)
    return cat(real.(field), imag.(field), abs2.(field); dims=3)
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
    phase, amplitude = filter!(sf, tel, src)
    return cat(Array(phase), Array(amplitude); dims=3)
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
    n = tel.params.resolution * zero_padding
    field = Matrix{ComplexF64}(undef, n, n)
    fill_telescope_field!(field, tel, src; zero_padding=zero_padding)
    flux = flux_map(tel, src)
    size(flux) == size(field) || throw(DimensionMismatchError(
        "Gate 0 radiometric fixture requires zero_padding=1"))
    psf = copy(compute_psf!(tel, src; zero_padding=zero_padding))
    seed = Int(get(case.config["compute"], "seed", 1))
    frames = [copy(capture!(detector, psf;
        rng=MersenneTwister(seed + index - 1)))
        for (index, detector) in enumerate(detectors)]
    return cat(Array(flux), abs2.(field), Array(psf), frames...; dims=3)
end

function gate0_spectral_psf(case::ReferenceCase)
    tel = gate0_telescope(case)
    src = build_reference_source(case.config["source"])
    cfg = case.config["spectrum"]
    bundle = SpectralBundle(Float64.(cfg["wavelengths"]),
        Float64.(cfg["weights"]); T=Float64)
    zero_padding = Int(get(case.config["compute"], "zero_padding", 1))
    n = tel.params.resolution * zero_padding
    stack = Array{Float64}(undef, n, n, length(bundle) + 1)
    combined = zeros(Float64, n, n)
    for (index, sample) in enumerate(bundle)
        sample_src = AdaptiveOpticsSim.source_with_wavelength_and_flux(src,
            sample.wavelength, photon_flux(src) * sample.weight)
        plane = compute_psf!(tel, sample_src; zero_padding=zero_padding)
        @views stack[:, :, index] .= plane
        combined .+= plane
    end
    @views stack[:, :, end] .= combined
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
    @views stack[:, :, 1] .= compute_psf!(tel, sources[1];
        zero_padding=zero_padding)
    @views stack[:, :, 2] .= compute_psf!(tel, sources[2];
        zero_padding=zero_padding)
    combined = copy(compute_psf!(tel, Asterism(sources);
        zero_padding=zero_padding))
    @views begin
        stack[:, :, 3] .= tel.state.psf_stack[:, :, 1]
        stack[:, :, 4] .= tel.state.psf_stack[:, :, 2]
        stack[:, :, 5] .= combined
    end
    return stack
end

function gate0_atmosphere_directions(case::ReferenceCase)
    tel = gate0_telescope(case)
    atmosphere = build_reference_atmosphere(case.config["atmosphere"], tel)
    seed = Int(get(case.config["atmosphere"], "advance_seed", 1))
    steps = Int(get(case.config["atmosphere"], "advance_steps", 1))
    rng = MersenneTwister(seed)
    for _ in 1:steps
        advance!(atmosphere, tel; rng=rng)
    end
    epoch = copy(Array(atmosphere.state.opd))
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

function validate_gate0_case(case::ReferenceCase)
    expected = load_reference_array(case)
    actual = compute_gate0_actual(case)
    if size(actual) != size(expected)
        return (ok=false, actual=actual, expected=expected, maxabs=Inf)
    end
    delta = abs.(actual .- expected)
    maxabs = isempty(delta) ? 0.0 : maximum(delta)
    ok = isapprox(actual, expected; atol=case.atol, rtol=case.rtol)
    return (ok=ok, actual=actual, expected=expected, maxabs=maxabs)
end
