using DelimitedFiles
using KernelAbstractions

struct ReferenceCase
    id::String
    kind::Symbol
    data_path::String
    shape::Tuple{Vararg{Int}}
    atol::Float64
    rtol::Float64
    config::Dict{String,Any}
end

struct ReferenceBundle
    root::String
    cases::Vector{ReferenceCase}
end

default_reference_root() = get(ENV, "ADAPTIVEOPTICS_REFERENCE_ROOT", joinpath(@__DIR__, "reference_data"))

reference_manifest_path(root::AbstractString) = joinpath(root, "manifest.toml")

has_reference_bundle(root::AbstractString=default_reference_root()) = isfile(reference_manifest_path(root))

function parse_reference_case(id::AbstractString, raw::AbstractDict{<:AbstractString,<:Any}, root::AbstractString)
    kind = Symbol(get(raw, "kind", ""))
    data = get(raw, "data", nothing)
    shape = Tuple(Int.(get(raw, "shape", Int[])))
    if kind === Symbol("")
        throw(InvalidConfiguration("reference case '$id' is missing kind"))
    end
    if data === nothing
        throw(InvalidConfiguration("reference case '$id' is missing data path"))
    end
    if isempty(shape)
        throw(InvalidConfiguration("reference case '$id' is missing shape"))
    end
    return ReferenceCase(
        String(id),
        kind,
        joinpath(root, String(data)),
        shape,
        Float64(get(raw, "atol", 0.0)),
        Float64(get(raw, "rtol", 0.0)),
        raw,
    )
end

function load_reference_bundle(root::AbstractString=default_reference_root())
    manifest_path = reference_manifest_path(root)
    if !isfile(manifest_path)
        throw(InvalidConfiguration("reference manifest not found at $(manifest_path)"))
    end
    raw = TOML.parsefile(manifest_path)
    version = Int(get(raw, "version", 1))
    if version != 1
        throw(InvalidConfiguration("unsupported reference manifest version $(version)"))
    end
    raw_cases = get(raw, "cases", Dict{String,Any}())
    case_ids = sort!(collect(keys(raw_cases)))
    cases = ReferenceCase[]
    for id in case_ids
        push!(cases, parse_reference_case(id, raw_cases[id], root))
    end
    return ReferenceBundle(String(root), cases)
end

function load_reference_array(case::ReferenceCase)
    flat = vec(readdlm(case.data_path, Float64))
    n_expected = prod(case.shape)
    if length(flat) != n_expected
        throw(DimensionMismatchError("reference data for '$(case.id)' has $(length(flat)) values, expected $(n_expected)"))
    end
    return reshape(copy(flat), case.shape)
end

function parse_sensing_mode(name)
    lname = lowercase(String(name))
    if lname == "geometric"
        return Geometric()
    elseif lname == "diffractive"
        return Diffractive()
    end
    throw(InvalidConfiguration("unknown sensing mode '$name'"))
end

function build_reference_telescope(cfg::AbstractDict{<:AbstractString,<:Any})
    return Telescope(
        resolution=Int(cfg["resolution"]),
        diameter=Float64(cfg["diameter"]),
        sampling_time=Float64(cfg["sampling_time"]),
        central_obstruction=Float64(get(cfg, "central_obstruction", 0.0)),
    )
end

function build_reference_source(cfg::AbstractDict{<:AbstractString,<:Any})
    kind = lowercase(String(get(cfg, "kind", "ngs")))
    coordinates_raw = get(cfg, "coordinates", [0.0, 0.0])
    coordinates = (Float64(coordinates_raw[1]), Float64(coordinates_raw[2]))
    if kind == "ngs"
        return Source(
            band=Symbol(cfg["band"]),
            magnitude=Float64(get(cfg, "magnitude", 0.0)),
            coordinates=coordinates,
        )
    elseif kind == "lgs"
        na_profile = nothing
        if haskey(cfg, "na_profile_altitudes") && haskey(cfg, "na_profile_weights")
            altitudes = Float64.(cfg["na_profile_altitudes"])
            weights = Float64.(cfg["na_profile_weights"])
            na_profile = permutedims(hcat(altitudes, weights))
        end
        laser_raw = get(cfg, "laser_coordinates", [0.0, 0.0])
        laser_coordinates = (Float64(laser_raw[1]), Float64(laser_raw[2]))
        return LGSSource(
            elongation_factor=Float64(get(cfg, "elongation_factor", 1.0)),
            magnitude=Float64(get(cfg, "magnitude", 0.0)),
            coordinates=coordinates,
            laser_coordinates=laser_coordinates,
            fwhm_spot_up=Float64(get(cfg, "fwhm_spot_up", 0.0)),
            na_profile=na_profile,
        )
    end
    throw(InvalidConfiguration("unknown source kind '$kind'"))
end

function apply_reference_opd!(tel::Telescope, cfg::AbstractDict{<:AbstractString,<:Any})
    kind = lowercase(String(get(cfg, "kind", "zeros")))
    if kind == "zeros"
        reset_opd!(tel)
        return tel
    elseif kind == "ramp"
        scale_x = Float64(get(cfg, "scale_x", 0.0))
        scale_y = Float64(get(cfg, "scale_y", 0.0))
        bias = Float64(get(cfg, "bias", 0.0))
        @inbounds for j in axes(tel.state.opd, 2), i in axes(tel.state.opd, 1)
            tel.state.opd[i, j] = bias + scale_x * i + scale_y * j
        end
        return tel
    elseif kind == "constant"
        fill!(tel.state.opd, Float64(get(cfg, "value", 0.0)))
        return tel
    end
    throw(InvalidConfiguration("unknown OPD initializer '$kind'"))
end

function build_reference_wfs(kind::Symbol, cfg::AbstractDict{<:AbstractString,<:Any}, tel::Telescope)
    n_subap = Int(cfg["n_subap"])
    threshold = Float64(get(cfg, "threshold", 0.1))
    mode = parse_sensing_mode(get(cfg, "mode", "geometric"))
    if kind === :shack_hartmann_slopes
        pixel_scale = get(cfg, "pixel_scale", nothing)
        n_pix_subap = get(cfg, "n_pix_subap", nothing)
        if pixel_scale === nothing && n_pix_subap === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode)
        elseif n_pix_subap === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
                pixel_scale=Float64(pixel_scale))
        elseif pixel_scale === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
                n_pix_subap=Int(n_pix_subap))
        end
        return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
            pixel_scale=Float64(pixel_scale), n_pix_subap=Int(n_pix_subap))
    elseif kind === :pyramid_slopes
        return PyramidWFS(tel;
            n_subap=n_subap,
            threshold=threshold,
            modulation=Float64(get(cfg, "modulation", 0.0)),
            modulation_points=Int(get(cfg, "modulation_points", 1)),
            diffraction_padding=Int(get(cfg, "diffraction_padding", 2)),
            psf_centering=Bool(get(cfg, "psf_centering", true)),
            n_pix_separation=get(cfg, "n_pix_separation", nothing),
            n_pix_edge=get(cfg, "n_pix_edge", nothing),
            binning=Int(get(cfg, "binning", 1)),
            mode=mode,
        )
    elseif kind === :bioedge_slopes
        return BioEdgeWFS(tel;
            n_subap=n_subap,
            threshold=threshold,
            diffraction_padding=Int(get(cfg, "diffraction_padding", 2)),
            psf_centering=Bool(get(cfg, "psf_centering", true)),
            n_pix_separation=get(cfg, "n_pix_separation", nothing),
            n_pix_edge=get(cfg, "n_pix_edge", nothing),
            binning=Int(get(cfg, "binning", 1)),
            mode=mode,
        )
    end
    throw(InvalidConfiguration("unsupported WFS reference kind '$(kind)'"))
end

function compute_reference_actual(case::ReferenceCase)
    tel = build_reference_telescope(case.config["telescope"])
    if haskey(case.config, "opd")
        apply_reference_opd!(tel, case.config["opd"])
    end
    if case.kind === :psf
        src = build_reference_source(case.config["source"])
        zero_padding = Int(get(case.config["compute"], "zero_padding", 2))
        return copy(compute_psf!(tel, src; zero_padding=zero_padding))
    elseif case.kind in (:shack_hartmann_slopes, :pyramid_slopes, :bioedge_slopes)
        src = build_reference_source(case.config["source"])
        wfs = build_reference_wfs(case.kind, case.config["wfs"], tel)
        return copy(measure!(wfs, tel, src))
    end
    throw(InvalidConfiguration("unsupported reference case kind '$(case.kind)'"))
end

function compute_reference_actual_ka_cpu(case::ReferenceCase)
    tel = build_reference_telescope(case.config["telescope"])
    if haskey(case.config, "opd")
        apply_reference_opd!(tel, case.config["opd"])
    end
    if case.kind === :shack_hartmann_slopes
        src = build_reference_source(case.config["source"])
        wfs = build_reference_wfs(case.kind, case.config["wfs"], tel)
        mode_name = lowercase(String(get(case.config["wfs"], "mode", "geometric")))
        if mode_name != "geometric"
            throw(InvalidConfiguration("KA CPU reference path currently supports only geometric Shack-Hartmann cases"))
        end
        update_valid_mask!(wfs, tel)
        n_sub = wfs.params.n_subap
        sub = div(tel.params.resolution, n_sub)
        offset = n_sub * n_sub
        slopes = similar(wfs.state.slopes)
        style = AdaptiveOptics.AcceleratorStyle(KernelAbstractions.CPU())
        AdaptiveOptics._geometric_slopes!(style, slopes, tel.state.opd, wfs.state.valid_mask, sub, n_sub, offset)
        return slopes
    end
    throw(InvalidConfiguration("KA CPU reference path not implemented for reference kind '$(case.kind)'"))
end

function adapt_reference_actual(case::ReferenceCase, actual)
    compare = get(case.config, "compare", nothing)
    compare === nothing && return actual

    adapted = copy(actual)
    if get(compare, "swap_halves", false)
        if ndims(adapted) != 1 || isodd(length(adapted))
            throw(InvalidConfiguration("swap_halves compare option requires an even-length vector"))
        end
        n = div(length(adapted), 2)
        adapted = vcat(@view(adapted[n+1:end]), @view(adapted[1:n]))
    end
    if haskey(compare, "scale")
        adapted .*= Float64(compare["scale"])
    end
    return adapted
end

function validate_reference_case(case::ReferenceCase)
    expected = load_reference_array(case)
    actual = adapt_reference_actual(case, compute_reference_actual(case))
    if size(actual) != size(expected)
        return (ok=false, actual=actual, expected=expected, maxabs=Inf)
    end
    maxabs = maximum(abs, actual .- expected)
    ok = isapprox(actual, expected; atol=case.atol, rtol=case.rtol)
    return (ok=ok, actual=actual, expected=expected, maxabs=maxabs)
end

function write_reference_array(path::AbstractString, data)
    writedlm(path, vec(data))
    return path
end

function create_reference_fixture(root::AbstractString)
    mkpath(root)
    manifest = Dict{String,Any}(
        "version" => 1,
        "cases" => Dict{String,Any}(),
    )

    psf_case = Dict{String,Any}(
        "kind" => "psf",
        "data" => "psf_baseline.txt",
        "shape" => [64, 64],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "telescope" => Dict(
            "resolution" => 32,
            "diameter" => 8.0,
            "sampling_time" => 1e-3,
            "central_obstruction" => 0.2,
        ),
        "source" => Dict(
            "kind" => "ngs",
            "band" => "I",
            "magnitude" => 0.0,
        ),
        "compute" => Dict(
            "zero_padding" => 2,
        ),
    )
    psf_ref = compute_reference_actual(parse_reference_case("psf_baseline", psf_case, root))
    write_reference_array(joinpath(root, psf_case["data"]), psf_ref)
    manifest["cases"]["psf_baseline"] = psf_case

    sh_case = Dict{String,Any}(
        "kind" => "shack_hartmann_slopes",
        "data" => "shack_hartmann_slopes.txt",
        "shape" => [32],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "telescope" => Dict(
            "resolution" => 32,
            "diameter" => 8.0,
            "sampling_time" => 1e-3,
            "central_obstruction" => 0.0,
        ),
        "source" => Dict(
            "kind" => "ngs",
            "band" => "I",
            "magnitude" => 0.0,
        ),
        "opd" => Dict(
            "kind" => "ramp",
            "scale_x" => 1.0,
            "scale_y" => 0.0,
            "bias" => 0.0,
        ),
        "wfs" => Dict(
            "n_subap" => 4,
            "mode" => "geometric",
            "threshold" => 0.1,
        ),
    )
    sh_ref = compute_reference_actual(parse_reference_case("shack_hartmann_slopes", sh_case, root))
    write_reference_array(joinpath(root, sh_case["data"]), sh_ref)
    manifest["cases"]["shack_hartmann_slopes"] = sh_case

    open(reference_manifest_path(root), "w") do io
        TOML.print(io, manifest)
    end
    return root
end
