using DelimitedFiles
using FFTW
using KernelAbstractions
using LinearAlgebra
using Random

struct ReferenceCase
    id::String
    kind::Symbol
    baseline::Symbol
    data_path::String
    shape::Tuple{Vararg{Int}}
    atol::Float64
    rtol::Float64
    config::Dict{String,Any}
end

struct ReferenceBundle
    root::String
    metadata::Dict{String,Any}
    cases::Vector{ReferenceCase}
end

abstract type ReferenceStorageConvention end

struct JuliaColumnMajorStorage <: ReferenceStorageConvention end
struct NumPyRowMajorStorage <: ReferenceStorageConvention end

abstract type ReferenceCompareConvention end

struct IdentityCompareConvention <: ReferenceCompareConvention end

struct OOPAOGeometricSHSignal2DConvention{T<:AbstractFloat} <: ReferenceCompareConvention
    slope_scale::T
end

const OOPAO_GEOMETRIC_SH_SLOPE_SCALE = 7.5949367088607559e6

function parse_reference_storage_convention(raw)
    order = lowercase(String(raw))
    if order in ("f", "julia_column_major", "column_major")
        return JuliaColumnMajorStorage()
    elseif order in ("c", "numpy_row_major", "row_major")
        return NumPyRowMajorStorage()
    end
    throw(InvalidConfiguration("unsupported storage order '$raw'"))
end

function reshape_reference_data(data, shape::Tuple, ::JuliaColumnMajorStorage)
    return reshape(data, shape)
end

function reshape_reference_data(data, shape::Tuple, ::NumPyRowMajorStorage)
    if length(shape) == 1
        return reshape(data, shape)
    end
    reshaped = reshape(data, reverse(shape))
    return permutedims(reshaped, reverse(1:length(shape)))
end

default_reference_root() = get(ENV, "ADAPTIVEOPTICS_REFERENCE_ROOT", joinpath(@__DIR__, "reference_data"))
default_specula_reference_root() = get(ENV, "ADAPTIVEOPTICS_SPECULA_REFERENCE_ROOT", joinpath(@__DIR__, "reference_data_specula"))

reference_manifest_path(root::AbstractString) = joinpath(root, "manifest.toml")

has_reference_bundle(root::AbstractString=default_reference_root()) = isfile(reference_manifest_path(root))
has_specula_reference_bundle(root::AbstractString=default_specula_reference_root()) = isfile(reference_manifest_path(root))

function parse_reference_baseline(kind::Symbol, raw::AbstractDict{<:AbstractString,<:Any})
    if haskey(raw, "baseline")
        return Symbol(lowercase(String(raw["baseline"])))
    elseif kind in (
        :tomography_model_gamma,
        :tomography_model_cxx,
        :tomography_model_cox,
        :tomography_model_cnz,
        :tomography_model_reconstructor,
        :tomography_model_wavefront,
        :tomography_im_reconstructor,
        :tomography_im_wavefront,
        :tomography_model_command_matrix,
        :tomography_model_dm_commands,
    )
        return :pytomoao
    end
    return :oopao
end

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
        parse_reference_baseline(kind, raw),
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
    metadata = Dict{String,Any}(get(raw, "metadata", Dict{String,Any}()))
    raw_cases = get(raw, "cases", Dict{String,Any}())
    case_ids = sort!(collect(keys(raw_cases)))
    cases = ReferenceCase[]
    for id in case_ids
        push!(cases, parse_reference_case(id, raw_cases[id], root))
    end
    return ReferenceBundle(String(root), metadata, cases)
end

reference_cases(bundle::ReferenceBundle, baseline::Symbol) = [case for case in bundle.cases if case.baseline === baseline]

function load_reference_array(case::ReferenceCase)
    flat = vec(readdlm(case.data_path, Float64))
    n_expected = prod(case.shape)
    if length(flat) != n_expected
        throw(DimensionMismatchError("reference data for '$(case.id)' has $(length(flat)) values, expected $(n_expected)"))
    end
    data = copy(flat)
    storage = parse_reference_storage_convention(get(case.config, "storage_order", "F"))
    return reshape_reference_data(data, case.shape, storage)
end

function load_reference_aux_array(
    case::ReferenceCase,
    relpath::AbstractString,
    shape,
    storage_order::AbstractString="F",
)
    path = joinpath(dirname(case.data_path), relpath)
    flat = vec(readdlm(path, Float64))
    dims = Tuple(Int.(shape))
    n_expected = prod(dims)
    if length(flat) != n_expected
        throw(DimensionMismatchError("reference aux data '$relpath' for '$(case.id)' has $(length(flat)) values, expected $(n_expected)"))
    end
    data = copy(flat)
    storage = parse_reference_storage_convention(storage_order)
    return reshape_reference_data(data, dims, storage)
end

function load_case_residual_opd(case::ReferenceCase)
    compute_cfg = get(case.config, "compute", nothing)
    compute_cfg === nothing && return nothing
    if !haskey(compute_cfg, "residual_opd_data")
        return nothing
    end
    return load_reference_aux_array(
        case,
        compute_cfg["residual_opd_data"],
        compute_cfg["residual_shape"],
        get(compute_cfg, "residual_storage_order", "C"),
    )
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

function parse_wfs_normalization(name)
    lname = lowercase(String(name))
    if lname in ("mean_valid_flux", "mean_flux", "slopesmaps")
        return MeanValidFluxNormalization()
    elseif lname in ("incidence_flux", "slopesmaps_incidence_flux")
        return IncidenceFluxNormalization()
    end
    throw(InvalidConfiguration("unknown WFS normalization '$name'"))
end

function parse_slope_order(name)
    lname = lowercase(String(name))
    if lname == "simu"
        return SimulationSlopes()
    elseif lname == "keck"
        return InterleavedSlopes()
    elseif lname == "inverted"
        return InvertedSlopes()
    end
    throw(InvalidConfiguration("unknown slope order '$name'"))
end

function build_reference_telescope(cfg::AbstractDict{<:AbstractString,<:Any})
    return Telescope(
        resolution=Int(cfg["resolution"]),
        diameter=Float64(cfg["diameter"]),
        sampling_time=Float64(cfg["sampling_time"]),
        central_obstruction=Float64(get(cfg, "central_obstruction", 0.0)),
        fov_arcsec=Float64(get(cfg, "fov_arcsec", 0.0)),
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

function build_reference_measurement_source(cfg::AbstractDict{<:AbstractString,<:Any})
    src = build_reference_source(cfg)
    spectrum_cfg = get(cfg, "spectrum", nothing)
    spectrum_cfg === nothing && return src
    bundle = SpectralBundle(
        Float64.(spectrum_cfg["wavelengths"]),
        Float64.(spectrum_cfg["weights"]),
    )
    return with_spectrum(src, bundle)
end

function build_reference_detector(cfg::AbstractDict{<:AbstractString,<:Any})
    noise_name = lowercase(String(get(cfg, "noise", "none")))
    integration_time = Float64(get(cfg, "integration_time", 1.0))
    qe = Float64(get(cfg, "qe", 1.0))
    psf_sampling = Int(get(cfg, "psf_sampling", 1))
    binning = Int(get(cfg, "binning", 1))
    noise = if noise_name == "none"
        NoiseNone()
    elseif noise_name == "photon"
        NoisePhoton()
    elseif noise_name == "readout"
        NoiseReadout(Float64(get(cfg, "readout_sigma", 1.0)))
    elseif noise_name == "photon_readout"
        NoisePhotonReadout(Float64(get(cfg, "readout_sigma", 1.0)))
    else
        throw(InvalidConfiguration("unknown detector noise kind '$noise_name'"))
    end
    return Detector(noise; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning)
end

function build_reference_atmosphere(cfg::AbstractDict{<:AbstractString,<:Any}, tel::Telescope)
    kind = lowercase(String(get(cfg, "kind", "finite")))
    kwargs = (
        r0=Float64(cfg["r0"]),
        L0=Float64(get(cfg, "L0", 25.0)),
        fractional_cn2=Float64.(cfg["fractional_cn2"]),
        wind_speed=Float64.(cfg["wind_speed"]),
        wind_direction=Float64.(cfg["wind_direction"]),
        altitude=Float64.(cfg["altitude"]),
        T=Float64,
    )
    if kind == "finite"
        return MultiLayerAtmosphere(tel; kwargs...)
    elseif kind == "infinite"
        return InfiniteMultiLayerAtmosphere(tel;
            kwargs...,
            screen_resolution=Int(get(cfg, "screen_resolution", default_infinite_screen_resolution(tel.params.resolution))),
            stencil_size=Int(get(cfg, "stencil_size", default_infinite_stencil_size(tel.params.resolution))),
        )
    end
    throw(InvalidConfiguration("unknown atmosphere kind '$kind'"))
end

function build_reference_atmospheric_model(cfg::AbstractDict{<:AbstractString,<:Any})
    kind = lowercase(String(get(cfg, "kind", "geometric")))
    chromatic_reference_wavelength = get(cfg, "chromatic_reference_wavelength", nothing)
    if kind == "geometric"
        return GeometricAtmosphericPropagation(
            chromatic_reference_wavelength=chromatic_reference_wavelength,
            T=Float64,
        )
    elseif kind in ("fresnel", "layered_fresnel")
        return LayeredFresnelAtmosphericPropagation(
            band_limit_factor=Float64(get(cfg, "band_limit_factor", 1.0)),
            chromatic_reference_wavelength=chromatic_reference_wavelength,
            T=Float64,
        )
    end
    throw(InvalidConfiguration("unknown atmospheric propagation model '$kind'"))
end

function reference_lift_img_resolution(
    tel::Telescope,
    det::AbstractDetector,
    compute_cfg::AbstractDict{<:AbstractString,<:Any},
)
    return Int(get(compute_cfg, "img_resolution", det.params.psf_sampling * tel.params.resolution))
end

function build_reference_basis(cfg::AbstractDict{<:AbstractString,<:Any}, tel::Telescope)
    kind = lowercase(String(get(cfg, "kind", "cartesian_polynomials")))
    n_modes = Int(get(cfg, "n_modes", 4))
    if kind == "cartesian_polynomials"
        n = tel.params.resolution
        basis = zeros(Float64, n, n, n_modes)
        x = collect(range(-1.0, 1.0; length=n + 1))[1:n]
        y = collect(range(-1.0, 1.0; length=n + 1))[1:n]
        @inbounds for j in 1:n, i in 1:n
            px = x[i]
            py = y[j]
            pupil = tel.state.pupil[i, j] ? 1.0 : 0.0
            if n_modes >= 1
                basis[i, j, 1] = pupil * px
            end
            if n_modes >= 2
                basis[i, j, 2] = pupil * py
            end
            if n_modes >= 3
                basis[i, j, 3] = pupil * px * py
            end
            if n_modes >= 4
                basis[i, j, 4] = pupil * (px * px - py * py)
            end
        end
        return basis
    end
    throw(InvalidConfiguration("unknown basis kind '$kind'"))
end

function build_reference_controllable_optic(cfg::AbstractDict{<:AbstractString,<:Any}, tel::Telescope)
    kind = lowercase(String(get(cfg, "kind", "modal")))
    label = Symbol(get(cfg, "label", "modal_optic"))
    T = Float64
    if kind == "tiptilt"
        scale = Float64(get(cfg, "scale", 1.0))
        return TipTiltMirror(tel; scale=scale, T=T, backend=backend(tel), label=label)
    elseif kind == "dm"
        n_act = Int(get(cfg, "n_act", 0))
        n_act > 0 || throw(InvalidConfiguration("dm controllable optic requires n_act"))
        influence_width = Float64(get(cfg, "influence_width", 0.2))
        return DeformableMirror(tel; n_act=n_act, influence_width=influence_width, T=T, backend=backend(tel))
    elseif kind == "modal"
        basis_name = lowercase(String(get(cfg, "basis", "cartesian_tilt")))
        if basis_name == "cartesian_tilt"
            scale = Float64(get(cfg, "scale", 1.0))
            return ModalControllableOptic(tel, CartesianTiltBasis(; scale=scale);
                labels=label, T=T, backend=backend(tel))
        elseif basis_name == "zernike"
            mode_ids = Int.(get(cfg, "mode_ids", Int[]))
            isempty(mode_ids) && throw(InvalidConfiguration("modal controllable optic zernike basis requires mode_ids"))
            scale = Float64(get(cfg, "scale", 1.0))
            return ModalControllableOptic(tel, ZernikeOpticBasis(mode_ids; scale=scale);
                labels=label, T=T, backend=backend(tel))
        end
        throw(InvalidConfiguration("unknown modal controllable optic basis '$basis_name'"))
    elseif kind == "composite"
        components = get(cfg, "components", nothing)
        components isa AbstractVector || throw(InvalidConfiguration("composite controllable optic requires components"))
        pairs = Pair{Symbol,AbstractControllableOptic}[]
        for raw_component in components
            component = raw_component isa AbstractDict{<:AbstractString,<:Any} ? raw_component :
                throw(InvalidConfiguration("composite controllable optic components must be tables"))
            component_label = Symbol(get(component, "label", "component"))
            push!(pairs, component_label => build_reference_controllable_optic(component, tel))
        end
        return CompositeControllableOptic(pairs...)
    end
    throw(InvalidConfiguration("unknown controllable optic kind '$kind'"))
end

function matrix_from_rows(rows)
    n_row = length(rows)
    n_col = n_row == 0 ? 0 : length(rows[1])
    mat = Matrix{Float64}(undef, n_row, n_col)
    @inbounds for i in 1:n_row
        row = rows[i]
        length(row) == n_col || throw(DimensionMismatchError("reference row lengths must match"))
        for j in 1:n_col
            mat[i, j] = Float64(row[j])
        end
    end
    return mat
end

function parse_reference_compare_convention(raw)
    raw === nothing && return IdentityCompareConvention()
    convention_name = lowercase(String(get(raw, "convention", "legacy")))
    if convention_name == "legacy"
        has_swap = get(raw, "swap_halves", false)
        has_scale = haskey(raw, "scale")
        if !has_swap && !has_scale
            return IdentityCompareConvention()
        elseif has_swap && has_scale
            return OOPAOGeometricSHSignal2DConvention(Float64(raw["scale"]))
        end
        throw(InvalidConfiguration("legacy compare adapters must specify both swap_halves=true and scale"))
    elseif convention_name == "oopao_geometric_sh_signal_2d"
        slope_scale = Float64(get(raw, "slope_scale", OOPAO_GEOMETRIC_SH_SLOPE_SCALE))
        return OOPAOGeometricSHSignal2DConvention(slope_scale)
    elseif convention_name == "identity"
        return IdentityCompareConvention()
    end
    throw(InvalidConfiguration("unknown compare convention '$convention_name'"))
end

adapt_compare_convention(::IdentityCompareConvention, actual) = actual

function adapt_compare_convention(
    convention::OOPAOGeometricSHSignal2DConvention,
    actual::AbstractVector,
)
    isodd(length(actual)) &&
        throw(InvalidConfiguration("OOPAO geometric SH compare convention requires an even-length slope vector"))
    n = div(length(actual), 2)
    adapted = similar(actual)
    @views adapted[1:n] .= actual[n+1:end]
    @views adapted[n+1:end] .= actual[1:n]
    adapted .*= convention.slope_scale
    return adapted
end

function adapt_compare_convention(
    ::OOPAOGeometricSHSignal2DConvention,
    actual,
)
    throw(InvalidConfiguration("OOPAO geometric SH compare convention requires a 1D slope vector"))
end

function parse_reference_real(value)
    if value isa AbstractString
        lowercase(String(value)) == "inf" && return Inf
    end
    return Float64(value)
end

function build_reference_tomography_atmosphere(cfg::AbstractDict{<:AbstractString,<:Any})
    return TomographyAtmosphereParams(
        zenith_angle_deg=Float64(cfg["zenith_angle_deg"]),
        altitude_km=Float64.(cfg["altitude_km"]),
        L0=Float64(cfg["L0"]),
        r0_zenith=Float64(cfg["r0_zenith"]),
        fractional_cn2=Float64.(cfg["fractional_cn2"]),
        wavelength=Float64(cfg["wavelength"]),
        wind_direction_deg=Float64.(cfg["wind_direction_deg"]),
        wind_speed=Float64.(cfg["wind_speed"]),
    )
end

function build_reference_tomography_asterism(cfg::AbstractDict{<:AbstractString,<:Any})
    return LGSAsterismParams(
        radius_arcsec=Float64(cfg["radius_arcsec"]),
        wavelength=Float64(cfg["wavelength"]),
        base_height_m=Float64(cfg["base_height_m"]),
        n_lgs=Int(cfg["n_lgs"]),
    )
end

function build_reference_tomography_wfs(cfg::AbstractDict{<:AbstractString,<:Any})
    valid_lenslet_map = convert.(Bool, matrix_from_rows(cfg["valid_lenslet_map"]))
    lenslet_offset = matrix_from_rows(cfg["lenslet_offset"])
    return LGSWFSParams(
        diameter=Float64(cfg["diameter"]),
        n_lenslet=Int(cfg["n_lenslet"]),
        n_px=Int(cfg["n_px"]),
        field_stop_size_arcsec=Float64(cfg["field_stop_size_arcsec"]),
        valid_lenslet_map=valid_lenslet_map,
        lenslet_rotation_rad=Float64.(cfg["lenslet_rotation_rad"]),
        lenslet_offset=lenslet_offset,
    )
end

function build_reference_tomography_params(cfg::AbstractDict{<:AbstractString,<:Any})
    return TomographyParams(
        n_fit_src=Int(cfg["n_fit_src"]),
        fov_optimization_arcsec=Float64(cfg["fov_optimization_arcsec"]),
        fit_src_height_m=parse_reference_real(get(cfg, "fit_src_height_m", Inf)),
    )
end

function build_reference_tomography_dm(cfg::AbstractDict{<:AbstractString,<:Any})
    valid_actuators = convert.(Bool, matrix_from_rows(cfg["valid_actuators"]))
    return TomographyDMParams(
        heights_m=Float64.(cfg["heights_m"]),
        pitch_m=Float64.(cfg["pitch_m"]),
        cross_coupling=Float64(cfg["cross_coupling"]),
        n_actuators=Int.(cfg["n_actuators"]),
        valid_actuators=valid_actuators,
    )
end

function row_major_mask_positions(mask::AbstractMatrix{Bool})
    positions = CartesianIndex{2}[]
    for i in axes(mask, 1), j in axes(mask, 2)
        mask[i, j] && push!(positions, CartesianIndex(i, j))
    end
    return positions
end

function pytomoao_phase_permutation(mask::AbstractMatrix{Bool})
    col_positions = findall(mask)
    col_index = Dict{CartesianIndex{2},Int}(idx => k for (k, idx) in enumerate(col_positions))
    return Int[col_index[idx] for idx in row_major_mask_positions(mask)]
end

function pytomoao_phase_block_permutation(mask::AbstractMatrix{Bool}, n_blocks::Integer)
    base = pytomoao_phase_permutation(mask)
    n_phase = length(base)
    perm = Vector{Int}(undef, n_blocks * n_phase)
    for block in 1:n_blocks
        offset = (block - 1) * n_phase
        @inbounds for k in eachindex(base)
            perm[offset + k] = offset + base[k]
        end
    end
    return perm
end

function pytomoao_model_mask(case::ReferenceCase)
    wfs = build_reference_tomography_wfs(case.config["wfs"])
    _, grid_mask = sparse_gradient_matrix(valid_lenslet_support(wfs); over_sampling=2)
    return grid_mask
end

function pytomoao_im_mask(case::ReferenceCase)
    dm = build_reference_tomography_dm(case.config["dm"])
    return dm.valid_actuators
end

function pytomoao_actuator_positions(mask::AbstractMatrix{Bool})
    positions = CartesianIndex{2}[]
    for i in axes(mask, 1), j in axes(mask, 2)
        mask[i, j] && push!(positions, CartesianIndex(i, j))
    end
    return positions
end

function pytomoao_actuator_permutation(dm::TomographyDMParams)
    support = dm_valid_support(dm)
    col_positions = findall(support)
    col_index = Dict{CartesianIndex{2},Int}(idx => k for (k, idx) in enumerate(col_positions))
    return Int[col_index[idx] for idx in pytomoao_actuator_positions(support)]
end

function adapt_actuator_rows_to_pytomoao(actual::AbstractMatrix, dm::TomographyDMParams)
    perm = pytomoao_actuator_permutation(dm)
    return actual[perm, :]
end

function adapt_actuator_vector_to_pytomoao(actual::AbstractVector, dm::TomographyDMParams)
    perm = pytomoao_actuator_permutation(dm)
    return actual[perm]
end

function build_reference_tomography_slopes(
    compute_cfg::AbstractDict{<:AbstractString,<:Any},
    wfs::LGSWFSParams,
    asterism::LGSAsterismParams,
)
    if haskey(compute_cfg, "slopes")
        return Float64.(compute_cfg["slopes"])
    end
    generator = lowercase(String(get(compute_cfg, "slopes_generator", "")))
    if generator == "pytomoao_tiptilt"
        n_valid = n_valid_subapertures(wfs)
        amp_pos = Float64(get(compute_cfg, "positive_amplitude", 4.0))
        amp_neg = Float64(get(compute_cfg, "negative_amplitude", -4.0))
        n_channels = Int(get(compute_cfg, "n_channels", asterism.n_lgs))
        base = zeros(Float64, 2 * n_valid)
        if n_valid > 1
            base[1:n_valid-1] .= amp_pos
        end
        base[n_valid+1:end] .= amp_neg
        return repeat(base, n_channels)
    end
    throw(InvalidConfiguration("reference tomography slopes are missing and no supported generator was provided"))
end

function adapt_phase_matrix_to_pytomoao(actual::AbstractMatrix, mask::AbstractMatrix{Bool})
    perm = pytomoao_phase_permutation(mask)
    return actual[perm, perm]
end

function adapt_phase_sensor_matrix_to_pytomoao(actual::AbstractMatrix, mask::AbstractMatrix{Bool}, n_blocks::Integer)
    row_perm = pytomoao_phase_permutation(mask)
    col_perm = pytomoao_phase_block_permutation(mask, n_blocks)
    return actual[row_perm, col_perm]
end

function adapt_reconstructor_rows_to_pytomoao(actual::AbstractMatrix, mask::AbstractMatrix{Bool})
    row_perm = pytomoao_phase_permutation(mask)
    return actual[row_perm, :]
end

function adapt_wavefront_map_to_pytomoao(actual::AbstractMatrix{T}, mask::AbstractMatrix{Bool}) where {T}
    values = actual[mask]
    row_positions = row_major_mask_positions(mask)
    adapted = fill(T(NaN), size(actual))
    @inbounds for (k, idx) in enumerate(row_positions)
        adapted[idx] = values[k]
    end
    return adapted
end

function combine_reference_modes!(dest::AbstractMatrix{T}, basis::AbstractArray{<:Real,3},
    coeffs::AbstractVector{<:Real}) where {T<:AbstractFloat}
    fill!(dest, zero(T))
    n_modes = min(size(basis, 3), length(coeffs))
    @inbounds for k in 1:n_modes
        coeff = T(coeffs[k])
        @views @. dest += coeff * basis[:, :, k]
    end
    return dest
end

function pupil_rms_nm(opd::AbstractMatrix{<:Real}, pupil::AbstractMatrix{Bool})
    acc = 0.0
    n = 0
    @inbounds for j in axes(opd, 2), i in axes(opd, 1)
        if pupil[i, j]
            val = Float64(opd[i, j])
            acc += val * val
            n += 1
        end
    end
    return 1e9 * sqrt(acc / n)
end

function center_crop(src::AbstractMatrix, size_out::Int)
    n, m = size(src)
    size_out <= min(n, m) || throw(DimensionMismatchError("crop size exceeds source size"))
    sx = fld(n - size_out, 2) + 1
    sy = fld(m - size_out, 2) + 1
    return @view src[sx:sx+size_out-1, sy:sy+size_out-1]
end

function strehl_ratio(psf::AbstractMatrix{T}, psf_ref::AbstractMatrix{T}) where {T<:AbstractFloat}
    size_min = min(size(psf, 1), size(psf_ref, 1))
    psf_crop = center_crop(psf, size_min)
    ref_crop = center_crop(psf_ref, size_min)
    otf = abs.(fft(psf_crop))
    otf_ref = abs.(fft(ref_crop))
    shifted = similar(otf)
    shifted_ref = similar(otf_ref)
    AdaptiveOpticsSim.fftshift2d!(shifted, otf)
    AdaptiveOpticsSim.fftshift2d!(shifted_ref, otf_ref)
    shifted ./= maximum(shifted)
    shifted_ref ./= maximum(shifted_ref)
    return 100.0 * sum(shifted) / sum(shifted_ref)
end

function reference_interaction_matrix(wfs::AbstractWFS, tel::Telescope, src::AbstractSource,
    basis::AbstractArray{<:Real,3}; amplitude::Real)
    n_modes = size(basis, 3)
    mat = nothing
    opd_base = copy(tel.state.opd)
    @inbounds for k in 1:n_modes
        @views @. tel.state.opd = Float64(amplitude) * basis[:, :, k]
        measure!(wfs, tel, src)
        if mat === nothing
            mat = Matrix{Float64}(undef, length(wfs.state.slopes), n_modes)
        end
        mat[:, k] .= wfs.state.slopes
    end
    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("reference interaction matrix requires at least one mode"))
    return mat
end

function closed_loop_trace(wfs::AbstractWFS, tel::Telescope, src::AbstractSource,
    control_basis::AbstractArray{<:Real,3}, forcing_coeffs::AbstractMatrix{<:Real};
    gain::Real, frame_delay::Int, calibration_amplitude::Real, psf_zero_padding::Int)
    n_iter, n_modes = size(forcing_coeffs)
    H = reference_interaction_matrix(wfs, tel, src, control_basis; amplitude=calibration_amplitude)
    recon = pinv(H)
    control_coeffs = zeros(Float64, n_modes)
    delayed_slopes = zeros(Float64, size(H, 1))
    forcing_opd = zeros(Float64, size(tel.state.opd))
    correction_opd = similar(forcing_opd)
    residual_opd = similar(forcing_opd)
    psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, src; zero_padding=psf_zero_padding))
    end
    trace = Matrix{Float64}(undef, n_iter, 5)

    for iter in 1:n_iter
        @views combine_reference_modes!(forcing_opd, control_basis, forcing_coeffs[iter, :])
        combine_reference_modes!(correction_opd, control_basis, control_coeffs)
        @. residual_opd = forcing_opd - correction_opd
        tel.state.opd .= residual_opd
        trace[iter, 1] = pupil_rms_nm(forcing_opd, tel.state.pupil)
        trace[iter, 2] = pupil_rms_nm(residual_opd, tel.state.pupil)
        psf = compute_psf!(tel, src; zero_padding=psf_zero_padding)
        trace[iter, 3] = strehl_ratio(psf, psf_ref)
        slopes = measure!(wfs, tel, src)
        trace[iter, 4] = norm(slopes)
        if frame_delay == 1
            delayed_slopes .= slopes
        end
        control_coeffs .+= gain .* (recon * delayed_slopes)
        trace[iter, 5] = norm(control_coeffs)
        if frame_delay == 2
            delayed_slopes .= slopes
        elseif frame_delay != 1
            throw(InvalidConfiguration("unsupported frame delay $(frame_delay)"))
        end
    end
    return trace
end

function gsc_closed_loop_trace(tel::Telescope, src::AbstractSource, wfs::PyramidWFS,
    control_basis::AbstractArray{<:Real,3}, forcing_coeffs::AbstractMatrix{<:Real};
    gain::Real, frame_delay::Int, calibration_amplitude::Real, psf_zero_padding::Int,
    og_floor::Real)
    n_iter, n_modes = size(forcing_coeffs)
    H = reference_interaction_matrix(wfs, tel, src, control_basis; amplitude=calibration_amplitude)
    recon = pinv(H)
    control_coeffs = zeros(Float64, n_modes)
    delayed_slopes = zeros(Float64, size(H, 1))
    forcing_opd = zeros(Float64, size(tel.state.opd))
    correction_opd = similar(forcing_opd)
    residual_opd = similar(forcing_opd)
    gsc = GainSensingCamera(wfs.state.pyramid_mask, control_basis)
    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, src)
    calibrate!(gsc, calibration_frame)
    frame = similar(calibration_frame)
    og_safe = similar(gsc.og)
    psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, src; zero_padding=psf_zero_padding))
    end
    trace = Matrix{Float64}(undef, n_iter, 6)

    for iter in 1:n_iter
        @views combine_reference_modes!(forcing_opd, control_basis, forcing_coeffs[iter, :])
        combine_reference_modes!(correction_opd, control_basis, control_coeffs)
        @. residual_opd = forcing_opd - correction_opd
        tel.state.opd .= residual_opd
        trace[iter, 1] = pupil_rms_nm(forcing_opd, tel.state.pupil)
        trace[iter, 2] = pupil_rms_nm(residual_opd, tel.state.pupil)
        psf = compute_psf!(tel, src; zero_padding=psf_zero_padding)
        trace[iter, 3] = strehl_ratio(psf, psf_ref)
        slopes = measure!(wfs, tel, src)
        trace[iter, 4] = norm(slopes)
        pyramid_modulation_frame!(frame, wfs, tel, src)
        og = compute_optical_gains!(gsc, frame)
        @. og_safe = max(abs(og), og_floor)
        trace[iter, 5] = sum(og_safe) / length(og_safe)
        if frame_delay == 1
            delayed_slopes .= slopes
        end
        control_coeffs .+= gain .* ((recon * delayed_slopes) ./ og_safe)
        trace[iter, 6] = norm(control_coeffs)
        if frame_delay == 2
            delayed_slopes .= slopes
        elseif frame_delay != 1
            throw(InvalidConfiguration("unsupported frame delay $(frame_delay)"))
        end
    end
    return trace
end

function gsc_atmosphere_replay_trace(
    tel::Telescope,
    ngs::AbstractSource,
    sci::AbstractSource,
    wfs::PyramidWFS,
    control_basis::AbstractArray{<:Real,3},
    forcing_ngs::AbstractArray{<:Real,3},
    forcing_src::AbstractArray{<:Real,3};
    gain::Real,
    frame_delay::Int,
    calibration_amplitude::Real,
    psf_zero_padding::Int,
    og_floor::Real,
)
    size(forcing_ngs) == size(forcing_src) ||
        throw(DimensionMismatchError("forcing_ngs and forcing_src must have matching shapes"))
    size(forcing_ngs, 1) == size(tel.state.opd, 1) == size(forcing_ngs, 2) == size(tel.state.opd, 2) ||
        throw(DimensionMismatchError("forcing OPD stack must match telescope resolution"))

    n_iter = size(forcing_ngs, 3)
    n_modes = size(control_basis, 3)
    H = reference_interaction_matrix(wfs, tel, ngs, control_basis; amplitude=calibration_amplitude)
    recon = pinv(H)
    control_coeffs = zeros(Float64, n_modes)
    delayed_slopes = zeros(Float64, size(H, 1))
    correction_opd = zeros(Float64, size(tel.state.opd))
    residual_ngs = similar(correction_opd)
    residual_src = similar(correction_opd)
    gsc = GainSensingCamera(wfs.state.pyramid_mask, control_basis)
    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, ngs)
    calibrate!(gsc, calibration_frame)
    frame = similar(calibration_frame)
    og_safe = similar(gsc.og)

    ngs_psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, ngs; zero_padding=psf_zero_padding))
    end
    sci_psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, sci; zero_padding=psf_zero_padding))
    end

    trace = Matrix{Float64}(undef, n_iter, 7)
    for iter in 1:n_iter
        forcing_ngs_i = @view forcing_ngs[:, :, iter]
        forcing_src_i = @view forcing_src[:, :, iter]
        combine_reference_modes!(correction_opd, control_basis, control_coeffs)
        @. residual_ngs = forcing_ngs_i - correction_opd
        @. residual_src = forcing_src_i - correction_opd

        apply_opd!(tel, residual_ngs)
        trace[iter, 1] = pupil_rms_nm(forcing_ngs_i, tel.state.pupil)
        trace[iter, 2] = pupil_rms_nm(residual_ngs, tel.state.pupil)
        ngs_psf = compute_psf!(tel, ngs; zero_padding=psf_zero_padding)
        trace[iter, 4] = strehl_ratio(ngs_psf, ngs_psf_ref)
        slopes = measure!(wfs, tel, ngs)
        trace[iter, 6] = norm(slopes)
        pyramid_modulation_frame!(frame, wfs, tel, ngs)
        og = compute_optical_gains!(gsc, frame)
        @. og_safe = max(abs(og), og_floor)

        apply_opd!(tel, residual_src)
        trace[iter, 3] = pupil_rms_nm(residual_src, tel.state.pupil)
        src_psf = compute_psf!(tel, sci; zero_padding=psf_zero_padding)
        trace[iter, 5] = strehl_ratio(src_psf, sci_psf_ref)

        if frame_delay == 1
            delayed_slopes .= slopes
        end
        control_coeffs .+= gain .* ((recon * delayed_slopes) ./ og_safe)
        trace[iter, 7] = sum(og_safe) / length(og_safe)
        if frame_delay == 2
            delayed_slopes .= slopes
        elseif frame_delay != 1
            throw(InvalidConfiguration("unsupported frame delay $(frame_delay)"))
        end
    end

    return trace
end

function transfer_functions(freq::AbstractVector{T}, loop_gain::T, Ti::T, Tau::T, Tdm::T) where {T<:AbstractFloat}
    S = complex.(zero(T), T(2π)) .* freq
    H_wfs = exp.(-Ti * S / 2)
    H_rtc = exp.(-Tau * S)
    H_dm = exp.(-Tdm * S)
    H_dac = (one(T) .- exp.(-S * Ti)) ./ (S * Ti)
    CC = loop_gain ./ (one(T) .- exp.(-S * Ti))
    H_ol = H_wfs .* H_rtc .* H_dac .* H_dm .* CC
    H_cl = H_ol ./ (one(T) .+ H_ol)
    H_er = one.(H_ol) ./ (one.(H_ol) .+ H_ol)
    H_n = H_cl ./ H_wfs
    return H_er, H_cl, H_ol, H_n
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
            tel.state.opd[i, j] = bias + scale_x * (i - 1) + scale_y * (j - 1)
        end
        return tel
    elseif kind == "constant"
        fill!(tel.state.opd, Float64(get(cfg, "value", 0.0)))
        return tel
    end
    throw(InvalidConfiguration("unknown OPD initializer '$kind'"))
end

function apply_reference_opd!(tel::Telescope, cfg::AbstractDict{<:AbstractString,<:Any},
    basis::AbstractArray{<:Real,3})
    kind = lowercase(String(get(cfg, "kind", "zeros")))
    if kind == "basis_mode"
        mode_index = Int(get(cfg, "mode_index", 1))
        amplitude = Float64(get(cfg, "amplitude", 0.0))
        if !(1 <= mode_index <= size(basis, 3))
            throw(InvalidConfiguration("basis mode index $(mode_index) is out of bounds"))
        end
        @views @. tel.state.opd = amplitude * basis[:, :, mode_index]
        return tel
    end
    return apply_reference_opd!(tel, cfg)
end

function build_reference_wfs(kind::Symbol, cfg::AbstractDict{<:AbstractString,<:Any}, tel::Telescope)
    n_subap = Int(cfg["n_subap"])
    threshold = Float64(get(cfg, "threshold", 0.1))
    mode = parse_sensing_mode(get(cfg, "mode", "geometric"))
    if kind in (:shack_hartmann_slopes, :shack_hartmann_frame)
        pixel_scale = get(cfg, "pixel_scale", nothing)
        n_pix_subap = get(cfg, "n_pix_subap", nothing)
        threshold_cog = Float64(get(cfg, "threshold_cog", 0.01))
        threshold_convolution = Float64(get(cfg, "threshold_convolution", 0.05))
        half_pixel_shift = Bool(get(cfg, "half_pixel_shift", false))
        if pixel_scale === nothing && n_pix_subap === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
                threshold_cog=threshold_cog, threshold_convolution=threshold_convolution,
                half_pixel_shift=half_pixel_shift)
        elseif n_pix_subap === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
                pixel_scale=Float64(pixel_scale), threshold_cog=threshold_cog,
                threshold_convolution=threshold_convolution, half_pixel_shift=half_pixel_shift)
        elseif pixel_scale === nothing
            return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
                n_pix_subap=Int(n_pix_subap), threshold_cog=threshold_cog,
                threshold_convolution=threshold_convolution, half_pixel_shift=half_pixel_shift)
        end
        return ShackHartmann(tel; n_subap=n_subap, threshold=threshold, mode=mode,
            pixel_scale=Float64(pixel_scale), n_pix_subap=Int(n_pix_subap),
            threshold_cog=threshold_cog, threshold_convolution=threshold_convolution,
            half_pixel_shift=half_pixel_shift)
    elseif kind in (:pyramid_slopes, :pyramid_frame)
        return PyramidWFS(tel;
            n_subap=n_subap,
            threshold=threshold,
            modulation=Float64(get(cfg, "modulation", 0.0)),
            light_ratio=Float64(get(cfg, "light_ratio", 0.0)),
            normalization=parse_wfs_normalization(get(cfg, "normalization", "mean_valid_flux")),
            modulation_points=get(cfg, "modulation_points", nothing),
            extra_modulation_factor=Int(get(cfg, "extra_modulation_factor", 0)),
            old_mask=Bool(get(cfg, "old_mask", false)),
            rooftop=Float64(get(cfg, "rooftop", 0.0)),
            theta_rotation=Float64(get(cfg, "theta_rotation", 0.0)),
            delta_theta=Float64(get(cfg, "delta_theta", 0.0)),
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
            modulation=Float64(get(cfg, "modulation", 0.0)),
            light_ratio=Float64(get(cfg, "light_ratio", 0.0)),
            normalization=parse_wfs_normalization(get(cfg, "normalization", "mean_valid_flux")),
            modulation_points=get(cfg, "modulation_points", nothing),
            extra_modulation_factor=Int(get(cfg, "extra_modulation_factor", 0)),
            delta_theta=Float64(get(cfg, "delta_theta", 0.0)),
            grey_width=Float64(get(cfg, "grey_width", 0.0)),
            grey_length=get(cfg, "grey_length", false),
            diffraction_padding=Int(get(cfg, "diffraction_padding", 2)),
            psf_centering=Bool(get(cfg, "psf_centering", true)),
            n_pix_separation=get(cfg, "n_pix_separation", nothing),
            n_pix_edge=get(cfg, "n_pix_edge", nothing),
            binning=Int(get(cfg, "binning", 1)),
            mode=mode,
        )
    elseif kind === :zernike_signal
        return ZernikeWFS(tel;
            n_subap=n_subap,
            threshold=threshold,
            phase_shift_pi=Float64(get(cfg, "phase_shift_pi", 0.5)),
            spot_radius_lambda_over_d=Float64(get(cfg, "spot_radius_lambda_over_d", 1.0)),
            normalization=parse_wfs_normalization(get(cfg, "normalization", "mean_valid_flux")),
            diffraction_padding=Int(get(cfg, "diffraction_padding", 2)),
            binning=Int(get(cfg, "binning", 1)),
        )
    elseif kind === :curvature_signal
        readout_name = lowercase(String(get(cfg, "readout_model", "frame")))
        readout_model = if readout_name == "frame"
            CurvatureFrameReadout()
        elseif readout_name == "counting"
            CurvatureCountingReadout()
        else
            throw(InvalidConfiguration("unknown curvature readout model '$readout_name'"))
        end
        branch_response = CurvatureBranchResponse(
            plus_throughput=Float64(get(cfg, "plus_throughput", 1.0)),
            minus_throughput=Float64(get(cfg, "minus_throughput", 1.0)),
            plus_background=Float64(get(cfg, "plus_background", 0.0)),
            minus_background=Float64(get(cfg, "minus_background", 0.0)),
        )
        return CurvatureWFS(tel;
            n_subap=n_subap,
            threshold=threshold,
            defocus_rms_nm=Float64(get(cfg, "defocus_rms_nm", 500.0)),
            diffraction_padding=Int(get(cfg, "diffraction_padding", 2)),
            readout_crop_resolution=Int(get(cfg, "readout_crop_resolution", tel.params.resolution)),
            readout_pixels_per_subap=Int(get(cfg, "readout_pixels_per_subap", 1)),
            readout_model=readout_model,
            branch_response=branch_response,
        )
    end
    throw(InvalidConfiguration("unsupported WFS reference kind '$(kind)'"))
end

function compute_reference_actual(case::ReferenceCase)
    if case.kind === :psf
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        if haskey(case.config, "opd")
            apply_reference_opd!(tel, case.config["opd"])
        end
        zero_padding = Int(get(case.config["compute"], "zero_padding", 2))
        return copy(compute_psf!(tel, src; zero_padding=zero_padding))
    elseif case.kind in (:shack_hartmann_slopes, :pyramid_slopes, :bioedge_slopes, :zernike_signal, :curvature_signal)
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_measurement_source(case.config["source"])
        wfs = build_reference_wfs(case.kind, case.config["wfs"], tel)
        residual = load_case_residual_opd(case)
        if residual !== nothing
            apply_opd!(tel, residual)
        elseif haskey(case.config, "controllable_optic")
            optic = build_reference_controllable_optic(case.config["controllable_optic"], tel)
            cmd = Float64.(get(case.config["controllable_optic"], "command", Float64[]))
            isempty(cmd) && throw(InvalidConfiguration("controllable_optic reference cases require a non-empty command"))
            set_command!(optic, cmd)
            apply!(optic, tel, DMReplace())
        elseif haskey(case.config, "opd")
            if haskey(case.config, "basis")
                basis = build_reference_basis(case.config["basis"], tel)
                apply_reference_opd!(tel, case.config["opd"], basis)
            else
                apply_reference_opd!(tel, case.config["opd"])
            end
        end
        return copy(measure!(wfs, tel, src))
    elseif case.kind === :shack_hartmann_frame
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_measurement_source(case.config["source"])
        wfs = build_reference_wfs(case.kind, case.config["wfs"], tel)
        if haskey(case.config, "opd")
            apply_reference_opd!(tel, case.config["opd"])
        end
        AdaptiveOpticsSim.prepare_sampling!(wfs, tel, src)
        AdaptiveOpticsSim.sampled_spots_peak!(wfs, tel, src)
        @views return copy(wfs.state.spot_cube[1, :, :])
    elseif case.kind === :pyramid_frame
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_measurement_source(case.config["source"])
        wfs = build_reference_wfs(case.kind, case.config["wfs"], tel)
        if haskey(case.config, "opd")
            apply_reference_opd!(tel, case.config["opd"])
        end
        if src isa SpectralSource
            AdaptiveOpticsSim.accumulate_pyramid_spectral_intensity!(AdaptiveOpticsSim.execution_style(wfs.state.intensity), wfs, tel, src)
        else
            AdaptiveOpticsSim.pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
        end
        return copy(AdaptiveOpticsSim.sample_pyramid_intensity!(wfs, tel, wfs.state.intensity))
    elseif case.kind === :atmospheric_intensity
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_measurement_source(case.config["source"])
        atm = build_reference_atmosphere(case.config["atmosphere"], tel)
        advance_steps = Int(get(case.config["atmosphere"], "advance_steps", 1))
        advance_seed = Int(get(case.config["atmosphere"], "advance_seed", 1))
        rng = MersenneTwister(advance_seed)
        for _ in 1:advance_steps
            advance!(atm, tel; rng=rng)
        end
        prop = AtmosphericFieldPropagation(
            atm,
            tel,
            src;
            model=build_reference_atmospheric_model(case.config["propagation"]),
            zero_padding=Int(get(case.config["propagation"], "zero_padding", 1)),
            T=Float64,
        )
        return copy(atmospheric_intensity!(prop, atm, tel, src))
    elseif case.kind === :gsc_optical_gains
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        wfs = build_reference_wfs(:pyramid_slopes, case.config["wfs"], tel)
        basis = build_reference_basis(case.config["basis"], tel)
        gsc = GainSensingCamera(wfs.state.pyramid_mask, basis)
        calibration_frame = similar(wfs.state.intensity)
        reset_opd!(tel)
        pyramid_modulation_frame!(calibration_frame, wfs, tel, src)
        calibrate!(gsc, calibration_frame)
        residual = load_case_residual_opd(case)
        if residual !== nothing
            apply_opd!(tel, residual)
        elseif haskey(case.config, "opd")
            apply_reference_opd!(tel, case.config["opd"], basis)
        end
        frame = similar(calibration_frame)
        pyramid_modulation_frame!(frame, wfs, tel, src)
        return copy(compute_optical_gains!(gsc, frame))
    elseif case.kind === :gsc_modulation_frame
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        wfs = build_reference_wfs(:pyramid_slopes, case.config["wfs"], tel)
        basis = build_reference_basis(case.config["basis"], tel)
        residual = load_case_residual_opd(case)
        if residual !== nothing
            apply_opd!(tel, residual)
        elseif haskey(case.config, "opd")
            apply_reference_opd!(tel, case.config["opd"], basis)
        end
        frame = similar(wfs.state.intensity)
        pyramid_modulation_frame!(frame, wfs, tel, src)
        return copy(frame)
    elseif case.kind === :transfer_function_rejection
        compute_cfg = case.config["compute"]
        fs = Float64(get(compute_cfg, "fs", 300.0))
        n_freq = Int(get(compute_cfg, "n_freq", 512))
        loop_gains = Float64.(get(compute_cfg, "loop_gains", [0.2, 0.6]))
        freq = collect(range(fs / n_freq, fs / 2; length=n_freq - 1))
        Ti = Float64(1 / fs)
        Tau = Ti / 2
        Tdm = Ti / 2
        rejection = Matrix{Float64}(undef, length(freq), length(loop_gains))
        for (idx, gain) in pairs(loop_gains)
            H_er, _, _, _ = transfer_functions(freq, gain, Ti, Tau, Tdm)
            rejection[:, idx] .= 20 .* log10.(abs.(H_er))
        end
        return rejection
    elseif case.kind in (
        :tomography_model_gamma,
        :tomography_model_cxx,
        :tomography_model_cox,
        :tomography_model_cnz,
        :tomography_model_reconstructor,
        :tomography_model_wavefront,
        :tomography_model_command_matrix,
        :tomography_model_dm_commands,
        :tomography_im_reconstructor,
        :tomography_im_wavefront,
    )
        atmosphere = build_reference_tomography_atmosphere(case.config["atmosphere"])
        asterism = build_reference_tomography_asterism(case.config["asterism"])
        wfs = build_reference_tomography_wfs(case.config["wfs"])
        tomography = build_reference_tomography_params(case.config["tomography"])
        dm = build_reference_tomography_dm(case.config["dm"])
        compute_cfg = get(case.config, "compute", Dict{String,Any}())
        if case.kind in (
            :tomography_model_gamma,
            :tomography_model_cxx,
            :tomography_model_cox,
            :tomography_model_cnz,
            :tomography_model_reconstructor,
            :tomography_model_wavefront,
            :tomography_model_command_matrix,
            :tomography_model_dm_commands,
        )
            recon = build_reconstructor(ModelBasedTomography(), atmosphere, asterism, wfs, tomography, dm)
            if case.kind === :tomography_model_gamma
                return Matrix(recon.operators.gamma)
            elseif case.kind === :tomography_model_cxx
                return Matrix(recon.operators.cxx)
            elseif case.kind === :tomography_model_cox
                return Matrix(recon.operators.cox)
            elseif case.kind === :tomography_model_cnz
                return Matrix(recon.operators.cnz)
            elseif case.kind === :tomography_model_reconstructor
                return Matrix(recon.reconstructor)
            elseif case.kind === :tomography_model_command_matrix || case.kind === :tomography_model_dm_commands
                assembly = get(compute_cfg, "assembly", Dict{String,Any}())
                command_recon = assemble_reconstructor_and_fitting(
                    recon,
                    dm;
                    n_channels=Int(get(assembly, "n_channels", asterism.n_lgs)),
                    slope_order=parse_slope_order(get(assembly, "slope_order", "simu")),
                    scaling_factor=Float64(get(assembly, "scaling_factor", 1.65e7)),
                )
                if haskey(assembly, "mask_actuators")
                    mask_actuators!(command_recon, Int.(assembly["mask_actuators"]))
                end
                if case.kind === :tomography_model_command_matrix
                    return Matrix(command_recon.matrix)
                end
                slopes = build_reference_tomography_slopes(compute_cfg, wfs, asterism)
                return dm_commands(command_recon, slopes)
            end
            slopes = build_reference_tomography_slopes(compute_cfg, wfs, asterism)
            return reconstruct_wavefront_map(recon, slopes)
        else
            imat = matrix_from_rows(compute_cfg["interaction_matrix"])
            recon = build_reconstructor(
                InteractionMatrixTomography(),
                imat,
                dm.valid_actuators,
                atmosphere,
                asterism,
                wfs,
                tomography,
                dm,
            )
            if case.kind === :tomography_im_reconstructor
                return Matrix(recon.reconstructor)
            end
            slopes = Float64.(compute_cfg["slopes"])
            return reconstruct_wavefront_map(recon, slopes)
        end
    elseif case.kind === :lift_interaction_matrix
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        basis = build_reference_basis(case.config["basis"], tel)
        det = build_reference_detector(case.config["detector"])
        diversity = zeros(Float64, size(tel.state.opd))
        if haskey(case.config, "diversity")
            apply_reference_opd!(tel, case.config["diversity"])
            diversity .= tel.state.opd
            reset_opd!(tel)
        end
        compute_cfg = case.config["compute"]
        img_resolution = reference_lift_img_resolution(tel, det, compute_cfg)
        numerical = Bool(get(compute_cfg, "numerical", false))
        mode_ids = Int.(get(compute_cfg, "mode_ids", collect(1:size(basis, 3))))
        flux_norm = Float64(get(compute_cfg, "flux_norm", 1.0))
        coeffs = Float64.(get(compute_cfg, "coefficients", zeros(length(mode_ids))))
        lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=Int(get(compute_cfg, "iterations", 3)),
            img_resolution=img_resolution, numerical=numerical)
        H = lift_interaction_matrix(lift, coeffs, mode_ids; flux_norm=flux_norm)
        stack = Array{Float64}(undef, img_resolution, img_resolution, length(mode_ids))
        for (idx, _) in pairs(mode_ids)
            @views stack[:, :, idx] .= reshape(H[:, idx], img_resolution, img_resolution)
        end
        return stack
    elseif case.kind === :closed_loop_trace
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        basis = build_reference_basis(case.config["basis"], tel)
        wfs = build_reference_wfs(Symbol(case.config["wfs"]["kind"]), case.config["wfs"], tel)
        compute_cfg = case.config["compute"]
        forcing_coeffs = matrix_from_rows(get(compute_cfg, "forcing_coefficients", Any[]))
        return closed_loop_trace(wfs, tel, src, basis, forcing_coeffs;
            gain=Float64(get(compute_cfg, "gain", 0.4)),
            frame_delay=Int(get(compute_cfg, "frame_delay", 2)),
            calibration_amplitude=Float64(get(compute_cfg, "calibration_amplitude", 1e-9)),
            psf_zero_padding=Int(get(compute_cfg, "psf_zero_padding", 2)))
    elseif case.kind === :gsc_closed_loop_trace
        tel = build_reference_telescope(case.config["telescope"])
        src = build_reference_source(case.config["source"])
        basis = build_reference_basis(case.config["basis"], tel)
        wfs = build_reference_wfs(Symbol(case.config["wfs"]["kind"]), case.config["wfs"], tel)
        wfs isa PyramidWFS || throw(InvalidConfiguration("GSC closed-loop trace requires a PyramidWFS"))
        compute_cfg = case.config["compute"]
        forcing_coeffs = matrix_from_rows(get(compute_cfg, "forcing_coefficients", Any[]))
        return gsc_closed_loop_trace(tel, src, wfs, basis, forcing_coeffs;
            gain=Float64(get(compute_cfg, "gain", 0.4)),
            frame_delay=Int(get(compute_cfg, "frame_delay", 2)),
            calibration_amplitude=Float64(get(compute_cfg, "calibration_amplitude", 1e-9)),
            psf_zero_padding=Int(get(compute_cfg, "psf_zero_padding", 2)),
            og_floor=Float64(get(compute_cfg, "og_floor", 0.05)))
    elseif case.kind === :gsc_atmosphere_replay_trace
        tel = build_reference_telescope(case.config["telescope"])
        ngs = build_reference_source(case.config["source"])
        sci = build_reference_source(case.config["science_source"])
        basis = build_reference_basis(case.config["basis"], tel)
        wfs = build_reference_wfs(Symbol(case.config["wfs"]["kind"]), case.config["wfs"], tel)
        wfs isa PyramidWFS || throw(InvalidConfiguration("GSC atmosphere replay trace requires a PyramidWFS"))
        compute_cfg = case.config["compute"]
        ngs_stack = load_reference_aux_array(
            case,
            compute_cfg["forcing_ngs_data"],
            compute_cfg["forcing_shape"],
            get(compute_cfg, "forcing_storage_order", "C"),
        )
        src_stack = load_reference_aux_array(
            case,
            compute_cfg["forcing_src_data"],
            compute_cfg["forcing_shape"],
            get(compute_cfg, "forcing_storage_order", "C"),
        )
        return gsc_atmosphere_replay_trace(tel, ngs, sci, wfs, basis, ngs_stack, src_stack;
            gain=Float64(get(compute_cfg, "gain", 0.2)),
            frame_delay=Int(get(compute_cfg, "frame_delay", 2)),
            calibration_amplitude=Float64(get(compute_cfg, "calibration_amplitude", 1e-9)),
            psf_zero_padding=Int(get(compute_cfg, "psf_zero_padding", 2)),
            og_floor=Float64(get(compute_cfg, "og_floor", 0.05)))
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
        style = AdaptiveOpticsSim.AcceleratorStyle(KernelAbstractions.CPU())
        AdaptiveOpticsSim._geometric_slopes!(style, slopes, tel.state.opd, wfs.state.valid_mask, sub, n_sub, offset)
        return slopes
    end
    throw(InvalidConfiguration("KA CPU reference path not implemented for reference kind '$(case.kind)'"))
end

function adapt_reference_actual(case::ReferenceCase, actual)
    if case.kind === :tomography_model_wavefront
        return adapt_wavefront_map_to_pytomoao(actual, pytomoao_model_mask(case))
    elseif case.kind === :tomography_im_wavefront
        return adapt_wavefront_map_to_pytomoao(actual, pytomoao_im_mask(case))
    elseif case.kind === :tomography_model_command_matrix
        return adapt_actuator_rows_to_pytomoao(actual, build_reference_tomography_dm(case.config["dm"]))
    elseif case.kind === :tomography_model_dm_commands
        return adapt_actuator_vector_to_pytomoao(actual, build_reference_tomography_dm(case.config["dm"]))
    end

    compare = get(case.config, "compare", nothing)
    compare === nothing && return actual
    return adapt_compare_convention(parse_reference_compare_convention(compare), actual)
end

function validate_reference_case(case::ReferenceCase)
    expected = load_reference_array(case)
    actual = adapt_reference_actual(case, compute_reference_actual(case))
    if size(actual) != size(expected)
        return (ok=false, actual=actual, expected=expected, maxabs=Inf)
    end
    expected_nan = isnan.(expected)
    actual_nan = isnan.(actual)
    if expected_nan != actual_nan
        return (ok=false, actual=actual, expected=expected, maxabs=Inf)
    end
    finite_mask = .!expected_nan
    maxabs = any(finite_mask) ? maximum(abs, actual[finite_mask] .- expected[finite_mask]) : 0.0
    ok = isapprox(actual[finite_mask], expected[finite_mask]; atol=case.atol, rtol=case.rtol)
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

    transfer_case = Dict{String,Any}(
        "kind" => "transfer_function_rejection",
        "data" => "transfer_function_rejection.txt",
        "shape" => [255, 2],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "compute" => Dict(
            "fs" => 300.0,
            "n_freq" => 256,
            "loop_gains" => [0.2, 0.6],
        ),
    )
    transfer_ref = compute_reference_actual(parse_reference_case("transfer_function_rejection", transfer_case, root))
    write_reference_array(joinpath(root, transfer_case["data"]), transfer_ref)
    manifest["cases"]["transfer_function_rejection"] = transfer_case

    lift_case = Dict{String,Any}(
        "kind" => "lift_interaction_matrix",
        "data" => "lift_interaction_matrix.txt",
        "shape" => [16, 16, 2],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "telescope" => Dict(
            "resolution" => 8,
            "diameter" => 8.0,
            "sampling_time" => 1e-3,
            "central_obstruction" => 0.0,
        ),
        "source" => Dict(
            "kind" => "ngs",
            "band" => "I",
            "magnitude" => 0.0,
        ),
        "basis" => Dict(
            "kind" => "cartesian_polynomials",
            "n_modes" => 2,
        ),
        "detector" => Dict(
            "noise" => "readout",
            "readout_sigma" => 1.0,
            "integration_time" => 1.0,
            "qe" => 1.0,
            "psf_sampling" => 2,
            "binning" => 1,
        ),
        "compute" => Dict(
            "img_resolution" => 16,
            "mode_ids" => [1, 2],
            "coefficients" => [0.0, 0.0],
            "iterations" => 3,
            "numerical" => false,
            "flux_norm" => 1.0,
        ),
    )
    lift_ref = compute_reference_actual(parse_reference_case("lift_interaction_matrix", lift_case, root))
    write_reference_array(joinpath(root, lift_case["data"]), lift_ref)
    manifest["cases"]["lift_interaction_matrix"] = lift_case

    closed_loop_case = Dict{String,Any}(
        "kind" => "closed_loop_trace",
        "data" => "closed_loop_trace.txt",
        "shape" => [4, 5],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "telescope" => Dict(
            "resolution" => 16,
            "diameter" => 8.0,
            "sampling_time" => 1e-3,
            "central_obstruction" => 0.0,
        ),
        "source" => Dict(
            "kind" => "ngs",
            "band" => "I",
            "magnitude" => 0.0,
        ),
        "basis" => Dict(
            "kind" => "cartesian_polynomials",
            "n_modes" => 2,
        ),
        "wfs" => Dict(
            "kind" => "shack_hartmann_slopes",
            "n_subap" => 4,
            "mode" => "geometric",
            "threshold" => 0.1,
        ),
        "compute" => Dict(
            "gain" => 0.4,
            "frame_delay" => 2,
            "calibration_amplitude" => 1e-9,
            "psf_zero_padding" => 2,
            "forcing_coefficients" => [
                [2e-8, -1e-8],
                [1.5e-8, 0.5e-8],
                [-1e-8, 1.25e-8],
                [0.5e-8, -0.75e-8],
            ],
        ),
    )
    closed_loop_ref = compute_reference_actual(parse_reference_case("closed_loop_trace", closed_loop_case, root))
    write_reference_array(joinpath(root, closed_loop_case["data"]), closed_loop_ref)
    manifest["cases"]["closed_loop_trace"] = closed_loop_case

    gsc_closed_loop_case = Dict{String,Any}(
        "kind" => "gsc_closed_loop_trace",
        "data" => "gsc_closed_loop_trace.txt",
        "shape" => [4, 6],
        "atol" => 1e-12,
        "rtol" => 1e-12,
        "telescope" => Dict(
            "resolution" => 16,
            "diameter" => 8.0,
            "sampling_time" => 1e-3,
            "central_obstruction" => 0.0,
        ),
        "source" => Dict(
            "kind" => "ngs",
            "band" => "R",
            "magnitude" => 8.0,
        ),
        "basis" => Dict(
            "kind" => "cartesian_polynomials",
            "n_modes" => 4,
        ),
        "wfs" => Dict(
            "kind" => "pyramid_slopes",
            "n_subap" => 4,
            "mode" => "diffractive",
            "threshold" => 0.0,
            "modulation" => 3.0,
            "modulation_points" => 8,
            "n_pix_separation" => 2,
            "n_pix_edge" => 1,
        ),
        "compute" => Dict(
            "gain" => 0.2,
            "frame_delay" => 2,
            "calibration_amplitude" => 1e-9,
            "psf_zero_padding" => 2,
            "og_floor" => 0.05,
            "forcing_coefficients" => [
                [2e-8, -1e-8, 5e-9, -2.5e-9],
                [1.5e-8, 5e-9, -7.5e-9, 5e-9],
                [-1e-8, 1.25e-8, 5e-9, -5e-9],
                [5e-9, -7.5e-9, 1e-8, 2.5e-9],
            ],
        ),
    )
    gsc_closed_loop_ref = compute_reference_actual(parse_reference_case("gsc_closed_loop_trace", gsc_closed_loop_case, root))
    write_reference_array(joinpath(root, gsc_closed_loop_case["data"]), gsc_closed_loop_ref)
    manifest["cases"]["gsc_closed_loop_trace"] = gsc_closed_loop_case

    open(reference_manifest_path(root), "w") do io
        TOML.print(io, manifest)
    end
    return root
end
