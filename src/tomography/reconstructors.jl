#
# Model-based and interaction-matrix tomography reconstructors
#
# This file implements covariance-based minimum-variance tomography used to
# recover layered or fitted wavefront estimates from guide-star slopes.
#
# Core operators:
# - `Gamma`: sparse gradient operator from pupil phase samples to x/y slopes
# - `Cxx`: guide-star phase auto-covariance
# - `Cox`: fit-phase/guide-star-phase cross-covariance
# - `Cnz`: measurement-noise covariance
# - `RecStatSA`: statistical Wiener-like reconstructor `Cox / (Cxx + Cnz)`
#
# The two main entry points differ only in how measured slopes are represented:
# - `ModelBasedTomography` builds `Gamma` explicitly and works on pupil samples
# - `InteractionMatrixTomography` assumes a measured interaction matrix already
#   maps the statistical covariances into slope space
#
abstract type AbstractTomographyMethod end
abstract type AbstractSlopeOrder end
abstract type TomographyNoiseModel end

const DEFAULT_TOMOGRAPHY_CROSS_SAMPLING = 49

struct ModelBasedTomography <: AbstractTomographyMethod end

struct InteractionMatrixTomography <: AbstractTomographyMethod end

struct SimulationSlopes <: AbstractSlopeOrder end

struct InterleavedSlopes <: AbstractSlopeOrder end

struct InvertedSlopes <: AbstractSlopeOrder end

struct RelativeSignalNoise{T<:AbstractFloat} <: TomographyNoiseModel
    fraction::T
end

struct ScalarMeasurementNoise{T<:AbstractFloat} <: TomographyNoiseModel
    variance::T
end

struct DiagonalMeasurementNoise{T<:AbstractFloat,V<:AbstractVector{T}} <: TomographyNoiseModel
    variances::V
end

struct PhotonReadoutSlopeNoise{T<:AbstractFloat} <: TomographyNoiseModel
    photons_per_subaperture::T
    readout_sigma::T
    qe::T
    n_pixels::Int
    excess_noise::T
end

RelativeSignalNoise(fraction::Real) = RelativeSignalNoise(float(fraction))
ScalarMeasurementNoise(variance::Real) = ScalarMeasurementNoise(float(variance))

function PhotonReadoutSlopeNoise(; photons_per_subaperture::Real, readout_sigma::Real=0.0,
    qe::Real=1.0, n_pixels::Integer=1, excess_noise::Real=1.0)
    T = promote_type(typeof(float(photons_per_subaperture)), typeof(float(readout_sigma)),
        typeof(float(qe)), typeof(float(excess_noise)))
    return PhotonReadoutSlopeNoise{T}(T(photons_per_subaperture), T(readout_sigma), T(qe), Int(n_pixels), T(excess_noise))
end

function PhotonReadoutSlopeNoise(det::Detector; photons_per_subaperture::Real, excess_noise::Real=1.0)
    T = eltype(det.products.frame)
    sigma = T(readout_noise(det))
    return PhotonReadoutSlopeNoise(
        photons_per_subaperture=T(photons_per_subaperture),
        readout_sigma=sigma,
        qe=det.params.qe,
        n_pixels=max(det.params.binning^2, 1),
        excess_noise=T(excess_noise),
    )
end

"""
    TomographyOperators

Cache the intermediate linear operators used to assemble a tomography
reconstructor.

- `gamma`: phase-to-slope projection for model-based tomography (`nothing`
  for a measured interaction-matrix projection)
- `cxx`: guide-star phase covariance
- `cox`: fit-phase/guide-star-phase cross-covariance
- `cnz`: slope measurement-noise covariance
- `recstat`: unscaled covariance reconstructor before physical conversion or
  DM fitting
"""
struct TomographyOperators{G,M,CX,CO,CN,RS,T}
    gamma::G
    grid_mask::M
    cxx::CX
    cox::CO
    cnz::CN
    recstat::RS
    wavefront_to_meter::T
end

"""
    TomographicReconstructor

Bundle a tomography reconstruction operator with the geometry and intermediate
operators that produced it.

`reconstructor` maps measured slopes to either masked wavefront samples or
interaction-matrix control outputs, depending on `method`.
"""
struct TomographicReconstructor{
    Method<:AbstractTomographyMethod,
    T<:AbstractFloat,
    R<:AbstractMatrix{T},
    G<:AbstractMatrix{Bool},
    AP<:TomographyAtmosphereParams,
    LP<:LGSAsterismParams,
    WP<:LGSWFSParams,
    TP<:TomographyParams,
    DP<:TomographyDMParams,
    F,
    O,
}
    method::Method
    reconstructor::R
    grid_mask::G
    atmosphere::AP
    asterism::LP
    wfs::WP
    tomography::TP
    dm::DP
    fitting::F
    operators::O
end

struct TomographyCommandReconstructor{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
    F<:TomographyFitting,
    R<:TomographicReconstructor,
    O<:AbstractSlopeOrder,
}
    matrix::M
    fitting::F
    reconstructor::R
    slope_order::O
    scaling_factor::T
end

@inline _fried_parameter(params::TomographyAtmosphereParams{T}) where {T<:AbstractFloat} =
    params.r0_zenith * cos(params.zenith_angle_rad)^(T(3) / T(5))

@inline _equal_fit_source_weights(params::TomographyParams{T}) where {T<:AbstractFloat} =
    fill(inv(T(params.n_fit_src^2)), params.n_fit_src^2)

function _active_guide_grid_params(
    rotations_rad::AbstractVector{T},
    offset_fractions_x::AbstractVector{T},
    offset_fractions_y::AbstractVector{T},
    n_gs::Integer,
) where {T<:AbstractFloat}
    n_gs == length(rotations_rad) == length(offset_fractions_x) ==
        length(offset_fractions_y) || throw(DimensionMismatchError(
        "lenslet-grid registration vectors must match the guide-star count"))
    return (
        view(rotations_rad, 1:n_gs),
        view(offset_fractions_x, 1:n_gs),
        view(offset_fractions_y, 1:n_gs),
    )
end

@kernel function fit_source_average_kernel!(out, cross, weights, n_fit::Int, n_row::Int, n_col::Int)
    i, j = @index(Global, NTuple)
    if i <= n_row && j <= n_col
        acc = zero(eltype(out))
        @inbounds for k in 1:n_fit
            acc += weights[k] * cross[k, i, j]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function diagonal_matrix_kernel!(out, variances, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds out[i, j] = i == j ? variances[i] : zero(eltype(out))
    end
end

@kernel function extract_diagonal_kernel!(out, matrix, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds out[i] = matrix[i, i]
    end
end

function _guide_star_grid(
    sampling::Integer,
    support_diameter_m::T,
    rotation_angle_rad::T,
    offset_fraction_x::T,
    offset_fraction_y::T,
) where {T<:AbstractFloat}
    coords = sampling == 1 ? range(zero(T), zero(T); length=1) :
        range(-support_diameter_m / 2, support_diameter_m / 2; length=sampling)
    s, c = sincos(rotation_angle_rad)
    xr = Matrix{T}(undef, sampling, sampling)
    yr = Matrix{T}(undef, sampling, sampling)
    @inbounds for j in 1:sampling
        y = coords[j]
        for i in 1:sampling
            x = coords[i]
            xr[j, i] = x * c - y * s - offset_fraction_x * support_diameter_m
            yr[j, i] = y * c + x * s - offset_fraction_y * support_diameter_m
        end
    end
    return xr, yr
end

function _guide_star_grid!(
    xr::AbstractMatrix{T},
    yr::AbstractMatrix{T},
    sampling::Integer,
    support_diameter_m::T,
    rotation_angle_rad::T,
    offset_fraction_x::T,
    offset_fraction_y::T,
) where {T<:AbstractFloat}
    size(xr) == (sampling, sampling) && size(yr) == (sampling, sampling) ||
        throw(DimensionMismatchError("guide-star grid workspaces must match sampling"))
    coords = sampling == 1 ? range(zero(T), zero(T); length=1) :
        range(-support_diameter_m / 2, support_diameter_m / 2; length=sampling)
    s, c = sincos(rotation_angle_rad)
    offset_x_m = offset_fraction_x * support_diameter_m
    offset_y_m = offset_fraction_y * support_diameter_m
    @inbounds for j in 1:sampling
        y = coords[j]
        for i in 1:sampling
            x = coords[i]
            xr[j, i] = x * c - y * s - offset_x_m
            yr[j, i] = y * c + x * s - offset_y_m
        end
    end
    return xr, yr
end

function _guide_star_grids!(
    xr::AbstractArray{T,3},
    yr::AbstractArray{T,3},
    sampling::Integer,
    support_diameter_m::T,
    rotations_rad::AbstractVector{T},
    offset_fractions_x::AbstractVector{T},
    offset_fractions_y::AbstractVector{T},
) where {T<:AbstractFloat}
    n_gs = length(rotations_rad)
    size(xr) == (sampling, sampling, n_gs) && size(yr) == (sampling, sampling, n_gs) ||
        throw(DimensionMismatchError("guide-star grid stack workspaces must match sampling and guide-star count"))
    length(offset_fractions_x) == n_gs == length(offset_fractions_y) ||
        throw(DimensionMismatchError("guide-star grid parameter vectors must have equal length"))
    @inbounds for gs in 1:n_gs
        _guide_star_grid!(
            @view(xr[:, :, gs]),
            @view(yr[:, :, gs]),
            sampling,
            support_diameter_m,
            rotations_rad[gs],
            offset_fractions_x[gs],
            offset_fractions_y[gs],
        )
    end
    return xr, yr
end

function _guide_star_grids(
    sampling::Integer,
    support_diameter_m::T,
    rotations_rad::AbstractVector{T},
    offset_fractions_x::AbstractVector{T},
    offset_fractions_y::AbstractVector{T},
) where {T<:AbstractFloat}
    n_gs = length(rotations_rad)
    length(offset_fractions_x) == n_gs == length(offset_fractions_y) ||
        throw(DimensionMismatchError("guide-star grid parameter vectors must have equal length"))
    xr = Array{T}(undef, sampling, sampling, n_gs)
    yr = similar(xr)
    _guide_star_grids!(xr, yr, sampling, support_diameter_m, rotations_rad,
        offset_fractions_x, offset_fractions_y)
    return xr, yr
end

struct TomographySourceLayerGeometry{T<:AbstractFloat,M<:AbstractMatrix{T}}
    beta_x::M
    beta_y::M
    scale::M
end

function _source_layer_geometry(
    direction_vectors::AbstractMatrix{T},
    layer_slant_ranges_m::AbstractVector{T},
    source_height_m::T,
) where {T<:AbstractFloat}
    n_src = size(direction_vectors, 2)
    n_layers = length(layer_slant_ranges_m)
    beta_x = Matrix{T}(undef, n_src, n_layers)
    beta_y = similar(beta_x)
    scale = similar(beta_x)
    finite_height = isfinite(source_height_m)
    @inbounds for layer in 1:n_layers
        slant_range_m = layer_slant_ranges_m[layer]
        layer_scale = finite_height ? one(T) - slant_range_m / source_height_m : one(T)
        for src in 1:n_src
            beta_x[src, layer] = direction_vectors[1, src] * slant_range_m
            beta_y[src, layer] = direction_vectors[2, src] * slant_range_m
            scale[src, layer] = layer_scale
        end
    end
    return TomographySourceLayerGeometry{T,typeof(beta_x)}(beta_x, beta_y, scale)
end

function _scaled_shifted_coords!(
    out::AbstractVector{Complex{T}},
    x::AbstractMatrix{T},
    y::AbstractMatrix{T},
    positions::AbstractVector{Int},
    beta_x::T,
    beta_y::T,
    scale::T,
) where {T<:AbstractFloat}
    length(out) == length(positions) ||
        throw(DimensionMismatchError("coordinate workspace length must match selected positions"))
    @inbounds for k in eachindex(out, positions)
        idx = positions[k]
        out[k] = complex(x[idx] * scale + beta_x, y[idx] * scale + beta_y)
    end
    return out
end

function _selected_transformed_coordinate_stack(
    x::AbstractArray{T,3},
    y::AbstractArray{T,3},
    positions::AbstractVector{Int},
    geometry::TomographySourceLayerGeometry{T},
) where {T<:AbstractFloat}
    size(x) == size(y) || throw(DimensionMismatchError(
        "coordinate grids must have matching axes",
    ))
    sample_count = length(positions)
    source_count = size(x, 3)
    layer_count = size(geometry.scale, 2)
    size(geometry.scale, 1) == source_count || throw(DimensionMismatchError(
        "source-layer geometry must match coordinate-grid source count",
    ))
    coordinates = Array{Complex{T}}(undef, sample_count, source_count, layer_count)
    @inbounds for layer in 1:layer_count, source in 1:source_count
        _scaled_shifted_coords!(
            @view(coordinates[:, source, layer]),
            @view(x[:, :, source]),
            @view(y[:, :, source]),
            positions,
            geometry.beta_x[source, layer],
            geometry.beta_y[source, layer],
            geometry.scale[source, layer],
        )
    end
    return coordinates
end

function _selected_transformed_coordinate_stack(
    x::AbstractMatrix{T},
    y::AbstractMatrix{T},
    positions::AbstractVector{Int},
    geometry::TomographySourceLayerGeometry{T},
) where {T<:AbstractFloat}
    size(x) == size(y) || throw(DimensionMismatchError(
        "coordinate grids must have matching axes",
    ))
    sample_count = length(positions)
    source_count, layer_count = size(geometry.scale)
    coordinates = Array{Complex{T}}(undef, sample_count, source_count, layer_count)
    @inbounds for layer in 1:layer_count, source in 1:source_count
        _scaled_shifted_coords!(
            @view(coordinates[:, source, layer]),
            x,
            y,
            positions,
            geometry.beta_x[source, layer],
            geometry.beta_y[source, layer],
            geometry.scale[source, layer],
        )
    end
    return coordinates
end

function _require_aoc_von_karman_precision(::Type{T}) where {T<:AbstractFloat}
    T <: Union{Float32,Float64} || throw(UnsupportedAlgorithm(
        "tomography von Kármán covariance requires Float32 or Float64 because AdaptiveOpticsCalibration supplies the numerical assembly",
    ))
    return nothing
end

function _aoc_auto_covariance(
    backend::GPUArrayBuildBackend{B},
    guide_coordinates::AbstractArray{Complex{T},3},
    atmosphere::TomographyAtmosphereParams{T},
) where {B,T<:AbstractFloat}
    _require_aoc_von_karman_precision(T)
    tomography = AdaptiveOpticsCalibration.Tomography
    specification = tomography.VonKarmanAutoCovarianceSpecification(
        guide_coordinates,
        _fried_parameter(atmosphere),
        atmosphere.L0,
        atmosphere.fractional_cn2,
    )
    ka_backend = KernelAbstractions.get_backend(_backend_array(B, T, 0))
    plan = AdaptiveOpticsCalibration.prepare(
        tomography.VonKarmanAutoCovariance(),
        AdaptiveOpticsCalibration.KernelExecution(specification, ka_backend),
    )
    return tomography.auto_covariance(AdaptiveOpticsCalibration.process(plan, nothing))
end

function _aoc_auto_covariance(
    guide_coordinates::AbstractArray{Complex{T},3},
    atmosphere::TomographyAtmosphereParams{T},
) where {T<:AbstractFloat}
    _require_aoc_von_karman_precision(T)
    tomography = AdaptiveOpticsCalibration.Tomography
    specification = tomography.VonKarmanAutoCovarianceSpecification(
        guide_coordinates,
        _fried_parameter(atmosphere),
        atmosphere.L0,
        atmosphere.fractional_cn2,
    )
    plan = AdaptiveOpticsCalibration.prepare(
        tomography.VonKarmanAutoCovariance(), specification)
    return tomography.auto_covariance(AdaptiveOpticsCalibration.process(plan, nothing))
end

function _aoc_cross_covariance(
    guide_coordinates::AbstractArray{Complex{T},3},
    fit_coordinates::AbstractArray{Complex{T},3},
    atmosphere::TomographyAtmosphereParams{T},
) where {T<:AbstractFloat}
    _require_aoc_von_karman_precision(T)
    tomography = AdaptiveOpticsCalibration.Tomography
    specification = tomography.VonKarmanCrossCovarianceSpecification(
        guide_coordinates,
        fit_coordinates,
        _fried_parameter(atmosphere),
        atmosphere.L0,
        atmosphere.fractional_cn2,
    )
    plan = AdaptiveOpticsCalibration.prepare(
        tomography.VonKarmanCrossCovariance(), specification)
    return tomography.cross_covariance(AdaptiveOpticsCalibration.process(plan, nothing))
end

function _aoc_cross_covariance(
    backend::GPUArrayBuildBackend{B},
    guide_coordinates::AbstractArray{Complex{T},3},
    fit_coordinates::AbstractArray{Complex{T},3},
    atmosphere::TomographyAtmosphereParams{T},
) where {B,T<:AbstractFloat}
    _require_aoc_von_karman_precision(T)
    tomography = AdaptiveOpticsCalibration.Tomography
    specification = tomography.VonKarmanCrossCovarianceSpecification(
        guide_coordinates,
        fit_coordinates,
        _fried_parameter(atmosphere),
        atmosphere.L0,
        atmosphere.fractional_cn2,
    )
    ka_backend = KernelAbstractions.get_backend(_backend_array(B, T, 0))
    plan = AdaptiveOpticsCalibration.prepare(
        tomography.VonKarmanCrossCovariance(),
        AdaptiveOpticsCalibration.KernelExecution(specification, ka_backend),
    )
    return tomography.cross_covariance(AdaptiveOpticsCalibration.process(plan, nothing))
end

function sparse_gradient_matrix(
    valid_lenslet::AbstractMatrix{Bool};
    amplitude_mask::Union{Nothing,AbstractMatrix}=nothing,
    over_sampling::Integer=2,
)
    over_sampling == 2 || throw(UnsupportedAlgorithm("only over_sampling=2 is implemented"))
    n_lenslets = size(valid_lenslet, 1)
    size(valid_lenslet, 2) == n_lenslets ||
        throw(DimensionMismatchError("valid_lenslet must be square"))
    n_map = over_sampling * n_lenslets + 1
    amp = amplitude_mask === nothing ? ones(Bool, n_map, n_map) : convert.(Bool, amplitude_mask)
    size(amp) == (n_map, n_map) ||
        throw(DimensionMismatchError("amplitude_mask size must match oversampled grid"))

    stencil_x = (-0.25, -0.5, -0.25, 0.0, 0.0, 0.0, 0.25, 0.5, 0.25)
    stencil_y = (-0.25, 0.0, 0.25, -0.5, 0.0, 0.5, -0.25, 0.0, 0.25)
    offsets_i = (0, 1, 2, 0, 1, 2, 0, 1, 2)
    offsets_j = (0, 0, 0, 1, 1, 1, 2, 2, 2)

    grid_mask = falses(n_map, n_map)
    valid_cells = Tuple{Int, Int, Int}[]
    row_id = 0
    @inbounds for j_lenslet in 1:n_lenslets
        j_offset = over_sampling * (j_lenslet - 1)
        for i_lenslet in 1:n_lenslets
            valid_lenslet[i_lenslet, j_lenslet] || continue
            i_offset = over_sampling * (i_lenslet - 1)
            patch = @view amp[i_offset+1:i_offset+over_sampling+1, j_offset+1:j_offset+over_sampling+1]
            all(patch) || continue
            row_id += 1
            push!(valid_cells, (row_id, i_offset, j_offset))
            for s in eachindex(offsets_i)
                grid_mask[i_offset + offsets_i[s] + 1, j_offset + offsets_j[s] + 1] = true
            end
        end
    end

    n_valid = length(valid_cells)
    col_map = zeros(Int, n_map, n_map)
    mask_positions = findall(grid_mask)
    for (k, idx) in enumerate(mask_positions)
        col_map[idx] = k
    end

    rows = Int[]
    cols = Int[]
    vals = Float64[]
    sizehint!(rows, 18 * n_valid)
    sizehint!(cols, 18 * n_valid)
    sizehint!(vals, 18 * n_valid)
    @inbounds for (local_row, i_offset, j_offset) in valid_cells
        xrow = local_row
        yrow = n_valid + local_row
        for s in eachindex(offsets_i)
            ii = i_offset + offsets_i[s] + 1
            jj = j_offset + offsets_j[s] + 1
            col = col_map[ii, jj]
            col == 0 && continue
            push!(rows, xrow)
            push!(cols, col)
            push!(vals, stencil_x[s])
            push!(rows, yrow)
            push!(cols, col)
            push!(vals, stencil_y[s])
        end
    end

    gamma = sparse(rows, cols, vals, 2 * n_valid, length(mask_positions))
    return gamma, grid_mask
end

"""
    auto_correlation(..., grid_mask)

Assemble the guide-star phase auto-covariance `Cxx` over the masked pupil grid.

Each block integrates the von Karman covariance across atmospheric layers after
shifting each guide-star pupil footprint by the layer geometry. AOS supplies
the physical coordinates; AdaptiveOpticsCalibration assembles the numerical
covariance on CPU or the selected accelerator for Float32 and Float64.
"""
function auto_correlation(
    backend::BuildBackend,
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    return auto_correlation(atmosphere, asterism, wfs, grid_mask)
end

function _auto_covariance_coordinates(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    sampling::Int,
    valid_positions::AbstractVector{Int},
) where {T<:AbstractFloat}
    n_gs = asterism.n_lgs
    slant_ranges_m = layer_slant_ranges_m(atmosphere)
    support_diameter_m = lenslet_grid_support_diameter_m(wfs)
    lgs_dir = lgs_directions(asterism)
    directions = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    source_height_m = lgs_height_m(asterism, atmosphere)
    geometry = _source_layer_geometry(directions, slant_ranges_m, source_height_m)
    rotations_rad, offset_fractions_x, offset_fractions_y = _active_guide_grid_params(
        wfs.lenslet_grid_rotations_rad,
        view(wfs.lenslet_grid_offsets_fraction, 1, :),
        view(wfs.lenslet_grid_offsets_fraction, 2, :),
        n_gs,
    )

    guide_x, guide_y = _guide_star_grids(
        sampling,
        support_diameter_m,
        rotations_rad,
        offset_fractions_x,
        offset_fractions_y,
    )
    return _selected_transformed_coordinate_stack(
        guide_x,
        guide_y,
        valid_positions,
        geometry,
    )
end

function auto_correlation(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    sampling = size(grid_mask, 1)
    size(grid_mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    valid_positions = findall(vec(grid_mask))
    n_valid = length(valid_positions)
    n_gs = asterism.n_lgs
    (iszero(n_valid) || iszero(n_gs)) && return zeros(T, n_gs * n_valid, n_gs * n_valid)
    guide_coordinates = _auto_covariance_coordinates(
        atmosphere, asterism, wfs, sampling, valid_positions)
    return _aoc_auto_covariance(guide_coordinates, atmosphere)
end

function auto_correlation(
    backend::GPUArrayBuildBackend{B},
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {B,T<:AbstractFloat}
    sampling = size(grid_mask, 1)
    size(grid_mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    valid_positions = findall(vec(grid_mask))
    n_valid = length(valid_positions)
    n_gs = asterism.n_lgs
    (iszero(n_valid) || iszero(n_gs)) &&
        return _backend_array(B, T, n_gs * n_valid, n_gs * n_valid)
    guide_coordinates = _auto_covariance_coordinates(
        atmosphere, asterism, wfs, sampling, valid_positions)
    return _aoc_auto_covariance(backend, guide_coordinates, atmosphere)
end

"""
    cross_correlation(...; grid_mask=nothing)

Assemble the phase cross-covariance `Cox` between fit directions and guide-star
pupil samples.

The result is stacked over fit sources. The model builder supplies its phase
grid mask and averages the selected fit sources for the statistical
reconstructor. AOS supplies transformed physical coordinates;
AdaptiveOpticsCalibration assembles the numerical covariance on CPU or the
selected accelerator for Float32 and Float64.
"""
function cross_correlation(
    backend::BuildBackend,
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {T<:AbstractFloat}
    return cross_correlation(atmosphere, asterism, wfs, tomography; grid_mask=grid_mask)
end

function _cross_covariance_coordinates(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T},
    sampling::Int,
    row_positions::AbstractVector{Int},
) where {T<:AbstractFloat}
    n_gs = asterism.n_lgs
    slant_ranges_m = layer_slant_ranges_m(atmosphere)
    support_diameter_m = lenslet_grid_support_diameter_m(wfs)
    lgs_dir = lgs_directions(asterism)
    lgs_directions_xyz = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    fit_zenith, fit_azimuth = optimization_geometry(tomography)
    fit_directions_xyz = direction_vectors(fit_zenith, fit_azimuth)
    source_height_m = lgs_height_m(asterism, atmosphere)
    lgs_geometry = _source_layer_geometry(lgs_directions_xyz, slant_ranges_m,
        source_height_m)
    fit_geometry = _source_layer_geometry(fit_directions_xyz, slant_ranges_m,
        tomography.fit_src_height_m)
    target_x, target_y = _guide_star_grid(sampling, support_diameter_m, zero(T),
        zero(T), zero(T))
    rotations_rad, offset_fractions_x, offset_fractions_y = _active_guide_grid_params(
        wfs.lenslet_grid_rotations_rad,
        view(wfs.lenslet_grid_offsets_fraction, 1, :),
        view(wfs.lenslet_grid_offsets_fraction, 2, :),
        n_gs,
    )
    guide_x, guide_y = _guide_star_grids(
        sampling,
        support_diameter_m,
        rotations_rad,
        offset_fractions_x,
        offset_fractions_y,
    )
    guide_coordinates = _selected_transformed_coordinate_stack(
        guide_x,
        guide_y,
        row_positions,
        lgs_geometry,
    )
    fit_coordinates = _selected_transformed_coordinate_stack(
        target_x,
        target_y,
        row_positions,
        fit_geometry,
    )
    return guide_coordinates, fit_coordinates
end

function cross_correlation(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {T<:AbstractFloat}
    sampling = isnothing(grid_mask) ? DEFAULT_TOMOGRAPHY_CROSS_SAMPLING : size(grid_mask, 1)
    mask = isnothing(grid_mask) ? trues(sampling, sampling) : grid_mask
    size(mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    row_positions = findall(vec(mask))
    n_row = length(row_positions)
    n_fit = tomography.n_fit_src^2
    n_gs = asterism.n_lgs
    (iszero(n_row) || iszero(n_gs)) && return Array{T}(undef, n_fit, n_row, n_gs * n_row)
    guide_coordinates, fit_coordinates = _cross_covariance_coordinates(
        atmosphere, asterism, wfs, tomography, sampling, row_positions)
    return _aoc_cross_covariance(guide_coordinates, fit_coordinates, atmosphere)
end

function cross_correlation(
    backend::GPUArrayBuildBackend{B},
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {B,T<:AbstractFloat}
    sampling = isnothing(grid_mask) ? DEFAULT_TOMOGRAPHY_CROSS_SAMPLING : size(grid_mask, 1)
    mask = isnothing(grid_mask) ? trues(sampling, sampling) : grid_mask
    size(mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    row_positions = findall(vec(mask))
    n_row = length(row_positions)
    n_fit = tomography.n_fit_src^2
    n_gs = asterism.n_lgs
    (iszero(n_row) || iszero(n_gs)) &&
        return _backend_array(B, T, n_fit, n_row, n_gs * n_row)
    guide_coordinates, fit_coordinates = _cross_covariance_coordinates(
        atmosphere, asterism, wfs, tomography, sampling, row_positions)
    return _aoc_cross_covariance(
        backend, guide_coordinates, fit_coordinates, atmosphere)
end

function _fit_source_average(cross::Array{T,3}, weights::AbstractVector{T}) where {T<:AbstractFloat}
    size(cross, 1) == length(weights) ||
        throw(DimensionMismatchError("fit-source weight length must match cross-correlation stack"))
    out = zeros(T, size(cross, 2), size(cross, 3))
    for i in axes(cross, 1)
        @views out .+= weights[i] .* cross[i, :, :]
    end
    return out
end

function _fit_source_average(cross::AbstractArray{T,3}, weights::AbstractVector{T}) where {T<:AbstractFloat}
    size(cross, 1) == length(weights) ||
        throw(DimensionMismatchError("fit-source weight length must match cross-correlation stack"))
    out = similar(cross, T, size(cross, 2), size(cross, 3))
    weights_native = similar(cross, T, size(cross, 1))
    copyto!(weights_native, weights)
    style = execution_style(out)
    launch_kernel_async!(style, fit_source_average_kernel!, out, cross, weights_native, size(cross, 1), size(cross, 2), size(cross, 3);
        ndrange=size(out))
    return out
end

_covariance_input_type(::Type{T}, matrix::AbstractMatrix{T}) where {T<:AbstractFloat} = matrix
_covariance_input_type(::Type{T}, matrix::AbstractMatrix) where {T<:AbstractFloat} =
    Matrix{T}(matrix)

function _tomographic_covariance_reconstructor(
    ::ScalarCPUStyle,
    ::BuildBackend,
    projection::AbstractMatrix,
    phase_covariance::AbstractMatrix{T},
    fit_phase_covariance::AbstractMatrix,
    measurement_noise_covariance::AbstractMatrix,
) where {T<:AbstractFloat}
    aoc_tomography = AdaptiveOpticsCalibration.Tomography
    S = promote_type(T, eltype(projection), eltype(fit_phase_covariance),
        eltype(measurement_noise_covariance))
    specification = aoc_tomography.CovarianceReconstructorSpecification(
        size(projection, 2), size(projection, 1), size(fit_phase_covariance, 1), S)
    plan = AdaptiveOpticsCalibration.prepare(
        aoc_tomography.CovarianceReconstructor(), specification)
    inputs = aoc_tomography.CovarianceReconstructorInputs(
        _covariance_input_type(S, projection),
        _covariance_input_type(S, phase_covariance),
        _covariance_input_type(S, fit_phase_covariance),
        _covariance_input_type(S, measurement_noise_covariance))
    return aoc_tomography.reconstructor(
        AdaptiveOpticsCalibration.process(plan, inputs))
end

function _tomographic_covariance_reconstructor(
    ::AcceleratorStyle,
    ::BuildBackend,
    projection::AbstractMatrix,
    phase_covariance::AbstractMatrix{T},
    fit_phase_covariance::AbstractMatrix,
    measurement_noise_covariance::AbstractMatrix,
) where {T<:AbstractFloat}
    aoc_tomography = AdaptiveOpticsCalibration.Tomography
    specification = aoc_tomography.CovarianceReconstructorSpecification(
        size(projection, 2), size(projection, 1), size(fit_phase_covariance, 1), T)
    plan = AdaptiveOpticsCalibration.prepare(
        aoc_tomography.CovarianceReconstructor(),
        AdaptiveOpticsCalibration.KernelExecution(
            specification, KernelAbstractions.get_backend(phase_covariance)))
    inputs = aoc_tomography.CovarianceReconstructorInputs(
        projection, phase_covariance, fit_phase_covariance,
        measurement_noise_covariance)
    return aoc_tomography.reconstructor(
        AdaptiveOpticsCalibration.process(plan, inputs))
end

tomography_noise_covariance(model::TomographyNoiseModel, reference_diag::AbstractVector) =
    tomography_noise_covariance(NativeBuildBackend(), model, reference_diag)

tomography_reference_diagonal(::BuildBackend, matrix::AbstractMatrix) = diag(matrix)

function tomography_reference_diagonal(
    ::GPUArrayBuildBackend{B},
    matrix::AbstractMatrix{T},
) where {B,T<:AbstractFloat}
    n = min(size(matrix)...)
    diagonal = _backend_array(B, T, n)
    launch_kernel_async!(execution_style(matrix), extract_diagonal_kernel!, diagonal, matrix, n;
        ndrange=n)
    return diagonal
end

function _build_diagonal_noise(::NativeBuildBackend, variances::AbstractVector{T}) where {T<:AbstractFloat}
    return Diagonal(variances)
end

function _build_diagonal_noise(::CPUBuildBackend, variances::AbstractVector{T}) where {T<:AbstractFloat}
    return Diagonal(Vector(variances))
end

function _build_diagonal_noise(backend::GPUArrayBuildBackend{B}, variances::AbstractVector{T}) where {B,T<:AbstractFloat}
    out = _backend_array(B, T, length(variances), length(variances))
    style = execution_style(out)
    launch_kernel_async!(style, diagonal_matrix_kernel!, out, variances, length(variances); ndrange=size(out))
    return out
end

function tomography_noise_covariance(::NativeBuildBackend, model::RelativeSignalNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    variances = similar(reference_diag, T)
    fraction = T(model.fraction)
    @. variances = max(fraction * max(reference_diag, eps(T)), eps(T))
    return _build_diagonal_noise(NativeBuildBackend(), variances)
end

function tomography_noise_covariance(backend::BuildBackend, model::RelativeSignalNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    variances = similar(reference_diag, T)
    fraction = T(model.fraction)
    @. variances = max(fraction * max(reference_diag, eps(T)), eps(T))
    return _build_diagonal_noise(backend, variances)
end

function tomography_noise_covariance(::NativeBuildBackend, model::ScalarMeasurementNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    n = length(reference_diag)
    return Diagonal(fill(max(T(model.variance), eps(T)), n))
end

function tomography_noise_covariance(backend::BuildBackend, model::ScalarMeasurementNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    variances = similar(reference_diag, T)
    fill!(variances, max(T(model.variance), eps(T)))
    return _build_diagonal_noise(backend, variances)
end

function tomography_noise_covariance(::NativeBuildBackend, model::DiagonalMeasurementNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    length(model.variances) == length(reference_diag) ||
        throw(DimensionMismatchError("tomography noise variances must match slope dimension"))
    variances = similar(reference_diag, T)
    @inbounds for i in eachindex(variances, model.variances)
        variances[i] = max(T(model.variances[i]), eps(T))
    end
    return Diagonal(variances)
end

function tomography_noise_covariance(backend::BuildBackend, model::DiagonalMeasurementNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    length(model.variances) == length(reference_diag) ||
        throw(DimensionMismatchError("tomography noise variances must match slope dimension"))
    variances = similar(reference_diag, T)
    @inbounds for i in eachindex(variances, model.variances)
        variances[i] = max(T(model.variances[i]), eps(T))
    end
    return _build_diagonal_noise(backend, variances)
end

function tomography_noise_covariance(::NativeBuildBackend, model::PhotonReadoutSlopeNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    n = length(reference_diag)
    photoelectrons = max(T(model.photons_per_subaperture) * T(model.qe), eps(T))
    shot = T(model.excess_noise)^2 / photoelectrons
    readout = T(model.n_pixels) * T(model.readout_sigma)^2 / (photoelectrons^2)
    σ2 = max(shot + readout, eps(T))
    return Diagonal(fill(σ2, n))
end

function tomography_noise_covariance(backend::BuildBackend, model::PhotonReadoutSlopeNoise{S},
    reference_diag::AbstractVector{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    photoelectrons = max(T(model.photons_per_subaperture) * T(model.qe), eps(T))
    shot = T(model.excess_noise)^2 / photoelectrons
    readout = T(model.n_pixels) * T(model.readout_sigma)^2 / (photoelectrons^2)
    σ2 = max(shot + readout, eps(T))
    variances = similar(reference_diag, T)
    fill!(variances, σ2)
    return _build_diagonal_noise(backend, variances)
end

"""
    build_reconstructor(InteractionMatrixTomography(), interaction_matrix, ...)

Build a statistical tomography reconstructor in measured slope space.

This path treats `interaction_matrix` as the slope-space projection operator,
forms `Cxx`, `Cox`, and `Cnz`, and computes the Wiener-like solve
`Cox * A' / (A * Cxx * A' + Cnz)`.
"""
function build_reconstructor(
    ::InteractionMatrixTomography,
    interaction_matrix::AbstractMatrix{T},
    grid_mask::AbstractMatrix{Bool},
    atmosphere::TomographyAtmosphereParams,
    asterism::LGSAsterismParams,
    wfs::LGSWFSParams,
    tomography::TomographyParams,
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
    α::Real=10,
    noise_model::TomographyNoiseModel=RelativeSignalNoise(1e-3 * α),
    build_backend::BuildBackend=default_build_backend(interaction_matrix),
) where {T<:AbstractFloat}
    size(interaction_matrix, 1) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one row"))
    size(interaction_matrix, 2) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one column"))

    interaction_native = materialize_build(build_backend, interaction_matrix, interaction_matrix)
    cxx = auto_correlation(build_backend, atmosphere, asterism, wfs, grid_mask)
    cross = cross_correlation(build_backend, atmosphere, asterism, wfs, tomography; grid_mask=grid_mask)
    cox = _fit_source_average(cross, _equal_fit_source_weights(tomography))
    cxx_native = materialize_build(build_backend, interaction_native, cxx)
    cox_native = materialize_build(build_backend, interaction_native, cox)
    css_signal = backend_symmetric_product(interaction_native, cxx_native)
    reference_diag = tomography_reference_diagonal(build_backend, css_signal)
    cnz = tomography_noise_covariance(build_backend, noise_model, reference_diag)
    recstat = _tomographic_covariance_reconstructor(execution_style(cxx_native),
        build_backend, interaction_native, cxx_native, cox_native, cnz)
    native_mask = materialize_build(build_backend, interaction_native, grid_mask)
    operators = TomographyOperators(
        nothing,
        native_mask,
        cxx_native,
        cox_native,
        cnz,
        recstat,
        one(T),
    )
    return TomographicReconstructor(
        InteractionMatrixTomography(),
        recstat,
        native_mask,
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm,
        fitting,
        operators,
    )
end

function swap_xy_blocks(
    matrix::AbstractMatrix{T},
    n_valid_subap::Integer;
    n_channels::Integer=1,
) where {T}
    cols_per_channel = 2 * Int(n_valid_subap)
    size(matrix, 2) == cols_per_channel * Int(n_channels) ||
        throw(DimensionMismatchError("matrix column count must equal 2*n_valid_subap*n_channels"))
    reordered = similar(matrix, T, size(matrix))
    @inbounds for channel in 1:Int(n_channels)
        channel_start = (channel - 1) * cols_per_channel
        src_x = channel_start + Int(n_valid_subap) + 1:channel_start + cols_per_channel
        src_y = channel_start + 1:channel_start + Int(n_valid_subap)
        dst = channel_start + 1:channel_start + cols_per_channel
        reordered[:, dst] .= matrix[:, vcat(src_x, src_y)]
    end
    return reordered
end

function interleave_xy_columns(
    matrix::AbstractMatrix{T},
    n_valid_subap::Integer;
    n_channels::Integer=1,
) where {T}
    cols_per_channel = 2 * Int(n_valid_subap)
    size(matrix, 2) == cols_per_channel * Int(n_channels) ||
        throw(DimensionMismatchError("matrix column count must equal 2*n_valid_subap*n_channels"))
    reordered = similar(matrix, T, size(matrix))
    @inbounds for channel in 1:Int(n_channels)
        channel_start = (channel - 1) * cols_per_channel
        for subap in 1:Int(n_valid_subap)
            reordered[:, channel_start + 2subap - 1] .= matrix[:, channel_start + subap]
            reordered[:, channel_start + 2subap] .= matrix[:, channel_start + Int(n_valid_subap) + subap]
        end
    end
    return reordered
end

prepare_slope_order(
    ::InvertedSlopes,
    matrix::AbstractMatrix,
    n_valid_subap::Integer;
    n_channels::Integer=1,
) = copy(matrix)

function prepare_slope_order(
    ::SimulationSlopes,
    matrix::AbstractMatrix,
    n_valid_subap::Integer;
    n_channels::Integer=1,
)
    return swap_xy_blocks(matrix, n_valid_subap; n_channels=n_channels)
end

function prepare_slope_order(
    ::InterleavedSlopes,
    matrix::AbstractMatrix,
    n_valid_subap::Integer;
    n_channels::Integer=1,
)
    swapped = swap_xy_blocks(matrix, n_valid_subap; n_channels=n_channels)
    return interleave_xy_columns(swapped, n_valid_subap; n_channels=n_channels)
end

"""
    assemble_reconstructor_and_fitting(reconstructor, dm; ...)

Project a model-based tomography reconstructor onto DM actuator commands.

This applies the requested slope-order convention, builds or reuses the DM
fitting operator, and composes the final slope-to-command matrix.
"""
function assemble_reconstructor_and_fitting(
    reconstructor::TomographicReconstructor{ModelBasedTomography,T},
    dm::TomographyDMParams{T};
    n_channels::Integer=reconstructor.asterism.n_lgs,
    slope_order::AbstractSlopeOrder=SimulationSlopes(),
    scaling_factor::Real=1.65e7,
    fitting::Union{Nothing,TomographyFitting}=nothing,
    regularization::Real=sqrt(eps(T)),
    w1::Real=2,
    w2::Real=-1,
    sigma1::Real=1.0,
    sigma2::Real=1.7,
    stretch_factor::Real=1.03,
    build_backend::BuildBackend=default_build_backend(reconstructor.reconstructor),
) where {T<:AbstractFloat}
    n_valid = n_valid_subapertures(reconstructor.wfs)
    n_channels >= 1 || throw(InvalidConfiguration("n_channels must be positive"))
    base = n_channels == 1 ?
        reconstructor.reconstructor[:, 1:2*n_valid] :
        reconstructor.reconstructor
    ordered = prepare_slope_order(slope_order, base, n_valid; n_channels=n_channels)
    fit = isnothing(fitting) ?
        TomographyFitting(
            dm;
            regularization=regularization,
            resolution=size(reconstructor.grid_mask, 1),
            w1=w1,
            w2=w2,
            sigma1=sigma1,
            sigma2=sigma2,
            stretch_factor=stretch_factor,
        ) :
        fitting
    modes = fit.influence_functions[vec(reconstructor.grid_mask), :]
    masked_fitting = TomographyFitting(modes; regularization=regularization, resolution=size(reconstructor.grid_mask, 1))
    fitting_matrix = materialize_build(build_backend, reconstructor.reconstructor, masked_fitting.fitting_matrix)
    matrix = -(fitting_matrix * ordered) * T(scaling_factor)
    return TomographyCommandReconstructor(
        matrix,
        masked_fitting,
        reconstructor,
        slope_order,
        T(scaling_factor),
    )
end

function mask_actuators!(
    reconstructor::TomographyCommandReconstructor,
    actuator_indices,
)
    reconstructor.matrix[actuator_indices, :] .= zero(eltype(reconstructor.matrix))
    return reconstructor
end

function build_reconstructor(
    ::InteractionMatrixTomography,
    imat::InteractionMatrix,
    grid_mask::AbstractMatrix{Bool},
    atmosphere::TomographyAtmosphereParams,
    asterism::LGSAsterismParams,
    wfs::LGSWFSParams,
    tomography::TomographyParams,
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
    α::Real=10,
    noise_model::TomographyNoiseModel=RelativeSignalNoise(1e-3 * α),
    build_backend::BuildBackend=default_build_backend(imat.matrix),
)
    return build_reconstructor(
        InteractionMatrixTomography(),
        imat.matrix,
        grid_mask,
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm;
        fitting=fitting,
        α=α,
        noise_model=noise_model,
        build_backend=build_backend,
    )
end

"""
    build_reconstructor(ModelBasedTomography(), atmosphere, asterism, wfs, tomography, dm; ...)

Build the full covariance-model tomography reconstructor.

This path constructs the phase-to-slope gradient operator `P`, forms the masked
covariance matrices `Cxx`, `Cox`, and `Cnz`, then asks
AdaptiveOpticsCalibration to evaluate the covariance reconstructor
`R = Cox * P' / (P * Cxx * P' + Cnz)` on the selected CPU or GPU backend.
"""
function build_reconstructor(
    ::ModelBasedTomography,
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T},
    dm::TomographyDMParams;
    fitting::Union{Nothing,TomographyFitting}=nothing,
    noise_model::TomographyNoiseModel=RelativeSignalNoise(0.1),
    build_backend::BuildBackend=NativeBuildBackend(),
) where {T<:AbstractFloat}
    gamma_single, grid_mask = sparse_gradient_matrix(valid_lenslet_support(wfs); over_sampling=2)
    gamma_t = SparseMatrixCSC{T, Int}(gamma_single)
    gamma = blockdiag(ntuple(_ -> gamma_t, asterism.n_lgs)...)
    cxx = auto_correlation(build_backend, atmosphere, asterism, wfs, grid_mask)
    cross = cross_correlation(build_backend, atmosphere, asterism, wfs, tomography;
        grid_mask=grid_mask)
    cox = _fit_source_average(cross, _equal_fit_source_weights(tomography))
    gamma_native = materialize_build(build_backend, gamma, gamma)
    cxx_native = materialize_build(build_backend, gamma_native, cxx)
    cox_native = materialize_build(build_backend, gamma_native, cox)
    native_mask = materialize_build(build_backend, gamma_native, grid_mask)
    css_signal = backend_symmetric_product(gamma_native, cxx_native)
    reference_diag = tomography_reference_diagonal(build_backend, css_signal)
    cnz = tomography_noise_covariance(build_backend, noise_model, reference_diag)
    recstat = _tomographic_covariance_reconstructor(execution_style(cxx_native),
        build_backend, gamma_native, cxx_native, cox_native, cnz)
    d = lenslet_grid_support_diameter_m(wfs) / size(valid_lenslet_support(wfs), 1)
    wavefront_to_meter = asterism.wavelength_m / d / 2
    recon = d * wavefront_to_meter .* recstat
    operators = TomographyOperators(
        gamma_native,
        native_mask,
        cxx_native,
        cox_native,
        cnz,
        recstat,
        wavefront_to_meter,
    )
    return TomographicReconstructor(
        ModelBasedTomography(),
        recon,
        native_mask,
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm,
        fitting,
        operators,
    )
end
