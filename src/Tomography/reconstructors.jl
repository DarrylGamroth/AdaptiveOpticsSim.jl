using LinearAlgebra
using SparseArrays

abstract type AbstractTomographyMethod end
abstract type AbstractSlopeOrder end
abstract type TomographyNoiseModel end

const PYTOMOAO_DEFAULT_CROSS_SAMPLING = 49

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
    T = eltype(det.state.frame)
    sigma = det.noise isa NoiseReadout ? T(det.noise.sigma) :
        det.noise isa NoisePhotonReadout ? T(det.noise.sigma) : zero(T)
    return PhotonReadoutSlopeNoise(
        photons_per_subaperture=T(photons_per_subaperture),
        readout_sigma=sigma,
        qe=det.params.qe,
        n_pixels=max(det.params.binning^2, 1),
        excess_noise=T(excess_noise),
    )
end

struct TomographyOperators{G,M,CX,CO,CN,RS,T}
    gamma::G
    grid_mask::M
    cxx::CX
    cox::CO
    cnz::CN
    recstat::RS
    wavefront_to_meter::T
end

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

mutable struct TomographyCommandReconstructor{
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
    params.r0_zenith * cosd(params.zenith_angle_deg)^(T(3) / T(5))

@inline _equal_fit_source_weights(params::TomographyParams{T}) where {T<:AbstractFloat} =
    fill(inv(T(params.n_fit_src^2)), params.n_fit_src^2)

@kernel function covariance_matrix_kernel!(out, rho1, rho2, cst, var_term, inv_L0, fractional_r0, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n1 && j <= n2
        rho = abs(@inbounds(rho1[i] - rho2[j]))
        if iszero(rho)
            @inbounds out[i, j] = var_term * fractional_r0
        else
            u = (2 * π) * rho * inv_L0
            @inbounds out[i, j] = cst * u^(5 / 6) * real(_kv56_scalar(u)) * fractional_r0
        end
    end
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

function _kv56_scalar(z::Complex{T}) where {T<:AbstractFloat}
    gamma_1_6 = T(5.56631600178)
    gamma_11_6 = T(0.94065585824)
    v = T(5) / T(6)
    z_abs = abs(z)
    if z_abs < T(2)
        sum_a = zero(Complex{T})
        sum_b = zero(Complex{T})
        term_a = (T(0.5) * z)^v / gamma_11_6
        term_b = (T(0.5) * z)^(-v) / gamma_1_6
        sum_a += term_a
        sum_b += term_b
        z_sq_over_4 = (T(0.5) * z)^2
        k = 1
        tol = T(1e-15)
        while k <= 1000
            factor_a = z_sq_over_4 / (T(k) * (T(k) + v))
            term_a *= factor_a
            sum_a += term_a
            factor_b = z_sq_over_4 / (T(k) * (T(k) - v))
            term_b *= factor_b
            sum_b += term_b
            if abs(term_a) < tol * abs(sum_a) && abs(term_b) < tol * abs(sum_b)
                break
            end
            k += 1
        end
        return T(π) * (sum_b - sum_a)
    end
    z_inv = inv(z)
    sum_terms = one(Complex{T}) +
        (T(2) / T(9)) * z_inv +
        (-T(7) / T(81)) * z_inv^2 +
        (T(175) / T(2187)) * z_inv^3 +
        (-T(2275) / T(19683)) * z_inv^4 +
        (T(5005) / T(177147)) * z_inv^5
    prefactor = sqrt(T(π) / (T(2) * z)) * exp(-z)
    return prefactor * sum_terms
end

@inline _kv56_scalar(z::T) where {T<:AbstractFloat} = _kv56_scalar(complex(z))

function _guide_star_grid(
    sampling::Integer,
    diameter::T,
    rotation_angle_rad::T,
    offset_x::T,
    offset_y::T,
) where {T<:AbstractFloat}
    coords = sampling == 1 ? range(zero(T), zero(T); length=1) : range(-diameter / 2, diameter / 2; length=sampling)
    s, c = sincos(rotation_angle_rad)
    xr = Matrix{T}(undef, sampling, sampling)
    yr = Matrix{T}(undef, sampling, sampling)
    @inbounds for j in 1:sampling
        y = coords[j]
        for i in 1:sampling
            x = coords[i]
            xr[j, i] = x * c - y * s - offset_x * diameter
            yr[j, i] = y * c + x * s - offset_y * diameter
        end
    end
    return xr, yr
end

function _scaled_shifted_coords(
    x::AbstractMatrix{T},
    y::AbstractMatrix{T},
    direction_vectors::AbstractMatrix{T},
    src_index::Int,
    altitude::AbstractVector{T},
    layer_index::Int,
    src_height::T,
) where {T<:AbstractFloat}
    beta_x = direction_vectors[1, src_index] * altitude[layer_index]
    beta_y = direction_vectors[2, src_index] * altitude[layer_index]
    scale = isfinite(src_height) ? one(T) - altitude[layer_index] / src_height : one(T)
    return @. complex(x * scale + beta_x, y * scale + beta_y)
end

function _covariance_matrix(
    ::BuildBackend,
    rho1::AbstractVector{Complex{T}},
    rho2::AbstractVector{Complex{T}},
    r0::T,
    L0::T,
    fractional_r0::T,
) where {T<:AbstractFloat}
    return _covariance_matrix(rho1, rho2, r0, L0, fractional_r0)
end

function _covariance_matrix(
    rho1::AbstractVector{Complex{T}},
    rho2::AbstractVector{Complex{T}},
    r0::T,
    L0::T,
    fractional_r0::T,
) where {T<:AbstractFloat}
    n1 = length(rho1)
    n2 = length(rho2)
    out = Matrix{T}(undef, n1, n2)

    gamma_6_5 = gamma(T(6) / T(5))
    gamma_11_6 = gamma(T(11) / T(6))
    gamma_5_6 = gamma(T(5) / T(6))
    base = (T(24) * gamma_6_5 / T(5))^(T(5) / T(6))
    cst = base * (gamma_11_6 / (T(2)^(T(5) / T(6)) * T(π)^(T(8) / T(3)))) * (L0 / r0)^(T(5) / T(3))
    var_term = base * gamma_11_6 * gamma_5_6 / (T(2) * T(π)^(T(8) / T(3))) * (L0 / r0)^(T(5) / T(3))

    @inbounds for j in 1:n2
        for i in 1:n1
            rho = abs(rho1[i] - rho2[j])
            if iszero(rho)
                out[i, j] = var_term * fractional_r0
            else
                u = T(2π) * rho / L0
                out[i, j] = cst * u^(T(5) / T(6)) * real(_kv56_scalar(u)) * fractional_r0
            end
        end
    end
    return out
end

function _covariance_matrix(
    backend::GPUArrayBuildBackend{B},
    rho1::AbstractVector{Complex{T}},
    rho2::AbstractVector{Complex{T}},
    r0::T,
    L0::T,
    fractional_r0::T,
) where {B,T<:AbstractFloat}
    rho1_native = materialize_build(backend, rho1)
    rho2_native = materialize_build(backend, rho2)
    out = _backend_array(B, T, length(rho1), length(rho2))
    gamma_6_5 = gamma(T(6) / T(5))
    gamma_11_6 = gamma(T(11) / T(6))
    gamma_5_6 = gamma(T(5) / T(6))
    base = (T(24) * gamma_6_5 / T(5))^(T(5) / T(6))
    cst = base * (gamma_11_6 / (T(2)^(T(5) / T(6)) * T(π)^(T(8) / T(3)))) * (L0 / r0)^(T(5) / T(3))
    var_term = base * gamma_11_6 * gamma_5_6 / (T(2) * T(π)^(T(8) / T(3))) * (L0 / r0)^(T(5) / T(3))
    style = execution_style(out)
    launch_kernel!(style, covariance_matrix_kernel!, out, rho1_native, rho2_native, cst, var_term, inv(L0), fractional_r0,
        length(rho1), length(rho2);
        ndrange=size(out))
    return out
end

function sparse_gradient_matrix(
    valid_lenslet::AbstractMatrix{Bool};
    amplitude_mask::Union{Nothing,AbstractMatrix}=nothing,
    over_sampling::Integer=2,
)
    over_sampling == 2 || throw(UnsupportedAlgorithm("only over_sampling=2 is implemented"))
    n_lenslet = size(valid_lenslet, 1)
    size(valid_lenslet, 2) == n_lenslet ||
        throw(DimensionMismatchError("valid_lenslet must be square"))
    n_map = over_sampling * n_lenslet + 1
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
    @inbounds for j_lenslet in 1:n_lenslet
        j_offset = over_sampling * (j_lenslet - 1)
        for i_lenslet in 1:n_lenslet
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

function auto_correlation(
    backend::BuildBackend,
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    return auto_correlation(atmosphere, asterism, wfs, grid_mask)
end

function auto_correlation(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    sampling = size(grid_mask, 1)
    size(grid_mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    mask_vec = vec(grid_mask)
    n_valid = count(mask_vec)
    n_gs = asterism.n_lgs
    result = zeros(T, n_gs * n_valid, n_gs * n_valid)

    altitude = layer_altitude_m(atmosphere)
    r0 = _fried_parameter(atmosphere)
    support_d = support_diameter(wfs)
    lgs_dir = lgs_directions(asterism)
    directions = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    source_height = lgs_height_m(asterism, atmosphere)

    guide_x = Vector{Matrix{T}}(undef, n_gs)
    guide_y = similar(guide_x)
    for gs in 1:n_gs
        guide_x[gs], guide_y[gs] = _guide_star_grid(
            sampling,
            support_d,
            wfs.lenslet_rotation_rad[gs],
            wfs.lenslet_offset[1, gs],
            wfs.lenslet_offset[2, gs],
        )
    end

    for jgs in 1:n_gs
        for igs in 1:jgs
            block = zeros(T, n_valid, n_valid)
            for layer in eachindex(altitude)
                iz = _scaled_shifted_coords(
                    guide_x[igs],
                    guide_y[igs],
                    directions,
                    igs,
                    altitude,
                    layer,
                    source_height,
                )
                jz = _scaled_shifted_coords(
                    guide_x[jgs],
                    guide_y[jgs],
                    directions,
                    jgs,
                    altitude,
                    layer,
                    source_height,
                )
                cov = _covariance_matrix(
                    vec(iz),
                    vec(jz),
                    r0,
                    atmosphere.L0,
                    atmosphere.fractional_r0[layer],
                )
                block .+= @view cov[mask_vec, mask_vec]
            end
            rows = (igs - 1) * n_valid + 1:igs * n_valid
            cols = (jgs - 1) * n_valid + 1:jgs * n_valid
            result[rows, cols] .= block
            if igs != jgs
                result[cols, rows] .= transpose(block)
            end
        end
    end
    return result
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
    mask_vec = vec(grid_mask)
    n_valid = count(mask_vec)
    n_gs = asterism.n_lgs
    result = _backend_array(B, T, n_gs * n_valid, n_gs * n_valid)
    fill!(result, zero(T))

    altitude = layer_altitude_m(atmosphere)
    r0 = _fried_parameter(atmosphere)
    support_d = support_diameter(wfs)
    lgs_dir = lgs_directions(asterism)
    directions = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    source_height = lgs_height_m(asterism, atmosphere)

    guide_x = Vector{Matrix{T}}(undef, n_gs)
    guide_y = similar(guide_x)
    for gs in 1:n_gs
        guide_x[gs], guide_y[gs] = _guide_star_grid(
            sampling,
            support_d,
            wfs.lenslet_rotation_rad[gs],
            wfs.lenslet_offset[1, gs],
            wfs.lenslet_offset[2, gs],
        )
    end

    for jgs in 1:n_gs
        for igs in 1:jgs
            block = _backend_array(B, T, n_valid, n_valid)
            fill!(block, zero(T))
            for layer in eachindex(altitude)
                iz = _scaled_shifted_coords(
                    guide_x[igs],
                    guide_y[igs],
                    directions,
                    igs,
                    altitude,
                    layer,
                    source_height,
                )
                jz = _scaled_shifted_coords(
                    guide_x[jgs],
                    guide_y[jgs],
                    directions,
                    jgs,
                    altitude,
                    layer,
                    source_height,
                )
                cov = _covariance_matrix(
                    backend,
                    vec(iz),
                    vec(jz),
                    r0,
                    atmosphere.L0,
                    atmosphere.fractional_r0[layer],
                )
                block .+= @view cov[mask_vec, mask_vec]
            end
            rows = (igs - 1) * n_valid + 1:igs * n_valid
            cols = (jgs - 1) * n_valid + 1:jgs * n_valid
            result[rows, cols] .= block
            if igs != jgs
                result[cols, rows] .= transpose(block)
            end
        end
    end
    return result
end

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

function cross_correlation(
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {T<:AbstractFloat}
    sampling = isnothing(grid_mask) ? PYTOMOAO_DEFAULT_CROSS_SAMPLING : size(grid_mask, 1)
    mask = isnothing(grid_mask) ? trues(sampling, sampling) : grid_mask
    size(mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    row_mask = vec(mask)
    n_row = count(row_mask)
    n_fit = tomography.n_fit_src^2
    n_gs = asterism.n_lgs
    result = Array{T}(undef, n_fit, n_row, n_gs * n_row)

    altitude = layer_altitude_m(atmosphere)
    r0 = _fried_parameter(atmosphere)
    support_d = support_diameter(wfs)
    lgs_dir = lgs_directions(asterism)
    lgs_directions_xyz = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    fit_zenith, fit_azimuth = optimization_geometry(tomography)
    fit_directions_xyz = direction_vectors(fit_zenith, fit_azimuth)
    source_height = lgs_height_m(asterism, atmosphere)
    target_x, target_y = _guide_star_grid(sampling, support_d, zero(T), zero(T), zero(T))

    guide_x = Vector{Matrix{T}}(undef, n_gs)
    guide_y = similar(guide_x)
    for gs in 1:n_gs
        guide_x[gs], guide_y[gs] = _guide_star_grid(
            sampling,
            support_d,
            wfs.lenslet_rotation_rad[gs],
            wfs.lenslet_offset[1, gs],
            wfs.lenslet_offset[2, gs],
        )
    end

    for fit_idx in 1:n_fit
        for gs in 1:n_gs
            block = zeros(T, n_row, n_row)
            for layer in eachindex(altitude)
                iz = _scaled_shifted_coords(
                    guide_x[gs],
                    guide_y[gs],
                    lgs_directions_xyz,
                    gs,
                    altitude,
                    layer,
                    source_height,
                )
                jz = _scaled_shifted_coords(
                    target_x,
                    target_y,
                    fit_directions_xyz,
                    fit_idx,
                    altitude,
                    layer,
                    tomography.fit_src_height_m,
                )
                cov = _covariance_matrix(
                    vec(iz),
                    vec(jz),
                    r0,
                    atmosphere.L0,
                    atmosphere.fractional_r0[layer],
                )
                block .+= transpose(@view cov[row_mask, row_mask])
            end
            cols = (gs - 1) * n_row + 1:gs * n_row
            result[fit_idx, :, cols] .= block
        end
    end
    return result
end

function cross_correlation(
    backend::GPUArrayBuildBackend{B},
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {B,T<:AbstractFloat}
    sampling = isnothing(grid_mask) ? PYTOMOAO_DEFAULT_CROSS_SAMPLING : size(grid_mask, 1)
    mask = isnothing(grid_mask) ? trues(sampling, sampling) : grid_mask
    size(mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    row_mask = vec(mask)
    n_row = count(row_mask)
    n_fit = tomography.n_fit_src^2
    n_gs = asterism.n_lgs
    result = _backend_array(B, T, n_fit, n_row, n_gs * n_row)

    altitude = layer_altitude_m(atmosphere)
    r0 = _fried_parameter(atmosphere)
    support_d = support_diameter(wfs)
    lgs_dir = lgs_directions(asterism)
    lgs_directions_xyz = direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    fit_zenith, fit_azimuth = optimization_geometry(tomography)
    fit_directions_xyz = direction_vectors(fit_zenith, fit_azimuth)
    source_height = lgs_height_m(asterism, atmosphere)
    target_x, target_y = _guide_star_grid(sampling, support_d, zero(T), zero(T), zero(T))

    guide_x = Vector{Matrix{T}}(undef, n_gs)
    guide_y = similar(guide_x)
    for gs in 1:n_gs
        guide_x[gs], guide_y[gs] = _guide_star_grid(
            sampling,
            support_d,
            wfs.lenslet_rotation_rad[gs],
            wfs.lenslet_offset[1, gs],
            wfs.lenslet_offset[2, gs],
        )
    end

    for fit_idx in 1:n_fit
        for gs in 1:n_gs
            block = _backend_array(B, T, n_row, n_row)
            fill!(block, zero(T))
            for layer in eachindex(altitude)
                iz = _scaled_shifted_coords(
                    guide_x[gs],
                    guide_y[gs],
                    lgs_directions_xyz,
                    gs,
                    altitude,
                    layer,
                    source_height,
                )
                jz = _scaled_shifted_coords(
                    target_x,
                    target_y,
                    fit_directions_xyz,
                    fit_idx,
                    altitude,
                    layer,
                    tomography.fit_src_height_m,
                )
                cov = _covariance_matrix(
                    backend,
                    vec(iz),
                    vec(jz),
                    r0,
                    atmosphere.L0,
                    atmosphere.fractional_r0[layer],
                )
                block .+= transpose(@view cov[row_mask, row_mask])
            end
            cols = (gs - 1) * n_row + 1:gs * n_row
            result[fit_idx, :, cols] .= block
        end
    end
    return result
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
    launch_kernel!(style, fit_source_average_kernel!, out, cross, weights_native, size(cross, 1), size(cross, 2), size(cross, 3);
        ndrange=size(out))
    return out
end

function stable_hermitian_right_division(rhs::AbstractMatrix{T}, css::AbstractMatrix{T}) where {T<:AbstractFloat}
    fact = cholesky(Hermitian(css), check=false)
    if issuccess(fact)
        return transpose(fact \ transpose(rhs))
    end
    return transpose(lu(css) \ transpose(rhs))
end

function stable_hermitian_right_division(
    backend::BuildBackend,
    rhs::AbstractMatrix{T},
    css::AbstractMatrix{T},
) where {T<:AbstractFloat}
    return materialize_build(backend, rhs, stable_hermitian_right_division(
        prepare_build_matrix(backend, rhs),
        prepare_build_matrix(backend, css),
    ))
end

tomography_noise_covariance(model::TomographyNoiseModel, reference_diag::AbstractVector) =
    tomography_noise_covariance(NativeBuildBackend(), model, reference_diag)

function _build_diagonal_noise(::NativeBuildBackend, variances::AbstractVector{T}) where {T<:AbstractFloat}
    return Diagonal(variances)
end

function _build_diagonal_noise(::CPUBuildBackend, variances::AbstractVector{T}) where {T<:AbstractFloat}
    return Diagonal(Vector(variances))
end

function _build_diagonal_noise(backend::GPUArrayBuildBackend{B}, variances::AbstractVector{T}) where {B,T<:AbstractFloat}
    out = _backend_array(B, T, length(variances), length(variances))
    style = execution_style(out)
    launch_kernel!(style, diagonal_matrix_kernel!, out, variances, length(variances); ndrange=size(out))
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
    css_signal = interaction_native * cxx_native * transpose(interaction_native)
    cnz = tomography_noise_covariance(build_backend, noise_model, diag(css_signal))
    css = css_signal .+ cnz
    recstat = stable_hermitian_right_division(build_backend, cox_native * transpose(interaction_native), css)
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

function dm_commands!(
    out::AbstractVector{T},
    reconstructor::TomographyCommandReconstructor{T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    size(reconstructor.matrix, 1) == length(out) ||
        throw(DimensionMismatchError("output length must match command matrix row count"))
    size(reconstructor.matrix, 2) == length(slopes) ||
        throw(DimensionMismatchError("slopes length must match command matrix column count"))
    mul!(out, reconstructor.matrix, slopes)
    return out
end

function dm_commands(
    reconstructor::TomographyCommandReconstructor{T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    out = similar(slopes, T, size(reconstructor.matrix, 1))
    return dm_commands!(out, reconstructor, slopes)
end

function dm_commands!(
    out::AbstractVector{T},
    reconstructor::TomographicReconstructor{InteractionMatrixTomography,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    size(reconstructor.reconstructor, 1) == length(out) ||
        throw(DimensionMismatchError("output length must match reconstructor row count"))
    size(reconstructor.reconstructor, 2) == length(slopes) ||
        throw(DimensionMismatchError("slopes length must match reconstructor column count"))
    mul!(out, reconstructor.reconstructor, slopes)
    return out
end

function dm_commands(
    reconstructor::TomographicReconstructor{InteractionMatrixTomography,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    out = similar(slopes, T, size(reconstructor.reconstructor, 1))
    return dm_commands!(out, reconstructor, slopes)
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
    gamma = blockdiag(ntuple(_ -> gamma_single, asterism.n_lgs)...)
    cxx = auto_correlation(build_backend, atmosphere, asterism, wfs, grid_mask)
    cross = cross_correlation(build_backend, atmosphere, asterism, wfs, tomography)
    cox_full = _fit_source_average(cross, _equal_fit_source_weights(tomography))
    row_positions = findall(vec(grid_mask))
    col_positions = findall(repeat(vec(grid_mask), asterism.n_lgs))
    cox = cox_full[row_positions, col_positions]
    gamma_native = materialize_build(build_backend, gamma, gamma)
    cxx_native = materialize_build(build_backend, gamma_native, cxx)
    cox_native = materialize_build(build_backend, gamma_native, cox)
    native_mask = materialize_build(build_backend, gamma_native, grid_mask)
    css_signal = gamma_native * cxx_native * transpose(gamma_native)
    cnz = tomography_noise_covariance(build_backend, noise_model, diag(css_signal))
    css = css_signal .+ cnz
    recstat = stable_hermitian_right_division(build_backend, cox_native * transpose(gamma_native), css)
    d = support_diameter(wfs) / size(valid_lenslet_support(wfs), 1)
    wavefront_to_meter = asterism.wavelength / d / 2
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

function reconstruct_wavefront!(
    out::AbstractVector{T},
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    size(reconstructor.reconstructor, 1) == length(out) ||
        throw(DimensionMismatchError("output length must match reconstructor row count"))
    size(reconstructor.reconstructor, 2) == length(slopes) ||
        throw(DimensionMismatchError("slopes length must match reconstructor column count"))
    mul!(out, reconstructor.reconstructor, slopes)
    return out
end

function reconstruct_wavefront(
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    out = similar(slopes, T, size(reconstructor.reconstructor, 1))
    return reconstruct_wavefront!(out, reconstructor, slopes)
end

function reconstruct_wavefront_map!(
    out::AbstractMatrix{T},
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T};
    masked_value::T=T(NaN),
) where {T<:AbstractFloat}
    size(out) == size(reconstructor.grid_mask) ||
        throw(DimensionMismatchError("output map size must match reconstructor grid mask"))
    wavefront = reconstruct_wavefront(reconstructor, slopes)
    fill!(out, masked_value)
    out[reconstructor.grid_mask] .= wavefront
    return out
end

function reconstruct_wavefront_map(
    reconstructor::TomographicReconstructor{<:AbstractTomographyMethod,T},
    slopes::AbstractVector{T};
    masked_value::T=T(NaN),
) where {T<:AbstractFloat}
    out = similar(reconstructor.reconstructor, T, size(reconstructor.grid_mask)...)
    return reconstruct_wavefront_map!(out, reconstructor, slopes; masked_value=masked_value)
end
