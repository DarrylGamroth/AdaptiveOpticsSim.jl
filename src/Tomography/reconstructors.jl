using LinearAlgebra
using SparseArrays

abstract type AbstractTomographyMethod end

struct ModelBasedTomography <: AbstractTomographyMethod end

struct InteractionMatrixTomography <: AbstractTomographyMethod end

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

@inline _fried_parameter(params::TomographyAtmosphereParams{T}) where {T<:AbstractFloat} =
    params.r0_zenith * cosd(params.zenith_angle_deg)^(T(3) / T(5))

@inline _equal_fit_source_weights(params::TomographyParams{T}) where {T<:AbstractFloat} =
    fill(inv(T(params.n_fit_src^2)), params.n_fit_src^2)

function _guide_star_grid(
    sampling::Integer,
    diameter::T,
    rotation_angle_rad::T,
    offset_x::T,
    offset_y::T,
) where {T<:AbstractFloat}
    coords = if sampling == 1
        T[zero(T)]
    else
        collect(range(-one(T), one(T); length=sampling) .* (diameter / 2))
    end
    x = repeat(reshape(coords, 1, :), sampling, 1)
    y = repeat(reshape(coords, :, 1), 1, sampling)
    c = cos(rotation_angle_rad)
    s = sin(rotation_angle_rad)
    xvec = vec(x)
    yvec = vec(y)
    xr = similar(xvec)
    yr = similar(yvec)
    @inbounds for k in eachindex(xvec, yvec)
        xr[k] = xvec[k] * c - yvec[k] * s - offset_x * diameter
        yr[k] = yvec[k] * c + xvec[k] * s - offset_y * diameter
    end
    return reshape(xr, sampling, sampling), reshape(yr, sampling, sampling)
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
                out[i, j] = cst * u^(T(5) / T(6)) * besselk(T(5) / T(6), u) * fractional_r0
            end
        end
    end
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
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    grid_mask::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    sampling = size(grid_mask, 1)
    size(grid_mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    mask_vec = vec(permutedims(grid_mask))
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
                block .+= transpose(@view cov[mask_vec, mask_vec])
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
    atmosphere::TomographyAtmosphereParams{T},
    asterism::LGSAsterismParams{T},
    wfs::LGSWFSParams{T},
    tomography::TomographyParams{T};
    grid_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
) where {T<:AbstractFloat}
    sampling = isnothing(grid_mask) ? 2 * size(valid_lenslet_support(wfs), 1) + 1 : size(grid_mask, 1)
    mask = isnothing(grid_mask) ? trues(sampling, sampling) : grid_mask
    size(mask, 2) == sampling || throw(DimensionMismatchError("grid_mask must be square"))
    row_mask = vec(permutedims(mask))
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

function _fit_source_average(cross::Array{T,3}, weights::AbstractVector{T}) where {T<:AbstractFloat}
    size(cross, 1) == length(weights) ||
        throw(DimensionMismatchError("fit-source weight length must match cross-correlation stack"))
    out = zeros(T, size(cross, 2), size(cross, 3))
    for i in axes(cross, 1)
        @views out .+= weights[i] .* cross[i, :, :]
    end
    return out
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
) where {T<:AbstractFloat}
    size(interaction_matrix, 1) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one row"))
    size(interaction_matrix, 2) > 0 ||
        throw(InvalidConfiguration("interaction_matrix must have at least one column"))

    cxx = auto_correlation(atmosphere, asterism, wfs, grid_mask)
    cross = cross_correlation(atmosphere, asterism, wfs, tomography; grid_mask=grid_mask)
    cox = _fit_source_average(cross, _equal_fit_source_weights(tomography))
    cnz = Diagonal(fill(T(1e-3 * α), size(interaction_matrix, 1)))
    css = Matrix(interaction_matrix * cxx * transpose(interaction_matrix) .+ cnz)
    recstat = (cox * transpose(interaction_matrix)) / css
    operators = TomographyOperators(
        nothing,
        Matrix(grid_mask),
        cxx,
        cox,
        cnz,
        recstat,
        one(T),
    )
    return TomographicReconstructor(
        InteractionMatrixTomography(),
        recstat,
        Matrix(grid_mask),
        atmosphere,
        asterism,
        wfs,
        tomography,
        dm,
        fitting,
        operators,
    )
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
) where {T<:AbstractFloat}
    gamma_single, grid_mask = sparse_gradient_matrix(valid_lenslet_support(wfs); over_sampling=2)
    gamma = blockdiag(ntuple(_ -> gamma_single, asterism.n_lgs)...)
    cxx = auto_correlation(atmosphere, asterism, wfs, grid_mask)
    cross = cross_correlation(atmosphere, asterism, wfs, tomography)
    cox_full = _fit_source_average(cross, _equal_fit_source_weights(tomography))
    row_mask = vec(permutedims(grid_mask))
    col_mask = repeat(row_mask, asterism.n_lgs)
    cox = cox_full[row_mask, col_mask]
    css_signal = Matrix(gamma * cxx * transpose(gamma))
    diagonal = diag(css_signal)
    mean_diag = sum(diagonal) / length(diagonal)
    cnz = Diagonal(fill(mean_diag / 10, size(gamma, 1)))
    css = css_signal .+ cnz
    recstat = (cox * transpose(gamma)) / css
    d = support_diameter(wfs) / size(valid_lenslet_support(wfs), 1)
    wavefront_to_meter = asterism.wavelength / d / 2
    recon = d * wavefront_to_meter .* recstat
    operators = TomographyOperators(
        gamma,
        grid_mask,
        cxx,
        cox,
        cnz,
        recstat,
        wavefront_to_meter,
    )
    return TomographicReconstructor(
        ModelBasedTomography(),
        recon,
        grid_mask,
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
    out = Vector{T}(undef, size(reconstructor.reconstructor, 1))
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
    out = Matrix{T}(undef, size(reconstructor.grid_mask)...)
    return reconstruct_wavefront_map!(out, reconstructor, slopes; masked_value=masked_value)
end
