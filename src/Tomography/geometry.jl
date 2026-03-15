const ARCSEC_TO_RAD = π / (180 * 3600)

function optimization_geometry!(
    zenith::AbstractVector{T},
    azimuth::AbstractVector{T},
    params::TomographyParams{T},
) where {T<:AbstractFloat}
    n = params.n_fit_src
    length(zenith) == n^2 ||
        throw(DimensionMismatchError("zenith output length must equal n_fit_src^2"))
    length(azimuth) == n^2 ||
        throw(DimensionMismatchError("azimuth output length must equal n_fit_src^2"))
    if n == 1
        zenith[1] = zero(T)
        azimuth[1] = zero(T)
        return zenith, azimuth
    end

    coords = range(-params.fov_optimization_arcsec / 2, params.fov_optimization_arcsec / 2; length=n)
    k = 1
    @inbounds for j in 1:n, i in 1:n
        x = coords[i]
        y = coords[j]
        zenith[k] = hypot(x, y) * ARCSEC_TO_RAD
        azimuth[k] = atan(y, x)
        k += 1
    end
    return zenith, azimuth
end

function optimization_geometry(params::TomographyParams{T}) where {T<:AbstractFloat}
    zenith = Vector{T}(undef, params.n_fit_src^2)
    azimuth = similar(zenith)
    optimization_geometry!(zenith, azimuth, params)
    return zenith, azimuth
end

function direction_vectors!(
    out::AbstractMatrix{T},
    zenith::AbstractVector{T},
    azimuth::AbstractVector{T},
) where {T<:AbstractFloat}
    size(out, 1) == 3 || throw(DimensionMismatchError("direction vector output must have 3 rows"))
    size(out, 2) == length(zenith) == length(azimuth) ||
        throw(DimensionMismatchError("direction vector output width must match angle vector lengths"))
    @inbounds for k in eachindex(zenith, azimuth)
        tangent = tan(zenith[k])
        out[1, k] = tangent * cos(azimuth[k])
        out[2, k] = tangent * sin(azimuth[k])
        out[3, k] = one(T)
    end
    return out
end

function direction_vectors(zenith::AbstractVector{T}, azimuth::AbstractVector{T}) where {T<:AbstractFloat}
    out = Matrix{T}(undef, 3, length(zenith))
    return direction_vectors!(out, zenith, azimuth)
end

function lgs_directions!(out::AbstractMatrix{T}, params::LGSAsterismParams{T}) where {T<:AbstractFloat}
    size(out) == (params.n_lgs, 2) ||
        throw(DimensionMismatchError("LGS directions output must have size (n_lgs, 2)"))
    radius = params.radius_arcsec * ARCSEC_TO_RAD
    @inbounds for k in 1:params.n_lgs
        out[k, 1] = radius
        out[k, 2] = T((k - 1) * 2π / params.n_lgs)
    end
    return out
end

function lgs_directions(params::LGSAsterismParams{T}) where {T<:AbstractFloat}
    out = Matrix{T}(undef, params.n_lgs, 2)
    return lgs_directions!(out, params)
end
