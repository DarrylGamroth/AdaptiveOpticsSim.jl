struct TomographyAtmosphereParams{
    T<:AbstractFloat,
    A<:AbstractVector{T},
    F<:AbstractVector{T},
    D<:AbstractVector{T},
    S<:AbstractVector{T},
}
    zenith_angle_deg::T
    altitude_km::A
    L0::T
    r0_zenith::T
    fractional_r0::F
    wavelength::T
    wind_direction_deg::D
    wind_speed::S
end

function TomographyAtmosphereParams(;
    zenith_angle_deg::Real,
    altitude_km::AbstractVector{<:Real},
    L0::Real,
    r0_zenith::Real,
    fractional_r0::AbstractVector{<:Real},
    wavelength::Real,
    wind_direction_deg::AbstractVector{<:Real},
    wind_speed::AbstractVector{<:Real},
)
    lengths = (
        length(altitude_km),
        length(fractional_r0),
        length(wind_direction_deg),
        length(wind_speed),
    )
    all(==(lengths[1]), lengths) ||
        throw(InvalidConfiguration("tomography atmosphere arrays must have matching lengths"))
    lengths[1] > 0 || throw(InvalidConfiguration("tomography atmosphere must have at least one layer"))
    0 <= zenith_angle_deg <= 90 ||
        throw(InvalidConfiguration("zenith_angle_deg must be between 0 and 90"))
    L0 > 0 || throw(InvalidConfiguration("L0 must be positive"))
    r0_zenith > 0 || throw(InvalidConfiguration("r0_zenith must be positive"))
    wavelength > 0 || throw(InvalidConfiguration("wavelength must be positive"))
    all(>=(0), altitude_km) || throw(InvalidConfiguration("altitude_km must be non-negative"))
    all(>=(0), wind_speed) || throw(InvalidConfiguration("wind_speed must be non-negative"))
    all(>=(0), fractional_r0) ||
        throw(InvalidConfiguration("fractional_r0 must be non-negative"))
    isapprox(sum(fractional_r0), 1; atol=1e-6, rtol=1e-6) ||
        throw(InvalidConfiguration("fractional_r0 must sum to 1"))

    T = promote_type(
        typeof(float(zenith_angle_deg)),
        eltype(float.(altitude_km)),
        typeof(float(L0)),
        typeof(float(r0_zenith)),
        eltype(float.(fractional_r0)),
        typeof(float(wavelength)),
        eltype(float.(wind_direction_deg)),
        eltype(float.(wind_speed)),
    )
    altitude = convert.(T, altitude_km)
    fractions = convert.(T, fractional_r0)
    directions = convert.(T, wind_direction_deg)
    speeds = convert.(T, wind_speed)
    return TomographyAtmosphereParams{
        T,
        typeof(altitude),
        typeof(fractions),
        typeof(directions),
        typeof(speeds),
    }(
        T(zenith_angle_deg),
        altitude,
        T(L0),
        T(r0_zenith),
        fractions,
        T(wavelength),
        directions,
        speeds,
    )
end

struct LGSAsterismParams{T<:AbstractFloat}
    radius_arcsec::T
    wavelength::T
    base_height_m::T
    n_lgs::Int
end

function LGSAsterismParams(;
    radius_arcsec::Real,
    wavelength::Real,
    base_height_m::Real,
    n_lgs::Integer,
)
    radius_arcsec >= 0 || throw(InvalidConfiguration("radius_arcsec must be non-negative"))
    wavelength > 0 || throw(InvalidConfiguration("wavelength must be positive"))
    base_height_m > 0 || throw(InvalidConfiguration("base_height_m must be positive"))
    n_lgs >= 0 || throw(InvalidConfiguration("n_lgs must be non-negative"))
    T = promote_type(
        typeof(float(radius_arcsec)),
        typeof(float(wavelength)),
        typeof(float(base_height_m)),
    )
    return LGSAsterismParams{T}(T(radius_arcsec), T(wavelength), T(base_height_m), Int(n_lgs))
end

struct LGSWFSParams{
    T<:AbstractFloat,
    M<:AbstractMatrix{Bool},
    R<:AbstractVector{T},
    O<:AbstractMatrix{T},
}
    diameter::T
    n_lenslet::Int
    n_px::Int
    field_stop_size_arcsec::T
    valid_lenslet_map::M
    lenslet_rotation_rad::R
    lenslet_offset::O
end

function LGSWFSParams(;
    diameter::Real,
    n_lenslet::Integer,
    n_px::Integer,
    field_stop_size_arcsec::Real,
    valid_lenslet_map::AbstractMatrix,
    lenslet_rotation_rad::AbstractVector{<:Real}=Float64[],
    lenslet_offset::AbstractMatrix{<:Real}=zeros(2, 0),
    n_lgs::Integer=max(length(lenslet_rotation_rad), size(lenslet_offset, 2)),
)
    diameter > 0 || throw(InvalidConfiguration("diameter must be positive"))
    n_lenslet > 0 || throw(InvalidConfiguration("n_lenslet must be positive"))
    n_px > 0 || throw(InvalidConfiguration("n_px must be positive"))
    field_stop_size_arcsec > 0 ||
        throw(InvalidConfiguration("field_stop_size_arcsec must be positive"))
    ndims(valid_lenslet_map) == 2 ||
        throw(InvalidConfiguration("valid_lenslet_map must be 2D"))

    T = promote_type(
        typeof(float(diameter)),
        typeof(float(field_stop_size_arcsec)),
        eltype(float.(lenslet_rotation_rad)),
        eltype(float.(lenslet_offset)),
    )
    lenslet_rotation = convert.(T, lenslet_rotation_rad)
    length(lenslet_rotation) == n_lgs ||
        throw(InvalidConfiguration("lenslet_rotation_rad length must match n_lgs"))
    size(lenslet_offset) == (2, n_lgs) ||
        throw(InvalidConfiguration("lenslet_offset must have size (2, n_lgs)"))
    offsets = convert.(T, lenslet_offset)
    lenslet_map = convert.(Bool, valid_lenslet_map)
    return LGSWFSParams{T, typeof(lenslet_map), typeof(lenslet_rotation), typeof(offsets)}(
        T(diameter),
        Int(n_lenslet),
        Int(n_px),
        T(field_stop_size_arcsec),
        lenslet_map,
        lenslet_rotation,
        offsets,
    )
end

struct TomographyParams{T<:AbstractFloat}
    n_fit_src::Int
    fov_optimization_arcsec::T
    fit_src_height_m::T
end

function TomographyParams(;
    n_fit_src::Integer,
    fov_optimization_arcsec::Real,
    fit_src_height_m::Real=Inf,
)
    n_fit_src > 0 || throw(InvalidConfiguration("n_fit_src must be positive"))
    fov_optimization_arcsec >= 0 ||
        throw(InvalidConfiguration("fov_optimization_arcsec must be non-negative"))
    if n_fit_src > 1 && iszero(fov_optimization_arcsec)
        throw(InvalidConfiguration("fov_optimization_arcsec must be positive when n_fit_src > 1"))
    end
    T = promote_type(typeof(float(fov_optimization_arcsec)), typeof(float(fit_src_height_m)))
    return TomographyParams{T}(Int(n_fit_src), T(fov_optimization_arcsec), T(fit_src_height_m))
end

struct TomographyDMParams{
    T<:AbstractFloat,
    H<:AbstractVector{T},
    P<:AbstractVector{T},
    N<:AbstractVector{Int},
    M<:AbstractMatrix{Bool},
}
    heights_m::H
    pitch_m::P
    cross_coupling::T
    n_actuators::N
    valid_actuators::M
end

function TomographyDMParams(;
    heights_m::AbstractVector{<:Real},
    pitch_m::AbstractVector{<:Real},
    cross_coupling::Real,
    n_actuators::AbstractVector{<:Integer},
    valid_actuators::AbstractMatrix,
)
    length(heights_m) > 0 || throw(InvalidConfiguration("heights_m cannot be empty"))
    length(heights_m) == length(pitch_m) == length(n_actuators) ||
        throw(InvalidConfiguration("DM arrays must have matching lengths"))
    0 <= cross_coupling <= 1 ||
        throw(InvalidConfiguration("cross_coupling must be between 0 and 1"))
    all(>=(0), heights_m) || throw(InvalidConfiguration("heights_m must be non-negative"))
    all(>(0), pitch_m) || throw(InvalidConfiguration("pitch_m must be positive"))
    all(>=(0), n_actuators) || throw(InvalidConfiguration("n_actuators must be non-negative"))
    ndims(valid_actuators) == 2 || throw(InvalidConfiguration("valid_actuators must be 2D"))

    T = promote_type(eltype(float.(heights_m)), eltype(float.(pitch_m)), typeof(float(cross_coupling)))
    heights = convert.(T, heights_m)
    pitch = convert.(T, pitch_m)
    nact = convert.(Int, n_actuators)
    valid = convert.(Bool, valid_actuators)
    return TomographyDMParams{T, typeof(heights), typeof(pitch), typeof(nact), typeof(valid)}(
        heights,
        pitch,
        T(cross_coupling),
        nact,
        valid,
    )
end

@inline airmass(params::TomographyAtmosphereParams) = inv(cosd(params.zenith_angle_deg))

@inline layer_altitude_m(params::TomographyAtmosphereParams) = params.altitude_km .* (1000 * airmass(params))

@inline function wind_direction_rad(params::TomographyAtmosphereParams)
    return deg2rad.(params.wind_direction_deg)
end

@inline function wind_velocity_components(params::TomographyAtmosphereParams)
    direction = wind_direction_rad(params)
    return params.wind_speed .* cos.(direction), params.wind_speed .* sin.(direction)
end

@inline lgs_height_m(params::LGSAsterismParams, atmosphere::TomographyAtmosphereParams) =
    params.base_height_m * airmass(atmosphere)

@inline n_valid_subapertures(params::LGSWFSParams) = count(params.valid_lenslet_map)

function valid_lenslet_support(params::LGSWFSParams)
    support = falses(size(params.valid_lenslet_map) .+ 4)
    @views support[3:end-2, 3:end-2] .= params.valid_lenslet_map
    return support
end

function dm_valid_support(params::TomographyDMParams)
    support = falses(size(params.valid_actuators) .+ 4)
    @views support[3:end-2, 3:end-2] .= params.valid_actuators
    return support
end
