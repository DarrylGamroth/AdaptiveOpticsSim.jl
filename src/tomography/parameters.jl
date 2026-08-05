struct TomographyAtmosphereParams{
    T<:AbstractFloat,
    A<:AbstractVector{T},
    F<:AbstractVector{T},
    D<:AbstractVector{T},
    S<:AbstractVector{T},
}
    zenith_angle_rad::T
    layer_altitudes_m::A
    L0::T
    r0_zenith::T
    fractional_cn2::F
    reference_wavelength_m::T
    wind_direction_rad::D
    wind_speed::S
end

function TomographyAtmosphereParams(;
    zenith_angle_deg::Real,
    layer_altitudes_m::AbstractVector{<:Real},
    L0::Real,
    r0_zenith::Real,
    fractional_cn2::AbstractVector{<:Real},
    reference_wavelength_m::Real,
    wind_direction_deg::AbstractVector{<:Real},
    wind_speed::AbstractVector{<:Real},
)
    lengths = (
        length(layer_altitudes_m),
        length(fractional_cn2),
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
    reference_wavelength_m > 0 ||
        throw(InvalidConfiguration("reference_wavelength_m must be positive"))
    all(>=(0), layer_altitudes_m) ||
        throw(InvalidConfiguration("layer_altitudes_m must be non-negative"))
    all(>=(0), wind_speed) || throw(InvalidConfiguration("wind_speed must be non-negative"))
    all(>=(0), fractional_cn2) ||
        throw(InvalidConfiguration("fractional_cn2 must be non-negative"))
    isapprox(sum(fractional_cn2), 1; atol=1e-6, rtol=1e-6) ||
        throw(InvalidConfiguration("fractional_cn2 must sum to 1"))

    T = promote_type(
        typeof(float(zenith_angle_deg)),
        eltype(float.(layer_altitudes_m)),
        typeof(float(L0)),
        typeof(float(r0_zenith)),
        eltype(float.(fractional_cn2)),
        typeof(float(reference_wavelength_m)),
        eltype(float.(wind_direction_deg)),
        eltype(float.(wind_speed)),
    )
    altitudes = convert.(T, layer_altitudes_m)
    fractions = convert.(T, fractional_cn2)
    directions = T.(deg2rad.(wind_direction_deg))
    speeds = convert.(T, wind_speed)
    return TomographyAtmosphereParams{
        T,
        typeof(altitudes),
        typeof(fractions),
        typeof(directions),
        typeof(speeds),
    }(
        T(deg2rad(zenith_angle_deg)),
        altitudes,
        T(L0),
        T(r0_zenith),
        fractions,
        T(reference_wavelength_m),
        directions,
        speeds,
    )
end

struct LGSAsterismParams{T<:AbstractFloat}
    radius_arcsec::T
    wavelength_m::T
    base_height_m::T
    n_lgs::Int
end

function LGSAsterismParams(;
    radius_arcsec::Real,
    wavelength_m::Real,
    base_height_m::Real,
    n_lgs::Integer,
)
    radius_arcsec >= 0 || throw(InvalidConfiguration("radius_arcsec must be non-negative"))
    wavelength_m > 0 || throw(InvalidConfiguration("wavelength_m must be positive"))
    base_height_m > 0 || throw(InvalidConfiguration("base_height_m must be positive"))
    n_lgs >= 0 || throw(InvalidConfiguration("n_lgs must be non-negative"))
    T = promote_type(
        typeof(float(radius_arcsec)),
        typeof(float(wavelength_m)),
        typeof(float(base_height_m)),
    )
    return LGSAsterismParams{T}(T(radius_arcsec), T(wavelength_m), T(base_height_m), Int(n_lgs))
end

struct LGSWFSParams{
    T<:AbstractFloat,
    M<:AbstractMatrix{Bool},
    R<:AbstractVector{T},
    O<:AbstractMatrix{T},
}
    pupil_diameter_m::T
    n_lenslets::Int
    n_px::Int
    field_stop_size_arcsec::T
    valid_lenslet_map::M
    lenslet_grid_rotations_rad::R
    lenslet_grid_offsets_fraction::O
end

function LGSWFSParams(;
    pupil_diameter_m::Real,
    n_lenslets::Integer,
    n_px::Integer,
    field_stop_size_arcsec::Real,
    valid_lenslet_map::AbstractMatrix,
    lenslet_grid_rotations_rad::AbstractVector{<:Real}=Float64[],
    lenslet_grid_offsets_fraction::AbstractMatrix{<:Real}=zeros(2, 0),
    n_lgs::Integer=max(length(lenslet_grid_rotations_rad),
        size(lenslet_grid_offsets_fraction, 2)),
)
    pupil_diameter_m > 0 ||
        throw(InvalidConfiguration("pupil_diameter_m must be positive"))
    n_lenslets > 0 || throw(InvalidConfiguration("n_lenslets must be positive"))
    n_px > 0 || throw(InvalidConfiguration("n_px must be positive"))
    field_stop_size_arcsec > 0 ||
        throw(InvalidConfiguration("field_stop_size_arcsec must be positive"))
    ndims(valid_lenslet_map) == 2 ||
        throw(InvalidConfiguration("valid_lenslet_map must be 2D"))

    T = promote_type(
        typeof(float(pupil_diameter_m)),
        typeof(float(field_stop_size_arcsec)),
        eltype(float.(lenslet_grid_rotations_rad)),
        eltype(float.(lenslet_grid_offsets_fraction)),
    )
    lenslet_grid_rotations = convert.(T, lenslet_grid_rotations_rad)
    length(lenslet_grid_rotations) == n_lgs || throw(InvalidConfiguration(
        "lenslet_grid_rotations_rad length must match n_lgs"))
    size(lenslet_grid_offsets_fraction) == (2, n_lgs) ||
        throw(InvalidConfiguration(
            "lenslet_grid_offsets_fraction must have size (2, n_lgs)"))
    lenslet_grid_offsets = convert.(T, lenslet_grid_offsets_fraction)
    lenslet_map = convert.(Bool, valid_lenslet_map)
    return LGSWFSParams{T, typeof(lenslet_map), typeof(lenslet_grid_rotations),
        typeof(lenslet_grid_offsets)}(
        T(pupil_diameter_m),
        Int(n_lenslets),
        Int(n_px),
        T(field_stop_size_arcsec),
        lenslet_map,
        lenslet_grid_rotations,
        lenslet_grid_offsets,
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

@inline zenith_angle_rad(params::TomographyAtmosphereParams) = params.zenith_angle_rad
@inline zenith_angle_deg(params::TomographyAtmosphereParams) = rad2deg(params.zenith_angle_rad)

@inline airmass(params::TomographyAtmosphereParams) = inv(cos(params.zenith_angle_rad))

@inline layer_slant_ranges_m(params::TomographyAtmosphereParams) =
    params.layer_altitudes_m .* airmass(params)

@inline wind_direction_rad(params::TomographyAtmosphereParams) = params.wind_direction_rad
@inline wind_direction_deg(params::TomographyAtmosphereParams) = rad2deg.(params.wind_direction_rad)

@inline function wind_velocity_components(params::TomographyAtmosphereParams)
    return params.wind_speed .* cos.(params.wind_direction_rad),
        params.wind_speed .* sin.(params.wind_direction_rad)
end

@inline lgs_height_m(params::LGSAsterismParams, atmosphere::TomographyAtmosphereParams) =
    params.base_height_m * airmass(atmosphere)

@inline n_valid_subapertures(params::LGSWFSParams) = count(params.valid_lenslet_map)

function valid_lenslet_support(params::LGSWFSParams)
    support = falses(size(params.valid_lenslet_map) .+ 4)
    @views support[3:end-2, 3:end-2] .= params.valid_lenslet_map
    return support
end

@inline function lenslet_grid_support_diameter_m(params::LGSWFSParams)
    return params.pupil_diameter_m * size(valid_lenslet_support(params), 1) /
        params.n_lenslets
end

function dm_valid_support(params::TomographyDMParams)
    support = falses(size(params.valid_actuators) .+ 4)
    @views support[3:end-2, 3:end-2] .= params.valid_actuators
    return support
end
