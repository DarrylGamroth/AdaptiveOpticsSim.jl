struct Misregistration{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    rotation_rad::T
    anamorphosis_angle_rad::T
    tangential_scaling::T
    radial_scaling::T
    transform::NTuple{4,T}
end

function Misregistration(; shift_x::Real=0.0, shift_y::Real=0.0, rotation_deg::Real=0.0,
    anamorphosis_angle::Real=0.0, tangential_scaling::Real=1.0, radial_scaling::Real=1.0, T::Type{<:AbstractFloat}=Float64)
    θ = T(deg2rad(anamorphosis_angle))
    sθ, cθ = sincos(θ)
    φ = T(deg2rad(rotation_deg))
    sφ, cφ = sincos(φ)
    tan_scale = T(tangential_scaling)
    rad_scale = T(radial_scaling)

    b11 = tan_scale * cθ * cθ + rad_scale * sθ * sθ
    b12 = (tan_scale - rad_scale) * sθ * cθ
    b21 = b12
    b22 = tan_scale * sθ * sθ + rad_scale * cθ * cθ
    transform = (
        cφ * b11 - sφ * b21,
        cφ * b12 - sφ * b22,
        sφ * b11 + cφ * b21,
        sφ * b12 + cφ * b22,
    )

    return Misregistration{T}(
        T(shift_x),
        T(shift_y),
        φ,
        θ,
        tan_scale,
        rad_scale,
        transform,
    )
end

@inline rotation_rad(mis::Misregistration) = mis.rotation_rad
@inline rotation_deg(mis::Misregistration) = rad2deg(mis.rotation_rad)
@inline anamorphosis_angle_rad(mis::Misregistration) = mis.anamorphosis_angle_rad
@inline anamorphosis_angle_deg(mis::Misregistration) = rad2deg(mis.anamorphosis_angle_rad)

@inline function misregistration_component(mis::Misregistration, field::Symbol)
    if field === :shift_x
        return mis.shift_x
    elseif field === :shift_y
        return mis.shift_y
    elseif field === :rotation_deg
        return rotation_deg(mis)
    elseif field === :anamorphosis_angle
        return anamorphosis_angle_deg(mis)
    elseif field === :tangential_scaling
        return mis.tangential_scaling
    elseif field === :radial_scaling
        return mis.radial_scaling
    end
    throw(InvalidConfiguration("unsupported misregistration field $(field)"))
end

@inline function apply_misregistration(mis::Misregistration, x::Real, y::Real)
    m11, m12, m21, m22 = mis.transform
    xr = muladd(m12, y, m11 * x)
    yr = muladd(m22, y, m21 * x)
    return xr - mis.shift_x, yr - mis.shift_y
end
