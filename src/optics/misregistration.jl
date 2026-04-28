struct Misregistration{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    rotation_rad::T
    anamorphosis_angle_rad::T
    tangential_scaling::T
    radial_scaling::T
end

function Misregistration(; shift_x::Real=0.0, shift_y::Real=0.0, rotation_deg::Real=0.0,
    anamorphosis_angle::Real=0.0, tangential_scaling::Real=1.0, radial_scaling::Real=1.0, T::Type{<:AbstractFloat}=Float64)

    return Misregistration{T}(
        T(shift_x),
        T(shift_y),
        T(deg2rad(rotation_deg)),
        T(deg2rad(anamorphosis_angle)),
        T(tangential_scaling),
        T(radial_scaling),
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

function apply_misregistration(mis::Misregistration, x::Real, y::Real)
    θ = mis.anamorphosis_angle_rad
    sθ, cθ = sincos(θ)
    x1 = cθ * x + sθ * y
    y1 = -sθ * x + cθ * y

    x1 *= mis.tangential_scaling
    y1 *= mis.radial_scaling

    x2 = cθ * x1 - sθ * y1
    y2 = sθ * x1 + cθ * y1

    φ = mis.rotation_rad
    sφ, cφ = sincos(φ)
    xr = cφ * x2 - sφ * y2
    yr = sφ * x2 + cφ * y2

    return xr - mis.shift_x, yr - mis.shift_y
end
