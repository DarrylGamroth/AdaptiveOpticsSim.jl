struct Misregistration{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    rotation_deg::T
    anamorphosis_angle::T
    tangential_scaling::T
    radial_scaling::T
end

function Misregistration(; shift_x::Real=0.0, shift_y::Real=0.0, rotation_deg::Real=0.0,
    anamorphosis_angle::Real=0.0, tangential_scaling::Real=1.0, radial_scaling::Real=1.0, T::Type{<:AbstractFloat}=Float64)

    return Misregistration{T}(
        T(shift_x),
        T(shift_y),
        T(rotation_deg),
        T(anamorphosis_angle),
        T(tangential_scaling),
        T(radial_scaling),
    )
end

function apply_misregistration(x::Real, y::Real, mis::Misregistration)
    θ = deg2rad(mis.anamorphosis_angle)
    x1 = cos(θ) * x + sin(θ) * y
    y1 = -sin(θ) * x + cos(θ) * y

    x1 *= mis.tangential_scaling
    y1 *= mis.radial_scaling

    x2 = cos(θ) * x1 - sin(θ) * y1
    y2 = sin(θ) * x1 + cos(θ) * y1

    φ = deg2rad(mis.rotation_deg)
    xr = cos(φ) * x2 - sin(φ) * y2
    yr = sin(φ) * x2 + cos(φ) * y2

    return xr - mis.shift_x, yr - mis.shift_y
end
