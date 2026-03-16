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

function apply_misregistration(mis::Misregistration, x::Real, y::Real)
    θ = deg2rad(mis.anamorphosis_angle)
    sθ, cθ = sincos(θ)
    x1 = cθ * x + sθ * y
    y1 = -sθ * x + cθ * y

    x1 *= mis.tangential_scaling
    y1 *= mis.radial_scaling

    x2 = cθ * x1 - sθ * y1
    y2 = sθ * x1 + cθ * y1

    φ = deg2rad(mis.rotation_deg)
    sφ, cφ = sincos(φ)
    xr = cφ * x2 - sφ * y2
    yr = sφ * x2 + cφ * y2

    return xr - mis.shift_x, yr - mis.shift_y
end

apply_misregistration(x::Real, y::Real, mis::Misregistration) = apply_misregistration(mis, x, y)
