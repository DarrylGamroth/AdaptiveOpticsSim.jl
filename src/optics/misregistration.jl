struct Misregistration{T<:AbstractFloat}
    shift_x::T
    shift_y::T
    rotation_rad::T
    anamorphosis_angle_rad::T
    tangential_scaling::T
    radial_scaling::T
    transform::SMatrix{2,2,T,4}
end

function Misregistration(; shift_x::Real=0.0, shift_y::Real=0.0, rotation_deg::Real=0.0,
    anamorphosis_angle::Real=0.0, tangential_scaling::Real=1.0, radial_scaling::Real=1.0, T::Type{<:AbstractFloat}=Float64)
    θ = T(deg2rad(anamorphosis_angle))
    sθ, cθ = sincos(θ)
    φ = T(deg2rad(rotation_deg))
    sφ, cφ = sincos(φ)
    tan_scale = T(tangential_scaling)
    rad_scale = T(radial_scaling)
    anamorph = @SMatrix [
        cθ sθ
        -sθ cθ
    ]
    inv_anamorph = @SMatrix [
        cθ -sθ
        sθ cθ
    ]
    scale = @SMatrix [
        tan_scale zero(T)
        zero(T) rad_scale
    ]
    rotation = @SMatrix [
        cφ -sφ
        sφ cφ
    ]

    return Misregistration{T}(
        T(shift_x),
        T(shift_y),
        φ,
        θ,
        tan_scale,
        rad_scale,
        rotation * inv_anamorph * scale * anamorph,
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
    v = mis.transform * SVector(x, y)
    return v[1] - mis.shift_x, v[2] - mis.shift_y
end
