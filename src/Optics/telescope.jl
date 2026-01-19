struct TelescopeParams{T<:AbstractFloat}
    resolution::Int
    diameter::T
    sampling_time::T
    central_obstruction::T
    fov_arcsec::T
end

mutable struct TelescopeState{T, Aopd<:AbstractMatrix{T}, Apsf<:AbstractMatrix{T}, Amask<:AbstractMatrix{Bool}}
    pupil::Amask
    opd::Aopd
    psf::Apsf
    psf_list::Vector{Apsf}
end

struct Telescope{P<:TelescopeParams,S<:TelescopeState} <: AbstractOpticalElement
    params::P
    state::S
end

function Telescope(; resolution::Int,
    diameter::Real,
    sampling_time::Real,
    central_obstruction::Real=0.0,
    fov_arcsec::Real=0.0,
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)

    params = TelescopeParams{T}(
        resolution,
        T(diameter),
        T(sampling_time),
        T(central_obstruction),
        T(fov_arcsec),
    )

    pupil = backend{Bool}(undef, resolution, resolution)
    generate_pupil!(pupil, params)

    opd = backend{T}(undef, resolution, resolution)
    fill!(opd, zero(T))

    psf = backend{T}(undef, resolution, resolution)
    fill!(psf, zero(T))

    state = TelescopeState{T, typeof(opd), typeof(psf), typeof(pupil)}(pupil, opd, psf, Vector{typeof(psf)}())
    return Telescope(params, state)
end

function generate_pupil!(pupil::AbstractMatrix{Bool}, params::TelescopeParams)
    Base.require_one_based_indexing(pupil)
    n = params.resolution
    r_outer = 1.0
    r_inner = params.central_obstruction
    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2

    @inbounds for i in 1:n, j in 1:n
        x = (i - cx) / scale
        y = (j - cy) / scale
        r = sqrt(x^2 + y^2)
        pupil[i, j] = (r <= r_outer) & (r >= r_inner)
    end
    return pupil
end

function reset_opd!(tel::Telescope)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    return tel
end

function apply_opd!(tel::Telescope, opd::AbstractMatrix)
    if size(opd) != size(tel.state.opd)
        throw(DimensionMismatchError("OPD size does not match telescope resolution"))
    end
    tel.state.opd .= opd
    return tel
end

function set_pupil!(tel::Telescope, pupil::AbstractMatrix{Bool})
    if size(pupil) != size(tel.state.pupil)
        throw(DimensionMismatchError("pupil mask size does not match telescope resolution"))
    end
    tel.state.pupil .= pupil
    return tel
end

function apply_spiders!(tel::Telescope; thickness::Real, angles::AbstractVector, offset_x::Real=0.0, offset_y::Real=0.0)
    Base.require_one_based_indexing(tel.state.pupil)
    n = tel.params.resolution
    radius = tel.params.diameter / 2
    thickness_norm = thickness / radius
    offset_x_norm = offset_x / radius
    offset_y_norm = offset_y / radius

    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2

    for angle in angles
        θ = deg2rad(angle)
        a = -sin(θ)
        b = cos(θ)
        @inbounds for i in 1:n, j in 1:n
            x = (i - cx) / scale - offset_x_norm
            y = (j - cy) / scale - offset_y_norm
            dist = abs(a * x + b * y)
            if dist <= thickness_norm
                tel.state.pupil[i, j] = false
            end
        end
    end
    return tel
end
