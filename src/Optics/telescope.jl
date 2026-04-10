struct TelescopeParams{T<:AbstractFloat}
    resolution::Int
    diameter::T
    sampling_time::T
    central_obstruction::T
    fov_arcsec::T
end

mutable struct TelescopeState{T,
    Aopd<:AbstractMatrix{T},
    Apsf<:AbstractMatrix{T},
    Amask<:AbstractMatrix{Bool},
    Aref<:AbstractMatrix{T},
    Spsf<:AbstractArray{T,3},
    W<:Workspace}
    pupil::Amask
    pupil_reflectivity::Aref
    opd::Aopd
    psf::Apsf
    psf_stack::Spsf
    psf_workspace::W
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
    pupil_reflectivity::Union{Real,AbstractMatrix}=1.0,
    T::Type{<:AbstractFloat}=Float64,
    backend=CPUBackend())

    backend = _resolve_array_backend(backend)

    params = TelescopeParams{T}(
        resolution,
        T(diameter),
        T(sampling_time),
        T(central_obstruction),
        T(fov_arcsec),
    )

    pupil = backend{Bool}(undef, resolution, resolution)
    generate_pupil!(pupil, params)
    reflectivity = initialize_reflectivity(pupil, pupil_reflectivity, T, backend)

    opd = backend{T}(undef, resolution, resolution)
    fill!(opd, zero(T))

    psf = backend{T}(undef, resolution, resolution)
    fill!(psf, zero(T))

    psf_stack = backend{T}(undef, resolution, resolution, 0)
    psf_workspace = Workspace(opd, resolution; T=T)
    state = TelescopeState{T, typeof(opd), typeof(psf), typeof(pupil), typeof(reflectivity), typeof(psf_stack), typeof(psf_workspace)}(
        pupil,
        reflectivity,
        opd,
        psf,
        psf_stack,
        psf_workspace,
    )
    return Telescope(params, state)
end

function generate_pupil!(pupil::AbstractMatrix{Bool}, params::TelescopeParams)
    Base.require_one_based_indexing(pupil)
    _generate_pupil!(execution_style(pupil), pupil, params)
    return pupil
end

function _generate_pupil!(::ScalarCPUStyle, pupil::AbstractMatrix{Bool}, params::TelescopeParams)
    build_mask!(pupil, AnnularAperture(inner_radius=params.central_obstruction, outer_radius=one(params.diameter), T=typeof(params.diameter));
        grid=default_mask_grid(pupil; T=typeof(params.diameter)))
    return pupil
end

function _generate_pupil!(style::AcceleratorStyle, pupil::AbstractMatrix{Bool}, params::TelescopeParams)
    build_mask!(pupil, AnnularAperture(inner_radius=params.central_obstruction, outer_radius=one(params.diameter), T=typeof(params.diameter));
        grid=default_mask_grid(pupil; T=typeof(params.diameter)))
    return pupil
end

function reset_opd!(tel::Telescope)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    return tel
end

function initialize_reflectivity(pupil::AbstractMatrix{Bool}, reflectivity::Real, ::Type{T}, backend) where {T<:AbstractFloat}
    out = backend{T}(undef, size(pupil)...)
    fill!(out, T(reflectivity))
    out .*= pupil
    return out
end

function initialize_reflectivity(pupil::AbstractMatrix{Bool}, reflectivity::AbstractMatrix, ::Type{T}, backend) where {T<:AbstractFloat}
    if size(reflectivity) != size(pupil)
        throw(DimensionMismatchError("pupil_reflectivity size does not match telescope resolution"))
    end
    out = backend{T}(undef, size(pupil)...)
    copyto!(out, T.(reflectivity))
    out .*= pupil
    return out
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
    tel.state.pupil_reflectivity .= pupil
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::Real)
    fill!(tel.state.pupil_reflectivity, eltype(tel.state.pupil_reflectivity)(reflectivity))
    tel.state.pupil_reflectivity .*= tel.state.pupil
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::AbstractMatrix)
    if size(reflectivity) != size(tel.state.pupil_reflectivity)
        throw(DimensionMismatchError("pupil_reflectivity size does not match telescope resolution"))
    end
    tel.state.pupil_reflectivity .= reflectivity
    tel.state.pupil_reflectivity .*= tel.state.pupil
    return tel
end

function flux_map(tel::Telescope, src::AbstractSource)
    scale = photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    out = similar(tel.state.pupil_reflectivity)
    @. out = scale * tel.state.pupil_reflectivity
    return out
end

function apply_spiders!(tel::Telescope; thickness::Real, angles::AbstractVector, offset_x::Real=0.0, offset_y::Real=0.0)
    Base.require_one_based_indexing(tel.state.pupil)
    radius = tel.params.diameter / 2
    thickness_norm = thickness / radius
    offset_x_norm = offset_x / radius
    offset_y_norm = offset_y / radius
    _apply_spiders!(tel.state.pupil, angles, thickness_norm, offset_x_norm, offset_y_norm)
    return tel
end

function _apply_spiders!(pupil::AbstractMatrix{Bool}, angles::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real)
    T = promote_type(typeof(thickness_norm), typeof(offset_x_norm), typeof(offset_y_norm))
    grid = default_mask_grid(pupil; T=T)
    for angle in angles
        apply_mask!(pupil, SpiderMask(thickness=thickness_norm, angle_rad=deg2rad(angle), offset_x=offset_x_norm, offset_y=offset_y_norm,
            T=T); grid=grid)
    end
    return pupil
end

function _apply_spiders!(::ScalarCPUStyle, pupil::AbstractMatrix{Bool}, angles::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real, cx::Real, cy::Real, scale::Real, n::Int)
    _apply_spiders!(pupil, angles, thickness_norm, offset_x_norm, offset_y_norm)
    return pupil
end

function _apply_spiders!(::AcceleratorStyle, pupil::AbstractMatrix{Bool}, angles::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real, cx::Real, cy::Real, scale::Real, n::Int)
    _apply_spiders!(pupil, angles, thickness_norm, offset_x_norm, offset_y_norm)
    return pupil
end
