struct TelescopeParams{T<:AbstractFloat}
    resolution::Int
    diameter::T
    sampling_time::T
    central_obstruction::T
    fov_arcsec::T
end

struct TelescopeAperture{T<:AbstractFloat,
    Amask<:AbstractMatrix{Bool},
    Aref<:AbstractMatrix{T}}
    pupil::Amask
    reflectivity::Aref
    sampling_m::NTuple{2,T}
    origin_m::NTuple{2,T}
end

# Transitional storage used by simulation paths that have not yet migrated to
# caller-owned optical products. Gate 0 removes these fields from `Telescope`
# as the atmosphere and WFS families move to the explicit plane API.
mutable struct LegacyTelescopePathState{T<:AbstractFloat,
    Aopd<:AbstractMatrix{T},
    Apsf<:AbstractMatrix{T},
    Spsf<:AbstractArray{T,3},
    W<:Workspace}
    opd::Aopd
    psf::Apsf
    psf_stack::Spsf
    psf_workspace::W
end

struct Telescope{P<:TelescopeParams,A<:TelescopeAperture,S<:LegacyTelescopePathState,
    B<:AbstractArrayBackend} <: AbstractTelescope
    params::P
    aperture::A
    state::S
end

@inline backend(::Telescope{<:Any,<:Any,<:Any,B}) where {B} = B()
@inline pupil_mask(tel::Telescope) = tel.aperture.pupil
@inline pupil_reflectivity(tel::Telescope) = tel.aperture.reflectivity
@inline opd_map(tel::Telescope) = tel.state.opd

function Telescope(; resolution::Int,
    diameter::Real,
    sampling_time::Real,
    central_obstruction::Real=0.0,
    fov_arcsec::Real=0.0,
    pupil_reflectivity::Union{Real,AbstractMatrix}=1.0,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())

    selector = _resolve_backend_selector(backend)
    backend = _resolve_array_backend(selector)

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
    sampling_m = T(diameter) / T(resolution)
    origin = -T(resolution - 1) * sampling_m / T(2)
    aperture = TelescopeAperture{T,typeof(pupil),typeof(reflectivity)}(
        pupil,
        reflectivity,
        (sampling_m, sampling_m),
        (origin, origin),
    )

    opd = backend{T}(undef, resolution, resolution)
    fill!(opd, zero(T))

    psf = backend{T}(undef, resolution, resolution)
    fill!(psf, zero(T))

    psf_stack = backend{T}(undef, resolution, resolution, 0)
    psf_workspace = Workspace(opd, resolution; T=T)
    state = LegacyTelescopePathState{T,typeof(opd),typeof(psf),
        typeof(psf_stack),typeof(psf_workspace)}(
        opd,
        psf,
        psf_stack,
        psf_workspace,
    )
    return Telescope{typeof(params),typeof(aperture),typeof(state),
        typeof(selector)}(params, aperture, state)
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
    if size(pupil) != size(pupil_mask(tel))
        throw(DimensionMismatchError("pupil mask size does not match telescope resolution"))
    end
    pupil_mask(tel) .= pupil
    pupil_reflectivity(tel) .= pupil
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::Real)
    fill!(pupil_reflectivity(tel), eltype(pupil_reflectivity(tel))(reflectivity))
    pupil_reflectivity(tel) .*= pupil_mask(tel)
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::AbstractMatrix)
    if size(reflectivity) != size(pupil_reflectivity(tel))
        throw(DimensionMismatchError("pupil_reflectivity size does not match telescope resolution"))
    end
    pupil_reflectivity(tel) .= reflectivity
    pupil_reflectivity(tel) .*= pupil_mask(tel)
    return tel
end

function pupil_expected_photon_map(tel::Telescope, src::AbstractSource)
    scale = photon_irradiance(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    out = similar(pupil_reflectivity(tel))
    reflectivity = pupil_reflectivity(tel)
    @. out = scale * reflectivity
    return out
end

function apply_spiders!(tel::Telescope; thickness::Real, angles::AbstractVector, offset_x::Real=0.0, offset_y::Real=0.0)
    Base.require_one_based_indexing(pupil_mask(tel))
    radius = tel.params.diameter / 2
    thickness_norm = thickness / radius
    offset_x_norm = offset_x / radius
    offset_y_norm = offset_y / radius
    _apply_spiders!(pupil_mask(tel), angles, thickness_norm, offset_x_norm,
        offset_y_norm)
    pupil_reflectivity(tel) .*= pupil_mask(tel)
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
