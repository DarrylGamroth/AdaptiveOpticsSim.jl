"""Configuration-only declaration of a telescope model."""
abstract type AbstractTelescopeDefinition end

"""
    TelescopeDefinition(; resolution, diameter, revision,
        central_obstruction=0, fov_arcsec=0, pupil_reflectivity=1,
        T=Float64)

Immutable, configuration-only declaration of the currently supported annular
telescope aperture. All fields are scalar values. The definition owns no
sampled aperture, mutable state, array backend, or compute device.

`pupil_reflectivity` is an intensity-throughput fraction. `revision` is the
explicit aperture-definition revision copied into the prepared telescope.
"""
struct TelescopeDefinition{T<:AbstractFloat} <: AbstractTelescopeDefinition
    resolution::Int
    diameter::T
    central_obstruction::T
    fov_arcsec::T
    pupil_reflectivity::T
    revision::UInt

    function TelescopeDefinition{T}(
        resolution::Int,
        diameter::T,
        central_obstruction::T,
        fov_arcsec::T,
        pupil_reflectivity::T,
        revision::UInt,
    ) where {T<:AbstractFloat}
        _require_valid_telescope_parameters(
            resolution,
            diameter,
            central_obstruction,
            fov_arcsec,
        )
        _require_valid_pupil_reflectivity(pupil_reflectivity)
        return new{T}(
            resolution,
            diameter,
            central_obstruction,
            fov_arcsec,
            pupil_reflectivity,
            revision,
        )
    end
end

function _checked_telescope_definition_resolution(::Bool)
    throw(InvalidConfiguration(
        "resolution must be an integer, not Bool"))
end

function _checked_telescope_definition_resolution(resolution::Integer)
    typemin(Int) <= resolution <= typemax(Int) || throw(
        InvalidConfiguration("resolution must be representable as Int"))
    return Int(resolution)
end

function _checked_telescope_definition_revision(::Bool)
    throw(InvalidConfiguration(
        "revision must be an integer, not Bool"))
end

function _checked_telescope_definition_revision(revision::Integer)
    zero(revision) <= revision <= typemax(UInt) || throw(
        InvalidConfiguration("revision must be representable as UInt"))
    return UInt(revision)
end

function TelescopeDefinition(;
    resolution::Integer,
    diameter::Real,
    revision::Integer,
    central_obstruction::Real=0.0,
    fov_arcsec::Real=0.0,
    pupil_reflectivity::Real=1.0,
    T::Type{<:AbstractFloat}=Float64,
)
    resolution_value = _checked_telescope_definition_resolution(resolution)
    revision_value = _checked_telescope_definition_revision(revision)
    diameter_value = _converted_positive_finite(
        diameter, T, "telescope diameter")
    central_obstruction_value = _converted_nonnegative_finite(
        central_obstruction, T, "telescope central obstruction")
    central_obstruction_value < one(T) || throw(InvalidConfiguration(
        "telescope central obstruction must lie in [0, 1)"))
    fov_value = _converted_nonnegative_finite(
        fov_arcsec, T, "telescope field of view")
    reflectivity_value = _converted_nonnegative_finite(
        pupil_reflectivity, T, "telescope pupil reflectivity")
    reflectivity_value <= one(T) || throw(InvalidConfiguration(
        "telescope pupil reflectivity must lie in [0, 1]"))
    return TelescopeDefinition{T}(
        resolution_value,
        diameter_value,
        central_obstruction_value,
        fov_value,
        reflectivity_value,
        revision_value,
    )
end

struct TelescopeParams{T<:AbstractFloat}
    resolution::Int
    diameter::T
    central_obstruction::T
    fov_arcsec::T
end

mutable struct TelescopeAperture{T<:AbstractFloat,
    Amask<:AbstractMatrix{Bool},
    Aref<:AbstractMatrix{T}}
    pupil::Amask
    reflectivity::Aref
    sampling_m::NTuple{2,T}
    origin_m::NTuple{2,T}
    revision::UInt
end

struct Telescope{P<:TelescopeParams,A<:TelescopeAperture,
    B<:AbstractArrayBackend} <: AbstractTelescope
    params::P
    aperture::A
end

@inline backend(::Telescope{<:Any,<:Any,B}) where {B} = B()
@inline pupil_mask(tel::Telescope) = tel.aperture.pupil
@inline pupil_reflectivity(tel::Telescope) = tel.aperture.reflectivity
@inline aperture_revision(tel::Telescope) = tel.aperture.revision

@inline function advance_aperture_revision!(tel::Telescope)
    tel.aperture.revision += one(UInt)
    return tel.aperture.revision
end

function Telescope(; resolution::Int,
    diameter::Real,
    central_obstruction::Real=0.0,
    fov_arcsec::Real=0.0,
    pupil_reflectivity::Union{Real,AbstractMatrix}=1.0,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())

    diameter_value = _converted_positive_finite(
        diameter, T, "telescope diameter")
    central_obstruction_value = _converted_nonnegative_finite(
        central_obstruction, T, "telescope central obstruction")
    central_obstruction_value < one(T) || throw(InvalidConfiguration(
        "telescope central obstruction must lie in [0, 1)"))
    fov_value = _converted_nonnegative_finite(
        fov_arcsec, T, "telescope field of view")
    params = TelescopeParams{T}(
        resolution,
        diameter_value,
        central_obstruction_value,
        fov_value,
    )
    _require_valid_telescope_parameters(
        params.resolution,
        params.diameter,
        params.central_obstruction,
        params.fov_arcsec,
    )
    _require_valid_pupil_reflectivity(pupil_reflectivity)
    selector = _resolve_backend_selector(backend)
    return _materialize_telescope(
        params,
        pupil_reflectivity,
        zero(UInt),
        selector,
    )
end

function _prepare_telescope(
    definition::TelescopeDefinition{T},
    target::AbstractComputeDevice,
) where {T<:AbstractFloat}
    _require_valid_telescope_definition(definition)
    selector = compute_device_backend(target)
    return _with_compute_device(target) do
        params = TelescopeParams{T}(
            definition.resolution,
            definition.diameter,
            definition.central_obstruction,
            definition.fov_arcsec,
        )
        telescope = _materialize_telescope(
            params,
            definition.pupil_reflectivity,
            definition.revision,
            selector,
        )
        validate_telescope_target(telescope, target)
        return telescope
    end
end

function _prepare_telescope(
    definition::AbstractTelescopeDefinition,
    ::AbstractComputeDevice,
)
    throw(InvalidConfiguration(
        "telescope definition $(typeof(definition)) does not implement " *
        "exact-target preparation"))
end

"""
    prepare_telescope(definition, target)

Validate and materialize a cold telescope definition on one exact compute
device. Target selection encloses all allocation and pupil construction, and
the caller's prior accelerator-device context is restored by the owning
backend extension.
"""
function prepare_telescope(
    definition::D,
    target::AbstractComputeDevice,
) where {D<:AbstractTelescopeDefinition}
    ismutabletype(D) && throw(InvalidConfiguration(
        "telescope definitions must be immutable"))
    return _prepare_telescope(definition, target)
end

@inline function _require_valid_telescope_definition(
    definition::TelescopeDefinition,
)
    _require_valid_telescope_parameters(
        definition.resolution,
        definition.diameter,
        definition.central_obstruction,
        definition.fov_arcsec,
    )
    _require_valid_pupil_reflectivity(definition.pupil_reflectivity)
    return nothing
end

@inline function _require_valid_telescope_parameters(
    resolution::Int,
    diameter::Real,
    central_obstruction::Real,
    fov_arcsec::Real,
)
    resolution > 0 || throw(InvalidConfiguration(
        "resolution must be positive"))
    isfinite(diameter) && diameter > zero(diameter) || throw(
        InvalidConfiguration("diameter must be finite and positive"))
    isfinite(central_obstruction) &&
        zero(central_obstruction) <= central_obstruction <
            one(central_obstruction) || throw(InvalidConfiguration(
        "central_obstruction must be finite and lie in [0, 1)"))
    isfinite(fov_arcsec) && fov_arcsec >= zero(fov_arcsec) || throw(
        InvalidConfiguration("fov_arcsec must be finite and nonnegative"))
    return nothing
end

function _materialize_telescope(
    params::TelescopeParams{T},
    pupil_reflectivity::Union{Real,AbstractMatrix},
    revision::UInt,
    selector::AbstractArrayBackend,
) where {T<:AbstractFloat}
    resolution = params.resolution
    pupil = allocate_array(selector, Bool, resolution, resolution)
    generate_pupil!(pupil, params)
    reflectivity = initialize_reflectivity(
        pupil,
        pupil_reflectivity,
        T,
        selector,
    )
    sampling_m = params.diameter / T(resolution)
    origin = -T(resolution - 1) * sampling_m / T(2)
    aperture = TelescopeAperture{T,typeof(pupil),typeof(reflectivity)}(
        pupil,
        reflectivity,
        (sampling_m, sampling_m),
        (origin, origin),
        revision,
    )
    return Telescope{typeof(params),typeof(aperture),typeof(selector)}(
        params, aperture)
end

"""Qualified fail-closed exact-target validation seam for prepared telescopes."""
@inline function validate_telescope_target(
    telescope::Telescope,
    target::AbstractComputeDevice,
)
    pupil_device = compute_device(pupil_mask(telescope))
    pupil_device == target || _throw_compute_device_error(
        :prepare_telescope,
        :wrong_device,
        target,
        "prepared telescope pupil occupies $(pupil_device)",
    )
    reflectivity_device = compute_device(pupil_reflectivity(telescope))
    reflectivity_device == target || _throw_compute_device_error(
        :prepare_telescope,
        :wrong_device,
        target,
        "prepared telescope reflectivity occupies $(reflectivity_device)",
    )
    return telescope
end

function validate_telescope_target(
    telescope::AbstractTelescope,
    ::AbstractComputeDevice,
)
    throw(InvalidConfiguration(
        "prepared telescope $(typeof(telescope)) does not implement " *
        "validate_telescope_target"))
end

function validate_telescope_target(telescope, ::AbstractComputeDevice)
    throw(InvalidConfiguration(
        "telescope preparation must return AbstractTelescope; got " *
        "$(typeof(telescope))"))
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

@inline function _require_valid_pupil_reflectivity(value::Real)
    isfinite(value) && zero(value) <= value <= one(value) || throw(
        InvalidConfiguration(
            "pupil_reflectivity must be finite and lie in [0, 1]"))
    return nothing
end

function _require_valid_pupil_reflectivity(values::AbstractMatrix)
    all(value -> isfinite(value) && zero(value) <= value <= one(value),
        values) || throw(InvalidConfiguration(
        "pupil_reflectivity values must be finite and lie in [0, 1]"))
    return nothing
end

function initialize_reflectivity(pupil::AbstractMatrix{Bool}, reflectivity::Real,
    ::Type{T}, selector::AbstractArrayBackend) where {T<:AbstractFloat}
    _require_valid_pupil_reflectivity(reflectivity)
    out = allocate_array(selector, T, size(pupil)...)
    fill!(out, T(reflectivity))
    out .*= pupil
    return out
end

function initialize_reflectivity(pupil::AbstractMatrix{Bool}, reflectivity::AbstractMatrix,
    ::Type{T}, selector::AbstractArrayBackend) where {T<:AbstractFloat}
    if size(reflectivity) != size(pupil)
        throw(DimensionMismatchError("pupil_reflectivity size does not match telescope resolution"))
    end
    _require_valid_pupil_reflectivity(reflectivity)
    out = allocate_array(selector, T, size(pupil)...)
    copyto!(out, T.(reflectivity))
    out .*= pupil
    return out
end

function set_pupil!(tel::Telescope, pupil::AbstractMatrix{Bool})
    if size(pupil) != size(pupil_mask(tel))
        throw(DimensionMismatchError("pupil mask size does not match telescope resolution"))
    end
    pupil_mask(tel) .= pupil
    pupil_reflectivity(tel) .= pupil
    advance_aperture_revision!(tel)
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::Real)
    _require_valid_pupil_reflectivity(reflectivity)
    fill!(pupil_reflectivity(tel), eltype(pupil_reflectivity(tel))(reflectivity))
    pupil_reflectivity(tel) .*= pupil_mask(tel)
    advance_aperture_revision!(tel)
    return tel
end

function set_pupil_reflectivity!(tel::Telescope, reflectivity::AbstractMatrix)
    if size(reflectivity) != size(pupil_reflectivity(tel))
        throw(DimensionMismatchError("pupil_reflectivity size does not match telescope resolution"))
    end
    _require_valid_pupil_reflectivity(reflectivity)
    pupil_reflectivity(tel) .= reflectivity
    pupil_reflectivity(tel) .*= pupil_mask(tel)
    advance_aperture_revision!(tel)
    return tel
end

function pupil_photon_rate_map(tel::Telescope, src::AbstractSource)
    irradiance = _require_physical_photon_irradiance(src,
        "pupil_photon_rate_map")
    scale = irradiance *
        (tel.params.diameter / tel.params.resolution)^2
    out = similar(pupil_reflectivity(tel))
    reflectivity = pupil_reflectivity(tel)
    @. out = scale * reflectivity
    return out
end

function apply_spiders!(tel::Telescope; thickness::Real,
    angles_deg::AbstractVector, offset_x::Real=0.0, offset_y::Real=0.0)
    Base.require_one_based_indexing(pupil_mask(tel))
    radius = tel.params.diameter / 2
    thickness_norm = thickness / radius
    offset_x_norm = offset_x / radius
    offset_y_norm = offset_y / radius
    _apply_spiders!(pupil_mask(tel), angles_deg, thickness_norm, offset_x_norm,
        offset_y_norm)
    pupil_reflectivity(tel) .*= pupil_mask(tel)
    advance_aperture_revision!(tel)
    return tel
end

function _apply_spiders!(pupil::AbstractMatrix{Bool}, angles_deg::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real)
    T = promote_type(typeof(thickness_norm), typeof(offset_x_norm), typeof(offset_y_norm))
    grid = default_mask_grid(pupil; T=T)
    for angle_deg in angles_deg
        apply_mask!(pupil, SpiderMask(thickness=thickness_norm, angle_rad=deg2rad(angle_deg), offset_x=offset_x_norm, offset_y=offset_y_norm,
            T=T); grid=grid)
    end
    return pupil
end

function _apply_spiders!(::ScalarCPUStyle, pupil::AbstractMatrix{Bool}, angles_deg::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real, cx::Real, cy::Real, scale::Real, n::Int)
    _apply_spiders!(pupil, angles_deg, thickness_norm, offset_x_norm, offset_y_norm)
    return pupil
end

function _apply_spiders!(::AcceleratorStyle, pupil::AbstractMatrix{Bool}, angles_deg::AbstractVector, thickness_norm::Real,
    offset_x_norm::Real, offset_y_norm::Real, cx::Real, cy::Real, scale::Real, n::Int)
    _apply_spiders!(pupil, angles_deg, thickness_norm, offset_x_norm, offset_y_norm)
    return pupil
end
