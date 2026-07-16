#
# Prepared focal-plane spatial filtering
#
# `SpatialFilter` owns only its fixed focal-plane mask definition.
# `SpatialFilterWorkspace` owns the single-writer FFT plans and scratch.
# The input `ElectricField` and output `PupilFunction` remain caller-owned.
#
abstract type SpatialFilterShape end
struct CircularFilter <: SpatialFilterShape end
struct SquareFilter <: SpatialFilterShape end
struct FoucaultFilter <: SpatialFilterShape end

struct SpatialFilterParams{T<:AbstractFloat}
    diameter::T
    zero_padding::Int
    pupil_resolution::Int
    padded_resolution::Int
end

struct SpatialFilter{
    S<:SpatialFilterShape,
    P<:SpatialFilterParams,
    C<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalElement
    params::P
    mask::C
    mask_shifted::C
end

@inline backend(::SpatialFilter{<:Any,<:Any,<:Any,B}) where {B} = B()

struct SpatialFilterWorkspace{
    C<:AbstractMatrix,
    Pf,
    Pi,
}
    fft_buffer::C
    filtered_field::C
    fft_plan::Pf
    ifft_plan::Pi
end

struct SpatialFilterPlan{
    T<:AbstractFloat,
    P<:SpatialFilterParams,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
    M<:AbstractMatrix{Bool},
}
    filter_params::P
    input_metadata::I
    output_metadata::O
    active_axes::NTuple{2,UnitRange{Int}}
    support::M
    wavelength_m::T
end

"""
    SpatialFilter(tel; shape=CircularFilter(), diameter=..., zero_padding=2)

Prepare an immutable focal-plane mask. Repeated filtering additionally uses a
`SpatialFilterWorkspace`, explicit input `ElectricField`, caller-owned output
`PupilFunction`, and a `SpatialFilterPlan` returned by
`prepare_spatial_filter`.
"""
function SpatialFilter(tel::Telescope;
    shape::SpatialFilterShape=CircularFilter(),
    diameter::Real=tel.params.resolution / 2,
    zero_padding::Int=2,
    T::Type{<:AbstractFloat}=eltype(opd_map(tel)),
    backend::AbstractArrayBackend=backend(tel))
    zero_padding >= 1 || throw(InvalidConfiguration(
        "spatial-filter zero_padding must be >= 1"))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    array_backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    n_pad = n * zero_padding
    params = SpatialFilterParams{T}(T(diameter), zero_padding, n, n_pad)
    mask = array_backend{Complex{T}}(undef, n_pad, n_pad)
    mask_shifted = similar(mask)
    spatial_filter = SpatialFilter{
        typeof(shape),typeof(params),typeof(mask),typeof(selector),
    }(params, mask, mask_shifted)
    set_spatial_filter!(spatial_filter, shape)
    return spatial_filter
end

function SpatialFilterWorkspace(spatial_filter::SpatialFilter)
    fft_buffer = similar(spatial_filter.mask_shifted)
    filtered_field = similar(spatial_filter.mask_shifted)
    fft_plan = plan_fft_backend!(fft_buffer)
    ifft_plan = plan_ifft_backend!(filtered_field)
    return SpatialFilterWorkspace(fft_buffer, filtered_field, fft_plan,
        ifft_plan)
end

set_spatial_filter!(spatial_filter::SpatialFilter{S}) where {S} =
    set_spatial_filter!(spatial_filter, S())

function set_spatial_filter!(spatial_filter::SpatialFilter,
    ::CircularFilter)
    diameter_padded = spatial_filter.params.diameter *
        spatial_filter.params.zero_padding
    T = real(eltype(spatial_filter.mask))
    inside = (one(eltype(spatial_filter.mask)) +
        im * one(eltype(spatial_filter.mask))) / sqrt(T(2))
    build_mask!(spatial_filter.mask,
        CircularAperture(radius=diameter_padded, T=T);
        grid=pixel_mask_grid(spatial_filter.mask; T=T), inside=inside)
    return finalize_spatial_filter_mask!(spatial_filter)
end

function set_spatial_filter!(spatial_filter::SpatialFilter, ::SquareFilter)
    n = spatial_filter.params.padded_resolution
    diameter_padded = spatial_filter.params.diameter *
        spatial_filter.params.zero_padding
    half = Int(round(diameter_padded / 2))
    center = Int(round((n + 1) / 2))
    range = max(1, center - half):min(n, center + half)
    T = real(eltype(spatial_filter.mask))
    inside = (one(eltype(spatial_filter.mask)) +
        im * one(eltype(spatial_filter.mask))) / sqrt(T(2))
    build_mask!(spatial_filter.mask, RectangularROI(range, range);
        inside=inside)
    return finalize_spatial_filter_mask!(spatial_filter)
end

function set_spatial_filter!(spatial_filter::SpatialFilter,
    ::FoucaultFilter)
    n = spatial_filter.params.padded_resolution
    T = real(eltype(spatial_filter.mask))
    inside = (one(eltype(spatial_filter.mask)) +
        im * one(eltype(spatial_filter.mask))) / sqrt(T(2))
    build_mask!(spatial_filter.mask,
        RectangularROI(1:floor(Int, n / 2), 1:n); inside=inside)
    return finalize_spatial_filter_mask!(spatial_filter)
end

function finalize_spatial_filter_mask!(spatial_filter::SpatialFilter)
    temporary = similar(spatial_filter.mask)
    fill!(temporary, zero(eltype(temporary)))
    @views temporary[1:end-1, 1:end-1] .=
        spatial_filter.mask[2:end, 2:end]
    copyto!(spatial_filter.mask, temporary)
    fftshift2d!(spatial_filter.mask_shifted, spatial_filter.mask)
    return spatial_filter
end

function prepare_spatial_filter(tel::Telescope,
    spatial_filter::SpatialFilter, input::ElectricField,
    output::PupilFunction)
    require_centered_plane_geometry(input.metadata;
        label="spatial-filter input ElectricField")
    require_centered_plane_geometry(output.metadata;
        label="spatial-filter output PupilFunction")
    require_metric_coordinates(input.metadata;
        label="spatial-filter input ElectricField")
    require_metric_coordinates(output.metadata;
        label="spatial-filter output PupilFunction")
    input.metadata.dimensions == (
        spatial_filter.params.padded_resolution,
        spatial_filter.params.padded_resolution,
    ) || throw(DimensionMismatchError(
        "ElectricField dimensions must match the spatial-filter mask"))
    output.metadata.dimensions == (
        spatial_filter.params.pupil_resolution,
        spatial_filter.params.pupil_resolution,
    ) || throw(DimensionMismatchError(
        "PupilFunction dimensions must match the spatial-filter pupil grid"))
    output.metadata.dimensions == size(pupil_mask(tel)) ||
        throw(DimensionMismatchError(
            "spatial-filter output must match the telescope aperture"))
    typeof(input.metadata.kind) === PupilPlane || throw(InvalidConfiguration(
        "spatial-filter input must be a pupil-plane ElectricField"))
    typeof(output.metadata.kind) === PupilPlane || throw(InvalidConfiguration(
        "spatial-filter output must be a pupil-plane PupilFunction"))
    input.metadata.sampling == output.metadata.sampling ||
        throw(InvalidConfiguration(
            "spatial-filter input and output sampling must match"))
    input.metadata.orientation == output.metadata.orientation ||
        throw(InvalidConfiguration(
            "spatial-filter input and output orientation must match"))
    selector = require_same_backend(tel, spatial_filter, input, output)
    input.metadata.device == output.metadata.device ||
        throw(InvalidConfiguration(
            "spatial-filter input and output must occupy the same device"))
    plane_device(spatial_filter.mask_shifted) == input.metadata.device ||
        throw(InvalidConfiguration(
            "spatial-filter mask and input must occupy the same device"))
    typeof(selector) === typeof(input.metadata.backend) ||
        throw(InvalidConfiguration(
            "spatial-filter metadata backend does not match storage"))
    ox, oy = field_embedding_offsets(
        spatial_filter.params.pupil_resolution,
        spatial_filter.params.padded_resolution)
    n = spatial_filter.params.pupil_resolution
    axes = (ox + 1:ox + n, oy + 1:oy + n)
    T = eltype(output.opd)
    return SpatialFilterPlan{
        T,typeof(spatial_filter.params),typeof(input.metadata),
        typeof(output.metadata),
        typeof(pupil_mask(tel)),
    }(
        spatial_filter.params,
        input.metadata,
        output.metadata,
        axes,
        pupil_mask(tel),
        T(electric_field_wavelength(input)),
    )
end

function filter!(output::PupilFunction, input::ElectricField,
    spatial_filter::SpatialFilter, plan::SpatialFilterPlan,
    workspace::SpatialFilterWorkspace)
    input.metadata == plan.input_metadata || throw(InvalidConfiguration(
        "ElectricField metadata does not match the prepared spatial-filter plan"))
    output.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "PupilFunction metadata does not match the prepared spatial-filter plan"))
    spatial_filter.params == plan.filter_params ||
        throw(InvalidConfiguration(
            "SpatialFilter does not match the prepared spatial-filter plan"))
    size(spatial_filter.mask_shifted) == input.metadata.dimensions ||
        throw(DimensionMismatchError(
            "spatial-filter mask does not match its prepared field"))
    size(workspace.fft_buffer) == input.metadata.dimensions ||
        throw(DimensionMismatchError(
            "spatial-filter workspace does not match its prepared field"))

    copyto!(workspace.fft_buffer, input.values)
    execute_fft_plan!(workspace.fft_buffer, workspace.fft_plan)
    @. workspace.filtered_field = workspace.fft_buffer *
        spatial_filter.mask_shifted
    execute_fft_plan!(workspace.filtered_field, workspace.ifft_plan)
    opd_per_radian = plan.wavelength_m / (2 * eltype(output.opd)(pi))
    @views begin
        region = workspace.filtered_field[plan.active_axes...]
        @. output.opd = angle(region) * opd_per_radian * plan.support
        @. output.amplitude = abs(region)
    end
    return output
end
