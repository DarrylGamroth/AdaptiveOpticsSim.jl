#
# Prepared focal-plane spatial filtering
#
# `SpatialFilter` owns only its fixed focal-plane mask definition.
# `SpatialFilterPlan` owns the reusable numerical and metadata contract.
# `SpatialFilterWorkspace` owns the single-writer FFT plans and scratch.
# `PreparedSpatialFilter` binds exact caller-visible input and output products.
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

"""Single-writer FFT resources and scratch for spatial filtering."""
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

"""Run-immutable spatial-filter numerical and optical-plane contract."""
struct SpatialFilterPlan{
    T<:AbstractFloat,
    P<:SpatialFilterParams,
    I<:OpticalPlaneMetadata,
    O<:OpticalPlaneMetadata,
}
    filter_params::P
    input_metadata::I
    output_metadata::O
    active_axes::NTuple{2,UnitRange{Int}}
    wavelength_m::T
    aperture_revision::UInt
end

struct _PreparedSpatialFilterToken end
const _PREPARED_SPATIAL_FILTER_TOKEN = _PreparedSpatialFilterToken()

"""
    PreparedSpatialFilter

Exact prepared owner for one focal-plane spatial-filter stage. The filter,
input field, output pupil function, and single-writer workspace are bound for
the lifetime of this owner. Call `filter!(prepared)` for repeated execution.
"""
struct PreparedSpatialFilter{
    F<:SpatialFilter,
    I<:ElectricField,
    O<:PupilFunction,
    P<:SpatialFilterPlan,
    W<:SpatialFilterWorkspace,
}
    spatial_filter::F
    input::I
    output::O
    plan::P
    workspace::W

    function PreparedSpatialFilter(
        ::_PreparedSpatialFilterToken,
        spatial_filter::F,
        input::I,
        output::O,
        plan::P,
        workspace::W,
    ) where {
        F<:SpatialFilter,
        I<:ElectricField,
        O<:PupilFunction,
        P<:SpatialFilterPlan,
        W<:SpatialFilterWorkspace,
    }
        return new{F,I,O,P,W}(
            spatial_filter, input, output, plan, workspace)
    end
end

"""Return the reusable plan from a prepared spatial-filter owner."""
@inline spatial_filter_plan(prepared::PreparedSpatialFilter) = prepared.plan

"""Return the single-writer workspace from a prepared spatial-filter owner."""
@inline spatial_filter_workspace(prepared::PreparedSpatialFilter) =
    prepared.workspace

"""Return the caller-visible output from a prepared spatial-filter owner."""
@inline spatial_filter_output(prepared::PreparedSpatialFilter) =
    prepared.output

"""
    SpatialFilter(tel; shape=CircularFilter(), diameter=..., zero_padding=2)

Prepare an immutable focal-plane mask. Repeated filtering additionally uses a
single-writer workspace, an explicit input `ElectricField`, and a caller-owned
output `PupilFunction` bound by `prepare_spatial_filter`.
"""
function SpatialFilter(tel::Telescope;
    shape::SpatialFilterShape=CircularFilter(),
    diameter::Real=tel.params.resolution / 2,
    zero_padding::Int=2,
    T::Type{<:AbstractFloat}=eltype(pupil_reflectivity(tel)),
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

function _require_spatial_filter_plan_workspace(
    plan::SpatialFilterPlan,
    workspace::SpatialFilterWorkspace,
)
    for (storage, label) in (
        (workspace.fft_buffer, "spatial-filter FFT scratch"),
        (workspace.filtered_field, "spatial-filter filtered-field scratch"),
    )
        size(storage) == plan.input_metadata.dimensions || throw(
            DimensionMismatchError(
                "$label dimensions do not match the prepared plan"))
        eltype(storage) === plan.input_metadata.numeric_type || throw(
            InvalidConfiguration(
                "$label numeric type does not match the prepared plan"))
        typeof(backend(storage)) === typeof(plan.input_metadata.backend) ||
            throw(InvalidConfiguration(
                "$label backend does not match the prepared plan"))
        compute_device(storage) == plan.input_metadata.device || throw(
            InvalidConfiguration(
                "$label device does not match the prepared plan"))
    end
    Base.mightalias(workspace.fft_buffer, workspace.filtered_field) && throw(
        InvalidConfiguration(
            "spatial-filter workspace arrays must not alias"))
    return workspace
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
    aperture_revision(output) == aperture_revision(tel) || throw(
        InvalidConfiguration(
            "spatial-filter output aperture revision does not match the telescope"))
    Array(output.support) == Array(pupil_mask(tel)) || throw(InvalidConfiguration(
        "spatial-filter output support does not match the telescope aperture"))
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
    compute_device(spatial_filter.mask_shifted) == input.metadata.device ||
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
    plan = SpatialFilterPlan{
        T,typeof(spatial_filter.params),typeof(input.metadata),
        typeof(output.metadata),
    }(
        spatial_filter.params,
        input.metadata,
        output.metadata,
        axes,
        T(electric_field_wavelength(input)),
        aperture_revision(tel),
    )
    workspace = SpatialFilterWorkspace(spatial_filter)
    _require_spatial_filter_plan_workspace(plan, workspace)
    return PreparedSpatialFilter(
        _PREPARED_SPATIAL_FILTER_TOKEN,
        spatial_filter,
        input,
        output,
        plan,
        workspace,
    )
end

function _require_prepared_spatial_filter_bindings(
    prepared::PreparedSpatialFilter,
    spatial_filter::SpatialFilter,
    input::ElectricField,
    output::PupilFunction,
    workspace::SpatialFilterWorkspace,
)
    spatial_filter === prepared.spatial_filter || throw(
        InvalidConfiguration(
            "spatial filter does not match its prepared owner"))
    input === prepared.input || throw(InvalidConfiguration(
        "spatial-filter input does not match its prepared owner"))
    output === prepared.output || throw(InvalidConfiguration(
        "spatial-filter output does not match its prepared owner"))
    workspace === prepared.workspace || throw(InvalidConfiguration(
        "spatial-filter workspace does not match its prepared owner"))
    return nothing
end

function _filter_prepared_spatial_filter!(
    prepared::PreparedSpatialFilter,
    spatial_filter::SpatialFilter,
    input::ElectricField,
    output::PupilFunction,
    workspace::SpatialFilterWorkspace,
)
    _require_prepared_spatial_filter_bindings(
        prepared, spatial_filter, input, output, workspace)
    plan = prepared.plan
    input.metadata == plan.input_metadata || throw(InvalidConfiguration(
        "ElectricField metadata does not match the prepared spatial-filter plan"))
    output.metadata == plan.output_metadata || throw(InvalidConfiguration(
        "PupilFunction metadata does not match the prepared spatial-filter plan"))
    aperture_revision(output) == plan.aperture_revision || throw(
        InvalidConfiguration(
            "PupilFunction aperture revision does not match the prepared spatial-filter plan"))
    spatial_filter.params == plan.filter_params ||
        throw(InvalidConfiguration(
            "SpatialFilter does not match the prepared spatial-filter plan"))
    size(spatial_filter.mask_shifted) == input.metadata.dimensions ||
        throw(DimensionMismatchError(
            "spatial-filter mask does not match its prepared field"))
    _require_spatial_filter_plan_workspace(plan, workspace)

    copyto!(workspace.fft_buffer, input.values)
    execute_fft_plan!(workspace.fft_buffer, workspace.fft_plan)
    @. workspace.filtered_field = workspace.fft_buffer *
        spatial_filter.mask_shifted
    execute_fft_plan!(workspace.filtered_field, workspace.ifft_plan)
    opd_per_radian = plan.wavelength_m / (2 * eltype(output.opd)(pi))
    @views begin
        region = workspace.filtered_field[plan.active_axes...]
        @. output.opd = angle(region) * opd_per_radian * output.support
        @. output.amplitude = abs(region)
    end
    return output
end

"""Execute one prepared focal-plane spatial-filter stage."""
function filter!(prepared::PreparedSpatialFilter)
    return _filter_prepared_spatial_filter!(
        prepared,
        prepared.spatial_filter,
        prepared.input,
        prepared.output,
        prepared.workspace,
    )
end
