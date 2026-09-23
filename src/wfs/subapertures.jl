abstract type AbstractValidSubaperturePolicy end

struct GeometryValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    threshold::T
end

struct RelativeIlluminationValidSubapertures{T<:AbstractFloat} <: AbstractValidSubaperturePolicy
    peak_fraction::T
end

function GeometryValidSubapertures(; threshold::Real=0.1,
    T::Type{<:AbstractFloat}=typeof(float(threshold)))
    value = T(threshold)
    isfinite(value) && zero(T) <= value <= one(T) ||
        throw(InvalidConfiguration(
            "GeometryValidSubapertures threshold must lie in [0, 1]"))
    return GeometryValidSubapertures{T}(value)
end

function RelativeIlluminationValidSubapertures(; peak_fraction::Real=0.5,
    T::Type{<:AbstractFloat}=typeof(float(peak_fraction)))
    value = T(peak_fraction)
    isfinite(value) && zero(T) <= value <= one(T) ||
        throw(InvalidConfiguration(
            "RelativeIlluminationValidSubapertures peak_fraction must lie in [0, 1]"))
    return RelativeIlluminationValidSubapertures{T}(value)
end

mutable struct SubapertureLayoutState
    revision::UInt
end

"""
Fixed subaperture geometry with caller-owned mutable mask storage.

The execution mask is authoritative at construction. Its host mirror and
column-major `CartesianIndex` list are synchronized before the layout is
returned. Maintained updates synchronize all three representations and advance
a cold-configuration revision; direct storage mutation is unsupported.
"""
struct SubapertureLayout{T<:AbstractFloat,A<:AbstractMatrix{Bool},M<:Matrix{Bool},V<:Vector{CartesianIndex{2}}}
    n_subap::Int
    subap_pixels::Int
    pitch_m::T
    threshold::T
    valid_mask::A
    valid_mask_host::M
    valid_indices_host::V
    state::SubapertureLayoutState
end

function SubapertureLayout(n_subap::Int, pupil_resolution::Int, diameter::Real, threshold::Real,
    valid_mask::AbstractMatrix{Bool}, valid_mask_host::Matrix{Bool})
    n_subap > 0 || throw(InvalidConfiguration(
        "SubapertureLayout n_subap must be positive"))
    pupil_resolution > 0 || throw(InvalidConfiguration(
        "SubapertureLayout pupil_resolution must be positive"))
    pupil_resolution % n_subap == 0 || throw(InvalidConfiguration(
        "SubapertureLayout pupil_resolution must be divisible by n_subap"))
    size(valid_mask) == (n_subap, n_subap) || throw(DimensionMismatchError(
        "SubapertureLayout valid_mask must have size (n_subap, n_subap)"))
    size(valid_mask_host) == size(valid_mask) || throw(DimensionMismatchError(
        "SubapertureLayout host and execution masks must have the same size"))
    T = promote_type(typeof(float(diameter)), typeof(float(threshold)))
    diameter_m = T(diameter)
    threshold_fraction = T(threshold)
    isfinite(diameter_m) && diameter_m > zero(T) ||
        throw(InvalidConfiguration(
            "SubapertureLayout diameter must be finite and positive"))
    isfinite(threshold_fraction) &&
        zero(T) <= threshold_fraction <= one(T) ||
        throw(InvalidConfiguration(
            "SubapertureLayout threshold must lie in [0, 1]"))
    subap_pixels = div(pupil_resolution, n_subap)
    pitch_m = diameter_m / T(n_subap)
    valid_indices_host = CartesianIndex{2}[]
    sizehint!(valid_indices_host, length(valid_mask_host))
    layout = SubapertureLayout{T, typeof(valid_mask), typeof(valid_mask_host), Vector{CartesianIndex{2}}}(
        n_subap,
        subap_pixels,
        pitch_m,
        threshold_fraction,
        valid_mask,
        valid_mask_host,
        valid_indices_host,
        SubapertureLayoutState(UInt(0)),
    )
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    return layout
end

@inline subaperture_layout_revision(layout::SubapertureLayout) =
    layout.state.revision

@inline function _advance_subaperture_layout_revision!(
    layout::SubapertureLayout)
    layout.state.revision += UInt(1)
    return layout.state.revision
end

@inline function _copy_valid_mask_to_host!(host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    return _copy_valid_mask_to_host!(execution_style(mask), host, mask)
end

@inline function _copy_valid_mask_to_host!(::ScalarCPUStyle, host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    copyto!(host, mask)
    return host
end

@inline function _copy_valid_mask_to_host!(::ExecutionStyle, host::Matrix{Bool}, mask::AbstractMatrix{Bool})
    copyto!(host, Array(mask))
    return host
end

@inline function _refresh_valid_indices_host!(layout::SubapertureLayout)
    indices = layout.valid_indices_host
    resize!(indices, 0)
    @inbounds for I in CartesianIndices(layout.valid_mask_host)
        layout.valid_mask_host[I] && push!(indices, I)
    end
    return indices
end

"""
    set_valid_subapertures!(layout, valid_subapertures)

Snapshot an explicit one-based square lenslet mask into a subaperture layout.
The execution mask, host mirror, and column-major valid-index list are updated
together before the cold-configuration revision advances. Prepared optics and
other owners that bind the prior revision must be prepared again.
"""
function set_valid_subapertures!(
    layout::SubapertureLayout,
    valid_subapertures::AbstractMatrix{Bool},
)
    Base.require_one_based_indexing(valid_subapertures)
    size(valid_subapertures) == (layout.n_subap, layout.n_subap) ||
        throw(DimensionMismatchError(
            "valid-subaperture mask must match the square lenslet layout",
        ))
    copyto!(layout.valid_mask, valid_subapertures)
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool})
    build_mask!(layout.valid_mask, SubapertureGridMask(threshold=layout.threshold, T=typeof(layout.threshold)), pupil)
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout, pupil::AbstractMatrix{Bool},
    policy::GeometryValidSubapertures)
    threshold = convert(typeof(layout.threshold), policy.threshold)
    build_mask!(layout.valid_mask,
        SubapertureGridMask(threshold=threshold,
            T=typeof(layout.threshold)), pupil)
    _copy_valid_mask_to_host!(layout.valid_mask_host, layout.valid_mask)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout!(layout::SubapertureLayout,
    support_map::AbstractMatrix{T},
    policy::GeometryValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(support_map) == (n_sub * sub, n_sub * sub) ||
        throw(DimensionMismatchError(
            "support map size must match subaperture layout"))
    support_host = _host_support_map(execution_style(support_map),
        support_map)
    threshold = convert(Float64, policy.threshold)
    denominator = sub * sub
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        illuminated = count(!iszero,
            @view support_host[xs:xe, ys:ye])
        layout.valid_mask_host[i, j] =
            illuminated / denominator >= threshold
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

@inline _host_support_map(::ScalarCPUStyle, support_map::AbstractMatrix) = support_map
@inline _host_support_map(::ExecutionStyle, support_map::AbstractMatrix) = Array(support_map)

function update_subaperture_layout!(layout::SubapertureLayout, support_map::AbstractMatrix{T},
    policy::RelativeIlluminationValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(support_map, 1) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    size(support_map, 2) == n_sub * sub || throw(DimensionMismatchError("support map size must match subaperture layout"))
    support_host = _host_support_map(execution_style(support_map), support_map)
    peak = zero(eltype(support_host))
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        total = sum(@view support_host[xs:xe, ys:ye])
        peak = max(peak, total)
    end
    cutoff = convert(eltype(support_host), policy.peak_fraction) * peak
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        layout.valid_mask_host[i, j] = sum(@view support_host[xs:xe, ys:ye]) >= cutoff
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

function update_subaperture_layout_from_amplitude!(
    layout::SubapertureLayout, amplitude::AbstractMatrix{T},
    policy::RelativeIlluminationValidSubapertures) where {T<:Real}
    n_sub = layout.n_subap
    sub = layout.subap_pixels
    size(amplitude) == (n_sub * sub, n_sub * sub) ||
        throw(DimensionMismatchError(
            "pupil amplitude size must match subaperture layout"))
    amplitude_host = _host_support_map(execution_style(amplitude),
        amplitude)
    peak = zero(eltype(amplitude_host))
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        total = sum(abs2, @view amplitude_host[xs:xe, ys:ye])
        peak = max(peak, total)
    end
    cutoff = convert(eltype(amplitude_host), policy.peak_fraction) * peak
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = i * sub
        ye = j * sub
        layout.valid_mask_host[i, j] =
            sum(abs2, @view amplitude_host[xs:xe, ys:ye]) >= cutoff
    end
    copyto!(layout.valid_mask, layout.valid_mask_host)
    _refresh_valid_indices_host!(layout)
    _advance_subaperture_layout_revision!(layout)
    return layout
end

@inline valid_subaperture_indices(layout::SubapertureLayout) = layout.valid_indices_host
@inline n_valid_subapertures(layout::SubapertureLayout) = length(layout.valid_indices_host)

function subaperture_layout end
