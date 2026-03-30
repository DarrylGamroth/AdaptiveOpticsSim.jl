abstract type AbstractMaskPrimitive end

struct MaskGrid{T<:AbstractFloat}
    center_i::T
    center_j::T
    scale_i::T
    scale_j::T
end

struct CircularAperture{T<:AbstractFloat} <: AbstractMaskPrimitive
    radius::T
    center_x::T
    center_y::T
end

struct AnnularAperture{T<:AbstractFloat} <: AbstractMaskPrimitive
    inner_radius::T
    outer_radius::T
    center_x::T
    center_y::T
end

struct SpiderMask{T<:AbstractFloat} <: AbstractMaskPrimitive
    thickness::T
    angle_rad::T
    offset_x::T
    offset_y::T
end

struct RectangularROI <: AbstractMaskPrimitive
    row_start::Int
    row_stop::Int
    col_start::Int
    col_stop::Int
end

struct SubapertureGridMask{T<:AbstractFloat} <: AbstractMaskPrimitive
    threshold::T
end

CircularAperture(; radius::Real=1.0, center_x::Real=0.0, center_y::Real=0.0, T::Type{<:AbstractFloat}=Float64) =
    CircularAperture{T}(T(radius), T(center_x), T(center_y))

AnnularAperture(; inner_radius::Real=0.0, outer_radius::Real=1.0, center_x::Real=0.0, center_y::Real=0.0, T::Type{<:AbstractFloat}=Float64) =
    AnnularAperture{T}(T(inner_radius), T(outer_radius), T(center_x), T(center_y))

SpiderMask(; thickness::Real, angle_rad::Real, offset_x::Real=0.0, offset_y::Real=0.0, T::Type{<:AbstractFloat}=Float64) =
    SpiderMask{T}(T(thickness), T(angle_rad), T(offset_x), T(offset_y))

RectangularROI(rows::UnitRange{Int}, cols::UnitRange{Int}) = RectangularROI(first(rows), last(rows), first(cols), last(cols))

SubapertureGridMask(; threshold::Real=0.1, T::Type{<:AbstractFloat}=Float64) =
    SubapertureGridMask{T}(T(threshold))

@inline _mask_inside_value(::Type{Bool}) = true
@inline _mask_inside_value(::Type{T}) where {T<:Number} = one(T)
@inline _mask_outside_value(::Type{Bool}) = false
@inline _mask_outside_value(::Type{T}) where {T<:Number} = zero(T)

@inline _mask_coord_type(::Type{Bool}) = Float64
@inline _mask_coord_type(::Type{T}) where {T<:AbstractFloat} = T
@inline _mask_coord_type(::Type{Complex{T}}) where {T<:AbstractFloat} = T
@inline _mask_coord_type(::Type{T}) where {T<:Real} = Float64
@inline _mask_execution_style(out::BitArray) = ScalarCPUStyle()
@inline _mask_execution_style(out::AbstractArray) = execution_style(out)

function default_mask_grid(out::AbstractMatrix; T::Type{<:AbstractFloat}=_mask_coord_type(eltype(out)))
    n_i, n_j = size(out)
    scale = T(min(n_i, n_j) / 2)
    return MaskGrid{T}(T((n_i + 1) / 2), T((n_j + 1) / 2), scale, scale)
end

function pixel_mask_grid(out::AbstractMatrix; T::Type{<:AbstractFloat}=_mask_coord_type(eltype(out)))
    n_i, n_j = size(out)
    return MaskGrid{T}(T((n_i + 1) / 2), T((n_j + 1) / 2), one(T), one(T))
end

@kernel function radial_mask_kernel!(out, inner2, outer2, center_x, center_y, grid_ci, grid_cj, scale_i, scale_j, inside, outside,
    n_i::Int, n_j::Int)
    i, j = @index(Global, NTuple)
    if i <= n_i && j <= n_j
        x = (i - grid_ci) / scale_i - center_x
        y = (j - grid_cj) / scale_j - center_y
        r2 = x^2 + y^2
        value = (r2 >= inner2) & (r2 <= outer2)
        @inbounds out[i, j] = value ? inside : outside
    end
end

@kernel function rectangular_roi_kernel!(out, row_start::Int, row_stop::Int, col_start::Int, col_stop::Int, inside, outside, n_i::Int, n_j::Int)
    i, j = @index(Global, NTuple)
    if i <= n_i && j <= n_j
        value = (row_start <= i <= row_stop) & (col_start <= j <= col_stop)
        @inbounds out[i, j] = value ? inside : outside
    end
end

@kernel function spider_mask_kernel!(out, thickness, a, b, offset_x, offset_y, grid_ci, grid_cj, scale_i, scale_j, masked, n_i::Int, n_j::Int)
    i, j = @index(Global, NTuple)
    if i <= n_i && j <= n_j
        x = (i - grid_ci) / scale_i - offset_x
        y = (j - grid_cj) / scale_j - offset_y
        dist = abs(a * x + b * y)
        if dist <= thickness
            @inbounds out[i, j] = masked
        end
    end
end

@kernel function subaperture_grid_mask_kernel!(valid_mask, pupil, threshold, sub_i::Int, sub_j::Int, n_sub_i::Int, n_sub_j::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub_i && j <= n_sub_j
        row_start = (i - 1) * sub_i + 1
        col_start = (j - 1) * sub_j + 1
        total = zero(promote_type(eltype(pupil), typeof(threshold)))
        @inbounds for di in 0:sub_i-1, dj in 0:sub_j-1
            total += pupil[row_start + di, col_start + dj]
        end
        @inbounds valid_mask[i, j] = total / (sub_i * sub_j) > threshold
    end
end

function build_mask!(out::AbstractMatrix, primitive::CircularAperture;
    grid::MaskGrid=default_mask_grid(out),
    inside=_mask_inside_value(eltype(out)),
    outside=_mask_outside_value(eltype(out)))
    _build_radial_mask!(_mask_execution_style(out), out, zero(primitive.radius), primitive.radius, primitive.center_x, primitive.center_y, grid, inside, outside)
    return out
end

function build_mask!(out::AbstractMatrix, primitive::AnnularAperture;
    grid::MaskGrid=default_mask_grid(out),
    inside=_mask_inside_value(eltype(out)),
    outside=_mask_outside_value(eltype(out)))
    _build_radial_mask!(_mask_execution_style(out), out, primitive.inner_radius, primitive.outer_radius, primitive.center_x, primitive.center_y, grid, inside, outside)
    return out
end

function build_mask!(out::AbstractMatrix, roi::RectangularROI;
    inside=_mask_inside_value(eltype(out)),
    outside=_mask_outside_value(eltype(out)))
    _build_rectangular_roi!(_mask_execution_style(out), out, roi, inside, outside)
    return out
end

function apply_mask!(out::AbstractMatrix, primitive::SpiderMask;
    grid::MaskGrid=default_mask_grid(out),
    masked=_mask_outside_value(eltype(out)))
    _apply_spider_mask!(_mask_execution_style(out), out, primitive, grid, masked)
    return out
end

function build_mask!(valid_mask::AbstractMatrix{Bool}, primitive::SubapertureGridMask, pupil::AbstractMatrix{Bool})
    Base.require_one_based_indexing(valid_mask, pupil)
    n_i, n_j = size(pupil)
    n_sub_i, n_sub_j = size(valid_mask)
    n_i == n_j || throw(DimensionMismatchError("pupil must be square"))
    n_sub_i == n_sub_j || throw(DimensionMismatchError("valid_mask must be square"))
    n_i % n_sub_i == 0 || throw(DimensionMismatchError("pupil size must be divisible by valid_mask size"))
    n_j % n_sub_j == 0 || throw(DimensionMismatchError("pupil size must be divisible by valid_mask size"))
    sub_i = div(n_i, n_sub_i)
    sub_j = div(n_j, n_sub_j)
    _build_subaperture_grid_mask!(_mask_execution_style(valid_mask), valid_mask, pupil, primitive.threshold, sub_i, sub_j, n_sub_i, n_sub_j)
    return valid_mask
end

function _build_radial_mask!(::ScalarCPUStyle, out::AbstractMatrix, inner_radius::Real, outer_radius::Real, center_x::Real, center_y::Real,
    grid::MaskGrid, inside, outside)
    inner2 = inner_radius^2
    outer2 = outer_radius^2
    n_i, n_j = size(out)
    @inbounds for i in 1:n_i, j in 1:n_j
        x = (i - grid.center_i) / grid.scale_i - center_x
        y = (j - grid.center_j) / grid.scale_j - center_y
        r2 = x^2 + y^2
        out[i, j] = (r2 >= inner2) & (r2 <= outer2) ? inside : outside
    end
    return out
end

function _build_radial_mask!(style::AcceleratorStyle, out::AbstractMatrix, inner_radius::Real, outer_radius::Real, center_x::Real, center_y::Real,
    grid::MaskGrid, inside, outside)
    n_i, n_j = size(out)
    launch_kernel!(style, radial_mask_kernel!, out, inner_radius^2, outer_radius^2, center_x, center_y,
        grid.center_i, grid.center_j, grid.scale_i, grid.scale_j, inside, outside, n_i, n_j; ndrange=size(out))
    return out
end

function _build_rectangular_roi!(::ScalarCPUStyle, out::AbstractMatrix, roi::RectangularROI, inside, outside)
    n_i, n_j = size(out)
    @inbounds for i in 1:n_i, j in 1:n_j
        out[i, j] = (roi.row_start <= i <= roi.row_stop) & (roi.col_start <= j <= roi.col_stop) ? inside : outside
    end
    return out
end

function _build_rectangular_roi!(style::AcceleratorStyle, out::AbstractMatrix, roi::RectangularROI, inside, outside)
    n_i, n_j = size(out)
    launch_kernel!(style, rectangular_roi_kernel!, out, roi.row_start, roi.row_stop, roi.col_start, roi.col_stop, inside, outside, n_i, n_j;
        ndrange=size(out))
    return out
end

function _apply_spider_mask!(::ScalarCPUStyle, out::AbstractMatrix, primitive::SpiderMask, grid::MaskGrid, masked)
    sθ, cθ = sincos(primitive.angle_rad)
    a = -sθ
    b = cθ
    n_i, n_j = size(out)
    @inbounds for i in 1:n_i, j in 1:n_j
        x = (i - grid.center_i) / grid.scale_i - primitive.offset_x
        y = (j - grid.center_j) / grid.scale_j - primitive.offset_y
        dist = abs(a * x + b * y)
        if dist <= primitive.thickness
            out[i, j] = masked
        end
    end
    return out
end

function _apply_spider_mask!(style::AcceleratorStyle, out::AbstractMatrix, primitive::SpiderMask, grid::MaskGrid, masked)
    sθ, cθ = sincos(primitive.angle_rad)
    a = -sθ
    b = cθ
    n_i, n_j = size(out)
    launch_kernel!(style, spider_mask_kernel!, out, primitive.thickness, a, b, primitive.offset_x, primitive.offset_y,
        grid.center_i, grid.center_j, grid.scale_i, grid.scale_j, masked, n_i, n_j; ndrange=size(out))
    return out
end

function _build_subaperture_grid_mask!(::ScalarCPUStyle, valid_mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool},
    threshold::Real, sub_i::Int, sub_j::Int, n_sub_i::Int, n_sub_j::Int)
    @inbounds for i in 1:n_sub_i, j in 1:n_sub_j
        row_start = (i - 1) * sub_i + 1
        col_start = (j - 1) * sub_j + 1
        row_stop = row_start + sub_i - 1
        col_stop = col_start + sub_j - 1
        valid_mask[i, j] = mean(@view pupil[row_start:row_stop, col_start:col_stop]) > threshold
    end
    return valid_mask
end

function _build_subaperture_grid_mask!(style::AcceleratorStyle, valid_mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool},
    threshold::Real, sub_i::Int, sub_j::Int, n_sub_i::Int, n_sub_j::Int)
    threshold_t = promote_type(eltype(pupil), typeof(threshold))(threshold)
    launch_kernel!(style, subaperture_grid_mask_kernel!, valid_mask, pupil, threshold_t, sub_i, sub_j, n_sub_i, n_sub_j;
        ndrange=size(valid_mask))
    return valid_mask
end
