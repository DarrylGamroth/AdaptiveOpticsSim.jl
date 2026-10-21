# Preserve the serial summation order for exceptional values and representations
# that cannot use the floating-point cooperative reduction.
@inline function _geometric_slope_sums(opd, xs, ys, sub, ::Type{T}) where {T}
    sx = zero(T)
    sy = zero(T)
    @inbounds for y in ys:(ys + sub - 1), x in xs:(xs + sub - 2)
        sx += opd[x + 1, y] - opd[x, y]
    end
    @inbounds for y in ys:(ys + sub - 2), x in xs:(xs + sub - 1)
        sy += opd[x, y + 1] - opd[x, y]
    end
    return sx, sy
end

# One workgroup owns one subaperture. Every adjacent difference still
# contributes; the reduction only changes the order of floating-point additions.
@kernel function cooperative_geometric_slopes_kernel!(slopes, opd, valid_mask,
    sub::Int, n_sub::Int, offset::Int, scale_x, scale_y)
    lane = @index(Local, Linear)
    idx = @index(Group, Linear)
    i = (idx - 1) % n_sub + 1
    j = (idx - 1) ÷ n_sub + 1
    xs = (i - 1) * sub + 1
    ys = (j - 1) * sub + 1
    sums_x = @localmem eltype(slopes) (32,)
    sums_y = @localmem eltype(slopes) (32,)
    serial = @localmem Bool (32,)
    sx = zero(eltype(slopes))
    sy = zero(eltype(slopes))
    count = sub * (sub - 1)
    # The factor two leaves headroom for rounding in either summation order.
    limit = floatmax(eltype(slopes)) / (2 * count) /
        max(one(eltype(slopes)), abs(scale_x), abs(scale_y))
    exceptional = false
    valid = @inbounds valid_mask[i, j]
    if valid
        @inbounds for k in (lane - 1):32:(count - 1)
            x = xs + k % (sub - 1)
            y = ys + k ÷ (sub - 1)
            dx = opd[x + 1, y] - opd[x, y]
            x = xs + k % sub
            y = ys + k ÷ sub
            dy = opd[x, y + 1] - opd[x, y]
            exceptional |= !isfinite(dx) || !isfinite(dy) ||
                abs(dx) > limit || abs(dy) > limit
            sx += dx
            sy += dy
        end
    end
    @inbounds sums_x[lane] = sx
    @inbounds sums_y[lane] = sy
    @inbounds serial[lane] = exceptional
    @synchronize
    for reduction_stride in (16, 8, 4, 2, 1)
        if lane <= reduction_stride
            @inbounds sums_x[lane] += sums_x[lane + reduction_stride]
            @inbounds sums_y[lane] += sums_y[lane + reduction_stride]
            @inbounds serial[lane] |= serial[lane + reduction_stride]
        end
        @synchronize
    end
    idx = @index(Group, Linear)
    i = (idx - 1) % n_sub + 1
    j = (idx - 1) ÷ n_sub + 1
    xs = (i - 1) * sub + 1
    ys = (j - 1) * sub + 1
    count = sub * (sub - 1)
    if lane == 1
        if @inbounds valid_mask[i, j]
            if @inbounds serial[1]
                sx, sy = _geometric_slope_sums(opd, xs, ys, sub, eltype(slopes))
            else
                sx = @inbounds sums_x[1]
                sy = @inbounds sums_y[1]
            end
            @inbounds slopes[idx] = scale_x * sx / count
            @inbounds slopes[idx + offset] = scale_y * sy / count
        else
            @inbounds slopes[idx] = zero(eltype(slopes))
            @inbounds slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

@kernel function geometric_slopes_kernel!(slopes, opd, valid_mask,
    sub::Int, n_sub::Int, offset::Int, scale_x, scale_y)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = i + (j - 1) * n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if @inbounds valid_mask[i, j]
            sx, sy = _geometric_slope_sums(opd, xs, ys, sub, eltype(slopes))
            count = max(sub * (sub - 1), 1)
            @inbounds slopes[idx] = scale_x * sx / count
            @inbounds slopes[idx + offset] = scale_y * sy / count
        else
            @inbounds slopes[idx] = zero(eltype(slopes))
            @inbounds slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

@kernel function edge_geometric_slopes_kernel!(
    slopes, opd, valid_mask, edge_mask,
    sub::Int, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = i + (j - 1) * n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        sx = zero(eltype(slopes))
        sy = zero(eltype(slopes))
        count_x = 0
        count_y = 0
        if @inbounds valid_mask[i, j]
            @inbounds for y in ys:ye, x in xs:(xe - 1)
                if edge_mask[x, y]
                    sx += opd[x + 1, y] - opd[x, y]
                    count_x += 1
                end
            end
            @inbounds for y in ys:(ye - 1), x in xs:xe
                if edge_mask[x, y]
                    sy += opd[x, y + 1] - opd[x, y]
                    count_y += 1
                end
            end
            @inbounds slopes[idx] = sx / max(count_x, 1)
            @inbounds slopes[idx + offset] = sy / max(count_y, 1)
        else
            @inbounds slopes[idx] = zero(eltype(slopes))
            @inbounds slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

function set_valid_subapertures!(
    valid_mask::AbstractMatrix{Bool},
    pupil::AbstractMatrix{Bool},
    threshold::Real,
)
    build_mask!(valid_mask, SubapertureGridMask(threshold), pupil)
    return valid_mask
end

function _set_valid_subapertures!(
    ::ScalarCPUStyle,
    valid_mask::AbstractMatrix{Bool},
    pupil::AbstractMatrix{Bool},
    threshold::Real,
    sub::Int,
    n_sub::Int,
)
    build_mask!(valid_mask, SubapertureGridMask(threshold), pupil)
    return valid_mask
end

function _set_valid_subapertures!(
    ::AcceleratorStyle,
    valid_mask::AbstractMatrix{Bool},
    pupil::AbstractMatrix{Bool},
    threshold::Real,
    sub::Int,
    n_sub::Int,
)
    build_mask!(valid_mask, SubapertureGridMask(threshold), pupil)
    return valid_mask
end

@inline function _geometric_slope_layout(
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
)
    Base.require_one_based_indexing(slopes, opd, valid_mask)
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    n > 0 || throw(DimensionMismatchError(
        "geometric wavefront OPD must be nonempty"))
    n_sub > 0 || throw(DimensionMismatchError(
        "geometric wavefront valid_mask must be nonempty"))
    size(opd, 2) == n || throw(DimensionMismatchError(
        "geometric wavefront OPD must be square"))
    size(valid_mask, 2) == n_sub || throw(DimensionMismatchError(
        "geometric wavefront valid_mask must be square"))
    n % n_sub == 0 || throw(DimensionMismatchError(
        "OPD size must be divisible by valid_mask size"))
    offset = n_sub * n_sub
    length(slopes) == 2 * offset ||
        throw(DimensionMismatchError(
            "slope vector length does not match valid_mask"))
    return div(n, n_sub), n_sub, offset
end

function geometric_slopes!(
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
)
    sub, n_sub, offset = _geometric_slope_layout(
        slopes, opd, valid_mask)
    unit_scale = one(eltype(slopes))
    _geometric_slopes!(
        execution_style(slopes), slopes, opd, valid_mask,
        sub, n_sub, offset, unit_scale, unit_scale)
    return slopes
end

"""
    geometric_wavefront_slopes!(slopes, opd, valid_mask, sampling_m)

Average the sampled OPD gradient over each valid subaperture. `opd` and
`sampling_m` are in metres, so the result is a dimensionless paraxial
wavefront angle conventionally reported in radians.
"""
function geometric_wavefront_slopes!(
    slopes::AbstractVector{T},
    opd::AbstractMatrix{S},
    valid_mask::AbstractMatrix{Bool},
    sampling_m::NTuple{2,<:Real},
) where {T<:AbstractFloat,S<:Real}
    all(value -> isfinite(value) && value > zero(value), sampling_m) ||
        throw(InvalidConfiguration(
            "geometric wavefront sampling must be finite and positive"))
    sub, n_sub, offset = _geometric_slope_layout(
        slopes, opd, valid_mask)
    scale_x = inv(T(sampling_m[1]))
    scale_y = inv(T(sampling_m[2]))
    _geometric_slopes!(
        execution_style(slopes), slopes, opd, valid_mask,
        sub, n_sub, offset, scale_x, scale_y)
    return slopes
end

function _geometric_slopes!(
    ::ScalarCPUStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
)
    unit_scale = one(eltype(slopes))
    return _geometric_slopes!(
        ScalarCPUStyle(), slopes, opd, valid_mask,
        sub, n_sub, offset, unit_scale, unit_scale)
end

function _geometric_slopes!(
    ::ScalarCPUStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
    scale_x,
    scale_y,
)
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        if valid_mask[i, j]
            sx = zero(eltype(slopes))
            sy = zero(eltype(slopes))
            count_x = 0
            count_y = 0
            for y in ys:ye, x in xs:(xe - 1)
                sx += opd[x + 1, y] - opd[x, y]
                count_x += 1
            end
            for y in ys:(ye - 1), x in xs:xe
                sy += opd[x, y + 1] - opd[x, y]
                count_y += 1
            end
            slopes[idx] = scale_x * sx / max(count_x, 1)
            slopes[idx + offset] = scale_y * sy / max(count_y, 1)
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
        end
        idx += 1
    end
    return slopes
end

function _geometric_slopes!(
    style::AcceleratorStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
)
    unit_scale = one(eltype(slopes))
    return _geometric_slopes!(
        style, slopes, opd, valid_mask,
        sub, n_sub, offset, unit_scale, unit_scale)
end

@inline _cooperative_geometric_slopes(slopes, opd, scale_x, scale_y) = false

@inline function _cooperative_geometric_slopes(
    ::AbstractVector{T}, ::AbstractMatrix{T}, scale_x::T, scale_y::T,
) where {T<:Union{Float32,Float64}}
    return isfinite(scale_x) && isfinite(scale_y)
end

function _geometric_slopes!(
    style::AcceleratorStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
    scale_x,
    scale_y,
)
    _enqueue_geometric_slopes!(style, slopes, opd, valid_mask,
        sub, n_sub, offset, scale_x, scale_y)
    synchronize_backend!(style)
    return slopes
end

function _enqueue_geometric_slopes!(
    style::AcceleratorStyle, slopes, opd, valid_mask,
    sub::Int, n_sub::Int, offset::Int, scale_x, scale_y,
)
    # Small subapertures do too little work to amortize one group per output.
    # Dense output grids already expose enough serial-kernel parallelism; keep
    # this reduction within the measured low-concurrency range.
    if sub >= 5 && n_sub <= 32 && _cooperative_geometric_slopes(slopes, opd, scale_x, scale_y)
        # The fixed workgroup size is part of this reduction's storage contract.
        cooperative_geometric_slopes_kernel!(style.backend, 32)(
            slopes, opd, valid_mask, sub, n_sub, offset, scale_x, scale_y;
            ndrange=32 * offset)
    else
        launch_kernel_async!(
            style, geometric_slopes_kernel!, slopes, opd, valid_mask,
            sub, n_sub, offset, scale_x, scale_y; ndrange=(n_sub, n_sub))
    end
    return slopes
end

function edge_geometric_slopes!(
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool},
)
    Base.require_one_based_indexing(edge_mask)
    sub, n_sub, offset = _geometric_slope_layout(
        slopes, opd, valid_mask)
    size(edge_mask) == size(opd) || throw(DimensionMismatchError(
        "edge_mask size must match OPD"))
    _edge_geometric_slopes!(
        execution_style(slopes), slopes, opd, valid_mask,
        edge_mask, sub, n_sub, offset)
    return slopes
end

function _edge_geometric_slopes!(
    ::ScalarCPUStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
)
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        if valid_mask[i, j]
            sx = zero(eltype(slopes))
            sy = zero(eltype(slopes))
            count_x = 0
            count_y = 0
            for y in ys:ye, x in xs:(xe - 1)
                if edge_mask[x, y]
                    sx += opd[x + 1, y] - opd[x, y]
                    count_x += 1
                end
            end
            for y in ys:(ye - 1), x in xs:xe
                if edge_mask[x, y]
                    sy += opd[x, y + 1] - opd[x, y]
                    count_y += 1
                end
            end
            slopes[idx] = sx / max(count_x, 1)
            slopes[idx + offset] = sy / max(count_y, 1)
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
        end
        idx += 1
    end
    return slopes
end

function _edge_geometric_slopes!(
    style::AcceleratorStyle,
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool},
    sub::Int,
    n_sub::Int,
    offset::Int,
)
    launch_kernel!(
        style, edge_geometric_slopes_kernel!, slopes, opd,
        valid_mask, edge_mask, sub, n_sub, offset;
        ndrange=(n_sub, n_sub))
    return slopes
end
