@kernel function geometric_slopes_kernel!(slopes, opd, valid_mask,
    sub::Int, n_sub::Int, offset::Int, scale_x, scale_y)
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
                sx += opd[x + 1, y] - opd[x, y]
                count_x += 1
            end
            @inbounds for y in ys:(ye - 1), x in xs:xe
                sy += opd[x, y + 1] - opd[x, y]
                count_y += 1
            end
            slopes[idx] = scale_x * sx / max(count_x, 1)
            slopes[idx + offset] = scale_y * sy / max(count_y, 1)
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
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
            slopes[idx] = sx / max(count_x, 1)
            slopes[idx + offset] = sy / max(count_y, 1)
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
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

function geometric_slopes!(
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
)
    Base.require_one_based_indexing(opd, valid_mask)
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    n % n_sub == 0 || throw(DimensionMismatchError(
        "OPD size must be divisible by valid_mask size"))
    length(slopes) == 2 * n_sub * n_sub ||
        throw(DimensionMismatchError(
            "slope vector length does not match valid_mask"))
    sub = div(n, n_sub)
    offset = n_sub * n_sub
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
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    sampling_m::NTuple{2,<:Real},
)
    Base.require_one_based_indexing(opd, valid_mask)
    all(value -> isfinite(value) && value > zero(value), sampling_m) ||
        throw(InvalidConfiguration(
            "geometric wavefront sampling must be finite and positive"))
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    size(opd, 2) == n || throw(DimensionMismatchError(
        "geometric wavefront OPD must be square"))
    size(valid_mask, 2) == n_sub || throw(DimensionMismatchError(
        "geometric wavefront valid_mask must be square"))
    n % n_sub == 0 || throw(DimensionMismatchError(
        "OPD size must be divisible by valid_mask size"))
    length(slopes) == 2 * n_sub * n_sub ||
        throw(DimensionMismatchError(
            "slope vector length does not match valid_mask"))
    sub = div(n, n_sub)
    offset = n_sub * n_sub
    T = eltype(slopes)
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
    launch_kernel!(
        style, geometric_slopes_kernel!, slopes, opd, valid_mask,
        sub, n_sub, offset, scale_x, scale_y; ndrange=(n_sub, n_sub))
    return slopes
end

function edge_geometric_slopes!(
    slopes::AbstractVector,
    opd::AbstractMatrix,
    valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool},
)
    Base.require_one_based_indexing(opd, valid_mask, edge_mask)
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    n % n_sub == 0 || throw(DimensionMismatchError(
        "OPD size must be divisible by valid_mask size"))
    size(edge_mask) == size(opd) || throw(DimensionMismatchError(
        "edge_mask size must match OPD"))
    length(slopes) == 2 * n_sub * n_sub ||
        throw(DimensionMismatchError(
            "slope vector length does not match valid_mask"))
    sub = div(n, n_sub)
    offset = n_sub * n_sub
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
