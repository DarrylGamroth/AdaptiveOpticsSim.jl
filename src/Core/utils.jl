@kernel function fftshift2d_kernel!(dest, src, sx::Int, sy::Int, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        ii = i + sx
        jj = j + sy
        ii = ii > n ? ii - n : ii
        jj = jj > m ? jj - m : jj
        @inbounds dest[i, j] = src[ii, jj]
    end
end

@kernel function circshift2d_kernel!(dest, src, dx::Int, dy::Int, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        ii = mod1(i - dx, n)
        jj = mod1(j - dy, m)
        @inbounds dest[i, j] = src[ii, jj]
    end
end

@kernel function bin2d_kernel!(out, input, binning::Int, n_out::Int, m_out::Int)
    i, j = @index(Global, NTuple)
    if i <= n_out && j <= m_out
        acc = zero(eltype(out))
        @inbounds for ii in 1:binning, jj in 1:binning
            acc += input[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function fftfreq_kernel!(dest, n::Int, val, offset)
    i = @index(Global, Linear)
    if i <= n
        k = i - 1
        freq = ifelse(k <= n ÷ 2, k * val, (k - n) * val)
        @inbounds dest[i] = freq + offset
    end
end

@kernel function valid_subaperture_mask_kernel!(valid_mask, pupil, threshold, sub::Int, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        occ = zero(typeof(threshold))
        @inbounds for x in xs:(xs + sub - 1), y in ys:(ys + sub - 1)
            occ += pupil[x, y]
        end
        @inbounds valid_mask[i, j] = occ / (sub * sub) > threshold
    end
end

@kernel function geometric_slopes_kernel!(slopes, opd, valid_mask, sub::Int, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        sx = zero(eltype(slopes))
        sy = zero(eltype(slopes))
        count_x = 0
        count_y = 0
        if @inbounds valid_mask[i, j]
            @inbounds for x in xs:(xe - 1), y in ys:ye
                sx += opd[x + 1, y] - opd[x, y]
                count_x += 1
            end
            @inbounds for x in xs:xe, y in ys:(ye - 1)
                sy += opd[x, y + 1] - opd[x, y]
                count_y += 1
            end
            slopes[idx] = sx / max(count_x, 1)
            slopes[idx + offset] = sy / max(count_y, 1)
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

@kernel function edge_geometric_slopes_kernel!(slopes, opd, valid_mask, edge_mask, sub::Int, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        sx = zero(eltype(slopes))
        sy = zero(eltype(slopes))
        count_x = 0
        count_y = 0
        if @inbounds valid_mask[i, j]
            @inbounds for x in xs:(xe - 1), y in ys:ye
                if edge_mask[x, y]
                    sx += opd[x + 1, y] - opd[x, y]
                    count_x += 1
                end
            end
            @inbounds for x in xs:xe, y in ys:(ye - 1)
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

@kernel function randn_fill_kernel!(out, seed::UInt64, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        u1 = uniform01(T, splitmix64(seed + UInt64(2 * i - 1)))
        u2 = uniform01(T, splitmix64(seed + UInt64(2 * i)))
        radius = sqrt(T(-2) * log(u1))
        phase = T(2 * pi) * u2
        @inbounds out[i] = radius * cos(phase)
    end
end

@inline function splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

@inline function uniform01(::Type{T}, x::UInt64) where {T<:AbstractFloat}
    u = T(ldexp(Float64(x >>> 11), -53))
    return clamp(u, eps(T), prevfloat(one(T)))
end

function pad_center!(dest::AbstractMatrix, src::AbstractMatrix)
    Base.require_one_based_indexing(dest, src)
    fill!(dest, zero(eltype(dest)))
    nx, ny = size(dest)
    sx, sy = size(src)
    if sx > nx || sy > ny
        throw(DimensionMismatchError("source larger than destination"))
    end
    ox = div(nx - sx, 2)
    oy = div(ny - sy, 2)
    @views dest[ox+1:ox+sx, oy+1:oy+sy] .= src
    return dest
end

function center_resize2d!(dest::AbstractMatrix, src::AbstractMatrix)
    Base.require_one_based_indexing(dest, src)
    fill!(dest, zero(eltype(dest)))
    nx, ny = size(dest)
    sx, sy = size(src)
    ox = div(nx - sx, 2)
    oy = div(ny - sy, 2)
    dx = max(1, 1 + ox)
    dy = max(1, 1 + oy)
    sx_start = max(1, 1 - ox)
    sy_start = max(1, 1 - oy)
    len_x = min(nx - dx + 1, sx - sx_start + 1)
    len_y = min(ny - dy + 1, sy - sy_start + 1)
    if len_x > 0 && len_y > 0
        @views dest[dx:dx+len_x-1, dy:dy+len_y-1] .= src[sx_start:sx_start+len_x-1, sy_start:sy_start+len_y-1]
    end
    return dest
end

function fftshift2d!(dest::AbstractMatrix, src::AbstractMatrix)
    Base.require_one_based_indexing(dest, src)
    n, m = size(src)
    if size(dest) != (n, m)
        throw(DimensionMismatchError("fftshift2d! destination size must match source"))
    end
    if dest === src
        tmp = similar(src)
        fftshift2d!(tmp, src)
        copyto!(dest, tmp)
        return dest
    end
    _fftshift2d!(execution_style(dest), dest, src)
    return dest
end

function _fftshift2d!(::ScalarCPUStyle, dest::AbstractMatrix, src::AbstractMatrix)
    n, m = size(src)
    sx = n ÷ 2
    sy = m ÷ 2
    @views begin
        dest[1:n-sx, 1:m-sy] .= src[sx+1:n, sy+1:m]
        dest[1:n-sx, m-sy+1:m] .= src[sx+1:n, 1:sy]
        dest[n-sx+1:n, 1:m-sy] .= src[1:sx, sy+1:m]
        dest[n-sx+1:n, m-sy+1:m] .= src[1:sx, 1:sy]
    end
    return dest
end

function _fftshift2d!(style::AcceleratorStyle, dest::AbstractMatrix, src::AbstractMatrix)
    launch_kernel!(style, fftshift2d_kernel!, dest, src, size(src, 1) ÷ 2, size(src, 2) ÷ 2, size(src, 1), size(src, 2);
        ndrange=size(dest))
    return dest
end

function circshift2d!(dest::AbstractMatrix, src::AbstractMatrix, shifts::NTuple{2,Int})
    Base.require_one_based_indexing(dest, src)
    if size(dest) != size(src)
        throw(DimensionMismatchError("circshift2d! destination size must match source"))
    end
    if dest === src
        tmp = similar(src)
        circshift2d!(tmp, src, shifts)
        copyto!(dest, tmp)
        return dest
    end
    _circshift2d!(execution_style(dest), dest, src, shifts)
    return dest
end

function _circshift2d!(::ScalarCPUStyle, dest::AbstractMatrix, src::AbstractMatrix, shifts::NTuple{2,Int})
    copyto!(dest, circshift(src, shifts))
    return dest
end

function _circshift2d!(style::AcceleratorStyle, dest::AbstractMatrix, src::AbstractMatrix, shifts::NTuple{2,Int})
    launch_kernel!(style, circshift2d_kernel!, dest, src, shifts[1], shifts[2], size(src, 1), size(src, 2);
        ndrange=size(dest))
    return dest
end

function bin2d(input::AbstractMatrix, binning::Int)
    Base.require_one_based_indexing(input)
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if binning == 1
        return copy(input)
    end
    n, m = size(input)
    n_out = div(n, binning)
    m_out = div(m, binning)
    out = similar(input, n_out, m_out)
    @inbounds for i in 1:n_out, j in 1:m_out
        acc = zero(eltype(out))
        for ii in 1:binning, jj in 1:binning
            acc += input[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = acc
    end
    return out
end

function bin2d!(out::AbstractMatrix, input::AbstractMatrix, binning::Int)
    Base.require_one_based_indexing(out, input)
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if binning == 1
        if size(out) != size(input)
            throw(DimensionMismatchError("output size must match input when binning=1"))
        end
        out .= input
        return out
    end
    n, m = size(input)
    n_out = div(n, binning)
    m_out = div(m, binning)
    if size(out) != (n_out, m_out)
        throw(DimensionMismatchError("output size does not match binned dimensions"))
    end
    _bin2d!(execution_style(out), out, input, binning)
    return out
end

function _bin2d!(::ScalarCPUStyle, out::AbstractMatrix, input::AbstractMatrix, binning::Int)
    n_out, m_out = size(out)
    @inbounds for i in 1:n_out, j in 1:m_out
        acc = zero(eltype(out))
        for ii in 1:binning, jj in 1:binning
            acc += input[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = acc
    end
    return out
end

function _bin2d!(style::AcceleratorStyle, out::AbstractMatrix, input::AbstractMatrix, binning::Int)
    launch_kernel!(style, bin2d_kernel!, out, input, binning, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function fftfreq!(dest::AbstractVector, n::Int; d::Real=1, offset::Real=0)
    Base.require_one_based_indexing(dest)
    if length(dest) != n
        throw(DimensionMismatchError("fftfreq! destination length must match n"))
    end
    _fftfreq!(execution_style(dest), dest, n, d, offset)
    return dest
end

function _fftfreq!(::ScalarCPUStyle, dest::AbstractVector, n::Int, d::Real, offset::Real)
    val = 1 / (n * d)
    @inbounds for i in 1:n
        k = i - 1
        if k <= n ÷ 2
            dest[i] = k * val
        else
            dest[i] = (k - n) * val
        end
    end
    if offset != 0
        dest .+= offset
    end
    return dest
end

function _fftfreq!(style::AcceleratorStyle, dest::AbstractVector, n::Int, d::Real, offset::Real)
    val = eltype(dest)(1 / (n * d))
    offset_t = eltype(dest)(offset)
    launch_kernel!(style, fftfreq_kernel!, dest, n, val, offset_t; ndrange=length(dest))
    return dest
end

function poisson_sample(rng::AbstractRNG, λ::Real)
    if λ <= 0
        return 0
    elseif λ < 30
        L = exp(-float(λ))
        k = 0
        p = 1.0
        while p > L
            k += 1
            p *= rand(rng)
        end
        return k - 1
    else
        return max(0, round(Int, float(λ) + sqrt(float(λ)) * randn(rng)))
    end
end

function poisson_noise!(rng::AbstractRNG, img::AbstractMatrix)
    _poisson_noise!(execution_style(img), rng, img)
    return img
end

function _poisson_noise!(::ScalarCPUStyle, rng::AbstractRNG, img::AbstractMatrix)
    @inbounds for i in eachindex(img)
        img[i] = poisson_sample(rng, img[i])
    end
    return img
end

function _poisson_noise!(::AcceleratorStyle, rng::AbstractRNG, img::AbstractMatrix{T}) where {T}
    host = Matrix{T}(undef, size(img))
    copyto!(host, img)
    _poisson_noise!(ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end

function randn_backend!(rng::AbstractRNG, out::AbstractArray)
    _randn_backend!(execution_style(out), rng, out)
    return out
end

function _randn_backend!(::ScalarCPUStyle, rng::AbstractRNG, out::AbstractArray)
    randn!(rng, out)
    return out
end

function _randn_backend!(style::AcceleratorStyle, rng::AbstractRNG, out::AbstractArray{T}) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel!(style, randn_fill_kernel!, out, seed, length(out); ndrange=length(out))
    return out
end

function set_valid_subapertures!(valid_mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool}, threshold::Real)
    Base.require_one_based_indexing(valid_mask, pupil)
    n = size(pupil, 1)
    n_sub = size(valid_mask, 1)
    if size(valid_mask, 2) != n_sub || size(pupil, 2) != n
        throw(DimensionMismatchError("valid_mask and pupil must be square"))
    end
    if n % n_sub != 0
        throw(DimensionMismatchError("pupil size must be divisible by valid_mask size"))
    end
    sub = div(n, n_sub)
    _set_valid_subapertures!(execution_style(valid_mask), valid_mask, pupil, threshold, sub, n_sub)
    return valid_mask
end

function _set_valid_subapertures!(::ScalarCPUStyle, valid_mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool},
    threshold::Real, sub::Int, n_sub::Int)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        valid_mask[i, j] = mean(@view pupil[xs:xe, ys:ye]) > threshold
    end
    return valid_mask
end

function _set_valid_subapertures!(style::AcceleratorStyle, valid_mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool},
    threshold::Real, sub::Int, n_sub::Int)
    threshold_t = promote_type(eltype(pupil), typeof(threshold))(threshold)
    launch_kernel!(style, valid_subaperture_mask_kernel!, valid_mask, pupil, threshold_t, sub, n_sub; ndrange=size(valid_mask))
    return valid_mask
end

function geometric_slopes!(slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool})
    Base.require_one_based_indexing(opd, valid_mask)
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    if n % n_sub != 0
        throw(DimensionMismatchError("OPD size must be divisible by valid_mask size"))
    end
    if length(slopes) != 2 * n_sub * n_sub
        throw(DimensionMismatchError("slope vector length does not match valid_mask"))
    end
    sub = div(n, n_sub)
    offset = n_sub * n_sub
    _geometric_slopes!(execution_style(slopes), slopes, opd, valid_mask, sub, n_sub, offset)
    return slopes
end

function _geometric_slopes!(::ScalarCPUStyle, slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, offset::Int)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        if valid_mask[i, j]
            sx = zero(eltype(slopes))
            sy = zero(eltype(slopes))
            count_x = 0
            count_y = 0
            for x in xs:(xe - 1), y in ys:ye
                sx += opd[x + 1, y] - opd[x, y]
                count_x += 1
            end
            for x in xs:xe, y in ys:(ye - 1)
                sy += opd[x, y + 1] - opd[x, y]
                count_y += 1
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

function _geometric_slopes!(style::AcceleratorStyle, slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, offset::Int)
    launch_kernel!(style, geometric_slopes_kernel!, slopes, opd, valid_mask, sub, n_sub, offset; ndrange=(n_sub, n_sub))
    return slopes
end

function edge_geometric_slopes!(slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool}, edge_mask::AbstractMatrix{Bool})
    Base.require_one_based_indexing(opd, valid_mask, edge_mask)
    n = size(opd, 1)
    n_sub = size(valid_mask, 1)
    if n % n_sub != 0
        throw(DimensionMismatchError("OPD size must be divisible by valid_mask size"))
    end
    if size(edge_mask) != size(opd)
        throw(DimensionMismatchError("edge_mask size must match OPD"))
    end
    if length(slopes) != 2 * n_sub * n_sub
        throw(DimensionMismatchError("slope vector length does not match valid_mask"))
    end
    sub = div(n, n_sub)
    offset = n_sub * n_sub
    _edge_geometric_slopes!(execution_style(slopes), slopes, opd, valid_mask, edge_mask, sub, n_sub, offset)
    return slopes
end

function _edge_geometric_slopes!(::ScalarCPUStyle, slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool}, sub::Int, n_sub::Int, offset::Int)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = xs + sub - 1
        ye = ys + sub - 1
        if valid_mask[i, j]
            sx = zero(eltype(slopes))
            sy = zero(eltype(slopes))
            count_x = 0
            count_y = 0
            for x in xs:(xe - 1), y in ys:ye
                if edge_mask[x, y]
                    sx += opd[x + 1, y] - opd[x, y]
                    count_x += 1
                end
            end
            for x in xs:xe, y in ys:(ye - 1)
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

function _edge_geometric_slopes!(style::AcceleratorStyle, slopes::AbstractVector, opd::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    edge_mask::AbstractMatrix{Bool}, sub::Int, n_sub::Int, offset::Int)
    launch_kernel!(style, edge_geometric_slopes_kernel!, slopes, opd, valid_mask, edge_mask, sub, n_sub, offset; ndrange=(n_sub, n_sub))
    return slopes
end
