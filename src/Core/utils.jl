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
    @inbounds for i in eachindex(img)
        img[i] = poisson_sample(rng, img[i])
    end
    return img
end
