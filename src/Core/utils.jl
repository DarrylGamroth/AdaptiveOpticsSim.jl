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
    @inbounds for i in 1:n_out, j in 1:m_out
        acc = zero(eltype(out))
        for ii in 1:binning, jj in 1:binning
            acc += input[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = acc
    end
    return out
end

function fftfreq!(dest::AbstractVector, n::Int; d::Real=1, offset::Real=0)
    Base.require_one_based_indexing(dest)
    if length(dest) != n
        throw(DimensionMismatchError("fftfreq! destination length must match n"))
    end
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
