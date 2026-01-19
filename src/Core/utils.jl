function pad_center!(dest::AbstractMatrix, src::AbstractMatrix)
    fill!(dest, zero(eltype(dest)))
    nx, ny = size(dest)
    sx, sy = size(src)
    if sx > nx || sy > ny
        throw(DimensionMismatchError("source larger than destination"))
    end
    ox = div(nx - sx, 2)
    oy = div(ny - sy, 2)
    dx0 = firstindex(dest, 1) + ox
    dy0 = firstindex(dest, 2) + oy
    dest_inds = CartesianIndices((dx0:dx0+sx-1, dy0:dy0+sy-1))
    copyto!(dest, dest_inds, src, CartesianIndices(src))
    return dest
end

function bin2d(input::AbstractMatrix, binning::Int)
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
    i0_in = firstindex(input, 1)
    j0_in = firstindex(input, 2)
    i0_out = firstindex(out, 1)
    j0_out = firstindex(out, 2)
    @inbounds for oi in 1:n_out, oj in 1:m_out
        acc = zero(eltype(out))
        base_i = i0_in + (oi - 1) * binning
        base_j = j0_in + (oj - 1) * binning
        for ii in 1:binning, jj in 1:binning
            acc += input[base_i + (ii - 1), base_j + (jj - 1)]
        end
        out[i0_out + (oi - 1), j0_out + (oj - 1)] = acc
    end
    return out
end

function bin2d!(out::AbstractMatrix, input::AbstractMatrix, binning::Int)
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if binning == 1
        if size(out) != size(input)
            throw(DimensionMismatchError("output size must match input when binning=1"))
        end
        copyto!(out, CartesianIndices(out), input, CartesianIndices(input))
        return out
    end
    n, m = size(input)
    n_out = div(n, binning)
    m_out = div(m, binning)
    if size(out) != (n_out, m_out)
        throw(DimensionMismatchError("output size does not match binned dimensions"))
    end
    i0_in = firstindex(input, 1)
    j0_in = firstindex(input, 2)
    i0_out = firstindex(out, 1)
    j0_out = firstindex(out, 2)
    @inbounds for oi in 1:n_out, oj in 1:m_out
        acc = zero(eltype(out))
        base_i = i0_in + (oi - 1) * binning
        base_j = j0_in + (oj - 1) * binning
        for ii in 1:binning, jj in 1:binning
            acc += input[base_i + (ii - 1), base_j + (jj - 1)]
        end
        out[i0_out + (oi - 1), j0_out + (oj - 1)] = acc
    end
    return out
end

function fftfreq!(dest::AbstractVector, n::Int; d::Real=1, offset::Real=0)
    if length(dest) != n
        throw(DimensionMismatchError("fftfreq! destination length must match n"))
    end
    val = 1 / (n * d)
    i0 = firstindex(dest)
    @inbounds for i in axes(dest, 1)
        k = i - i0
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
