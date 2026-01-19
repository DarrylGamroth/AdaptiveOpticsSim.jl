function pad_center!(dest::AbstractMatrix, src::AbstractMatrix)
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
    @inbounds for i in 1:n_out, j in 1:m_out
        acc = zero(eltype(out))
        for ii in 1:binning, jj in 1:binning
            acc += input[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = acc
    end
    return out
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
