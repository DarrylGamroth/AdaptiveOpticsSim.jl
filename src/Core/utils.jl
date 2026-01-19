function fftshift(a::AbstractMatrix)
    return circshift(a, (div(size(a, 1), 2), div(size(a, 2), 2)))
end

function ifftshift(a::AbstractMatrix)
    return circshift(a, (-div(size(a, 1), 2), -div(size(a, 2), 2)))
end

function fftshift!(dest::AbstractMatrix, src::AbstractMatrix)
    n1, n2 = size(src)
    if size(dest) != (n1, n2)
        throw(DimensionMismatchError("fftshift! requires matching sizes"))
    end
    s1 = div(n1, 2)
    s2 = div(n2, 2)
    @inbounds for i in 1:n1, j in 1:n2
        ii = i + s1
        jj = j + s2
        ii = ii > n1 ? ii - n1 : ii
        jj = jj > n2 ? jj - n2 : jj
        dest[i, j] = src[ii, jj]
    end
    return dest
end

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

backend_type(A::AbstractArray) = Base.typename(typeof(A)).wrapper
