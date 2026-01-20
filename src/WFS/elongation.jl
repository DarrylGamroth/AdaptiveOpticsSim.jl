@inline function ensure_kernel(kernel::Vector{T}, needed::Int) where {T}
    if length(kernel) < needed
        resize!(kernel, needed)
    end
    return kernel
end

@inline function ensure_kernel(kernel::AbstractVector{T}, needed::Int) where {T}
    if length(kernel) < needed
        kernel = similar(kernel, needed)
    end
    return kernel
end

function apply_elongation!(intensity::AbstractMatrix{T}, factor::Real, tmp::AbstractMatrix{T},
    kernel::AbstractVector{T}) where {T<:AbstractFloat}
    if factor <= 1
        return kernel
    end
    sigma = T(0.5) * (T(factor) - one(T))
    if sigma <= 0
        return kernel
    end
    half = max(1, ceil(Int, 2 * sigma))
    needed = 2 * half + 1
    kernel = ensure_kernel(kernel, needed)
    @inbounds for k in -half:half
        kernel[k + half + 1] = exp(-T(0.5) * (T(k) / sigma)^2)
    end
    view_kernel = @view kernel[1:needed]
    view_kernel ./= sum(view_kernel)

    n1, n2 = size(intensity)
    @inbounds for i in 1:n1, j in 1:n2
        acc = zero(T)
        for k in -half:half
            jj = clamp(j + k, 1, n2)
            acc += intensity[i, jj] * kernel[k + half + 1]
        end
        tmp[i, j] = acc
    end
    copyto!(intensity, tmp)
    return kernel
end
