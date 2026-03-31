@kernel function randn_fill_kernel!(out, seed::UInt64, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        @inbounds out[i] = normal01(T, seed + UInt64(2 * i - 1), seed + UInt64(2 * i))
    end
end

@kernel function poisson_noise_kernel!(img, seed::UInt64, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(img)
        λ = @inbounds img[i]
        if λ <= zero(T)
            @inbounds img[i] = zero(T)
        elseif λ < T(30)
            limit = exp(-λ)
            prod = one(T)
            count = 0
            key = seed + UInt64(0x9e3779b97f4a7c15) * UInt64(i)
            while prod > limit
                count += 1
                prod *= uniform01(T, splitmix64(key + UInt64(count)))
            end
            @inbounds img[i] = T(count - 1)
        else
            z = normal01(T, seed + UInt64(2 * i - 1), seed + UInt64(2 * i))
            sample = floor(λ + sqrt(λ) * z + T(0.5))
            @inbounds img[i] = max(zero(T), sample)
        end
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

@inline function normal01(::Type{T}, x1::UInt64, x2::UInt64) where {T<:AbstractFloat}
    u1 = uniform01(T, splitmix64(x1))
    u2 = uniform01(T, splitmix64(x2))
    radius = sqrt(T(-2) * log(u1))
    phase = T(2 * pi) * u2
    return radius * cos(phase)
end

function poisson_noise!(rng::AbstractRNG, img::AbstractArray)
    _poisson_noise!(execution_style(img), rng, img)
    return img
end

function _poisson_noise!(::ScalarCPUStyle, rng::AbstractRNG, img::AbstractArray)
    @inbounds for i in eachindex(img)
        img[i] = poisson_sample(rng, img[i])
    end
    return img
end

function _poisson_noise!(style::AcceleratorStyle, rng::AbstractRNG, img::AbstractArray{T}) where {T<:AbstractFloat}
    poisson_noise_async!(style, rng, img)
    synchronize_backend!(style)
    return img
end

function poisson_noise_async!(style::AcceleratorStyle, rng::AbstractRNG, img::AbstractArray{T}) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel_async!(style, poisson_noise_kernel!, img, seed, length(img); ndrange=length(img))
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
    randn_backend_async!(style, rng, out)
    synchronize_backend!(style)
    return out
end

function randn_backend_async!(style::AcceleratorStyle, rng::AbstractRNG, out::AbstractArray{T}) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel_async!(style, randn_fill_kernel!, out, seed, length(out); ndrange=length(out))
    return out
end
