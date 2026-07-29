"""
    runtime_rng(seed=0)

Return the preferred RNG for new runtime, benchmark, and RTC/HIL simulations.
The stream is repeatable for a fixed software stack, but is not intended to
preserve historical regression fixtures.
"""
runtime_rng(seed::Integer=0) = Xoshiro(seed)

"""
    deterministic_reference_rng(seed=0)

Return the preferred RNG for reference data, regression fixtures, and examples
where preserving the historical stream matters.
"""
deterministic_reference_rng(seed::Integer=0) = MersenneTwister(seed)

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

function poisson_sample(rng::AbstractRNG, λ::Real)
    if λ <= 0
        return 0
    elseif λ < 30
        limit = exp(-float(λ))
        count = 0
        product = 1.0
        while product > limit
            count += 1
            product *= rand(rng)
        end
        return count - 1
    else
        return max(
            0,
            round(Int, float(λ) + sqrt(float(λ)) * randn(rng)),
        )
    end
end
