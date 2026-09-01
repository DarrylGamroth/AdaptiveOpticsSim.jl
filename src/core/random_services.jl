const _SPLITMIX64_INCREMENT = UInt64(0x9e3779b97f4a7c15)

@inline function splitmix64(x::UInt64)
    z = x + _SPLITMIX64_INCREMENT
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

@noinline function _invalid_splitmix64_seed(seed::Integer)
    throw(InvalidConfiguration(
        "SplitMix64 seeds must be integers from 0 through $(typemax(UInt64)); got $seed",
    ))
end

@inline function _splitmix64_seed(seed::Integer)
    seed isa Bool && _invalid_splitmix64_seed(seed)
    (seed < 0 || seed > typemax(UInt64)) && _invalid_splitmix64_seed(seed)
    return UInt64(seed)
end

"""
    SplitMix64RNG(seed=0)

Explicit state owner for the package's host-side SplitMix64 random stream.
Retain one instance per simulation concern and mutate it from one execution
owner. The generator is deterministic for a fixed seed and does not depend on
Julia's task-local default RNG.
"""
mutable struct SplitMix64RNG <: AbstractRNG
    state::UInt64

    SplitMix64RNG(seed::Integer=0) = new(_splitmix64_seed(seed))
end

function Random.seed!(rng::SplitMix64RNG, seed::Integer)
    rng.state = _splitmix64_seed(seed)
    return rng
end

Base.copy(rng::SplitMix64RNG) = SplitMix64RNG(rng.state)

function Base.copy!(destination::SplitMix64RNG, source::SplitMix64RNG)
    destination.state = source.state
    return destination
end

Base.:(==)(left::SplitMix64RNG, right::SplitMix64RNG) =
    left.state == right.state

function Base.show(io::IO, rng::SplitMix64RNG)
    print(io, "SplitMix64RNG(state=0x", string(rng.state; base=16, pad=16), ")")
end

@inline function Random.rand(
    rng::SplitMix64RNG,
    ::Random.SamplerType{UInt64},
)
    value = splitmix64(rng.state)
    rng.state += _SPLITMIX64_INCREMENT
    return value
end

for T in (Bool, Base.BitInteger64_types...)
    T === UInt64 && continue
    @eval @inline Random.rand(
        rng::SplitMix64RNG,
        ::Random.SamplerType{$T},
    ) = rand(rng, UInt64) % $T
end

@inline Random.rand(
    rng::SplitMix64RNG,
    ::Random.SamplerType{UInt128},
) = rand(rng, UInt64) | ((rand(rng, UInt64) % UInt128) << 64)

@inline Random.rand(
    rng::SplitMix64RNG,
    ::Random.SamplerType{Int128},
) = rand(rng, UInt128) % Int128

Random.rng_native_52(::SplitMix64RNG) = UInt64

"""
    runtime_rng(seed=0)

Return a `SplitMix64RNG` for new runtime, benchmark, and RTC/HIL simulations.
The stream is stable for a fixed seed across CPU and host-driven accelerator
execution. Retain the returned object as explicit single-writer state.
"""
runtime_rng(seed::Integer=0) = SplitMix64RNG(seed)

"""
    deterministic_reference_rng(seed=0)

Return the preferred RNG for reference data, regression fixtures, and examples
where preserving the historical stream matters.
"""
deterministic_reference_rng(seed::Integer=0) = MersenneTwister(seed)

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
