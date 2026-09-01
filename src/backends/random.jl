const _RANDOM_DOMAIN_NORMAL = UInt64(0x6e6f726d616c5f31)
const _RANDOM_DOMAIN_UNIFORM = UInt64(0x756e69666f726d31)
const _RANDOM_DOMAIN_POISSON = UInt64(0x706f6973736f6e31)
const _RANDOM_ELEMENT_STRIDE = UInt64(0x9e3779b97f4a7c15)

"""Run-immutable seed for an accelerator graph owner's random stream."""
struct _CounterRandomPlan
    seed::UInt64
end

"""Persistent, single-writer draw sequence stored on the target device."""
mutable struct _CounterRandomState{A<:AbstractVector{UInt64}}
    draw_sequence::A
end

"""
Internal array-only RNG for prepared accelerator graph owners.

The scientific seed is an immutable plan value. The evolving draw sequence is
device-resident state, so replaying a captured accelerator graph advances the
stream instead of replaying a host seed captured by value. This type does not
implement scalar `rand`; stochastic kernels address samples by draw sequence,
domain, element, and sub-draw.
"""
struct _PreparedCounterRNG{S<:_CounterRandomState} <: AbstractRNG
    plan::_CounterRandomPlan
    state::S
end

@kernel function _advance_random_draw_kernel!(draw_sequence)
    i = @index(Global, Linear)
    if i == 1
        @inbounds draw_sequence[1] += UInt64(1)
    end
end

@kernel function _reset_random_draw_kernel!(draw_sequence)
    i = @index(Global, Linear)
    if i == 1
        @inbounds draw_sequence[1] = UInt64(0)
    end
end

@inline function _counter_draw_key(
    seed::UInt64,
    draw_sequence::UInt64,
    domain::UInt64,
)
    return splitmix64(
        xor(xor(seed, splitmix64(draw_sequence)), splitmix64(domain)),
    )
end

@inline function _counter_element_key(draw_key::UInt64, index::Int)
    return splitmix64(draw_key + _RANDOM_ELEMENT_STRIDE * UInt64(index))
end

@kernel function _counter_randn_fill_kernel!(
    out,
    draw_sequence,
    seed::UInt64,
    n::Int,
)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        draw_key = _counter_draw_key(
            seed,
            @inbounds(draw_sequence[1]),
            _RANDOM_DOMAIN_NORMAL,
        )
        element_key = _counter_element_key(draw_key, i)
        @inbounds out[i] = normal01(
            T,
            element_key,
            element_key + _RANDOM_ELEMENT_STRIDE,
        )
    end
end

@kernel function _counter_uniform_fill_kernel!(
    out,
    draw_sequence,
    seed::UInt64,
    n::Int,
)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        draw_key = _counter_draw_key(
            seed,
            @inbounds(draw_sequence[1]),
            _RANDOM_DOMAIN_UNIFORM,
        )
        @inbounds out[i] = uniform01(
            T,
            _counter_element_key(draw_key, i),
        )
    end
end

@kernel function _counter_poisson_noise_kernel!(
    img,
    draw_sequence,
    seed::UInt64,
    n::Int,
)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(img)
        λ = @inbounds img[i]
        if λ <= zero(T)
            @inbounds img[i] = zero(T)
        else
            draw_key = _counter_draw_key(
                seed,
                @inbounds(draw_sequence[1]),
                _RANDOM_DOMAIN_POISSON,
            )
            element_key = _counter_element_key(draw_key, i)
            if λ < T(30)
                limit = exp(-λ)
                product = one(T)
                count = 0
                while product > limit
                    count += 1
                    product *= uniform01(
                        T,
                        splitmix64(element_key + UInt64(count)),
                    )
                end
                @inbounds img[i] = T(count - 1)
            else
                z = normal01(
                    T,
                    element_key,
                    element_key + _RANDOM_ELEMENT_STRIDE,
                )
                sample = floor(λ + sqrt(λ) * z + T(0.5))
                @inbounds img[i] = max(zero(T), sample)
            end
        end
    end
end

function _prepare_graph_rng(device::HostComputeDevice, seed::UInt64)
    return Xoshiro(seed)
end

function _prepare_graph_rng(
    ::AcceleratorComputeDevice,
    seed::UInt64,
)
    return Xoshiro(seed)
end

function _prepare_counter_rng(
    device::AcceleratorComputeDevice,
    seed::UInt64,
)
    draw_sequence = allocate_device_array(device, UInt64, 1)
    fill!(draw_sequence, UInt64(0))
    state = _CounterRandomState(draw_sequence)
    return _PreparedCounterRNG(_CounterRandomPlan(seed), state)
end

function _reset_graph_rng!(rng::AbstractRNG, seed::UInt64)
    Random.seed!(rng, seed)
    return rng
end

function _reset_graph_rng!(rng::_PreparedCounterRNG, seed::UInt64)
    seed == rng.plan.seed || throw(InvalidConfiguration(
        "a prepared counter RNG can only be reset to its construction seed",
    ))
    draw_sequence = rng.state.draw_sequence
    style = execution_style(draw_sequence)
    _reset_random_draw!(style, draw_sequence)
    return rng
end

function _reset_random_draw!(
    ::ScalarCPUStyle,
    draw_sequence::AbstractVector{UInt64},
)
    fill!(draw_sequence, UInt64(0))
    return draw_sequence
end

function _reset_random_draw!(
    style::AcceleratorStyle,
    draw_sequence::AbstractVector{UInt64},
)
    launch_kernel_async!(
        style,
        _reset_random_draw_kernel!,
        draw_sequence;
        ndrange=1,
    )
    return draw_sequence
end

@inline function _queue_counter_draw!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
)
    launch_kernel_async!(
        style,
        _advance_random_draw_kernel!,
        rng.state.draw_sequence;
        ndrange=1,
    )
    return nothing
end

@kernel function randn_fill_kernel!(out, seed::UInt64, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        @inbounds out[i] = normal01(
            T,
            seed + UInt64(2 * i - 1),
            seed + UInt64(2 * i),
        )
    end
end

function _counter_poisson_noise_async!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    img::AbstractArray{T},
) where {T<:AbstractFloat}
    isempty(img) && return img
    _queue_counter_draw!(style, rng)
    launch_kernel_async!(
        style,
        _counter_poisson_noise_kernel!,
        img,
        rng.state.draw_sequence,
        rng.plan.seed,
        length(img);
        ndrange=length(img),
    )
    return img
end

function _counter_randn_backend_async!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    isempty(out) && return out
    _queue_counter_draw!(style, rng)
    launch_kernel_async!(
        style,
        _counter_randn_fill_kernel!,
        out,
        rng.state.draw_sequence,
        rng.plan.seed,
        length(out);
        ndrange=length(out),
    )
    return out
end

function _counter_uniform_backend_async!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    isempty(out) && return out
    _queue_counter_draw!(style, rng)
    launch_kernel_async!(
        style,
        _counter_uniform_fill_kernel!,
        out,
        rng.state.draw_sequence,
        rng.plan.seed,
        length(out);
        ndrange=length(out),
    )
    return out
end

@kernel function uniform_fill_kernel!(out, seed::UInt64, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(out)
        @inbounds out[i] = uniform01(T, splitmix64(seed + UInt64(i)))
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
            product = one(T)
            count = 0
            key = seed + UInt64(0x9e3779b97f4a7c15) * UInt64(i)
            while product > limit
                count += 1
                product *= uniform01(T, splitmix64(key + UInt64(count)))
            end
            @inbounds img[i] = T(count - 1)
        else
            z = normal01(
                T,
                seed + UInt64(2 * i - 1),
                seed + UInt64(2 * i),
            )
            sample = floor(λ + sqrt(λ) * z + T(0.5))
            @inbounds img[i] = max(zero(T), sample)
        end
    end
end

function poisson_noise!(rng::AbstractRNG, img::AbstractArray)
    _poisson_noise!(execution_style(img), rng, img)
    return img
end

function _poisson_noise!(
    ::ScalarCPUStyle,
    rng::AbstractRNG,
    img::AbstractArray,
)
    @inbounds for i in eachindex(img)
        img[i] = poisson_sample(rng, img[i])
    end
    return img
end

function _poisson_noise!(
    style::AcceleratorStyle,
    rng::AbstractRNG,
    img::AbstractArray{T},
) where {T<:AbstractFloat}
    poisson_noise_async!(style, rng, img)
    synchronize_backend!(style)
    return img
end

function poisson_noise_async!(
    style::AcceleratorStyle,
    rng::AbstractRNG,
    img::AbstractArray{T},
) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel_async!(
        style,
        poisson_noise_kernel!,
        img,
        seed,
        length(img);
        ndrange=length(img),
    )
    return img
end

function poisson_noise_async!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    img::AbstractArray{T},
) where {T<:AbstractFloat}
    return _counter_poisson_noise_async!(style, rng, img)
end

function randn_backend!(rng::AbstractRNG, out::AbstractArray)
    _randn_backend!(execution_style(out), rng, out)
    return out
end

function _randn_backend!(
    ::ScalarCPUStyle,
    rng::AbstractRNG,
    out::AbstractArray,
)
    randn!(rng, out)
    return out
end

function _randn_backend!(
    style::AcceleratorStyle,
    rng::AbstractRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    randn_backend_async!(style, rng, out)
    synchronize_backend!(style)
    return out
end

function randn_backend_async!(
    style::AcceleratorStyle,
    rng::AbstractRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel_async!(
        style,
        randn_fill_kernel!,
        out,
        seed,
        length(out);
        ndrange=length(out),
    )
    return out
end

function randn_backend_async!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    return _counter_randn_backend_async!(style, rng, out)
end

function rand_uniform_backend!(rng::AbstractRNG, out::AbstractArray)
    _rand_uniform_backend!(execution_style(out), rng, out)
    return out
end

function _rand_uniform_backend!(
    ::ScalarCPUStyle,
    rng::AbstractRNG,
    out::AbstractArray,
)
    rand!(rng, out)
    return out
end


function _rand_uniform_backend!(
    style::AcceleratorStyle,
    rng::_PreparedCounterRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    _counter_uniform_backend_async!(style, rng, out)
    synchronize_backend!(style)
    return out
end

function _rand_uniform_backend!(
    style::AcceleratorStyle,
    rng::AbstractRNG,
    out::AbstractArray{T},
) where {T<:AbstractFloat}
    seed = rand(rng, UInt64)
    launch_kernel!(
        style,
        uniform_fill_kernel!,
        out,
        seed,
        length(out);
        ndrange=length(out),
    )
    return out
end
