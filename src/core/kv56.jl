@inline _kv56_order(::Type{T}) where {T<:AbstractFloat} = T(5) / T(6)
@inline _kv56_small_cutoff(::Type{T}) where {T<:AbstractFloat} = T(0.1)
@inline _kv56_large_cutoff(::Type{T}) where {T<:AbstractFloat} = T(20)
@inline _kv56_integral_tmax(::Type{T}) where {T<:AbstractFloat} = T(9)
@inline _kv56_integral_bins() = 64

@inline function _kv56_small_coefficients(::Type{T}) where {T<:AbstractFloat}
    nu = _kv56_order(T)
    a0 = T(0.5) * T(gamma(T(5) / T(6))) * T(2)^nu
    b0 = T(0.5) * T(gamma(-T(5) / T(6))) * T(2)^(-nu)
    c1 = inv(T(4) * (one(T) - nu))
    c2 = inv(T(32) * (one(T) - nu) * (T(2) - nu))
    d1 = inv(T(4) * (one(T) + nu))
    d2 = inv(T(32) * (one(T) + nu) * (T(2) + nu))
    return (; nu, a0, b0, c1, c2, d1, d2)
end

@inline function _scaled_kv56_small(x::T, nu::T, a0::T, b0::T, c1::T, c2::T, d1::T, d2::T) where {T<:AbstractFloat}
    x2 = x * x
    lead = a0 * (one(T) + x2 * c1 + (x2 * x2) * c2)
    tail = b0 * x^(T(2) * nu) * (one(T) + x2 * d1 + (x2 * x2) * d2)
    return lead + tail
end

@inline function _scaled_kv56_integral(x::T, nu::T, dt::T, n_bins::Int) where {T<:AbstractFloat}
    acc = zero(T)
    @inbounds for k in 0:(n_bins - 1)
        t = (T(k) + T(0.5)) * dt
        acc += exp(-x * cosh(t)) * cosh(nu * t)
    end
    return x^nu * dt * acc
end

@inline function _scaled_kv56_asymptotic(x::T, nu::T) where {T<:AbstractFloat}
    xinv = inv(x)
    sum_terms = one(T) +
        (T(2) / T(9)) * xinv +
        (-T(7) / T(81)) * xinv^2 +
        (T(175) / T(2187)) * xinv^3 +
        (-T(2275) / T(19683)) * xinv^4 +
        (T(5005) / T(177147)) * xinv^5
    return x^nu * sqrt(T(pi) / (T(2) * x)) * exp(-x) * sum_terms
end

@inline function _scaled_kv56_scalar(x::T) where {T<:AbstractFloat}
    coeffs = _kv56_small_coefficients(T)
    if x < _kv56_small_cutoff(T)
        return _scaled_kv56_small(x, coeffs.nu, coeffs.a0, coeffs.b0, coeffs.c1, coeffs.c2, coeffs.d1, coeffs.d2)
    end
    if x >= _kv56_large_cutoff(T)
        return _scaled_kv56_asymptotic(x, coeffs.nu)
    end
    dt = _kv56_integral_tmax(T) / T(_kv56_integral_bins())
    return _scaled_kv56_integral(x, coeffs.nu, dt, _kv56_integral_bins())
end

@inline function _kv56_scalar(x::T) where {T<:AbstractFloat}
    x == zero(T) && return T(Inf)
    return _scaled_kv56_scalar(x) / x^_kv56_order(T)
end
