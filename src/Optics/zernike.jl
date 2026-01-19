struct ZernikeBasis{T<:AbstractFloat}
    n_modes::Int
    modes::Array{T,3}
end

function ZernikeBasis(tel::Telescope, n_modes::Int; T::Type{<:AbstractFloat}=Float64)
    n = tel.params.resolution
    modes = zeros(T, n, n, n_modes)
    return ZernikeBasis{T}(n_modes, modes)
end

function noll_to_nm(j::Int)
    if j < 1
        throw(InvalidConfiguration("Noll index must be >= 1"))
    end
    count = 0
    n = 0
    while true
        for m in -n:2:n
            count += 1
            if count == j
                return n, m
            end
        end
        n += 1
    end
end

function zernike_radial(n::Int, m::Int, r::Real)
    m = abs(m)
    if isodd(n - m)
        return zero(r)
    end
    radial = zero(r)
    for k in 0:((n - m) ÷ 2)
        num = (-1)^k * factorial(n - k)
        den = factorial(k) * factorial((n + m) ÷ 2 - k) * factorial((n - m) ÷ 2 - k)
        radial += (num / den) * r^(n - 2k)
    end
    return radial
end

function compute_zernike!(zb::ZernikeBasis, tel::Telescope)
    ax1 = axes(zb.modes, 1)
    ax2 = axes(zb.modes, 2)
    if axes(tel.state.pupil, 1) != ax1 || axes(tel.state.pupil, 2) != ax2
        throw(DimensionMismatchError("pupil axes must match Zernike modes"))
    end
    n = length(ax1)
    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2
    i0 = first(ax1)
    j0 = first(ax2)
    pupil = tel.state.pupil
    modes = zb.modes

    for j in 1:zb.n_modes
        n_mode, m_mode = noll_to_nm(j)
        @inbounds for i in ax1, k in ax2
            if pupil[i, k]
                ii = i - i0 + 1
                kk = k - j0 + 1
                x = (ii - cx) / scale
                y = (kk - cy) / scale
                r = sqrt(x^2 + y^2)
                theta = atan(y, x)
                radial = zernike_radial(n_mode, m_mode, r)
                if m_mode >= 0
                    modes[i, k, j] = radial * cos(m_mode * theta)
                else
                    modes[i, k, j] = radial * sin(abs(m_mode) * theta)
                end
            else
                modes[i, k, j] = zero(eltype(modes))
            end
        end
    end
    return zb
end
