using Statistics

struct ShackHartmannParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
end

mutable struct ShackHartmannState{T<:AbstractFloat,A<:AbstractMatrix{Bool},V<:AbstractVector{T}}
    valid_mask::A
    slopes::V
end

struct ShackHartmann{M<:SensingMode,P<:ShackHartmannParams,S<:ShackHartmannState} <: AbstractWFS
    params::P
    state::S
end

function ShackHartmann(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    params = ShackHartmannParams{T}(n_subap, T(threshold))
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    state = ShackHartmannState{T, typeof(valid_mask), typeof(slopes)}(valid_mask, slopes)
    wfs = ShackHartmann{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ShackHartmann{M}) where {M} = M()

function update_valid_mask!(wfs::ShackHartmann, tel::Telescope)
    pupil = tel.state.pupil
    ax1 = axes(pupil, 1)
    ax2 = axes(pupil, 2)
    n = length(ax1)
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    x0 = firstindex(pupil, 1)
    y0 = firstindex(pupil, 2)
    x_last = lastindex(pupil, 1)
    y_last = lastindex(pupil, 2)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = x0 + (i - 1) * sub
        ys = y0 + (j - 1) * sub
        xe = min(x0 + i * sub - 1, x_last)
        ye = min(y0 + j * sub - 1, y_last)
        subap_pupil = @view pupil[xs:xe, ys:ye]
        wfs.state.valid_mask[i, j] = mean(subap_pupil) > wfs.params.threshold
    end
    return wfs
end

function measure!(mode::Geometric, wfs::ShackHartmann, tel::Telescope)
    opd = tel.state.opd
    ax1 = axes(opd, 1)
    ax2 = axes(opd, 2)
    n = length(ax1)
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    idx = 1
    x0 = firstindex(opd, 1)
    y0 = firstindex(opd, 2)
    x_last = lastindex(opd, 1)
    y_last = lastindex(opd, 2)

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = x0 + (i - 1) * sub
        ys = y0 + (j - 1) * sub
        xe = min(x0 + i * sub - 1, x_last)
        ye = min(y0 + j * sub - 1, y_last)
        if wfs.state.valid_mask[i, j]
            sx = 0.0
            sy = 0.0
            count_x = 0
            count_y = 0
            for x in xs:(xe - 1), y in ys:ye
                sx += opd[x + 1, y] - opd[x, y]
                count_x += 1
            end
            for x in xs:xe, y in ys:(ye - 1)
                sy += opd[x, y + 1] - opd[x, y]
                count_y += 1
            end
            wfs.state.slopes[idx] = sx / max(count_x, 1)
            wfs.state.slopes[idx + n_sub * n_sub] = sy / max(count_y, 1)
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end

    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope)
    return measure!(Geometric(), wfs, tel)
end

function measure!(wfs::ShackHartmann, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    slopes = measure!(sensing_mode(wfs), wfs, tel)
    n_sub = wfs.params.n_subap
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end
