using Statistics

struct ShackHartmannParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
end

mutable struct ShackHartmannState{T<:AbstractFloat,A<:AbstractMatrix{Bool},V<:AbstractVector{T}}
    valid_mask::A
    slopes::V
end

struct ShackHartmann{P<:ShackHartmannParams,S<:ShackHartmannState} <: AbstractWFS
    params::P
    state::S
end

function ShackHartmann(tel::Telescope; n_subap::Int, threshold::Real=0.1, T::Type{<:AbstractFloat}=Float64, backend=Array)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    params = ShackHartmannParams{T}(n_subap, T(threshold))
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    state = ShackHartmannState{T, typeof(valid_mask), typeof(slopes)}(valid_mask, slopes)
    wfs = ShackHartmann(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

function update_valid_mask!(wfs::ShackHartmann, tel::Telescope)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        subap_pupil = @view tel.state.pupil[xs:xe, ys:ye]
        wfs.state.valid_mask[i, j] = mean(subap_pupil) > wfs.params.threshold
    end
    return wfs
end

function measure!(wfs::ShackHartmann, tel::Telescope)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            sx = 0.0
            sy = 0.0
            count_x = 0
            count_y = 0
            @views subap = tel.state.opd[xs:xe, ys:ye]
            for x in 1:(size(subap, 1) - 1), y in 1:size(subap, 2)
                sx += subap[x + 1, y] - subap[x, y]
                count_x += 1
            end
            for x in 1:size(subap, 1), y in 1:(size(subap, 2) - 1)
                sy += subap[x, y + 1] - subap[x, y]
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
