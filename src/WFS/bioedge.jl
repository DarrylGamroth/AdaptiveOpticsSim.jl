using Statistics

struct BioEdgeParams{T<:AbstractFloat,M<:SensingMode}
    n_subap::Int
    threshold::T
    mode::M
end

mutable struct BioEdgeState{T<:AbstractFloat,A<:AbstractMatrix{Bool},V<:AbstractVector{T}}
    valid_mask::A
    edge_mask::A
    slopes::V
end

struct BioEdgeWFS{P<:BioEdgeParams,S<:BioEdgeState} <: AbstractWFS
    params::P
    state::S
end

function BioEdgeWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    params = BioEdgeParams{T, typeof(mode)}(n_subap, T(threshold), mode)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    edge_mask = backend{Bool}(undef, size(tel.state.pupil))
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    state = BioEdgeState{T, typeof(valid_mask), typeof(slopes)}(valid_mask, edge_mask, slopes)
    wfs = BioEdgeWFS(params, state)
    update_valid_mask!(wfs, tel)
    update_edge_mask!(wfs, tel)
    return wfs
end

sensing_mode(wfs::BioEdgeWFS) = wfs.params.mode

function update_valid_mask!(wfs::BioEdgeWFS, tel::Telescope)
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

function update_edge_mask!(wfs::BioEdgeWFS, tel::Telescope)
    n = tel.params.resolution
    mask = wfs.state.edge_mask
    @inbounds for i in 1:n, j in 1:n
        if tel.state.pupil[i, j]
            neighbor = false
            for di in -1:1, dj in -1:1
                ii = i + di
                jj = j + dj
                if ii < 1 || ii > n || jj < 1 || jj > n || !tel.state.pupil[ii, jj]
                    neighbor = true
                end
            end
            mask[i, j] = neighbor
        else
            mask[i, j] = false
        end
    end
    return wfs
end

function measure!(mode::Geometric, wfs::BioEdgeWFS, tel::Telescope)
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
            for x in xs:(xe - 1), y in ys:ye
                if wfs.state.edge_mask[x, y]
                    sx += tel.state.opd[x + 1, y] - tel.state.opd[x, y]
                    count_x += 1
                end
            end
            for x in xs:xe, y in ys:(ye - 1)
                if wfs.state.edge_mask[x, y]
                    sy += tel.state.opd[x, y + 1] - tel.state.opd[x, y]
                    count_y += 1
                end
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

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope)
    return measure!(Geometric(), wfs, tel)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, mode::SensingMode=sensing_mode(wfs))
    return measure!(mode, wfs, tel)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource, mode::SensingMode=sensing_mode(wfs))
    slopes = measure!(mode, wfs, tel)
    if is_lgs(src)
        n_sub = wfs.params.n_subap
        factor = lgs_elongation_factor(src)
        @views slopes[n_sub * n_sub + 1:end] .*= factor
    end
    return slopes
end
