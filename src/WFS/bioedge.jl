using Statistics

struct BioEdgeParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    diffraction_padding::Int
    binning::Int
end

mutable struct BioEdgeState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    SF<:SpatialFilter,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}}}
    valid_mask::A
    edge_mask::A
    slopes::V
    spatial_filter::SF
    fft_buffer::C
    binned_phase::R
    edge_mask_binned::A
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::Kf
    lgs_kernel_tag::UInt
    optical_gain::V
    binned_resolution::Int
end

struct BioEdgeWFS{M<:SensingMode,P<:BioEdgeParams,S<:BioEdgeState} <: AbstractWFS
    params::P
    state::S
end

function BioEdgeWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    diffraction_padding::Int=2, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    params = BioEdgeParams{T}(n_subap, T(threshold), diffraction_padding, binning)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    edge_mask = backend{Bool}(undef, size(tel.state.pupil))
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    sf = SpatialFilter(tel; shape=FoucaultFilter(), zero_padding=diffraction_padding, T=T, backend=backend)
    fft_buffer = backend{Complex{T}}(undef, tel.params.resolution, tel.params.resolution)
    binned_phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    edge_mask_binned = similar(edge_mask)
    fft_plan = plan_fft!(fft_buffer)
    ifft_plan = plan_ifft!(fft_buffer)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    state = BioEdgeState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(sf),
        typeof(fft_buffer),
        typeof(binned_phase),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
    }(
        valid_mask,
        edge_mask,
        slopes,
        sf,
        fft_buffer,
        binned_phase,
        edge_mask_binned,
        fft_plan,
        ifft_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
        optical_gain,
        tel.params.resolution,
    )
    wfs = BioEdgeWFS{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    update_edge_mask!(wfs, tel)
    return wfs
end

sensing_mode(::BioEdgeWFS{M}) where {M} = M()

function update_valid_mask!(wfs::BioEdgeWFS, tel::Telescope)
    Base.require_one_based_indexing(tel.state.pupil)
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
    Base.require_one_based_indexing(wfs.state.edge_mask, tel.state.pupil)
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

function prepare_bioedge_sampling!(wfs::BioEdgeWFS, tel::Telescope)
    binning = wfs.params.binning
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n = tel.params.resolution
    if n % binning != 0
        throw(InvalidConfiguration("binning must evenly divide telescope resolution"))
    end
    n_binned = div(n, binning)
    if n_binned != wfs.state.binned_resolution
        wfs.state.binned_phase = similar(wfs.state.binned_phase, n_binned, n_binned)
        wfs.state.edge_mask_binned = similar(wfs.state.edge_mask_binned, n_binned, n_binned)
        wfs.state.binned_resolution = n_binned
    end
    if binning > 1
        bin_edge_mask!(wfs.state.edge_mask_binned, wfs.state.edge_mask, binning)
    end
    return wfs
end

function bin_edge_mask!(out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int)
    Base.require_one_based_indexing(out, mask)
    n, m = size(mask)
    n_out = div(n, binning)
    m_out = div(m, binning)
    if size(out) != (n_out, m_out)
        throw(DimensionMismatchError("edge_mask_binned size mismatch"))
    end
    @inbounds for i in 1:n_out, j in 1:m_out
        val = false
        for ii in 1:binning, jj in 1:binning
            val |= mask[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = val
    end
    return out
end

function sample_bioedge_phase!(wfs::BioEdgeWFS, phase::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.params.binning
    if binning == 1
        return phase, wfs.state.edge_mask
    end
    bin2d!(wfs.state.binned_phase, phase, binning)
    wfs.state.binned_phase ./= binning * binning
    return wfs.state.binned_phase, wfs.state.edge_mask_binned
end

function measure!(mode::Geometric, wfs::BioEdgeWFS, tel::Telescope)
    Base.require_one_based_indexing(tel.state.opd)
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
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.n_subap
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive BioEdgeWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::BioEdgeWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    phase, _ = filter!(wfs.state.spatial_filter, tel, src)
    prepare_bioedge_sampling!(wfs, tel)
    phase_sampled, mask = sample_bioedge_phase!(wfs, phase)
    return bioedge_slopes!(wfs, phase_sampled, mask)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    phase, _ = filter!(wfs.state.spatial_filter, tel, src)
    apply_lgs_elongation!(lgs_profile(src), phase, wfs, tel, src)
    prepare_bioedge_sampling!(wfs, tel)
    phase_sampled, mask = sample_bioedge_phase!(wfs, phase)
    return bioedge_slopes!(wfs, phase_sampled, mask)
end

function bioedge_slopes!(wfs::BioEdgeWFS, phase::AbstractMatrix, edge_mask::AbstractMatrix{Bool})
    Base.require_one_based_indexing(phase, edge_mask)
    n = size(phase, 1)
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
                if edge_mask[x, y]
                    sx += phase[x + 1, y] - phase[x, y]
                    count_x += 1
                end
            end
            for x in xs:xe, y in ys:(ye - 1)
                if edge_mask[x, y]
                    sy += phase[x, y + 1] - phase[x, y]
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
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function apply_lgs_elongation!(::LGSProfileNone, phase::AbstractMatrix{T}, wfs::BioEdgeWFS,
    ::Telescope, src::LGSSource) where {T<:AbstractFloat}
    tmp = wfs.state.spatial_filter.state.amplitude
    wfs.state.elongation_kernel = apply_elongation!(
        phase,
        lgs_elongation_factor(src),
        tmp,
        wfs.state.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, phase::AbstractMatrix{T}, wfs::BioEdgeWFS,
    tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, tel, src)
    apply_lgs_convolution!(
        phase,
        wfs.state.lgs_kernel_fft,
        wfs.state.fft_buffer,
        wfs.state.fft_plan,
        wfs.state.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    if !(wfs.state.fft_buffer isa Array)
        throw(InvalidConfiguration("LGS Na-profile kernels currently require Array backend"))
    end
    pad = size(wfs.state.fft_buffer, 1)
    tag = objectid(na_profile) ⊻ hash(src.params.laser_coordinates) ⊻ hash(src.params.fwhm_spot_up) ⊻
        hash(pad) ⊻ hash(wfs.params.n_subap)
    if size(wfs.state.lgs_kernel_fft, 1) == pad && wfs.state.lgs_kernel_tag == tag
        return wfs
    end
    pixel_scale = lgs_pixel_scale(tel.params.diameter, wfs.params.diffraction_padding, wavelength(src))
    wfs.state.lgs_kernel_fft = lgs_average_kernel_fft(
        tel,
        src,
        pad,
        wfs.params.n_subap,
        pixel_scale,
        wfs.state.fft_buffer,
        wfs.state.fft_plan,
    )
    wfs.state.lgs_kernel_tag = tag
    return wfs
end

function set_optical_gain!(wfs::BioEdgeWFS, gain::Real)
    fill!(wfs.state.optical_gain, gain)
    return wfs
end

function set_optical_gain!(wfs::BioEdgeWFS, gain::AbstractVector)
    if length(gain) != length(wfs.state.optical_gain)
        throw(InvalidConfiguration("optical_gain length must match slope vector"))
    end
    copyto!(wfs.state.optical_gain, gain)
    return wfs
end
