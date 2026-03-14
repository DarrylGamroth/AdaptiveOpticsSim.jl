using Statistics

struct BioEdgeParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
    binning::Int
end

mutable struct BioEdgeState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    SF<:SpatialFilter,
    C<:AbstractMatrix{Complex{T}},
    C3<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}}}
    valid_mask::A
    edge_mask::A
    slopes::V
    spatial_filter::SF
    field::C
    focal_field::C
    pupil_field::C
    bioedge_masks::C3
    phasor::C
    intensity::R
    temp::R
    scratch::R
    binned_intensity::R
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
    effective_resolution::Int
end

struct BioEdgeWFS{M<:SensingMode,P<:BioEdgeParams,S<:BioEdgeState} <: AbstractWFS
    params::P
    state::S
end

function BioEdgeWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    diffraction_padding::Int=2, psf_centering::Bool=true, n_pix_separation=nothing,
    n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    params = BioEdgeParams{T}(n_subap, T(threshold), diffraction_padding, psf_centering,
        n_pix_separation, n_pix_edge, binning)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    edge_mask = backend{Bool}(undef, size(tel.state.pupil))
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    sf = SpatialFilter(tel; shape=FoucaultFilter(), zero_padding=diffraction_padding, T=T, backend=backend)
    pad = tel.params.resolution * diffraction_padding
    if n_pix_separation !== nothing
        edge = n_pix_edge === nothing ? div(n_pix_separation, 2) : n_pix_edge
        pad = Int(round((n_subap * 2 + n_pix_separation + 2 * edge) * tel.params.resolution / n_subap))
    end
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    bioedge_masks = backend{Complex{T}}(undef, pad, pad, 4)
    phasor = similar(field)
    intensity = backend{T}(undef, 2 * pad, 2 * pad)
    temp = backend{T}(undef, pad, pad)
    scratch = similar(temp)
    binned_intensity = similar(intensity)
    fft_buffer = similar(field)
    binned_phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    edge_mask_binned = similar(edge_mask)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    state = BioEdgeState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(sf),
        typeof(field),
        typeof(bioedge_masks),
        typeof(intensity),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
    }(
        valid_mask,
        edge_mask,
        slopes,
        sf,
        field,
        focal_field,
        pupil_field,
        bioedge_masks,
        phasor,
        intensity,
        temp,
        scratch,
        binned_intensity,
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
        pad,
    )
    wfs = BioEdgeWFS{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    update_edge_mask!(wfs, tel)
    build_bioedge_phasor!(wfs.state.phasor)
    build_bioedge_masks!(wfs)
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

function build_bioedge_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function build_bioedge_masks!(wfs::BioEdgeWFS)
    masks = wfs.state.bioedge_masks
    T = eltype(wfs.state.slopes)
    one_c = complex(one(T), zero(T))
    zero_c = complex(zero(T), zero(T))
    n = size(masks, 1)
    half = n ÷ 2
    @inbounds for i in 1:n, j in 1:n
        left = j <= half
        top = i <= half
        masks[i, j, 1] = left ? one_c : zero_c
        masks[i, j, 2] = left ? zero_c : one_c
        masks[i, j, 3] = top ? one_c : zero_c
        masks[i, j, 4] = top ? zero_c : one_c
    end
    return masks
end

function ensure_bioedge_buffers!(wfs::BioEdgeWFS, pad::Int, tel::Telescope)
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.focal_field = similar(wfs.state.focal_field, pad, pad)
        wfs.state.pupil_field = similar(wfs.state.pupil_field, pad, pad)
        wfs.state.bioedge_masks = similar(wfs.state.bioedge_masks, pad, pad, 4)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.intensity = similar(wfs.state.intensity, 2 * pad, 2 * pad)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.scratch = similar(wfs.state.scratch, pad, pad)
        wfs.state.binned_intensity = similar(wfs.state.binned_intensity, 2 * pad, 2 * pad)
        wfs.state.fft_buffer = similar(wfs.state.fft_buffer, pad, pad)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.focal_field)
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.pupil_field)
        wfs.state.lgs_kernel_fft = similar(wfs.state.focal_field, Complex{eltype(wfs.state.focal_field)}, 0, 0)
        wfs.state.lgs_kernel_tag = UInt(0)
        wfs.state.effective_resolution = pad
        build_bioedge_phasor!(wfs.state.phasor)
        build_bioedge_masks!(wfs)
    end
    return wfs
end

function prepare_bioedge_sampling!(wfs::BioEdgeWFS, tel::Telescope)
    binning = wfs.params.binning
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_sub = wfs.params.n_subap
    pad = tel.params.resolution * wfs.params.diffraction_padding
    if wfs.params.n_pix_separation !== nothing
        edge = wfs.params.n_pix_edge === nothing ? div(wfs.params.n_pix_separation, 2) : wfs.params.n_pix_edge
        pad = Int(round((n_sub * 2 + wfs.params.n_pix_separation + 2 * edge) * tel.params.resolution / n_sub))
    end
    if pad < tel.params.resolution
        throw(InvalidConfiguration("bioedge padding must be >= telescope resolution"))
    end
    if pad % binning != 0
        throw(InvalidConfiguration("bioedge binning must evenly divide padded resolution"))
    end
    n = tel.params.resolution
    if n % binning != 0
        throw(InvalidConfiguration("binning must evenly divide telescope resolution"))
    end
    ensure_bioedge_buffers!(wfs, pad, tel)
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

function sample_bioedge_intensity!(wfs::BioEdgeWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.params.binning
    if binning == 1
        return intensity
    end
    pad = size(intensity, 1)
    if pad % binning != 0
        throw(InvalidConfiguration("bioedge binning must evenly divide padded resolution"))
    end
    n_binned = div(pad, binning)
    if size(wfs.state.binned_intensity) != (n_binned, n_binned)
        wfs.state.binned_intensity = similar(wfs.state.binned_intensity, n_binned, n_binned)
    end
    bin2d!(wfs.state.binned_intensity, intensity, binning)
    return wfs.state.binned_intensity
end

function bioedge_intensity_core!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource; apply_lgs::Bool=false) where {T<:AbstractFloat}
    prepare_bioedge_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    phase_scale = (2 * pi) / wavelength(src)

    fill!(out, zero(T))
    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = tel.state.pupil * cis(phase_scale * tel.state.opd)
    copyto!(wfs.state.focal_field, wfs.state.field)
    if wfs.params.psf_centering
        @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
        mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
    else
        mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
        fftshift2d!(wfs.state.fft_buffer, wfs.state.focal_field)
    end

    focal_source = wfs.params.psf_centering ? wfs.state.focal_field : wfs.state.fft_buffer
    lgs_fft_buffer = wfs.params.psf_centering ? wfs.state.fft_buffer : wfs.state.focal_field
    lgs_ifft_buffer = wfs.state.pupil_field
    profile = apply_lgs ? lgs_profile(src) : LGSProfileNone()
    if apply_lgs && profile isa LGSProfileNaProfile
        ensure_lgs_kernel!(wfs, tel, src)
    end

    @inbounds for k in 1:4
        @views @. wfs.state.pupil_field = focal_source * wfs.state.bioedge_masks[:, :, k]
        mul!(wfs.state.pupil_field, wfs.state.ifft_plan, wfs.state.pupil_field)
        @. wfs.state.temp = abs2(wfs.state.pupil_field)
        if apply_lgs
            if profile isa LGSProfileNone
                wfs.state.elongation_kernel = apply_elongation!(
                    wfs.state.temp,
                    lgs_elongation_factor(src),
                    wfs.state.scratch,
                    wfs.state.elongation_kernel,
                )
            else
                apply_lgs_convolution!(
                    wfs.state.temp,
                    wfs.state.lgs_kernel_fft,
                    lgs_fft_buffer,
                    wfs.state.fft_plan,
                    lgs_ifft_buffer,
                    wfs.state.ifft_plan,
                )
            end
        end
        oxq = k in (3, 4) ? pad : 0
        oyq = k in (2, 4) ? pad : 0
        @views out[oxq+1:oxq+pad, oyq+1:oyq+pad] .= wfs.state.temp
    end
    return out
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, tel, src; apply_lgs=false)
end

function bioedge_intensity!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::LGSSource) where {T<:AbstractFloat}
    return bioedge_intensity_core!(out, wfs, tel, src; apply_lgs=true)
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
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, wfs.state.intensity)
    return bioedge_slopes_intensity!(wfs, tel, intensity)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, wfs.state.intensity)
    return bioedge_slopes_intensity!(wfs, tel, intensity)
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

function bioedge_slopes_intensity!(wfs::BioEdgeWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    pad = div(size(intensity, 1), 2)
    n_sub = wfs.params.n_subap
    binning = wfs.params.binning
    pupil_size = div(tel.params.resolution, binning)
    edge_pix = div(pad - pupil_size, 2)
    if edge_pix < 0
        throw(InvalidConfiguration("bioedge padding is too small for pupil extraction"))
    end
    sub = div(pupil_size, n_sub)
    if sub < 1
        throw(InvalidConfiguration("bioedge pupil size is too small for n_subap"))
    end
    ox1 = edge_pix
    oy1 = edge_pix
    ox2 = edge_pix
    oy2 = edge_pix + pad
    ox3 = edge_pix + pad
    oy3 = edge_pix
    ox4 = edge_pix + pad
    oy4 = edge_pix + pad
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if wfs.state.valid_mask[i, j]
            q1 = quadrant_sum(intensity, xs, ys, sub, ox1, oy1, 0, 0)
            q2 = quadrant_sum(intensity, xs, ys, sub, ox2, oy2, 0, 0)
            q3 = quadrant_sum(intensity, xs, ys, sub, ox3, oy3, 0, 0)
            q4 = quadrant_sum(intensity, xs, ys, sub, ox4, oy4, 0, 0)
            left = q1 + q3
            right = q2 + q4
            bottom = q1 + q2
            top = q3 + q4
            total = left + right
            if total <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                wfs.state.slopes[idx] = (right - left) / total
                wfs.state.slopes[idx + n_sub * n_sub] = (top - bottom) / total
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function apply_lgs_elongation!(::LGSProfileNone, intensity::AbstractMatrix{T}, wfs::BioEdgeWFS,
    ::Telescope, src::LGSSource) where {T<:AbstractFloat}
    wfs.state.elongation_kernel = apply_elongation!(
        intensity,
        lgs_elongation_factor(src),
        wfs.state.scratch,
        wfs.state.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, intensity::AbstractMatrix{T}, wfs::BioEdgeWFS,
    tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, tel, src)
    apply_lgs_convolution!(
        intensity,
        wfs.state.lgs_kernel_fft,
        wfs.state.fft_buffer,
        wfs.state.fft_plan,
        wfs.state.pupil_field,
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
    padding = wfs.state.effective_resolution / tel.params.resolution
    pixel_scale = lgs_pixel_scale(tel.params.diameter, padding, wavelength(src))
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
