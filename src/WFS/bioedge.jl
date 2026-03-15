using Statistics

@kernel function edge_mask_kernel!(mask, pupil, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        if @inbounds pupil[i, j]
            neighbor = false
            @inbounds for di in -1:1, dj in -1:1
                ii = i + di
                jj = j + dj
                if ii < 1 || ii > n || jj < 1 || jj > n || !pupil[ii, jj]
                    neighbor = true
                end
            end
            mask[i, j] = neighbor
        else
            @inbounds mask[i, j] = false
        end
    end
end

@kernel function bioedge_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end

@kernel function bioedge_masks_kernel!(masks, one_c, zero_c, half::Int, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        left = j <= half
        top = i <= half
        @inbounds masks[i, j, 1] = left ? one_c : zero_c
        @inbounds masks[i, j, 2] = left ? zero_c : one_c
        @inbounds masks[i, j, 3] = top ? one_c : zero_c
        @inbounds masks[i, j, 4] = top ? zero_c : one_c
    end
end

@kernel function bin_edge_mask_kernel!(out, mask, binning::Int, n_out::Int, m_out::Int)
    i, j = @index(Global, NTuple)
    if i <= n_out && j <= m_out
        val = false
        @inbounds for ii in 1:binning, jj in 1:binning
            val |= mask[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        @inbounds out[i, j] = val
    end
end

struct BioEdgeParams{T<:AbstractFloat,M,N<:WFSNormalization}
    n_subap::Int
    threshold::T
    light_ratio::T
    normalization::N
    calib_modulation::T
    modulation::T
    modulation_points::Int
    extra_modulation_factor::Int
    delta_theta::T
    user_modulation_path::M
    grey_width::T
    grey_length::Union{Bool,T}
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
    modulation_phases::C3
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
    valid_i4q::A
    valid_signal::A
    signal_2d::R
    reference_signal_2d::R
    camera_frame::R
    binned_resolution::Int
    effective_resolution::Int
    nominal_detector_resolution::Int
    calibrated::Bool
    calibration_wavelength::T
end

struct BioEdgeWFS{M<:SensingMode,P<:BioEdgeParams,S<:BioEdgeState} <: AbstractWFS
    params::P
    state::S
end

function BioEdgeWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    light_ratio::Real=0.0, modulation::Real=0.0, modulation_points::Union{Int,Nothing}=nothing,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    extra_modulation_factor::Int=0, delta_theta::Real=0.0, user_modulation_path=nothing,
    grey_width::Real=0.0, grey_length=false,
    diffraction_padding::Int=2, psf_centering::Bool=true, n_pix_separation=nothing,
    n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    grey_length_val = grey_length === false ? false : T(grey_length)
    n_mod = resolve_modulation_points(T(modulation), modulation_points, extra_modulation_factor, user_modulation_path)
    params = BioEdgeParams{T,typeof(user_modulation_path),typeof(normalization)}(
        n_subap,
        T(threshold),
        T(light_ratio),
        normalization,
        T(calib_modulation),
        T(modulation),
        n_mod,
        extra_modulation_factor,
        T(delta_theta),
        user_modulation_path,
        T(grey_width),
        grey_length_val,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge,
        binning,
    )
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
    modulation_phases = backend{Complex{T}}(undef, tel.params.resolution, tel.params.resolution, n_mod)
    intensity = backend{T}(undef, 2 * pad, 2 * pad)
    temp = backend{T}(undef, pad, pad)
    scratch = similar(temp)
    binned_intensity = similar(intensity)
    nominal_detector_resolution = max(1, round(Int, n_subap * pad / tel.params.resolution))
    camera_frame = backend{T}(undef, 2 * nominal_detector_resolution, 2 * nominal_detector_resolution)
    valid_i4q = backend{Bool}(undef, nominal_detector_resolution ÷ 2, nominal_detector_resolution ÷ 2)
    valid_signal = backend{Bool}(undef, nominal_detector_resolution, nominal_detector_resolution ÷ 2)
    signal_2d = backend{T}(undef, nominal_detector_resolution, nominal_detector_resolution ÷ 2)
    reference_signal_2d = similar(signal_2d)
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
        modulation_phases,
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
        valid_i4q,
        valid_signal,
        signal_2d,
        reference_signal_2d,
        camera_frame,
        tel.params.resolution,
        pad,
        nominal_detector_resolution,
        false,
        zero(T),
    )
    wfs = BioEdgeWFS{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    update_edge_mask!(wfs, tel)
    build_bioedge_phasor!(wfs.state.phasor)
    build_bioedge_masks!(wfs)
    build_modulation_phases!(wfs, tel)
    return wfs
end

sensing_mode(::BioEdgeWFS{M}) where {M} = M()

function update_valid_mask!(wfs::BioEdgeWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function update_edge_mask!(wfs::BioEdgeWFS, tel::Telescope)
    Base.require_one_based_indexing(wfs.state.edge_mask, tel.state.pupil)
    _update_edge_mask!(execution_style(wfs.state.edge_mask), wfs.state.edge_mask, tel.state.pupil, tel.params.resolution)
    return wfs
end

function _update_edge_mask!(::ScalarCPUStyle, mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool}, n::Int)
    @inbounds for i in 1:n, j in 1:n
        if pupil[i, j]
            neighbor = false
            for di in -1:1, dj in -1:1
                ii = i + di
                jj = j + dj
                if ii < 1 || ii > n || jj < 1 || jj > n || !pupil[ii, jj]
                    neighbor = true
                end
            end
            mask[i, j] = neighbor
        else
            mask[i, j] = false
        end
    end
    return mask
end

function _update_edge_mask!(style::AcceleratorStyle, mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool}, n::Int)
    launch_kernel!(style, edge_mask_kernel!, mask, pupil, n; ndrange=size(mask))
    return mask
end

function build_bioedge_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_bioedge_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_bioedge_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_bioedge_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, bioedge_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_bioedge_masks!(wfs::BioEdgeWFS)
    masks = wfs.state.bioedge_masks
    copyto!(masks, host_bioedge_masks(wfs))
    return masks
end

function build_modulation_phases!(wfs::BioEdgeWFS, tel::Telescope)
    copyto!(wfs.state.modulation_phases, host_modulation_phases(eltype(wfs.state.slopes), tel,
        wfs.params.modulation, wfs.params.modulation_points, wfs.params.delta_theta,
        wfs.params.user_modulation_path))
    return wfs.state.modulation_phases
end

function host_bioedge_masks(wfs::BioEdgeWFS)
    T = eltype(wfs.state.slopes)
    n = size(wfs.state.bioedge_masks, 1)
    host = Array{Complex{T}}(undef, n, n, 4)
    build_bioedge_masks_host!(host, wfs)
    return host
end

function build_bioedge_masks_host!(masks::AbstractArray{Complex{T},3}, wfs::BioEdgeWFS) where {T<:AbstractFloat}
    n = size(masks, 1)
    half = n ÷ 2
    bw = zeros(T, n)
    bw[1:half] .= one(T)
    r = round(Int, wfs.params.diffraction_padding * wfs.params.grey_width)
    if r > 0
        gradient = vcat(segment_values(one(T), T(0.5), r), segment_values(T(0.5), zero(T), r))
        lo = max(1, half - r + 1)
        hi = min(n, half + r)
        bw[lo:hi] .= gradient[1:(hi - lo + 1)]
    end
    X = repeat(reshape(bw, 1, :), n, 1)
    A = sqrt.(X)
    if wfs.params.grey_length !== false
        r_grey = wfs.params.diffraction_padding
        r_length = round(Int, r_grey * wfs.params.grey_length)
        top_stop = max(1, half - r_length)
        bot_start = min(n + 1, half + r_length + 1)
        if top_stop >= 1
            A[1:top_stop, 1:half] .= one(T)
            A[1:top_stop, half+1:end] .= zero(T)
        end
        if bot_start <= n
            A[bot_start:end, 1:half] .= one(T)
            A[bot_start:end, half+1:end] .= zero(T)
        end
    end
    B = sqrt.(max.(zero(T), one(T) .- A .^ 2))
    C = permutedims(A)
    D = permutedims(B)
    @views begin
        masks[:, :, 1] .= complex.(A, zero(T))
        masks[:, :, 2] .= complex.(B, zero(T))
        masks[:, :, 3] .= complex.(C, zero(T))
        masks[:, :, 4] .= complex.(D, zero(T))
    end
    return masks
end

segment_values(a::T, b::T, n::Int) where {T<:AbstractFloat} = n == 1 ? T[a] : collect(range(a, b; length=n))

function _build_bioedge_masks!(::ScalarCPUStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
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

function _build_bioedge_masks!(style::AcceleratorStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
    n = size(masks, 1)
    half = n ÷ 2
    one_c = complex(one(T), zero(T))
    zero_c = complex(zero(T), zero(T))
    launch_kernel!(style, bioedge_masks_kernel!, masks, one_c, zero_c, half, n; ndrange=(n, n))
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
        wfs.state.calibrated = false
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
    _bin_edge_mask!(execution_style(out), out, mask, binning, n_out, m_out)
    return out
end

function _bin_edge_mask!(::ScalarCPUStyle, out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int, n_out::Int, m_out::Int)
    @inbounds for i in 1:n_out, j in 1:m_out
        val = false
        for ii in 1:binning, jj in 1:binning
            val |= mask[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = val
    end
    return out
end

function _bin_edge_mask!(style::AcceleratorStyle, out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int, n_out::Int, m_out::Int)
    launch_kernel!(style, bin_edge_mask_kernel!, out, mask, binning, n_out, m_out; ndrange=size(out))
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

function resize_bioedge_signal_buffers!(wfs::BioEdgeWFS, frame_rows::Int)
    nominal = wfs.state.nominal_detector_resolution
    n_pixels = max(1, round(Int, nominal / (2 * wfs.params.binning)))
    if size(wfs.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.state.valid_i4q = similar(wfs.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.state.valid_i4q, false)
    end
    if size(wfs.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.state.valid_signal = similar(wfs.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.signal_2d = similar(wfs.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    elseif size(wfs.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    end
    if size(wfs.state.camera_frame) != (frame_rows, frame_rows)
        wfs.state.camera_frame = similar(wfs.state.camera_frame, frame_rows, frame_rows)
    end
    update_bioedge_valid_signal!(wfs)
    return wfs
end

function sample_bioedge_intensity!(wfs::BioEdgeWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    sub = div(tel.params.resolution, wfs.params.n_subap)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("bioedge intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    if size(wfs.state.camera_frame) != (n_camera, n_camera)
        wfs.state.camera_frame = similar(wfs.state.camera_frame, n_camera, n_camera)
    end
    bin2d!(wfs.state.camera_frame, intensity, sub)
    frame = wfs.state.camera_frame
    wfs.state.nominal_detector_resolution = round(Int, wfs.params.n_subap * wfs.state.effective_resolution / tel.params.resolution)
    if wfs.params.binning != 1
        target = div(wfs.state.nominal_detector_resolution, wfs.params.binning)
        factor = div(size(frame, 1), target)
        if factor < 1 || size(frame, 1) % target != 0
            throw(InvalidConfiguration("bioedge detector binning is not compatible with the sampled frame"))
        end
        if size(wfs.state.binned_intensity) != (target, target)
            wfs.state.binned_intensity = similar(wfs.state.binned_intensity, target, target)
        end
        bin2d!(wfs.state.binned_intensity, frame, factor)
        frame = wfs.state.binned_intensity
    end
    resize_bioedge_signal_buffers!(wfs, size(frame, 1))
    return frame
end

function bioedge_intensity_core!(out::AbstractMatrix{T}, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource; apply_lgs::Bool=false) where {T<:AbstractFloat}
    prepare_bioedge_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    phase_scale = (2 * pi) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2 /
        wfs.params.modulation_points))

    fill!(out, zero(T))
    profile = apply_lgs ? lgs_profile(src) : LGSProfileNone()
    if apply_lgs && profile isa LGSProfileNaProfile
        ensure_lgs_kernel!(wfs, tel, src)
    end

    @inbounds for p in 1:wfs.params.modulation_points
        fill!(wfs.state.field, zero(eltype(wfs.state.field)))
        @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
            wfs.state.modulation_phases[:, :, p] * cis(phase_scale * tel.state.opd)
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
            @views out[oxq+1:oxq+pad, oyq+1:oyq+pad] .+= wfs.state.temp
        end
    end
    out ./= wfs.params.modulation_points
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
    edge_geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask, wfs.state.edge_mask)
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
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    return bioedge_signal!(wfs, tel, intensity, src)
end

function measure!(::Diffractive, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource)
    ensure_bioedge_calibration!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    return bioedge_signal!(wfs, tel, intensity, src)
end

function bioedge_slopes!(wfs::BioEdgeWFS, phase::AbstractMatrix, edge_mask::AbstractMatrix{Bool})
    edge_geometric_slopes!(wfs.state.slopes, phase, wfs.state.valid_mask, edge_mask)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function bioedge_slopes_intensity!(wfs::BioEdgeWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    return bioedge_signal!(wfs, tel, intensity)
end

function bioedge_signal!(wfs::BioEdgeWFS, tel::Telescope, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    return bioedge_signal!(wfs, tel, frame, nothing)
end

function bioedge_signal!(wfs::BioEdgeWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
    n_pixels = center
    count = 0
    norma = zero(T)
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        if wfs.state.valid_i4q[i, j]
            count += 1
            norma += q1 + q2 + q3 + q4
        end
    end
    norma = bioedge_normalization(wfs.params.normalization, wfs, tel, src, count, norma)
    idx = 1
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        sx = (q1 - q2 + q4 - q3) / norma
        sy = (q1 - q4 + q2 - q3) / norma
        wfs.state.signal_2d[i, j] = sx - wfs.state.reference_signal_2d[i, j]
        wfs.state.signal_2d[i + n_pixels, j] = sy - wfs.state.reference_signal_2d[i + n_pixels, j]
        if wfs.state.valid_i4q[i, j]
            wfs.state.slopes[idx] = wfs.state.signal_2d[i, j]
            wfs.state.slopes[idx + count] = wfs.state.signal_2d[i + n_pixels, j]
            idx += 1
        end
    end
    if idx <= count
        fill!(@view(wfs.state.slopes[idx:count]), zero(T))
        fill!(@view(wfs.state.slopes[count + idx:end]), zero(T))
    end
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function bioedge_normalization(::MeanValidFluxNormalization, ::BioEdgeWFS, ::Telescope,
    ::Union{Nothing,AbstractSource}, count::Int, summed_i4q)
    T = typeof(summed_i4q)
    return count == 0 ? one(T) : summed_i4q / count
end

function bioedge_normalization(::IncidenceFluxNormalization, wfs::BioEdgeWFS, tel::Telescope,
    src::AbstractSource, ::Int, summed_i4q)
    T = typeof(summed_i4q)
    sub_area = (tel.params.diameter / wfs.params.n_subap)^2
    return T(photon_flux(src) * tel.params.sampling_time * sub_area)
end

function bioedge_normalization(::IncidenceFluxNormalization, ::BioEdgeWFS, ::Telescope,
    ::Nothing, ::Int, summed_i4q)
    return one(typeof(summed_i4q))
end

function update_bioedge_valid_signal!(wfs::BioEdgeWFS)
    n_pixels = size(wfs.state.valid_i4q, 1)
    fill!(wfs.state.valid_signal, false)
    @views begin
        wfs.state.valid_signal[1:n_pixels, :] .= wfs.state.valid_i4q
        wfs.state.valid_signal[n_pixels+1:end, :] .= wfs.state.valid_i4q
    end
    return wfs
end

function select_bioedge_valid_i4q!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    n_pixels = max(1, round(Int, wfs.state.nominal_detector_resolution / (2 * wfs.params.binning)))
    if size(wfs.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.state.valid_i4q = similar(wfs.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.state.valid_i4q, false)
    end
    if size(wfs.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.state.valid_signal = similar(wfs.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.signal_2d = similar(wfs.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    elseif size(wfs.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    end
    if iszero(wfs.params.light_ratio)
        fill!(wfs.state.valid_i4q, true)
        update_bioedge_valid_signal!(wfs)
        return wfs
    end

    copyto!(wfs.state.modulation_phases, host_modulation_phases(
        eltype(wfs.state.slopes),
        tel,
        wfs.params.calib_modulation,
        wfs.params.modulation_points,
        wfs.params.delta_theta,
        wfs.params.user_modulation_path,
    ))
    bioedge_intensity!(wfs.state.temp, wfs, tel, src)
    frame = sample_bioedge_intensity!(wfs, tel, wfs.state.temp)
    build_modulation_phases!(wfs, tel)

    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
    max_i4q = zero(eltype(frame))
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        i4q = q1 + q2 + q3 + q4
        if i4q > max_i4q
            max_i4q = i4q
        end
    end
    cutoff = wfs.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        wfs.state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
    end
    update_bioedge_valid_signal!(wfs)
    return wfs
end

function ensure_bioedge_calibration!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    λ = eltype(wfs.state.slopes)(wavelength(src))
    if wfs.state.calibrated && wfs.state.calibration_wavelength == λ
        return wfs
    end
    opd_saved = copy(tel.state.opd)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    select_bioedge_valid_i4q!(wfs, tel, src)
    bioedge_intensity!(wfs.state.intensity, wfs, tel, src)
    frame = sample_bioedge_intensity!(wfs, tel, wfs.state.intensity)
    fill!(wfs.state.reference_signal_2d, zero(eltype(wfs.state.reference_signal_2d)))
    bioedge_signal!(wfs, tel, frame, src)
    copyto!(wfs.state.reference_signal_2d, wfs.state.signal_2d)
    fill!(wfs.state.slopes, zero(eltype(wfs.state.slopes)))
    copyto!(tel.state.opd, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    return wfs
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
