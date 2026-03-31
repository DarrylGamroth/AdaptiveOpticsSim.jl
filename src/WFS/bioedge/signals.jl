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
    return bioedge_signal!(execution_style(frame), wfs, tel, frame, src)
end

function bioedge_signal!(::ScalarCPUStyle, wfs::BioEdgeWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
    n_pixels = center
    count = div(length(wfs.state.slopes), 2)
    norma = zero(T)
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        if wfs.state.valid_i4q[i, j]
            norma += q1 + q2 + q3 + q4
        end
    end
    norma = bioedge_normalization(wfs.params.normalization, wfs, tel, src, count, norma)
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
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
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function bioedge_signal!(::AcceleratorStyle, wfs::BioEdgeWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    count = wfs.state.valid_signal_count
    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
    n_pixels = center
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    sx = @view wfs.state.signal_2d[1:n_pixels, :]
    sy = @view wfs.state.signal_2d[n_pixels+1:2*n_pixels, :]
    refx = @view wfs.state.reference_signal_2d[1:n_pixels, :]
    refy = @view wfs.state.reference_signal_2d[n_pixels+1:2*n_pixels, :]
    i4q = @view wfs.state.temp[1:n_pixels, 1:n_pixels]
    @. i4q = q1 + q2 + q3 + q4
    summed_i4q = bioedge_valid_flux_sum!(execution_style(frame), wfs, i4q)
    norma = bioedge_normalization(wfs.params.normalization, wfs, tel, src, count, summed_i4q)
    @. sx = (q1 - q2 + q4 - q3) / norma - refx
    @. sy = (q1 - q4 + q2 - q3) / norma - refy
    launch_kernel!(execution_style(frame), gather_bioedge_slopes_kernel!, wfs.state.slopes,
        wfs.state.signal_2d, wfs.state.valid_signal_indices, count, n_pixels; ndrange=count)
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

function update_bioedge_valid_signal_indices!(wfs::BioEdgeWFS)
    valid_host = wfs.state.valid_i4q_host
    if size(valid_host) != size(wfs.state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(wfs.state.valid_i4q)...)
        wfs.state.valid_i4q_host = valid_host
    end
    copyto!(valid_host, wfs.state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(wfs.state.valid_signal_indices) < n_valid
        wfs.state.valid_signal_indices = similar(wfs.state.valid_signal_indices, n_valid)
    end
    if length(wfs.state.valid_signal_indices_host) < n_valid
        wfs.state.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = wfs.state.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(wfs.state.valid_signal_indices, 1, host_indices, 1, n_valid)
    wfs.state.valid_signal_count = n_valid
    return n_valid
end

function resize_bioedge_slope_buffers!(wfs::BioEdgeWFS)
    n_valid = wfs.state.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("bioedge valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(wfs.state.slopes) != n_slopes
        wfs.state.slopes = similar(wfs.state.slopes, n_slopes)
    end
    if length(wfs.state.optical_gain) != n_slopes
        wfs.state.optical_gain = similar(wfs.state.optical_gain, n_slopes)
        fill!(wfs.state.optical_gain, one(eltype(wfs.state.optical_gain)))
    end
    return wfs
end

function bioedge_valid_flux_sum!(::ScalarCPUStyle, wfs::BioEdgeWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    return masked_sum2d(ScalarCPUStyle(), i4q, wfs.state.valid_i4q_host)
end

function bioedge_valid_flux_sum!(style::AcceleratorStyle, wfs::BioEdgeWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.state.valid_i4q,
        wfs.state.valid_i4q_host,
        wfs.state.valid_flux_sum_buffer,
        wfs.state.valid_flux_sum_host,
        wfs.state.valid_flux_i4q_host,
    )
    wfs.state.valid_flux_i4q_host = host_parent
    return summed
end

function select_bioedge_valid_i4q!(wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
    return select_bioedge_valid_i4q!(execution_style(wfs.state.valid_i4q), wfs, tel, src)
end

function select_bioedge_valid_i4q!(::ScalarCPUStyle, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
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
        update_bioedge_valid_signal_indices!(wfs)
        resize_bioedge_slope_buffers!(wfs)
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
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
    return wfs
end

function select_bioedge_valid_i4q!(::AcceleratorStyle, wfs::BioEdgeWFS, tel::Telescope, src::AbstractSource)
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
        update_bioedge_valid_signal_indices!(wfs)
        resize_bioedge_slope_buffers!(wfs)
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
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view wfs.state.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.params.light_ratio * maximum(i4q)
    @. wfs.state.valid_i4q = i4q >= cutoff
    update_bioedge_valid_signal!(wfs)
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
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
