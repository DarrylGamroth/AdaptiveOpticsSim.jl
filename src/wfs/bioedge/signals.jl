function bioedge_slopes!(wfs::BioEdgeWFS, phase::AbstractMatrix, edge_mask::AbstractMatrix{Bool})
    edge_geometric_slopes!(wfs.estimator.state.slopes, phase, wfs.estimator.state.valid_mask, edge_mask)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function bioedge_slopes_intensity!(wfs::BioEdgeWFS, pupil::PupilFunction,
    intensity::AbstractMatrix{F}) where {F<:Real}
    return bioedge_signal!(wfs, pupil, intensity)
end

function bioedge_signal!(wfs::BioEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F}) where {F<:Real}
    return bioedge_signal!(wfs, pupil, frame, nothing)
end

function bioedge_signal!(wfs::BioEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(wfs.estimator.state.slopes)
    return bioedge_signal!(wfs, pupil, frame, src, one(S))
end

function bioedge_signal!(wfs::BioEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::Real) where {F<:Real}
    S = eltype(wfs.estimator.state.slopes)
    return bioedge_signal!(execution_style(frame), wfs, pupil, frame, src,
        S(normalization_scale))
end

function bioedge_signal!(::ScalarCPUStyle, wfs::BioEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(wfs.estimator.state.slopes)
    return bioedge_signal!(ScalarCPUStyle(), wfs, pupil, frame, src, one(S))
end

function bioedge_signal!(::ScalarCPUStyle, wfs::BioEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    return bioedge_signal!(ScalarCPUStyle(), wfs, frame, src,
        normalization_scale)
end

function bioedge_signal!(::ScalarCPUStyle, wfs::BioEdgeWFS,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    n_pixels = size(wfs.estimator.state.signal_2d, 2)
    center = require_bioedge_frame_geometry(wfs, frame)
    count = div(length(wfs.estimator.state.slopes), 2)
    norma = zero(S)
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = S(frame[center - n_pixels + i, center - n_pixels + j])
        q2 = S(frame[center - n_pixels + i, center + j])
        q3 = S(frame[center + i, center + j])
        q4 = S(frame[center + i, center - n_pixels + j])
        if wfs.estimator.state.valid_i4q[i, j]
            norma += q1 + q2 + q3 + q4
        end
    end
    norma = bioedge_normalization(wfs.estimator.params.normalization, wfs, src,
        count, norma, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(wfs.estimator.state.signal_2d, zero(S))
        fill!(wfs.estimator.state.slopes, zero(S))
        return wfs.estimator.state.slopes
    end
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        q1 = S(frame[center - n_pixels + i, center - n_pixels + j])
        q2 = S(frame[center - n_pixels + i, center + j])
        q3 = S(frame[center + i, center + j])
        q4 = S(frame[center + i, center - n_pixels + j])
        sx = (q1 - q2 + q4 - q3) / norma
        sy = (q1 - q4 + q2 - q3) / norma
        wfs.estimator.state.signal_2d[i, j] = sx - wfs.estimator.state.reference_signal_2d[i, j]
        wfs.estimator.state.signal_2d[i + n_pixels, j] = sy - wfs.estimator.state.reference_signal_2d[i + n_pixels, j]
        if wfs.estimator.state.valid_i4q[i, j]
            wfs.estimator.state.slopes[idx] = wfs.estimator.state.signal_2d[i, j]
            wfs.estimator.state.slopes[idx + count] = wfs.estimator.state.signal_2d[i + n_pixels, j]
            idx += 1
        end
    end
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function bioedge_signal!(style::AcceleratorStyle, wfs::BioEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(wfs.estimator.state.slopes)
    return bioedge_signal!(style, wfs, pupil, frame, src, one(S))
end

function bioedge_signal!(style::AcceleratorStyle, wfs::BioEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    return bioedge_signal!(style, wfs, frame, src, normalization_scale)
end

function bioedge_signal!(style::AcceleratorStyle, wfs::BioEdgeWFS,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    count = wfs.estimator.state.valid_signal_count
    n_pixels = size(wfs.estimator.state.signal_2d, 2)
    center = require_bioedge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    sx = @view wfs.estimator.state.signal_2d[1:n_pixels, :]
    sy = @view wfs.estimator.state.signal_2d[n_pixels+1:2*n_pixels, :]
    refx = @view wfs.estimator.state.reference_signal_2d[1:n_pixels, :]
    refy = @view wfs.estimator.state.reference_signal_2d[n_pixels+1:2*n_pixels, :]
    i4q = @view wfs.estimator.state.flux_i4q[1:n_pixels, 1:n_pixels]
    @. i4q = S(q1) + S(q2) + S(q3) + S(q4)
    summed_i4q = bioedge_valid_flux_sum!(style, wfs, i4q)
    norma = bioedge_normalization(wfs.estimator.params.normalization, wfs, src,
        count, summed_i4q, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(wfs.estimator.state.signal_2d, zero(S))
        fill!(wfs.estimator.state.slopes, zero(S))
        return wfs.estimator.state.slopes
    end
    @. sx = (S(q1) - S(q2) + S(q4) - S(q3)) / norma - refx
    @. sy = (S(q1) - S(q4) + S(q2) - S(q3)) / norma - refy
    launch_kernel!(style, gather_bioedge_slopes_kernel!, wfs.estimator.state.slopes,
        wfs.estimator.state.signal_2d, wfs.estimator.state.valid_signal_indices, count, n_pixels; ndrange=count)
    @. wfs.estimator.state.slopes *= wfs.estimator.state.optical_gain
    return wfs.estimator.state.slopes
end

function bioedge_normalization(normalization::MeanValidFluxNormalization,
    wfs::BioEdgeWFS, pupil::PupilFunction, src::Union{Nothing,AbstractSource},
    count::Int, summed_i4q)
    return bioedge_normalization(normalization, wfs, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function bioedge_normalization(normalization::MeanValidFluxNormalization,
    wfs::BioEdgeWFS, ::PupilFunction, src::Union{Nothing,AbstractSource},
    count::Int, summed_i4q, normalization_scale::Real)
    return bioedge_normalization(normalization, wfs, src, count,
        summed_i4q, normalization_scale)
end

function bioedge_normalization(::MeanValidFluxNormalization, ::BioEdgeWFS,
    ::Union{Nothing,AbstractSource}, count::Int, summed_i4q,
    ::Real)
    T = typeof(summed_i4q)
    return count == 0 ? one(T) : summed_i4q / count
end

function bioedge_normalization(normalization::IncidenceFluxNormalization,
    wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource, count::Int,
    summed_i4q)
    return bioedge_normalization(normalization, wfs, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function bioedge_normalization(normalization::IncidenceFluxNormalization,
    wfs::BioEdgeWFS, ::PupilFunction, src::AbstractSource, count::Int,
    summed_i4q, normalization_scale::Real)
    return bioedge_normalization(normalization, wfs, src, count, summed_i4q,
        normalization_scale)
end

function bioedge_normalization(::IncidenceFluxNormalization,
    wfs::BioEdgeWFS, src::AbstractSource, ::Int,
    summed_i4q, normalization_scale::Real)
    T = typeof(summed_i4q)
    sub_area = (wfs.estimator.params.pupil_diameter_m /
        wfs.estimator.params.pupil_samples)^2
    return wfs_incident_photon_irradiance(src, T) * T(sub_area) *
        T(normalization_scale)
end

function bioedge_normalization(::IncidenceFluxNormalization, ::BioEdgeWFS, ::PupilFunction,
    ::Nothing, ::Int, summed_i4q)
    return one(typeof(summed_i4q))
end


function bioedge_normalization(::IncidenceFluxNormalization, ::BioEdgeWFS,
    ::PupilFunction, ::Nothing, ::Int, summed_i4q, ::Real)
    return one(typeof(summed_i4q))
end

function bioedge_normalization(::IncidenceFluxNormalization, ::BioEdgeWFS,
    ::Nothing, ::Int, summed_i4q, ::Real)
    return one(typeof(summed_i4q))
end

function update_bioedge_valid_signal!(wfs::BioEdgeWFS)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    fill!(wfs.estimator.state.valid_signal, false)
    @views begin
        wfs.estimator.state.valid_signal[1:n_pixels, :] .= wfs.estimator.state.valid_i4q
        wfs.estimator.state.valid_signal[n_pixels+1:end, :] .= wfs.estimator.state.valid_i4q
    end
    return wfs
end

function update_bioedge_valid_signal_indices!(wfs::BioEdgeWFS)
    valid_host = wfs.estimator.state.valid_i4q_host
    if size(valid_host) != size(wfs.estimator.state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(wfs.estimator.state.valid_i4q)...)
        wfs.estimator.state.valid_i4q_host = valid_host
    end
    copyto!(valid_host, wfs.estimator.state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(wfs.estimator.state.valid_signal_indices) < n_valid
        wfs.estimator.state.valid_signal_indices = similar(wfs.estimator.state.valid_signal_indices, n_valid)
    end
    if length(wfs.estimator.state.valid_signal_indices_host) < n_valid
        wfs.estimator.state.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = wfs.estimator.state.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(wfs.estimator.state.valid_signal_indices, 1, host_indices, 1, n_valid)
    wfs.estimator.state.valid_signal_count = n_valid
    return n_valid
end

function resize_bioedge_slope_buffers!(wfs::BioEdgeWFS)
    n_valid = wfs.estimator.state.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("bioedge valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(wfs.estimator.state.slopes) != n_slopes
        wfs.estimator.state.slopes = similar(wfs.estimator.state.slopes, n_slopes)
    end
    if length(wfs.estimator.state.optical_gain) != n_slopes
        wfs.estimator.state.optical_gain = similar(wfs.estimator.state.optical_gain, n_slopes)
        fill!(wfs.estimator.state.optical_gain, one(eltype(wfs.estimator.state.optical_gain)))
    end
    return wfs
end

function bioedge_valid_flux_sum!(::ScalarCPUStyle, wfs::BioEdgeWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    return masked_sum2d(ScalarCPUStyle(), i4q, wfs.estimator.state.valid_i4q_host)
end

function bioedge_valid_flux_sum!(style::AcceleratorStyle, wfs::BioEdgeWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.estimator.state.valid_i4q,
        wfs.estimator.state.valid_i4q_host,
        wfs.estimator.state.valid_flux_sum_buffer,
        wfs.estimator.state.valid_flux_sum_host,
        wfs.estimator.state.valid_flux_i4q_host,
    )
    wfs.estimator.state.valid_flux_i4q_host = host_parent
    return summed
end

function select_bioedge_valid_i4q_from_frame!(::ScalarCPUStyle,
    wfs::BioEdgeWFS, frame::AbstractMatrix)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center = require_bioedge_frame_geometry(wfs, frame)
    max_i4q = zero(eltype(frame))
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        max_i4q = max(max_i4q, q1 + q2 + q3 + q4)
    end
    cutoff = wfs.estimator.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        wfs.estimator.state.valid_i4q[i, j] =
            (q1 + q2 + q3 + q4) >= cutoff
    end
    update_bioedge_valid_signal!(wfs)
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
    return wfs
end

function select_bioedge_valid_i4q_from_frame!(::AcceleratorStyle,
    wfs::BioEdgeWFS, frame::AbstractMatrix)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center = require_bioedge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view wfs.estimator.state.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. wfs.estimator.state.valid_i4q = i4q >= cutoff
    update_bioedge_valid_signal!(wfs)
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
    return wfs
end

function select_bioedge_valid_i4q!(wfs::BioEdgeWFS, pupil::PupilFunction,
    src::AbstractSource, det::Detector)
    bioedge_intensity_core!(wfs.front_end.propagation.intensity, wfs, pupil,
        src, bioedge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    sampled = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
    frame = detector_calibration_frame!(det, sampled, src)
    resize_bioedge_signal_buffers!(wfs, size(frame, 1), det)
    return select_bioedge_valid_i4q_from_frame!(execution_style(frame), wfs,
        frame)
end

function select_bioedge_valid_i4q!(wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    return select_bioedge_valid_i4q!(execution_style(wfs.estimator.state.valid_i4q), wfs, pupil, src)
end

function select_bioedge_valid_i4q!(::ScalarCPUStyle, wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    n_pixels = max(1, round(Int, wfs.acquisition.state.nominal_detector_resolution / (2 * wfs.acquisition.binning)))
    if size(wfs.estimator.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.valid_i4q = similar(wfs.estimator.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.estimator.state.valid_i4q, false)
    end
    if size(wfs.estimator.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.valid_signal = similar(wfs.estimator.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.flux_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.flux_i4q = similar(
            wfs.estimator.state.flux_i4q, n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.signal_2d = similar(wfs.estimator.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    elseif size(wfs.estimator.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(wfs.estimator.state.valid_i4q, true)
        update_bioedge_valid_signal!(wfs)
        update_bioedge_valid_signal_indices!(wfs)
        resize_bioedge_slope_buffers!(wfs)
        return wfs
    end

    bioedge_intensity_core!(wfs.front_end.propagation.intensity, wfs, pupil,
        src, bioedge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    frame = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)

    center = require_bioedge_frame_geometry(wfs, frame)
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
    cutoff = wfs.estimator.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_pixels + i, center - n_pixels + j]
        q2 = frame[center - n_pixels + i, center + j]
        q3 = frame[center + i, center + j]
        q4 = frame[center + i, center - n_pixels + j]
        wfs.estimator.state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
    end
    update_bioedge_valid_signal!(wfs)
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
    return wfs
end

function select_bioedge_valid_i4q!(::AcceleratorStyle, wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    n_pixels = max(1, round(Int, wfs.acquisition.state.nominal_detector_resolution / (2 * wfs.acquisition.binning)))
    if size(wfs.estimator.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.valid_i4q = similar(wfs.estimator.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.estimator.state.valid_i4q, false)
    end
    if size(wfs.estimator.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.valid_signal = similar(wfs.estimator.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.flux_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.flux_i4q = similar(
            wfs.estimator.state.flux_i4q, n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.signal_2d = similar(wfs.estimator.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    elseif size(wfs.estimator.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(wfs.estimator.state.valid_i4q, true)
        update_bioedge_valid_signal!(wfs)
        update_bioedge_valid_signal_indices!(wfs)
        resize_bioedge_slope_buffers!(wfs)
        return wfs
    end

    bioedge_intensity_core!(wfs.front_end.propagation.intensity, wfs, pupil,
        src, bioedge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    frame = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)

    center = require_bioedge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view wfs.estimator.state.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. wfs.estimator.state.valid_i4q = i4q >= cutoff
    update_bioedge_valid_signal!(wfs)
    update_bioedge_valid_signal_indices!(wfs)
    resize_bioedge_slope_buffers!(wfs)
    return wfs
end

function ensure_bioedge_calibration!(wfs::BioEdgeWFS, pupil::PupilFunction, src::AbstractSource)
    λ = calibration_wavelength(src, eltype(wfs.estimator.state.slopes))
    sig = pupil_aperture_calibration_signature(pupil,
        calibration_signature(src))
    if calibration_matches(wfs.estimator.state.calibrated,
        wfs.estimator.state.calibration_wavelength, λ,
        wfs.estimator.state.calibration_signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        select_bioedge_valid_i4q!(wfs, pupil, src)
        bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
        frame = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        bioedge_signal!(wfs, pupil, frame, src)
        store_reference_signal!(wfs.estimator.state.reference_signal_2d,
            wfs.estimator.state.signal_2d, wfs.estimator.state.slopes)
    finally
        restore_opd!(pupil, opd_saved)
    end
    wfs.estimator.state.calibrated = true
    wfs.estimator.state.calibration_wavelength = λ
    wfs.estimator.state.calibration_signature = sig
    wfs.estimator.state.calibration_revision += UInt(1)
    return wfs
end

@inline function ensure_bioedge_calibration!(wfs::BioEdgeWFS,
    pupil::PupilFunction, src::AbstractSource, ::AbstractDetector)
    return ensure_bioedge_calibration!(wfs, pupil, src)
end

function ensure_bioedge_calibration!(wfs::BioEdgeWFS, pupil::PupilFunction,
    src::AbstractSource, det::Detector)
    T = eltype(wfs.estimator.state.slopes)
    λ = calibration_wavelength(src, T)
    sig = detector_calibration_signature(det,
        pupil_aperture_calibration_signature(pupil,
            calibration_signature(src)))
    if calibration_matches(wfs.estimator.state.calibrated,
        wfs.estimator.state.calibration_wavelength, λ,
        wfs.estimator.state.calibration_signature, sig)
        return wfs
    end

    require_whole_capture_idle(det)
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        if !iszero(wfs.estimator.params.light_ratio)
            select_bioedge_valid_i4q!(wfs, pupil, src, det)
        end
        bioedge_intensity!(wfs.front_end.propagation.intensity, wfs, pupil, src)
        sampled = sample_bioedge_intensity!(wfs, pupil, wfs.front_end.propagation.intensity)
        frame = detector_calibration_frame!(det, sampled, src)
        resize_bioedge_signal_buffers!(wfs, size(frame, 1), det)
        if iszero(wfs.estimator.params.light_ratio)
            fill!(wfs.estimator.state.valid_i4q, true)
            update_bioedge_valid_signal!(wfs)
            update_bioedge_valid_signal_indices!(wfs)
            resize_bioedge_slope_buffers!(wfs)
        end
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        normalization_scale = wfs_detector_incidence_scale(det, src,
            eltype(frame))
        bioedge_signal!(wfs, pupil, frame, src, normalization_scale)
        store_reference_signal!(wfs.estimator.state.reference_signal_2d,
            wfs.estimator.state.signal_2d, wfs.estimator.state.slopes)
    finally
        restore_opd!(pupil, opd_saved)
    end
    wfs.estimator.state.calibrated = true
    wfs.estimator.state.calibration_wavelength = λ
    wfs.estimator.state.calibration_signature = sig
    wfs.estimator.state.calibration_revision += UInt(1)
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNone, intensity::AbstractMatrix{T}, wfs::BioEdgeWFS,
    ::PupilFunction, src::LGSSource) where {T<:AbstractFloat}
    wfs.front_end.propagation.elongation_kernel = apply_elongation!(
        intensity,
        lgs_elongation_factor(src),
        wfs.front_end.propagation.scratch,
        wfs.front_end.propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, intensity::AbstractMatrix{T}, wfs::BioEdgeWFS,
    pupil::PupilFunction, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, pupil, src)
    apply_lgs_convolution!(
        intensity,
        wfs.front_end.propagation.lgs_kernel_fft,
        wfs.front_end.propagation.fft_buffer,
        wfs.front_end.propagation.fft_plan,
        wfs.front_end.propagation.pupil_field,
        wfs.front_end.propagation.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::BioEdgeWFS, pupil::PupilFunction, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    pad = size(wfs.front_end.propagation.fft_buffer, 1)
    padding = wfs.front_end.propagation.effective_resolution / _pupil_resolution(pupil)
    pixel_scale = lgs_pixel_scale(_pupil_diameter_m(pupil), padding,
        wavelength(src))
    tag = lgs_kernel_signature(
        pupil,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        eltype(wfs.front_end.propagation.intensity);
        model=:subaperture_average,
    )
    if size(wfs.front_end.propagation.lgs_kernel_fft, 1) == pad && wfs.front_end.propagation.lgs_kernel_tag == tag
        return wfs
    end
    wfs.front_end.propagation.lgs_kernel_fft = lgs_average_kernel_fft(
        pupil,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        wfs.front_end.propagation.fft_buffer,
        wfs.front_end.propagation.fft_plan,
    )
    wfs.front_end.propagation.lgs_kernel_tag = tag
    return wfs
end

function set_optical_gain!(wfs::BioEdgeWFS, gain::Real)
    fill!(wfs.estimator.state.optical_gain, gain)
    return wfs
end

function set_optical_gain!(wfs::BioEdgeWFS, gain::AbstractVector)
    if length(gain) != length(wfs.estimator.state.optical_gain)
        throw(InvalidConfiguration("optical_gain length must match slope vector"))
    end
    copyto!(wfs.estimator.state.optical_gain, gain)
    return wfs
end
