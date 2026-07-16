
function pyramid_slopes!(wfs::PyramidWFS, tel::Telescope)
    return pyramid_signal!(wfs, tel, wfs.state.camera_frame)
end

function pyramid_slopes!(wfs::PyramidWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    return pyramid_signal!(wfs, tel, intensity)
end

function pyramid_signal!(wfs::PyramidWFS, tel::Telescope, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    return pyramid_signal!(wfs, tel, frame, nothing)
end

function pyramid_signal!(wfs::PyramidWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    return pyramid_signal!(wfs, tel, frame, src, one(T))
end

function pyramid_signal!(wfs::PyramidWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}, normalization_scale::Real) where {T<:AbstractFloat}
    return pyramid_signal!(execution_style(frame), wfs, tel, frame, src,
        T(normalization_scale))
end

function pyramid_signal!(::ScalarCPUStyle, wfs::PyramidWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    return pyramid_signal!(ScalarCPUStyle(), wfs, tel, frame, src, one(T))
end

function pyramid_signal!(::ScalarCPUStyle, wfs::PyramidWFS, tel::Telescope,
    frame::AbstractMatrix{T}, src::Union{Nothing,AbstractSource},
    normalization_scale::T) where {T<:AbstractFloat}
    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    n_pixels = size(wfs.state.signal_2d, 2)
    count = div(length(wfs.state.slopes), 2)
    norma = zero(T)
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i, center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i, center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i, center - n_extra - n_pixels + j]
        if wfs.state.valid_i4q[i, j]
            norma += q1 + q2 + q3 + q4
        end
    end
    norma = pyramid_normalization(wfs.params.normalization, wfs, tel, src,
        count, norma, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(wfs.state.signal_2d, zero(T))
        fill!(wfs.state.slopes, zero(T))
        return wfs.state.slopes
    end
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i, center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i, center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i, center - n_extra - n_pixels + j]
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
    return wfs.state.slopes
end

function pyramid_signal!(style::AcceleratorStyle, wfs::PyramidWFS, tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}) where {T<:AbstractFloat}
    return pyramid_signal!(style, wfs, tel, frame, src, one(T))
end

function pyramid_signal!(style::AcceleratorStyle, wfs::PyramidWFS,
    tel::Telescope, frame::AbstractMatrix{T},
    src::Union{Nothing,AbstractSource}, normalization_scale::T) where {T<:AbstractFloat}
    count = wfs.state.valid_signal_count
    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    n_pixels = size(wfs.state.signal_2d, 2)
    rows_lo = center - n_extra - n_pixels + 1:center - n_extra
    rows_hi = center + n_extra + 1:center + n_extra + n_pixels
    cols_lo = center - n_extra - n_pixels + 1:center - n_extra
    cols_hi = center + n_extra + 1:center + n_extra + n_pixels
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
    summed_i4q = pyramid_valid_flux_sum!(style, wfs, i4q)
    norma = pyramid_normalization(wfs.params.normalization, wfs, tel, src,
        count, summed_i4q, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(wfs.state.signal_2d, zero(T))
        fill!(wfs.state.slopes, zero(T))
        return wfs.state.slopes
    end
    @. sx = (q1 - q2 + q4 - q3) / norma - refx
    @. sy = (q1 - q4 + q2 - q3) / norma - refy
    launch_kernel!(style, gather_pyramid_slopes_kernel!, wfs.state.slopes,
        wfs.state.signal_2d, wfs.state.valid_signal_indices, count, n_pixels; ndrange=count)
    return wfs.state.slopes
end

function pyramid_normalization(normalization::MeanValidFluxNormalization,
    wfs::PyramidWFS, tel::Telescope, src::Union{Nothing,AbstractSource},
    count::Int, summed_i4q)
    return pyramid_normalization(normalization, wfs, tel, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function pyramid_normalization(::MeanValidFluxNormalization, ::PyramidWFS,
    ::Telescope, ::Union{Nothing,AbstractSource}, count::Int, summed_i4q,
    ::Real)
    T = typeof(summed_i4q)
    return count == 0 ? one(T) : summed_i4q / count
end

function pyramid_normalization(normalization::IncidenceFluxNormalization,
    wfs::PyramidWFS, tel::Telescope, src::AbstractSource, count::Int,
    summed_i4q)
    return pyramid_normalization(normalization, wfs, tel, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function pyramid_normalization(::IncidenceFluxNormalization, wfs::PyramidWFS,
    tel::Telescope, src::AbstractSource, ::Int, summed_i4q,
    normalization_scale::Real)
    T = typeof(summed_i4q)
    sub_area = (tel.params.diameter / wfs.params.pupil_samples)^2
    return wfs_incident_photon_irradiance(src, T) * T(sub_area) *
        T(normalization_scale)
end

function pyramid_normalization(::IncidenceFluxNormalization, ::PyramidWFS, ::Telescope,
    ::Nothing, ::Int, summed_i4q)
    return one(typeof(summed_i4q))
end

function pyramid_normalization(::IncidenceFluxNormalization, ::PyramidWFS,
    ::Telescope, ::Nothing, ::Int, summed_i4q, ::Real)
    return one(typeof(summed_i4q))
end

@inline pyramid_calibration_qe_scale(::AbstractSource, ::Nothing,
    ::Type{T}) where {T<:AbstractFloat} = one(T)

@inline function pyramid_calibration_qe_scale(src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel,
    ::Type{T}) where {T<:AbstractFloat}
    return effective_qe(qe_model, src, T)
end

@inline function pyramid_calibration_signature(src::AbstractSource, ::Nothing)
    return calibration_signature(src)
end

function pyramid_calibration_signature(src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel)
    sig = calibration_signature(src)
    @inbounds for sample in spectral_bundle(src)
        sig = hash(qe_at(qe_model, sample.wavelength), sig)
    end
    return sig
end

@inline pyramid_detector_qe_model(::AbstractSource, ::Detector) = nothing
@inline pyramid_detector_qe_model(::SpectralSource, det::Detector) =
    quantum_efficiency_model(det)

@inline function pyramid_detector_calibration_qe(src::AbstractSource,
    det::Detector, ::Type{T}) where {T<:AbstractFloat}
    return effective_qe(det, src, T)
end

@inline pyramid_detector_calibration_qe(::SpectralSource, ::Detector,
    ::Type{T}) where {T<:AbstractFloat} = one(T)

function ensure_pyramid_calibration!(wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource)
    return ensure_pyramid_calibration!(wfs, tel, src, nothing)
end

@inline function ensure_pyramid_calibration!(wfs::PyramidWFS,
    tel::Telescope, src::AbstractSource, ::AbstractDetector)
    return ensure_pyramid_calibration!(wfs, tel, src)
end

function ensure_pyramid_calibration!(wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel})
    λ = calibration_wavelength(src, eltype(wfs.state.slopes))
    sig = telescope_aperture_calibration_signature(tel,
        pyramid_calibration_signature(src, qe_model))
    if calibration_matches(wfs.state.calibrated, wfs.state.calibration_wavelength, λ,
        wfs.state.calibration_signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, tel)
    opd_saved = save_zero_opd!(tel)
    try
        select_pyramid_valid_i4q!(wfs, tel, src, qe_model)
        pyramid_calibration_intensity!(wfs.state.intensity, wfs, tel, src,
            qe_model)
        frame = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
        fill!(wfs.state.reference_signal_2d,
            zero(eltype(wfs.state.reference_signal_2d)))
        normalization_scale = pyramid_calibration_qe_scale(src, qe_model,
            eltype(frame))
        pyramid_signal!(wfs, tel, frame, src, normalization_scale)
        store_reference_signal!(wfs.state.reference_signal_2d,
            wfs.state.signal_2d, wfs.state.slopes)
    finally
        restore_opd!(tel, opd_saved)
    end
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    wfs.state.calibration_signature = sig
    return wfs
end

function ensure_pyramid_calibration!(wfs::PyramidWFS, tel::Telescope,
    src::AbstractSource, det::Detector)
    T = eltype(wfs.state.slopes)
    λ = calibration_wavelength(src, T)
    qe_model = pyramid_detector_qe_model(src, det)
    sig = detector_calibration_signature(det,
        telescope_aperture_calibration_signature(tel,
            pyramid_calibration_signature(src, qe_model)))
    if calibration_matches(wfs.state.calibrated,
        wfs.state.calibration_wavelength, λ,
        wfs.state.calibration_signature, sig)
        return wfs
    end

    require_whole_capture_idle(det)
    update_valid_mask!(wfs, tel)
    opd_saved = save_zero_opd!(tel)
    try
        if !iszero(wfs.params.light_ratio)
            select_pyramid_valid_i4q!(wfs, tel, src, qe_model, det)
        end
        pyramid_calibration_intensity!(wfs.state.intensity, wfs, tel, src,
            qe_model)
        sampled = sample_pyramid_intensity!(wfs, tel, wfs.state.intensity)
        frame = detector_calibration_frame!(det, sampled,
            pyramid_detector_calibration_qe(src, det,
                eltype(det.state.frame)))
        resize_pyramid_signal_buffers!(wfs, size(frame, 1))
        if iszero(wfs.params.light_ratio)
            fill!(wfs.state.valid_i4q, true)
            update_pyramid_valid_signal!(wfs)
            update_pyramid_valid_signal_indices!(wfs)
            resize_pyramid_slope_buffers!(wfs)
        end
        fill!(wfs.state.reference_signal_2d,
            zero(eltype(wfs.state.reference_signal_2d)))
        normalization_scale = wfs_detector_incidence_scale(det, src,
            eltype(frame))
        pyramid_signal!(wfs, tel, frame, src, normalization_scale)
        store_reference_signal!(wfs.state.reference_signal_2d,
            wfs.state.signal_2d, wfs.state.slopes)
    finally
        restore_opd!(tel, opd_saved)
    end
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    wfs.state.calibration_signature = sig
    return wfs
end

function _pyramid_slopes!(::ScalarCPUStyle, slopes::AbstractVector, intensity::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, ::Int, offset::Int, ox1::Int, oy1::Int, ox2::Int, oy2::Int, ox3::Int, oy3::Int, ox4::Int, oy4::Int,
    shift_x::NTuple{4,Int}, shift_y::NTuple{4,Int}) where {T<:AbstractFloat}
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if valid_mask[i, j]
            q1 = quadrant_sum(intensity, xs, ys, sub, ox1, oy1, shift_x[1], shift_y[1])
            q2 = quadrant_sum(intensity, xs, ys, sub, ox2, oy2, shift_x[2], shift_y[2])
            q3 = quadrant_sum(intensity, xs, ys, sub, ox3, oy3, shift_x[3], shift_y[3])
            q4 = quadrant_sum(intensity, xs, ys, sub, ox4, oy4, shift_x[4], shift_y[4])
            left = q1 + q3
            right = q2 + q4
            bottom = q1 + q2
            top = q3 + q4
            total = left + right
            if total <= 0
                slopes[idx] = zero(eltype(slopes))
                slopes[idx + offset] = zero(eltype(slopes))
            else
                slopes[idx] = (right - left) / total
                slopes[idx + offset] = (top - bottom) / total
            end
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
        end
        idx += 1
    end
    return slopes
end

function _pyramid_slopes!(style::AcceleratorStyle, slopes::AbstractVector, intensity::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, pad::Int, offset::Int, ox1::Int, oy1::Int, ox2::Int, oy2::Int, ox3::Int, oy3::Int, ox4::Int, oy4::Int,
    shift_x::NTuple{4,Int}, shift_y::NTuple{4,Int})
    launch_kernel!(style, pyramid_slopes_kernel!, slopes, intensity, valid_mask, sub, n_sub, pad, offset,
        ox1, oy1, ox2, oy2, ox3, oy3, ox4, oy4,
        shift_x[1], shift_y[1], shift_x[2], shift_y[2], shift_x[3], shift_y[3], shift_x[4], shift_y[4];
        ndrange=(n_sub, n_sub))
    return slopes
end

function quadrant_sum(intensity::AbstractMatrix{T}, xs::Int, ys::Int, sub::Int,
    ox::Int, oy::Int, sx::Int, sy::Int) where {T<:AbstractFloat}
    acc = zero(T)
    n = size(intensity, 1)
    @inbounds for i in 0:(sub - 1), j in 0:(sub - 1)
        xi = ox + xs + i - 1 + sx
        yj = oy + ys + j - 1 + sy
        if 1 <= xi <= n && 1 <= yj <= n
            acc += intensity[xi, yj]
        end
    end
    return acc
end

function apply_shift_wfs!(wfs::PyramidWFS; sx, sy)
    shift_x, shift_y = pyramid_shift_components(sx, sy)
    wfs.state.shift_x = shift_x
    wfs.state.shift_y = shift_y
    return wfs
end

@inline function pyramid_shift_components(sx::Real, sy::Real)
    shift_x = (round(Int, sx), round(Int, sx), round(Int, sx), round(Int, sx))
    shift_y = (round(Int, sy), round(Int, sy), round(Int, sy), round(Int, sy))
    return shift_x, shift_y
end

@inline function pyramid_shift_components(sx, sy)
    if length(sx) != 4 || length(sy) != 4
        throw(InvalidConfiguration("pyramid shift must have 4 elements"))
    end
    shift_x = (round(Int, sx[1]), round(Int, sx[2]), round(Int, sx[3]), round(Int, sx[4]))
    shift_y = (round(Int, sy[1]), round(Int, sy[2]), round(Int, sy[3]), round(Int, sy[4]))
    return shift_x, shift_y
end

function set_optical_gain!(wfs::PyramidWFS, gain::Real)
    fill!(wfs.state.optical_gain, gain)
    return wfs
end

function set_optical_gain!(wfs::PyramidWFS, gain::AbstractVector)
    if length(gain) != length(wfs.state.optical_gain)
        throw(InvalidConfiguration("optical_gain length must match slope vector"))
    end
    copyto!(wfs.state.optical_gain, gain)
    return wfs
end
