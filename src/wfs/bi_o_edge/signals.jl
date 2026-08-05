function bi_o_edge_slopes!(wfs::BiOEdgeWFS, phase::AbstractMatrix, edge_mask::AbstractMatrix{Bool})
    slopes = bi_o_edge_estimator_products(wfs).slopes
    edge_geometric_slopes!(slopes, phase, wfs.estimator.state.valid_mask,
        edge_mask)
    return apply_bi_o_edge_optical_gain!(wfs)
end

@inline function apply_bi_o_edge_optical_gain!(wfs::BiOEdgeWFS)
    slopes = bi_o_edge_estimator_products(wfs).slopes
    @. slopes *= wfs.estimator.state.optical_gain
    return slopes
end

function bi_o_edge_slopes_intensity!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    intensity::AbstractMatrix{F}) where {F<:Real}
    return bi_o_edge_signal!(wfs, pupil, intensity)
end

function bi_o_edge_signal!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F}) where {F<:Real}
    return bi_o_edge_signal!(wfs, pupil, frame, nothing)
end

function bi_o_edge_signal!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(bi_o_edge_estimator_products(wfs).slopes)
    return bi_o_edge_signal!(wfs, pupil, frame, src, one(S))
end

function bi_o_edge_signal!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::Real) where {F<:Real}
    S = eltype(bi_o_edge_estimator_products(wfs).slopes)
    return bi_o_edge_signal!(execution_style(frame), wfs, pupil, frame, src,
        S(normalization_scale))
end

function bi_o_edge_signal!(::ScalarCPUStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(bi_o_edge_estimator_products(wfs).slopes)
    return bi_o_edge_signal!(ScalarCPUStyle(), wfs, pupil, frame, src, one(S))
end

function bi_o_edge_signal!(::ScalarCPUStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    return bi_o_edge_signal!(ScalarCPUStyle(), wfs, frame, src,
        normalization_scale)
end

function bi_o_edge_signal!(::ScalarCPUStyle, wfs::BiOEdgeWFS,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    products = bi_o_edge_estimator_products(wfs)
    n_pixels = size(workspace.signal_2d, 2)
    center = require_bi_o_edge_frame_geometry(wfs, frame)
    count = div(length(products.slopes), 2)
    norma = zero(S)
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = S(frame[center - n_pixels + i, center - n_pixels + j])
        q2 = S(frame[center - n_pixels + i, center + j])
        q3 = S(frame[center + i, center + j])
        q4 = S(frame[center + i, center - n_pixels + j])
        if state.valid_i4q[i, j]
            norma += q1 + q2 + q3 + q4
        end
    end
    norma = bi_o_edge_normalization(wfs.estimator.params.normalization, wfs, src,
        count, norma, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(workspace.signal_2d, zero(S))
        fill!(products.slopes, zero(S))
        return products.slopes
    end
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        q1 = S(frame[center - n_pixels + i, center - n_pixels + j])
        q2 = S(frame[center - n_pixels + i, center + j])
        q3 = S(frame[center + i, center + j])
        q4 = S(frame[center + i, center - n_pixels + j])
        sx = (q1 - q2 + q4 - q3) / norma
        sy = (q1 - q4 + q2 - q3) / norma
        workspace.signal_2d[i, j] = sx - state.reference_signal_2d[i, j]
        workspace.signal_2d[i + n_pixels, j] =
            sy - state.reference_signal_2d[i + n_pixels, j]
        if state.valid_i4q[i, j]
            products.slopes[idx] = workspace.signal_2d[i, j]
            products.slopes[idx + count] =
                workspace.signal_2d[i + n_pixels, j]
            idx += 1
        end
    end
    return apply_bi_o_edge_optical_gain!(wfs)
end

function bi_o_edge_signal!(style::AcceleratorStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource}) where {F<:Real}
    S = eltype(bi_o_edge_estimator_products(wfs).slopes)
    return bi_o_edge_signal!(style, wfs, pupil, frame, src, one(S))
end

function bi_o_edge_signal!(style::AcceleratorStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{F},
    src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    return bi_o_edge_signal!(style, wfs, frame, src, normalization_scale)
end

function bi_o_edge_signal!(style::AcceleratorStyle, wfs::BiOEdgeWFS,
    frame::AbstractMatrix{F}, src::Union{Nothing,AbstractSource},
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    products = bi_o_edge_estimator_products(wfs)
    count = workspace.valid_signal_count
    n_pixels = size(workspace.signal_2d, 2)
    center = require_bi_o_edge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    sx = @view workspace.signal_2d[1:n_pixels, :]
    sy = @view workspace.signal_2d[n_pixels+1:2*n_pixels, :]
    refx = @view state.reference_signal_2d[1:n_pixels, :]
    refy = @view state.reference_signal_2d[n_pixels+1:2*n_pixels, :]
    i4q = @view workspace.flux_i4q[1:n_pixels, 1:n_pixels]
    @. i4q = S(q1) + S(q2) + S(q3) + S(q4)
    summed_i4q = bi_o_edge_valid_flux_sum!(style, wfs, i4q)
    norma = bi_o_edge_normalization(wfs.estimator.params.normalization, wfs, src,
        count, summed_i4q, normalization_scale)
    if !usable_wfs_normalization(norma)
        fill!(workspace.signal_2d, zero(S))
        fill!(products.slopes, zero(S))
        return products.slopes
    end
    @. sx = (S(q1) - S(q2) + S(q4) - S(q3)) / norma - refx
    @. sy = (S(q1) - S(q4) + S(q2) - S(q3)) / norma - refy
    launch_kernel!(style, gather_bi_o_edge_slopes_kernel!, products.slopes,
        workspace.signal_2d, workspace.valid_signal_indices, count, n_pixels;
        ndrange=count)
    return apply_bi_o_edge_optical_gain!(wfs)
end

function bi_o_edge_normalization(normalization::MeanValidFluxNormalization,
    wfs::BiOEdgeWFS, pupil::PupilFunction, src::Union{Nothing,AbstractSource},
    count::Int, summed_i4q)
    return bi_o_edge_normalization(normalization, wfs, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function bi_o_edge_normalization(normalization::MeanValidFluxNormalization,
    wfs::BiOEdgeWFS, ::PupilFunction, src::Union{Nothing,AbstractSource},
    count::Int, summed_i4q, normalization_scale::Real)
    return bi_o_edge_normalization(normalization, wfs, src, count,
        summed_i4q, normalization_scale)
end

function bi_o_edge_normalization(::MeanValidFluxNormalization, ::BiOEdgeWFS,
    ::Union{Nothing,AbstractSource}, count::Int, summed_i4q,
    ::Real)
    T = typeof(summed_i4q)
    return count == 0 ? one(T) : summed_i4q / count
end

function bi_o_edge_normalization(normalization::IncidenceFluxNormalization,
    wfs::BiOEdgeWFS, pupil::PupilFunction, src::AbstractSource, count::Int,
    summed_i4q)
    return bi_o_edge_normalization(normalization, wfs, src, count,
        summed_i4q, one(typeof(summed_i4q)))
end

function bi_o_edge_normalization(normalization::IncidenceFluxNormalization,
    wfs::BiOEdgeWFS, ::PupilFunction, src::AbstractSource, count::Int,
    summed_i4q, normalization_scale::Real)
    return bi_o_edge_normalization(normalization, wfs, src, count, summed_i4q,
        normalization_scale)
end

function bi_o_edge_normalization(::IncidenceFluxNormalization,
    wfs::BiOEdgeWFS, src::AbstractSource, ::Int,
    summed_i4q, normalization_scale::Real)
    T = typeof(summed_i4q)
    sub_area = (wfs.estimator.params.pupil_diameter_m /
        wfs.estimator.params.pupil_samples)^2
    return wfs_incident_photon_irradiance(src, T) * T(sub_area) *
        T(normalization_scale)
end

function bi_o_edge_normalization(::IncidenceFluxNormalization, ::BiOEdgeWFS, ::PupilFunction,
    ::Nothing, ::Int, summed_i4q)
    return one(typeof(summed_i4q))
end


function bi_o_edge_normalization(::IncidenceFluxNormalization, ::BiOEdgeWFS,
    ::PupilFunction, ::Nothing, ::Int, summed_i4q, ::Real)
    return one(typeof(summed_i4q))
end

function bi_o_edge_normalization(::IncidenceFluxNormalization, ::BiOEdgeWFS,
    ::Nothing, ::Int, summed_i4q, ::Real)
    return one(typeof(summed_i4q))
end

function update_bi_o_edge_valid_signal!(wfs::BiOEdgeWFS)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    n_pixels = size(state.valid_i4q, 1)
    fill!(workspace.valid_signal, false)
    @views begin
        workspace.valid_signal[1:n_pixels, :] .= state.valid_i4q
        workspace.valid_signal[n_pixels+1:end, :] .= state.valid_i4q
    end
    return wfs
end

function update_bi_o_edge_valid_signal_indices!(wfs::BiOEdgeWFS)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    valid_host = workspace.valid_i4q_host
    if size(valid_host) != size(state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(state.valid_i4q)...)
        workspace.valid_i4q_host = valid_host
    end
    copyto!(valid_host, state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(workspace.valid_signal_indices) < n_valid
        workspace.valid_signal_indices = similar(
            workspace.valid_signal_indices, n_valid)
    end
    if length(workspace.valid_signal_indices_host) < n_valid
        workspace.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = workspace.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(workspace.valid_signal_indices, 1, host_indices, 1, n_valid)
    workspace.valid_signal_count = n_valid
    return n_valid
end

function resize_bi_o_edge_slope_buffers!(wfs::BiOEdgeWFS)
    workspace = bi_o_edge_estimator_workspace(wfs)
    products = bi_o_edge_estimator_products(wfs)
    state = bi_o_edge_estimator_state(wfs)
    n_valid = workspace.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("bi_o_edge valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(products.slopes) != n_slopes
        products.slopes = similar(products.slopes, n_slopes)
    end
    if length(state.optical_gain) != n_slopes
        state.optical_gain = similar(state.optical_gain, n_slopes)
        fill!(state.optical_gain, one(eltype(state.optical_gain)))
    end
    return wfs
end

function bi_o_edge_valid_flux_sum!(::ScalarCPUStyle, wfs::BiOEdgeWFS,
    i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    workspace = bi_o_edge_estimator_workspace(wfs)
    return masked_sum2d(ScalarCPUStyle(), i4q, workspace.valid_i4q_host)
end

function bi_o_edge_valid_flux_sum!(style::AcceleratorStyle,
    wfs::BiOEdgeWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    workspace = bi_o_edge_estimator_workspace(wfs)
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.estimator.state.valid_i4q,
        workspace.valid_i4q_host,
        workspace.valid_flux_sum_buffer,
        workspace.valid_flux_sum_host,
        workspace.valid_flux_i4q_host,
    )
    workspace.valid_flux_i4q_host = host_parent
    return summed
end

function select_bi_o_edge_valid_i4q_from_frame!(::ScalarCPUStyle,
    wfs::BiOEdgeWFS, frame::AbstractMatrix)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center = require_bi_o_edge_frame_geometry(wfs, frame)
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
    update_bi_o_edge_valid_signal!(wfs)
    update_bi_o_edge_valid_signal_indices!(wfs)
    resize_bi_o_edge_slope_buffers!(wfs)
    return wfs
end

function select_bi_o_edge_valid_i4q_from_frame!(::AcceleratorStyle,
    wfs::BiOEdgeWFS, frame::AbstractMatrix)
    workspace = bi_o_edge_estimator_workspace(wfs)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center = require_bi_o_edge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view workspace.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. wfs.estimator.state.valid_i4q = i4q >= cutoff
    update_bi_o_edge_valid_signal!(wfs)
    update_bi_o_edge_valid_signal_indices!(wfs)
    resize_bi_o_edge_slope_buffers!(wfs)
    return wfs
end

function select_bi_o_edge_valid_i4q!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::AbstractSource, det::Detector)
    propagation = bi_o_edge_propagation_workspace(wfs)
    bi_o_edge_intensity_core!(propagation.intensity, wfs, pupil, src,
        bi_o_edge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    sampled = sample_bi_o_edge_intensity!(wfs, pupil, propagation.intensity)
    frame = detector_calibration_frame!(det, sampled, src)
    resize_bi_o_edge_signal_buffers!(wfs, size(frame, 1), det)
    return select_bi_o_edge_valid_i4q_from_frame!(execution_style(frame), wfs,
        frame)
end

function select_bi_o_edge_valid_i4q!(wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource)
    style = execution_style(wfs.estimator.state.valid_i4q)
    return select_bi_o_edge_valid_i4q!(style, wfs, pupil, src)
end

function select_bi_o_edge_valid_i4q!(::ScalarCPUStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    acquisition = bi_o_edge_acquisition_workspace(wfs)
    binning = bi_o_edge_acquisition_plan(wfs).binning
    propagation = bi_o_edge_propagation_workspace(wfs)
    n_pixels = max(1, round(Int,
        acquisition.nominal_detector_resolution / (2 * binning)))
    if size(state.valid_i4q) != (n_pixels, n_pixels)
        state.valid_i4q = similar(state.valid_i4q, n_pixels, n_pixels)
        fill!(state.valid_i4q, false)
    end
    if size(workspace.valid_signal) != (2 * n_pixels, n_pixels)
        workspace.valid_signal = similar(workspace.valid_signal,
            2 * n_pixels, n_pixels)
    end
    if size(workspace.flux_i4q) != (n_pixels, n_pixels)
        workspace.flux_i4q = similar(workspace.flux_i4q,
            n_pixels, n_pixels)
    end
    if size(workspace.signal_2d) != (2 * n_pixels, n_pixels)
        workspace.signal_2d = similar(workspace.signal_2d,
            2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
    elseif size(state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(state.valid_i4q, true)
        update_bi_o_edge_valid_signal!(wfs)
        update_bi_o_edge_valid_signal_indices!(wfs)
        resize_bi_o_edge_slope_buffers!(wfs)
        return wfs
    end

    bi_o_edge_intensity_core!(propagation.intensity, wfs, pupil, src,
        bi_o_edge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    frame = sample_bi_o_edge_intensity!(wfs, pupil, propagation.intensity)

    center = require_bi_o_edge_frame_geometry(wfs, frame)
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
        state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
    end
    update_bi_o_edge_valid_signal!(wfs)
    update_bi_o_edge_valid_signal_indices!(wfs)
    resize_bi_o_edge_slope_buffers!(wfs)
    return wfs
end

function select_bi_o_edge_valid_i4q!(::AcceleratorStyle, wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    acquisition = bi_o_edge_acquisition_workspace(wfs)
    binning = bi_o_edge_acquisition_plan(wfs).binning
    propagation = bi_o_edge_propagation_workspace(wfs)
    n_pixels = max(1, round(Int,
        acquisition.nominal_detector_resolution / (2 * binning)))
    if size(state.valid_i4q) != (n_pixels, n_pixels)
        state.valid_i4q = similar(state.valid_i4q, n_pixels, n_pixels)
        fill!(state.valid_i4q, false)
    end
    if size(workspace.valid_signal) != (2 * n_pixels, n_pixels)
        workspace.valid_signal = similar(workspace.valid_signal,
            2 * n_pixels, n_pixels)
    end
    if size(workspace.flux_i4q) != (n_pixels, n_pixels)
        workspace.flux_i4q = similar(workspace.flux_i4q,
            n_pixels, n_pixels)
    end
    if size(workspace.signal_2d) != (2 * n_pixels, n_pixels)
        workspace.signal_2d = similar(workspace.signal_2d,
            2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
    elseif size(state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(state.valid_i4q, true)
        update_bi_o_edge_valid_signal!(wfs)
        update_bi_o_edge_valid_signal_indices!(wfs)
        resize_bi_o_edge_slope_buffers!(wfs)
        return wfs
    end

    bi_o_edge_intensity_core!(propagation.intensity, wfs, pupil, src,
        bi_o_edge_calibration_modulation(wfs);
        apply_lgs=src isa LGSSource)
    frame = sample_bi_o_edge_intensity!(wfs, pupil, propagation.intensity)

    center = require_bi_o_edge_frame_geometry(wfs, frame)
    rows_lo = center - n_pixels + 1:center
    rows_hi = center + 1:center + n_pixels
    cols_lo = center - n_pixels + 1:center
    cols_hi = center + 1:center + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view workspace.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. state.valid_i4q = i4q >= cutoff
    update_bi_o_edge_valid_signal!(wfs)
    update_bi_o_edge_valid_signal_indices!(wfs)
    resize_bi_o_edge_slope_buffers!(wfs)
    return wfs
end

function ensure_bi_o_edge_calibration!(wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    products = bi_o_edge_estimator_products(wfs)
    propagation = bi_o_edge_propagation_workspace(wfs)
    λ = calibration_wavelength(src, eltype(products.slopes))
    sig = pupil_aperture_calibration_signature(pupil,
        calibration_signature(src))
    if calibration_matches(state.calibrated,
        state.calibration_wavelength, λ, state.calibration_signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        select_bi_o_edge_valid_i4q!(wfs, pupil, src)
        bi_o_edge_intensity!(propagation.intensity, wfs, pupil, src)
        frame = sample_bi_o_edge_intensity!(wfs, pupil,
            propagation.intensity)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        bi_o_edge_signal!(wfs, pupil, frame, src)
        store_reference_signal!(state.reference_signal_2d,
            workspace.signal_2d, products.slopes)
    finally
        restore_opd!(pupil, opd_saved)
    end
    state.calibrated = true
    state.calibration_wavelength = λ
    state.calibration_signature = sig
    state.calibration_revision += UInt(1)
    return wfs
end

@inline function ensure_bi_o_edge_calibration!(wfs::BiOEdgeWFS,
    pupil::PupilFunction, src::AbstractSource, ::AbstractDetector)
    return ensure_bi_o_edge_calibration!(wfs, pupil, src)
end

function ensure_bi_o_edge_calibration!(wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::AbstractSource, det::Detector)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    products = bi_o_edge_estimator_products(wfs)
    propagation = bi_o_edge_propagation_workspace(wfs)
    T = eltype(products.slopes)
    λ = calibration_wavelength(src, T)
    sig = detector_calibration_signature(det,
        pupil_aperture_calibration_signature(pupil,
            calibration_signature(src)))
    if calibration_matches(state.calibrated,
        state.calibration_wavelength, λ, state.calibration_signature, sig)
        return wfs
    end

    require_whole_capture_idle(det)
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        if !iszero(wfs.estimator.params.light_ratio)
            select_bi_o_edge_valid_i4q!(wfs, pupil, src, det)
        end
        bi_o_edge_intensity!(propagation.intensity, wfs, pupil, src)
        sampled = sample_bi_o_edge_intensity!(wfs, pupil,
            propagation.intensity)
        frame = detector_calibration_frame!(det, sampled, src)
        resize_bi_o_edge_signal_buffers!(wfs, size(frame, 1), det)
        if iszero(wfs.estimator.params.light_ratio)
            fill!(state.valid_i4q, true)
            update_bi_o_edge_valid_signal!(wfs)
            update_bi_o_edge_valid_signal_indices!(wfs)
            resize_bi_o_edge_slope_buffers!(wfs)
        end
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        normalization_scale = wfs_detector_incidence_scale(det, src,
            eltype(frame))
        bi_o_edge_signal!(wfs, pupil, frame, src, normalization_scale)
        store_reference_signal!(state.reference_signal_2d,
            workspace.signal_2d, products.slopes)
    finally
        restore_opd!(pupil, opd_saved)
    end
    state.calibrated = true
    state.calibration_wavelength = λ
    state.calibration_signature = sig
    state.calibration_revision += UInt(1)
    return wfs
end

function apply_lgs_elongation!(::NoSodiumLayerProfileStyle,
    intensity::AbstractMatrix{T}, wfs::BiOEdgeWFS, ::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    propagation = bi_o_edge_propagation_workspace(wfs)
    propagation.elongation_kernel = apply_elongation!(
        intensity,
        lgs_elongation_factor(src),
        propagation.scratch,
        propagation.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::SampledSodiumLayerProfileStyle,
    intensity::AbstractMatrix{T}, wfs::BiOEdgeWFS, pupil::PupilFunction,
    src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, pupil, src)
    propagation = bi_o_edge_propagation_workspace(wfs)
    apply_lgs_convolution!(
        intensity,
        propagation.lgs_kernel_fft,
        propagation.fft_buffer,
        propagation.fft_plan,
        propagation.pupil_field,
        propagation.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::BiOEdgeWFS, pupil::PupilFunction, src::LGSSource)
    profile = src.params.sodium_layer_profile
    if profile === nothing
        return wfs
    end
    propagation = bi_o_edge_propagation_workspace(wfs)
    pad = size(propagation.fft_buffer, 1)
    padding = propagation.effective_resolution / _pupil_resolution(pupil)
    pixel_scale = lgs_pixel_scale(_pupil_diameter_m(pupil), padding,
        wavelength(src))
    tag = lgs_kernel_signature(
        pupil,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        eltype(propagation.intensity);
        model=:subaperture_average,
    )
    if size(propagation.lgs_kernel_fft, 1) == pad &&
        propagation.lgs_kernel_tag == tag
        return wfs
    end
    propagation.lgs_kernel_fft = lgs_average_kernel_fft(
        pupil,
        src,
        pad,
        wfs.estimator.params.pupil_samples,
        pixel_scale,
        propagation.fft_buffer,
        propagation.fft_plan,
    )
    propagation.lgs_kernel_tag = tag
    return wfs
end

function set_optical_gain!(wfs::BiOEdgeWFS, gain::Real)
    fill!(wfs.estimator.state.optical_gain, gain)
    return wfs
end

function set_optical_gain!(wfs::BiOEdgeWFS, gain::AbstractVector)
    if length(gain) != length(wfs.estimator.state.optical_gain)
        throw(InvalidConfiguration("optical_gain length must match slope vector"))
    end
    copyto!(wfs.estimator.state.optical_gain, gain)
    return wfs
end
