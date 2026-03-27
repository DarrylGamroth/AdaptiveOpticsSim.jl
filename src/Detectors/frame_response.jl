apply_frame_response!(::NullFrameResponse, det::Detector) = det.state.frame

function apply_frame_response!(model::SeparableGaussianPixelResponse, det::Detector)
    return apply_frame_response!(execution_style(det.state.frame), model, det)
end

function apply_frame_response!(::ScalarCPUStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * frame[i, jj]
        end
        scratch[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * scratch[ii, j]
        end
        frame[i, j] = acc
    end
    return frame
end

function apply_frame_response!(style::AcceleratorStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    launch_kernel!(style, separable_response_rows_kernel!, scratch, frame, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    launch_kernel!(style, separable_response_cols_kernel!, frame, scratch, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    return frame
end

function ensure_buffers!(det::Detector, n_mid::Int, m_mid::Int, n_out::Int, m_out::Int)
    if size(det.state.frame) != (n_out, m_out)
        det.state.frame = similar(det.state.frame, n_out, m_out)
    end
    if size(det.state.response_buffer) != (n_out, m_out)
        det.state.response_buffer = similar(det.state.response_buffer, n_out, m_out)
    end
    if size(det.state.bin_buffer) != (n_mid, m_mid)
        det.state.bin_buffer = similar(det.state.bin_buffer, n_mid, m_mid)
    end
    if size(det.state.noise_buffer) != (n_out, m_out)
        det.state.noise_buffer = similar(det.state.noise_buffer, n_out, m_out)
    end
    if size(det.state.accum_buffer) != (n_out, m_out)
        det.state.accum_buffer = similar(det.state.accum_buffer, n_out, m_out)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    window = det.params.readout_window === nothing ? nothing : validate_readout_window(det.params.readout_window, n_out, m_out)
    out_rows = window === nothing ? n_out : length(window.rows)
    out_cols = window === nothing ? m_out : length(window.cols)
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != (out_rows, out_cols)
        det.state.output_buffer = similar(det.state.output_buffer, out_rows, out_cols)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}, exposure_time::Real) where {T}
    n_in, m_in = size(psf)
    sampling = det.params.psf_sampling
    binning = det.params.binning
    if sampling < 1 || binning < 1
        throw(InvalidConfiguration("psf_sampling and binning must be >= 1"))
    end
    if n_in % sampling != 0 || m_in % sampling != 0
        throw(DimensionMismatchError("psf_sampling must evenly divide input dimensions"))
    end
    n_mid = div(n_in, sampling)
    m_mid = div(m_in, sampling)
    if n_mid % binning != 0 || m_mid % binning != 0
        throw(DimensionMismatchError("binning must evenly divide sampled dimensions"))
    end
    n_out = div(n_mid, binning)
    m_out = div(m_mid, binning)
    ensure_buffers!(det, n_mid, m_mid, n_out, m_out)

    if sampling > 1
        if binning > 1
            bin2d!(det.state.frame, psf, sampling * binning)
        else
            bin2d!(det.state.bin_buffer, psf, sampling)
            det.state.frame .= det.state.bin_buffer
        end
    else
        if binning > 1
            bin2d!(det.state.frame, psf, binning)
        else
            det.state.frame .= psf
        end
    end
    @. det.state.frame *= det.params.qe * exposure_time
    apply_frame_response!(det.params.response_model, det)
    return det.state.frame
end

fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T} = fill_frame!(det, psf, det.params.integration_time)
