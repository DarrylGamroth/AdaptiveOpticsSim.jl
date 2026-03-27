apply_response!(::NullFrameResponse, det::Detector) = det.state.frame

function apply_response!(model::AbstractFrameResponse, det::Detector)
    return apply_response!(execution_style(det.state.frame), model, det.state.frame, det.state.response_buffer)
end

function _apply_separable_response!(::ScalarCPUStyle, frame::AbstractMatrix, scratch::AbstractMatrix,
    kernel_y::AbstractVector, kernel_x::AbstractVector)
    n, m = size(frame)
    radius_x = fld(length(kernel_x), 2)
    radius_y = fld(length(kernel_y), 2)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel_x)
            jj = clamp(j + kk - radius_x - 1, 1, m)
            acc += kernel_x[kk] * frame[i, jj]
        end
        scratch[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel_y)
            ii = clamp(i + kk - radius_y - 1, 1, n)
            acc += kernel_y[kk] * scratch[ii, j]
        end
        frame[i, j] = acc
    end
    return frame
end

function _apply_separable_response!(style::AcceleratorStyle, frame::AbstractMatrix, scratch::AbstractMatrix,
    kernel_y::AbstractVector, kernel_x::AbstractVector)
    n, m = size(frame)
    radius_x = fld(length(kernel_x), 2)
    radius_y = fld(length(kernel_y), 2)
    launch_kernel_async!(style, separable_response_rows_kernel!, scratch, frame, kernel_x, radius_x, n, m, length(kernel_x); ndrange=(n, m))
    launch_kernel_async!(style, separable_response_cols_kernel!, frame, scratch, kernel_y, radius_y, n, m, length(kernel_y); ndrange=(n, m))
    return frame
end

function apply_response!(style::ExecutionStyle, model::GaussianPixelResponse, frame::AbstractMatrix, scratch::AbstractMatrix)
    return _apply_separable_response!(style, frame, scratch, model.kernel, model.kernel)
end

function apply_response!(::ScalarCPUStyle, model::SampledFrameResponse, frame::AbstractMatrix, scratch::AbstractMatrix)
    n, m = size(frame)
    kn, km = size(model.kernel)
    radius_i = fld(kn, 2)
    radius_j = fld(km, 2)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for ki in 1:kn
            ii = clamp(i + ki - radius_i - 1, 1, n)
            for kj in 1:km
                jj = clamp(j + kj - radius_j - 1, 1, m)
                acc += model.kernel[ki, kj] * frame[ii, jj]
            end
        end
        scratch[i, j] = acc
    end
    copyto!(frame, scratch)
    return frame
end

function apply_response!(style::AcceleratorStyle, model::SampledFrameResponse, frame::AbstractMatrix, scratch::AbstractMatrix)
    n, m = size(frame)
    kn, km = size(model.kernel)
    radius_i = fld(kn, 2)
    radius_j = fld(km, 2)
    launch_kernel_async!(style, sampled_response_kernel!, scratch, frame, model.kernel, radius_i, radius_j, n, m, kn, km; ndrange=(n, m))
    copyto!(frame, scratch)
    return frame
end

function apply_response!(style::ExecutionStyle, model::RectangularPixelAperture, frame::AbstractMatrix, scratch::AbstractMatrix)
    return _apply_separable_response!(style, frame, scratch, model.kernel_y, model.kernel_x)
end

function apply_response!(style::ExecutionStyle, model::SeparablePixelMTF, frame::AbstractMatrix, scratch::AbstractMatrix)
    return _apply_separable_response!(style, frame, scratch, model.kernel_y, model.kernel_x)
end

const apply_frame_response! = apply_response!

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
    apply_response!(det.params.response_model, det)
    return det.state.frame
end

fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T} = fill_frame!(det, psf, det.params.integration_time)
