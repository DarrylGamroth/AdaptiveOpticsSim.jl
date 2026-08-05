@inline _apply_lift_mapping!(::LiFTIdentityMapping,
    ::LiFTForwardWorkspace, optical_rate::AbstractMatrix) = optical_rate

@inline function _apply_lift_response!(::NullFrameResponse,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    copyto!(workspace.response_buffer, optical_rate)
    return workspace.response_buffer
end

@inline function _apply_lift_response!(response::AbstractFrameResponse,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    copyto!(workspace.response_buffer, optical_rate)
    apply_response!(execution_style(workspace.response_buffer), response,
        workspace.response_buffer, workspace.response_scratch)
    return workspace.response_buffer
end

function _apply_lift_mapping!(mapping::LiFTFrameMapping,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    response_rate = _apply_lift_response!(mapping.response, workspace,
        optical_rate)
    if mapping.sampling > 1
        bin2d!(workspace.sampling_buffer, response_rate, mapping.sampling)
    else
        copyto!(workspace.sampling_buffer, response_rate)
    end
    if mapping.binning > 1
        bin2d!(workspace.mapped_rate_buffer, workspace.sampling_buffer,
            mapping.binning)
    else
        copyto!(workspace.mapped_rate_buffer, workspace.sampling_buffer)
    end
    return workspace.mapped_rate_buffer
end

function _lift_rate_values_from_opd!(forward::PreparedLiFTForward,
    opd::AbstractMatrix; rate_scale::Real=1.0)
    model = forward.plan
    workspace = forward.workspace
    size(opd) == size(model.diversity_opd) || throw(DimensionMismatchError(
        "LiFT forward OPD must match the prepared pupil dimensions"))
    T = eltype(workspace.optical_rate_buffer)
    scale = T(rate_scale)
    isfinite(scale) && scale >= zero(T) || throw(InvalidConfiguration(
        "LiFT forward rate_scale must be finite and nonnegative"))
    amplitude_scale = sqrt(model.photon_irradiance * scale *
        model.pupil_cell_area_m2)
    @. workspace.amplitude_buffer = model.pupil_amplitude * amplitude_scale
    oversampling = focal_field_from_opd!(workspace.focal_buffer, forward,
        workspace.amplitude_buffer, opd)
    field_intensity!(workspace.optical_rate_buffer, workspace.focal_buffer,
        oversampling, workspace.field_scratch)
    maybe_object_convolve!(forward, workspace.optical_rate_buffer)
    return _apply_lift_mapping!(model.mapping, workspace,
        workspace.optical_rate_buffer)
end

"""
    evaluate_lift_forward!(forward; rate_scale=1)

Evaluate the prepared LiFT forward model from its bound pupil-plane OPD input.
The returned `IntensityMap` contains cell-integrated photon-arrival rate on the
prepared observation grid after object convolution and deterministic spatial
preprocessing. No exposure, QE, noise, or readout operation is applied.
"""
function evaluate_lift_forward!(forward::PreparedLiFTForward;
    rate_scale::Real=1.0)
    _require_lift_forward_owner(forward)
    values = _lift_rate_values_from_opd!(forward, forward.input;
        rate_scale=rate_scale)
    copyto!(forward.output.values, values)
    return forward.output
end

"""
    predict_lift_observation!(dest, forward, domain)

Evaluate `forward` at its bound OPD input and write values in the explicitly
selected native observation `domain`. The caller owns `dest`; no detector is
invoked.
"""
function predict_lift_observation!(dest::AbstractMatrix,
    forward::PreparedLiFTForward, domain::AbstractLiFTObservationDomain)
    _require_lift_forward_owner(forward)
    size(dest) == size(forward.output.values) || throw(DimensionMismatchError(
        "LiFT prediction destination must match the observation dimensions"))
    typeof(backend(dest)) === typeof(backend(forward.output.values)) || throw(
        InvalidConfiguration(
            "LiFT prediction destination must use the prepared array backend"))
    compute_device(dest) == compute_device(forward.output.values) || throw(
        InvalidConfiguration(
            "LiFT prediction destination must occupy the prepared compute device"))
    _lift_mightalias_any(dest,
        (forward.input, forward.output.values,
            _lift_forward_workspace_arrays(forward.workspace)...)) && throw(
        InvalidConfiguration(
            "LiFT prediction destination must not alias its prepared owner"))
    rate = _lift_rate_values_from_opd!(forward, forward.input)
    T = eltype(rate)
    native_scale = inv(lift_observation_to_rate_scale(domain, T))
    @. dest = rate * native_scale
    return dest
end

function center_crop!(dest::AbstractMatrix, src::AbstractMatrix)
    if size(dest) == size(src)
        copyto!(dest, src)
        return dest
    end
    n = size(dest, 1)
    cx = div(size(src, 1) - n, 2)
    cy = div(size(src, 2) - n, 2)
    @views copyto!(dest, src[cx+1:cx+n, cy+1:cy+n])
    return dest
end

@inline function maybe_object_convolve!(
    forward::PreparedLiFTForward, matrix::AbstractMatrix)
    return _maybe_object_convolve!(forward.plan.object_kernel,
        forward.workspace, matrix)
end

@inline _maybe_object_convolve!(::Nothing, ::LiFTForwardWorkspace,
    matrix::AbstractMatrix) = matrix

function _maybe_object_convolve!(kernel::LiFTDenseObjectKernel,
    workspace::LiFTForwardWorkspace, matrix::AbstractMatrix)
    conv2d_same!(workspace.convolution_buffer, matrix, kernel.kernel,
        kernel.inv_norm)
    copyto!(matrix, workspace.convolution_buffer)
    return matrix
end

function _maybe_object_convolve!(kernel::LiFTSeparableObjectKernel,
    workspace::LiFTForwardWorkspace, matrix::AbstractMatrix)
    conv2d_same_separable!(workspace.convolution_buffer,
        workspace.convolution_scratch, matrix, kernel.row, kernel.col,
        kernel.inv_norm)
    copyto!(matrix, workspace.convolution_buffer)
    return matrix
end

function _lift_object_kernel(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    inv_norm = inv(T(sum(Array(kernel))))
    row, col = _separable_kernel_factors(kernel)
    if row === nothing
        return LiFTDenseObjectKernel{T,typeof(kernel)}(kernel, inv_norm)
    end
    return LiFTSeparableObjectKernel{T,typeof(row)}(row, col, inv_norm)
end

function _separable_kernel_factors(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    F = svd(Matrix(kernel); full=false)
    isempty(F.S) && return nothing, nothing
    σ1 = F.S[1]
    σ2 = length(F.S) >= 2 ? F.S[2] : zero(T)
    tol = sqrt(eps(T)) * max(one(T), σ1)
    σ2 <= tol || return nothing, nothing
    scale = sqrt(σ1)
    host_row = T.(F.U[:, 1] .* scale)
    host_col = T.(F.V[:, 1] .* scale)
    row = similar(kernel, T, length(host_row))
    col = similar(kernel, T, length(host_col))
    copyto!(row, host_row)
    copyto!(col, host_col)
    return row, col
end

function conv2d_same!(dest::AbstractMatrix{T}, src::AbstractMatrix{T},
    kernel::AbstractMatrix) where {T<:AbstractFloat}
    norm = backend_sum_value(kernel)
    inv_norm = iszero(norm) ? one(T) : inv(T(norm))
    return conv2d_same!(dest, src, kernel, inv_norm)
end

@inline function conv2d_same!(dest::AbstractMatrix{T},
    src::AbstractMatrix{T}, kernel::AbstractMatrix,
    inv_norm::T) where {T<:AbstractFloat}
    return conv2d_same!(execution_style(dest), dest, src, kernel, inv_norm)
end

function conv2d_same!(::ScalarCPUStyle, dest::AbstractMatrix{T},
    src::AbstractMatrix{T}, kernel::AbstractMatrix,
    inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    cx = div(kh, 2)
    cy = div(kw, 2)
    fill!(dest, zero(T))
    @inbounds for ki in 1:kh, kj in 1:kw
        offset_i = ki - cx - 1
        offset_j = kj - cy - 1
        weight = kernel[ki, kj]
        for j in 1:m
            jj = symm_index(j + offset_j, m)
            @simd for i in 1:n
                ii = symm_index(i + offset_i, n)
                dest[i, j] += src[ii, jj] * weight
            end
        end
    end
    @inbounds @simd for index in eachindex(dest)
        dest[index] *= inv_norm
    end
    return dest
end

function conv2d_same!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, src::AbstractMatrix{T},
    kernel::AbstractMatrix, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    launch_kernel!(style, lift_dense_convolution_kernel!, dest, src,
        kernel, inv_norm, n, m, kh, kw; ndrange=(n, m))
    return dest
end

function conv2d_same_separable!(dest::AbstractMatrix{T},
    tmp::AbstractMatrix{T}, src::AbstractMatrix{T},
    row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}) where {T<:AbstractFloat}
    row_norm = backend_sum_value(row_kernel)
    col_norm = backend_sum_value(col_kernel)
    inv_norm = (iszero(row_norm) || iszero(col_norm)) ? one(T) :
        inv(T(row_norm * col_norm))
    return conv2d_same_separable!(dest, tmp, src, row_kernel,
        col_kernel, inv_norm)
end

@inline function conv2d_same_separable!(dest::AbstractMatrix{T},
    tmp::AbstractMatrix{T}, src::AbstractMatrix{T},
    row_kernel::AbstractVector{T}, col_kernel::AbstractVector{T},
    inv_norm::T) where {T<:AbstractFloat}
    return conv2d_same_separable!(execution_style(dest), dest, tmp, src,
        row_kernel, col_kernel, inv_norm)
end

function conv2d_same_separable!(::ScalarCPUStyle,
    dest::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    src::AbstractMatrix{T}, row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kr = length(row_kernel)
    kc = length(col_kernel)
    cx = div(kr, 2)
    cy = div(kc, 2)
    @inbounds for j in 1:m
        for i in 1:n
            acc = zero(T)
            for ki in 1:kr
                ii = symm_index(i + ki - cx - 1, n)
                acc += src[ii, j] * row_kernel[ki]
            end
            tmp[i, j] = acc
        end
    end
    fill!(dest, zero(T))
    @inbounds for kj in 1:kc
        offset_j = kj - cy - 1
        weight = col_kernel[kj]
        for j in 1:m
            jj = symm_index(j + offset_j, m)
            @simd for i in 1:n
                dest[i, j] += tmp[i, jj] * weight
            end
        end
    end
    @inbounds @simd for index in eachindex(dest)
        dest[index] *= inv_norm
    end
    return dest
end

function conv2d_same_separable!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    src::AbstractMatrix{T}, row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, lift_row_convolution_kernel!, tmp, src,
        row_kernel, n, m, length(row_kernel); ndrange=(n, m))
    queue_kernel!(phase, lift_column_convolution_kernel!, dest, tmp,
        col_kernel, inv_norm, n, m, length(col_kernel); ndrange=(n, m))
    finish_kernel_phase!(phase)
    return dest
end

@inline function symm_index(i::Int, n::Int)
    n == 1 && return 1
    while i < 1 || i > n
        if i < 1
            i = 2 - i
        else
            i = 2 * n - i
        end
    end
    return i
end

@inline function lift_oversampling(zero_padding::Int)
    zero_padding < 1 && throw(InvalidConfiguration("LiFT zero_padding must be >= 1"))
    return zero_padding < 2 ? cld(2, zero_padding) : 1
end

@inline function lift_pad_size(resolution::Int, zero_padding::Int)
    oversampling = lift_oversampling(zero_padding)
    nominal = zero_padding * oversampling * resolution
    pad_width = cld(nominal - resolution, 2)
    return resolution + 2 * pad_width
end

function focal_field_from_opd!(dest::AbstractMatrix{Complex{T}},
    forward::PreparedLiFTForward,
    amplitude::AbstractMatrix{T}, opd::AbstractMatrix) where {T<:AbstractFloat}
    model = forward.plan
    n = size(model.pupil_mask, 1)
    oversampling = lift_oversampling(model.zero_padding)
    n_pad = lift_pad_size(n, model.zero_padding)
    image_size = model.focal_resolution * oversampling
    ws = forward.workspace.propagation
    ensure_psf_buffers!(ws, n_pad)
    if size(dest) != (image_size, image_size)
        throw(DimensionMismatchError("LiFT focal field buffer size must match oversampled image size"))
    end

    opd_to_cycles = T(2) / model.wavelength_m
    ox = cld(n_pad - n, 2)
    oy = cld(n_pad - n, 2)
    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = amplitude * cispi(opd_to_cycles * opd)
    if iseven(model.focal_resolution)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(execution_style(ws.pupil_field), ws.pupil_field, phase_shift)
    end

    copyto!(ws.fft_buffer, ws.pupil_field)
    execute_fft_plan!(ws.fft_buffer, ws.fft_plan)
    ws.fft_buffer ./= T(n_pad)

    shift_pix = if n_pad % 2 == image_size % 2
        0
    elseif iseven(n_pad)
        1
    else
        -1
    end
    start = Int(ceil(n_pad / 2)) - div(image_size, 2) + (1 - (n_pad % 2))
    stop = Int(ceil(n_pad / 2)) + div(image_size, 2) + shift_pix
    @views copyto!(dest, ws.fft_buffer[start:stop, start:stop])
    return oversampling
end

function field_intensity!(dest::AbstractMatrix{T}, field::AbstractMatrix{Complex{T}}, oversampling::Int,
    scratch::AbstractMatrix{T}) where {T<:AbstractFloat}
    if oversampling == 1
        @. dest = abs2(field)
        return dest
    end
    n_out, m_out = size(dest)
    if size(field) != (n_out * oversampling, m_out * oversampling)
        throw(DimensionMismatchError("LiFT field size does not match oversampled rate dimensions"))
    end
    if size(scratch) != size(field)
        throw(DimensionMismatchError("LiFT scratch buffer size must match oversampled field size"))
    end
    @. scratch = abs2(field)
    bin2d!(dest, scratch, oversampling)
    return dest
end

function field_derivative!(dest::AbstractMatrix{T}, buf::AbstractMatrix{Complex{T}},
    Pd::AbstractMatrix{Complex{T}}, oversampling::Int, scale::T, scratch::AbstractMatrix{T}) where {T<:AbstractFloat}
    if size(buf) != size(Pd)
        throw(DimensionMismatchError("LiFT focal fields must have matching sizes"))
    end
    if oversampling == 1
        @. dest = scale * real(im * buf * Pd)
        return dest
    end
    n_out, m_out = size(dest)
    if size(buf) != (n_out * oversampling, m_out * oversampling)
        throw(DimensionMismatchError("LiFT derivative field size does not match oversampled image size"))
    end
    if size(scratch) != size(buf)
        throw(DimensionMismatchError("LiFT scratch buffer size must match oversampled derivative field size"))
    end
    @. scratch = real(im * buf * Pd)
    bin2d!(dest, scratch, oversampling)
    @. dest = scale * dest
    return dest
end

function conjugate_field!(dest::AbstractMatrix{Complex{T}}, src::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. dest = conj(src)
    return dest
end
