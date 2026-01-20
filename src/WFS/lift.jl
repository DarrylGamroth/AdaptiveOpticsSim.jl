using LinearAlgebra
using Logging

abstract type LiFTMode end
struct LiFTAnalytic <: LiFTMode end
struct LiFTNumerical <: LiFTMode end

struct LiFTParams{T<:AbstractFloat,A<:AbstractMatrix{T},K}
    diversity_opd::A
    iterations::Int
    img_resolution::Int
    zero_padding::Int
    object_kernel::K
end

mutable struct LiFTState{T<:AbstractFloat,
    W<:Workspace,
    B<:AbstractMatrix{T},
    C<:AbstractMatrix{Complex{T}},
    V<:AbstractVector{T}}
    workspace::W
    psf_buffer::B
    amp_buffer::B
    focal_buffer::C
    mode_buffer::C
    conv_buffer::B
    opd_buffer::B
    residual_buffer::V
    weight_buffer::V
    H_buffer::B
end

struct LiFT{M<:LiFTMode,P<:LiFTParams,S<:LiFTState,B<:AbstractArray{<:AbstractFloat,3},SRC<:AbstractSource,D<:AbstractDetector}
    tel::Telescope
    src::SRC
    det::D
    basis::B
    params::P
    state::S
end

function LiFT(tel::Telescope, src::AbstractSource, basis::AbstractArray, det::AbstractDetector;
    diversity_opd::AbstractMatrix, iterations::Int=5, img_resolution::Int=0,
    numerical::Bool=false, ang_pixel_arcsec=nothing, object_kernel=nothing)

    if size(basis, 1) != tel.params.resolution || size(basis, 2) != tel.params.resolution
        throw(InvalidConfiguration("basis resolution must match telescope resolution"))
    end
    if size(diversity_opd) != size(tel.state.opd)
        throw(InvalidConfiguration("diversity_opd must match telescope resolution"))
    end

    zero_padding = det.params.psf_sampling
    if ang_pixel_arcsec !== nothing
        scale = psf_pixel_scale_arcsec(tel, src, 1)
        zero_padding = max(1, round(Int, scale / ang_pixel_arcsec))
        @info "LiFT using angular sampling override", zero_padding
    end
    if img_resolution <= 0
        img_resolution = tel.params.resolution * zero_padding
    end
    kernel = object_kernel === nothing ? nothing : eltype(tel.state.opd).(object_kernel)
    params = LiFTParams(float.(diversity_opd), iterations, img_resolution, zero_padding, kernel)
    ws = Workspace(tel.state.opd, tel.params.resolution * zero_padding; T=eltype(tel.state.opd))
    psf_buffer = similar(tel.state.opd, eltype(tel.state.opd), img_resolution, img_resolution)
    amp_buffer = similar(tel.state.opd, eltype(tel.state.opd), tel.params.resolution, tel.params.resolution)
    focal_buffer = similar(psf_buffer, Complex{eltype(psf_buffer)})
    mode_buffer = similar(focal_buffer)
    conv_buffer = similar(psf_buffer)
    opd_buffer = similar(amp_buffer)
    residual_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution)
    weight_buffer = similar(residual_buffer)
    H_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution, size(basis, 3))
    state = LiFTState(ws, psf_buffer, amp_buffer, focal_buffer, mode_buffer, conv_buffer,
        opd_buffer, residual_buffer, weight_buffer, H_buffer)
    mode = numerical ? LiFTNumerical() : LiFTAnalytic()
    return LiFT{typeof(mode), typeof(params), typeof(state), typeof(basis), typeof(src), typeof(det)}(
        tel,
        src,
        det,
        basis,
        params,
        state,
    )
end

@inline function prepare_opd!(lift::LiFT, coeffs::AbstractVector)
    combine_basis!(lift.state.opd_buffer, lift.basis, coeffs, lift.tel.state.pupil)
    @. lift.state.opd_buffer += lift.params.diversity_opd
    return lift.state.opd_buffer
end

function lift_interaction_matrix!(H::AbstractMatrix, lift::LiFT{LiFTNumerical}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    psf_size = lift.params.img_resolution
    if size(H, 1) < psf_size * psf_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end
    delta = T(1e-9)

    initial_opd = prepare_opd!(lift, coefficients)
    opd_work = lift.state.amp_buffer
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        @views mode = lift.basis[:, :, mode_id]
        @. opd_work = initial_opd + delta * mode
        psf_p = psf_from_opd!(lift, opd_work; flux_norm=flux_norm)
        copyto!(lift.state.conv_buffer, psf_p)
        @. opd_work = initial_opd - delta * mode
        psf_m = psf_from_opd!(lift, opd_work; flux_norm=flux_norm)
        @. lift.state.psf_buffer = (lift.state.conv_buffer - psf_m) / (2 * delta)
        maybe_object_convolve!(lift, lift.state.psf_buffer)
        @views H[:, idx] .= reshape(lift.state.psf_buffer, :)
    end
    return H
end

function lift_interaction_matrix!(H::AbstractMatrix, lift::LiFT{LiFTAnalytic}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    psf_size = lift.params.img_resolution
    if size(H, 1) < psf_size * psf_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end

    initial_opd = prepare_opd!(lift, coefficients)

    amp_scale = sqrt(T(flux_norm))
    @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale, zero(T))

    focal_field_from_opd!(lift.state.focal_buffer, lift, lift.state.amp_buffer, initial_opd)
    Pd = conj.(lift.state.focal_buffer)

    k = T(2 * pi) / wavelength(lift.src)
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        @views mode = lift.basis[:, :, mode_id]
        @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale * mode, zero(T))
        focal_field_from_opd!(lift.state.mode_buffer, lift, lift.state.amp_buffer, initial_opd)
        @. lift.state.psf_buffer = real(im * lift.state.mode_buffer * Pd)
        lift.state.psf_buffer .*= 2 * k
        maybe_object_convolve!(lift, lift.state.psf_buffer)
        @views H[:, idx] .= reshape(lift.state.psf_buffer, :)
    end
    return H
end

function lift_interaction_matrix(lift::LiFT{LiFTNumerical}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    psf_size = lift.params.img_resolution
    H = zeros(T, psf_size * psf_size, length(mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients, mode_ids; flux_norm=flux_norm)
end

function lift_interaction_matrix(lift::LiFT{LiFTAnalytic}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    psf_size = lift.params.img_resolution
    H = zeros(T, psf_size * psf_size, length(mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients, mode_ids; flux_norm=flux_norm)
end

function reconstruct(lift::LiFT, psf_in::AbstractMatrix, mode_ids::AbstractVector;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)

    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    coeffs = coeffs0 === nothing ? zeros(T, maximum(mode_ids)) : copy(coeffs0)

    residual = lift.state.residual_buffer
    weight = lift.state.weight_buffer
    use_iter_weight = R_n === :model || R_n === :iterative
    if !use_iter_weight
        weight_vector!(weight, psf_in, R_n, lift.det)
    end
    H = @view lift.state.H_buffer[:, 1:n_modes]
    for iter in 1:lift.params.iterations
        current_opd = prepare_opd!(lift, coeffs)
        model_psf = psf_from_opd!(lift, current_opd)
        scale = one(T)
        if optimize_norm == :sum
            denom = sum(model_psf)
            if denom > 0
                scale = sum(psf_in) / denom
            end
        elseif optimize_norm == :max
            denom = maximum(model_psf)
            if denom > 0
                scale = maximum(psf_in) / denom
            end
        end
        if scale != one(T)
            model_psf .*= scale
        end
        @inbounds for idx in eachindex(model_psf)
            residual[idx] = psf_in[idx] - model_psf[idx]
        end
        lift_interaction_matrix!(H, lift, coeffs, mode_ids; flux_norm=scale)

        if use_iter_weight
            weight_vector!(weight, model_psf, R_n, lift.det)
        end
        W = Diagonal(weight)
        normal = H' * W * H
        rhs = H' * W * residual
        delta = try
            normal \ rhs
        catch err
            if err isa SingularException
                pinv(normal) * rhs
            else
                rethrow()
            end
        end
        for (j, mode_id) in enumerate(mode_ids)
            coeffs[mode_id] += delta[j]
        end
        if check_convergence && norm(delta) / max(norm(coeffs), eps(T)) < 1e-3
            break
        end
    end

    return coeffs[mode_ids]
end

function weight_vector!(out::AbstractVector{T}, psf::AbstractMatrix{T}, R_n,
    det::AbstractDetector) where {T<:AbstractFloat}
    σ = readout_noise(det)
    σ2 = T(σ * σ)
    if R_n === :model || R_n === :iterative
        vals = vec(psf)
        @inbounds for idx in eachindex(out)
            denom = vals[idx] + σ2
            out[idx] = inv(max(denom, eps(T)))
        end
        return out
    elseif R_n === nothing
        fill!(out, T(1) / max(σ2, eps(T)))
        return out
    elseif R_n isa AbstractMatrix
        vals = vec(R_n)
        @inbounds for idx in eachindex(out)
            out[idx] = inv(max(vals[idx], eps(T)))
        end
        return out
    end
    throw(InvalidConfiguration("R_n must be :model, :iterative, a matrix, or nothing"))
end

function weight_vector(psf::AbstractMatrix{T}, R_n, det::AbstractDetector) where {T<:AbstractFloat}
    out = similar(psf, T, length(psf))
    return weight_vector!(out, psf, R_n, det)
end

function psf_from_opd!(lift::LiFT, opd::AbstractMatrix; flux_norm::Real=1.0)
    lift.tel.state.opd .= opd
    psf = compute_psf!(lift.tel, lift.src, lift.state.workspace, lift.params.zero_padding)
    psf .*= flux_norm
    center_crop!(lift.state.psf_buffer, psf)
    maybe_object_convolve!(lift, lift.state.psf_buffer)
    return lift.state.psf_buffer
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

function readout_noise(det::Detector{NoiseNone})
    return 0.0
end

function readout_noise(det::Detector{NoisePhoton})
    return 0.0
end

function readout_noise(det::Detector{NoiseReadout})
    return det.noise.sigma
end

function readout_noise(det::Detector{NoisePhotonReadout})
    return det.noise.sigma
end

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,Nothing}}, mat::AbstractMatrix)
    return mat
end

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,<:AbstractMatrix}}, mat::AbstractMatrix)
    conv2d_same!(lift.state.conv_buffer, mat, lift.params.object_kernel)
    copyto!(mat, lift.state.conv_buffer)
    return mat
end

function conv2d_same!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, kernel::AbstractMatrix) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    cx = div(kh, 2)
    cy = div(kw, 2)
    norm = sum(kernel)
    inv_norm = norm == 0 ? one(T) : T(1) / T(norm)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for ki in 1:kh, kj in 1:kw
            ii = symm_index(i + ki - cx - 1, n)
            jj = symm_index(j + kj - cy - 1, m)
            acc += src[ii, jj] * kernel[ki, kj]
        end
        dest[i, j] = acc * inv_norm
    end
    return dest
end

@inline function symm_index(i::Int, n::Int)
    while i < 1 || i > n
        if i < 1
            i = 2 - i
        else
            i = 2 * n - i
        end
    end
    return i
end

function focal_field_from_opd!(dest::AbstractMatrix{Complex{T}}, lift::LiFT,
    amplitude::AbstractMatrix{T}, opd::AbstractMatrix) where {T<:AbstractFloat}
    n = lift.tel.params.resolution
    n_pad = n * lift.params.zero_padding
    ws = lift.state.workspace
    ensure_psf_buffers!(ws, n_pad)

    phase_scale = T(2 * pi) / wavelength(lift.src)
    ox = div(n_pad - n, 2)
    oy = div(n_pad - n, 2)
    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = amplitude * cis(phase_scale * opd)

    copyto!(ws.fft_buffer, ws.pupil_field)
    mul!(ws.fft_buffer, ws.fft_plan, ws.fft_buffer)
    FFTW.fftshift!(ws.pupil_field, ws.fft_buffer)
    return center_crop!(dest, ws.pupil_field)
end
