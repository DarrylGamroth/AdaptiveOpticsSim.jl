using LinearAlgebra
using Logging

struct LiFTParams{T<:AbstractFloat,A<:AbstractMatrix{T}}
    diversity_opd::A
    iterations::Int
    img_resolution::Int
    zero_padding::Int
    numerical::Bool
end

mutable struct LiFTState{T<:AbstractFloat,
    W<:Workspace,
    B<:AbstractMatrix{T},
    C<:AbstractMatrix{Complex{T}}}
    workspace::W
    psf_buffer::B
    amp_buffer::B
    focal_buffer::C
    mode_buffer::C
end

struct LiFT{P<:LiFTParams,S<:LiFTState,B<:AbstractArray{<:AbstractFloat,3},SRC<:AbstractSource,D<:AbstractDetector}
    tel::Telescope
    src::SRC
    det::D
    basis::B
    params::P
    state::S
end

function LiFT(tel::Telescope, src::AbstractSource, basis::AbstractArray, det::AbstractDetector;
    diversity_opd::AbstractMatrix, iterations::Int=5, img_resolution::Int=0,
    numerical::Bool=false, ang_pixel_arcsec=nothing)

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
    params = LiFTParams(float.(diversity_opd), iterations, img_resolution, zero_padding, numerical)
    ws = Workspace(tel.state.opd, tel.params.resolution * zero_padding; T=eltype(tel.state.opd))
    psf_buffer = similar(tel.state.opd, eltype(tel.state.opd), img_resolution, img_resolution)
    amp_buffer = similar(tel.state.opd, eltype(tel.state.opd), tel.params.resolution, tel.params.resolution)
    focal_buffer = similar(psf_buffer, Complex{eltype(psf_buffer)})
    mode_buffer = similar(focal_buffer)
    state = LiFTState(ws, psf_buffer, amp_buffer, focal_buffer, mode_buffer)
    return LiFT(tel, src, det, basis, params, state)
end

function lift_interaction_matrix(lift::LiFT, coefficients::AbstractVector, mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    psf_size = lift.params.img_resolution
    H = zeros(T, psf_size * psf_size, n_modes)
    delta = T(1e-9)

    initial_opd = combine_basis(lift.basis, coefficients, lift.tel.state.pupil) .+ lift.params.diversity_opd
    if lift.params.numerical
        for (idx, mode_id) in enumerate(mode_ids)
            mode = lift.basis[:, :, mode_id]
            psf_p = psf_from_opd!(lift, initial_opd .+ delta .* mode; flux_norm=flux_norm)
            psf_m = psf_from_opd!(lift, initial_opd .- delta .* mode; flux_norm=flux_norm)
            deriv = (psf_p .- psf_m) ./ (2 * delta)
            H[:, idx] .= reshape(deriv, :)
        end
        return H
    end

    fill!(lift.state.amp_buffer, zero(T))
    @inbounds for i in 1:size(lift.state.amp_buffer, 1), j in 1:size(lift.state.amp_buffer, 2)
        lift.state.amp_buffer[i, j] = lift.tel.state.pupil[i, j] ? sqrt(T(flux_norm)) : zero(T)
    end

    focal_field_from_opd!(lift.state.focal_buffer, lift, lift.state.amp_buffer, initial_opd)
    Pd = conj.(lift.state.focal_buffer)

    k = T(2 * pi) / wavelength(lift.src)
    for (idx, mode_id) in enumerate(mode_ids)
        mode = lift.basis[:, :, mode_id]
        @inbounds for i in 1:size(lift.state.amp_buffer, 1), j in 1:size(lift.state.amp_buffer, 2)
            lift.state.amp_buffer[i, j] = lift.tel.state.pupil[i, j] ? sqrt(T(flux_norm)) * mode[i, j] : zero(T)
        end
        focal_field_from_opd!(lift.state.mode_buffer, lift, lift.state.amp_buffer, initial_opd)
        @. lift.state.psf_buffer = real(im * lift.state.mode_buffer * Pd)
        lift.state.psf_buffer .*= 2 * k
        H[:, idx] .= reshape(lift.state.psf_buffer, :)
    end
    return H
end

function reconstruct(lift::LiFT, psf_in::AbstractMatrix, mode_ids::AbstractVector;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)

    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    coeffs = coeffs0 === nothing ? zeros(T, maximum(mode_ids)) : copy(coeffs0)

    for iter in 1:lift.params.iterations
        current_opd = combine_basis(lift.basis, coeffs, lift.tel.state.pupil) .+ lift.params.diversity_opd
        model_psf = psf_from_opd!(lift, current_opd)
        residual = reshape(psf_in .- model_psf, :)
        H = lift_interaction_matrix(lift, coeffs, mode_ids)

        weight = weight_vector(model_psf, R_n)
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

function weight_vector(psf::AbstractMatrix{T}, R_n) where {T<:AbstractFloat}
    if R_n === :model || R_n === :iterative
        return 1 ./ max.(reshape(psf, :), eps(T))
    end
    return ones(T, length(psf))
end

function psf_from_opd!(lift::LiFT, opd::AbstractMatrix; flux_norm::Real=1.0)
    lift.tel.state.opd .= opd
    psf = compute_psf!(lift.tel, lift.src, lift.state.workspace, lift.params.zero_padding)
    psf .*= flux_norm
    return center_crop!(lift.state.psf_buffer, psf)
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
