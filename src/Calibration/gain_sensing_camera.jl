using FFTW
using Logging

mutable struct GainSensingCamera{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    B<:AbstractArray{T,3}}
    mask::C
    basis::B
    n_modes::Int
    calibration_ready::Bool
    basis_product::Union{Nothing,Array{Complex{T},3}}
    ir_calib::Union{Nothing,Matrix{T}}
    sensi_calib::Union{Nothing,Vector{Complex{T}}}
    og::Vector{T}
end

function GainSensingCamera(mask::AbstractMatrix, basis::AbstractArray;
    n_jobs::Int=10, T::Type{<:AbstractFloat}=Float64)
    mask_c = Complex{T}.(mask)
    basis_t = T.(basis)
    if size(basis_t, 3) == 1 && ndims(basis_t) == 2
        basis_t = reshape(basis_t, size(basis_t, 1), size(basis_t, 2), 1)
    end
    if size(basis_t, 1) != size(mask_c, 1) || size(basis_t, 2) != size(mask_c, 2)
        basis_t = pad_basis(basis_t, size(mask_c, 1))
    end
    n_modes = size(basis_t, 3)
    og = ones(T, n_modes)
    gsc = GainSensingCamera{T, typeof(mask_c), typeof(basis_t)}(
        mask_c,
        basis_t,
        n_modes,
        false,
        nothing,
        nothing,
        nothing,
        og,
    )
    if n_jobs < 1
        @warn "n_jobs must be >= 1; using 1 instead."
    end
    return gsc
end

function calibrate!(gsc::GainSensingCamera, frame::AbstractMatrix; n_jobs::Int=10)
    total = sum(frame)
    if total <= 0
        throw(InvalidConfiguration("frame sum must be > 0 for calibration"))
    end
    @info "GainSensingCamera calibration started"
    frame_norm = frame ./ total
    gsc.basis_product = split_basis_product(gsc.basis; n_jobs=max(n_jobs, 1))
    gsc.ir_calib = impulse_response(gsc.mask, frame_norm)
    gsc.sensi_calib = sensitivity(gsc.ir_calib, gsc.ir_calib, gsc.basis_product)
    gsc.calibration_ready = true
    fill!(gsc.og, one(eltype(gsc.og)))
    @info "GainSensingCamera calibration complete"
    return gsc
end

function reset_calibration!(gsc::GainSensingCamera)
    gsc.calibration_ready = false
    fill!(gsc.og, one(eltype(gsc.og)))
    return gsc
end

function compute_optical_gains!(gsc::GainSensingCamera, frame::AbstractMatrix)
    if !gsc.calibration_ready
        throw(InvalidConfiguration("GainSensingCamera must be calibrated before computing optical gains"))
    end
    total = sum(frame)
    if total <= 0
        throw(InvalidConfiguration("frame sum must be > 0 for optical gain computation"))
    end
    frame_norm = frame ./ total
    ir_sky = impulse_response(gsc.mask, frame_norm)
    sensi_sky = sensitivity(ir_sky, gsc.ir_calib, gsc.basis_product)
    @. gsc.og = real(sensi_sky / gsc.sensi_calib)
    return gsc.og
end

function pad_basis(basis::AbstractArray{T,3}, size_out::Int) where {T}
    if size_out <= size(basis, 1)
        throw(InvalidConfiguration("size_out must be larger than basis size"))
    end
    pad = div(size_out - size(basis, 1), 2)
    out = zeros(T, size_out, size_out, size(basis, 3))
    @views out[pad+1:pad+size(basis, 1), pad+1:pad+size(basis, 2), :] .= basis
    return out
end

function fft2c(mat::AbstractMatrix)
    shifted = FFTW.fftshift(mat, (1, 2))
    return FFTW.fftshift(fft(shifted), (1, 2)) / size(mat, 1)
end

function ifft2c(mat::AbstractMatrix)
    shifted = FFTW.fftshift(mat, (1, 2))
    return FFTW.fftshift(ifft(shifted), (1, 2)) * size(mat, 1)
end

function split_basis_product(basis::AbstractArray{T,3}; n_jobs::Int=10) where {T<:AbstractFloat}
    n_modes = size(basis, 3)
    n_chunks = ceil(Int, n_modes / max(n_jobs, 1))
    out = Array{Complex{T}}(undef, size(basis, 1), size(basis, 2), n_modes)
    idx = 1
    for _ in 1:n_chunks
        last = min(idx + n_jobs - 1, n_modes)
        for k in idx:last
            @views begin
                fft_b = fft2c(basis[:, :, k])
                ifft_b = ifft2c(basis[:, :, k])
                out[:, :, k] .= fft_b .* ifft_b
            end
        end
        idx = last + 1
    end
    return out
end

function impulse_response(mask::AbstractMatrix, frame::AbstractMatrix)
    fft_mask = fft2c(mask)
    fft_field = fft2c(mask .* frame)
    return 2 .* imag.(conj.(fft_mask) .* fft_field)
end

function sensitivity(IR_1::AbstractMatrix, IR_2::AbstractMatrix, basis_product::AbstractArray)
    fft_IR_1 = vec(fft2c(IR_1))
    ifft_IR_2 = vec(ifft2c(IR_2))
    product = fft_IR_1 .* ifft_IR_2
    basis_mat = reshape(basis_product, size(IR_1, 1) * size(IR_1, 2), size(basis_product, 3))
    return vec(transpose(product) * basis_mat)
end
