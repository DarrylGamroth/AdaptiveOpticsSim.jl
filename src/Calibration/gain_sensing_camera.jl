using Logging

mutable struct GSCFFTWorkspace{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    Pf,
    Pi}
    buffer::C
    fft_out::C
    ifft_out::C
    fft_plan::Pf
    ifft_plan::Pi
end

function GSCFFTWorkspace(n::Int; T::Type{<:AbstractFloat}=Float64, backend=Array)
    buffer = backend{Complex{T}}(undef, n, n)
    fft_out = similar(buffer)
    ifft_out = similar(buffer)
    fft_plan = plan_fft_backend!(buffer)
    ifft_plan = plan_ifft_backend!(buffer)
    return GSCFFTWorkspace{T, typeof(buffer), typeof(fft_plan), typeof(ifft_plan)}(
        buffer,
        fft_out,
        ifft_out,
        fft_plan,
        ifft_plan,
    )
end

mutable struct GainSensingCamera{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    B<:AbstractArray{T,3},
    BP<:AbstractArray{Complex{T},3},
    IR<:AbstractMatrix{T},
    SV<:AbstractVector{Complex{T}},
    OG<:AbstractVector{T},
    W<:GSCFFTWorkspace{T}}
    mask::C
    basis::B
    n_modes::Int
    calibration_ready::Bool
    basis_product::Union{Nothing,BP}
    ir_calib::Union{Nothing,IR}
    sensi_calib::Union{Nothing,SV}
    og::OG
    fftws::W
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
    fftws = GSCFFTWorkspace(size(mask_c, 1); T=T)
    bp_probe = similar(fftws.buffer, Complex{T}, size(mask_c, 1), size(mask_c, 2), 1)
    ir_probe = similar(mask_c, T, size(mask_c, 1), size(mask_c, 2))
    sensi_probe = similar(og, Complex{T})
    gsc = GainSensingCamera{T, typeof(mask_c), typeof(basis_t), typeof(bp_probe), typeof(ir_probe), typeof(sensi_probe), typeof(og), typeof(fftws)}(
        mask_c,
        basis_t,
        n_modes,
        false,
        nothing,
        nothing,
        nothing,
        og,
        fftws,
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
    gsc.basis_product = split_basis_product(gsc.basis, gsc.fftws; n_jobs=max(n_jobs, 1))
    gsc.ir_calib = impulse_response(gsc.mask, frame_norm, gsc.fftws)
    gsc.sensi_calib = sensitivity(gsc.ir_calib, gsc.ir_calib, gsc.basis_product, gsc.fftws)
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
    ir_sky = impulse_response(gsc.mask, frame_norm, gsc.fftws)
    sensi_sky = sensitivity(ir_sky, gsc.ir_calib, gsc.basis_product, gsc.fftws)
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

function fft2c!(out::AbstractMatrix{Complex{T}}, ws::GSCFFTWorkspace{T}, input::AbstractMatrix) where {T}
    fftshift2d!(ws.buffer, input)
    mul!(ws.buffer, ws.fft_plan, ws.buffer)
    fftshift2d!(out, ws.buffer)
    out ./= size(out, 1)
    return out
end

function ifft2c!(out::AbstractMatrix{Complex{T}}, ws::GSCFFTWorkspace{T}, input::AbstractMatrix) where {T}
    fftshift2d!(ws.buffer, input)
    mul!(ws.buffer, ws.ifft_plan, ws.buffer)
    fftshift2d!(out, ws.buffer)
    out .*= size(out, 1)
    return out
end

function split_basis_product(basis::AbstractArray{T,3}, ws::GSCFFTWorkspace{T}; n_jobs::Int=10) where {T<:AbstractFloat}
    n_modes = size(basis, 3)
    out = similar(ws.buffer, Complex{T}, size(basis, 1), size(basis, 2), n_modes)
    jobs = max(1, n_jobs)
    if jobs == 1 || Base.Threads.nthreads() == 1
        for k in 1:n_modes
            @views begin
                fft2c!(ws.fft_out, ws, basis[:, :, k])
                ifft2c!(ws.ifft_out, ws, basis[:, :, k])
                out[:, :, k] .= ws.fft_out .* ws.ifft_out
            end
        end
        return out
    end

    ws_list = [GSCFFTWorkspace(size(basis, 1); T=T, backend=typeof(ws.buffer).name.wrapper) for _ in 1:Base.Threads.nthreads()]
    Base.Threads.@threads for k in 1:n_modes
        ws_local = ws_list[Base.Threads.threadid()]
        @views begin
            fft2c!(ws_local.fft_out, ws_local, basis[:, :, k])
            ifft2c!(ws_local.ifft_out, ws_local, basis[:, :, k])
            out[:, :, k] .= ws_local.fft_out .* ws_local.ifft_out
        end
    end
    return out
end

function impulse_response(mask::AbstractMatrix, frame::AbstractMatrix, ws::GSCFFTWorkspace)
    fft2c!(ws.fft_out, ws, mask)
    fft2c!(ws.ifft_out, ws, mask .* frame)
    return 2 .* imag.(conj.(ws.fft_out) .* ws.ifft_out)
end

function sensitivity(IR_1::AbstractMatrix, IR_2::AbstractMatrix, basis_product::AbstractArray, ws::GSCFFTWorkspace)
    fft2c!(ws.fft_out, ws, IR_1)
    fft_IR_1 = vec(ws.fft_out)
    ifft2c!(ws.ifft_out, ws, IR_2)
    ifft_IR_2 = vec(ws.ifft_out)
    product = fft_IR_1 .* ifft_IR_2
    basis_mat = reshape(basis_product, size(IR_1, 1) * size(IR_1, 2), size(basis_product, 3))
    return vec(transpose(product) * basis_mat)
end

function sensitivity(IR_1::AbstractMatrix, IR_2::AbstractMatrix, basis_product::AbstractArray)
    ws = GSCFFTWorkspace(size(IR_1, 1); T=eltype(IR_1))
    return sensitivity(IR_1, IR_2, basis_product, ws)
end
