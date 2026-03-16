using Logging

struct GSCDetectorMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    noise::Symbol
    readout_sigma::Union{Nothing,T}
end

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

function GSCFFTWorkspace(ref::AbstractMatrix; T::Type{<:AbstractFloat}=real(eltype(ref)))
    n = size(ref, 1)
    buffer = similar(ref, Complex{T}, n, n)
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
    mask_fft::C
    basis::B
    n_modes::Int
    calibration_ready::Bool
    basis_product::Union{Nothing,BP}
    ir_calib::Union{Nothing,IR}
    sensi_calib::Union{Nothing,SV}
    frame_buffer::IR
    ir_buffer::IR
    sensi_buffer::SV
    og::OG
    fftws::W
    detector_metadata::Union{Nothing,GSCDetectorMetadata{T}}
end

function GainSensingCamera(mask::AbstractMatrix, basis::AbstractArray;
    n_jobs::Int=10, T::Type{<:AbstractFloat}=Float64, detector::Union{Nothing,Detector}=nothing)
    mask_c = Complex{T}.(mask)
    basis_t = T.(basis)
    if size(basis_t, 3) == 1 && ndims(basis_t) == 2
        basis_t = reshape(basis_t, size(basis_t, 1), size(basis_t, 2), 1)
    end
    if size(basis_t, 1) != size(mask_c, 1) || size(basis_t, 2) != size(mask_c, 2)
        basis_t = pad_basis(basis_t, size(mask_c, 1))
    end
    n_modes = size(basis_t, 3)
    og = similar(vec(mask_c), T, n_modes)
    fill!(og, one(T))
    fftws = GSCFFTWorkspace(mask_c; T=T)
    mask_fft = similar(mask_c)
    fft2c!(mask_fft, fftws, mask_c)
    bp_probe = similar(fftws.buffer, Complex{T}, size(mask_c, 1), size(mask_c, 2), 1)
    ir_probe = similar(mask_c, T, size(mask_c, 1), size(mask_c, 2))
    sensi_probe = similar(og, Complex{T})
    metadata = detector === nothing ? nothing : GSCDetectorMetadata(detector; T=T)
    gsc = GainSensingCamera{T, typeof(mask_c), typeof(basis_t), typeof(bp_probe), typeof(ir_probe), typeof(sensi_probe), typeof(og), typeof(fftws)}(
        mask_c,
        mask_fft,
        basis_t,
        n_modes,
        false,
        nothing,
        nothing,
        nothing,
        similar(ir_probe),
        similar(ir_probe),
        similar(sensi_probe),
        og,
        fftws,
        metadata,
    )
    if n_jobs < 1
        @warn "n_jobs must be >= 1; using 1 instead."
    end
    return gsc
end

noise_symbol(::NoiseNone) = :none
noise_symbol(::NoisePhoton) = :photon
noise_symbol(::NoiseReadout) = :readout
noise_symbol(::NoisePhotonReadout) = :photon_readout

readout_sigma(::NoiseNone, ::Type{T}) where {T<:AbstractFloat} = nothing
readout_sigma(::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} = nothing
readout_sigma(noise::NoiseReadout, ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)
readout_sigma(noise::NoisePhotonReadout, ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)

function GSCDetectorMetadata(det::Detector; T::Type{<:AbstractFloat}=eltype(det.state.frame))
    return GSCDetectorMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        det.params.psf_sampling,
        det.params.binning,
        noise_symbol(det.noise),
        readout_sigma(det.noise, T),
    )
end

detector_metadata(gsc::GainSensingCamera) = gsc.detector_metadata

function attach_detector!(gsc::GainSensingCamera{T}, det::Detector) where {T<:AbstractFloat}
    gsc.detector_metadata = GSCDetectorMetadata(det; T=T)
    return gsc
end

function detach_detector!(gsc::GainSensingCamera)
    gsc.detector_metadata = nothing
    return gsc
end

function Base.show(io::IO, ::MIME"text/plain", gsc::GainSensingCamera)
    print(io, "GainSensingCamera(")
    print(io, "calibrated=", gsc.calibration_ready)
    print(io, ", n_modes=", gsc.n_modes)
    metadata = gsc.detector_metadata
    if metadata !== nothing
        print(io, ", detector=(")
        print(io, "integration_time=", metadata.integration_time)
        print(io, ", qe=", metadata.qe)
        print(io, ", psf_sampling=", metadata.psf_sampling)
        print(io, ", binning=", metadata.binning)
        print(io, ", noise=", metadata.noise)
        if metadata.readout_sigma !== nothing
            print(io, ", readout_sigma=", metadata.readout_sigma)
        end
        print(io, ")")
    end
    print(io, ")")
end

function calibrate!(gsc::GainSensingCamera, frame::AbstractMatrix; n_jobs::Int=10)
    total = sum(frame)
    if total <= 0
        throw(InvalidConfiguration("frame sum must be > 0 for calibration"))
    end
    @info "GainSensingCamera calibration started"
    normalize_frame!(gsc.frame_buffer, frame, total)
    gsc.basis_product = split_basis_product(gsc.basis, gsc.fftws; n_jobs=max(n_jobs, 1))
    ir_calib = gsc.ir_calib === nothing ? similar(gsc.ir_buffer) : gsc.ir_calib
    sensi_calib = gsc.sensi_calib === nothing ? similar(gsc.sensi_buffer) : gsc.sensi_calib
    impulse_response!(ir_calib, gsc, gsc.frame_buffer)
    sensitivity!(sensi_calib, gsc, ir_calib, ir_calib)
    gsc.ir_calib = ir_calib
    gsc.sensi_calib = sensi_calib
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
    ir_calib = gsc.ir_calib
    ir_calib === nothing && throw(InvalidConfiguration("GainSensingCamera calibration IR is not available"))
    sensi_calib = gsc.sensi_calib
    sensi_calib === nothing && throw(InvalidConfiguration("GainSensingCamera calibration sensitivity is not available"))
    total = sum(frame)
    if total <= 0
        throw(InvalidConfiguration("frame sum must be > 0 for optical gain computation"))
    end
    normalize_frame!(gsc.frame_buffer, frame, total)
    impulse_response!(gsc.ir_buffer, gsc, gsc.frame_buffer)
    sensitivity!(gsc.sensi_buffer, gsc, gsc.ir_buffer, ir_calib)
    @. gsc.og = real(gsc.sensi_buffer / sensi_calib)
    return gsc.og
end

function normalize_frame!(out::AbstractMatrix{T}, frame::AbstractMatrix, total::Real) where {T<:AbstractFloat}
    inv_total = inv(T(total))
    @. out = frame * inv_total
    return out
end

function pad_basis(basis::AbstractArray{T,3}, size_out::Int) where {T}
    if size_out <= size(basis, 1)
        throw(InvalidConfiguration("size_out must be larger than basis size"))
    end
    pad = div(size_out - size(basis, 1), 2)
    out = similar(basis, T, size_out, size_out, size(basis, 3))
    fill!(out, zero(T))
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
    if jobs == 1 || Base.Threads.nthreads() == 1 || !(ws.buffer isa Array)
        for k in 1:n_modes
            @views begin
                fft2c!(ws.fft_out, ws, basis[:, :, k])
                ifft2c!(ws.ifft_out, ws, basis[:, :, k])
                out[:, :, k] .= ws.fft_out .* ws.ifft_out
            end
        end
        return out
    end

    ws_list = [GSCFFTWorkspace(ws.buffer; T=T) for _ in 1:Base.Threads.nthreads()]
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

function impulse_response!(out::AbstractMatrix{T}, gsc::GainSensingCamera{T}, frame::AbstractMatrix) where {T<:AbstractFloat}
    @. gsc.fftws.fft_out = gsc.mask * frame
    fft2c!(gsc.fftws.ifft_out, gsc.fftws, gsc.fftws.fft_out)
    @. out = 2 * imag(conj(gsc.mask_fft) * gsc.fftws.ifft_out)
    return out
end

function sensitivity!(out::AbstractVector{Complex{T}}, gsc::GainSensingCamera{T},
    IR_1::AbstractMatrix, IR_2::AbstractMatrix) where {T<:AbstractFloat}
    fft2c!(gsc.fftws.fft_out, gsc.fftws, IR_1)
    ifft2c!(gsc.fftws.ifft_out, gsc.fftws, IR_2)
    basis_product = gsc.basis_product
    basis_product === nothing && throw(InvalidConfiguration("GainSensingCamera basis product is not available"))
    basis_mat = reshape(basis_product, :, size(basis_product, 3))
    @. gsc.fftws.buffer = gsc.fftws.fft_out * gsc.fftws.ifft_out
    mul!(out, adjoint(basis_mat), vec(gsc.fftws.buffer))
    return out
end

function impulse_response(mask::AbstractMatrix, frame::AbstractMatrix, ws::GSCFFTWorkspace)
    mask_fft = similar(frame, Complex{real(eltype(ws.buffer))}, size(frame)...)
    fft2c!(mask_fft, ws, mask)
    tmp = similar(ws.buffer)
    @. tmp = mask * frame
    fft2c!(ws.ifft_out, ws, tmp)
    return 2 .* imag.(conj.(mask_fft) .* ws.ifft_out)
end

function sensitivity(IR_1::AbstractMatrix, IR_2::AbstractMatrix, basis_product::AbstractArray, ws::GSCFFTWorkspace)
    fft2c!(ws.fft_out, ws, IR_1)
    ifft2c!(ws.ifft_out, ws, IR_2)
    basis_mat = reshape(basis_product, :, size(basis_product, 3))
    out = similar(vec(IR_1), Complex{real(eltype(ws.buffer))}, size(basis_product, 3))
    @. ws.buffer = ws.fft_out * ws.ifft_out
    mul!(out, adjoint(basis_mat), vec(ws.buffer))
    return out
end

function sensitivity(IR_1::AbstractMatrix, IR_2::AbstractMatrix, basis_product::AbstractArray)
    ws = GSCFFTWorkspace(size(IR_1, 1); T=eltype(IR_1))
    return sensitivity(IR_1, IR_2, basis_product, ws)
end
