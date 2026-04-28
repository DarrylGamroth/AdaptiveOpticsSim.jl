using Logging

#
# Gain-sensing camera optical-gain estimation
#
# This file implements the focal-plane gain-sensing camera method used to
# estimate modal optical gains from detector images.
#
# The calibration/evaluation flow is:
# 1. normalize a detector frame by its total flux
# 2. compute the impulse response by mixing the frame with the focal mask and
#    transforming it in Fourier space
# 3. project pairs of impulse responses through precomputed basis products to
#    obtain complex modal sensitivities
# 4. compare on-sky sensitivities to calibration sensitivities to recover
#    per-mode optical gains
#
# Weak modes are explicitly masked so poorly conditioned sensitivity ratios do
# not inject spurious gain estimates.
#
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

"""
    GSCFFTWorkspace(n; T=Float64, backend::AbstractArrayBackend=CPUBackend())
    GSCFFTWorkspace(ref; T=real(eltype(ref)))

Preallocated FFT workspace used by gain-sensing camera calibration and
evaluation.
"""
function GSCFFTWorkspace(n::Int; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
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

"""
    GainSensingCamera(mask, basis; ...)

Construct the gain-sensing camera estimator.

`mask` is the focal-plane complex mask, and `basis[:, :, k]` defines the modal
phase basis used to project the impulse-response sensitivities into modal
optical gains.
"""
mutable struct GainSensingCamera{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    B<:AbstractArray{T,3},
    BP<:AbstractArray{Complex{T},3},
    IR<:AbstractMatrix{T},
    SV<:AbstractVector{Complex{T}},
    OG<:AbstractVector{T},
    WM<:AbstractVector{Bool},
    W<:GSCFFTWorkspace{T}}
    mask::C
    mask_fft::C
    basis::B
    n_modes::Int
    calibration_ready::Bool
    sensitivity_floor::T
    basis_product::Union{Nothing,BP}
    ir_calib::Union{Nothing,IR}
    sensi_calib::Union{Nothing,SV}
    frame_buffer::IR
    ir_buffer::IR
    sensi_buffer::SV
    og::OG
    weak_mode_mask::WM
    fftws::W
    detector_metadata::Union{Nothing,GSCDetectorMetadata{T}}
end

function GainSensingCamera(mask::AbstractMatrix, basis::AbstractArray;
    n_jobs::Int=10, T::Type{<:AbstractFloat}=Float64, detector::Union{Nothing,Detector}=nothing,
    sensitivity_floor::Real=sqrt(eps(T)))
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
    weak_mask = similar(og, Bool)
    fill!(weak_mask, false)
    fftws = GSCFFTWorkspace(mask_c; T=T)
    mask_fft = similar(mask_c)
    fft2c!(mask_fft, fftws, mask_c)
    bp_probe = similar(fftws.buffer, Complex{T}, size(mask_c, 1), size(mask_c, 2), 1)
    ir_probe = similar(mask_c, T, size(mask_c, 1), size(mask_c, 2))
    sensi_probe = similar(og, Complex{T})
    metadata = detector === nothing ? nothing : GSCDetectorMetadata(detector; T=T)
    gsc = GainSensingCamera{T, typeof(mask_c), typeof(basis_t), typeof(bp_probe), typeof(ir_probe), typeof(sensi_probe), typeof(og), typeof(weak_mask), typeof(fftws)}(
        mask_c,
        mask_fft,
        basis_t,
        n_modes,
        false,
        T(sensitivity_floor),
        nothing,
        nothing,
        nothing,
        similar(ir_probe),
        similar(ir_probe),
        similar(sensi_probe),
        og,
        weak_mask,
        fftws,
        metadata,
    )
    if n_jobs < 1
        @warn "n_jobs must be >= 1; using 1 instead."
    end
    return gsc
end

function GSCDetectorMetadata(det::Detector; T::Type{<:AbstractFloat}=eltype(det.state.frame))
    return GSCDetectorMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        det.params.psf_sampling,
        det.params.binning,
        detector_noise_symbol(det.noise),
        detector_readout_sigma(det.noise, det.params.sensor, T),
    )
end

detector_metadata(gsc::GainSensingCamera) = gsc.detector_metadata
weak_mode_mask(gsc::GainSensingCamera) = gsc.weak_mode_mask

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
    print(io, ", sensitivity_floor=", gsc.sensitivity_floor)
    print(io, ")")
end

"""
    calibrate!(gsc, frame; n_jobs=10)

Calibrate the gain-sensing camera from a reference detector frame.

This computes the calibration impulse response, the reference modal
sensitivities, and the weak-mode mask used later when converting sensitivities
into optical-gain estimates.
"""
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
    @. gsc.weak_mode_mask = abs(sensi_calib) <= gsc.sensitivity_floor
    gsc.calibration_ready = true
    fill!(gsc.og, one(eltype(gsc.og)))
    @info "GainSensingCamera calibration complete"
    return gsc
end

function reset_calibration!(gsc::GainSensingCamera)
    gsc.calibration_ready = false
    fill!(gsc.og, one(eltype(gsc.og)))
    fill!(gsc.weak_mode_mask, false)
    return gsc
end

"""
    compute_optical_gains!(gsc, frame)

Estimate the current modal optical gains from a detector frame.

The frame is normalized, converted into an impulse response, projected into
modal sensitivities, and finally divided by the stored calibration
sensitivities.
"""
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
    compute_optical_gains!(execution_style(gsc.og), gsc.og, gsc.sensi_buffer, sensi_calib, gsc.weak_mode_mask)
    return gsc.og
end

function compute_optical_gains!(out::AbstractVector{T}, sensi_sky::AbstractVector, sensi_calib::AbstractVector,
    weak_mode_mask::AbstractVector{Bool}) where {T<:AbstractFloat}
    return compute_optical_gains!(execution_style(out), out, sensi_sky, sensi_calib, weak_mode_mask)
end

function compute_optical_gains!(::ScalarCPUStyle, out::AbstractVector{T}, sensi_sky::AbstractVector, sensi_calib::AbstractVector,
    weak_mode_mask::AbstractVector{Bool}) where {T<:AbstractFloat}
    @inbounds for i in eachindex(out, sensi_sky, sensi_calib, weak_mode_mask)
        out[i] = weak_mode_mask[i] ? one(T) : real(sensi_sky[i] / sensi_calib[i])
    end
    return out
end

function compute_optical_gains!(::AcceleratorStyle, out::AbstractVector{T}, sensi_sky::AbstractVector, sensi_calib::AbstractVector,
    weak_mode_mask::AbstractVector{Bool}) where {T<:AbstractFloat}
    @. out = ifelse(weak_mode_mask, one(T), real(sensi_sky / sensi_calib))
    return out
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
    execute_fft_plan!(ws.buffer, ws.fft_plan)
    fftshift2d!(out, ws.buffer)
    out ./= size(out, 1)
    return out
end

function ifft2c!(out::AbstractMatrix{Complex{T}}, ws::GSCFFTWorkspace{T}, input::AbstractMatrix) where {T}
    fftshift2d!(ws.buffer, input)
    execute_fft_plan!(ws.buffer, ws.ifft_plan)
    fftshift2d!(out, ws.buffer)
    out .*= size(out, 1)
    return out
end

supports_threaded_split_basis(::Array) = true
supports_threaded_split_basis(::AbstractArray) = false

function split_basis_product(basis::AbstractArray{T,3}, ws::GSCFFTWorkspace{T}; n_jobs::Int=10) where {T<:AbstractFloat}
    n_modes = size(basis, 3)
    out = similar(ws.buffer, Complex{T}, size(basis, 1), size(basis, 2), n_modes)
    jobs = max(1, n_jobs)
    if jobs == 1 || Base.Threads.nthreads() == 1 || !supports_threaded_split_basis(ws.buffer)
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

"""
    impulse_response!(out, gsc, frame)

Compute the gain-sensing camera impulse response for a normalized detector
frame.

The implementation mixes the frame with the focal mask and evaluates the
imaginary part of the corresponding Fourier-domain product.
"""
function impulse_response!(out::AbstractMatrix{T}, gsc::GainSensingCamera{T}, frame::AbstractMatrix) where {T<:AbstractFloat}
    @. gsc.fftws.fft_out = gsc.mask * frame
    fft2c!(gsc.fftws.ifft_out, gsc.fftws, gsc.fftws.fft_out)
    @. out = 2 * imag(conj(gsc.mask_fft) * gsc.fftws.ifft_out)
    return out
end

"""
    sensitivity!(out, gsc, IR_1, IR_2)

Project two impulse responses into the complex modal sensitivity domain.

This forms the Fourier-domain product of the two impulse responses and then
projects that product through the precomputed basis products for each mode.
"""
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
