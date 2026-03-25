#
# Curvature wavefront sensing
#
# This sensor is modeled as a twin-defocus propagated pupil measurement:
#
# 1. embed the complex pupil field on a padded FFT grid
# 2. apply equal-and-opposite defocus masks to form two branches
# 3. propagate both branches with a batched centered FFT
# 4. crop back to the pupil support and bin to the exported readout
# 5. form the normalized signal `(I⁺ - I⁻) / (I⁺ + I⁻)`
#
# The exported `state.slopes` remains one scalar signal per valid subaperture
# so the runtime and reconstructor interfaces stay compatible with the rest of
# the package.
#

@kernel function curvature_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds phasor[i, j] = cis(scale * (i + j - 2))
    end
end

@kernel function curvature_signal_from_frame_kernel!(signal_2d, slopes, frame, reference_signal_2d, valid_mask,
    epsval, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        if @inbounds valid_mask[i, j]
            plus = @inbounds frame[i, j]
            minus = @inbounds frame[i + n_sub, j]
            signal = (plus - minus) / (plus + minus + epsval)
            corrected = signal - @inbounds(reference_signal_2d[i, j])
            @inbounds begin
                signal_2d[i, j] = corrected
                slopes[idx] = corrected
            end
        else
            @inbounds begin
                signal_2d[i, j] = zero(eltype(signal_2d))
                slopes[idx] = zero(eltype(slopes))
            end
        end
    end
end

@kernel function curvature_branch_field_stack_kernel!(field_stack, pupil, opd, defocus_stack, phasor,
    amp_scale, opd_to_cycles, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= size(field_stack, 3)
        xi = x - ox
        yi = y - oy
        val = zero(eltype(field_stack))
        if 1 <= xi <= n && 1 <= yi <= n && @inbounds(pupil[xi, yi])
            val = amp_scale * cispi(opd_to_cycles * @inbounds(opd[xi, yi]))
        end
        @inbounds field_stack[x, y, branch] = val * defocus_stack[x, y, branch] * phasor[x, y]
    end
end

@kernel function curvature_abs2_stack_kernel!(intensity_stack, field_stack, pad::Int, n_branches::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= n_branches
        @inbounds intensity_stack[x, y, branch] = abs2(field_stack[x, y, branch])
    end
end

"""
    CurvatureWFS

Curvature wavefront sensor using propagated intra-/extra-focal detector planes.

The sensor generates two defocused propagated intensity maps and exports them as
an `(2n_subap, n_subap)` detector-like frame. The 1-D `state.slopes` vector is
the normalized curvature signal sampled on the valid subaperture grid.
"""
struct CurvatureWFSParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    response_gain::T
    diffraction_padding::Int
end

mutable struct CurvatureWFSState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    CS<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pfs,
}
    valid_mask::A
    slopes::V
    phasor::C
    field_stack::CS
    defocus_stack::CS
    intensity_stack::RS
    cropped_plus::R
    cropped_minus::R
    frame_plus::R
    frame_minus::R
    signal_2d::R
    reference_signal_2d::R
    camera_frame::R
    fft_stack_plan::Pfs
    effective_padding::Int
    calibrated::Bool
    calibration_wavelength::T
end

struct CurvatureWFS{P<:CurvatureWFSParams,S<:CurvatureWFSState} <: AbstractWFS
    params::P
    state::S
end

"""
    CurvatureWFS(tel; n_subap, threshold=0.1, response_gain=0.5, diffraction_padding=2, ...)

Construct a curvature WFS using two propagated defocus branches.

`response_gain` scales the RMS defocus phase amplitude. In the current optical
model, `1.0` corresponds to an approximately `1000 nm` RMS defocus term on the
pupil.
"""
function CurvatureWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1, response_gain::Real=0.5,
    diffraction_padding::Int=2, T::Type{<:AbstractFloat}=Float64, backend=Array)
    tel.params.resolution % n_subap == 0 ||
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    diffraction_padding >= 1 ||
        throw(InvalidConfiguration("diffraction_padding must be >= 1"))
    n = tel.params.resolution
    pad = n * diffraction_padding
    params = CurvatureWFSParams{T}(n_subap, T(threshold), T(response_gain), diffraction_padding)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, n_subap * n_subap)
    n_branches = 2
    phasor = backend{Complex{T}}(undef, pad, pad)
    field_stack = backend{Complex{T}}(undef, pad, pad, n_branches)
    defocus_stack = similar(field_stack)
    intensity_stack = backend{T}(undef, pad, pad, n_branches)
    cropped_plus = backend{T}(undef, n, n)
    cropped_minus = backend{T}(undef, n, n)
    frame_plus = backend{T}(undef, n_subap, n_subap)
    frame_minus = backend{T}(undef, n_subap, n_subap)
    signal_2d = backend{T}(undef, n_subap, n_subap)
    reference_signal_2d = similar(signal_2d)
    camera_frame = backend{T}(undef, 2 * n_subap, n_subap)
    fill!(slopes, zero(T))
    fill!(field_stack, zero(eltype(field_stack)))
    fill!(intensity_stack, zero(T))
    fill!(cropped_plus, zero(T))
    fill!(cropped_minus, zero(T))
    fill!(frame_plus, zero(T))
    fill!(frame_minus, zero(T))
    fill!(signal_2d, zero(T))
    fill!(reference_signal_2d, zero(T))
    fill!(camera_frame, zero(T))
    fft_stack_plan = plan_fft_backend!(field_stack, (1, 2))
    state = CurvatureWFSState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(phasor),
        typeof(field_stack),
        typeof(frame_plus),
        typeof(intensity_stack),
        typeof(fft_stack_plan),
    }(
        valid_mask,
        slopes,
        phasor,
        field_stack,
        defocus_stack,
        intensity_stack,
        cropped_plus,
        cropped_minus,
        frame_plus,
        frame_minus,
        signal_2d,
        reference_signal_2d,
        camera_frame,
        fft_stack_plan,
        diffraction_padding,
        false,
        zero(T),
    )
    wfs = CurvatureWFS{typeof(params),typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    build_curvature_phasor!(wfs.state.phasor)
    build_curvature_defocus_masks!(wfs, tel)
    return wfs
end

sensing_mode(::CurvatureWFS) = Diffractive()

function update_valid_mask!(wfs::CurvatureWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function ensure_curvature_buffers!(wfs::CurvatureWFS, tel::Telescope)
    n = tel.params.resolution
    pad = n * wfs.params.diffraction_padding
    if size(wfs.state.field_stack, 1) != pad || size(wfs.state.field_stack, 2) != pad
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        n_branches = size(wfs.state.field_stack, 3)
        wfs.state.field_stack = similar(wfs.state.field_stack, pad, pad, n_branches)
        wfs.state.defocus_stack = similar(wfs.state.defocus_stack, pad, pad, n_branches)
        wfs.state.intensity_stack = similar(wfs.state.intensity_stack, pad, pad, n_branches)
        wfs.state.cropped_plus = similar(wfs.state.cropped_plus, n, n)
        wfs.state.cropped_minus = similar(wfs.state.cropped_minus, n, n)
        wfs.state.fft_stack_plan = plan_fft_backend!(wfs.state.field_stack, (1, 2))
        wfs.state.effective_padding = wfs.params.diffraction_padding
        wfs.state.calibrated = false
        build_curvature_phasor!(wfs.state.phasor)
        build_curvature_defocus_masks!(wfs, tel)
    end
    return wfs
end

function build_curvature_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    return build_curvature_phasor!(execution_style(phasor), phasor)
end

function build_curvature_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for j in 1:n, i in 1:n
        phasor[i, j] = cis(scale * (i + j - 2))
    end
    return phasor
end

function build_curvature_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, curvature_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function host_curvature_defocus_stack(wfs::CurvatureWFS, tel::Telescope)
    T = eltype(wfs.state.frame_plus)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    cx = (pad + 1) / 2
    radius = T(n) / 2
    defocus_rms_nm = T(1000) * wfs.params.response_gain
    out = Array{Complex{T}}(undef, pad, pad, 2)
    inv_waves = defocus_rms_nm / T(1e9)
    @inbounds for j in 1:pad, i in 1:pad
        x = (T(i) - T(cx)) / radius
        y = (T(j) - T(cx)) / radius
        r2 = x * x + y * y
        z4 = r2 <= one(T) ? sqrt(T(3)) * (T(2) * r2 - one(T)) : zero(T)
        phase_cycles = inv_waves * z4
        out[i, j, 1] = cispi(T(2) * phase_cycles)
        out[i, j, 2] = cispi(-T(2) * phase_cycles)
    end
    return out
end

function build_curvature_defocus_masks!(wfs::CurvatureWFS, tel::Telescope)
    copyto!(wfs.state.defocus_stack, host_curvature_defocus_stack(wfs, tel))
    return wfs
end

function sample_curvature_frames!(wfs::CurvatureWFS, tel::Telescope)
    n = tel.params.resolution
    pad = size(wfs.state.intensity_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    sub = div(n, wfs.params.n_subap)
    copyto!(wfs.state.cropped_plus, @view(wfs.state.intensity_stack[ox+1:ox+n, oy+1:oy+n, 1]))
    copyto!(wfs.state.cropped_minus, @view(wfs.state.intensity_stack[ox+1:ox+n, oy+1:oy+n, 2]))
    bin2d!(wfs.state.frame_plus, wfs.state.cropped_plus, sub)
    bin2d!(wfs.state.frame_minus, wfs.state.cropped_minus, sub)
    @views copyto!(wfs.state.camera_frame[1:wfs.params.n_subap, :], wfs.state.frame_plus)
    @views copyto!(wfs.state.camera_frame[wfs.params.n_subap+1:2*wfs.params.n_subap, :], wfs.state.frame_minus)
    return wfs.state.camera_frame
end

function curvature_intensity!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    return curvature_intensity!(execution_style(wfs.state.camera_frame), wfs, tel, src)
end

function curvature_intensity!(::ScalarCPUStyle, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.frame_plus)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))
    fill!(wfs.state.field_stack, zero(eltype(wfs.state.field_stack)))
    @inbounds for y in 1:n, x in 1:n
        if tel.state.pupil[x, y]
            val = amp_scale * cispi(opd_to_cycles * tel.state.opd[x, y])
            xx = ox + x
            yy = oy + y
            common = val * wfs.state.phasor[xx, yy]
            wfs.state.field_stack[xx, yy, 1] = common * wfs.state.defocus_stack[xx, yy, 1]
            wfs.state.field_stack[xx, yy, 2] = common * wfs.state.defocus_stack[xx, yy, 2]
        end
    end
    execute_fft_plan!(wfs.state.field_stack, wfs.state.fft_stack_plan)
    @. wfs.state.intensity_stack = abs2(wfs.state.field_stack)
    return sample_curvature_frames!(wfs, tel)
end

function curvature_intensity!(style::AcceleratorStyle, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.frame_plus)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))
    launch_kernel_async!(style, curvature_branch_field_stack_kernel!, wfs.state.field_stack, tel.state.pupil,
        tel.state.opd, wfs.state.defocus_stack, wfs.state.phasor, amp_scale, opd_to_cycles, ox, oy, n, pad;
        ndrange=size(wfs.state.field_stack))
    synchronize_backend!(style)
    execute_fft_plan!(wfs.state.field_stack, wfs.state.fft_stack_plan)
    launch_kernel!(style, curvature_abs2_stack_kernel!, wfs.state.intensity_stack, wfs.state.field_stack,
        pad, size(wfs.state.field_stack, 3); ndrange=size(wfs.state.intensity_stack))
    return sample_curvature_frames!(wfs, tel)
end

function curvature_signal!(wfs::CurvatureWFS, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(wfs.state.camera_frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled camera frame"))
    n_sub = wfs.params.n_subap
    epsval = eps(eltype(wfs.state.signal_2d))
    style = execution_style(frame)
    if style isa ScalarCPUStyle
        @inbounds for i in 1:n_sub, j in 1:n_sub
            idx = (i - 1) * n_sub + j
            if wfs.state.valid_mask[i, j]
                plus = frame[i, j]
                minus = frame[i + n_sub, j]
                signal = (plus - minus) / (plus + minus + epsval)
                corrected = signal - wfs.state.reference_signal_2d[i, j]
                wfs.state.signal_2d[i, j] = corrected
                wfs.state.slopes[idx] = corrected
            else
                wfs.state.signal_2d[i, j] = zero(eltype(wfs.state.signal_2d))
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            end
        end
    else
        launch_kernel!(style, curvature_signal_from_frame_kernel!, wfs.state.signal_2d, wfs.state.slopes,
            frame, wfs.state.reference_signal_2d, wfs.state.valid_mask, epsval, n_sub;
            ndrange=size(wfs.state.signal_2d))
    end
    return wfs.state.slopes
end

function ensure_curvature_calibration!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    if wfs.state.calibrated && wfs.state.calibration_wavelength == src.params.wavelength
        return wfs
    end
    opd_saved = copy(tel.state.opd)
    reset_opd!(tel)
    update_valid_mask!(wfs, tel)
    curvature_intensity!(wfs, tel, src)
    curvature_signal!(wfs, wfs.state.camera_frame)
    copyto!(wfs.state.reference_signal_2d, wfs.state.signal_2d)
    copyto!(tel.state.opd, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = src.params.wavelength
    return wfs
end

function measure!(::Diffractive, wfs::CurvatureWFS, tel::Telescope)
    throw(InvalidConfiguration("CurvatureWFS requires a source; call measure!(wfs, tel, src)."))
end

measure!(wfs::CurvatureWFS, tel::Telescope) = measure!(sensing_mode(wfs), wfs, tel)
measure!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource) = measure!(sensing_mode(wfs), wfs, tel, src)

function measure!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::CurvatureWFS, tel::Telescope, ast::Asterism)
    throw(InvalidConfiguration("CurvatureWFS asterism support is not implemented"))
end

function measure!(wfs::CurvatureWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration("CurvatureWFS asterism support is not implemented"))
end

function measure!(::Diffractive, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src)
    return curvature_signal!(wfs, wfs.state.camera_frame)
end

function measure!(::Diffractive, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src)
    capture!(det, wfs.state.camera_frame; rng=rng)
    size(output_frame(det)) == size(wfs.state.camera_frame) ||
        throw(InvalidConfiguration("CurvatureWFS detector output size must match the sampled camera frame"))
    return curvature_signal!(wfs, output_frame(det))
end
