#
# Curvature wavefront sensing
#
# This implementation follows the same broad optical structure used by
# SPECULA's curvature sensor:
#
# 1. embed the pupil field on a padded FFT grid
# 2. apply equal-and-opposite defocus phase masks in the pupil plane
# 3. propagate each branch to a detector plane with a centered FFT
# 4. crop the central detector region back to the pupil scale
# 5. bin to the exported `n_subap x n_subap` readout
# 6. form the normalized curvature signal `(I⁺ - I⁻) / (I⁺ + I⁻)`
#
# The current exported `state.slopes` remains one scalar signal per valid
# subaperture so the runtime/reconstructor surface stays compatible with the
# rest of the package, but the signal now comes from propagated twin-defocus
# images rather than a discrete Laplacian surrogate.
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
    R<:AbstractMatrix{T},
    Pf,
}
    valid_mask::A
    slopes::V
    field::C
    focal_field::C
    phasor::C
    defocus_plus::C
    defocus_minus::C
    intensity_plus::R
    intensity_minus::R
    cropped_plus::R
    cropped_minus::R
    frame_plus::R
    frame_minus::R
    signal_2d::R
    reference_signal_2d::R
    camera_frame::R
    fft_plan::Pf
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

`response_gain` is a legacy high-level control on the defocus strength. In the
current optical model it scales the RMS defocus phase amplitude, with `1.0`
corresponding to an approximately `1000 nm` RMS defocus term on the pupil.
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
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    phasor = similar(field)
    defocus_plus = similar(field)
    defocus_minus = similar(field)
    intensity_plus = backend{T}(undef, pad, pad)
    intensity_minus = backend{T}(undef, pad, pad)
    cropped_plus = backend{T}(undef, n, n)
    cropped_minus = backend{T}(undef, n, n)
    frame_plus = backend{T}(undef, n_subap, n_subap)
    frame_minus = backend{T}(undef, n_subap, n_subap)
    signal_2d = backend{T}(undef, n_subap, n_subap)
    reference_signal_2d = similar(signal_2d)
    camera_frame = backend{T}(undef, 2 * n_subap, n_subap)
    fill!(slopes, zero(T))
    fill!(field, zero(eltype(field)))
    fill!(focal_field, zero(eltype(focal_field)))
    fill!(intensity_plus, zero(T))
    fill!(intensity_minus, zero(T))
    fill!(cropped_plus, zero(T))
    fill!(cropped_minus, zero(T))
    fill!(frame_plus, zero(T))
    fill!(frame_minus, zero(T))
    fill!(signal_2d, zero(T))
    fill!(reference_signal_2d, zero(T))
    fill!(camera_frame, zero(T))
    fft_plan = plan_fft_backend!(focal_field)
    state = CurvatureWFSState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(field),
        typeof(frame_plus),
        typeof(fft_plan),
    }(
        valid_mask,
        slopes,
        field,
        focal_field,
        phasor,
        defocus_plus,
        defocus_minus,
        intensity_plus,
        intensity_minus,
        cropped_plus,
        cropped_minus,
        frame_plus,
        frame_minus,
        signal_2d,
        reference_signal_2d,
        camera_frame,
        fft_plan,
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
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.focal_field = similar(wfs.state.focal_field, pad, pad)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.defocus_plus = similar(wfs.state.defocus_plus, pad, pad)
        wfs.state.defocus_minus = similar(wfs.state.defocus_minus, pad, pad)
        wfs.state.intensity_plus = similar(wfs.state.intensity_plus, pad, pad)
        wfs.state.intensity_minus = similar(wfs.state.intensity_minus, pad, pad)
        wfs.state.cropped_plus = similar(wfs.state.cropped_plus, n, n)
        wfs.state.cropped_minus = similar(wfs.state.cropped_minus, n, n)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.focal_field)
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

function host_curvature_defocus_masks(wfs::CurvatureWFS, tel::Telescope)
    T = eltype(wfs.state.frame_plus)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    cx = (pad + 1) / 2
    radius = T(n) / 2
    # SPECULA uses a defocus RMS in nanometers; keep the current public surface
    # stable by mapping the legacy response gain to a similar phase amplitude.
    defocus_rms_nm = T(1000) * wfs.params.response_gain
    plus = Matrix{Complex{T}}(undef, pad, pad)
    minus = Matrix{Complex{T}}(undef, pad, pad)
    inv_waves = defocus_rms_nm / T(1e9)
    @inbounds for j in 1:pad, i in 1:pad
        x = (T(i) - T(cx)) / radius
        y = (T(j) - T(cx)) / radius
        r2 = x * x + y * y
        z4 = r2 <= one(T) ? sqrt(T(3)) * (T(2) * r2 - one(T)) : zero(T)
        phase_cycles = inv_waves * z4
        plus[i, j] = cispi(T(2) * phase_cycles)
        minus[i, j] = cispi(-T(2) * phase_cycles)
    end
    return plus, minus
end

function build_curvature_defocus_masks!(wfs::CurvatureWFS, tel::Telescope)
    plus, minus = host_curvature_defocus_masks(wfs, tel)
    copyto!(wfs.state.defocus_plus, plus)
    copyto!(wfs.state.defocus_minus, minus)
    return wfs
end

function sample_curvature_frames!(wfs::CurvatureWFS, tel::Telescope)
    n = tel.params.resolution
    pad = size(wfs.state.intensity_plus, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    sub = div(n, wfs.params.n_subap)
    copyto!(wfs.state.cropped_plus, @view(wfs.state.intensity_plus[ox+1:ox+n, oy+1:oy+n]))
    copyto!(wfs.state.cropped_minus, @view(wfs.state.intensity_minus[ox+1:ox+n, oy+1:oy+n]))
    bin2d!(wfs.state.frame_plus, wfs.state.cropped_plus, sub)
    bin2d!(wfs.state.frame_minus, wfs.state.cropped_minus, sub)
    @views copyto!(wfs.state.camera_frame[1:wfs.params.n_subap, :], wfs.state.frame_plus)
    @views copyto!(wfs.state.camera_frame[wfs.params.n_subap+1:2*wfs.params.n_subap, :], wfs.state.frame_minus)
    return wfs.state.camera_frame
end

function curvature_intensity!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.frame_plus)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))

    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
        cispi(opd_to_cycles * tel.state.opd)

    copyto!(wfs.state.focal_field, wfs.state.field)
    @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.defocus_plus * wfs.state.phasor
    execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
    @. wfs.state.intensity_plus = abs2(wfs.state.focal_field)

    copyto!(wfs.state.focal_field, wfs.state.field)
    @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.defocus_minus * wfs.state.phasor
    execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
    @. wfs.state.intensity_minus = abs2(wfs.state.focal_field)

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
