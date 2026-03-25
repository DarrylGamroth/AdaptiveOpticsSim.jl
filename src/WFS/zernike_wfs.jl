#
# Zernike wavefront sensing
#
# The maintained MVP follows the standard phase-shifting Zernike-wavefront-
# sensor optical chain:
#
# 1. embed the pupil field on a padded FFT grid
# 2. propagate to the focal plane
# 3. apply a circular phase-shifting spot of radius `spot_radius_lambda_over_d`
# 4. propagate back to the pupil plane
# 5. sample the re-imaged pupil intensity on a compact detector grid
# 6. normalize and subtract a stored reference signal
#
# The exported 1-D `slopes` vector is therefore not a centroid slope field. It
# is the normalized valid-pupil intensity signal used by the current runtime and
# controller surfaces for compatibility with other WFS types.
#

struct ZernikeWFSParams{T<:AbstractFloat,N<:WFSNormalization}
    n_subap::Int
    threshold::T
    phase_shift_pi::T
    spot_radius_lambda_over_d::T
    normalization::N
    diffraction_padding::Int
    binning::Int
end

mutable struct ZernikeWFSState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    I<:AbstractVector{Int},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
}
    valid_mask::A
    valid_signal_indices::I
    slopes::V
    field::C
    focal_field::C
    pupil_field::C
    phasor::C
    phase_mask::C
    pupil_intensity::R
    nominal_frame::R
    camera_frame::R
    signal_2d::R
    reference_signal_2d::R
    reference_frame::R
    normalization_frame::R
    fft_plan::Pf
    ifft_plan::Pi
    effective_padding::Int
    calibrated::Bool
    calibration_wavelength::T
end

"""
    ZernikeWFS

Diffractive Zernike wavefront sensor with a circular focal-plane phase spot.

The sensor stores the sampled pupil-intensity signal in `state.slopes` for
runtime compatibility. For this sensor the vector represents a normalized pupil
intensity signal, not geometric or centroid slopes.
"""
struct ZernikeWFS{P<:ZernikeWFSParams,S<:ZernikeWFSState} <: AbstractWFS
    params::P
    state::S
end

"""
    ZernikeWFS(tel; ...)

Construct a Zernike WFS using a focal-plane circular phase-shifting spot.

`n_subap` defines the nominal sampled pupil grid before optional `binning`
coarsens the final exported camera/signal frame.
"""
function ZernikeWFS(tel::Telescope; n_subap::Int,
    phase_shift_pi::Real=0.5,
    spot_radius_lambda_over_d::Real=1.0,
    threshold::Real=0.0,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    diffraction_padding::Int=2,
    binning::Int=1,
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if n_subap % binning != 0
        throw(InvalidConfiguration("n_subap must be divisible by binning"))
    end
    if diffraction_padding < 1
        throw(InvalidConfiguration("diffraction_padding must be >= 1"))
    end
    n_signal = div(n_subap, binning)
    pad = tel.params.resolution * diffraction_padding
    params = ZernikeWFSParams{T,typeof(normalization)}(
        n_subap,
        T(threshold),
        T(phase_shift_pi),
        T(spot_radius_lambda_over_d),
        normalization,
        diffraction_padding,
        binning,
    )
    valid_mask = backend{Bool}(undef, n_signal, n_signal)
    fill!(valid_mask, false)
    valid_signal_indices = backend{Int}(undef, max(1, n_signal * n_signal))
    slopes = backend{T}(undef, max(1, n_signal * n_signal))
    fill!(slopes, zero(T))
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    phasor = similar(field)
    phase_mask = similar(field)
    pupil_intensity = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    nominal_frame = backend{T}(undef, n_subap, n_subap)
    camera_frame = backend{T}(undef, n_signal, n_signal)
    signal_2d = backend{T}(undef, n_signal, n_signal)
    reference_signal_2d = similar(signal_2d)
    reference_frame = similar(camera_frame)
    normalization_frame = similar(camera_frame)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    state = ZernikeWFSState{
        T,
        typeof(valid_mask),
        typeof(valid_signal_indices),
        typeof(slopes),
        typeof(field),
        typeof(camera_frame),
        typeof(fft_plan),
        typeof(ifft_plan),
    }(
        valid_mask,
        valid_signal_indices,
        slopes,
        field,
        focal_field,
        pupil_field,
        phasor,
        phase_mask,
        pupil_intensity,
        nominal_frame,
        camera_frame,
        signal_2d,
        reference_signal_2d,
        reference_frame,
        normalization_frame,
        fft_plan,
        ifft_plan,
        diffraction_padding,
        false,
        zero(T),
    )
    wfs = ZernikeWFS{typeof(params),typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    build_zernike_phasor!(wfs.state.phasor)
    build_zernike_phase_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ZernikeWFS) = Diffractive()

function zernike_signal_resolution(wfs::ZernikeWFS)
    return div(wfs.params.n_subap, wfs.params.binning)
end

function update_zernike_valid_indices!(wfs::ZernikeWFS)
    valid_host = execution_style(wfs.state.valid_mask) isa ScalarCPUStyle ?
        wfs.state.valid_mask : Array(wfs.state.valid_mask)
    n_valid = count(valid_host)
    if length(wfs.state.valid_signal_indices) != n_valid
        wfs.state.valid_signal_indices = similar(wfs.state.valid_signal_indices, n_valid)
    end
    if length(wfs.state.slopes) != n_valid
        wfs.state.slopes = similar(wfs.state.slopes, n_valid)
    end
    host_indices = Vector{Int}(undef, n_valid)
    idx = 1
    @inbounds for j in axes(valid_host, 2), i in axes(valid_host, 1)
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * size(valid_host, 1)
            idx += 1
        end
    end
    copyto!(wfs.state.valid_signal_indices, host_indices)
    fill!(wfs.state.slopes, zero(eltype(wfs.state.slopes)))
    return wfs
end

function update_valid_mask!(wfs::ZernikeWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    update_zernike_valid_indices!(wfs)
    return wfs
end

function ensure_zernike_buffers!(wfs::ZernikeWFS, tel::Telescope)
    n = tel.params.resolution
    pad = n * wfs.params.diffraction_padding
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.focal_field = similar(wfs.state.focal_field, pad, pad)
        wfs.state.pupil_field = similar(wfs.state.pupil_field, pad, pad)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.phase_mask = similar(wfs.state.phase_mask, pad, pad)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.focal_field)
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.pupil_field)
        wfs.state.effective_padding = wfs.params.diffraction_padding
        wfs.state.calibrated = false
        build_zernike_phasor!(wfs.state.phasor)
        build_zernike_phase_mask!(wfs, tel)
    end
    return wfs
end

function build_zernike_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for j in 1:n, i in 1:n
        phasor[i, j] = cis(scale * (i + j - 2))
    end
    return phasor
end

function host_zernike_phase_mask(wfs::ZernikeWFS, tel::Telescope)
    n = tel.params.resolution
    pad = size(wfs.state.phase_mask, 1)
    T = eltype(wfs.state.camera_frame)
    host = Matrix{Complex{T}}(undef, pad, pad)
    center = T(pad) / 2
    radius = wfs.params.spot_radius_lambda_over_d * (T(pad) / T(n))
    phase = T(pi) * wfs.params.phase_shift_pi
    shifted = cis(phase)
    @inbounds for j in 1:pad, i in 1:pad
        x = T(i) - center - T(0.5)
        y = T(j) - center - T(0.5)
        host[i, j] = hypot(x, y) <= radius ? shifted : one(shifted)
    end
    return host
end

function build_zernike_phase_mask!(wfs::ZernikeWFS, tel::Telescope)
    copyto!(wfs.state.phase_mask, host_zernike_phase_mask(wfs, tel))
    return wfs.state.phase_mask
end

function sample_zernike_frame!(out::AbstractMatrix{T}, nominal::AbstractMatrix{T}, wfs::ZernikeWFS,
    input::AbstractMatrix{T}, tel::Telescope) where {T<:AbstractFloat}
    sub = div(tel.params.resolution, wfs.params.n_subap)
    bin2d!(nominal, input, sub)
    if wfs.params.binning == 1
        copyto!(out, nominal)
    else
        bin2d!(out, nominal, wfs.params.binning)
    end
    return out
end

function zernike_pupil_intensity!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    ensure_zernike_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.camera_frame)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.camera_frame)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))
    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * tel.state.pupil *
        cispi(opd_to_cycles * tel.state.opd)
    copyto!(wfs.state.focal_field, wfs.state.field)
    @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
    execute_fft_plan!(wfs.state.focal_field, wfs.state.fft_plan)
    @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phase_mask
    copyto!(wfs.state.pupil_field, wfs.state.focal_field)
    execute_fft_plan!(wfs.state.pupil_field, wfs.state.ifft_plan)
    @views @. wfs.state.pupil_intensity = abs2(wfs.state.pupil_field[ox+1:ox+n, oy+1:oy+n])
    return wfs.state.pupil_intensity
end

function zernike_normalization(::MeanValidFluxNormalization, wfs::ZernikeWFS, tel::Telescope,
    src::AbstractSource, frame::AbstractMatrix)
    vals = frame[wfs.state.valid_mask]
    return isempty(vals) ? one(eltype(frame)) : max(mean(vals), eps(eltype(frame)))
end

function zernike_normalization(::IncidenceFluxNormalization, wfs::ZernikeWFS, tel::Telescope,
    src::AbstractSource, frame::AbstractMatrix)
    fmap = flux_map(tel, src)
    sample_zernike_frame!(wfs.state.normalization_frame, wfs.state.nominal_frame, wfs, fmap, tel)
    vals = wfs.state.normalization_frame[wfs.state.valid_mask]
    return isempty(vals) ? one(eltype(frame)) : max(mean(vals), eps(eltype(frame)))
end

function zernike_signal!(wfs::ZernikeWFS, tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    if size(frame) != size(wfs.state.signal_2d)
        throw(DimensionMismatchError("ZernikeWFS frame size must match sampled camera frame"))
    end
    norm_factor = zernike_normalization(wfs.params.normalization, wfs, tel, src, frame)
    fill!(wfs.state.signal_2d, zero(T))
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        if wfs.state.valid_mask[i, j]
            wfs.state.signal_2d[i, j] = frame[i, j] / norm_factor - wfs.state.reference_signal_2d[i, j]
        end
    end
    @inbounds for idx in eachindex(wfs.state.valid_signal_indices)
        wfs.state.slopes[idx] = wfs.state.signal_2d[wfs.state.valid_signal_indices[idx]]
    end
    return wfs.state.slopes
end

function ensure_zernike_calibration!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    λ = eltype(wfs.state.slopes)(wavelength(src))
    if wfs.state.calibrated && wfs.state.calibration_wavelength == λ
        return wfs
    end
    opd_saved = copy(tel.state.opd)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    zernike_pupil_intensity!(wfs, tel, src)
    sample_zernike_frame!(wfs.state.camera_frame, wfs.state.nominal_frame, wfs, wfs.state.pupil_intensity, tel)
    fill!(wfs.state.reference_signal_2d, zero(eltype(wfs.state.reference_signal_2d)))
    zernike_signal!(wfs, tel, wfs.state.camera_frame, src)
    copyto!(wfs.state.reference_signal_2d, wfs.state.signal_2d)
    copyto!(wfs.state.reference_frame, wfs.state.camera_frame)
    fill!(wfs.state.signal_2d, zero(eltype(wfs.state.signal_2d)))
    fill!(wfs.state.slopes, zero(eltype(wfs.state.slopes)))
    copyto!(tel.state.opd, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    return wfs
end

function measure!(::Diffractive, wfs::ZernikeWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive ZernikeWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::ZernikeWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::ZernikeWFS, tel::Telescope, ast::Asterism)
    throw(InvalidConfiguration("ZernikeWFS asterism support is not implemented in the Phase 1 MVP"))
end

function measure!(wfs::ZernikeWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration("ZernikeWFS asterism support is not implemented in the Phase 1 MVP"))
end

function measure!(::Diffractive, wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    ensure_zernike_calibration!(wfs, tel, src)
    zernike_pupil_intensity!(wfs, tel, src)
    sample_zernike_frame!(wfs.state.camera_frame, wfs.state.nominal_frame, wfs, wfs.state.pupil_intensity, tel)
    return zernike_signal!(wfs, tel, wfs.state.camera_frame, src)
end

function measure!(::Diffractive, wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_zernike_calibration!(wfs, tel, src)
    zernike_pupil_intensity!(wfs, tel, src)
    sample_zernike_frame!(wfs.state.camera_frame, wfs.state.nominal_frame, wfs, wfs.state.pupil_intensity, tel)
    capture!(det, wfs.state.camera_frame; rng=rng)
    size(output_frame(det)) == size(wfs.state.camera_frame) ||
        throw(InvalidConfiguration("ZernikeWFS detector output size must match the sampled camera frame"))
    return zernike_signal!(wfs, tel, output_frame(det), src)
end
