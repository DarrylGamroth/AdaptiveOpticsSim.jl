#
# Curvature wavefront sensing
#
# This file implements a compact transport-of-intensity-style curvature sensor.
# The pupil is partitioned into subapertures, the mean OPD in each subaperture
# is converted to a discrete curvature signal through a 5-point Laplacian on
# the subaperture grid, and the resulting signal is represented as a pair of
# intra-/extra-focal pupil images.
#
# This is intentionally a simpler optical model than the diffractive SH/Pyramid
# paths, but it preserves the important AO-facing structure:
#
# 1. a scalar curvature measurement per valid subaperture,
# 2. a two-plane detector-like readout surface,
# 3. detector-coupled measurement when an explicit detector is attached, and
# 4. interaction-matrix / reconstructor compatibility through `state.slopes`.
#

@kernel function curvature_subaperture_kernel!(subap_phase, subap_flux, opd, pupil, sub::Int, n_sub::Int, inv_area)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        acc = zero(eltype(subap_phase))
        cnt = zero(eltype(subap_phase))
        @inbounds for di in 0:(sub - 1), dj in 0:(sub - 1)
            x = xs + di
            y = ys + dj
            if pupil[x, y]
                acc += opd[x, y]
                cnt += one(eltype(subap_phase))
            end
        end
        @inbounds begin
            subap_phase[i, j] = cnt > 0 ? acc / cnt : zero(eltype(subap_phase))
            subap_flux[i, j] = cnt * inv_area
        end
    end
end

@kernel function curvature_frame_kernel!(signal_2d, frame_plus, frame_minus, camera_frame, valid_mask,
    subap_phase, subap_flux, inv_wavelength, response_gain, epsval, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        valid = @inbounds valid_mask[i, j]
        center = @inbounds subap_phase[i, j]
        left = j > 1 && @inbounds(valid_mask[i, j - 1]) ? @inbounds(subap_phase[i, j - 1]) : center
        right = j < n_sub && @inbounds(valid_mask[i, j + 1]) ? @inbounds(subap_phase[i, j + 1]) : center
        down = i > 1 && @inbounds(valid_mask[i - 1, j]) ? @inbounds(subap_phase[i - 1, j]) : center
        up = i < n_sub && @inbounds(valid_mask[i + 1, j]) ? @inbounds(subap_phase[i + 1, j]) : center
        lap = left + right + up + down - 4 * center
        normalized = response_gain * lap * inv_wavelength
        flux = @inbounds subap_flux[i, j]
        plus = valid ? max(zero(eltype(signal_2d)), flux * (one(eltype(signal_2d)) + normalized)) : zero(eltype(signal_2d))
        minus = valid ? max(zero(eltype(signal_2d)), flux * (one(eltype(signal_2d)) - normalized)) : zero(eltype(signal_2d))
        denom = plus + minus + epsval
        signal = valid ? (plus - minus) / denom : zero(eltype(signal_2d))
        @inbounds begin
            signal_2d[i, j] = signal
            frame_plus[i, j] = plus
            frame_minus[i, j] = minus
            camera_frame[i, j] = plus
            camera_frame[i + n_sub, j] = minus
        end
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

Curvature wavefront sensor using a two-plane pupil-intensity representation.

The measurement is a scalar curvature signal per valid subaperture. The sensor
exports a detector-like frame with the intra-focal image stacked above the
extra-focal image, so it can be used in detector-coupled runtime paths.
"""
struct CurvatureWFSParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    response_gain::T
end

mutable struct CurvatureWFSState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    M<:AbstractMatrix{T},
}
    valid_mask::A
    slopes::V
    subap_phase::M
    subap_flux::M
    signal_2d::M
    reference_signal_2d::M
    frame_plus::M
    frame_minus::M
    camera_frame::M
    calibrated::Bool
    calibration_wavelength::T
end

struct CurvatureWFS{P<:CurvatureWFSParams,S<:CurvatureWFSState} <: AbstractWFS
    params::P
    state::S
end

"""
    CurvatureWFS(tel; n_subap, threshold=0.1, response_gain=0.5, ...)

Construct a curvature WFS on the telescope pupil grid.

The sensor partitions the pupil into `n_subap × n_subap` cells, computes a
transport-of-intensity-style curvature signal on that grid, and exports a
`(2n_subap, n_subap)` camera frame containing the intra-/extra-focal images.
"""
function CurvatureWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1, response_gain::Real=0.5,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    tel.params.resolution % n_subap == 0 ||
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    params = CurvatureWFSParams{T}(n_subap, T(threshold), T(response_gain))
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, n_subap * n_subap)
    subap_phase = backend{T}(undef, n_subap, n_subap)
    subap_flux = backend{T}(undef, n_subap, n_subap)
    signal_2d = backend{T}(undef, n_subap, n_subap)
    reference_signal_2d = backend{T}(undef, n_subap, n_subap)
    frame_plus = backend{T}(undef, n_subap, n_subap)
    frame_minus = backend{T}(undef, n_subap, n_subap)
    camera_frame = backend{T}(undef, 2 * n_subap, n_subap)
    fill!(slopes, zero(T))
    fill!(subap_phase, zero(T))
    fill!(subap_flux, zero(T))
    fill!(signal_2d, zero(T))
    fill!(reference_signal_2d, zero(T))
    fill!(frame_plus, zero(T))
    fill!(frame_minus, zero(T))
    fill!(camera_frame, zero(T))
    state = CurvatureWFSState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(subap_phase),
    }(
        valid_mask,
        slopes,
        subap_phase,
        subap_flux,
        signal_2d,
        reference_signal_2d,
        frame_plus,
        frame_minus,
        camera_frame,
        false,
        zero(T),
    )
    wfs = CurvatureWFS{typeof(params),typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

sensing_mode(::CurvatureWFS) = Diffractive()

function update_valid_mask!(wfs::CurvatureWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function curvature_subapertures!(wfs::CurvatureWFS, tel::Telescope)
    n_sub = wfs.params.n_subap
    sub = div(tel.params.resolution, n_sub)
    inv_area = inv(eltype(wfs.state.subap_phase)(sub * sub))
    style = execution_style(wfs.state.subap_phase)
    if style isa ScalarCPUStyle
        @inbounds for i in 1:n_sub, j in 1:n_sub
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            acc = zero(eltype(wfs.state.subap_phase))
            cnt = zero(eltype(wfs.state.subap_phase))
            for di in 0:(sub - 1), dj in 0:(sub - 1)
                x = xs + di
                y = ys + dj
                if tel.state.pupil[x, y]
                    acc += tel.state.opd[x, y]
                    cnt += one(eltype(wfs.state.subap_phase))
                end
            end
            wfs.state.subap_phase[i, j] = cnt > 0 ? acc / cnt : zero(eltype(wfs.state.subap_phase))
            wfs.state.subap_flux[i, j] = cnt * inv_area
        end
    else
        launch_kernel!(style, curvature_subaperture_kernel!, wfs.state.subap_phase, wfs.state.subap_flux,
            tel.state.opd, tel.state.pupil, sub, n_sub, inv_area; ndrange=size(wfs.state.subap_phase))
    end
    return wfs
end

function curvature_frame!(wfs::CurvatureWFS, src::AbstractSource)
    n_sub = wfs.params.n_subap
    epsval = eps(eltype(wfs.state.signal_2d))
    style = execution_style(wfs.state.signal_2d)
    if style isa ScalarCPUStyle
        @inbounds for i in 1:n_sub, j in 1:n_sub
            if wfs.state.valid_mask[i, j]
                center = wfs.state.subap_phase[i, j]
                left = j > 1 && wfs.state.valid_mask[i, j - 1] ? wfs.state.subap_phase[i, j - 1] : center
                right = j < n_sub && wfs.state.valid_mask[i, j + 1] ? wfs.state.subap_phase[i, j + 1] : center
                down = i > 1 && wfs.state.valid_mask[i - 1, j] ? wfs.state.subap_phase[i - 1, j] : center
                up = i < n_sub && wfs.state.valid_mask[i + 1, j] ? wfs.state.subap_phase[i + 1, j] : center
                lap = left + right + up + down - 4 * center
                normalized = wfs.params.response_gain * lap / src.params.wavelength
                flux = wfs.state.subap_flux[i, j]
                plus = max(zero(eltype(wfs.state.signal_2d)), flux * (one(eltype(wfs.state.signal_2d)) + normalized))
                minus = max(zero(eltype(wfs.state.signal_2d)), flux * (one(eltype(wfs.state.signal_2d)) - normalized))
                signal = (plus - minus) / (plus + minus + epsval)
                wfs.state.signal_2d[i, j] = signal
                wfs.state.frame_plus[i, j] = plus
                wfs.state.frame_minus[i, j] = minus
                wfs.state.camera_frame[i, j] = plus
                wfs.state.camera_frame[i + n_sub, j] = minus
            else
                wfs.state.signal_2d[i, j] = zero(eltype(wfs.state.signal_2d))
                wfs.state.frame_plus[i, j] = zero(eltype(wfs.state.signal_2d))
                wfs.state.frame_minus[i, j] = zero(eltype(wfs.state.signal_2d))
                wfs.state.camera_frame[i, j] = zero(eltype(wfs.state.signal_2d))
                wfs.state.camera_frame[i + n_sub, j] = zero(eltype(wfs.state.signal_2d))
            end
        end
    else
        launch_kernel!(style, curvature_frame_kernel!, wfs.state.signal_2d, wfs.state.frame_plus,
            wfs.state.frame_minus, wfs.state.camera_frame, wfs.state.valid_mask, wfs.state.subap_phase,
            wfs.state.subap_flux, inv(src.params.wavelength), wfs.params.response_gain, epsval, n_sub;
            ndrange=size(wfs.state.signal_2d))
    end
    return wfs
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
    curvature_subapertures!(wfs, tel)
    curvature_frame!(wfs, src)
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
    curvature_subapertures!(wfs, tel)
    curvature_frame!(wfs, src)
    return curvature_signal!(wfs, wfs.state.camera_frame)
end

function measure!(::Diffractive, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_subapertures!(wfs, tel)
    curvature_frame!(wfs, src)
    capture!(det, wfs.state.camera_frame; rng=rng)
    size(output_frame(det)) == size(wfs.state.camera_frame) ||
        throw(InvalidConfiguration("CurvatureWFS detector output size must match the sampled camera frame"))
    return curvature_signal!(wfs, output_frame(det))
end
