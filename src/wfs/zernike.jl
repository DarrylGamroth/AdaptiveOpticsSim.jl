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
# controller's family-neutral WFS signal contract.
#

@kernel function zernike_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end

@kernel function zernike_masked_row_sum_kernel!(partials, frame, valid_mask,
    n1::Int, n2::Int)
    i = @index(Global, Linear)
    if i <= n1
        summed = zero(eltype(partials))
        @inbounds for j in 1:n2
            if valid_mask[i, j]
                summed += eltype(partials)(frame[i, j])
            end
        end
        @inbounds partials[i] = summed
    end
end

@kernel function zernike_finalize_normalization_sum_kernel!(
    normalization_sum, partials, count::Int)
    idx = @index(Global, Linear)
    if idx == 1
        summed = zero(eltype(normalization_sum))
        @inbounds for i in 1:count
            summed += partials[i]
        end
        @inbounds normalization_sum[1] = summed
    end
end

@kernel function zernike_signal_kernel!(signal_2d, frame, valid_mask,
    reference_signal_2d, normalization_sum, normalization_multiplier,
    clamp_to_epsilon::Bool, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n1 && j <= n2
        T = eltype(signal_2d)
        normalization = @inbounds(normalization_sum[1]) *
            T(normalization_multiplier)
        if clamp_to_epsilon
            normalization = max(normalization, eps(T))
        end
        usable = isfinite(normalization) && normalization > eps(T)
        @inbounds signal_2d[i, j] = valid_mask[i, j] && usable ?
            T(frame[i, j]) / normalization - reference_signal_2d[i, j] : zero(T)
    end
end

@kernel function gather_zernike_signal_kernel!(slopes, signal_2d, valid_signal_indices, count::Int)
    idx = @index(Global, Linear)
    if idx <= count
        src = @inbounds valid_signal_indices[idx]
        @inbounds slopes[idx] = signal_2d[src]
    end
end

struct ZernikeWFSParams{T<:AbstractFloat,N<:WFSNormalization}
    pupil_samples::Int
    pupil_resolution::Int
    pupil_diameter_m::T
    threshold::T
    phase_shift_pi::T
    spot_radius_lambda_over_d::T
    normalization::N
    diffraction_padding::Int
    binning::Int
end

"""Immutable focal-plane phase-shifting spot for a Zernike WFS."""
struct ZernikePhaseSpot{T<:AbstractFloat}
    phase_shift_pi::T
    radius_lambda_over_d::T
    diffraction_padding::Int
end

"""Single-writer FFT and pupil-relay workspace."""
mutable struct PreparedZernikePropagation{
    T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
}
    field::C
    focal_field::C
    pupil_field::C
    phasor::C
    phase_mask::C
    pupil_intensity::R
    nominal_frame::R
    fft_plan::Pf
    ifft_plan::Pi
    effective_padding::Int
    revision::UInt
end

"""Prepared Zernike phase spot and re-imaged-pupil optical front end."""
struct ZernikeOpticalFrontEnd{T<:AbstractFloat,M,P,S}
    phase_spot::M
    propagation::P
    pupil_resolution::Int
    pupil_diameter_m::T
    pupil_samples::Int
    binning::Int
    source::S
end

"""Mutable detector-facing frame storage owned by acquisition."""
mutable struct ZernikeAcquisitionState{R<:AbstractMatrix}
    camera_frame::R
end

struct ZernikeDetectorAcquisition{S}
    state::S
end

"""Mutable support, reference, and output storage owned by estimation."""
mutable struct ZernikeEstimatorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    I<:AbstractVector{Int},
    V<:AbstractVector{T},
    R<:AbstractMatrix{T},
    H<:AbstractVector{T},
}
    valid_mask::A
    valid_signal_indices::I
    slopes::V
    signal_2d::R
    reference_signal_2d::R
    reference_frame::R
    normalization_frame::R
    normalization_partials::V
    normalization_sum::V
    normalization_sum_host::H
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

struct ZernikePupilEstimator{S}
    state::S
end

"""
    ZernikeWFS

Diffractive Zernike wavefront sensor with a circular focal-plane phase spot.

The estimator stores a normalized pupil-intensity signal, not geometric or
centroid slopes. Access it through `slopes(sensor)` or
`sensor.estimator.state.slopes`.
"""
struct ZernikeWFS{P<:ZernikeWFSParams,F,A,E,B<:AbstractArrayBackend} <:
    AbstractWFS
    params::P
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::ZernikeWFS{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

"""
    ZernikeWFS(tel; ...)

Construct a Zernike WFS using a focal-plane circular phase-shifting spot.

`pupil_samples` defines the nominal sampled pupil grid before optional `binning`
coarsens the final exported camera/signal frame.
"""
function ZernikeWFS(tel::Telescope; pupil_samples::Int,
    phase_shift_pi::Real=0.5,
    spot_radius_lambda_over_d::Real=1.0,
    threshold::Real=0.0,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    diffraction_padding::Int=2,
    binning::Int=1,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    if tel.params.resolution % pupil_samples != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by pupil_samples"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if pupil_samples % binning != 0
        throw(InvalidConfiguration("pupil_samples must be divisible by binning"))
    end
    if diffraction_padding < 1
        throw(InvalidConfiguration("diffraction_padding must be >= 1"))
    end
    n_signal = div(pupil_samples, binning)
    pad = tel.params.resolution * diffraction_padding
    params = ZernikeWFSParams{T,typeof(normalization)}(
        pupil_samples,
        tel.params.resolution,
        T(tel.params.diameter),
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
    nominal_frame = backend{T}(undef, pupil_samples, pupil_samples)
    camera_frame = backend{T}(undef, n_signal, n_signal)
    signal_2d = backend{T}(undef, n_signal, n_signal)
    reference_signal_2d = similar(signal_2d)
    reference_frame = similar(camera_frame)
    normalization_frame = similar(camera_frame)
    normalization_partials = backend{T}(undef, size(camera_frame, 1))
    fill!(normalization_partials, zero(T))
    normalization_sum = backend{T}(undef, 1)
    fill!(normalization_sum, zero(T))
    normalization_sum_host = zeros(T, 1)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    propagation = PreparedZernikePropagation{
        T,
        typeof(field),
        typeof(pupil_intensity),
        typeof(fft_plan),
        typeof(ifft_plan),
    }(
        field,
        focal_field,
        pupil_field,
        phasor,
        phase_mask,
        pupil_intensity,
        nominal_frame,
        fft_plan,
        ifft_plan,
        diffraction_padding,
        UInt(0),
    )
    acquisition = ZernikeDetectorAcquisition(
        ZernikeAcquisitionState(camera_frame))
    estimator_state = ZernikeEstimatorState{
        T,
        typeof(valid_mask),
        typeof(valid_signal_indices),
        typeof(slopes),
        typeof(signal_2d),
        typeof(normalization_sum_host),
    }(
        valid_mask,
        valid_signal_indices,
        slopes,
        signal_2d,
        reference_signal_2d,
        reference_frame,
        normalization_frame,
        normalization_partials,
        normalization_sum,
        normalization_sum_host,
        false,
        zero(T),
        UInt(0),
        UInt(0),
    )
    estimator = ZernikePupilEstimator(estimator_state)
    spot = ZernikePhaseSpot(params.phase_shift_pi,
        params.spot_radius_lambda_over_d, params.diffraction_padding)
    front_end = ZernikeOpticalFrontEnd(spot, propagation,
        params.pupil_resolution, params.pupil_diameter_m,
        params.pupil_samples, params.binning, nothing)
    wfs = ZernikeWFS{
        typeof(params),typeof(front_end),typeof(acquisition),
        typeof(estimator),typeof(selector),
    }(params, front_end, acquisition, estimator)
    update_valid_mask!(wfs, tel)
    build_zernike_phasor!(wfs.front_end.propagation.phasor)
    build_zernike_phase_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ZernikeWFS) = Diffractive()

function zernike_signal_resolution(wfs::ZernikeWFS)
    return div(wfs.params.pupil_samples, wfs.params.binning)
end

function update_zernike_valid_indices!(wfs::ZernikeWFS)
    valid_host = host_array(wfs.estimator.state.valid_mask)
    n_valid = count(valid_host)
    if length(wfs.estimator.state.valid_signal_indices) != n_valid
        wfs.estimator.state.valid_signal_indices = similar(wfs.estimator.state.valid_signal_indices, n_valid)
    end
    if length(wfs.estimator.state.slopes) != n_valid
        wfs.estimator.state.slopes = similar(wfs.estimator.state.slopes, n_valid)
    end
    host_indices = Vector{Int}(undef, n_valid)
    idx = 1
    @inbounds for j in axes(valid_host, 2), i in axes(valid_host, 1)
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * size(valid_host, 1)
            idx += 1
        end
    end
    copyto!(wfs.estimator.state.valid_signal_indices, host_indices)
    fill!(wfs.estimator.state.slopes, zero(eltype(wfs.estimator.state.slopes)))
    return wfs
end

function update_valid_mask!(wfs::ZernikeWFS, tel::Telescope)
    set_valid_subapertures!(wfs.estimator.state.valid_mask, pupil_mask(tel), wfs.params.threshold)
    sample_zernike_frame!(wfs.estimator.state.normalization_frame,
        wfs.front_end.propagation.nominal_frame, wfs, pupil_reflectivity(tel), tel)
    update_zernike_valid_indices!(wfs)
    return wfs
end

function ensure_zernike_buffers!(wfs::ZernikeWFS, tel::Telescope)
    n = tel.params.resolution
    pad = n * wfs.params.diffraction_padding
    if size(wfs.front_end.propagation.field) != (pad, pad)
        wfs.front_end.propagation.field = similar(wfs.front_end.propagation.field, pad, pad)
        wfs.front_end.propagation.focal_field = similar(wfs.front_end.propagation.focal_field, pad, pad)
        wfs.front_end.propagation.pupil_field = similar(wfs.front_end.propagation.pupil_field, pad, pad)
        wfs.front_end.propagation.phasor = similar(wfs.front_end.propagation.phasor, pad, pad)
        wfs.front_end.propagation.phase_mask = similar(wfs.front_end.propagation.phase_mask, pad, pad)
        wfs.front_end.propagation.fft_plan = plan_fft_backend!(wfs.front_end.propagation.focal_field)
        wfs.front_end.propagation.ifft_plan = plan_ifft_backend!(wfs.front_end.propagation.pupil_field)
        wfs.front_end.propagation.effective_padding = wfs.params.diffraction_padding
        wfs.front_end.propagation.revision += UInt(1)
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
        build_zernike_phasor!(wfs.front_end.propagation.phasor)
        build_zernike_phase_mask!(wfs, tel)
    end
    return wfs
end

function build_zernike_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    return build_zernike_phasor!(execution_style(phasor), phasor)
end

function build_zernike_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for j in 1:n, i in 1:n
        phasor[i, j] = cis(scale * (i + j - 2))
    end
    return phasor
end

function build_zernike_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, zernike_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function host_zernike_phase_mask(wfs::ZernikeWFS, tel::Telescope)
    n = tel.params.resolution
    pad = size(wfs.front_end.propagation.phase_mask, 1)
    T = eltype(wfs.acquisition.state.camera_frame)
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
    copyto!(wfs.front_end.propagation.phase_mask, host_zernike_phase_mask(wfs, tel))
    return wfs.front_end.propagation.phase_mask
end

function sample_zernike_frame!(out::AbstractMatrix{T}, nominal::AbstractMatrix{T}, wfs::ZernikeWFS,
    input::AbstractMatrix{T}, tel::Telescope) where {T<:AbstractFloat}
    sub = div(tel.params.resolution, wfs.params.pupil_samples)
    bin2d!(nominal, input, sub)
    if wfs.params.binning == 1
        copyto!(out, nominal)
    else
        bin2d!(out, nominal, wfs.params.binning)
    end
    return out
end

function zernike_pupil_intensity!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS")
    ensure_zernike_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.front_end.propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.acquisition.state.camera_frame)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.acquisition.state.camera_frame)(
        photon_irradiance(src) * (tel.params.diameter / tel.params.resolution)^2
    ))
    fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
    @views @. wfs.front_end.propagation.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * sqrt($(pupil_reflectivity(tel))) *
        cispi(opd_to_cycles * tel.state.opd)
    copyto!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.field)
    @. wfs.front_end.propagation.focal_field = wfs.front_end.propagation.focal_field * wfs.front_end.propagation.phasor
    execute_fft_plan!(wfs.front_end.propagation.focal_field, wfs.front_end.propagation.fft_plan)
    @. wfs.front_end.propagation.focal_field = wfs.front_end.propagation.focal_field * wfs.front_end.propagation.phase_mask
    copyto!(wfs.front_end.propagation.pupil_field, wfs.front_end.propagation.focal_field)
    execute_fft_plan!(wfs.front_end.propagation.pupil_field, wfs.front_end.propagation.ifft_plan)
    @views @. wfs.front_end.propagation.pupil_intensity = abs2(wfs.front_end.propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    return wfs.front_end.propagation.pupil_intensity
end

function zernike_normalization(normalization::MeanValidFluxNormalization,
    wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    frame::AbstractMatrix)
    return zernike_normalization(normalization, wfs, tel, src, frame,
        one(eltype(frame)))
end

function zernike_normalization(normalization::WFSNormalization,
    wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    frame::AbstractMatrix{T}, normalization_scale::Real) where {T<:AbstractFloat}
    return zernike_normalization(execution_style(frame), normalization,
        wfs, tel, src, frame, T(normalization_scale))
end

function zernike_normalization(normalization::IncidenceFluxNormalization,
    wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    frame::AbstractMatrix)
    return zernike_normalization(normalization, wfs, tel, src, frame,
        one(eltype(frame)))
end

@inline zernike_normalization_count(wfs::ZernikeWFS) =
    length(wfs.estimator.state.valid_signal_indices)

@inline function zernike_incidence_sample_scale(tel::Telescope,
    src::AbstractSource, normalization_scale::T) where {T<:AbstractFloat}
    irradiance = T(_require_physical_photon_irradiance(src,
        "ZernikeWFS incidence normalization"))
    pupil_sample = T(tel.params.diameter) / T(tel.params.resolution)
    return irradiance * abs2(pupil_sample) * normalization_scale
end

@inline function finalize_zernike_normalization(
    ::MeanValidFluxNormalization, summed::T, multiplier::T) where {T<:AbstractFloat}
    return max(summed * multiplier, eps(T))
end

@inline function finalize_zernike_normalization(
    ::IncidenceFluxNormalization, summed::T, multiplier::T) where {T<:AbstractFloat}
    return summed * multiplier
end

function zernike_normalization(::ScalarCPUStyle,
    normalization::MeanValidFluxNormalization,
    wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    frame::AbstractMatrix{T}, ::T) where {T<:AbstractFloat}
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    summed = masked_sum2d(ScalarCPUStyle(), frame, wfs.estimator.state.valid_mask)
    return finalize_zernike_normalization(normalization, summed,
        inv(T(count)))
end

function zernike_normalization(::ScalarCPUStyle,
    normalization::IncidenceFluxNormalization,
    wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    frame::AbstractMatrix{T}, normalization_scale::T) where {T<:AbstractFloat}
    scale = zernike_incidence_sample_scale(tel, src, normalization_scale)
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    sample_binning = div(tel.params.resolution,
        size(wfs.estimator.state.normalization_frame, 1))
    bin2d!(wfs.estimator.state.normalization_frame, pupil_reflectivity(tel),
        sample_binning)
    summed = masked_sum2d(ScalarCPUStyle(),
        wfs.estimator.state.normalization_frame, wfs.estimator.state.valid_mask)
    return finalize_zernike_normalization(normalization, summed,
        scale / T(count))
end

@inline function queue_zernike_masked_sum!(phase::KernelLaunchPhase,
    wfs::ZernikeWFS, frame::AbstractMatrix)
    n1, n2 = size(wfs.estimator.state.valid_mask)
    queue_kernel!(phase, zernike_masked_row_sum_kernel!,
        wfs.estimator.state.normalization_partials, frame, wfs.estimator.state.valid_mask,
        n1, n2; ndrange=n1)
    queue_kernel!(phase, zernike_finalize_normalization_sum_kernel!,
        wfs.estimator.state.normalization_sum, wfs.estimator.state.normalization_partials,
        n1; ndrange=1)
    return nothing
end

function queue_zernike_normalization!(phase::KernelLaunchPhase,
    ::MeanValidFluxNormalization, wfs::ZernikeWFS, ::Telescope,
    ::AbstractSource, frame::AbstractMatrix{T}, ::T) where {T<:AbstractFloat}
    queue_zernike_masked_sum!(phase, wfs, frame)
    return inv(T(zernike_normalization_count(wfs))), true
end

function queue_zernike_normalization!(phase::KernelLaunchPhase,
    ::IncidenceFluxNormalization, wfs::ZernikeWFS, tel::Telescope,
    src::AbstractSource, ::AbstractMatrix{T},
    normalization_scale::T) where {T<:AbstractFloat}
    sample_binning = div(tel.params.resolution,
        size(wfs.estimator.state.normalization_frame, 1))
    n1, n2 = size(wfs.estimator.state.normalization_frame)
    queue_kernel!(phase, bin2d_kernel!, wfs.estimator.state.normalization_frame,
        pupil_reflectivity(tel), sample_binning, n1, n2;
        ndrange=(n1, n2))
    queue_zernike_masked_sum!(phase, wfs,
        wfs.estimator.state.normalization_frame)
    scale = zernike_incidence_sample_scale(tel, src,
        normalization_scale)
    return scale / T(zernike_normalization_count(wfs)), false
end

function zernike_normalization(style::AcceleratorStyle,
    normalization::WFSNormalization, wfs::ZernikeWFS, tel::Telescope,
    src::AbstractSource, frame::AbstractMatrix{T},
    normalization_scale::T) where {T<:AbstractFloat}
    # A host scalar result necessarily synchronizes. The measurement hot path
    # below queues this reduction with signal formation and never calls this
    # scalar inspection seam.
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    phase = begin_kernel_phase(style)
    multiplier, _ = queue_zernike_normalization!(phase, normalization,
        wfs, tel, src, frame, normalization_scale)
    finish_kernel_phase!(phase)
    copyto!(wfs.estimator.state.normalization_sum_host,
        wfs.estimator.state.normalization_sum)
    return finalize_zernike_normalization(normalization,
        wfs.estimator.state.normalization_sum_host[1], multiplier)
end

function zernike_signal!(wfs::ZernikeWFS, tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(wfs, tel, frame, src, one(T))
end

function zernike_signal!(wfs::ZernikeWFS, tel::Telescope,
    frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::Real) where {T<:AbstractFloat}
    return zernike_signal!(execution_style(frame), wfs, tel, frame, src,
        T(normalization_scale))
end

function zernike_signal!(::ScalarCPUStyle, wfs::ZernikeWFS, tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(ScalarCPUStyle(), wfs, tel, frame, src, one(T))
end

function zernike_signal!(::ScalarCPUStyle, wfs::ZernikeWFS,
    tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::T) where {T<:AbstractFloat}
    if size(frame) != size(wfs.estimator.state.signal_2d)
        throw(DimensionMismatchError("ZernikeWFS frame size must match sampled camera frame"))
    end
    norm_factor = zernike_normalization(wfs.params.normalization, wfs, tel,
        src, frame, normalization_scale)
    fill!(wfs.estimator.state.signal_2d, zero(T))
    if !usable_wfs_normalization(norm_factor)
        fill!(wfs.estimator.state.slopes, zero(T))
        return wfs.estimator.state.slopes
    end
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        if wfs.estimator.state.valid_mask[i, j]
            wfs.estimator.state.signal_2d[i, j] = frame[i, j] / norm_factor - wfs.estimator.state.reference_signal_2d[i, j]
        end
    end
    @inbounds for idx in eachindex(wfs.estimator.state.valid_signal_indices)
        wfs.estimator.state.slopes[idx] = wfs.estimator.state.signal_2d[wfs.estimator.state.valid_signal_indices[idx]]
    end
    return wfs.estimator.state.slopes
end

function zernike_signal!(style::AcceleratorStyle, wfs::ZernikeWFS, tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(style, wfs, tel, frame, src, one(T))
end

function zernike_signal!(style::AcceleratorStyle, wfs::ZernikeWFS,
    tel::Telescope, frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::T) where {T<:AbstractFloat}
    if size(frame) != size(wfs.estimator.state.signal_2d)
        throw(DimensionMismatchError("ZernikeWFS frame size must match sampled camera frame"))
    end
    count = zernike_normalization_count(wfs)
    if count == 0
        fill!(wfs.estimator.state.signal_2d, zero(T))
        fill!(wfs.estimator.state.slopes, zero(T))
        return wfs.estimator.state.slopes
    end
    phase = begin_kernel_phase(style)
    normalization_multiplier, clamp_to_epsilon =
        queue_zernike_normalization!(phase, wfs.params.normalization,
            wfs, tel, src, frame, normalization_scale)
    queue_kernel!(phase, zernike_signal_kernel!, wfs.estimator.state.signal_2d,
        frame, wfs.estimator.state.valid_mask, wfs.estimator.state.reference_signal_2d,
        wfs.estimator.state.normalization_sum, normalization_multiplier,
        clamp_to_epsilon, size(frame, 1), size(frame, 2);
        ndrange=size(frame))
    queue_kernel!(phase, gather_zernike_signal_kernel!, wfs.estimator.state.slopes,
        wfs.estimator.state.signal_2d, wfs.estimator.state.valid_signal_indices, count;
        ndrange=count)
    finish_kernel_phase!(phase)
    return wfs.estimator.state.slopes
end

function ensure_zernike_calibration!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    λ = calibration_wavelength(src, eltype(wfs.estimator.state.slopes))
    sig = telescope_aperture_calibration_signature(tel,
        calibration_signature(src))
    if calibration_matches(wfs.estimator.state.calibrated,
        wfs.estimator.state.calibration_wavelength, λ,
        wfs.estimator.state.calibration_signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, tel)
    opd_saved = save_zero_opd!(tel)
    try
        zernike_pupil_intensity!(wfs, tel, src)
        sample_zernike_frame!(wfs.acquisition.state.camera_frame,
            wfs.front_end.propagation.nominal_frame, wfs, wfs.front_end.propagation.pupil_intensity, tel)
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        zernike_signal!(wfs, tel, wfs.acquisition.state.camera_frame, src)
        copyto!(wfs.estimator.state.reference_signal_2d, wfs.estimator.state.signal_2d)
        copyto!(wfs.estimator.state.reference_frame, wfs.acquisition.state.camera_frame)
        fill!(wfs.estimator.state.signal_2d, zero(eltype(wfs.estimator.state.signal_2d)))
        fill!(wfs.estimator.state.slopes, zero(eltype(wfs.estimator.state.slopes)))
    finally
        restore_opd!(tel, opd_saved)
    end
    wfs.estimator.state.calibrated = true
    wfs.estimator.state.calibration_wavelength = λ
    wfs.estimator.state.calibration_signature = sig
    wfs.estimator.state.calibration_revision += UInt(1)
    return wfs
end

@inline function ensure_zernike_calibration!(wfs::ZernikeWFS,
    tel::Telescope, src::AbstractSource, ::AbstractDetector)
    return ensure_zernike_calibration!(wfs, tel, src)
end

function ensure_zernike_calibration!(wfs::ZernikeWFS, tel::Telescope,
    src::AbstractSource, det::Detector)
    T = eltype(wfs.estimator.state.slopes)
    λ = calibration_wavelength(src, T)
    sig = detector_calibration_signature(det,
        telescope_aperture_calibration_signature(tel,
            calibration_signature(src)))
    if calibration_matches(wfs.estimator.state.calibrated,
        wfs.estimator.state.calibration_wavelength, λ,
        wfs.estimator.state.calibration_signature, sig)
        return wfs
    end

    require_whole_capture_idle(det)
    update_valid_mask!(wfs, tel)
    opd_saved = save_zero_opd!(tel)
    try
        zernike_pupil_intensity!(wfs, tel, src)
        sample_zernike_frame!(wfs.acquisition.state.camera_frame,
            wfs.front_end.propagation.nominal_frame, wfs, wfs.front_end.propagation.pupil_intensity, tel)
        frame = detector_calibration_frame!(det, wfs.acquisition.state.camera_frame, src)
        size(frame) == size(wfs.estimator.state.signal_2d) || throw(
            InvalidConfiguration(
                "ZernikeWFS detector sampling and binning must preserve " *
                "the sampled camera frame size"))
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        normalization_scale = wfs_detector_incidence_scale(det, src,
            eltype(frame))
        zernike_signal!(wfs, tel, frame, src, normalization_scale)
        copyto!(wfs.estimator.state.reference_signal_2d, wfs.estimator.state.signal_2d)
        copyto!(wfs.estimator.state.reference_frame, frame)
        fill!(wfs.estimator.state.signal_2d, zero(eltype(wfs.estimator.state.signal_2d)))
        fill!(wfs.estimator.state.slopes, zero(eltype(wfs.estimator.state.slopes)))
    finally
        restore_opd!(tel, opd_saved)
    end
    wfs.estimator.state.calibrated = true
    wfs.estimator.state.calibration_wavelength = λ
    wfs.estimator.state.calibration_signature = sig
    wfs.estimator.state.calibration_revision += UInt(1)
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
    sample_zernike_frame!(wfs.acquisition.state.camera_frame, wfs.front_end.propagation.nominal_frame, wfs, wfs.front_end.propagation.pupil_intensity, tel)
    return zernike_signal!(wfs, tel, wfs.acquisition.state.camera_frame, src)
end

function measure!(::Diffractive, wfs::ZernikeWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_zernike_calibration!(wfs, tel, src, det)
    zernike_pupil_intensity!(wfs, tel, src)
    sample_zernike_frame!(wfs.acquisition.state.camera_frame, wfs.front_end.propagation.nominal_frame, wfs, wfs.front_end.propagation.pupil_intensity, tel)
    capture!(det, wfs.acquisition.state.camera_frame, src; rng=rng)
    size(output_frame(det)) == size(wfs.acquisition.state.camera_frame) ||
        throw(InvalidConfiguration("ZernikeWFS detector output size must match the sampled camera frame"))
    frame = output_frame(det)
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    return zernike_signal!(wfs, tel, frame, src, normalization_scale)
end

@inline slopes(wfs::ZernikeWFS) = wfs.estimator.state.slopes
@inline valid_subaperture_mask(wfs::ZernikeWFS) =
    wfs.estimator.state.valid_mask
@inline reference_signal(wfs::ZernikeWFS) =
    wfs.estimator.state.reference_signal_2d
@inline camera_frame(wfs::ZernikeWFS) =
    wfs.acquisition.state.camera_frame

@inline wfs_output_frame(wfs::ZernikeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame(wfs::ZernikeWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::ZernikeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::ZernikeWFS, det::AbstractDetector) = camera_frame(wfs)
@inline supports_prepared_runtime(::ZernikeWFS, src::AbstractSource) =
    is_leaf_source(src)
@inline supports_detector_output(::ZernikeWFS, ::AbstractDetector) = true

@inline function prepare_runtime_wfs!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS runtime preparation")
    ensure_zernike_calibration!(wfs, tel, src)
    return wfs
end

include("zernike/stages.jl")
