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

@kernel function gather_zernike_signal_kernel!(signal, signal_2d,
    valid_signal_indices, count::Int)
    idx = @index(Global, Linear)
    if idx <= count
        src = @inbounds valid_signal_indices[idx]
        @inbounds signal[idx] = signal_2d[src]
    end
end

"""Immutable configuration for normalized-pupil signal estimation."""
struct ZernikeEstimatorParams{T<:AbstractFloat,N<:WFSNormalization}
    pupil_samples::Int
    pupil_resolution::Int
    pupil_diameter_m::T
    threshold::T
    normalization::N
end

"""Run-immutable numerical contract for phase-spot pupil-relay propagation."""
struct ZernikePropagationPlan{M<:ZernikePhaseSpot,T<:AbstractFloat}
    phase_spot::M
    pupil_resolution::Int
    pupil_samples::Int
    numeric_type::Type{T}
end

"""
Backend-bound FFT handles and replaceable single-writer scratch for phase-spot
pupil-relay propagation. No field is a caller-visible optical product.
"""
mutable struct ZernikePropagationWorkspace{
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

"""Exact plan/workspace owner for one Zernike propagation execution."""
struct PreparedZernikePropagation{
    P<:ZernikePropagationPlan,W<:ZernikePropagationWorkspace}
    plan::P
    workspace::W
end

@inline zernike_propagation_plan(
    propagation::PreparedZernikePropagation) = propagation.plan
@inline zernike_propagation_workspace(
    propagation::PreparedZernikePropagation) = propagation.workspace

"""Zernike phase spot and prepared re-imaged-pupil optical front end."""
struct ZernikeOpticalFrontEnd{M,P,S}
    phase_spot::M
    propagation::P
    binning::Int
    source::S
end

"""Run-immutable binning contract for the internal detector-facing frame."""
struct ZernikeAcquisitionPlan
    binning::Int
end

"""Caller-visible detector-facing photon-arrival-rate frame product."""
mutable struct ZernikeAcquisitionProducts{T<:AbstractFloat,
    R<:AbstractMatrix{T}}
    frame::R
end

"""
Internal convenience owner for frame binning and its rate product. Detector
response, exposure integration, and readout belong to the generic acquisition
stage.
"""
struct ZernikeDetectorAcquisition{P,PR}
    plan::P
    products::PR
end

"""Persistent support and reference-calibration state."""
mutable struct ZernikeEstimatorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    R<:AbstractMatrix{T},
}
    valid_mask::A
    reference_signal_2d::R
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

"""Replaceable single-writer scratch for normalized-pupil estimation."""
mutable struct ZernikeEstimatorWorkspace{
    T<:AbstractFloat,
    I<:AbstractVector{Int},
    R<:AbstractMatrix{T},
    V<:AbstractVector{T},
    H<:AbstractVector{T},
}
    valid_signal_indices::I
    signal_2d::R
    normalization_frame::R
    normalization_partials::V
    normalization_sum::V
    normalization_sum_host::H
end

"""Caller-visible normalized-pupil signal product."""
mutable struct ZernikeEstimatorProducts{T<:AbstractFloat,
    V<:AbstractVector{T}}
    signal::V
end

struct ZernikePupilEstimator{P<:ZernikeEstimatorParams,S,W,PR}
    params::P
    state::S
    workspace::W
    products::PR
end

"""
    ZernikeWFS

Diffractive Zernike wavefront sensor with a circular focal-plane phase spot.

The estimator stores a normalized pupil-intensity signal, not geometric or
centroid slopes. Access it through `slopes(sensor)` or
`sensor.estimator.products.signal`.
"""
struct ZernikeWFS{F,A,E,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::ZernikeWFS{F,A,E,B}) where {F,A,E,B} = B()

@inline zernike_estimator_params(wfs::ZernikeWFS) = wfs.estimator.params
@inline zernike_estimator_state(wfs::ZernikeWFS) = wfs.estimator.state
@inline zernike_estimator_workspace(wfs::ZernikeWFS) =
    wfs.estimator.workspace
@inline zernike_estimator_products(wfs::ZernikeWFS) =
    wfs.estimator.products
@inline zernike_acquisition_plan(wfs::ZernikeWFS) = wfs.acquisition.plan
@inline zernike_acquisition_products(wfs::ZernikeWFS) =
    wfs.acquisition.products
@inline zernike_propagation(wfs::ZernikeWFS) = wfs.front_end.propagation
@inline zernike_propagation_plan(wfs::ZernikeWFS) =
    zernike_propagation_plan(zernike_propagation(wfs))
@inline zernike_propagation_workspace(wfs::ZernikeWFS) =
    zernike_propagation_workspace(zernike_propagation(wfs))
@inline zernike_propagation_workspace(front_end::ZernikeOpticalFrontEnd) =
    zernike_propagation_workspace(front_end.propagation)

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
    estimator_params = ZernikeEstimatorParams{T,typeof(normalization)}(
        pupil_samples,
        tel.params.resolution,
        T(tel.params.diameter),
        T(threshold),
        normalization,
    )
    spot = ZernikePhaseSpot(T(phase_shift_pi),
        T(spot_radius_lambda_over_d), diffraction_padding)
    valid_mask = backend{Bool}(undef, n_signal, n_signal)
    fill!(valid_mask, false)
    valid_signal_indices = backend{Int}(undef, max(1, n_signal * n_signal))
    signal = backend{T}(undef, max(1, n_signal * n_signal))
    fill!(signal, zero(T))
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
    normalization_frame = similar(camera_frame)
    normalization_partials = backend{T}(undef, size(camera_frame, 1))
    fill!(normalization_partials, zero(T))
    normalization_sum = backend{T}(undef, 1)
    fill!(normalization_sum, zero(T))
    normalization_sum_host = zeros(T, 1)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    propagation_workspace = ZernikePropagationWorkspace{
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
    propagation_plan = ZernikePropagationPlan(spot,
        tel.params.resolution, pupil_samples, T)
    propagation = PreparedZernikePropagation(propagation_plan,
        propagation_workspace)
    acquisition = ZernikeDetectorAcquisition(
        ZernikeAcquisitionPlan(binning),
        ZernikeAcquisitionProducts(camera_frame))
    estimator_state = ZernikeEstimatorState{
        T,
        typeof(valid_mask),
        typeof(signal_2d),
    }(
        valid_mask,
        reference_signal_2d,
        false,
        zero(T),
        UInt(0),
        UInt(0),
    )
    estimator_workspace = ZernikeEstimatorWorkspace{
        T,
        typeof(valid_signal_indices),
        typeof(signal_2d),
        typeof(normalization_partials),
        typeof(normalization_sum_host),
    }(
        valid_signal_indices,
        signal_2d,
        normalization_frame,
        normalization_partials,
        normalization_sum,
        normalization_sum_host,
    )
    estimator_products = ZernikeEstimatorProducts(signal)
    estimator = ZernikePupilEstimator(estimator_params, estimator_state,
        estimator_workspace, estimator_products)
    front_end = ZernikeOpticalFrontEnd(spot, propagation, binning, nothing)
    wfs = ZernikeWFS{
        typeof(front_end),typeof(acquisition),typeof(estimator),typeof(selector),
    }(front_end, acquisition, estimator)
    initial_pupil = PupilFunction(tel)
    update_valid_mask!(wfs, initial_pupil)
    build_zernike_phasor!(zernike_propagation_workspace(wfs).phasor)
    build_zernike_phase_mask!(wfs, initial_pupil)
    return wfs
end

sensing_mode(::ZernikeWFS) = Diffractive()

function zernike_signal_resolution(wfs::ZernikeWFS)
    plan = zernike_propagation_plan(wfs)
    return div(plan.pupil_samples, zernike_acquisition_plan(wfs).binning)
end

function update_zernike_valid_indices!(wfs::ZernikeWFS)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    products = zernike_estimator_products(wfs)
    valid_host = host_array(state.valid_mask)
    n_valid = count(valid_host)
    if length(workspace.valid_signal_indices) != n_valid
        workspace.valid_signal_indices = similar(
            workspace.valid_signal_indices, n_valid)
    end
    if length(products.signal) != n_valid
        products.signal = similar(products.signal, n_valid)
    end
    host_indices = Vector{Int}(undef, n_valid)
    idx = 1
    @inbounds for j in axes(valid_host, 2), i in axes(valid_host, 1)
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * size(valid_host, 1)
            idx += 1
        end
    end
    copyto!(workspace.valid_signal_indices, host_indices)
    fill!(products.signal, zero(eltype(products.signal)))
    return wfs
end

function update_valid_mask!(wfs::ZernikeWFS, pupil::PupilFunction)
    params = zernike_estimator_params(wfs)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    propagation = zernike_propagation_workspace(wfs)
    set_valid_subapertures!(state.valid_mask, pupil.support, params.threshold)
    sample_zernike_amplitude_squared!(
        workspace.normalization_frame,
        propagation.nominal_frame, wfs,
        pupil.amplitude, pupil)
    update_zernike_valid_indices!(wfs)
    return wfs
end

function ensure_zernike_buffers!(wfs::ZernikeWFS, pupil::PupilFunction)
    plan = zernike_propagation_plan(wfs)
    workspace = zernike_propagation_workspace(wfs)
    state = zernike_estimator_state(wfs)
    n = _pupil_resolution(pupil)
    diffraction_padding = plan.phase_spot.diffraction_padding
    pad = n * diffraction_padding
    if size(workspace.field) != (pad, pad)
        workspace.field = similar(workspace.field, pad, pad)
        workspace.focal_field = similar(workspace.focal_field, pad, pad)
        workspace.pupil_field = similar(workspace.pupil_field, pad, pad)
        workspace.phasor = similar(workspace.phasor, pad, pad)
        workspace.phase_mask = similar(workspace.phase_mask, pad, pad)
        workspace.fft_plan = plan_fft_backend!(workspace.focal_field)
        workspace.ifft_plan = plan_ifft_backend!(workspace.pupil_field)
        workspace.effective_padding = diffraction_padding
        workspace.revision += UInt(1)
        state.calibrated = false
        state.calibration_revision += UInt(1)
        build_zernike_phasor!(workspace.phasor)
        build_zernike_phase_mask!(wfs, pupil)
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

function host_zernike_phase_mask(wfs::ZernikeWFS, pupil::PupilFunction)
    phase_spot = wfs.front_end.phase_spot
    workspace = zernike_propagation_workspace(wfs)
    n = _pupil_resolution(pupil)
    pad = size(workspace.phase_mask, 1)
    T = eltype(zernike_acquisition_products(wfs).frame)
    host = Matrix{Complex{T}}(undef, pad, pad)
    center = T(pad) / 2
    radius = phase_spot.radius_lambda_over_d * (T(pad) / T(n))
    phase = T(pi) * phase_spot.phase_shift_pi
    shifted = cis(phase)
    @inbounds for j in 1:pad, i in 1:pad
        x = T(i) - center - T(0.5)
        y = T(j) - center - T(0.5)
        host[i, j] = hypot(x, y) <= radius ? shifted : one(shifted)
    end
    return host
end

function build_zernike_phase_mask!(wfs::ZernikeWFS, pupil::PupilFunction)
    phase_mask = zernike_propagation_workspace(wfs).phase_mask
    copyto!(phase_mask, host_zernike_phase_mask(wfs, pupil))
    return phase_mask
end

function sample_zernike_frame!(out::AbstractMatrix{T}, nominal::AbstractMatrix{T}, wfs::ZernikeWFS,
    input::AbstractMatrix{T}, pupil::PupilFunction) where {T<:AbstractFloat}
    plan = zernike_propagation_plan(wfs)
    sub = div(_pupil_resolution(pupil), plan.pupil_samples)
    bin2d!(nominal, input, sub)
    binning = zernike_acquisition_plan(wfs).binning
    if binning == 1
        copyto!(out, nominal)
    else
        bin2d!(out, nominal, binning)
    end
    return out
end

function sample_zernike_amplitude_squared!(out::AbstractMatrix{T},
    nominal::AbstractMatrix{T}, wfs::ZernikeWFS,
    amplitude::AbstractMatrix{T}, pupil::PupilFunction) where {
    T<:AbstractFloat,
}
    plan = zernike_propagation_plan(wfs)
    sub = div(_pupil_resolution(pupil), plan.pupil_samples)
    bin2d_abs2!(nominal, amplitude, sub)
    binning = zernike_acquisition_plan(wfs).binning
    if binning == 1
        copyto!(out, nominal)
    else
        bin2d!(out, nominal, binning)
    end
    return out
end

function zernike_pupil_intensity!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS")
    ensure_zernike_buffers!(wfs, pupil)
    propagation = zernike_propagation_workspace(wfs)
    T = eltype(zernike_acquisition_products(wfs).frame)
    n = _pupil_resolution(pupil)
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(
        photon_irradiance(src) * (_pupil_diameter_m(pupil) / _pupil_resolution(pupil))^2
    ))
    amplitude = pupil.amplitude
    fill!(propagation.field, zero(eltype(propagation.field)))
    @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] = amp_scale * amplitude *
        cispi(opd_to_cycles * pupil.opd)
    copyto!(propagation.focal_field, propagation.field)
    @. propagation.focal_field *= propagation.phasor
    execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
    @. propagation.focal_field *= propagation.phase_mask
    copyto!(propagation.pupil_field, propagation.focal_field)
    execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    @views @. propagation.pupil_intensity =
        abs2(propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    return propagation.pupil_intensity
end

function zernike_normalization(normalization::MeanValidFluxNormalization,
    wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    frame::AbstractMatrix)
    return zernike_normalization(normalization, wfs, pupil, src, frame,
        one(eltype(frame)))
end

function zernike_normalization(normalization::WFSNormalization,
    wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    frame::AbstractMatrix{T}, normalization_scale::Real) where {T<:AbstractFloat}
    return zernike_normalization(execution_style(frame), normalization,
        wfs, pupil, src, frame, T(normalization_scale))
end

function zernike_normalization(normalization::IncidenceFluxNormalization,
    wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    frame::AbstractMatrix)
    return zernike_normalization(normalization, wfs, pupil, src, frame,
        one(eltype(frame)))
end

@inline zernike_normalization_count(wfs::ZernikeWFS) =
    length(zernike_estimator_workspace(wfs).valid_signal_indices)

@inline function zernike_incidence_sample_scale(pupil::PupilFunction,
    src::AbstractSource, normalization_scale::T) where {T<:AbstractFloat}
    irradiance = T(_require_physical_photon_irradiance(src,
        "ZernikeWFS incidence normalization"))
    pupil_sample = T(_pupil_diameter_m(pupil)) / T(_pupil_resolution(pupil))
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
    wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    frame::AbstractMatrix{T}, ::T) where {T<:AbstractFloat}
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    state = zernike_estimator_state(wfs)
    summed = masked_sum2d(ScalarCPUStyle(), frame, state.valid_mask)
    return finalize_zernike_normalization(normalization, summed,
        inv(T(count)))
end

function zernike_normalization(::ScalarCPUStyle,
    normalization::IncidenceFluxNormalization,
    wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    frame::AbstractMatrix{T}, normalization_scale::T) where {T<:AbstractFloat}
    scale = zernike_incidence_sample_scale(pupil, src, normalization_scale)
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    sample_binning = div(_pupil_resolution(pupil),
        size(workspace.normalization_frame, 1))
    bin2d_abs2!(workspace.normalization_frame, pupil.amplitude,
        sample_binning)
    summed = masked_sum2d(ScalarCPUStyle(),
        workspace.normalization_frame, state.valid_mask)
    return finalize_zernike_normalization(normalization, summed,
        scale / T(count))
end

@inline function queue_zernike_masked_sum!(phase::KernelLaunchPhase,
    wfs::ZernikeWFS, frame::AbstractMatrix)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    n1, n2 = size(state.valid_mask)
    queue_kernel!(phase, zernike_masked_row_sum_kernel!,
        workspace.normalization_partials, frame, state.valid_mask,
        n1, n2; ndrange=n1)
    queue_kernel!(phase, zernike_finalize_normalization_sum_kernel!,
        workspace.normalization_sum, workspace.normalization_partials,
        n1; ndrange=1)
    return nothing
end

function queue_zernike_normalization!(phase::KernelLaunchPhase,
    ::MeanValidFluxNormalization, wfs::ZernikeWFS, ::PupilFunction,
    ::AbstractSource, frame::AbstractMatrix{T}, ::T) where {T<:AbstractFloat}
    queue_zernike_masked_sum!(phase, wfs, frame)
    return inv(T(zernike_normalization_count(wfs))), true
end

function queue_zernike_normalization!(phase::KernelLaunchPhase,
    ::IncidenceFluxNormalization, wfs::ZernikeWFS, pupil::PupilFunction,
    src::AbstractSource, ::AbstractMatrix{T},
    normalization_scale::T) where {T<:AbstractFloat}
    workspace = zernike_estimator_workspace(wfs)
    sample_binning = div(_pupil_resolution(pupil),
        size(workspace.normalization_frame, 1))
    n1, n2 = size(workspace.normalization_frame)
    queue_kernel!(phase, bin2d_abs2_kernel!,
        workspace.normalization_frame,
        pupil.amplitude, sample_binning, n1, n2;
        ndrange=(n1, n2))
    queue_zernike_masked_sum!(phase, wfs,
        workspace.normalization_frame)
    scale = zernike_incidence_sample_scale(pupil, src,
        normalization_scale)
    return scale / T(zernike_normalization_count(wfs)), false
end

function zernike_normalization(style::AcceleratorStyle,
    normalization::WFSNormalization, wfs::ZernikeWFS, pupil::PupilFunction,
    src::AbstractSource, frame::AbstractMatrix{T},
    normalization_scale::T) where {T<:AbstractFloat}
    # A host scalar result necessarily synchronizes. The measurement hot path
    # below queues this reduction with signal formation and never calls this
    # scalar inspection seam.
    count = zernike_normalization_count(wfs)
    count == 0 && return one(T)
    phase = begin_kernel_phase(style)
    multiplier, _ = queue_zernike_normalization!(phase, normalization,
        wfs, pupil, src, frame, normalization_scale)
    finish_kernel_phase!(phase)
    workspace = zernike_estimator_workspace(wfs)
    copyto!(workspace.normalization_sum_host,
        workspace.normalization_sum)
    return finalize_zernike_normalization(normalization,
        workspace.normalization_sum_host[1], multiplier)
end

function zernike_signal!(wfs::ZernikeWFS, pupil::PupilFunction, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(wfs, pupil, frame, src, one(T))
end

function zernike_signal!(wfs::ZernikeWFS, pupil::PupilFunction,
    frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::Real) where {T<:AbstractFloat}
    return zernike_signal!(execution_style(frame), wfs, pupil, frame, src,
        T(normalization_scale))
end

function zernike_signal!(::ScalarCPUStyle, wfs::ZernikeWFS, pupil::PupilFunction, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(ScalarCPUStyle(), wfs, pupil, frame, src, one(T))
end

function zernike_signal!(::ScalarCPUStyle, wfs::ZernikeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::T) where {T<:AbstractFloat}
    params = zernike_estimator_params(wfs)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    products = zernike_estimator_products(wfs)
    if size(frame) != size(workspace.signal_2d)
        throw(DimensionMismatchError("ZernikeWFS frame size must match sampled camera frame"))
    end
    norm_factor = zernike_normalization(params.normalization, wfs, pupil,
        src, frame, normalization_scale)
    fill!(workspace.signal_2d, zero(T))
    if !usable_wfs_normalization(norm_factor)
        fill!(products.signal, zero(T))
        return products.signal
    end
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        if state.valid_mask[i, j]
            workspace.signal_2d[i, j] = frame[i, j] / norm_factor -
                state.reference_signal_2d[i, j]
        end
    end
    @inbounds for idx in eachindex(workspace.valid_signal_indices)
        products.signal[idx] =
            workspace.signal_2d[workspace.valid_signal_indices[idx]]
    end
    return products.signal
end

function zernike_signal!(style::AcceleratorStyle, wfs::ZernikeWFS, pupil::PupilFunction, frame::AbstractMatrix{T}, src::AbstractSource) where {T<:AbstractFloat}
    return zernike_signal!(style, wfs, pupil, frame, src, one(T))
end

function zernike_signal!(style::AcceleratorStyle, wfs::ZernikeWFS,
    pupil::PupilFunction, frame::AbstractMatrix{T}, src::AbstractSource,
    normalization_scale::T) where {T<:AbstractFloat}
    params = zernike_estimator_params(wfs)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    products = zernike_estimator_products(wfs)
    if size(frame) != size(workspace.signal_2d)
        throw(DimensionMismatchError("ZernikeWFS frame size must match sampled camera frame"))
    end
    count = zernike_normalization_count(wfs)
    if count == 0
        fill!(workspace.signal_2d, zero(T))
        fill!(products.signal, zero(T))
        return products.signal
    end
    phase = begin_kernel_phase(style)
    normalization_multiplier, clamp_to_epsilon =
        queue_zernike_normalization!(phase, params.normalization,
            wfs, pupil, src, frame, normalization_scale)
    queue_kernel!(phase, zernike_signal_kernel!, workspace.signal_2d,
        frame, state.valid_mask, state.reference_signal_2d,
        workspace.normalization_sum, normalization_multiplier,
        clamp_to_epsilon, size(frame, 1), size(frame, 2);
        ndrange=size(frame))
    queue_kernel!(phase, gather_zernike_signal_kernel!, products.signal,
        workspace.signal_2d, workspace.valid_signal_indices, count;
        ndrange=count)
    finish_kernel_phase!(phase)
    return products.signal
end

function ensure_zernike_calibration!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    products = zernike_estimator_products(wfs)
    propagation = zernike_propagation_workspace(wfs)
    frame = zernike_acquisition_products(wfs).frame
    λ = calibration_wavelength(src, eltype(products.signal))
    sig = pupil_aperture_calibration_signature(pupil,
        calibration_signature(src))
    if calibration_matches(state.calibrated,
        state.calibration_wavelength, λ,
        state.calibration_signature, sig)
        return wfs
    end
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        zernike_pupil_intensity!(wfs, pupil, src)
        sample_zernike_frame!(frame, propagation.nominal_frame, wfs,
            propagation.pupil_intensity, pupil)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        zernike_signal!(wfs, pupil, frame, src)
        copyto!(state.reference_signal_2d, workspace.signal_2d)
        fill!(workspace.signal_2d, zero(eltype(workspace.signal_2d)))
        fill!(products.signal, zero(eltype(products.signal)))
    finally
        restore_opd!(pupil, opd_saved)
    end
    state.calibrated = true
    state.calibration_wavelength = λ
    state.calibration_signature = sig
    state.calibration_revision += UInt(1)
    return wfs
end

@inline function ensure_zernike_calibration!(wfs::ZernikeWFS,
    pupil::PupilFunction, src::AbstractSource, ::AbstractDetector)
    return ensure_zernike_calibration!(wfs, pupil, src)
end

function ensure_zernike_calibration!(wfs::ZernikeWFS, pupil::PupilFunction,
    src::AbstractSource, det::Detector)
    state = zernike_estimator_state(wfs)
    workspace = zernike_estimator_workspace(wfs)
    products = zernike_estimator_products(wfs)
    propagation = zernike_propagation_workspace(wfs)
    camera_frame = zernike_acquisition_products(wfs).frame
    T = eltype(products.signal)
    λ = calibration_wavelength(src, T)
    sig = detector_calibration_signature(det,
        pupil_aperture_calibration_signature(pupil,
            calibration_signature(src)))
    if calibration_matches(state.calibrated,
        state.calibration_wavelength, λ,
        state.calibration_signature, sig)
        return wfs
    end

    require_whole_capture_idle(det)
    update_valid_mask!(wfs, pupil)
    opd_saved = save_zero_opd!(pupil)
    try
        zernike_pupil_intensity!(wfs, pupil, src)
        sample_zernike_frame!(camera_frame, propagation.nominal_frame, wfs,
            propagation.pupil_intensity, pupil)
        frame = detector_calibration_frame!(det, camera_frame, src)
        size(frame) == size(workspace.signal_2d) || throw(
            InvalidConfiguration(
                "ZernikeWFS detector sampling and binning must preserve " *
                "the sampled camera frame size"))
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        normalization_scale = wfs_detector_incidence_scale(det, src,
            eltype(frame))
        zernike_signal!(wfs, pupil, frame, src, normalization_scale)
        copyto!(state.reference_signal_2d, workspace.signal_2d)
        fill!(workspace.signal_2d, zero(eltype(workspace.signal_2d)))
        fill!(products.signal, zero(eltype(products.signal)))
    finally
        restore_opd!(pupil, opd_saved)
    end
    state.calibrated = true
    state.calibration_wavelength = λ
    state.calibration_signature = sig
    state.calibration_revision += UInt(1)
    return wfs
end

function measure!(::Diffractive, wfs::ZernikeWFS, pupil::PupilFunction)
    throw(InvalidConfiguration("Diffractive ZernikeWFS requires a source; call measure!(wfs, pupil, src)."))
end

function measure!(wfs::ZernikeWFS, pupil::PupilFunction)
    return measure!(sensing_mode(wfs), wfs, pupil)
end

function measure!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, pupil, src)
end

function measure!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::ZernikeWFS, pupil::PupilFunction, ast::Asterism)
    throw(InvalidConfiguration("ZernikeWFS asterism support is not implemented in the Phase 1 MVP"))
end

function measure!(wfs::ZernikeWFS, pupil::PupilFunction, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration("ZernikeWFS asterism support is not implemented in the Phase 1 MVP"))
end

function measure!(::Diffractive, wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    ensure_zernike_calibration!(wfs, pupil, src)
    propagation = zernike_propagation_workspace(wfs)
    frame = zernike_acquisition_products(wfs).frame
    zernike_pupil_intensity!(wfs, pupil, src)
    sample_zernike_frame!(frame, propagation.nominal_frame, wfs,
        propagation.pupil_intensity, pupil)
    return zernike_signal!(wfs, pupil, frame, src)
end

function measure!(::Diffractive, wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_zernike_calibration!(wfs, pupil, src, det)
    propagation = zernike_propagation_workspace(wfs)
    camera_frame = zernike_acquisition_products(wfs).frame
    zernike_pupil_intensity!(wfs, pupil, src)
    sample_zernike_frame!(camera_frame, propagation.nominal_frame, wfs,
        propagation.pupil_intensity, pupil)
    capture!(det, camera_frame, src; rng=rng)
    size(output_frame(det)) == size(camera_frame) ||
        throw(InvalidConfiguration("ZernikeWFS detector output size must match the sampled camera frame"))
    frame = output_frame(det)
    normalization_scale = wfs_detector_incidence_scale(det, src,
        eltype(frame))
    return zernike_signal!(wfs, pupil, frame, src, normalization_scale)
end

@inline slopes(wfs::ZernikeWFS) = zernike_estimator_products(wfs).signal
@inline valid_subaperture_mask(wfs::ZernikeWFS) =
    zernike_estimator_state(wfs).valid_mask
@inline reference_signal(wfs::ZernikeWFS) =
    zernike_estimator_state(wfs).reference_signal_2d
@inline wfs_calibration_signature(wfs::ZernikeWFS) =
    zernike_estimator_state(wfs).calibration_signature

@inline supports_prepared_runtime(::ZernikeWFS, src::AbstractSource) =
    is_leaf_source(src)
@inline supports_detector_output(::ZernikeWFS, ::AbstractDetector) = true

@inline function prepare_runtime_wfs!(wfs::ZernikeWFS, pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "ZernikeWFS runtime preparation")
    ensure_zernike_calibration!(wfs, pupil, src)
    return wfs
end

include("zernike/stages.jl")
