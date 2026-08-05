#
# Bi-O-edge wavefront sensing
#
# Bi-O-edge is implemented here as a four-edge/Foucault-style pupil-plane sensor.
# The diffractive path:
#
# 1. propagates the pupil field to the focal plane
# 2. applies four complementary knife-edge masks
# 3. propagates back to pupil intensity images
# 4. combines those edge images into x/y differential signals
#
# Modulation, detector binning, and asterism batching follow the same pattern
# as the Pyramid implementation, but the focal-plane filtering is performed
# with the Bi-O-edge mask family rather than a pyramid phase ramp.
#
@kernel function edge_mask_kernel!(mask, pupil, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        if @inbounds pupil[i, j]
            neighbor = false
            @inbounds for di in -1:1, dj in -1:1
                ii = i + di
                jj = j + dj
                if ii < 1 || ii > n || jj < 1 || jj > n || !pupil[ii, jj]
                    neighbor = true
                end
            end
            mask[i, j] = neighbor
        else
            @inbounds mask[i, j] = false
        end
    end
end

@kernel function bi_o_edge_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end

@kernel function bi_o_edge_masks_kernel!(masks, one_c, zero_c, half::Int, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        left = j <= half
        top = i <= half
        @inbounds masks[i, j, 1] = left ? one_c : zero_c
        @inbounds masks[i, j, 2] = left ? zero_c : one_c
        @inbounds masks[i, j, 3] = top ? one_c : zero_c
        @inbounds masks[i, j, 4] = top ? zero_c : one_c
    end
end

@kernel function bin_edge_mask_kernel!(out, mask, binning::Int, n_out::Int, m_out::Int)
    i, j = @index(Global, NTuple)
    if i <= n_out && j <= m_out
        val = false
        @inbounds for ii in 1:binning, jj in 1:binning
            val |= mask[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        @inbounds out[i, j] = val
    end
end

@kernel function gather_bi_o_edge_slopes_kernel!(slopes, signal_2d, valid_signal_indices, count::Int, y_offset::Int)
    idx = @index(Global, Linear)
    if idx <= count
        src = @inbounds valid_signal_indices[idx]
        @inbounds begin
            slopes[idx] = signal_2d[src]
            slopes[idx + count] = signal_2d[src + y_offset]
        end
    end
end

"""Immutable differential-estimator configuration."""
struct BiOEdgeEstimatorParams{T<:AbstractFloat,N<:WFSNormalization}
    pupil_samples::Int
    pupil_resolution::Int
    pupil_diameter_m::T
    threshold::T
    light_ratio::T
    normalization::N
end

"""Run-immutable numerical contract for Bi-O-edge propagation."""
struct BiOEdgePropagationPlan{M<:BiOEdgeAmplitudeMask,T<:AbstractFloat}
    amplitude_mask::M
    pupil_samples::Int
    binning::Int
    numeric_type::Type{T}
end

"""
Backend-bound FFT handles, caches, and replaceable single-writer scratch for
Bi-O-edge propagation. No field is a caller-visible optical product.
"""
mutable struct BiOEdgePropagationWorkspace{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    C3<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}}}
    field::C
    focal_field::C
    pupil_field::C
    bi_o_edge_masks::C3
    phasor::C
    intensity::R
    temp::R
    scratch::R
    asterism_stack::RS
    fft_buffer::C
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::Kf
    lgs_kernel_tag::UInt
    effective_resolution::Int
    asterism_capacity::Int
    revision::UInt
end

"""Exact plan/workspace owner for one Bi-O-edge propagation execution."""
struct PreparedBiOEdgePropagation{
    P<:BiOEdgePropagationPlan,W<:BiOEdgePropagationWorkspace}
    plan::P
    workspace::W
end

@inline bi_o_edge_propagation_plan(
    propagation::PreparedBiOEdgePropagation) = propagation.plan
@inline bi_o_edge_propagation_workspace(
    propagation::PreparedBiOEdgePropagation) = propagation.workspace

"""A physically distinct Bi-O-edge front end with prepared modulation."""
struct BiOEdgeOpticalFrontEnd{O<:BiOEdgeAmplitudeMask,M,C,P,S}
    amplitude_mask::O
    modulation::M
    calibration_modulation::C
    propagation::P
    pupil_samples::Int
    binning::Int
    source::S
end

"""Run-immutable family acquisition contract."""
struct BiOEdgeAcquisitionPlan
    binning::Int
end

"""Derived native sampling metadata for convenience-frame acquisition."""
mutable struct BiOEdgeAcquisitionWorkspace
    nominal_detector_resolution::Int
end

"""Caller-visible convenience-frame product."""
mutable struct BiOEdgeAcquisitionProducts{T<:AbstractFloat,
    R<:AbstractMatrix{T}}
    frame::R
end

struct BiOEdgeDetectorAcquisition{P,W,PR}
    plan::P
    workspace::W
    products::PR
end

"""Persistent support and calibration state for Bi-O-edge estimation."""
mutable struct BiOEdgeEstimatorState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},V<:AbstractVector{T},
    R<:AbstractMatrix{T}}
    valid_mask::A
    edge_mask::A
    optical_gain::V
    valid_i4q::A
    reference_signal_2d::R
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

"""Replaceable single-writer scratch for Bi-O-edge estimation."""
mutable struct BiOEdgeEstimatorWorkspace{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},V<:AbstractVector{T},I<:AbstractVector{Int},
    R<:AbstractMatrix{T}}
    valid_i4q_host::Matrix{Bool}
    valid_signal::A
    valid_signal_indices::I
    valid_signal_indices_host::Vector{Int}
    valid_signal_count::Int
    valid_flux_sum_buffer::V
    valid_flux_sum_host::Vector{T}
    valid_flux_i4q_host::Matrix{T}
    flux_i4q::R
    signal_2d::R
    binned_phase::R
    edge_mask_binned::A
    binned_resolution::Int
end

"""Caller-visible Bi-O-edge differential-slope product."""
mutable struct BiOEdgeEstimatorProducts{T<:AbstractFloat,
    V<:AbstractVector{T}}
    slopes::V
end

struct BiOEdgeDifferentialEstimator{P<:BiOEdgeEstimatorParams,S,W,PR}
    params::P
    state::S
    workspace::W
    products::PR
end

struct BiOEdgeWFS{M<:SensingMode,F,A,E,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::BiOEdgeWFS{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

@inline bi_o_edge_propagation(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.front_end.propagation
@inline bi_o_edge_propagation_plan(wfs::BiOEdgeWFS{<:Diffractive}) =
    bi_o_edge_propagation_plan(bi_o_edge_propagation(wfs))
@inline bi_o_edge_propagation_workspace(wfs::BiOEdgeWFS{<:Diffractive}) =
    bi_o_edge_propagation_workspace(bi_o_edge_propagation(wfs))
@inline bi_o_edge_propagation_workspace(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end.propagation)
@inline four_pupil_propagation_workspace(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end)
@inline bi_o_edge_acquisition_plan(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.acquisition.plan
@inline bi_o_edge_acquisition_workspace(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.acquisition.workspace
@inline bi_o_edge_acquisition_products(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.acquisition.products
@inline bi_o_edge_estimator_state(wfs::BiOEdgeWFS) = wfs.estimator.state
@inline bi_o_edge_estimator_workspace(wfs::BiOEdgeWFS) =
    wfs.estimator.workspace
@inline bi_o_edge_estimator_products(wfs::BiOEdgeWFS) =
    wfs.estimator.products
@inline bi_o_edge_amplitude_mask(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.front_end.amplitude_mask
@inline bi_o_edge_operating_modulation(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.front_end.modulation
@inline bi_o_edge_calibration_modulation(wfs::BiOEdgeWFS{<:Diffractive}) =
    wfs.front_end.calibration_modulation

"""
    BiOEdgeWFS(tel; ...)

Construct a Bi-O-edge wavefront sensor.

The diffractive model forms four edge-filtered pupil images using complementary
focal-plane Bi-O-edge masks. Slopes are then built from the resulting
left/right and top/bottom differential signals after optional binning and
modulation averaging.
`modulation_phase_offset_rad` selects the circular modulation quadrature
origin in radians.
"""
function BiOEdgeWFS(tel::Telescope; pupil_samples::Int, threshold::Real=0.1,
    light_ratio::Real=0.0, modulation::Real=0.0, modulation_points::Union{Int,Nothing}=nothing,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    extra_modulation_factor::Int=0,
    modulation_phase_offset_rad::Real=0.0,
    user_modulation_path=nothing,
    grey_width::Real=0.0, grey_length=false,
    diffraction_padding::Int=2, psf_centering::Bool=true, n_pix_separation=nothing,
    n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel))

    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    pupil_samples >= 1 || throw(InvalidConfiguration(
        "pupil_samples must be >= 1"))
    if tel.params.resolution % pupil_samples != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by pupil_samples"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    if pupil_samples % binning != 0
        throw(InvalidConfiguration(
            "Bi-O-edge binning must evenly divide pupil_samples"))
    end
    grey_length_val = grey_length === false ? false : T(grey_length)
    typed_modulation_phase_offset_rad = T(modulation_phase_offset_rad)
    isfinite(typed_modulation_phase_offset_rad) || throw(
        InvalidConfiguration(
            "Bi-O-edge modulation_phase_offset_rad must be finite"))
    estimator_params = BiOEdgeEstimatorParams{T,typeof(normalization)}(
        pupil_samples,
        tel.params.resolution,
        T(tel.params.diameter),
        T(threshold),
        T(light_ratio),
        normalization)
    amplitude_mask = BiOEdgeAmplitudeMask{T}(
        T(grey_width),
        grey_length_val,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge)
    operating_policy = legacy_modulation_policy(T(modulation),
        modulation_points, extra_modulation_factor,
        typed_modulation_phase_offset_rad,
        user_modulation_path)
    calibration_policy = calibration_modulation_policy(operating_policy,
        T(calib_modulation), typed_modulation_phase_offset_rad)
    valid_mask = backend{Bool}(undef, pupil_samples, pupil_samples)
    edge_mask = backend{Bool}(undef, size(pupil_mask(tel)))
    slopes = backend{T}(undef, 2 * pupil_samples * pupil_samples)
    fill!(slopes, zero(T))
    n_pix_signal = div(pupil_samples, binning)
    valid_i4q = backend{Bool}(undef, n_pix_signal, n_pix_signal)
    valid_i4q_host = Matrix{Bool}(undef, size(valid_i4q)...)
    valid_signal = backend{Bool}(undef, 2 * n_pix_signal, n_pix_signal)
    valid_signal_indices = backend{Int}(undef, length(valid_i4q))
    valid_signal_indices_host = Vector{Int}(undef, length(valid_signal_indices))
    signal_2d = backend{T}(undef, 2 * n_pix_signal, n_pix_signal)
    reference_signal_2d = similar(signal_2d)
    fill!(reference_signal_2d, zero(T))
    valid_flux_sum_buffer = backend{T}(undef, 1)
    valid_flux_sum_host = Vector{T}(undef, 1)
    valid_flux_i4q_host = Matrix{T}(undef, n_pix_signal, n_pix_signal)
    flux_i4q = backend{T}(undef, n_pix_signal, n_pix_signal)
    binned_phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    edge_mask_binned = similar(edge_mask)
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    estimator_state = BiOEdgeEstimatorState(valid_mask, edge_mask,
        optical_gain, valid_i4q, reference_signal_2d, false, zero(T), UInt(0),
        UInt(0))
    estimator_workspace = BiOEdgeEstimatorWorkspace(valid_i4q_host,
        valid_signal, valid_signal_indices, valid_signal_indices_host, 0,
        valid_flux_sum_buffer, valid_flux_sum_host, valid_flux_i4q_host,
        flux_i4q, signal_2d, binned_phase, edge_mask_binned,
        tel.params.resolution)
    estimator_products = BiOEdgeEstimatorProducts(slopes)
    estimator = BiOEdgeDifferentialEstimator(estimator_params,
        estimator_state, estimator_workspace, estimator_products)
    front_end, acquisition = prepare_bi_o_edge_mode(mode, backend, T, tel,
        amplitude_mask, operating_policy, calibration_policy, pupil_samples,
        binning)
    wfs = BiOEdgeWFS{
        typeof(mode),typeof(front_end),typeof(acquisition),typeof(estimator),
        typeof(selector),
    }(front_end, acquisition, estimator)
    initialize_bi_o_edge_masks!(wfs, tel)
    prepare_bi_o_edge_front_end!(mode, wfs)
    return wfs
end

@inline prepare_bi_o_edge_mode(::Geometric, backend, ::Type{T}, tel,
    amplitude_mask, operating_policy, calibration_policy, pupil_samples,
    binning) where {T} = (nothing, nothing)

function prepare_bi_o_edge_mode(::Diffractive, backend, ::Type{T}, tel,
    amplitude_mask, operating_policy, calibration_policy, pupil_samples,
    binning) where {T<:AbstractFloat}
    pad = tel.params.resolution * amplitude_mask.diffraction_padding
    if amplitude_mask.n_pix_separation !== nothing
        edge = amplitude_mask.n_pix_edge === nothing ?
            div(amplitude_mask.n_pix_separation, 2) : amplitude_mask.n_pix_edge
        pad = Int(round((2 * pupil_samples +
            amplitude_mask.n_pix_separation + 2 * edge) *
            tel.params.resolution / pupil_samples))
    end
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    masks = backend{Complex{T}}(undef, pad, pad, 4)
    phasor = similar(field)
    intensity = backend{T}(undef, 2 * pad, 2 * pad)
    temp = backend{T}(undef, pad, pad)
    scratch = similar(temp)
    asterism_stack = backend{T}(undef, 2 * pad, 2 * pad, 1)
    fft_buffer = similar(field)
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    propagation_plan = BiOEdgePropagationPlan(
        amplitude_mask, pupil_samples, binning, T)
    propagation_workspace = BiOEdgePropagationWorkspace(
        field, focal_field, pupil_field,
        masks, phasor, intensity, temp, scratch, asterism_stack, fft_buffer,
        fft_plan, ifft_plan, elongation_kernel, lgs_kernel_fft, UInt(0), pad,
        1, UInt(0))
    propagation = PreparedBiOEdgePropagation(
        propagation_plan, propagation_workspace)
    prepared_modulation = prepare_focal_plane_modulation(operating_policy,
        tel.params.resolution, field, T)
    prepared_calibration = prepare_focal_plane_modulation(
        calibration_policy, tel.params.resolution, field, T)
    front_end = BiOEdgeOpticalFrontEnd(amplitude_mask, prepared_modulation,
        prepared_calibration, propagation, pupil_samples, binning, nothing)
    nominal = max(1,
        round(Int, pupil_samples * pad / tel.params.resolution))
    camera_frame = backend{T}(undef, 2 * nominal, 2 * nominal)
    acquisition = BiOEdgeDetectorAcquisition(BiOEdgeAcquisitionPlan(binning),
        BiOEdgeAcquisitionWorkspace(nominal),
        BiOEdgeAcquisitionProducts(camera_frame))
    return front_end, acquisition
end

@inline prepare_bi_o_edge_front_end!(::Geometric, ::BiOEdgeWFS) = nothing

function prepare_bi_o_edge_front_end!(::Diffractive, wfs::BiOEdgeWFS)
    build_bi_o_edge_phasor!(bi_o_edge_propagation_workspace(wfs).phasor)
    build_bi_o_edge_masks!(wfs)
    return nothing
end

function BiOEdgeOpticalFrontEnd(sensor::BiOEdgeWFS{<:Diffractive},
    source=nothing)
    front_end = sensor.front_end
    return BiOEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.calibration_modulation,
        front_end.propagation, front_end.pupil_samples, front_end.binning,
        source)
end

@inline function bi_o_edge_front_end_with_source(
    front_end::BiOEdgeOpticalFrontEnd, source)
    return BiOEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.calibration_modulation,
        front_end.propagation, front_end.pupil_samples, front_end.binning,
        source)
end

function BiOEdgeOpticalFrontEnd(::BiOEdgeWFS{<:Geometric}, source=nothing)
    throw(WFSPreparationError(:wfs_optics, :unsupported,
        "geometric Bi-O-edge sensing uses DirectMeasurementPath and has no optical front end"))
end

sensing_mode(::BiOEdgeWFS{M}) where {M} = M()

function initialize_bi_o_edge_masks!(wfs::BiOEdgeWFS,
    tel::Telescope)
    set_valid_subapertures!(wfs.estimator.state.valid_mask,
        pupil_mask(tel), wfs.estimator.params.threshold)
    Base.require_one_based_indexing(wfs.estimator.state.edge_mask,
        pupil_mask(tel))
    _update_edge_mask!(execution_style(wfs.estimator.state.edge_mask),
        wfs.estimator.state.edge_mask, pupil_mask(tel),
        tel.params.resolution)
    return wfs
end

function update_valid_mask!(wfs::BiOEdgeWFS, pupil::PupilFunction)
    set_valid_subapertures!(wfs.estimator.state.valid_mask,
        pupil.support, wfs.estimator.params.threshold)
    return wfs
end

function update_edge_mask!(wfs::BiOEdgeWFS, pupil::PupilFunction)
    Base.require_one_based_indexing(wfs.estimator.state.edge_mask, pupil.support)
    _update_edge_mask!(execution_style(wfs.estimator.state.edge_mask),
        wfs.estimator.state.edge_mask, pupil.support,
        _pupil_resolution(pupil))
    return wfs
end

function _update_edge_mask!(::ScalarCPUStyle, mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool}, n::Int)
    @inbounds for i in 1:n, j in 1:n
        if pupil[i, j]
            neighbor = false
            for di in -1:1, dj in -1:1
                ii = i + di
                jj = j + dj
                if ii < 1 || ii > n || jj < 1 || jj > n || !pupil[ii, jj]
                    neighbor = true
                end
            end
            mask[i, j] = neighbor
        else
            mask[i, j] = false
        end
    end
    return mask
end

function _update_edge_mask!(style::AcceleratorStyle, mask::AbstractMatrix{Bool}, pupil::AbstractMatrix{Bool}, n::Int)
    if gpu_backend_name(typeof(mask)) === :amdgpu
        host_pupil = Array(pupil)
        host_mask = Matrix{Bool}(undef, size(mask))
        _update_edge_mask!(ScalarCPUStyle(), host_mask, host_pupil, n)
        copyto!(mask, host_mask)
        return mask
    end
    launch_kernel!(style, edge_mask_kernel!, mask, pupil, n; ndrange=size(mask))
    return mask
end

@inline ensure_bi_o_edge_lgs_kernel!(::NoSodiumLayerProfileStyle, wfs::BiOEdgeWFS, pupil::PupilFunction, src::LGSSource) = wfs
@inline ensure_bi_o_edge_lgs_kernel!(::SampledSodiumLayerProfileStyle, wfs::BiOEdgeWFS, pupil::PupilFunction, src::LGSSource) =
    ensure_lgs_kernel!(wfs, pupil, src)

@inline function apply_bi_o_edge_sodium_layer_profile!(::NoSodiumLayerProfileStyle, wfs::BiOEdgeWFS, src::LGSSource,
    lgs_fft_buffer, lgs_ifft_buffer)
    propagation = bi_o_edge_propagation_workspace(wfs)
    propagation.elongation_kernel = apply_elongation!(
        propagation.temp,
        lgs_elongation_factor(src),
        propagation.scratch,
        propagation.elongation_kernel,
    )
    return propagation.temp
end

@inline function apply_bi_o_edge_sodium_layer_profile!(::SampledSodiumLayerProfileStyle, wfs::BiOEdgeWFS, src::LGSSource,
    lgs_fft_buffer, lgs_ifft_buffer)
    propagation = bi_o_edge_propagation_workspace(wfs)
    apply_lgs_convolution!(
        propagation.temp,
        propagation.lgs_kernel_fft,
        lgs_fft_buffer,
        propagation.fft_plan,
        lgs_ifft_buffer,
        propagation.ifft_plan,
    )
    return propagation.temp
end

function build_bi_o_edge_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_bi_o_edge_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_bi_o_edge_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_bi_o_edge_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, bi_o_edge_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_bi_o_edge_masks!(wfs::BiOEdgeWFS)
    masks = bi_o_edge_propagation_workspace(wfs).bi_o_edge_masks
    copyto!(masks, host_bi_o_edge_masks(wfs))
    return masks
end

function host_bi_o_edge_masks(wfs::BiOEdgeWFS)
    T = eltype(bi_o_edge_estimator_products(wfs).slopes)
    n = size(bi_o_edge_propagation_workspace(wfs).bi_o_edge_masks, 1)
    host = Array{Complex{T}}(undef, n, n, 4)
    build_bi_o_edge_masks_host!(host, wfs)
    return host
end

function build_bi_o_edge_masks_host!(masks::AbstractArray{Complex{T},3}, wfs::BiOEdgeWFS) where {T<:AbstractFloat}
    n = size(masks, 1)
    half = n ÷ 2
    bw = zeros(T, n)
    bw[1:half] .= one(T)
    r = round(Int, wfs.front_end.amplitude_mask.diffraction_padding * wfs.front_end.amplitude_mask.grey_width)
    if r > 0
        gradient = vcat(segment_values(one(T), T(0.5), r), segment_values(T(0.5), zero(T), r))
        lo = max(1, half - r + 1)
        hi = min(n, half + r)
        bw[lo:hi] .= gradient[1:(hi - lo + 1)]
    end
    X = repeat(reshape(bw, 1, :), n, 1)
    A = sqrt.(X)
    if wfs.front_end.amplitude_mask.grey_length !== false
        r_grey = wfs.front_end.amplitude_mask.diffraction_padding
        r_length = round(Int, r_grey * wfs.front_end.amplitude_mask.grey_length)
        top_stop = max(1, half - r_length)
        bot_start = min(n + 1, half + r_length + 1)
        if top_stop >= 1
            A[1:top_stop, 1:half] .= one(T)
            A[1:top_stop, half+1:end] .= zero(T)
        end
        if bot_start <= n
            A[bot_start:end, 1:half] .= one(T)
            A[bot_start:end, half+1:end] .= zero(T)
        end
    end
    B = sqrt.(max.(zero(T), one(T) .- A .^ 2))
    C = permutedims(A)
    D = permutedims(B)
    @views begin
        masks[:, :, 1] .= complex.(A, zero(T))
        masks[:, :, 2] .= complex.(B, zero(T))
        masks[:, :, 3] .= complex.(C, zero(T))
        masks[:, :, 4] .= complex.(D, zero(T))
    end
    return masks
end

segment_values(a::T, b::T, n::Int) where {T<:AbstractFloat} =
    n == 1 ? reshape(fill(a, 1), :) : range(a, b; length=n)

function _build_bi_o_edge_masks!(::ScalarCPUStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
    one_c = complex(one(T), zero(T))
    zero_c = complex(zero(T), zero(T))
    n = size(masks, 1)
    half = n ÷ 2
    @inbounds for i in 1:n, j in 1:n
        left = j <= half
        top = i <= half
        masks[i, j, 1] = left ? one_c : zero_c
        masks[i, j, 2] = left ? zero_c : one_c
        masks[i, j, 3] = top ? one_c : zero_c
        masks[i, j, 4] = top ? zero_c : one_c
    end
    return masks
end

function _build_bi_o_edge_masks!(style::AcceleratorStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
    n = size(masks, 1)
    half = n ÷ 2
    one_c = complex(one(T), zero(T))
    zero_c = complex(zero(T), zero(T))
    launch_kernel!(style, bi_o_edge_masks_kernel!, masks, one_c, zero_c, half, n; ndrange=(n, n))
    return masks
end

function ensure_bi_o_edge_buffers!(wfs::BiOEdgeWFS, pad::Int, pupil::PupilFunction)
    propagation = bi_o_edge_propagation_workspace(wfs)
    acquisition = bi_o_edge_acquisition_products(wfs)
    if size(propagation.field) != (pad, pad)
        propagation.revision += UInt(1)
        propagation.field = similar(propagation.field, pad, pad)
        propagation.focal_field = similar(propagation.focal_field, pad, pad)
        propagation.pupil_field = similar(propagation.pupil_field, pad, pad)
        propagation.bi_o_edge_masks = similar(propagation.bi_o_edge_masks,
            pad, pad, 4)
        propagation.phasor = similar(propagation.phasor, pad, pad)
        propagation.intensity = similar(propagation.intensity,
            2 * pad, 2 * pad)
        propagation.temp = similar(propagation.temp, pad, pad)
        propagation.scratch = similar(propagation.scratch, pad, pad)
        acquisition.frame = similar(acquisition.frame, 2 * pad, 2 * pad)
        propagation.asterism_stack = similar(propagation.asterism_stack,
            2 * pad, 2 * pad, propagation.asterism_capacity)
        propagation.fft_buffer = similar(propagation.fft_buffer, pad, pad)
        propagation.fft_plan = plan_fft_backend!(propagation.focal_field)
        propagation.ifft_plan = plan_ifft_backend!(propagation.pupil_field)
        propagation.lgs_kernel_fft = similar(propagation.focal_field,
            eltype(propagation.focal_field), 0, 0)
        propagation.lgs_kernel_tag = UInt(0)
        propagation.effective_resolution = pad
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
        build_bi_o_edge_phasor!(propagation.phasor)
        build_bi_o_edge_masks!(wfs)
    end
    return wfs
end

function ensure_bi_o_edge_asterism_stack!(wfs::BiOEdgeWFS, n_src::Int)
    n_src >= 1 || throw(InvalidConfiguration("asterism source count must be >= 1"))
    propagation = bi_o_edge_propagation_workspace(wfs)
    dims = size(propagation.intensity)
    if size(propagation.asterism_stack, 1) != dims[1] ||
            size(propagation.asterism_stack, 2) != dims[2] ||
            size(propagation.asterism_stack, 3) < n_src
        capacity = max(n_src, propagation.asterism_capacity)
        propagation.asterism_stack = similar(propagation.asterism_stack,
            dims[1], dims[2], capacity)
        propagation.asterism_capacity = capacity
    end
    return propagation.asterism_stack
end

@inline grouped_staging_buffer(wfs::BiOEdgeWFS, out::AbstractMatrix) = bi_o_edge_propagation_workspace(wfs).intensity

function accumulate_bi_o_edge_asterism_intensity!(::ScalarCPUStyle, wfs::BiOEdgeWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_bi_o_edge_asterism_stack!(wfs, count), count)
    intensity = bi_o_edge_propagation_workspace(wfs).intensity
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, intensity,
        stack, ast.sources, bi_o_edge_intensity!, wfs, pupil)
end

function accumulate_bi_o_edge_asterism_intensity!(style::AcceleratorStyle, wfs::BiOEdgeWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_bi_o_edge_asterism_stack!(wfs, count), count)
    intensity = bi_o_edge_propagation_workspace(wfs).intensity
    return accumulate_grouped_sources!(style, wfs, intensity, stack,
        ast.sources, bi_o_edge_intensity!, wfs, pupil)
end

function prepare_bi_o_edge_sampling!(wfs::BiOEdgeWFS, pupil::PupilFunction)
    binning = bi_o_edge_acquisition_plan(wfs).binning
    workspace = bi_o_edge_estimator_workspace(wfs)
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_sub = wfs.estimator.params.pupil_samples
    pad = _pupil_resolution(pupil) * wfs.front_end.amplitude_mask.diffraction_padding
    if wfs.front_end.amplitude_mask.n_pix_separation !== nothing
        edge = wfs.front_end.amplitude_mask.n_pix_edge === nothing ? div(wfs.front_end.amplitude_mask.n_pix_separation, 2) : wfs.front_end.amplitude_mask.n_pix_edge
        pad = Int(round((n_sub * 2 + wfs.front_end.amplitude_mask.n_pix_separation + 2 * edge) * _pupil_resolution(pupil) / n_sub))
    end
    if pad < _pupil_resolution(pupil)
        throw(InvalidConfiguration("bi_o_edge padding must be >= telescope resolution"))
    end
    if pad % binning != 0
        throw(InvalidConfiguration("bi_o_edge binning must evenly divide padded resolution"))
    end
    n = _pupil_resolution(pupil)
    if n % binning != 0
        throw(InvalidConfiguration("binning must evenly divide telescope resolution"))
    end
    ensure_bi_o_edge_buffers!(wfs, pad, pupil)
    n_binned = div(n, binning)
    if n_binned != workspace.binned_resolution
        workspace.binned_phase = similar(workspace.binned_phase,
            n_binned, n_binned)
        workspace.edge_mask_binned = similar(workspace.edge_mask_binned,
            n_binned, n_binned)
        workspace.binned_resolution = n_binned
    end
    if binning > 1
        bin_edge_mask!(workspace.edge_mask_binned,
            wfs.estimator.state.edge_mask, binning)
    end
    return wfs
end

function bin_edge_mask!(out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int)
    Base.require_one_based_indexing(out, mask)
    n, m = size(mask)
    n_out = div(n, binning)
    m_out = div(m, binning)
    if size(out) != (n_out, m_out)
        throw(DimensionMismatchError("edge_mask_binned size mismatch"))
    end
    _bin_edge_mask!(execution_style(out), out, mask, binning, n_out, m_out)
    return out
end

function _bin_edge_mask!(::ScalarCPUStyle, out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int, n_out::Int, m_out::Int)
    @inbounds for i in 1:n_out, j in 1:m_out
        val = false
        for ii in 1:binning, jj in 1:binning
            val |= mask[(i - 1) * binning + ii, (j - 1) * binning + jj]
        end
        out[i, j] = val
    end
    return out
end

function _bin_edge_mask!(style::AcceleratorStyle, out::AbstractMatrix{Bool}, mask::AbstractMatrix{Bool}, binning::Int, n_out::Int, m_out::Int)
    launch_kernel!(style, bin_edge_mask_kernel!, out, mask, binning, n_out, m_out; ndrange=size(out))
    return out
end

function sample_bi_o_edge_phase!(wfs::BiOEdgeWFS, phase::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = bi_o_edge_acquisition_plan(wfs).binning
    if binning == 1
        return phase, wfs.estimator.state.edge_mask
    end
    workspace = bi_o_edge_estimator_workspace(wfs)
    bin2d!(workspace.binned_phase, phase, binning)
    workspace.binned_phase ./= binning * binning
    return workspace.binned_phase, workspace.edge_mask_binned
end

@inline resize_bi_o_edge_signal_buffers!(wfs::BiOEdgeWFS,
    frame_rows::Int) = resize_bi_o_edge_signal_buffers!(wfs, frame_rows, 1)

@inline resize_bi_o_edge_signal_buffers!(wfs::BiOEdgeWFS,
    frame_rows::Int, ::AbstractDetector) =
    resize_bi_o_edge_signal_buffers!(wfs, frame_rows)

@inline resize_bi_o_edge_signal_buffers!(wfs::BiOEdgeWFS,
    frame_rows::Int, det::Detector) = resize_bi_o_edge_signal_buffers!(wfs,
    frame_rows, det.params.psf_sampling * det.params.binning)

function resize_bi_o_edge_signal_buffers!(wfs::BiOEdgeWFS, frame_rows::Int,
    detector_reduction::Int)
    acquisition_plan = bi_o_edge_acquisition_plan(wfs)
    acquisition_workspace = bi_o_edge_acquisition_workspace(wfs)
    acquisition_products = bi_o_edge_acquisition_products(wfs)
    state = bi_o_edge_estimator_state(wfs)
    workspace = bi_o_edge_estimator_workspace(wfs)
    nominal = acquisition_workspace.nominal_detector_resolution
    detector_reduction >= 1 || throw(InvalidConfiguration(
        "Bi-O-edge detector sampling reduction must be >= 1"))
    nominal_pixels = max(1,
        round(Int, nominal / (2 * acquisition_plan.binning)))
    nominal_pixels % detector_reduction == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must preserve an integer Bi-O-edge pupil image"))
    n_pixels = div(nominal_pixels, detector_reduction)
    n_pixels >= 1 || throw(InvalidConfiguration(
        "detector sampling and binning removed every Bi-O-edge pupil sample"))
    iseven(frame_rows) || throw(InvalidConfiguration(
        "Bi-O-edge camera frame must have even dimensions for symmetric pupil extraction"))
    frame_rows >= 2 * n_pixels || throw(InvalidConfiguration(
        "Bi-O-edge camera frame does not contain four complete pupil images"))
    calibration_storage_changed = false
    if size(state.valid_i4q) != (n_pixels, n_pixels)
        state.valid_i4q = similar(state.valid_i4q, n_pixels, n_pixels)
        fill!(state.valid_i4q, false)
        calibration_storage_changed = true
    end
    if size(workspace.valid_signal) != (2 * n_pixels, n_pixels)
        workspace.valid_signal = similar(workspace.valid_signal,
            2 * n_pixels, n_pixels)
    end
    if size(workspace.flux_i4q) != (n_pixels, n_pixels)
        workspace.flux_i4q = similar(workspace.flux_i4q,
            n_pixels, n_pixels)
    end
    if size(workspace.signal_2d) != (2 * n_pixels, n_pixels)
        workspace.signal_2d = similar(workspace.signal_2d,
            2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        calibration_storage_changed = true
    elseif size(state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        state.reference_signal_2d = similar(state.reference_signal_2d,
            2 * n_pixels, n_pixels)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        calibration_storage_changed = true
    end
    if size(acquisition_products.frame) != (frame_rows, frame_rows)
        acquisition_products.frame = similar(acquisition_products.frame,
            frame_rows, frame_rows)
    end
    update_bi_o_edge_valid_signal!(wfs)
    if calibration_storage_changed
        state.calibrated = false
        state.calibration_revision += UInt(1)
    end
    return wfs
end

@inline function require_bi_o_edge_frame_geometry(wfs::BiOEdgeWFS,
    frame::AbstractMatrix)
    n_rows, n_cols = size(frame)
    n_rows == n_cols || throw(DimensionMismatchError(
        "Bi-O-edge camera frame must be square"))
    iseven(n_rows) || throw(InvalidConfiguration(
        "Bi-O-edge camera frame must have even dimensions for symmetric pupil extraction"))
    n_pixels = size(bi_o_edge_estimator_workspace(wfs).signal_2d, 2)
    n_rows >= 2 * n_pixels || throw(DimensionMismatchError(
        "Bi-O-edge camera frame does not contain four complete pupil images"))
    return div(n_rows, 2)
end

function sample_bi_o_edge_intensity!(wfs::BiOEdgeWFS, pupil::PupilFunction, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    acquisition_plan = bi_o_edge_acquisition_plan(wfs)
    acquisition_workspace = bi_o_edge_acquisition_workspace(wfs)
    acquisition_products = bi_o_edge_acquisition_products(wfs)
    propagation = bi_o_edge_propagation_workspace(wfs)
    sub = div(_pupil_resolution(pupil), wfs.estimator.params.pupil_samples)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("bi_o_edge intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    if size(acquisition_products.frame) != (n_camera, n_camera)
        acquisition_products.frame = similar(acquisition_products.frame,
            n_camera, n_camera)
    end
    frame = acquisition_products.frame
    acquisition_workspace.nominal_detector_resolution = round(Int,
        wfs.estimator.params.pupil_samples * propagation.effective_resolution /
        _pupil_resolution(pupil))
    if acquisition_plan.binning != 1
        target = div(acquisition_workspace.nominal_detector_resolution,
            acquisition_plan.binning)
        factor = div(size(frame, 1), target)
        if factor < 1 || size(frame, 1) % target != 0
            throw(InvalidConfiguration("bi_o_edge detector binning is not compatible with the sampled frame"))
        end
        if size(acquisition_products.frame) != (target, target)
            acquisition_products.frame = similar(acquisition_products.frame,
                target, target)
        end
        bin2d!(acquisition_products.frame, intensity, sub * factor)
        frame = acquisition_products.frame
    else
        bin2d!(acquisition_products.frame, intensity, sub)
    end
    resize_bi_o_edge_signal_buffers!(wfs, size(frame, 1))
    return frame
end
