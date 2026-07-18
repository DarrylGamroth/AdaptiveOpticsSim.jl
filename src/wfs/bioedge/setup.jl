#
# BioEdge wavefront sensing
#
# BioEdge is implemented here as a four-edge/Foucault-style pupil-plane sensor.
# The diffractive path:
#
# 1. propagates the pupil field to the focal plane
# 2. applies four complementary knife-edge masks
# 3. propagates back to pupil intensity images
# 4. combines those edge images into x/y differential signals
#
# Modulation, detector binning, and asterism batching follow the same pattern
# as the Pyramid implementation, but the focal-plane filtering is performed
# with the BioEdge mask family rather than a pyramid phase ramp.
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

@kernel function bioedge_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end

@kernel function bioedge_masks_kernel!(masks, one_c, zero_c, half::Int, n::Int)
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

@kernel function gather_bioedge_slopes_kernel!(slopes, signal_2d, valid_signal_indices, count::Int, y_offset::Int)
    idx = @index(Global, Linear)
    if idx <= count
        src = @inbounds valid_signal_indices[idx]
        @inbounds begin
            slopes[idx] = signal_2d[src]
            slopes[idx + count] = signal_2d[src + y_offset]
        end
    end
end

"""Immutable physical definition of the four BioEdge amplitude masks."""
struct BioEdgeAmplitudeMask{T<:AbstractFloat}
    grey_width::T
    grey_length::Union{Bool,T}
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
end

"""Immutable differential-estimator configuration."""
struct BioEdgeEstimatorParams{T<:AbstractFloat,N<:WFSNormalization}
    pupil_samples::Int
    pupil_resolution::Int
    pupil_diameter_m::T
    threshold::T
    light_ratio::T
    normalization::N
end

"""Single-writer FFT and scratch state for BioEdge propagation."""
mutable struct PreparedBioEdgePropagation{T<:AbstractFloat,
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
    bioedge_masks::C3
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
end

"""A physically distinct BioEdge front end with prepared modulation."""
struct BioEdgeOpticalFrontEnd{O<:BioEdgeAmplitudeMask,M,C,P,S}
    amplitude_mask::O
    modulation::M
    calibration_modulation::C
    propagation::P
    pupil_samples::Int
    binning::Int
    source::S
end

mutable struct BioEdgeAcquisitionState{T<:AbstractFloat,R<:AbstractMatrix{T}}
    binned_intensity::R
    camera_frame::R
    nominal_detector_resolution::Int
end

struct BioEdgeDetectorAcquisition{S}
    binning::Int
    state::S
end

"""Mutable support, calibration, and output storage for BioEdge estimation."""
mutable struct BioEdgeEstimatorState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},V<:AbstractVector{T},I<:AbstractVector{Int},
    R<:AbstractMatrix{T}}
    valid_mask::A
    edge_mask::A
    slopes::V
    optical_gain::V
    valid_i4q::A
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
    reference_signal_2d::R
    binned_phase::R
    edge_mask_binned::A
    binned_resolution::Int
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

struct BioEdgeDifferentialEstimator{P<:BioEdgeEstimatorParams,S}
    params::P
    state::S
end

struct BioEdgeWFS{M<:SensingMode,F,A,E,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::BioEdgeWFS{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

@inline bioedge_propagation(wfs::BioEdgeWFS{<:Diffractive}) =
    wfs.front_end.propagation
@inline bioedge_amplitude_mask(wfs::BioEdgeWFS{<:Diffractive}) =
    wfs.front_end.amplitude_mask
@inline bioedge_operating_modulation(wfs::BioEdgeWFS{<:Diffractive}) =
    wfs.front_end.modulation
@inline bioedge_calibration_modulation(wfs::BioEdgeWFS{<:Diffractive}) =
    wfs.front_end.calibration_modulation

"""
    BioEdgeWFS(tel; ...)

Construct a BioEdge wavefront sensor.

The diffractive model forms four edge-filtered pupil images using complementary
focal-plane BioEdge masks. Slopes are then built from the resulting
left/right and top/bottom differential signals after optional binning and
modulation averaging.
"""
function BioEdgeWFS(tel::Telescope; pupil_samples::Int, threshold::Real=0.1,
    light_ratio::Real=0.0, modulation::Real=0.0, modulation_points::Union{Int,Nothing}=nothing,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    extra_modulation_factor::Int=0, delta_theta::Real=0.0, user_modulation_path=nothing,
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
            "BioEdge binning must evenly divide pupil_samples"))
    end
    grey_length_val = grey_length === false ? false : T(grey_length)
    estimator_params = BioEdgeEstimatorParams{T,typeof(normalization)}(
        pupil_samples,
        tel.params.resolution,
        T(tel.params.diameter),
        T(threshold),
        T(light_ratio),
        normalization)
    amplitude_mask = BioEdgeAmplitudeMask{T}(
        T(grey_width),
        grey_length_val,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge)
    operating_policy = legacy_modulation_policy(T(modulation),
        modulation_points, extra_modulation_factor, T(delta_theta),
        user_modulation_path)
    calibration_policy = calibration_modulation_policy(operating_policy,
        T(calib_modulation), T(delta_theta))
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
    valid_flux_sum_buffer = backend{T}(undef, 1)
    valid_flux_sum_host = Vector{T}(undef, 1)
    valid_flux_i4q_host = Matrix{T}(undef, n_pix_signal, n_pix_signal)
    flux_i4q = backend{T}(undef, n_pix_signal, n_pix_signal)
    binned_phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    edge_mask_binned = similar(edge_mask)
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    estimator_state = BioEdgeEstimatorState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(valid_signal_indices),
        typeof(signal_2d),
    }(
        valid_mask,
        edge_mask,
        slopes,
        optical_gain,
        valid_i4q,
        valid_i4q_host,
        valid_signal,
        valid_signal_indices,
        valid_signal_indices_host,
        0,
        valid_flux_sum_buffer,
        valid_flux_sum_host,
        valid_flux_i4q_host,
        flux_i4q,
        signal_2d,
        reference_signal_2d,
        binned_phase,
        edge_mask_binned,
        tel.params.resolution,
        false,
        zero(T),
        UInt(0),
        UInt(0),
    )
    estimator = BioEdgeDifferentialEstimator(estimator_params,
        estimator_state)
    front_end, acquisition = prepare_bioedge_mode(mode, backend, T, tel,
        amplitude_mask, operating_policy, calibration_policy, pupil_samples,
        binning)
    wfs = BioEdgeWFS{
        typeof(mode),typeof(front_end),typeof(acquisition),typeof(estimator),
        typeof(selector),
    }(front_end, acquisition, estimator)
    update_valid_mask!(wfs, tel)
    update_edge_mask!(wfs, tel)
    prepare_bioedge_front_end!(mode, wfs)
    return wfs
end

@inline prepare_bioedge_mode(::Geometric, backend, ::Type{T}, tel,
    amplitude_mask, operating_policy, calibration_policy, pupil_samples,
    binning) where {T} = (nothing, nothing)

function prepare_bioedge_mode(::Diffractive, backend, ::Type{T}, tel,
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
    propagation = PreparedBioEdgePropagation(field, focal_field, pupil_field,
        masks, phasor, intensity, temp, scratch, asterism_stack, fft_buffer,
        fft_plan, ifft_plan, elongation_kernel, lgs_kernel_fft, UInt(0), pad,
        1)
    prepared_modulation = prepare_focal_plane_modulation(operating_policy,
        tel.params.resolution, field, T)
    prepared_calibration = prepare_focal_plane_modulation(
        calibration_policy, tel.params.resolution, field, T)
    front_end = BioEdgeOpticalFrontEnd(amplitude_mask, prepared_modulation,
        prepared_calibration, propagation, pupil_samples, binning, nothing)
    binned_intensity = similar(intensity)
    nominal = max(1,
        round(Int, pupil_samples * pad / tel.params.resolution))
    camera_frame = backend{T}(undef, 2 * nominal, 2 * nominal)
    acquisition = BioEdgeDetectorAcquisition(binning,
        BioEdgeAcquisitionState(binned_intensity, camera_frame, nominal))
    return front_end, acquisition
end

@inline prepare_bioedge_front_end!(::Geometric, ::BioEdgeWFS) = nothing

function prepare_bioedge_front_end!(::Diffractive, wfs::BioEdgeWFS)
    build_bioedge_phasor!(wfs.front_end.propagation.phasor)
    build_bioedge_masks!(wfs)
    return nothing
end

function BioEdgeOpticalFrontEnd(sensor::BioEdgeWFS{<:Diffractive},
    source=nothing)
    front_end = sensor.front_end
    return BioEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.calibration_modulation,
        front_end.propagation, front_end.pupil_samples, front_end.binning,
        source)
end

@inline function bioedge_front_end_with_source(
    front_end::BioEdgeOpticalFrontEnd, source)
    return BioEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.calibration_modulation,
        front_end.propagation, front_end.pupil_samples, front_end.binning,
        source)
end

function BioEdgeOpticalFrontEnd(::BioEdgeWFS{<:Geometric}, source=nothing)
    throw(WFSPreparationError(:optical_formation, :unsupported,
        "geometric BioEdge sensing uses DirectMeasurementPath and has no optical front end"))
end

sensing_mode(::BioEdgeWFS{M}) where {M} = M()

function update_valid_mask!(wfs::BioEdgeWFS, tel::Telescope)
    set_valid_subapertures!(wfs.estimator.state.valid_mask, pupil_mask(tel), wfs.estimator.params.threshold)
    return wfs
end

function update_edge_mask!(wfs::BioEdgeWFS, tel::Telescope)
    Base.require_one_based_indexing(wfs.estimator.state.edge_mask, pupil_mask(tel))
    _update_edge_mask!(execution_style(wfs.estimator.state.edge_mask), wfs.estimator.state.edge_mask, pupil_mask(tel), tel.params.resolution)
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

@inline ensure_bioedge_lgs_kernel!(::LGSProfileNone, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource) = wfs
@inline ensure_bioedge_lgs_kernel!(::LGSProfileNaProfile, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource) =
    ensure_lgs_kernel!(wfs, tel, src)
@inline ensure_bioedge_lgs_kernel!(profile, wfs::BioEdgeWFS, tel::Telescope, src::LGSSource) = wfs

@inline function apply_bioedge_lgs_profile!(::LGSProfileNone, wfs::BioEdgeWFS, src::LGSSource,
    lgs_fft_buffer, lgs_ifft_buffer)
    wfs.front_end.propagation.elongation_kernel = apply_elongation!(
        wfs.front_end.propagation.temp,
        lgs_elongation_factor(src),
        wfs.front_end.propagation.scratch,
        wfs.front_end.propagation.elongation_kernel,
    )
    return wfs.front_end.propagation.temp
end

@inline function apply_bioedge_lgs_profile!(::LGSProfileNaProfile, wfs::BioEdgeWFS, src::LGSSource,
    lgs_fft_buffer, lgs_ifft_buffer)
    apply_lgs_convolution!(
        wfs.front_end.propagation.temp,
        wfs.front_end.propagation.lgs_kernel_fft,
        lgs_fft_buffer,
        wfs.front_end.propagation.fft_plan,
        lgs_ifft_buffer,
        wfs.front_end.propagation.ifft_plan,
    )
    return wfs.front_end.propagation.temp
end

@inline apply_bioedge_lgs_profile!(profile, wfs::BioEdgeWFS, src::LGSSource, lgs_fft_buffer, lgs_ifft_buffer) =
    wfs.front_end.propagation.temp

function build_bioedge_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_bioedge_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_bioedge_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_bioedge_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, bioedge_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_bioedge_masks!(wfs::BioEdgeWFS)
    masks = wfs.front_end.propagation.bioedge_masks
    copyto!(masks, host_bioedge_masks(wfs))
    return masks
end

function host_bioedge_masks(wfs::BioEdgeWFS)
    T = eltype(wfs.estimator.state.slopes)
    n = size(wfs.front_end.propagation.bioedge_masks, 1)
    host = Array{Complex{T}}(undef, n, n, 4)
    build_bioedge_masks_host!(host, wfs)
    return host
end

function build_bioedge_masks_host!(masks::AbstractArray{Complex{T},3}, wfs::BioEdgeWFS) where {T<:AbstractFloat}
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

function _build_bioedge_masks!(::ScalarCPUStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
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

function _build_bioedge_masks!(style::AcceleratorStyle, masks::AbstractArray{Complex{T},3}, ::Type{T}) where {T<:AbstractFloat}
    n = size(masks, 1)
    half = n ÷ 2
    one_c = complex(one(T), zero(T))
    zero_c = complex(zero(T), zero(T))
    launch_kernel!(style, bioedge_masks_kernel!, masks, one_c, zero_c, half, n; ndrange=(n, n))
    return masks
end

function ensure_bioedge_buffers!(wfs::BioEdgeWFS, pad::Int, tel::Telescope)
    if size(wfs.front_end.propagation.field) != (pad, pad)
        wfs.front_end.propagation.field = similar(wfs.front_end.propagation.field, pad, pad)
        wfs.front_end.propagation.focal_field = similar(wfs.front_end.propagation.focal_field, pad, pad)
        wfs.front_end.propagation.pupil_field = similar(wfs.front_end.propagation.pupil_field, pad, pad)
        wfs.front_end.propagation.bioedge_masks = similar(wfs.front_end.propagation.bioedge_masks, pad, pad, 4)
        wfs.front_end.propagation.phasor = similar(wfs.front_end.propagation.phasor, pad, pad)
        wfs.front_end.propagation.intensity = similar(wfs.front_end.propagation.intensity, 2 * pad, 2 * pad)
        wfs.front_end.propagation.temp = similar(wfs.front_end.propagation.temp, pad, pad)
        wfs.front_end.propagation.scratch = similar(wfs.front_end.propagation.scratch, pad, pad)
        wfs.acquisition.state.binned_intensity = similar(wfs.acquisition.state.binned_intensity, 2 * pad, 2 * pad)
        wfs.front_end.propagation.asterism_stack = similar(wfs.front_end.propagation.asterism_stack, 2 * pad, 2 * pad, wfs.front_end.propagation.asterism_capacity)
        wfs.front_end.propagation.fft_buffer = similar(wfs.front_end.propagation.fft_buffer, pad, pad)
        wfs.front_end.propagation.fft_plan = plan_fft_backend!(wfs.front_end.propagation.focal_field)
        wfs.front_end.propagation.ifft_plan = plan_ifft_backend!(wfs.front_end.propagation.pupil_field)
        wfs.front_end.propagation.lgs_kernel_fft = similar(wfs.front_end.propagation.focal_field, Complex{eltype(wfs.front_end.propagation.focal_field)}, 0, 0)
        wfs.front_end.propagation.lgs_kernel_tag = UInt(0)
        wfs.front_end.propagation.effective_resolution = pad
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
        build_bioedge_phasor!(wfs.front_end.propagation.phasor)
        build_bioedge_masks!(wfs)
    end
    return wfs
end

function ensure_bioedge_asterism_stack!(wfs::BioEdgeWFS, n_src::Int)
    n_src >= 1 || throw(InvalidConfiguration("asterism source count must be >= 1"))
    dims = size(wfs.front_end.propagation.intensity)
    if size(wfs.front_end.propagation.asterism_stack, 1) != dims[1] || size(wfs.front_end.propagation.asterism_stack, 2) != dims[2] ||
            size(wfs.front_end.propagation.asterism_stack, 3) < n_src
        capacity = max(n_src, wfs.front_end.propagation.asterism_capacity)
        wfs.front_end.propagation.asterism_stack = similar(wfs.front_end.propagation.asterism_stack, dims[1], dims[2], capacity)
        wfs.front_end.propagation.asterism_capacity = capacity
    end
    return wfs.front_end.propagation.asterism_stack
end

@inline grouped_staging_buffer(wfs::BioEdgeWFS, out::AbstractMatrix) = wfs.front_end.propagation.intensity

function accumulate_bioedge_asterism_intensity!(::ScalarCPUStyle, wfs::BioEdgeWFS, tel::Telescope, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_bioedge_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, wfs.front_end.propagation.intensity, stack, ast.sources, bioedge_intensity!, wfs, tel)
end

function accumulate_bioedge_asterism_intensity!(style::AcceleratorStyle, wfs::BioEdgeWFS, tel::Telescope, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_bioedge_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(style, wfs, wfs.front_end.propagation.intensity, stack, ast.sources, bioedge_intensity!, wfs, tel)
end

function prepare_bioedge_sampling!(wfs::BioEdgeWFS, tel::Telescope)
    binning = wfs.acquisition.binning
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_sub = wfs.estimator.params.pupil_samples
    pad = tel.params.resolution * wfs.front_end.amplitude_mask.diffraction_padding
    if wfs.front_end.amplitude_mask.n_pix_separation !== nothing
        edge = wfs.front_end.amplitude_mask.n_pix_edge === nothing ? div(wfs.front_end.amplitude_mask.n_pix_separation, 2) : wfs.front_end.amplitude_mask.n_pix_edge
        pad = Int(round((n_sub * 2 + wfs.front_end.amplitude_mask.n_pix_separation + 2 * edge) * tel.params.resolution / n_sub))
    end
    if pad < tel.params.resolution
        throw(InvalidConfiguration("bioedge padding must be >= telescope resolution"))
    end
    if pad % binning != 0
        throw(InvalidConfiguration("bioedge binning must evenly divide padded resolution"))
    end
    n = tel.params.resolution
    if n % binning != 0
        throw(InvalidConfiguration("binning must evenly divide telescope resolution"))
    end
    ensure_bioedge_buffers!(wfs, pad, tel)
    n_binned = div(n, binning)
    if n_binned != wfs.estimator.state.binned_resolution
        wfs.estimator.state.binned_phase = similar(wfs.estimator.state.binned_phase, n_binned, n_binned)
        wfs.estimator.state.edge_mask_binned = similar(wfs.estimator.state.edge_mask_binned, n_binned, n_binned)
        wfs.estimator.state.binned_resolution = n_binned
    end
    if binning > 1
        bin_edge_mask!(wfs.estimator.state.edge_mask_binned, wfs.estimator.state.edge_mask, binning)
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

function sample_bioedge_phase!(wfs::BioEdgeWFS, phase::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.acquisition.binning
    if binning == 1
        return phase, wfs.estimator.state.edge_mask
    end
    bin2d!(wfs.estimator.state.binned_phase, phase, binning)
    wfs.estimator.state.binned_phase ./= binning * binning
    return wfs.estimator.state.binned_phase, wfs.estimator.state.edge_mask_binned
end

@inline resize_bioedge_signal_buffers!(wfs::BioEdgeWFS,
    frame_rows::Int) = resize_bioedge_signal_buffers!(wfs, frame_rows, 1)

@inline resize_bioedge_signal_buffers!(wfs::BioEdgeWFS,
    frame_rows::Int, ::AbstractDetector) =
    resize_bioedge_signal_buffers!(wfs, frame_rows)

@inline resize_bioedge_signal_buffers!(wfs::BioEdgeWFS,
    frame_rows::Int, det::Detector) = resize_bioedge_signal_buffers!(wfs,
    frame_rows, det.params.psf_sampling * det.params.binning)

function resize_bioedge_signal_buffers!(wfs::BioEdgeWFS, frame_rows::Int,
    detector_reduction::Int)
    nominal = wfs.acquisition.state.nominal_detector_resolution
    detector_reduction >= 1 || throw(InvalidConfiguration(
        "BioEdge detector sampling reduction must be >= 1"))
    nominal_pixels = max(1,
        round(Int, nominal / (2 * wfs.acquisition.binning)))
    nominal_pixels % detector_reduction == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must preserve an integer BioEdge pupil image"))
    n_pixels = div(nominal_pixels, detector_reduction)
    n_pixels >= 1 || throw(InvalidConfiguration(
        "detector sampling and binning removed every BioEdge pupil sample"))
    iseven(frame_rows) || throw(InvalidConfiguration(
        "BioEdge camera frame must have even dimensions for symmetric pupil extraction"))
    frame_rows >= 2 * n_pixels || throw(InvalidConfiguration(
        "BioEdge camera frame does not contain four complete pupil images"))
    calibration_storage_changed = false
    if size(wfs.estimator.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.valid_i4q = similar(wfs.estimator.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.estimator.state.valid_i4q, false)
        calibration_storage_changed = true
    end
    if size(wfs.estimator.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.valid_signal = similar(wfs.estimator.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.flux_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.flux_i4q = similar(
            wfs.estimator.state.flux_i4q, n_pixels, n_pixels)
    end
    if size(wfs.estimator.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.signal_2d = similar(wfs.estimator.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
        calibration_storage_changed = true
    elseif size(wfs.estimator.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
        calibration_storage_changed = true
    end
    if size(wfs.acquisition.state.camera_frame) != (frame_rows, frame_rows)
        wfs.acquisition.state.camera_frame = similar(wfs.acquisition.state.camera_frame, frame_rows, frame_rows)
    end
    update_bioedge_valid_signal!(wfs)
    if calibration_storage_changed
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
    end
    return wfs
end

@inline function require_bioedge_frame_geometry(wfs::BioEdgeWFS,
    frame::AbstractMatrix)
    n_rows, n_cols = size(frame)
    n_rows == n_cols || throw(DimensionMismatchError(
        "BioEdge camera frame must be square"))
    iseven(n_rows) || throw(InvalidConfiguration(
        "BioEdge camera frame must have even dimensions for symmetric pupil extraction"))
    n_pixels = size(wfs.estimator.state.signal_2d, 2)
    n_rows >= 2 * n_pixels || throw(DimensionMismatchError(
        "BioEdge camera frame does not contain four complete pupil images"))
    return div(n_rows, 2)
end

function sample_bioedge_intensity!(wfs::BioEdgeWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    sub = div(tel.params.resolution, wfs.estimator.params.pupil_samples)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("bioedge intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    if size(wfs.acquisition.state.camera_frame) != (n_camera, n_camera)
        wfs.acquisition.state.camera_frame = similar(wfs.acquisition.state.camera_frame, n_camera, n_camera)
    end
    frame = wfs.acquisition.state.camera_frame
    wfs.acquisition.state.nominal_detector_resolution = round(Int, wfs.estimator.params.pupil_samples * wfs.front_end.propagation.effective_resolution / tel.params.resolution)
    if wfs.acquisition.binning != 1
        target = div(wfs.acquisition.state.nominal_detector_resolution, wfs.acquisition.binning)
        factor = div(size(frame, 1), target)
        if factor < 1 || size(frame, 1) % target != 0
            throw(InvalidConfiguration("bioedge detector binning is not compatible with the sampled frame"))
        end
        if size(wfs.acquisition.state.binned_intensity) != (target, target)
            wfs.acquisition.state.binned_intensity = similar(wfs.acquisition.state.binned_intensity, target, target)
        end
        bin2d!(wfs.acquisition.state.binned_intensity, intensity, sub * factor)
        frame = wfs.acquisition.state.binned_intensity
    else
        bin2d!(wfs.acquisition.state.camera_frame, intensity, sub)
    end
    resize_bioedge_signal_buffers!(wfs, size(frame, 1))
    return frame
end
