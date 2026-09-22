#
# Bi-O-edge wavefront sensing
#
# Bi-O-edge is implemented here as a four-edge/Foucault-style pupil-plane sensor.
# The diffractive path:
#
# 1. propagates the pupil field to the focal plane
# 2. applies four complementary knife-edge masks
# 3. propagates back to pupil intensity images
# 4. publishes the four-pupil detector-plane photon-rate map
#
# Modulation, detector binning, and asterism batching follow the same pattern
# as the Pyramid implementation, but the focal-plane filtering is performed
# with the Bi-O-edge mask family rather than a pyramid phase ramp.
#
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
struct BiOEdgeOpticalFrontEnd{O<:BiOEdgeAmplitudeMask,M,P,S}
    amplitude_mask::O
    modulation::M
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

struct BiOEdgeWFS{F,A,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
end

@inline backend(::BiOEdgeWFS{<:Any,<:Any,B}) where {B} = B()

@inline bi_o_edge_propagation(wfs::BiOEdgeWFS) =
    wfs.front_end.propagation
@inline bi_o_edge_propagation_plan(wfs::BiOEdgeWFS) =
    bi_o_edge_propagation_plan(bi_o_edge_propagation(wfs))
@inline bi_o_edge_propagation_workspace(wfs::BiOEdgeWFS) =
    bi_o_edge_propagation_workspace(bi_o_edge_propagation(wfs))
@inline bi_o_edge_propagation_workspace(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end.propagation)
@inline four_pupil_propagation_workspace(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end)
@inline bi_o_edge_acquisition_plan(wfs::BiOEdgeWFS) =
    wfs.acquisition.plan
@inline bi_o_edge_acquisition_workspace(wfs::BiOEdgeWFS) =
    wfs.acquisition.workspace
@inline bi_o_edge_acquisition_products(wfs::BiOEdgeWFS) =
    wfs.acquisition.products
@inline bi_o_edge_amplitude_mask(wfs::BiOEdgeWFS) =
    wfs.front_end.amplitude_mask
@inline bi_o_edge_operating_modulation(wfs::BiOEdgeWFS) =
    wfs.front_end.modulation

"""
    BiOEdgeWFS(tel; ...)

Construct a Bi-O-edge wavefront sensor.

The model forms four edge-filtered pupil images using complementary focal-plane
Bi-O-edge masks and publishes the physical detector-plane photon-rate map.
Detector acquisition and RTC signal estimation are separate owners.
`modulation_phase_offset_rad` selects the circular modulation quadrature
origin in radians.
"""
function BiOEdgeWFS(tel::Telescope; pupil_samples::Int,
    modulation::Real=0.0, modulation_points::Union{Int,Nothing}=nothing,
    extra_modulation_factor::Int=0,
    modulation_phase_offset_rad::Real=0.0,
    user_modulation_path=nothing,
    grey_width::Real=0.0, grey_length=false,
    diffraction_padding::Int=2, psf_centering::Bool=true, n_pix_separation=nothing,
    n_pix_edge=nothing, binning::Int=1,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))

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
    front_end, acquisition = _prepare_bi_o_edge_diffractive_storage(
        backend, T, tel, amplitude_mask, operating_policy, pupil_samples,
        binning)
    wfs = BiOEdgeWFS{
        typeof(front_end),typeof(acquisition),typeof(selector),
    }(front_end, acquisition)
    prepare_bi_o_edge_front_end!(wfs)
    return wfs
end

function _prepare_bi_o_edge_diffractive_storage(backend, ::Type{T}, tel,
    amplitude_mask, operating_policy, pupil_samples,
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
    front_end = BiOEdgeOpticalFrontEnd(amplitude_mask, prepared_modulation,
        propagation, pupil_samples, binning, nothing)
    nominal = max(1,
        round(Int, pupil_samples * pad / tel.params.resolution))
    camera_frame = backend{T}(undef, 2 * nominal, 2 * nominal)
    acquisition = BiOEdgeDetectorAcquisition(BiOEdgeAcquisitionPlan(binning),
        BiOEdgeAcquisitionWorkspace(nominal),
        BiOEdgeAcquisitionProducts(camera_frame))
    return front_end, acquisition
end

function prepare_bi_o_edge_front_end!(wfs::BiOEdgeWFS)
    build_bi_o_edge_phasor!(bi_o_edge_propagation_workspace(wfs).phasor)
    build_bi_o_edge_masks!(wfs)
    return nothing
end

function BiOEdgeOpticalFrontEnd(sensor::BiOEdgeWFS,
    source=nothing)
    front_end = sensor.front_end
    return BiOEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.propagation,
        front_end.pupil_samples, front_end.binning,
        source)
end

@inline function bi_o_edge_front_end_with_source(
    front_end::BiOEdgeOpticalFrontEnd, source)
    return BiOEdgeOpticalFrontEnd(front_end.amplitude_mask,
        front_end.modulation, front_end.propagation,
        front_end.pupil_samples, front_end.binning,
        source)
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
    propagation = bi_o_edge_propagation_workspace(wfs)
    T = eltype(propagation.intensity)
    n = size(propagation.bi_o_edge_masks, 1)
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
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_sub = wfs.front_end.pupil_samples
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
    return wfs
end
