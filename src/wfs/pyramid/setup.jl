#
# Pyramid wavefront sensing
#
# The diffractive pyramid model follows the standard optical sequence:
#
# 1. propagate the pupil field to the focal plane
# 2. apply the pyramid phase mask
# 3. propagate back to the re-imaged pupil plane
# 4. combine the four pupil images into differential slope signals
#
# Modulation is represented explicitly by averaging across a discrete set of
# focal-plane phase tilts. GPU/runtime optimizations keep the same optical model
# but batch modulation points and compatible asterism sources where possible.
#
@kernel function pyramid_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        phase = scale * (i + j - 2)
        @inbounds phasor[i, j] = cis(phase)
    end
end

@kernel function pyramid_mask_kernel!(mask, r, norma, start, step, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        x = start + (i - 1) * step
        y = start + (j - 1) * step
        p1 = x * r + y * r
        p2 = -x * r + y * r
        p3 = -x * r - y * r
        p4 = x * r - y * r
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        @inbounds mask[i, j] = cis(phase)
    end
end

@kernel function pyramid_pupil_modulation_batch_kernel!(stack, amplitude,
    opd, phases, weights, phasor, amplitude_scale, opd_to_cycles,
    offset::Int, resolution::Int, first_point::Int, pad::Int,
    batch_size::Int)
    i, j, batch_index = @index(Global, NTuple)
    if i <= pad && j <= pad && batch_index <= batch_size
        pupil_i = i - offset
        pupil_j = j - offset
        if 1 <= pupil_i <= resolution && 1 <= pupil_j <= resolution
            point = first_point + batch_index - 1
            @inbounds begin
                value = amplitude_scale * weights[point] *
                    amplitude[pupil_i, pupil_j] *
                    phases[pupil_i, pupil_j, point] *
                    cispi(opd_to_cycles * opd[pupil_i, pupil_j])
                stack[i, j, batch_index] = value * phasor[i, j]
            end
        else
            @inbounds stack[i, j, batch_index] = zero(eltype(stack))
        end
    end
end

@kernel function pyramid_electric_field_modulation_batch_kernel!(stack,
    field, phases, weights, phasor, offset::Int, resolution::Int,
    first_point::Int, pad::Int, batch_size::Int)
    i, j, batch_index = @index(Global, NTuple)
    if i <= pad && j <= pad && batch_index <= batch_size
        pupil_i = i - offset
        pupil_j = j - offset
        if 1 <= pupil_i <= resolution && 1 <= pupil_j <= resolution
            point = first_point + batch_index - 1
            @inbounds begin
                value = weights[point] * field[pupil_i, pupil_j] *
                    phases[pupil_i, pupil_j, point]
                stack[i, j, batch_index] = value * phasor[i, j]
            end
        else
            @inbounds stack[i, j, batch_index] = zero(eltype(stack))
        end
    end
end

@kernel function pyramid_modulation_batch_mask_kernel!(stack, mask,
    pad::Int, batch_size::Int)
    i, j, batch_index = @index(Global, NTuple)
    if i <= pad && j <= pad && batch_index <= batch_size
        @inbounds stack[i, j, batch_index] *= mask[i, j]
    end
end

@kernel function pyramid_modulation_batch_intensity_kernel!(out, stack,
    pad::Int, batch_size::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        value = @inbounds out[i, j]
        @inbounds for batch_index in 1:batch_size
            value += abs2(stack[i, j, batch_index])
        end
        @inbounds out[i, j] = value
    end
end

@kernel function pyramid_shifted_mask_stack_kernel!(masks,
    axis_1_shifts_rad, axis_2_shifts_rad, r, norma, rooftop_pixels,
    coordinate_start, coordinate_step, rotation_cos, rotation_sin,
    shift_x_1, shift_y_1, shift_x_2, shift_y_2,
    shift_x_3, shift_y_3, shift_x_4, shift_y_4,
    pad::Int, point_count::Int)
    axis_1, axis_2, point = @index(Global, NTuple)
    if axis_1 <= pad && axis_2 <= pad && point <= point_count
        unrotated_x = coordinate_start + (axis_1 - 1) * coordinate_step +
            @inbounds(axis_1_shifts_rad[point])
        unrotated_y = coordinate_start + (axis_2 - 1) * coordinate_step +
            @inbounds(axis_2_shifts_rad[point])
        x = unrotated_x * rotation_cos - unrotated_y * rotation_sin
        y = unrotated_y * rotation_cos + unrotated_x * rotation_sin
        phase_1 = x * r + unrotated_x * shift_x_1 + y * r -
            unrotated_y * shift_y_1 + rooftop_pixels
        phase_2 = -x * r + unrotated_x * shift_x_2 + y * r -
            unrotated_y * shift_y_2
        phase_3 = -x * r + unrotated_x * shift_x_3 - y * r -
            unrotated_y * shift_y_3 + rooftop_pixels
        phase_4 = x * r + unrotated_x * shift_x_4 - y * r -
            unrotated_y * shift_y_4
        phase = -max(max(phase_1, phase_2), max(phase_3, phase_4)) * norma
        @inbounds masks[axis_1, axis_2, point] = cis(phase)
    end
end

@kernel function pyramid_unmodulated_pupil_field_kernel!(focal_field,
    amplitude, opd, phasor, amplitude_scale, opd_to_cycles,
    offset::Int, resolution::Int, pad::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        pupil_i = i - offset
        pupil_j = j - offset
        if 1 <= pupil_i <= resolution && 1 <= pupil_j <= resolution
            @inbounds focal_field[i, j] = amplitude_scale *
                amplitude[pupil_i, pupil_j] *
                cispi(opd_to_cycles * opd[pupil_i, pupil_j]) * phasor[i, j]
        else
            @inbounds focal_field[i, j] = zero(eltype(focal_field))
        end
    end
end

@kernel function pyramid_unmodulated_electric_field_kernel!(focal_field,
    field, phasor, offset::Int, resolution::Int, pad::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        pupil_i = i - offset
        pupil_j = j - offset
        if 1 <= pupil_i <= resolution && 1 <= pupil_j <= resolution
            @inbounds focal_field[i, j] =
                field[pupil_i, pupil_j] * phasor[i, j]
        else
            @inbounds focal_field[i, j] = zero(eltype(focal_field))
        end
    end
end

@kernel function pyramid_shifted_mask_batch_kernel!(stack, focal_field,
    shifted_masks, weights, first_point::Int, pad::Int, batch_size::Int)
    i, j, batch_index = @index(Global, NTuple)
    if i <= pad && j <= pad && batch_index <= batch_size
        point = first_point + batch_index - 1
        @inbounds stack[i, j, batch_index] = weights[point] *
            focal_field[i, j] * shifted_masks[i, j, point]
    end
end

@kernel function pyramid_slopes_kernel!(slopes, intensity, valid_mask, sub::Int, n_sub::Int, pad::Int, offset::Int,
    ox1::Int, oy1::Int, ox2::Int, oy2::Int, ox3::Int, oy3::Int, ox4::Int, oy4::Int,
    sx1::Int, sy1::Int, sx2::Int, sy2::Int, sx3::Int, sy3::Int, sx4::Int, sy4::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if @inbounds valid_mask[i, j]
            q1 = zero(eltype(slopes))
            q2 = zero(eltype(slopes))
            q3 = zero(eltype(slopes))
            q4 = zero(eltype(slopes))
            @inbounds for di in 0:(sub - 1), dj in 0:(sub - 1)
                x = xs + di - 1
                y = ys + dj - 1
                x1 = ox1 + x + sx1
                y1 = oy1 + y + sy1
                x2 = ox2 + x + sx2
                y2 = oy2 + y + sy2
                x3 = ox3 + x + sx3
                y3 = oy3 + y + sy3
                x4 = ox4 + x + sx4
                y4 = oy4 + y + sy4
                if 1 <= x1 <= pad && 1 <= y1 <= pad
                    q1 += intensity[x1, y1]
                end
                if 1 <= x2 <= pad && 1 <= y2 <= pad
                    q2 += intensity[x2, y2]
                end
                if 1 <= x3 <= pad && 1 <= y3 <= pad
                    q3 += intensity[x3, y3]
                end
                if 1 <= x4 <= pad && 1 <= y4 <= pad
                    q4 += intensity[x4, y4]
                end
            end
            left = q1 + q3
            right = q2 + q4
            bottom = q1 + q2
            top = q3 + q4
            total = left + right
            if total <= 0
                slopes[idx] = zero(eltype(slopes))
                slopes[idx + offset] = zero(eltype(slopes))
            else
                slopes[idx] = (right - left) / total
                slopes[idx + offset] = (top - bottom) / total
            end
        else
            slopes[idx] = zero(eltype(slopes))
            slopes[idx + offset] = zero(eltype(slopes))
        end
    end
end

@kernel function gather_pyramid_slopes_kernel!(slopes, signal_2d, valid_signal_indices, count::Int, y_offset::Int)
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
struct PyramidEstimatorParams{T<:AbstractFloat,N<:WFSNormalization}
    pupil_samples::Int
    pupil_resolution::Int
    pupil_diameter_m::T
    threshold::T
    light_ratio::T
    normalization::N
    geometric_modulation_radius::T
end

"""Dispatch family for prepared Pyramid modulation propagation."""
abstract type AbstractPyramidModulationPropagationStrategy end

"""Reference modulation propagation using one pupil tilt per point."""
struct PyramidPupilTiltStrategy <:
       AbstractPyramidModulationPropagationStrategy end

"""Approximate fixed-modulation propagation using shifted focal masks."""
struct PyramidShiftedMaskStrategy <:
       AbstractPyramidModulationPropagationStrategy end

"""Run-immutable numerical contract for pyramid-mask propagation."""
struct PyramidPropagationPlan{M<:PyramidPhaseMask,T<:AbstractFloat,
    S<:AbstractPyramidModulationPropagationStrategy}
    phase_mask::M
    pupil_samples::Int
    binning::Int
    numeric_type::Type{T}
    modulation_propagation_strategy::S
end

"""Marker for scalar or uncentered propagation without modulation batching."""
struct NoPyramidModulationBatchWorkspace end

"""Accelerator scratch and plans for one bounded modulation-point tile."""
struct PyramidModulationBatchWorkspace{C,V,Pf,Pi}
    field_stack::C
    operating_weights::V
    calibration_weights::V
    fft_plan::Pf
    ifft_plan::Pi
    batch_size::Int
end

"""
Prepared shifted focal-plane masks and one bounded inverse-propagation tile.

The masks are derived cache, not persistent scientific state. They are valid
only for the fixed operating modulation used during preparation.
"""
struct PyramidShiftedMaskModulationWorkspace{C,M,V,Pi}
    field_stack::C
    shifted_masks::M
    operating_weights::V
    axis_1_shifts_rad::V
    axis_2_shifts_rad::V
    ifft_plan::Pi
    batch_size::Int
end

"""
Backend-bound FFT handles, caches, and replaceable single-writer scratch for
pyramid-mask propagation. No field is a caller-visible optical product.
"""
mutable struct PyramidPropagationWorkspace{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}},
    MB}
    field::C
    focal_field::C
    pupil_field::C
    pyramid_mask::C
    phasor::C
    intensity::R
    temp::R
    scratch::R
    asterism_stack::RS
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::Kf
    lgs_kernel_tag::UInt
    effective_resolution::Int
    asterism_capacity::Int
    revision::UInt
    modulation_batch::MB
end

"""Exact plan/workspace owner for one pyramid propagation execution."""
struct PreparedPyramidPropagation{
    P<:PyramidPropagationPlan,W<:PyramidPropagationWorkspace}
    plan::P
    workspace::W
end

const _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT = 8

function _pyramid_modulation_batch_size(point_count::Int)
    for candidate in min(
        point_count, _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT):-1:1
        point_count % candidate == 0 && return candidate
    end
    return 1
end

@inline _pyramid_shifted_mask_batch_size(
    ::ScalarCPUStyle, ::Int) = 1
@inline _pyramid_shifted_mask_batch_size(
    ::AcceleratorStyle, point_count::Int) =
    _pyramid_modulation_batch_size(point_count)

@inline function _prepare_pyramid_modulation_batch(
    ::ScalarCPUStyle, field, phase_mask, modulation,
    calibration_modulation, ::PyramidPupilTiltStrategy)
    return NoPyramidModulationBatchWorkspace()
end

function _prepare_pyramid_modulation_batch(
    ::AcceleratorStyle, field, phase_mask, modulation,
    calibration_modulation, ::PyramidPupilTiltStrategy)
    phase_mask.psf_centering ||
        return NoPyramidModulationBatchWorkspace()
    point_count = modulation_point_count(modulation)
    modulation_point_count(calibration_modulation) == point_count || throw(
        InvalidConfiguration(
            "Pyramid operating and calibration modulation counts must match",
        ),
    )
    batch_size = _pyramid_modulation_batch_size(point_count)
    pad = size(field, 1)
    field_stack = similar(field, eltype(field), pad, pad, batch_size)
    operating_weights = similar(field, real(eltype(field)), point_count)
    calibration_weights = similar(operating_weights)
    copyto!(operating_weights, modulation.amplitude_weights)
    copyto!(calibration_weights, calibration_modulation.amplitude_weights)
    fft_plan = plan_repeated_fft_backend!(field_stack, (1, 2))
    ifft_plan = plan_repeated_ifft_backend!(field_stack, (1, 2))
    return PyramidModulationBatchWorkspace(
        field_stack,
        operating_weights,
        calibration_weights,
        fft_plan,
        ifft_plan,
        batch_size,
    )
end


function _prepare_pyramid_modulation_batch(
    style::ExecutionStyle, field, phase_mask, modulation,
    calibration_modulation, ::PyramidShiftedMaskStrategy)
    phase_mask.psf_centering || throw(InvalidConfiguration(
        "Pyramid shifted-mask modulation requires psf_centering=true",
    ))
    phase_mask.old_mask && throw(InvalidConfiguration(
        "Pyramid shifted-mask modulation does not support old_mask=true",
    ))
    resolution = size(modulation.phases, 1)
    resolution > 1 || throw(InvalidConfiguration(
        "Pyramid shifted-mask modulation requires a pupil resolution greater than one",
    ))
    point_count = modulation_point_count(modulation)
    batch_size = _pyramid_shifted_mask_batch_size(style, point_count)
    pad = size(field, 1)
    field_stack = similar(field, eltype(field), pad, pad, batch_size)
    shifted_masks = similar(field, eltype(field), pad, pad, point_count)
    T = real(eltype(field))
    operating_weights = similar(field, T, point_count)
    axis_1_shifts_rad = similar(operating_weights)
    axis_2_shifts_rad = similar(operating_weights)
    host_axis_1_shifts_rad = Vector{T}(undef, point_count)
    host_axis_2_shifts_rad = Vector{T}(undef, point_count)
    coordinate_scale = T(2pi) / T(resolution - 1)
    @inbounds for point in 1:point_count
        offset_x, offset_y = modulation_offset(modulation.policy, point, T)
        host_axis_1_shifts_rad[point] = coordinate_scale * offset_y
        host_axis_2_shifts_rad[point] = coordinate_scale * offset_x
    end
    copyto!(operating_weights, modulation.amplitude_weights)
    copyto!(axis_1_shifts_rad, host_axis_1_shifts_rad)
    copyto!(axis_2_shifts_rad, host_axis_2_shifts_rad)
    ifft_plan = plan_repeated_ifft_backend!(field_stack, (1, 2))
    return PyramidShiftedMaskModulationWorkspace(
        field_stack,
        shifted_masks,
        operating_weights,
        axis_1_shifts_rad,
        axis_2_shifts_rad,
        ifft_plan,
        batch_size,
    )
end

@inline function _resize_pyramid_modulation_batch(
    batch::NoPyramidModulationBatchWorkspace, field)
    return batch
end

function _resize_pyramid_modulation_batch(
    batch::PyramidModulationBatchWorkspace, field)
    pad = size(field, 1)
    field_stack = similar(
        batch.field_stack,
        eltype(batch.field_stack),
        pad,
        pad,
        batch.batch_size,
    )
    fft_plan = plan_repeated_fft_backend!(field_stack, (1, 2))
    ifft_plan = plan_repeated_ifft_backend!(field_stack, (1, 2))
    return PyramidModulationBatchWorkspace(
        field_stack,
        batch.operating_weights,
        batch.calibration_weights,
        fft_plan,
        ifft_plan,
        batch.batch_size,
    )
end


function _resize_pyramid_modulation_batch(
    batch::PyramidShiftedMaskModulationWorkspace, field)
    pad = size(field, 1)
    point_count = size(batch.shifted_masks, 3)
    field_stack = similar(
        batch.field_stack,
        eltype(batch.field_stack),
        pad,
        pad,
        batch.batch_size,
    )
    shifted_masks = similar(
        batch.shifted_masks,
        eltype(batch.shifted_masks),
        pad,
        pad,
        point_count,
    )
    ifft_plan = plan_repeated_ifft_backend!(field_stack, (1, 2))
    return PyramidShiftedMaskModulationWorkspace(
        field_stack,
        shifted_masks,
        batch.operating_weights,
        batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad,
        ifft_plan,
        batch.batch_size,
    )
end

@inline pyramid_propagation_plan(
    propagation::PreparedPyramidPropagation) = propagation.plan
@inline pyramid_propagation_workspace(
    propagation::PreparedPyramidPropagation) = propagation.workspace

"""A physically distinct pyramid front end with prepared modulation."""
struct PyramidOpticalFrontEnd{O<:PyramidPhaseMask,M,C,P,S}
    phase_mask::O
    modulation::M
    calibration_modulation::C
    propagation::P
    pupil_samples::Int
    binning::Int
    source::S
end

"""Run-immutable family acquisition contract."""
struct PyramidAcquisitionPlan
    binning::Int
end

"""Derived native sampling metadata for convenience-frame acquisition."""
mutable struct PyramidAcquisitionWorkspace
    nominal_detector_resolution::Int
end

"""Caller-visible convenience-frame product."""
mutable struct PyramidAcquisitionProducts{T<:AbstractFloat,
    R<:AbstractMatrix{T}}
    frame::R
end

struct PyramidDetectorAcquisition{P,W,PR}
    plan::P
    workspace::W
    products::PR
end

"""Persistent support and calibration state for slope estimation."""
mutable struct PyramidEstimatorState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},V<:AbstractVector{T},
    R<:AbstractMatrix{T}}
    valid_mask::A
    optical_gain::V
    valid_i4q::A
    reference_signal_2d::R
    shift_x::NTuple{4,Int}
    shift_y::NTuple{4,Int}
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

"""Replaceable single-writer scratch for differential estimation."""
mutable struct PyramidEstimatorWorkspace{T<:AbstractFloat,
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
end

"""Caller-visible differential-slope product."""
mutable struct PyramidEstimatorProducts{T<:AbstractFloat,
    V<:AbstractVector{T}}
    slopes::V
end

struct PyramidDifferentialEstimator{P<:PyramidEstimatorParams,S,W,PR}
    params::P
    state::S
    workspace::W
    products::PR
end

struct PyramidWFS{M<:SensingMode,F,A,E,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::PyramidWFS{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

@inline pyramid_estimator_params(wfs::PyramidWFS) = wfs.estimator.params
@inline pyramid_estimator_state(wfs::PyramidWFS) = wfs.estimator.state
@inline pyramid_estimator_workspace(wfs::PyramidWFS) =
    wfs.estimator.workspace
@inline pyramid_estimator_products(wfs::PyramidWFS) =
    wfs.estimator.products
@inline pyramid_acquisition_plan(wfs::PyramidWFS{<:Diffractive}) =
    wfs.acquisition.plan
@inline pyramid_acquisition_workspace(wfs::PyramidWFS{<:Diffractive}) =
    wfs.acquisition.workspace
@inline pyramid_acquisition_products(wfs::PyramidWFS{<:Diffractive}) =
    wfs.acquisition.products
@inline pyramid_propagation(wfs::PyramidWFS{<:Diffractive}) =
    wfs.front_end.propagation
@inline pyramid_propagation_plan(wfs::PyramidWFS{<:Diffractive}) =
    pyramid_propagation_plan(pyramid_propagation(wfs))
@inline pyramid_propagation_workspace(wfs::PyramidWFS{<:Diffractive}) =
    pyramid_propagation_workspace(pyramid_propagation(wfs))
@inline pyramid_propagation_workspace(
    front_end::PyramidOpticalFrontEnd) =
    pyramid_propagation_workspace(front_end.propagation)
@inline four_pupil_propagation_workspace(
    front_end::PyramidOpticalFrontEnd) =
    pyramid_propagation_workspace(front_end)
@inline pyramid_phase_mask(wfs::PyramidWFS{<:Diffractive}) =
    wfs.front_end.phase_mask
@inline pyramid_operating_modulation(wfs::PyramidWFS{<:Diffractive}) =
    wfs.front_end.modulation
@inline pyramid_calibration_modulation(wfs::PyramidWFS{<:Diffractive}) =
    wfs.front_end.calibration_modulation

@inline function pyramid_sampled_geometry(pupil_samples::Int,
    n_pix_separation::Union{Int,Nothing}, n_pix_edge::Union{Int,Nothing},
    sampling::Int)
    sampling >= 1 || throw(InvalidConfiguration(
        "pyramid geometry sampling must be >= 1"))
    pupil_samples % sampling == 0 || throw(InvalidConfiguration(
        "pyramid sampling must preserve an integer pupil image"))
    n_pixels = div(pupil_samples, sampling)
    if n_pix_separation === nothing
        n_pix_edge === nothing || throw(InvalidConfiguration(
            "pyramid pupil-image edge padding requires an explicit separation"))
        return n_pixels, 0, 0
    end

    n_pix_separation >= 0 || throw(InvalidConfiguration(
        "pyramid pupil-image separation must be nonnegative"))
    edge = n_pix_edge === nothing ? div(n_pix_separation, 2) : n_pix_edge
    edge >= 0 || throw(InvalidConfiguration(
        "pyramid pupil-image edge padding must be nonnegative"))
    n_pix_separation % (2 * sampling) == 0 || throw(InvalidConfiguration(
        "pyramid pupil-image separation must remain an even integer after sampling"))
    edge % sampling == 0 || throw(InvalidConfiguration(
        "pyramid pupil-image edge padding must remain an integer after sampling"))
    return n_pixels, div(n_pix_separation, 2 * sampling), div(edge, sampling)
end

@inline function pyramid_native_frame_size(pupil_samples::Int,
    n_pix_separation::Int, n_pix_edge::Union{Int,Nothing})
    edge = n_pix_edge === nothing ? div(n_pix_separation, 2) : n_pix_edge
    return 2 * pupil_samples + n_pix_separation + 2 * edge
end

"""
    PyramidWFS(tel; ...)

Construct a pyramid wavefront sensor.

The diffractive model forms four re-imaged pupil intensities through a
focal-plane pyramid mask. Slopes are obtained from left/right and top/bottom
intensity differences after optional modulation averaging and binning.
`phase_mask_rotation_rad` sets the rotation angle applied to the physical
mask-coordinate transform, while
`modulation_phase_offset_rad` selects the circular modulation quadrature
origin. Both values are in radians. `modulation_propagation_strategy` selects
the exact pupil-tilt reference formulation or the approximate shifted-mask
formulation for fixed prepared modulation.
"""
function PyramidWFS(tel::Telescope; pupil_samples::Int, threshold::Real=0.1, modulation::Real=2.0,
    light_ratio::Real=0.0,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    modulation_points::Union{Int,Nothing}=nothing, extra_modulation_factor::Int=0,
    old_mask::Bool=false, rooftop::Real=0.0,
    phase_mask_rotation_rad::Real=0.0,
    modulation_phase_offset_rad::Real=0.0,
    user_modulation_path=nothing, mask_scale::Real=1.0, diffraction_padding::Int=2,
    psf_centering::Bool=true, n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
    modulation_propagation_strategy::AbstractPyramidModulationPropagationStrategy=
        PyramidPupilTiltStrategy(),
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
            "pyramid binning must evenly divide pupil_samples"))
    end
    typed_mask_scale = T(mask_scale)
    isfinite(typed_mask_scale) && typed_mask_scale > zero(T) ||
        throw(InvalidConfiguration(
            "pyramid mask_scale must be finite and > 0"))
    typed_phase_mask_rotation_rad = T(phase_mask_rotation_rad)
    isfinite(typed_phase_mask_rotation_rad) || throw(InvalidConfiguration(
        "pyramid phase_mask_rotation_rad must be finite"))
    typed_modulation_phase_offset_rad = T(modulation_phase_offset_rad)
    isfinite(typed_modulation_phase_offset_rad) || throw(
        InvalidConfiguration(
            "pyramid modulation_phase_offset_rad must be finite"))
    pyramid_sampled_geometry(pupil_samples, n_pix_separation, n_pix_edge,
        binning)
    estimator_params = PyramidEstimatorParams{T,typeof(normalization)}(
        pupil_samples,
        tel.params.resolution,
        T(tel.params.diameter),
        T(threshold),
        T(light_ratio),
        normalization,
        T(modulation))
    phase_mask = PyramidPhaseMask{T}(
        old_mask,
        T(rooftop),
        typed_phase_mask_rotation_rad,
        typed_mask_scale,
        diffraction_padding, psf_centering, n_pix_separation, n_pix_edge)
    operating_policy = legacy_modulation_policy(T(modulation),
        modulation_points, extra_modulation_factor,
        typed_modulation_phase_offset_rad,
        user_modulation_path)
    calibration_policy = calibration_modulation_policy(operating_policy,
        T(calib_modulation), typed_modulation_phase_offset_rad)
    valid_mask = backend{Bool}(undef, pupil_samples, pupil_samples)
    slopes = backend{T}(undef, 2 * pupil_samples * pupil_samples)
    fill!(slopes, zero(T))
    n_pix_signal = div(pupil_samples, binning)
    valid_i4q = backend{Bool}(undef, n_pix_signal, n_pix_signal)
    valid_i4q_host = Matrix{Bool}(undef, n_pix_signal, n_pix_signal)
    valid_signal = backend{Bool}(undef, 2 * n_pix_signal, n_pix_signal)
    valid_signal_indices = backend{Int}(undef, n_pix_signal * n_pix_signal)
    valid_signal_indices_host = Vector{Int}(undef, n_pix_signal * n_pix_signal)
    valid_flux_sum_buffer = backend{T}(undef, 1)
    valid_flux_sum_host = Vector{T}(undef, 1)
    valid_flux_i4q_host = Matrix{T}(undef, n_pix_signal, n_pix_signal)
    flux_i4q = backend{T}(undef, n_pix_signal, n_pix_signal)
    signal_2d = backend{T}(undef, 2 * n_pix_signal, n_pix_signal)
    reference_signal_2d = similar(signal_2d)
    fill!(reference_signal_2d, zero(T))
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    estimator_state = PyramidEstimatorState(valid_mask, optical_gain,
        valid_i4q, reference_signal_2d, (0, 0, 0, 0), (0, 0, 0, 0),
        false, zero(T), UInt(0), UInt(0))
    estimator_workspace = PyramidEstimatorWorkspace(valid_i4q_host,
        valid_signal, valid_signal_indices, valid_signal_indices_host, 0,
        valid_flux_sum_buffer, valid_flux_sum_host, valid_flux_i4q_host,
        flux_i4q, signal_2d)
    estimator_products = PyramidEstimatorProducts(slopes)
    estimator = PyramidDifferentialEstimator(estimator_params,
        estimator_state, estimator_workspace, estimator_products)
    front_end, acquisition = prepare_pyramid_mode(mode, backend, T, tel,
        phase_mask, operating_policy, calibration_policy, pupil_samples,
        binning, modulation_propagation_strategy)
    wfs = PyramidWFS{
        typeof(mode),typeof(front_end),typeof(acquisition),typeof(estimator),
        typeof(selector),
    }(front_end, acquisition, estimator)
    initialize_pyramid_valid_mask!(wfs, tel)
    prepare_pyramid_front_end!(mode, wfs, tel)
    return wfs
end

@inline prepare_pyramid_mode(::Geometric, backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning, ::PyramidPupilTiltStrategy) where {T} =
    (nothing, nothing)

function prepare_pyramid_mode(::Geometric, backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning, ::PyramidShiftedMaskStrategy) where {T}
    throw(InvalidConfiguration(
        "Pyramid shifted-mask modulation applies only to diffractive sensing",
    ))
end

function prepare_pyramid_mode(::Diffractive, backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning, modulation_propagation_strategy) where {
    T<:AbstractFloat,
}
    return _prepare_pyramid_diffractive_storage(backend, T, tel, phase_mask,
        operating_policy, calibration_policy, pupil_samples, binning,
        modulation_propagation_strategy)
end

function _prepare_pyramid_diffractive_storage(backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning, modulation_propagation_strategy) where {
    T<:AbstractFloat,
}
    pad = tel.params.resolution * phase_mask.diffraction_padding
    if phase_mask.n_pix_separation !== nothing
        edge = phase_mask.n_pix_edge === nothing ?
            div(phase_mask.n_pix_separation, 2) : phase_mask.n_pix_edge
        pad = Int(round((2 * pupil_samples + phase_mask.n_pix_separation +
            2 * edge) * tel.params.resolution / pupil_samples))
    end
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    mask = similar(field)
    phasor = similar(field)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    scratch = similar(intensity)
    asterism_stack = backend{T}(undef, pad, pad, 1)
    style = execution_style(field)
    fft_plan = plan_repeated_fft_backend!(focal_field)
    ifft_plan = plan_repeated_ifft_backend!(pupil_field)
    modulation = prepare_focal_plane_modulation(operating_policy,
        tel.params.resolution, field, T)
    calibration_modulation = prepare_focal_plane_modulation(
        calibration_policy, tel.params.resolution, field, T)
    modulation_batch = _prepare_pyramid_modulation_batch(
        style, field, phase_mask, modulation, calibration_modulation,
        modulation_propagation_strategy)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    propagation_plan = PyramidPropagationPlan(
        phase_mask, pupil_samples, binning, T,
        modulation_propagation_strategy)
    propagation_workspace = PyramidPropagationWorkspace(
        field, focal_field, pupil_field, mask, phasor, intensity, temp,
        scratch, asterism_stack, fft_plan, ifft_plan, elongation_kernel,
        lgs_kernel_fft, UInt(0), pad, 1, UInt(0), modulation_batch)
    propagation = PreparedPyramidPropagation(
        propagation_plan, propagation_workspace)
    front_end = PyramidOpticalFrontEnd(phase_mask, modulation,
        calibration_modulation, propagation, pupil_samples, binning, nothing)
    subaperture_pixels = div(tel.params.resolution, pupil_samples)
    nominal = div(pad, subaperture_pixels)
    camera_frame = backend{T}(undef, nominal, nominal)
    acquisition = PyramidDetectorAcquisition(PyramidAcquisitionPlan(binning),
        PyramidAcquisitionWorkspace(nominal),
        PyramidAcquisitionProducts(camera_frame))
    return front_end, acquisition
end

@inline prepare_pyramid_front_end!(::Geometric, ::PyramidWFS,
    ::Telescope) = nothing

function prepare_pyramid_front_end!(::Diffractive, wfs::PyramidWFS,
    tel::Telescope)
    build_pyramid_phasor!(pyramid_propagation_workspace(wfs).phasor)
    pupil = PupilFunction(tel)
    build_pyramid_mask!(wfs, pupil)
    build_pyramid_shifted_masks!(wfs, pupil)
    return nothing
end

function PyramidOpticalFrontEnd(sensor::PyramidWFS{<:Diffractive},
    source=nothing)
    front_end = sensor.front_end
    return PyramidOpticalFrontEnd(front_end.phase_mask, front_end.modulation,
        front_end.calibration_modulation, front_end.propagation,
        front_end.pupil_samples, front_end.binning, source)
end

@inline function pyramid_front_end_with_source(
    front_end::PyramidOpticalFrontEnd, source)
    return PyramidOpticalFrontEnd(front_end.phase_mask,
        front_end.modulation, front_end.calibration_modulation,
        front_end.propagation, front_end.pupil_samples, front_end.binning,
        source)
end

function PyramidOpticalFrontEnd(::PyramidWFS{<:Geometric}, source=nothing)
    throw(WFSPreparationError(:wfs_optics, :unsupported,
        "geometric pyramid sensing uses DirectMeasurementPath and has no optical front end"))
end

sensing_mode(::PyramidWFS{M}) where {M} = M()

function initialize_pyramid_valid_mask!(wfs::PyramidWFS,
    tel::Telescope)
    set_valid_subapertures!(wfs.estimator.state.valid_mask,
        pupil_mask(tel), wfs.estimator.params.threshold)
    return wfs
end

function update_valid_mask!(wfs::PyramidWFS, pupil::PupilFunction)
    set_valid_subapertures!(wfs.estimator.state.valid_mask,
        pupil.support, wfs.estimator.params.threshold)
    return wfs
end

function ensure_pyramid_buffers!(wfs::PyramidWFS, pad::Int, pupil::PupilFunction)
    propagation = pyramid_propagation_workspace(wfs)
    acquisition = pyramid_acquisition_products(wfs)
    if size(propagation.field) != (pad, pad)
        propagation.revision += UInt(1)
        propagation.field = similar(propagation.field, pad, pad)
        propagation.focal_field = similar(propagation.focal_field, pad, pad)
        propagation.pupil_field = similar(propagation.pupil_field, pad, pad)
        propagation.pyramid_mask = similar(propagation.pyramid_mask, pad, pad)
        propagation.phasor = similar(propagation.phasor, pad, pad)
        propagation.intensity = similar(propagation.intensity, pad, pad)
        propagation.temp = similar(propagation.temp, pad, pad)
        propagation.scratch = similar(propagation.scratch, pad, pad)
        acquisition.frame = similar(acquisition.frame, pad, pad)
        propagation.asterism_stack = similar(propagation.asterism_stack,
            pad, pad, propagation.asterism_capacity)
        propagation.fft_plan = plan_repeated_fft_backend!(
            propagation.focal_field)
        propagation.ifft_plan = plan_repeated_ifft_backend!(
            propagation.pupil_field)
        propagation.modulation_batch = _resize_pyramid_modulation_batch(
            propagation.modulation_batch, propagation.field)
        propagation.lgs_kernel_fft = similar(propagation.focal_field,
            eltype(propagation.focal_field), 0, 0)
        propagation.lgs_kernel_tag = UInt(0)
        propagation.effective_resolution = pad
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
        build_pyramid_phasor!(propagation.phasor)
        build_pyramid_mask!(wfs, pupil)
        build_pyramid_shifted_masks!(wfs, pupil)
    end
    return wfs
end

function ensure_pyramid_asterism_stack!(wfs::PyramidWFS, n_src::Int)
    n_src >= 1 || throw(InvalidConfiguration("asterism source count must be >= 1"))
    propagation = pyramid_propagation_workspace(wfs)
    pad = size(propagation.intensity, 1)
    if size(propagation.asterism_stack, 1) != pad ||
            size(propagation.asterism_stack, 2) != pad ||
            size(propagation.asterism_stack, 3) < n_src
        capacity = max(n_src, propagation.asterism_capacity)
        propagation.asterism_stack = similar(propagation.asterism_stack,
            pad, pad, capacity)
        propagation.asterism_capacity = capacity
    end
    return propagation.asterism_stack
end

@inline grouped_staging_buffer(wfs::PyramidWFS, out::AbstractMatrix) = pyramid_propagation_workspace(wfs).intensity

function accumulate_pyramid_asterism_intensity!(::ScalarCPUStyle, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    intensity = pyramid_propagation_workspace(wfs).intensity
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, intensity,
        stack, ast.sources, pyramid_intensity!, wfs, pupil)
end

function accumulate_pyramid_asterism_intensity!(style::AcceleratorStyle, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    intensity = pyramid_propagation_workspace(wfs).intensity
    return accumulate_grouped_sources!(style, wfs, intensity, stack,
        ast.sources, pyramid_intensity!, wfs, pupil)
end

@inline pyramid_spectral_component_qe(::Nothing, sample,
    ::Type{T}) where {T<:AbstractFloat} = one(T)

@inline pyramid_spectral_component_qe(model::AbstractQuantumEfficiencyModel,
    sample, ::Type{T}) where {T<:AbstractFloat} =
    T(qe_at(model, sample.wavelength))

function accumulate_pyramid_spectral_intensity!(style::ExecutionStyle,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel},
    modulation=pyramid_operating_modulation(wfs))
    count = length(src.bundle.samples)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    total_irradiance = photon_irradiance(src)
    T = eltype(pyramid_propagation_workspace(wfs).intensity)
    @inbounds for (sample_idx, sample) in pairs(src.bundle.samples)
        channel_qe = pyramid_spectral_component_qe(qe_model, sample, T)
        variant = source_with_wavelength_and_radiometric_value(src, sample.wavelength,
            T(total_irradiance * sample.weight * channel_qe))
        pyramid_intensity_core!(@view(stack[:, :, sample_idx]), wfs, pupil,
            variant, modulation)
    end
    return reduce_grouped_stack!(style, pyramid_propagation_workspace(wfs).intensity, stack, count)
end

accumulate_pyramid_spectral_intensity!(style::ExecutionStyle,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource) =
    accumulate_pyramid_spectral_intensity!(style, wfs, pupil, src, nothing)

@inline function pyramid_support_selection_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource, ::Nothing)
    return pyramid_intensity_core!(out, wfs, pupil, src,
        pyramid_calibration_modulation(wfs))
end

@inline function pyramid_support_selection_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource, ::Nothing)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil,
        src, nothing, pyramid_calibration_modulation(wfs))
    out === pyramid_propagation_workspace(wfs).intensity || copyto!(out, pyramid_propagation_workspace(wfs).intensity)
    return out
end

@inline function pyramid_support_selection_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil, src,
        qe_model, pyramid_calibration_modulation(wfs))
    out === pyramid_propagation_workspace(wfs).intensity || copyto!(out, pyramid_propagation_workspace(wfs).intensity)
    return out
end

# The legacy `calib_modulation` path broadens illumination only while selecting
# valid support. A zero-aberration reference must use the operating modulation
# so that subtracting it removes the actual sensor's static response.
@inline function pyramid_calibration_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource, ::Nothing)
    return pyramid_intensity_core!(out, wfs, pupil, src,
        pyramid_operating_modulation(wfs))
end

@inline function pyramid_calibration_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource, ::Nothing)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil,
        src, nothing, pyramid_operating_modulation(wfs))
    out === pyramid_propagation_workspace(wfs).intensity ||
        copyto!(out, pyramid_propagation_workspace(wfs).intensity)
    return out
end

@inline function pyramid_calibration_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil,
        src, qe_model, pyramid_operating_modulation(wfs))
    out === pyramid_propagation_workspace(wfs).intensity ||
        copyto!(out, pyramid_propagation_workspace(wfs).intensity)
    return out
end

function accumulate_pyramid_extended_intensity!(::ScalarCPUStyle, out::AbstractMatrix, wfs::PyramidWFS,
    pupil::PupilFunction, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return pyramid_intensity!(out, wfs, pupil, ast.sources[1])
    end
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, out, stack, ast.sources, pyramid_intensity!, wfs, pupil)
end

function accumulate_pyramid_extended_intensity!(style::AcceleratorStyle, out::AbstractMatrix, wfs::PyramidWFS,
    pupil::PupilFunction, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return pyramid_intensity!(out, wfs, pupil, ast.sources[1])
    end
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(style, wfs, out, stack, ast.sources, pyramid_intensity!, wfs, pupil)
end

function prepare_pyramid_sampling!(wfs::PyramidWFS, pupil::PupilFunction)
    n_sub = wfs.estimator.params.pupil_samples
    pad = _pupil_resolution(pupil) * wfs.front_end.phase_mask.diffraction_padding
    if wfs.front_end.phase_mask.n_pix_separation !== nothing
        pixels_per_pupil_sample = div(_pupil_resolution(pupil), n_sub)
        pad = pyramid_native_frame_size(n_sub,
            wfs.front_end.phase_mask.n_pix_separation, wfs.front_end.phase_mask.n_pix_edge) *
              pixels_per_pupil_sample
    end
    if pad < _pupil_resolution(pupil)
        throw(InvalidConfiguration("pyramid padding must be >= telescope resolution"))
    end
    if pad % pyramid_acquisition_plan(wfs).binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide padded resolution"))
    end
    if _pupil_resolution(pupil) % pyramid_acquisition_plan(wfs).binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide telescope resolution"))
    end
    ensure_pyramid_buffers!(wfs, pad, pupil)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration("sample_pyramid_intensity! requires telescope context"))
end

function resize_pyramid_signal_buffers!(wfs::PyramidWFS, frame_size::Int)
    acquisition_plan = pyramid_acquisition_plan(wfs)
    acquisition_workspace = pyramid_acquisition_workspace(wfs)
    acquisition_products = pyramid_acquisition_products(wfs)
    state = pyramid_estimator_state(wfs)
    workspace = pyramid_estimator_workspace(wfs)
    nominal = acquisition_workspace.nominal_detector_resolution
    frame_size >= 1 || throw(InvalidConfiguration(
        "pyramid camera frame size must be >= 1"))
    nominal > 0 || throw(InvalidConfiguration(
        "pyramid nominal detector resolution must be prepared before signal resizing"))
    nominal % acquisition_plan.binning == 0 || throw(InvalidConfiguration(
        "pyramid binning must evenly divide the nominal detector resolution"))
    sampled_size = div(nominal, acquisition_plan.binning)
    sampled_size % frame_size == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must evenly divide the pyramid camera frame"))
    reduction = div(sampled_size, frame_size)
    total_sampling = acquisition_plan.binning * reduction
    n_pixels, half_separation, edge_padding = pyramid_sampled_geometry(
        wfs.estimator.params.pupil_samples, wfs.front_end.phase_mask.n_pix_separation,
        wfs.front_end.phase_mask.n_pix_edge, total_sampling)
    n_pixels >= 1 || throw(InvalidConfiguration(
        "detector sampling and binning removed every pyramid pupil sample"))
    iseven(frame_size) || throw(InvalidConfiguration(
        "pyramid camera frame must have even dimensions for symmetric pupil extraction"))
    if wfs.front_end.phase_mask.n_pix_separation === nothing
        frame_size >= 2 * n_pixels || throw(InvalidConfiguration(
            "pyramid camera frame does not contain four complete pupil images"))
    else
        frame_size == 2 * (n_pixels + half_separation + edge_padding) ||
            throw(InvalidConfiguration(
                "pyramid camera frame does not exactly preserve the configured pupil-image geometry"))
    end
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
    if size(acquisition_products.frame) != (frame_size, frame_size)
        acquisition_products.frame = similar(acquisition_products.frame,
            frame_size, frame_size)
    end
    update_pyramid_valid_signal!(wfs)
    if calibration_storage_changed
        state.calibrated = false
        state.calibration_revision += UInt(1)
    end
    return wfs
end

@inline function require_pyramid_frame_geometry(wfs::PyramidWFS,
    frame::AbstractMatrix)
    n_rows, n_cols = size(frame)
    n_rows == n_cols || throw(DimensionMismatchError(
        "pyramid camera frame must be square"))
    n_rows >= 1 || throw(DimensionMismatchError(
        "pyramid camera frame must be nonempty"))
    iseven(n_rows) || throw(InvalidConfiguration(
        "pyramid camera frame must have even dimensions for symmetric pupil extraction"))
    n_pixels = size(pyramid_estimator_workspace(wfs).signal_2d, 2)
    nominal = pyramid_acquisition_workspace(wfs).nominal_detector_resolution
    nominal > 0 || throw(InvalidConfiguration(
        "pyramid nominal detector resolution must be prepared before signal extraction"))
    nominal % pyramid_acquisition_plan(wfs).binning == 0 || throw(InvalidConfiguration(
        "pyramid binning must evenly divide the nominal detector resolution"))
    sampled_size = div(nominal, pyramid_acquisition_plan(wfs).binning)
    sampled_size % n_rows == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must evenly divide the pyramid camera frame"))
    reduction = div(sampled_size, n_rows)
    total_sampling = pyramid_acquisition_plan(wfs).binning * reduction
    geometry_pixels, half_separation, edge_padding = pyramid_sampled_geometry(
        wfs.estimator.params.pupil_samples, wfs.front_end.phase_mask.n_pix_separation,
        wfs.front_end.phase_mask.n_pix_edge, total_sampling)
    geometry_pixels == n_pixels || throw(InvalidConfiguration(
        "pyramid signal buffers do not match the sampled pupil-image geometry"))
    center = div(n_rows, 2)
    if wfs.front_end.phase_mask.n_pix_separation === nothing
        center >= n_pixels || throw(DimensionMismatchError(
            "pyramid camera frame does not contain four complete pupil images"))
    else
        n_rows == 2 * (n_pixels + half_separation + edge_padding) ||
            throw(DimensionMismatchError(
                "pyramid camera frame does not exactly preserve the configured pupil-image geometry"))
    end
    return center, half_separation
end

function update_pyramid_valid_signal!(wfs::PyramidWFS)
    state = pyramid_estimator_state(wfs)
    workspace = pyramid_estimator_workspace(wfs)
    n_pixels = size(state.valid_i4q, 1)
    fill!(workspace.valid_signal, false)
    @views begin
        workspace.valid_signal[1:n_pixels, :] .= state.valid_i4q
        workspace.valid_signal[n_pixels+1:end, :] .= state.valid_i4q
    end
    return wfs
end

function update_pyramid_valid_signal_indices!(wfs::PyramidWFS)
    state = pyramid_estimator_state(wfs)
    workspace = pyramid_estimator_workspace(wfs)
    valid_host = workspace.valid_i4q_host
    if size(valid_host) != size(state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(state.valid_i4q)...)
        workspace.valid_i4q_host = valid_host
    end
    copyto!(valid_host, state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(workspace.valid_signal_indices) < n_valid
        workspace.valid_signal_indices = similar(
            workspace.valid_signal_indices, n_valid)
    end
    if length(workspace.valid_signal_indices_host) < n_valid
        workspace.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = workspace.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(workspace.valid_signal_indices, 1, host_indices, 1, n_valid)
    workspace.valid_signal_count = n_valid
    return n_valid
end

function resize_pyramid_slope_buffers!(wfs::PyramidWFS)
    workspace = pyramid_estimator_workspace(wfs)
    products = pyramid_estimator_products(wfs)
    state = pyramid_estimator_state(wfs)
    n_valid = workspace.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("pyramid valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(products.slopes) != n_slopes
        products.slopes = similar(products.slopes, n_slopes)
    end
    if length(state.optical_gain) != n_slopes
        state.optical_gain = similar(state.optical_gain, n_slopes)
        fill!(state.optical_gain, one(eltype(state.optical_gain)))
    end
    return wfs
end

function pyramid_valid_flux_sum!(::ScalarCPUStyle, wfs::PyramidWFS,
    i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    workspace = pyramid_estimator_workspace(wfs)
    return masked_sum2d(ScalarCPUStyle(), i4q, workspace.valid_i4q_host)
end

function pyramid_valid_flux_sum!(style::AcceleratorStyle, wfs::PyramidWFS,
    i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    workspace = pyramid_estimator_workspace(wfs)
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.estimator.state.valid_i4q,
        workspace.valid_i4q_host,
        workspace.valid_flux_sum_buffer,
        workspace.valid_flux_sum_host,
        workspace.valid_flux_i4q_host,
    )
    workspace.valid_flux_i4q_host = host_parent
    return summed
end

function select_pyramid_valid_i4q_from_frame!(::ScalarCPUStyle,
    wfs::PyramidWFS, frame::AbstractMatrix)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    max_i4q = zero(eltype(frame))
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i,
            center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i,
            center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i,
            center - n_extra - n_pixels + j]
        max_i4q = max(max_i4q, q1 + q2 + q3 + q4)
    end
    cutoff = wfs.estimator.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i,
            center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i,
            center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i,
            center - n_extra - n_pixels + j]
        wfs.estimator.state.valid_i4q[i, j] =
            (q1 + q2 + q3 + q4) >= cutoff
    end
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function select_pyramid_valid_i4q_from_frame!(::AcceleratorStyle,
    wfs::PyramidWFS, frame::AbstractMatrix)
    workspace = pyramid_estimator_workspace(wfs)
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    rows_lo = center - n_extra - n_pixels + 1:center - n_extra
    rows_hi = center + n_extra + 1:center + n_extra + n_pixels
    cols_lo = center - n_extra - n_pixels + 1:center - n_extra
    cols_hi = center + n_extra + 1:center + n_extra + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view workspace.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. wfs.estimator.state.valid_i4q = i4q >= cutoff
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function select_pyramid_valid_i4q!(wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel}, det::Detector)
    propagation = pyramid_propagation_workspace(wfs)
    pyramid_support_selection_intensity!(propagation.temp, wfs, pupil, src,
        qe_model)
    sampled = sample_pyramid_intensity!(wfs, pupil, propagation.temp)
    frame = detector_calibration_frame!(det, sampled,
        pyramid_detector_calibration_qe(src, det, eltype(det.products.frame)))
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    return select_pyramid_valid_i4q_from_frame!(execution_style(frame), wfs,
        frame)
end

function select_pyramid_valid_i4q!(wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource)
    return select_pyramid_valid_i4q!(wfs, pupil, src, nothing)
end

function select_pyramid_valid_i4q!(wfs::PyramidWFS, pupil::PupilFunction,
    src::AbstractSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel})
    return select_pyramid_valid_i4q!(execution_style(wfs.estimator.state.valid_i4q),
        wfs, pupil, src, qe_model)
end

function select_pyramid_valid_i4q!(style::ScalarCPUStyle,
    wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    return select_pyramid_valid_i4q!(style, wfs, pupil, src, nothing)
end

function select_pyramid_valid_i4q!(::ScalarCPUStyle, wfs::PyramidWFS,
    pupil::PupilFunction, src::AbstractSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel})
    state = pyramid_estimator_state(wfs)
    propagation = pyramid_propagation_workspace(wfs)
    binning = pyramid_acquisition_plan(wfs).binning
    n_pixels = div(wfs.estimator.params.pupil_samples, binning)
    if size(state.valid_i4q) != (n_pixels, n_pixels)
        state.valid_i4q = similar(state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    pyramid_support_selection_intensity!(propagation.temp, wfs, pupil, src,
        qe_model)
    frame = sample_pyramid_intensity!(wfs, pupil, propagation.temp)

    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    max_i4q = zero(eltype(frame))
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i, center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i, center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i, center - n_extra - n_pixels + j]
        i4q = q1 + q2 + q3 + q4
        if i4q > max_i4q
            max_i4q = i4q
        end
    end
    cutoff = wfs.estimator.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i, center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i, center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i, center - n_extra - n_pixels + j]
        state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
    end
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function select_pyramid_valid_i4q!(style::AcceleratorStyle,
    wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource)
    return select_pyramid_valid_i4q!(style, wfs, pupil, src, nothing)
end

function select_pyramid_valid_i4q!(::AcceleratorStyle,
    wfs::PyramidWFS, pupil::PupilFunction, src::AbstractSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel})
    state = pyramid_estimator_state(wfs)
    workspace = pyramid_estimator_workspace(wfs)
    propagation = pyramid_propagation_workspace(wfs)
    binning = pyramid_acquisition_plan(wfs).binning
    n_pixels = div(wfs.estimator.params.pupil_samples, binning)
    if size(state.valid_i4q) != (n_pixels, n_pixels)
        state.valid_i4q = similar(state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    pyramid_support_selection_intensity!(propagation.temp, wfs, pupil, src,
        qe_model)
    frame = sample_pyramid_intensity!(wfs, pupil, propagation.temp)

    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    rows_lo = center - n_extra - n_pixels + 1:center - n_extra
    rows_hi = center + n_extra + 1:center + n_extra + n_pixels
    cols_lo = center - n_extra - n_pixels + 1:center - n_extra
    cols_hi = center + n_extra + 1:center + n_extra + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view workspace.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. state.valid_i4q = i4q >= cutoff
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, pupil::PupilFunction, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    acquisition_workspace = pyramid_acquisition_workspace(wfs)
    acquisition_products = pyramid_acquisition_products(wfs)
    binning = pyramid_acquisition_plan(wfs).binning
    sub = div(_pupil_resolution(pupil), wfs.estimator.params.pupil_samples)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("pyramid intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    acquisition_workspace.nominal_detector_resolution = n_camera
    frame = acquisition_products.frame
    if binning != 1
        if n_camera % binning != 0
            throw(InvalidConfiguration("pyramid binning must evenly divide detector resolution"))
        end
        n_binned = div(n_camera, binning)
        if size(frame) != (n_binned, n_binned)
            acquisition_products.frame = similar(frame, n_binned, n_binned)
            frame = acquisition_products.frame
        end
        bin2d!(frame, intensity, sub * binning)
    else
        if size(frame) != (n_camera, n_camera)
            acquisition_products.frame = similar(frame, n_camera, n_camera)
            frame = acquisition_products.frame
        end
        bin2d!(frame, intensity, sub)
    end
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    return frame
end
