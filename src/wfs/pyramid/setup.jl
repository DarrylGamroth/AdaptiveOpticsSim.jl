#
# Pyramid wavefront sensing
#
# The diffractive pyramid model follows the standard optical sequence:
#
# 1. propagate the pupil field to the focal plane
# 2. apply the pyramid phase mask
# 3. propagate back to the re-imaged pupil plane
# 4. publish the four-pupil detector-plane photon-rate map
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
    intensity_scale, pad::Int, batch_size::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        value = @inbounds out[i, j]
        batch_intensity = zero(eltype(out))
        @inbounds for batch_index in 1:batch_size
            batch_intensity += abs2(stack[i, j, batch_index])
        end
        @inbounds out[i, j] = value + intensity_scale * batch_intensity
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

# One work item owns every prepared modulation plane at a focal coordinate.
# The 2-D launch therefore loads the shared focal field once per coordinate;
# writes remain coalesced because all work items advance through the same plane.
@kernel function pyramid_shifted_mask_batch_kernel!(stack, focal_field,
    shifted_masks, weights, first_point::Int, pad::Int, batch_size::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        focal_value = @inbounds focal_field[i, j]
        @inbounds for batch_index in 1:batch_size
            point = first_point + batch_index - 1
            stack[i, j, batch_index] = weights[point] *
                focal_value * shifted_masks[i, j, point]
        end
    end
end

@kernel function pyramid_separable_shifted_mask_factors_kernel!(
    axis_1_factors, axis_2_factors, axis_1_shifts_rad, axis_2_shifts_rad,
    r, norma, coordinate_start, coordinate_step,
    shift_x_positive, shift_x_negative,
    shift_y_positive, shift_y_negative,
    pad::Int, point_count::Int)
    axis, point = @index(Global, NTuple)
    if axis <= pad && point <= point_count
        coordinate = coordinate_start + (axis - 1) * coordinate_step
        axis_1_coordinate = coordinate + @inbounds(axis_1_shifts_rad[point])
        axis_2_coordinate = coordinate + @inbounds(axis_2_shifts_rad[point])
        axis_1_phase = -max(
            axis_1_coordinate * (r + shift_x_positive),
            axis_1_coordinate * (-r + shift_x_negative),
        ) * norma
        axis_2_phase = -max(
            axis_2_coordinate * (r - shift_y_positive),
            axis_2_coordinate * (-r - shift_y_negative),
        ) * norma
        @inbounds begin
            axis_1_factors[axis, point] = cis(axis_1_phase)
            axis_2_factors[axis, point] = cis(axis_2_phase)
        end
    end
end

@kernel function pyramid_separable_shifted_mask_batch_kernel!(stack,
    focal_field, axis_1_factors, axis_2_factors, weights,
    first_point::Int, pad::Int, batch_size::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        focal_value = @inbounds focal_field[i, j]
        @inbounds for batch_index in 1:batch_size
            point = first_point + batch_index - 1
            stack[i, j, batch_index] = weights[point] *
                focal_value * axis_1_factors[i, point] *
                axis_2_factors[j, point]
        end
    end
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
struct PyramidModulationBatchWorkspace{C,V,Pf,Pb}
    field_stack::C
    operating_weights::V
    fft_plan::Pf
    bfft_plan::Pb
    batch_size::Int
end

"""Dispatch family for prepared shifted-mask modulation storage."""
abstract type AbstractPyramidShiftedMaskModulationWorkspace end

"""
Prepared full shifted focal-plane masks and one bounded inverse-propagation tile.

The masks are derived cache, not persistent scientific state. They are valid
only for the fixed operating modulation used during preparation.
"""
struct PyramidShiftedMaskModulationWorkspace{C,M,V,Pb} <:
       AbstractPyramidShiftedMaskModulationWorkspace
    field_stack::C
    shifted_masks::M
    operating_weights::V
    axis_1_shifts_rad::V
    axis_2_shifts_rad::V
    bfft_plan::Pb
    batch_size::Int
end

"""
Prepared separable factors for an ideal shifted focal-plane Pyramid mask.

This representation is exact for the shifted-mask formulation when physical
mask rotation and rooftop terms are both zero. It replaces each `pad × pad`
mask with two `pad`-element factors while retaining the same bounded inverse
propagation tile.
"""
struct PyramidSeparableShiftedMaskModulationWorkspace{C,M,V,Pb} <:
       AbstractPyramidShiftedMaskModulationWorkspace
    field_stack::C
    axis_1_factors::M
    axis_2_factors::M
    operating_weights::V
    axis_1_shifts_rad::V
    axis_2_shifts_rad::V
    bfft_plan::Pb
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

# Larger batches amortize FFT submission and intermediate kernels, while this
# scratch budget prevents the prepared field stack from scaling without bound
# on larger propagation grids. The largest divisor within both limits is fixed
# during preparation; the repeated path performs no device query or tuning.
const _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT = 32
const _PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES = 256 * 1024 * 1024

function _pyramid_modulation_batch_size(
    point_count::Int,
    plane_bytes::Int,
)
    point_count >= 1 || throw(InvalidConfiguration(
        "Pyramid modulation requires at least one point",
    ))
    plane_bytes >= 1 || throw(InvalidConfiguration(
        "Pyramid modulation field planes must occupy at least one byte",
    ))
    memory_limit = max(1, div(
        _PYRAMID_ACCELERATOR_MODULATION_BATCH_BUDGET_BYTES,
        plane_bytes,
    ))
    candidate_limit = min(
        point_count,
        _PYRAMID_ACCELERATOR_MODULATION_BATCH_LIMIT,
        memory_limit,
    )
    for candidate in candidate_limit:-1:1
        point_count % candidate == 0 && return candidate
    end
    return 1
end

@inline function _pyramid_modulation_batch_size(field, point_count::Int)
    plane_bytes = Base.checked_mul(sizeof(eltype(field)), length(field))
    return _pyramid_modulation_batch_size(point_count, plane_bytes)
end

@inline _pyramid_shifted_mask_batch_size(
    ::ScalarCPUStyle, field, ::Int) = 1
@inline _pyramid_shifted_mask_batch_size(
    ::AcceleratorStyle, field, point_count::Int) =
    _pyramid_modulation_batch_size(field, point_count)

@inline _pyramid_shifted_mask_is_separable(phase_mask) =
    iszero(phase_mask.rotation_rad) && iszero(phase_mask.rooftop)

@inline function _prepare_pyramid_modulation_batch(
    ::ScalarCPUStyle, field, phase_mask, modulation,
    ::PyramidPupilTiltStrategy)
    return NoPyramidModulationBatchWorkspace()
end

function _prepare_pyramid_modulation_batch(
    ::AcceleratorStyle, field, phase_mask, modulation,
    ::PyramidPupilTiltStrategy)
    phase_mask.psf_centering ||
        return NoPyramidModulationBatchWorkspace()
    point_count = modulation_point_count(modulation)
    pad = size(field, 1)
    batch_size = _pyramid_modulation_batch_size(field, point_count)
    field_stack = similar(field, eltype(field), pad, pad, batch_size)
    operating_weights = similar(field, real(eltype(field)), point_count)
    copyto!(operating_weights, modulation.amplitude_weights)
    fft_plan = plan_repeated_fft_backend!(field_stack, (1, 2))
    bfft_plan = plan_repeated_bfft_backend!(field_stack, (1, 2))
    return PyramidModulationBatchWorkspace(
        field_stack,
        operating_weights,
        fft_plan,
        bfft_plan,
        batch_size,
    )
end


function _prepare_pyramid_modulation_batch(
    style::ExecutionStyle, field, phase_mask, modulation,
    ::PyramidShiftedMaskStrategy)
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
    batch_size = _pyramid_shifted_mask_batch_size(style, field, point_count)
    pad = size(field, 1)
    field_stack = similar(field, eltype(field), pad, pad, batch_size)
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
    bfft_plan = plan_repeated_bfft_backend!(field_stack, (1, 2))
    if _pyramid_shifted_mask_is_separable(phase_mask)
        axis_1_factors = similar(
            field, eltype(field), pad, point_count)
        axis_2_factors = similar(axis_1_factors)
        return PyramidSeparableShiftedMaskModulationWorkspace(
            field_stack,
            axis_1_factors,
            axis_2_factors,
            operating_weights,
            axis_1_shifts_rad,
            axis_2_shifts_rad,
            bfft_plan,
            batch_size,
        )
    end
    shifted_masks = similar(field, eltype(field), pad, pad, point_count)
    return PyramidShiftedMaskModulationWorkspace(
        field_stack,
        shifted_masks,
        operating_weights,
        axis_1_shifts_rad,
        axis_2_shifts_rad,
        bfft_plan,
        batch_size,
    )
end

function _resize_pyramid_modulation_batch(
    batch::PyramidSeparableShiftedMaskModulationWorkspace, field)
    pad = size(field, 1)
    point_count = size(batch.axis_1_factors, 2)
    field_stack = similar(
        batch.field_stack,
        eltype(batch.field_stack),
        pad,
        pad,
        batch.batch_size,
    )
    axis_1_factors = similar(
        batch.axis_1_factors,
        eltype(batch.axis_1_factors),
        pad,
        point_count,
    )
    axis_2_factors = similar(axis_1_factors)
    bfft_plan = plan_repeated_bfft_backend!(field_stack, (1, 2))
    return PyramidSeparableShiftedMaskModulationWorkspace(
        field_stack,
        axis_1_factors,
        axis_2_factors,
        batch.operating_weights,
        batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad,
        bfft_plan,
        batch.batch_size,
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
    bfft_plan = plan_repeated_bfft_backend!(field_stack, (1, 2))
    return PyramidModulationBatchWorkspace(
        field_stack,
        batch.operating_weights,
        fft_plan,
        bfft_plan,
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
    bfft_plan = plan_repeated_bfft_backend!(field_stack, (1, 2))
    return PyramidShiftedMaskModulationWorkspace(
        field_stack,
        shifted_masks,
        batch.operating_weights,
        batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad,
        bfft_plan,
        batch.batch_size,
    )
end

@inline pyramid_propagation_plan(
    propagation::PreparedPyramidPropagation) = propagation.plan
@inline pyramid_propagation_workspace(
    propagation::PreparedPyramidPropagation) = propagation.workspace

"""A physically distinct pyramid front end with prepared modulation."""
struct PyramidOpticalFrontEnd{O<:PyramidPhaseMask,M,P,S}
    phase_mask::O
    modulation::M
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

struct PyramidWFS{F,A,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
end

@inline backend(::PyramidWFS{<:Any,<:Any,B}) where {B} = B()

@inline pyramid_acquisition_plan(wfs::PyramidWFS) =
    wfs.acquisition.plan
@inline pyramid_acquisition_workspace(wfs::PyramidWFS) =
    wfs.acquisition.workspace
@inline pyramid_acquisition_products(wfs::PyramidWFS) =
    wfs.acquisition.products
@inline pyramid_propagation(wfs::PyramidWFS) =
    wfs.front_end.propagation
@inline pyramid_propagation_plan(wfs::PyramidWFS) =
    pyramid_propagation_plan(pyramid_propagation(wfs))
@inline pyramid_propagation_workspace(wfs::PyramidWFS) =
    pyramid_propagation_workspace(pyramid_propagation(wfs))
@inline pyramid_propagation_workspace(
    front_end::PyramidOpticalFrontEnd) =
    pyramid_propagation_workspace(front_end.propagation)
@inline four_pupil_propagation_workspace(
    front_end::PyramidOpticalFrontEnd) =
    pyramid_propagation_workspace(front_end)
@inline pyramid_phase_mask(wfs::PyramidWFS) =
    wfs.front_end.phase_mask
"""
    pyramid_focal_mask(wfs::PyramidWFS)

Borrow the prepared sampled complex focal-plane transmission used by this
Pyramid sensor. The returned array is plant-owned storage; copy or prepare a
calibration plan from it before changing the sensor's sampling geometry.
"""
@inline pyramid_focal_mask(wfs::PyramidWFS) =
    pyramid_propagation_workspace(wfs).pyramid_mask
@inline pyramid_operating_modulation(wfs::PyramidWFS) =
    wfs.front_end.modulation

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

@inline function _pyramid_shift_pixels(value)
    value isa Real && isfinite(value) || throw(InvalidConfiguration(
        "pyramid pupil shifts must be finite real values"))
    try
        return round(Int, value)
    catch error
        error isa InexactError || rethrow()
        throw(InvalidConfiguration(
            "pyramid pupil shifts must round to Int"))
    end
end

@inline _pyramid_face_shifts_pixels(value::Real) =
    ntuple(_ -> _pyramid_shift_pixels(value), 4)

function _pyramid_face_shifts_pixels(values)
    (values isa Tuple || values isa AbstractVector) || throw(
        InvalidConfiguration(
            "pyramid per-face pupil shifts must be a scalar or four-element tuple/vector"))
    length(values) == 4 || throw(InvalidConfiguration(
        "pyramid per-face pupil shifts must have four elements"))
    return ntuple(index -> _pyramid_shift_pixels(values[index]), 4)
end

"""
    PyramidWFS(tel; ...)

Construct a pyramid wavefront sensor.

The diffractive model forms four re-imaged pupil intensities through a
focal-plane pyramid mask. It publishes the physical four-pupil photon-rate
map; detector acquisition and RTC signal estimation are separate owners.
`phase_mask_rotation_rad` sets the rotation angle applied to the physical
mask-coordinate transform, while
`modulation_phase_offset_rad` selects the circular modulation quadrature
origin. Both values are in radians. `modulation_propagation_strategy` selects
the exact pupil-tilt reference formulation or the approximate shifted-mask
formulation for fixed prepared modulation. `pupil_shift_x_pixels` and
`pupil_shift_y_pixels` specify immutable per-face q1–q4 detector-pixel shifts;
a scalar value applies to all four faces.
"""
function PyramidWFS(tel::Telescope; pupil_samples::Int, modulation::Real=2.0,
    modulation_points::Union{Int,Nothing}=nothing, extra_modulation_factor::Int=0,
    old_mask::Bool=false, rooftop::Real=0.0,
    phase_mask_rotation_rad::Real=0.0,
    modulation_phase_offset_rad::Real=0.0,
    pupil_shift_x_pixels=0,
    pupil_shift_y_pixels=0,
    user_modulation_path=nothing, mask_scale::Real=1.0, diffraction_padding::Int=2,
    psf_centering::Bool=true, n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
    modulation_propagation_strategy::AbstractPyramidModulationPropagationStrategy=
        PyramidPupilTiltStrategy(),
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel))

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
    shift_x = _pyramid_face_shifts_pixels(pupil_shift_x_pixels)
    shift_y = _pyramid_face_shifts_pixels(pupil_shift_y_pixels)
    phase_mask = PyramidPhaseMask{T}(
        old_mask,
        T(rooftop),
        typed_phase_mask_rotation_rad,
        typed_mask_scale,
        diffraction_padding, psf_centering, n_pix_separation, n_pix_edge,
        shift_x, shift_y)
    operating_policy = legacy_modulation_policy(T(modulation),
        modulation_points, extra_modulation_factor,
        typed_modulation_phase_offset_rad,
        user_modulation_path)
    front_end, acquisition = _prepare_pyramid_diffractive_storage(
        backend, T, tel, phase_mask, operating_policy, pupil_samples,
        binning, modulation_propagation_strategy)
    wfs = PyramidWFS{typeof(front_end),typeof(acquisition),typeof(selector)}(
        front_end, acquisition)
    prepare_pyramid_front_end!(wfs, tel)
    return wfs
end

function _prepare_pyramid_diffractive_storage(backend, ::Type{T}, tel,
    phase_mask, operating_policy, pupil_samples,
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
    modulation_batch = _prepare_pyramid_modulation_batch(
        style, field, phase_mask, modulation,
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
    front_end = PyramidOpticalFrontEnd(phase_mask, modulation, propagation,
        pupil_samples, binning, nothing)
    subaperture_pixels = div(tel.params.resolution, pupil_samples)
    nominal = div(pad, subaperture_pixels)
    camera_frame = backend{T}(undef, nominal, nominal)
    acquisition = PyramidDetectorAcquisition(PyramidAcquisitionPlan(binning),
        PyramidAcquisitionWorkspace(nominal),
        PyramidAcquisitionProducts(camera_frame))
    return front_end, acquisition
end

function prepare_pyramid_front_end!(wfs::PyramidWFS, tel::Telescope)
    build_pyramid_phasor!(pyramid_propagation_workspace(wfs).phasor)
    pupil = PupilFunction(tel)
    build_pyramid_mask!(wfs, pupil)
    build_pyramid_shifted_masks!(wfs, pupil)
    return nothing
end

function PyramidOpticalFrontEnd(sensor::PyramidWFS,
    source=nothing)
    front_end = sensor.front_end
    return PyramidOpticalFrontEnd(front_end.phase_mask, front_end.modulation,
        front_end.propagation,
        front_end.pupil_samples, front_end.binning, source)
end

@inline function pyramid_front_end_with_source(
    front_end::PyramidOpticalFrontEnd, source)
    return PyramidOpticalFrontEnd(front_end.phase_mask,
        front_end.modulation, front_end.propagation,
        front_end.pupil_samples, front_end.binning,
        source)
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
    pupil_samples = wfs.front_end.pupil_samples
    pad = _pupil_resolution(pupil) * wfs.front_end.phase_mask.diffraction_padding
    if wfs.front_end.phase_mask.n_pix_separation !== nothing
        pixels_per_pupil_sample = div(_pupil_resolution(pupil), pupil_samples)
        pad = pyramid_native_frame_size(pupil_samples,
            wfs.front_end.phase_mask.n_pix_separation,
            wfs.front_end.phase_mask.n_pix_edge) * pixels_per_pupil_sample
    end
    pad >= _pupil_resolution(pupil) || throw(InvalidConfiguration(
        "pyramid padding must be >= telescope resolution"))
    pad % pyramid_acquisition_plan(wfs).binning == 0 || throw(
        InvalidConfiguration("pyramid binning must evenly divide padded resolution"))
    _pupil_resolution(pupil) % pyramid_acquisition_plan(wfs).binning == 0 ||
        throw(InvalidConfiguration(
            "pyramid binning must evenly divide telescope resolution"))
    ensure_pyramid_buffers!(wfs, pad, pupil)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, pupil::PupilFunction, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    acquisition_workspace = pyramid_acquisition_workspace(wfs)
    acquisition_products = pyramid_acquisition_products(wfs)
    binning = pyramid_acquisition_plan(wfs).binning
    sub = div(_pupil_resolution(pupil), wfs.front_end.pupil_samples)
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
    return frame
end
