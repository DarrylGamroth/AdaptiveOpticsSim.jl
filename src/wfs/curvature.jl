#
# Curvature wavefront sensing
#
# This sensor is modeled as a twin-defocus propagated pupil measurement:
#
# 1. embed the complex pupil field on a padded FFT grid
# 2. apply equal-and-opposite defocus masks to form two branches
# 3. propagate both branches with a batched centered FFT
# 4. crop back to the pupil support and bin to the exported readout
# 5. apply branch response terms and form the normalized signal
#    `(I⁺ - I⁻) / (I⁺ + I⁻)`
#
# The exported `slopes(wfs)` product contains one scalar signal per valid
# subaperture and implements the controller's family-neutral WFS signal
# contract.
#

@kernel function curvature_phasor_kernel!(phasor, scale, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds phasor[i, j] = cis(scale * (i + j - 2))
    end
end

@kernel function curvature_signal_from_frame_kernel!(signal_2d, slopes, frame_plus, frame_minus, reference_signal_2d, valid_mask,
    epsval, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        if @inbounds valid_mask[i, j]
            plus = @inbounds frame_plus[i, j]
            minus = @inbounds frame_minus[i, j]
            total = plus + minus
            corrected = total > epsval ?
                (plus - minus) / (total + epsval) - @inbounds(reference_signal_2d[i, j]) :
                zero(eltype(signal_2d))
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

@kernel function curvature_branch_field_stack_kernel!(field_stack, amplitude, opd, defocus_stack, phasor,
    amp_scale, opd_to_cycles, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= size(field_stack, 3)
        xi = x - ox
        yi = y - oy
        val = zero(eltype(field_stack))
        if 1 <= xi <= n && 1 <= yi <= n
            val = amp_scale * @inbounds(amplitude[xi, yi]) *
                cispi(opd_to_cycles * @inbounds(opd[xi, yi]))
        end
        @inbounds field_stack[x, y, branch] = val * defocus_stack[x, y, branch] * phasor[x, y]
    end
end

@kernel function curvature_branch_field_from_input_kernel!(field_stack, input_field, defocus_stack, phasor, pad::Int, n_branches::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= n_branches
        @inbounds field_stack[x, y, branch] = input_field[x, y] * defocus_stack[x, y, branch] * phasor[x, y]
    end
end

@kernel function curvature_abs2_stack_kernel!(intensity_stack, field_stack, pad::Int, n_branches::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= n_branches
        @inbounds intensity_stack[x, y, branch] = abs2(field_stack[x, y, branch])
    end
end

@kernel function curvature_sample_branch_kernel!(out, intensity_stack, ox::Int, oy::Int, binning::Int,
    branch::Int, throughput, background, n_out::Int, m_out::Int)
    i, j = @index(Global, NTuple)
    if i <= n_out && j <= m_out
        acc = zero(eltype(out))
        base_i = ox + (i - 1) * binning
        base_j = oy + (j - 1) * binning
        @inbounds for ii in 1:binning, jj in 1:binning
            acc += intensity_stack[base_i + ii, base_j + jj, branch]
        end
        @inbounds out[i, j] = throughput * acc + background
    end
end

@kernel function curvature_frame_pack_kernel!(camera_frame, frame_plus, frame_minus, side::Int)
    i, j = @index(Global, NTuple)
    if i <= side && j <= side
        @inbounds begin
            camera_frame[i, j] = frame_plus[i, j]
            camera_frame[i + side, j] = frame_minus[i, j]
        end
    end
end

@kernel function curvature_frame_unpack_kernel!(frame_plus, frame_minus, camera_frame, side::Int)
    i, j = @index(Global, NTuple)
    if i <= side && j <= size(camera_frame, 2)
        @inbounds begin
            frame_plus[i, j] = camera_frame[i, j]
            frame_minus[i, j] = camera_frame[i + side, j]
        end
    end
end

@kernel function curvature_channel_pack_kernel!(camera_frame, frame_plus, frame_minus, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        @inbounds begin
            camera_frame[1, idx] = frame_plus[i, j]
            camera_frame[2, idx] = frame_minus[i, j]
        end
    end
end

@kernel function curvature_signal_from_channels_kernel!(signal_2d, slopes, frame, reference_signal_2d, valid_mask,
    epsval, n_sub::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = (i - 1) * n_sub + j
        if @inbounds valid_mask[i, j]
            plus = @inbounds frame[1, idx]
            minus = @inbounds frame[2, idx]
            total = plus + minus
            corrected = total > epsval ?
                (plus - minus) / (total + epsval) - @inbounds(reference_signal_2d[i, j]) :
                zero(eltype(signal_2d))
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
a configurable readout surface. The 1-D `state.slopes` vector is the normalized
curvature signal sampled on the valid subaperture grid.
"""
abstract type CurvatureReadoutModel end
struct CurvatureFrameReadout <: CurvatureReadoutModel end

"""
    CurvatureChannelReadout()

Map the positive- and negative-defocus Curvature WFS branches to two values per
pupil sample for a compatible channel or counting detector.
"""
struct CurvatureChannelReadout <: CurvatureReadoutModel end

"""
    CurvatureBranchResponse(; plus_throughput=1, minus_throughput=1, plus_background=0, minus_background=0, T=Float64)

Per-branch optical relay rate adjustment for `CurvatureWFS`.

The throughputs scale the intra-/extra-focal branch photon-rate planes. The
backgrounds are added after optical sampling and represent fixed relay
stray-light or rate-offset approximations. Detector presampling response,
physical-pixel integration, quantum efficiency, exposure, stochastic response,
and readout belong to the downstream detector-acquisition stage instead.
"""
struct CurvatureBranchResponse{T<:AbstractFloat}
    plus_throughput::T
    minus_throughput::T
    plus_background::T
    minus_background::T
end

function CurvatureBranchResponse(; plus_throughput::Real=1.0, minus_throughput::Real=1.0,
    plus_background::Real=0.0, minus_background::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    plus_throughput >= 0 || throw(InvalidConfiguration("plus_throughput must be >= 0"))
    minus_throughput >= 0 || throw(InvalidConfiguration("minus_throughput must be >= 0"))
    plus_background >= 0 || throw(InvalidConfiguration("plus_background must be >= 0"))
    minus_background >= 0 || throw(InvalidConfiguration("minus_background must be >= 0"))
    return CurvatureBranchResponse{T}(T(plus_throughput), T(minus_throughput), T(plus_background), T(minus_background))
end

"""Immutable configuration for differential Curvature signal estimation."""
struct CurvatureEstimatorParams{T<:AbstractFloat}
    pupil_samples::Int
    threshold::T
end

"""Run-immutable contract for two-branch defocus propagation and sampling."""
struct CurvaturePropagationPlan{
    D<:CurvatureDefocusPair,T<:AbstractFloat,R<:CurvatureBranchResponse{T}}
    defocus_pair::D
    pupil_resolution::Int
    pupil_samples::Int
    readout_crop_resolution::Int
    readout_pixels_per_sample::Int
    branch_response::R
end

"""
Replaceable single-writer diffraction and sampling workspace for one Curvature
propagation plan.
"""
mutable struct CurvaturePropagationWorkspace{
    T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    CS<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pfs,
}
    phasor::C
    field_stack::CS
    defocus_stack::CS
    intensity_stack::RS
    cropped_plus::R
    cropped_minus::R
    frame_plus::R
    frame_minus::R
    fft_stack_plan::Pfs
end

"""Exact plan/workspace owner for one Curvature propagation execution."""
struct PreparedCurvaturePropagation{
    P<:CurvaturePropagationPlan,W<:CurvaturePropagationWorkspace}
    plan::P
    workspace::W
end

@inline curvature_propagation_plan(
    propagation::PreparedCurvaturePropagation) = propagation.plan
@inline curvature_propagation_workspace(
    propagation::PreparedCurvaturePropagation) = propagation.workspace

"""Prepared intra-/extra-focal Curvature optical front end."""
struct CurvatureOpticalFrontEnd{D,P,S}
    defocus_pair::D
    propagation::P
    source::S
end

"""Run-immutable detector-facing Curvature readout mapping."""
struct CurvatureAcquisitionPlan{R<:CurvatureReadoutModel}
    readout_model::R
end

"""Caller-visible packed detector-facing Curvature rate product."""
mutable struct CurvatureAcquisitionProducts{R<:AbstractMatrix}
    frame::R
end

struct CurvatureDetectorAcquisition{P,PR}
    plan::P
    products::PR
end

"""Persistent support and reference-calibration state."""
mutable struct CurvatureEstimatorState{
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

"""Replaceable single-writer reduction and signal-image scratch."""
mutable struct CurvatureEstimatorWorkspace{
    T<:AbstractFloat,R<:AbstractMatrix{T}}
    reduced_plus::R
    reduced_minus::R
    signal_2d::R
end

"""Caller-visible differential Curvature signal product."""
mutable struct CurvatureEstimatorProducts{
    T<:AbstractFloat,V<:AbstractVector{T}}
    signal::V
end

struct CurvatureDifferentialEstimator{P<:CurvatureEstimatorParams,S,W,PR}
    params::P
    state::S
    workspace::W
    products::PR
end

struct CurvatureWFS{F,A,E,B<:AbstractArrayBackend} <:
    AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::CurvatureWFS{F,A,E,B}) where {F,A,E,B} = B()

@inline curvature_estimator_params(wfs::CurvatureWFS) = wfs.estimator.params
@inline curvature_estimator_state(wfs::CurvatureWFS) = wfs.estimator.state
@inline curvature_estimator_workspace(wfs::CurvatureWFS) =
    wfs.estimator.workspace
@inline curvature_estimator_products(wfs::CurvatureWFS) =
    wfs.estimator.products
@inline curvature_acquisition_plan(wfs::CurvatureWFS) = wfs.acquisition.plan
@inline curvature_acquisition_products(wfs::CurvatureWFS) =
    wfs.acquisition.products
@inline curvature_propagation(wfs::CurvatureWFS) = wfs.front_end.propagation
@inline curvature_propagation_plan(wfs::CurvatureWFS) =
    curvature_propagation_plan(curvature_propagation(wfs))
@inline curvature_propagation_workspace(wfs::CurvatureWFS) =
    curvature_propagation_workspace(curvature_propagation(wfs))
@inline curvature_propagation_workspace(front_end::CurvatureOpticalFrontEnd) =
    curvature_propagation_workspace(front_end.propagation)

"""
    CurvatureWFS(tel; pupil_samples, threshold=0.1, defocus_rms_nm=500.0, diffraction_padding=2,
                 readout_crop_resolution=tel.params.resolution, readout_pixels_per_sample=1,
                 readout_model=CurvatureFrameReadout(), branch_response=CurvatureBranchResponse(), ...)

Construct a curvature WFS using two propagated defocus branches.

`defocus_rms_nm` sets the RMS amplitude of the defocus term in nanometers.
`branch_response` models intra-/extra-focal throughput and background imbalance.
`readout_crop_resolution` selects the square detector-plane crop before export and
`readout_pixels_per_sample` controls the exported frame sampling for frame-style
readout.
"""
function CurvatureWFS(tel::Telescope; pupil_samples::Int, threshold::Real=0.1, defocus_rms_nm::Real=500.0,
    diffraction_padding::Int=2, readout_model::CurvatureReadoutModel=CurvatureFrameReadout(),
    readout_crop_resolution::Integer=tel.params.resolution, readout_pixels_per_sample::Integer=1,
    branch_response::CurvatureBranchResponse=CurvatureBranchResponse(),
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    tel.params.resolution % pupil_samples == 0 ||
        throw(InvalidConfiguration("telescope resolution must be divisible by pupil_samples"))
    diffraction_padding >= 1 ||
        throw(InvalidConfiguration("diffraction_padding must be >= 1"))
    defocus_rms_nm >= 0 ||
        throw(InvalidConfiguration("defocus_rms_nm must be >= 0"))
    readout_crop_resolution > 0 ||
        throw(InvalidConfiguration("readout_crop_resolution must be > 0"))
    readout_crop_resolution <= tel.params.resolution ||
        throw(InvalidConfiguration("readout_crop_resolution must be <= telescope resolution"))
    readout_pixels_per_sample >= 1 ||
        throw(InvalidConfiguration("readout_pixels_per_sample must be >= 1"))
    readout_side = pupil_samples * Int(readout_pixels_per_sample)
    readout_crop_resolution % readout_side == 0 ||
        throw(InvalidConfiguration("readout_crop_resolution must be divisible by pupil_samples * readout_pixels_per_sample"))
    validate_curvature_readout_geometry(readout_model, Int(readout_pixels_per_sample))
    n = tel.params.resolution
    pad = n * diffraction_padding
    response = convert_curvature_branch_response(branch_response, T)
    defocus_pair = CurvatureDefocusPair(T(defocus_rms_nm),
        diffraction_padding)
    propagation_plan = CurvaturePropagationPlan(defocus_pair, n,
        pupil_samples, Int(readout_crop_resolution),
        Int(readout_pixels_per_sample), response)
    estimator_params = CurvatureEstimatorParams(pupil_samples, T(threshold))
    valid_mask = backend{Bool}(undef, pupil_samples, pupil_samples)
    slopes = backend{T}(undef, pupil_samples * pupil_samples)
    n_branches = 2
    phasor = backend{Complex{T}}(undef, pad, pad)
    field_stack = backend{Complex{T}}(undef, pad, pad, n_branches)
    defocus_stack = similar(field_stack)
    intensity_stack = backend{T}(undef, pad, pad, n_branches)
    cropped_plus = backend{T}(undef, propagation_plan.readout_crop_resolution,
        propagation_plan.readout_crop_resolution)
    cropped_minus = similar(cropped_plus)
    frame_side = propagation_plan.pupil_samples *
        propagation_plan.readout_pixels_per_sample
    frame_plus = backend{T}(undef, frame_side, frame_side)
    frame_minus = backend{T}(undef, frame_side, frame_side)
    reduced_plus = backend{T}(undef, pupil_samples, pupil_samples)
    reduced_minus = similar(reduced_plus)
    signal_2d = backend{T}(undef, pupil_samples, pupil_samples)
    reference_signal_2d = similar(signal_2d)
    camera_frame = curvature_camera_frame(backend, T, pupil_samples, readout_model;
        readout_pixels_per_sample=Int(readout_pixels_per_sample))
    fill!(slopes, zero(T))
    fill!(field_stack, zero(eltype(field_stack)))
    fill!(intensity_stack, zero(T))
    fill!(cropped_plus, zero(T))
    fill!(cropped_minus, zero(T))
    fill!(frame_plus, zero(T))
    fill!(frame_minus, zero(T))
    fill!(reduced_plus, zero(T))
    fill!(reduced_minus, zero(T))
    fill!(signal_2d, zero(T))
    fill!(reference_signal_2d, zero(T))
    fill!(camera_frame, zero(T))
    fft_stack_plan = plan_fft_backend!(field_stack, (1, 2))
    propagation_workspace = CurvaturePropagationWorkspace{
        T,
        typeof(phasor),
        typeof(field_stack),
        typeof(frame_plus),
        typeof(intensity_stack),
        typeof(fft_stack_plan),
    }(
        phasor,
        field_stack,
        defocus_stack,
        intensity_stack,
        cropped_plus,
        cropped_minus,
        frame_plus,
        frame_minus,
        fft_stack_plan,
    )
    propagation = PreparedCurvaturePropagation(propagation_plan,
        propagation_workspace)
    acquisition = CurvatureDetectorAcquisition(
        CurvatureAcquisitionPlan(readout_model),
        CurvatureAcquisitionProducts(camera_frame))
    estimator_state = CurvatureEstimatorState{
        T,
        typeof(valid_mask),
        typeof(signal_2d),
    }(
        valid_mask,
        reference_signal_2d,
        false,
        zero(T),
        zero(UInt),
        zero(UInt),
    )
    estimator_workspace = CurvatureEstimatorWorkspace(reduced_plus,
        reduced_minus, signal_2d)
    estimator_products = CurvatureEstimatorProducts(slopes)
    estimator = CurvatureDifferentialEstimator(estimator_params,
        estimator_state, estimator_workspace, estimator_products)
    front_end = CurvatureOpticalFrontEnd(defocus_pair, propagation, nothing)
    wfs = CurvatureWFS{
        typeof(front_end),typeof(acquisition),typeof(estimator),typeof(selector),
    }(front_end, acquisition, estimator)
    initial_pupil = PupilFunction(tel; T=T, backend=selector)
    update_valid_mask!(wfs, initial_pupil)
    build_curvature_phasor!(curvature_propagation_workspace(wfs).phasor)
    build_curvature_defocus_masks!(wfs)
    return wfs
end

sensing_mode(::CurvatureWFS) = Diffractive()
validate_curvature_readout_geometry(::CurvatureFrameReadout, readout_pixels_per_sample::Int) = nothing
function validate_curvature_readout_geometry(::CurvatureChannelReadout, readout_pixels_per_sample::Int)
    readout_pixels_per_sample == 1 ||
        throw(InvalidConfiguration("CurvatureChannelReadout requires readout_pixels_per_sample == 1"))
    return nothing
end
convert_curvature_branch_response(response::CurvatureBranchResponse, ::Type{T}) where {T<:AbstractFloat} =
    CurvatureBranchResponse(T=T, plus_throughput=response.plus_throughput, minus_throughput=response.minus_throughput,
        plus_background=response.plus_background, minus_background=response.minus_background)
curvature_camera_dims(pupil_samples::Int, ::CurvatureFrameReadout) = (2 * pupil_samples, pupil_samples)
curvature_camera_dims(pupil_samples::Int, ::CurvatureChannelReadout) = (2, pupil_samples * pupil_samples)
curvature_camera_dims(pupil_samples::Int, readout_pixels_per_sample::Int, ::CurvatureFrameReadout) =
    (2 * pupil_samples * readout_pixels_per_sample, pupil_samples * readout_pixels_per_sample)
curvature_camera_dims(pupil_samples::Int, readout_pixels_per_sample::Int, ::CurvatureChannelReadout) = (2, pupil_samples * pupil_samples)
curvature_camera_frame(backend, ::Type{T}, pupil_samples::Int, readout_model::CurvatureReadoutModel;
    readout_pixels_per_sample::Int=1) where {T<:AbstractFloat} =
    backend{T}(undef, curvature_camera_dims(pupil_samples, readout_pixels_per_sample, readout_model)...)

function update_valid_mask!(wfs::CurvatureWFS, pupil::PupilFunction)
    _require_curvature_pupil_geometry(wfs, pupil)
    params = curvature_estimator_params(wfs)
    state = curvature_estimator_state(wfs)
    set_valid_subapertures!(state.valid_mask, pupil.support,
        params.threshold)
    return wfs
end

@inline function _require_curvature_pupil_geometry(wfs::CurvatureWFS,
    pupil::PupilFunction)
    resolution = curvature_propagation_plan(wfs).pupil_resolution
    pupil.metadata.dimensions == (resolution, resolution) || throw(
        DimensionMismatchError(
            "CurvatureWFS PupilFunction dimensions do not match its prepared pupil grid"))
    return nothing
end

function _prepare_curvature_atmospheric_propagation(wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource,
    atm::AbstractAtmosphere, model::AbstractAtmosphericFieldModel)
    plan = curvature_propagation_plan(wfs)
    workspace = curvature_propagation_workspace(wfs)
    return AtmosphericFieldPropagation(atm, pupil, src;
        model=model,
        zero_padding=plan.defocus_pair.diffraction_padding,
        T=eltype(workspace.frame_plus))
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

function host_curvature_defocus_stack(wfs::CurvatureWFS)
    plan = curvature_propagation_plan(wfs)
    workspace = curvature_propagation_workspace(wfs)
    T = eltype(workspace.frame_plus)
    n = plan.pupil_resolution
    pad = size(workspace.field_stack, 1)
    cx = (pad + 1) / 2
    radius = T(n) / 2
    out = Array{Complex{T}}(undef, pad, pad, 2)
    inv_waves = plan.defocus_pair.rms_nm / T(1e9)
    @inbounds for j in 1:pad, i in 1:pad
        x = (T(i) - T(cx)) / radius
        y = (T(j) - T(cx)) / radius
        r2 = x * x + y * y
        z4 = r2 <= one(T) ? sqrt(T(3)) * (T(2) * r2 - one(T)) : zero(T)
        phase_cycles = inv_waves * z4
        out[i, j, 1] = cispi(T(2) * phase_cycles)
        out[i, j, 2] = cispi(-T(2) * phase_cycles)
    end
    return out
end

function build_curvature_defocus_masks!(wfs::CurvatureWFS)
    copyto!(curvature_propagation_workspace(wfs).defocus_stack,
        host_curvature_defocus_stack(wfs))
    return wfs
end

function sample_curvature_frames!(wfs::CurvatureWFS)
    frame = curvature_acquisition_products(wfs).frame
    return sample_curvature_frames!(
        execution_style(frame), wfs)
end

function sample_curvature_frames!(::ScalarCPUStyle, wfs::CurvatureWFS)
    plan = curvature_propagation_plan(wfs)
    workspace = curvature_propagation_workspace(wfs)
    crop_n = plan.readout_crop_resolution
    pad = size(workspace.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    sub = div(crop_n,
        plan.pupil_samples * plan.readout_pixels_per_sample)
    copyto!(workspace.cropped_plus,
        @view(workspace.intensity_stack[ox+1:ox+crop_n,
            oy+1:oy+crop_n, 1]))
    copyto!(workspace.cropped_minus,
        @view(workspace.intensity_stack[ox+1:ox+crop_n,
            oy+1:oy+crop_n, 2]))
    bin2d!(workspace.frame_plus, workspace.cropped_plus, sub)
    bin2d!(workspace.frame_minus, workspace.cropped_minus, sub)
    apply_curvature_branch_response!(wfs)
    pack_curvature_readout!(wfs)
    return curvature_acquisition_products(wfs).frame
end

function sample_curvature_frames!(style::AcceleratorStyle,
    wfs::CurvatureWFS)
    plan = curvature_propagation_plan(wfs)
    workspace = curvature_propagation_workspace(wfs)
    crop_n = plan.readout_crop_resolution
    pad = size(workspace.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    sub = div(crop_n,
        plan.pupil_samples * plan.readout_pixels_per_sample)
    response = plan.branch_response
    n_out, m_out = size(workspace.frame_plus)
    launch_kernel_async!(style, curvature_sample_branch_kernel!,
        workspace.frame_plus, workspace.intensity_stack,
        ox, oy, sub, 1, response.plus_throughput, response.plus_background, n_out, m_out;
        ndrange=size(workspace.frame_plus))
    launch_kernel_async!(style, curvature_sample_branch_kernel!,
        workspace.frame_minus, workspace.intensity_stack,
        ox, oy, sub, 2, response.minus_throughput, response.minus_background, n_out, m_out;
        ndrange=size(workspace.frame_minus))
    pack_curvature_readout!(wfs)
    return curvature_acquisition_products(wfs).frame
end

function apply_curvature_branch_response!(wfs::CurvatureWFS)
    plan = curvature_propagation_plan(wfs)
    workspace = curvature_propagation_workspace(wfs)
    response = plan.branch_response
    @. workspace.frame_plus = response.plus_throughput *
        workspace.frame_plus + response.plus_background
    @. workspace.frame_minus = response.minus_throughput *
        workspace.frame_minus + response.minus_background
    return wfs
end

function pack_curvature_readout!(wfs::CurvatureWFS)
    frame = curvature_acquisition_products(wfs).frame
    return pack_curvature_readout!(execution_style(frame),
        curvature_acquisition_plan(wfs).readout_model, wfs)
end

function pack_curvature_readout!(::ScalarCPUStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS)
    workspace = curvature_propagation_workspace(wfs)
    frame = curvature_acquisition_products(wfs).frame
    side = size(workspace.frame_plus, 1)
    @views copyto!(frame[1:side, :], workspace.frame_plus)
    @views copyto!(frame[side+1:2*side, :], workspace.frame_minus)
    return frame
end

function pack_curvature_readout!(::ScalarCPUStyle, ::CurvatureChannelReadout, wfs::CurvatureWFS)
    workspace = curvature_propagation_workspace(wfs)
    frame = curvature_acquisition_products(wfs).frame
    copyto!(@view(frame[1, :]), vec(workspace.frame_plus))
    copyto!(@view(frame[2, :]), vec(workspace.frame_minus))
    return frame
end

function pack_curvature_readout!(style::AcceleratorStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS)
    workspace = curvature_propagation_workspace(wfs)
    frame = curvature_acquisition_products(wfs).frame
    side = size(workspace.frame_plus, 1)
    launch_kernel!(style, curvature_frame_pack_kernel!, frame,
        workspace.frame_plus, workspace.frame_minus, side;
        ndrange=size(workspace.frame_plus))
    return frame
end

function pack_curvature_readout!(style::AcceleratorStyle, ::CurvatureChannelReadout, wfs::CurvatureWFS)
    workspace = curvature_propagation_workspace(wfs)
    frame = curvature_acquisition_products(wfs).frame
    n_sub = curvature_estimator_params(wfs).pupil_samples
    launch_kernel!(style, curvature_channel_pack_kernel!, frame,
        workspace.frame_plus, workspace.frame_minus, n_sub;
        ndrange=size(workspace.frame_plus))
    return frame
end

function curvature_intensity!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource)
    require_leaf_source(src, "CurvatureWFS")
    return curvature_intensity!(execution_style(
        curvature_acquisition_products(wfs).frame), wfs, pupil, src)
end

function curvature_intensity!(::ScalarCPUStyle, wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource)
    _require_curvature_pupil_geometry(wfs, pupil)
    workspace = curvature_propagation_workspace(wfs)
    n = pupil.metadata.dimensions[1]
    pad = size(workspace.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(workspace.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(workspace.frame_plus)(
        photon_irradiance(src) * _pupil_cell_area(pupil)
    ))
    amplitude = pupil.amplitude
    fill!(workspace.field_stack, zero(eltype(workspace.field_stack)))
    @inbounds for y in 1:n, x in 1:n
        if amplitude[x, y] > zero(eltype(amplitude))
            val = amp_scale * amplitude[x, y] *
                cispi(opd_to_cycles * pupil.opd[x, y])
            xx = ox + x
            yy = oy + y
            common = val * workspace.phasor[xx, yy]
            workspace.field_stack[xx, yy, 1] = common * workspace.defocus_stack[xx, yy, 1]
            workspace.field_stack[xx, yy, 2] = common * workspace.defocus_stack[xx, yy, 2]
        end
    end
    fraunhofer_intensity_stack!(workspace.intensity_stack,
        workspace.field_stack, workspace.fft_stack_plan)
    return sample_curvature_frames!(wfs)
end

function curvature_intensity!(style::AcceleratorStyle, wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource)
    _require_curvature_pupil_geometry(wfs, pupil)
    workspace = curvature_propagation_workspace(wfs)
    n = pupil.metadata.dimensions[1]
    pad = size(workspace.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(workspace.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(workspace.frame_plus)(
        photon_irradiance(src) * _pupil_cell_area(pupil)
    ))
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_branch_field_stack_kernel!,
        workspace.field_stack, pupil.amplitude, pupil.opd,
        workspace.defocus_stack, workspace.phasor, amp_scale, opd_to_cycles, ox, oy,
        n, pad;
        ndrange=size(workspace.field_stack))
    finish_kernel_phase!(phase)
    fraunhofer_intensity_stack!(workspace.intensity_stack,
        workspace.field_stack, workspace.fft_stack_plan)
    return sample_curvature_frames!(wfs)
end

function curvature_branch_stack_from_field!(wfs::CurvatureWFS, field::ElectricField)
    return curvature_branch_stack_from_field!(
        execution_style(curvature_propagation_workspace(wfs).field_stack),
        wfs, field)
end

function curvature_branch_stack_from_field!(::ScalarCPUStyle, wfs::CurvatureWFS, field::ElectricField)
    workspace = curvature_propagation_workspace(wfs)
    size(field.values) == (size(workspace.field_stack, 1),
        size(workspace.field_stack, 2)) ||
        throw(DimensionMismatchError("ElectricField padded resolution must match CurvatureWFS diffraction grid"))
    n_branches = size(workspace.field_stack, 3)
    pad = size(workspace.field_stack, 1)
    @inbounds for branch in 1:n_branches, y in 1:pad, x in 1:pad
        workspace.field_stack[x, y, branch] = field.values[x, y] *
            workspace.defocus_stack[x, y, branch] * workspace.phasor[x, y]
    end
    return workspace.field_stack
end

function curvature_branch_stack_from_field!(style::AcceleratorStyle, wfs::CurvatureWFS, field::ElectricField)
    workspace = curvature_propagation_workspace(wfs)
    size(field.values) == (size(workspace.field_stack, 1),
        size(workspace.field_stack, 2)) ||
        throw(DimensionMismatchError("ElectricField padded resolution must match CurvatureWFS diffraction grid"))
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_branch_field_from_input_kernel!, workspace.field_stack,
        field.values, workspace.defocus_stack, workspace.phasor,
        size(workspace.field_stack, 1), size(workspace.field_stack, 3);
        ndrange=size(workspace.field_stack))
    finish_kernel_phase!(phase)
    return workspace.field_stack
end

function curvature_intensity_from_field!(wfs::CurvatureWFS,
    pupil::PupilFunction, field::ElectricField)
    _require_curvature_pupil_geometry(wfs, pupil)
    workspace = curvature_propagation_workspace(wfs)
    curvature_branch_stack_from_field!(wfs, field)
    fraunhofer_intensity_stack!(workspace.intensity_stack,
        workspace.field_stack, workspace.fft_stack_plan)
    return sample_curvature_frames!(wfs)
end

function curvature_intensity!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, atm::AbstractAtmosphere;
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(
        T=eltype(curvature_propagation_workspace(wfs).frame_plus)))
    require_leaf_source(src, "atmosphere-aware CurvatureWFS")
    prop = isnothing(propagation) ?
        _prepare_curvature_atmospheric_propagation(wfs, pupil, src, atm,
            model) : propagation
    field = propagate_atmosphere_field!(prop, atm, current_epoch(atm))
    return curvature_intensity_from_field!(wfs, pupil, field)
end

function curvature_signal!(wfs::CurvatureWFS, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    return curvature_signal!(execution_style(frame),
        curvature_acquisition_plan(wfs).readout_model, wfs, frame)
end

function curvature_signal!(::ScalarCPUStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(curvature_acquisition_products(wfs).frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled camera frame"))
    unpack_curvature_frame!(wfs, frame)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(wfs)
end

function curvature_signal!(style::AcceleratorStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(curvature_acquisition_products(wfs).frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled camera frame"))
    unpack_curvature_frame!(wfs, frame)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(style, wfs)
end

function unpack_curvature_frame!(wfs::CurvatureWFS, frame::AbstractMatrix)
    return _unpack_curvature_frame!(execution_style(frame), wfs, frame)
end

function _unpack_curvature_frame!(::ScalarCPUStyle, wfs::CurvatureWFS, frame::AbstractMatrix)
    workspace = curvature_propagation_workspace(wfs)
    side = size(workspace.frame_plus, 1)
    @views copyto!(workspace.frame_plus, frame[1:side, :])
    @views copyto!(workspace.frame_minus, frame[side+1:2*side, :])
    return wfs
end

function _unpack_curvature_frame!(style::AcceleratorStyle, wfs::CurvatureWFS,
    frame::AbstractMatrix)
    workspace = curvature_propagation_workspace(wfs)
    side = size(workspace.frame_plus, 1)
    launch_kernel!(style, curvature_frame_unpack_kernel!, workspace.frame_plus,
        workspace.frame_minus, frame, side; ndrange=size(workspace.frame_plus))
    return wfs
end

function reduce_curvature_frame_signal!(wfs::CurvatureWFS)
    factor = curvature_propagation_plan(wfs).readout_pixels_per_sample
    propagation_workspace = curvature_propagation_workspace(wfs)
    estimator_workspace = curvature_estimator_workspace(wfs)
    if factor == 1
        copyto!(estimator_workspace.reduced_plus,
            propagation_workspace.frame_plus)
        copyto!(estimator_workspace.reduced_minus,
            propagation_workspace.frame_minus)
    else
        bin2d!(estimator_workspace.reduced_plus,
            propagation_workspace.frame_plus, factor)
        bin2d!(estimator_workspace.reduced_minus,
            propagation_workspace.frame_minus, factor)
    end
    return wfs
end

function curvature_signal_from_planes!(wfs::CurvatureWFS)
    n_sub = curvature_estimator_params(wfs).pupil_samples
    state = curvature_estimator_state(wfs)
    workspace = curvature_estimator_workspace(wfs)
    products = curvature_estimator_products(wfs)
    epsval = eps(eltype(workspace.signal_2d))
    @inbounds for i in 1:n_sub, j in 1:n_sub
        idx = (i - 1) * n_sub + j
        if state.valid_mask[i, j]
            plus = workspace.reduced_plus[i, j]
            minus = workspace.reduced_minus[i, j]
            total = plus + minus
            corrected = total > epsval ?
                (plus - minus) / (total + epsval) - state.reference_signal_2d[i, j] :
                zero(eltype(workspace.signal_2d))
            workspace.signal_2d[i, j] = corrected
            products.signal[idx] = corrected
        else
            workspace.signal_2d[i, j] = zero(eltype(workspace.signal_2d))
            products.signal[idx] = zero(eltype(products.signal))
        end
    end
    return products.signal
end

function curvature_signal_from_planes!(style::AcceleratorStyle, wfs::CurvatureWFS)
    n_sub = curvature_estimator_params(wfs).pupil_samples
    state = curvature_estimator_state(wfs)
    workspace = curvature_estimator_workspace(wfs)
    products = curvature_estimator_products(wfs)
    epsval = eps(eltype(workspace.signal_2d))
    launch_kernel!(style, curvature_signal_from_frame_kernel!,
        workspace.signal_2d, products.signal, workspace.reduced_plus,
        workspace.reduced_minus, state.reference_signal_2d,
        state.valid_mask, epsval, n_sub;
        ndrange=size(workspace.signal_2d))
    return products.signal
end

function curvature_signal_from_current_frames!(wfs::CurvatureWFS)
    return curvature_signal_from_current_frames!(execution_style(
        curvature_acquisition_products(wfs).frame), wfs)
end

function curvature_signal_from_current_frames!(::ScalarCPUStyle, wfs::CurvatureWFS)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(wfs)
end

function curvature_signal_from_current_frames!(style::AcceleratorStyle, wfs::CurvatureWFS)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(style, wfs)
end

function curvature_signal!(::ScalarCPUStyle, ::CurvatureChannelReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(curvature_acquisition_products(wfs).frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled channel readout"))
    n_sub = curvature_estimator_params(wfs).pupil_samples
    state = curvature_estimator_state(wfs)
    workspace = curvature_estimator_workspace(wfs)
    products = curvature_estimator_products(wfs)
    epsval = eps(eltype(workspace.signal_2d))
    @inbounds for i in 1:n_sub, j in 1:n_sub
        idx = (i - 1) * n_sub + j
        if state.valid_mask[i, j]
            plus = frame[1, idx]
            minus = frame[2, idx]
            total = plus + minus
            corrected = total > epsval ?
                (plus - minus) / (total + epsval) - state.reference_signal_2d[i, j] :
                zero(eltype(workspace.signal_2d))
            workspace.signal_2d[i, j] = corrected
            products.signal[idx] = corrected
        else
            workspace.signal_2d[i, j] = zero(eltype(workspace.signal_2d))
            products.signal[idx] = zero(eltype(products.signal))
        end
    end
    return products.signal
end

function curvature_signal!(style::AcceleratorStyle, ::CurvatureChannelReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(curvature_acquisition_products(wfs).frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled channel readout"))
    n_sub = curvature_estimator_params(wfs).pupil_samples
    state = curvature_estimator_state(wfs)
    workspace = curvature_estimator_workspace(wfs)
    products = curvature_estimator_products(wfs)
    epsval = eps(eltype(workspace.signal_2d))
    launch_kernel!(style, curvature_signal_from_channels_kernel!,
        workspace.signal_2d, products.signal, frame,
        state.reference_signal_2d, state.valid_mask, epsval, n_sub;
        ndrange=size(workspace.signal_2d))
    return products.signal
end

function ensure_curvature_calibration!(wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "CurvatureWFS calibration")
    state = curvature_estimator_state(wfs)
    workspace = curvature_estimator_workspace(wfs)
    products = curvature_estimator_products(wfs)
    frame = curvature_acquisition_products(wfs).frame
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
        curvature_intensity!(wfs, pupil, src)
        fill!(state.reference_signal_2d,
            zero(eltype(state.reference_signal_2d)))
        curvature_signal!(wfs, frame)
        copyto!(state.reference_signal_2d, workspace.signal_2d)
    finally
        restore_opd!(pupil, opd_saved)
    end
    state.calibrated = true
    state.calibration_wavelength = λ
    state.calibration_signature = sig
    state.calibration_revision += UInt(1)
    return wfs
end


include("curvature/stages.jl")

function measure!(::Diffractive, wfs::CurvatureWFS,
    pupil::PupilFunction)
    throw(InvalidConfiguration(
        "CurvatureWFS requires a source; call measure!(wfs, pupil, src)."))
end

measure!(wfs::CurvatureWFS, pupil::PupilFunction) =
    measure!(sensing_mode(wfs), wfs, pupil)
measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource) = measure!(sensing_mode(wfs), wfs, pupil, src)

function measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, pupil, src, det; rng=rng)
end

function measure!(wfs::CurvatureWFS, ::PupilFunction, ::Asterism)
    throw(InvalidConfiguration("CurvatureWFS asterism support is not implemented"))
end

function measure!(wfs::CurvatureWFS, ::PupilFunction, ::Asterism,
    ::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration("CurvatureWFS asterism support is not implemented"))
end

function measure!(::Diffractive, wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource)
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, atm::AbstractAtmosphere;
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(curvature_propagation_workspace(wfs).frame_plus)))
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src, atm; propagation=propagation,
        model=model)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, atm::AbstractAtmosphere, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng(),
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(curvature_propagation_workspace(wfs).frame_plus)))
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src, atm; propagation=propagation,
        model=model)
    capture!(det, curvature_acquisition_products(wfs).frame, src; rng=rng)
    return curvature_signal!(wfs, output_frame(det))
end

function measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    ast::Asterism, atm::AbstractAtmosphere;
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(curvature_propagation_workspace(wfs).frame_plus)))
    common_source = common_wfs_calibration_source(ast, "CurvatureWFS")
    ensure_curvature_calibration!(wfs, pupil, common_source)
    propagation_workspace = curvature_propagation_workspace(wfs)
    acc_plus = similar(propagation_workspace.frame_plus)
    acc_minus = similar(propagation_workspace.frame_minus)
    fill!(acc_plus, zero(eltype(acc_plus)))
    fill!(acc_minus, zero(eltype(acc_minus)))
    @inbounds for src in ast.sources
        curvature_intensity!(wfs, pupil, src, atm; model=model)
        acc_plus .+= propagation_workspace.frame_plus
        acc_minus .+= propagation_workspace.frame_minus
    end
    copyto!(propagation_workspace.frame_plus, acc_plus)
    copyto!(propagation_workspace.frame_minus, acc_minus)
    pack_curvature_readout!(wfs)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, pupil::PupilFunction,
    ast::Asterism, atm::AbstractAtmosphere, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng(),
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(curvature_propagation_workspace(wfs).frame_plus)))
    common_source = common_wfs_calibration_source(ast, "CurvatureWFS")
    measure!(wfs, pupil, ast, atm; model=model)
    capture!(det, curvature_acquisition_products(wfs).frame, common_source;
        rng=rng)
    return curvature_signal!(wfs, output_frame(det))
end

function measure!(::Diffractive, wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    return measure_detector_coupled!(
        curvature_acquisition_plan(wfs).readout_model, wfs, pupil,
        src, det; rng=rng)
end

function measure_detector_coupled!(::CurvatureChannelReadout,
    wfs::CurvatureWFS, ::PupilFunction,
    src::AbstractSource, det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration(
        "CurvatureChannelReadout requires a linear-mode APD channel bank or a counting detector"))
end

function measure_detector_coupled!(::CurvatureChannelReadout,
    wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, det::AbstractCountingDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src)
    frame = curvature_acquisition_products(wfs).frame
    capture!(det, frame, src; rng=rng)
    size(output_frame(det)) == size(frame) ||
        throw(InvalidConfiguration("CurvatureWFS counting-detector output size must match the sampled channel readout"))
    return curvature_signal!(wfs, output_frame(det))
end

function measure_detector_coupled!(::CurvatureChannelReadout,
    wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, det::LinearAPDDetector;
    rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src)
    frame = curvature_acquisition_products(wfs).frame
    input = vec(frame)
    length(input) == length(channel_output(det)) || throw(
        InvalidConfiguration(
            "CurvatureWFS linear-APD channel count must match the sampled channel readout"))
    capture!(det, input; rng=rng)
    copyto!(frame, channel_output(det))
    return curvature_signal!(wfs, frame)
end

function measure_detector_coupled!(::CurvatureFrameReadout,
    wfs::CurvatureWFS, pupil::PupilFunction,
    src::AbstractSource, det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, pupil, src)
    curvature_intensity!(wfs, pupil, src)
    frame = curvature_acquisition_products(wfs).frame
    capture!(det, frame, src; rng=rng)
    size(output_frame(det)) == size(frame) ||
        throw(InvalidConfiguration("CurvatureWFS detector output size must match the sampled camera frame"))
    return curvature_signal!(wfs, output_frame(det))
end

@inline slopes(wfs::CurvatureWFS) = curvature_estimator_products(wfs).signal
@inline valid_subaperture_mask(wfs::CurvatureWFS) =
    curvature_estimator_state(wfs).valid_mask
@inline reference_signal(wfs::CurvatureWFS) =
    curvature_estimator_state(wfs).reference_signal_2d
@inline wfs_calibration_signature(wfs::CurvatureWFS) =
    curvature_estimator_state(wfs).calibration_signature

@inline wfs_output_metadata(wfs::CurvatureWFS) = wfs_output_metadata(
    curvature_acquisition_plan(wfs).readout_model, wfs)
@inline wfs_output_metadata(::CurvatureFrameReadout, wfs::CurvatureWFS) = nothing
@inline wfs_output_metadata(::CurvatureChannelReadout, wfs::CurvatureWFS) =
    ChannelReadoutMetadata(:branch_by_channel,
        size(curvature_acquisition_products(wfs).frame),
        length(curvature_acquisition_products(wfs).frame))

@inline supports_prepared_runtime(::CurvatureWFS, src::AbstractSource) =
    is_leaf_source(src)
@inline supports_detector_output(wfs::CurvatureWFS,
    det::AbstractDetector) = supports_detector_output(
    curvature_acquisition_plan(wfs).readout_model, det)
@inline supports_detector_output(::CurvatureFrameReadout, ::AbstractDetector) = true
@inline supports_detector_output(::CurvatureChannelReadout, ::AbstractDetector) = false
@inline supports_detector_output(::CurvatureChannelReadout, ::AbstractCountingDetector) = true
@inline supports_detector_output(::CurvatureChannelReadout,
    ::LinearAPDDetector) = true

@inline function prepare_runtime_wfs!(wfs::CurvatureWFS,
    pupil::PupilFunction, src::AbstractSource)
    require_leaf_source(src, "CurvatureWFS runtime preparation")
    ensure_curvature_calibration!(wfs, pupil, src)
    return wfs
end
