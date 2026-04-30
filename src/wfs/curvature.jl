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
# The exported `state.slopes` remains one scalar signal per valid subaperture
# so the runtime and reconstructor interfaces stay compatible with the rest of
# the package.
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

@kernel function curvature_branch_field_stack_kernel!(field_stack, pupil, opd, defocus_stack, phasor,
    amp_scale, opd_to_cycles, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, branch = @index(Global, NTuple)
    if x <= pad && y <= pad && branch <= size(field_stack, 3)
        xi = x - ox
        yi = y - oy
        val = zero(eltype(field_stack))
        if 1 <= xi <= n && 1 <= yi <= n && @inbounds(pupil[xi, yi])
            val = amp_scale * cispi(opd_to_cycles * @inbounds(opd[xi, yi]))
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

Curvature wavefront sensor using propagated intra-/extra-focal detector planes.

The sensor generates two defocused propagated intensity maps and exports them as
a configurable readout surface. The 1-D `state.slopes` vector is the normalized
curvature signal sampled on the valid subaperture grid.
"""
abstract type CurvatureReadoutModel end
struct CurvatureFrameReadout <: CurvatureReadoutModel end
struct CurvatureCountingReadout <: CurvatureReadoutModel end

"""
    CurvatureBranchResponse(; plus_throughput=1, minus_throughput=1, plus_background=0, minus_background=0, T=Float64)

Per-branch detector-plane response model for `CurvatureWFS`.

The throughputs scale the intra-/extra-focal branch intensities before signal
formation. The backgrounds are added after detector-plane sampling and can be
used to represent fixed branch offsets or asymmetric relay background.
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

struct CurvatureWFSParams{T<:AbstractFloat,R<:CurvatureReadoutModel,B<:CurvatureBranchResponse{T}}
    pupil_samples::Int
    threshold::T
    defocus_rms_nm::T
    diffraction_padding::Int
    readout_crop_resolution::Int
    readout_pixels_per_sample::Int
    readout_model::R
    branch_response::B
end

mutable struct CurvatureWFSState{
    T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    CS<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pfs,
    AP,
}
    valid_mask::A
    slopes::V
    phasor::C
    field_stack::CS
    defocus_stack::CS
    intensity_stack::RS
    cropped_plus::R
    cropped_minus::R
    frame_plus::R
    frame_minus::R
    reduced_plus::R
    reduced_minus::R
    signal_2d::R
    reference_signal_2d::R
    camera_frame::R
    fft_stack_plan::Pfs
    effective_padding::Int
    calibrated::Bool
    calibration_wavelength::T
    atmospheric_propagation::AP
    atmospheric_signature::UInt
end

struct CurvatureWFS{P<:CurvatureWFSParams,S<:CurvatureWFSState,B<:AbstractArrayBackend} <: AbstractWFS
    params::P
    state::S
end

@inline backend(::CurvatureWFS{<:Any,<:Any,B}) where {B} = B()

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
    params = CurvatureWFSParams{T, typeof(readout_model), typeof(response)}(
        pupil_samples, T(threshold), T(defocus_rms_nm), diffraction_padding,
        Int(readout_crop_resolution), Int(readout_pixels_per_sample), readout_model, response)
    valid_mask = backend{Bool}(undef, pupil_samples, pupil_samples)
    slopes = backend{T}(undef, pupil_samples * pupil_samples)
    n_branches = 2
    phasor = backend{Complex{T}}(undef, pad, pad)
    field_stack = backend{Complex{T}}(undef, pad, pad, n_branches)
    defocus_stack = similar(field_stack)
    intensity_stack = backend{T}(undef, pad, pad, n_branches)
    cropped_plus = backend{T}(undef, params.readout_crop_resolution, params.readout_crop_resolution)
    cropped_minus = backend{T}(undef, params.readout_crop_resolution, params.readout_crop_resolution)
    frame_side = params.pupil_samples * params.readout_pixels_per_sample
    frame_plus = backend{T}(undef, frame_side, frame_side)
    frame_minus = backend{T}(undef, frame_side, frame_side)
    reduced_plus = backend{T}(undef, params.pupil_samples, params.pupil_samples)
    reduced_minus = backend{T}(undef, params.pupil_samples, params.pupil_samples)
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
    state = CurvatureWFSState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(phasor),
        typeof(field_stack),
        typeof(frame_plus),
        typeof(intensity_stack),
        typeof(fft_stack_plan),
        Union{Nothing,AtmosphericFieldPropagation},
    }(
        valid_mask,
        slopes,
        phasor,
        field_stack,
        defocus_stack,
        intensity_stack,
        cropped_plus,
        cropped_minus,
        frame_plus,
        frame_minus,
        reduced_plus,
        reduced_minus,
        signal_2d,
        reference_signal_2d,
        camera_frame,
        fft_stack_plan,
        diffraction_padding,
        false,
        zero(T),
        nothing,
        zero(UInt),
    )
    wfs = CurvatureWFS{typeof(params),typeof(state),typeof(selector)}(params, state)
    update_valid_mask!(wfs, tel)
    build_curvature_phasor!(wfs.state.phasor)
    build_curvature_defocus_masks!(wfs, tel)
    return wfs
end

sensing_mode(::CurvatureWFS) = Diffractive()
validate_curvature_readout_geometry(::CurvatureFrameReadout, readout_pixels_per_sample::Int) = nothing
function validate_curvature_readout_geometry(::CurvatureCountingReadout, readout_pixels_per_sample::Int)
    readout_pixels_per_sample == 1 ||
        throw(InvalidConfiguration("CurvatureCountingReadout requires readout_pixels_per_sample == 1"))
    return nothing
end
convert_curvature_branch_response(response::CurvatureBranchResponse, ::Type{T}) where {T<:AbstractFloat} =
    CurvatureBranchResponse(T=T, plus_throughput=response.plus_throughput, minus_throughput=response.minus_throughput,
        plus_background=response.plus_background, minus_background=response.minus_background)
curvature_camera_dims(pupil_samples::Int, ::CurvatureFrameReadout) = (2 * pupil_samples, pupil_samples)
curvature_camera_dims(pupil_samples::Int, ::CurvatureCountingReadout) = (2, pupil_samples * pupil_samples)
curvature_camera_dims(pupil_samples::Int, readout_pixels_per_sample::Int, ::CurvatureFrameReadout) =
    (2 * pupil_samples * readout_pixels_per_sample, pupil_samples * readout_pixels_per_sample)
curvature_camera_dims(pupil_samples::Int, readout_pixels_per_sample::Int, ::CurvatureCountingReadout) = (2, pupil_samples * pupil_samples)
curvature_camera_frame(backend, ::Type{T}, pupil_samples::Int, readout_model::CurvatureReadoutModel;
    readout_pixels_per_sample::Int=1) where {T<:AbstractFloat} =
    backend{T}(undef, curvature_camera_dims(pupil_samples, readout_pixels_per_sample, readout_model)...)

function update_valid_mask!(wfs::CurvatureWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function ensure_curvature_buffers!(wfs::CurvatureWFS, tel::Telescope)
    n = tel.params.resolution
    pad = n * wfs.params.diffraction_padding
    crop_n = wfs.params.readout_crop_resolution
    frame_side = wfs.params.pupil_samples * wfs.params.readout_pixels_per_sample
    if size(wfs.state.field_stack, 1) != pad || size(wfs.state.field_stack, 2) != pad
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        n_branches = size(wfs.state.field_stack, 3)
        wfs.state.field_stack = similar(wfs.state.field_stack, pad, pad, n_branches)
        wfs.state.defocus_stack = similar(wfs.state.defocus_stack, pad, pad, n_branches)
        wfs.state.intensity_stack = similar(wfs.state.intensity_stack, pad, pad, n_branches)
        wfs.state.cropped_plus = similar(wfs.state.cropped_plus, crop_n, crop_n)
        wfs.state.cropped_minus = similar(wfs.state.cropped_minus, crop_n, crop_n)
        wfs.state.frame_plus = similar(wfs.state.frame_plus, frame_side, frame_side)
        wfs.state.frame_minus = similar(wfs.state.frame_minus, frame_side, frame_side)
        wfs.state.reduced_plus = similar(wfs.state.reduced_plus, wfs.params.pupil_samples, wfs.params.pupil_samples)
        wfs.state.reduced_minus = similar(wfs.state.reduced_minus, wfs.params.pupil_samples, wfs.params.pupil_samples)
        wfs.state.camera_frame = similar(wfs.state.camera_frame,
            curvature_camera_dims(wfs.params.pupil_samples, wfs.params.readout_pixels_per_sample, wfs.params.readout_model)...)
        wfs.state.fft_stack_plan = plan_fft_backend!(wfs.state.field_stack, (1, 2))
        wfs.state.effective_padding = wfs.params.diffraction_padding
        wfs.state.calibrated = false
        wfs.state.atmospheric_propagation = nothing
        wfs.state.atmospheric_signature = zero(UInt)
        build_curvature_phasor!(wfs.state.phasor)
        build_curvature_defocus_masks!(wfs, tel)
    end
    return wfs
end

function curvature_atmospheric_signature(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource,
    atm::AbstractAtmosphere, model::AbstractAtmosphericFieldModel)
    return hash((
        typeof(atm),
        typeof(src),
        source_measurement_signature(src),
        typeof(model),
        tel.params.resolution,
        tel.params.diameter,
        wfs.params.diffraction_padding,
        length(atm.layers),
    ))
end

@inline _curvature_cached_propagation(::Nothing, ::UInt, ::UInt) = nothing
@inline _curvature_cached_propagation(prop::AtmosphericFieldPropagation, signature::UInt, cached_signature::UInt) =
    cached_signature == signature ? prop : nothing

function ensure_curvature_atmospheric_propagation!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource,
    atm::AbstractAtmosphere, model::AbstractAtmosphericFieldModel)
    sig = curvature_atmospheric_signature(wfs, tel, src, atm, model)
    cache = _curvature_cached_propagation(wfs.state.atmospheric_propagation, sig, wfs.state.atmospheric_signature)
    if !isnothing(cache)
        return cache
    end
    prop = AtmosphericFieldPropagation(atm, tel, src;
        model=model,
        zero_padding=wfs.params.diffraction_padding,
        T=eltype(wfs.state.frame_plus))
    wfs.state.atmospheric_propagation = prop
    wfs.state.atmospheric_signature = sig
    return prop
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

function host_curvature_defocus_stack(wfs::CurvatureWFS, tel::Telescope)
    T = eltype(wfs.state.frame_plus)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    cx = (pad + 1) / 2
    radius = T(n) / 2
    out = Array{Complex{T}}(undef, pad, pad, 2)
    inv_waves = wfs.params.defocus_rms_nm / T(1e9)
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

function build_curvature_defocus_masks!(wfs::CurvatureWFS, tel::Telescope)
    copyto!(wfs.state.defocus_stack, host_curvature_defocus_stack(wfs, tel))
    return wfs
end

function sample_curvature_frames!(wfs::CurvatureWFS, tel::Telescope)
    return sample_curvature_frames!(execution_style(wfs.state.camera_frame), wfs, tel)
end

function sample_curvature_frames!(::ScalarCPUStyle, wfs::CurvatureWFS, tel::Telescope)
    crop_n = wfs.params.readout_crop_resolution
    pad = size(wfs.state.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    sub = div(crop_n, wfs.params.pupil_samples * wfs.params.readout_pixels_per_sample)
    copyto!(wfs.state.cropped_plus, @view(wfs.state.intensity_stack[ox+1:ox+crop_n, oy+1:oy+crop_n, 1]))
    copyto!(wfs.state.cropped_minus, @view(wfs.state.intensity_stack[ox+1:ox+crop_n, oy+1:oy+crop_n, 2]))
    bin2d!(wfs.state.frame_plus, wfs.state.cropped_plus, sub)
    bin2d!(wfs.state.frame_minus, wfs.state.cropped_minus, sub)
    apply_curvature_branch_response!(wfs)
    pack_curvature_readout!(wfs)
    return wfs.state.camera_frame
end

function sample_curvature_frames!(style::AcceleratorStyle, wfs::CurvatureWFS, tel::Telescope)
    crop_n = wfs.params.readout_crop_resolution
    pad = size(wfs.state.intensity_stack, 1)
    ox = div(pad - crop_n, 2)
    oy = div(pad - crop_n, 2)
    sub = div(crop_n, wfs.params.pupil_samples * wfs.params.readout_pixels_per_sample)
    response = wfs.params.branch_response
    n_out, m_out = size(wfs.state.frame_plus)
    launch_kernel_async!(style, curvature_sample_branch_kernel!, wfs.state.frame_plus, wfs.state.intensity_stack,
        ox, oy, sub, 1, response.plus_throughput, response.plus_background, n_out, m_out;
        ndrange=size(wfs.state.frame_plus))
    launch_kernel_async!(style, curvature_sample_branch_kernel!, wfs.state.frame_minus, wfs.state.intensity_stack,
        ox, oy, sub, 2, response.minus_throughput, response.minus_background, n_out, m_out;
        ndrange=size(wfs.state.frame_minus))
    pack_curvature_readout!(wfs)
    return wfs.state.camera_frame
end

function apply_curvature_branch_response!(wfs::CurvatureWFS)
    response = wfs.params.branch_response
    @. wfs.state.frame_plus = response.plus_throughput * wfs.state.frame_plus + response.plus_background
    @. wfs.state.frame_minus = response.minus_throughput * wfs.state.frame_minus + response.minus_background
    return wfs
end

function pack_curvature_readout!(wfs::CurvatureWFS)
    return pack_curvature_readout!(execution_style(wfs.state.camera_frame), wfs.params.readout_model, wfs)
end

function pack_curvature_readout!(::ScalarCPUStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS)
    side = size(wfs.state.frame_plus, 1)
    @views copyto!(wfs.state.camera_frame[1:side, :], wfs.state.frame_plus)
    @views copyto!(wfs.state.camera_frame[side+1:2*side, :], wfs.state.frame_minus)
    return wfs.state.camera_frame
end

function pack_curvature_readout!(::ScalarCPUStyle, ::CurvatureCountingReadout, wfs::CurvatureWFS)
    copyto!(@view(wfs.state.camera_frame[1, :]), vec(wfs.state.frame_plus))
    copyto!(@view(wfs.state.camera_frame[2, :]), vec(wfs.state.frame_minus))
    return wfs.state.camera_frame
end

function pack_curvature_readout!(style::AcceleratorStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS)
    side = size(wfs.state.frame_plus, 1)
    launch_kernel!(style, curvature_frame_pack_kernel!, wfs.state.camera_frame, wfs.state.frame_plus,
        wfs.state.frame_minus, side; ndrange=size(wfs.state.frame_plus))
    return wfs.state.camera_frame
end

function pack_curvature_readout!(style::AcceleratorStyle, ::CurvatureCountingReadout, wfs::CurvatureWFS)
    launch_kernel!(style, curvature_channel_pack_kernel!, wfs.state.camera_frame, wfs.state.frame_plus,
        wfs.state.frame_minus, wfs.params.pupil_samples; ndrange=size(wfs.state.frame_plus))
    return wfs.state.camera_frame
end

function curvature_intensity!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    return curvature_intensity!(execution_style(wfs.state.camera_frame), wfs, tel, src)
end

function curvature_intensity!(::ScalarCPUStyle, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.frame_plus)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))
    fill!(wfs.state.field_stack, zero(eltype(wfs.state.field_stack)))
    @inbounds for y in 1:n, x in 1:n
        if tel.state.pupil[x, y]
            val = amp_scale * cispi(opd_to_cycles * tel.state.opd[x, y])
            xx = ox + x
            yy = oy + y
            common = val * wfs.state.phasor[xx, yy]
            wfs.state.field_stack[xx, yy, 1] = common * wfs.state.defocus_stack[xx, yy, 1]
            wfs.state.field_stack[xx, yy, 2] = common * wfs.state.defocus_stack[xx, yy, 2]
        end
    end
    fraunhofer_intensity_stack!(wfs.state.intensity_stack, wfs.state.field_stack, wfs.state.fft_stack_plan)
    return sample_curvature_frames!(wfs, tel)
end

function curvature_intensity!(style::AcceleratorStyle, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_buffers!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field_stack, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    opd_to_cycles = eltype(wfs.state.frame_plus)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.frame_plus)(
        photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2
    ))
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_branch_field_stack_kernel!, wfs.state.field_stack, tel.state.pupil,
        tel.state.opd, wfs.state.defocus_stack, wfs.state.phasor, amp_scale, opd_to_cycles, ox, oy, n, pad;
        ndrange=size(wfs.state.field_stack))
    finish_kernel_phase!(phase)
    fraunhofer_intensity_stack!(wfs.state.intensity_stack, wfs.state.field_stack, wfs.state.fft_stack_plan)
    return sample_curvature_frames!(wfs, tel)
end

function curvature_branch_stack_from_field!(wfs::CurvatureWFS, field::ElectricField)
    return curvature_branch_stack_from_field!(execution_style(wfs.state.field_stack), wfs, field)
end

function curvature_branch_stack_from_field!(::ScalarCPUStyle, wfs::CurvatureWFS, field::ElectricField)
    size(field.state.field) == (size(wfs.state.field_stack, 1), size(wfs.state.field_stack, 2)) ||
        throw(DimensionMismatchError("ElectricField padded resolution must match CurvatureWFS diffraction grid"))
    n_branches = size(wfs.state.field_stack, 3)
    pad = size(wfs.state.field_stack, 1)
    @inbounds for branch in 1:n_branches, y in 1:pad, x in 1:pad
        wfs.state.field_stack[x, y, branch] = field.state.field[x, y] * wfs.state.defocus_stack[x, y, branch] * wfs.state.phasor[x, y]
    end
    return wfs.state.field_stack
end

function curvature_branch_stack_from_field!(style::AcceleratorStyle, wfs::CurvatureWFS, field::ElectricField)
    size(field.state.field) == (size(wfs.state.field_stack, 1), size(wfs.state.field_stack, 2)) ||
        throw(DimensionMismatchError("ElectricField padded resolution must match CurvatureWFS diffraction grid"))
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, curvature_branch_field_from_input_kernel!, wfs.state.field_stack,
        field.state.field, wfs.state.defocus_stack, wfs.state.phasor, size(wfs.state.field_stack, 1), size(wfs.state.field_stack, 3);
        ndrange=size(wfs.state.field_stack))
    finish_kernel_phase!(phase)
    return wfs.state.field_stack
end

function curvature_intensity_from_field!(wfs::CurvatureWFS, tel::Telescope, field::ElectricField)
    ensure_curvature_buffers!(wfs, tel)
    curvature_branch_stack_from_field!(wfs, field)
    fraunhofer_intensity_stack!(wfs.state.intensity_stack, wfs.state.field_stack, wfs.state.fft_stack_plan)
    return sample_curvature_frames!(wfs, tel)
end

function curvature_intensity!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource, atm::AbstractAtmosphere;
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(wfs.state.frame_plus)))
    prop = isnothing(propagation) ? ensure_curvature_atmospheric_propagation!(wfs, tel, src, atm, model) : propagation
    field = propagate_atmosphere_field!(prop, atm, tel, src)
    return curvature_intensity_from_field!(wfs, tel, field)
end

function curvature_signal!(wfs::CurvatureWFS, frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    return curvature_signal!(execution_style(frame), wfs.params.readout_model, wfs, frame)
end

function curvature_signal!(::ScalarCPUStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(wfs.state.camera_frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled camera frame"))
    unpack_curvature_frame!(wfs, frame)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(wfs)
end

function curvature_signal!(style::AcceleratorStyle, ::CurvatureFrameReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(wfs.state.camera_frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled camera frame"))
    unpack_curvature_frame!(wfs, frame)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(style, wfs)
end

function unpack_curvature_frame!(wfs::CurvatureWFS, frame::AbstractMatrix)
    side = size(wfs.state.frame_plus, 1)
    @views copyto!(wfs.state.frame_plus, frame[1:side, :])
    @views copyto!(wfs.state.frame_minus, frame[side+1:2*side, :])
    return wfs
end

function reduce_curvature_frame_signal!(wfs::CurvatureWFS)
    factor = wfs.params.readout_pixels_per_sample
    if factor == 1
        copyto!(wfs.state.reduced_plus, wfs.state.frame_plus)
        copyto!(wfs.state.reduced_minus, wfs.state.frame_minus)
    else
        bin2d!(wfs.state.reduced_plus, wfs.state.frame_plus, factor)
        bin2d!(wfs.state.reduced_minus, wfs.state.frame_minus, factor)
    end
    return wfs
end

function curvature_signal_from_planes!(wfs::CurvatureWFS)
    n_sub = wfs.params.pupil_samples
    epsval = eps(eltype(wfs.state.signal_2d))
    @inbounds for i in 1:n_sub, j in 1:n_sub
        idx = (i - 1) * n_sub + j
        if wfs.state.valid_mask[i, j]
            plus = wfs.state.reduced_plus[i, j]
            minus = wfs.state.reduced_minus[i, j]
            signal = (plus - minus) / (plus + minus + epsval)
            corrected = signal - wfs.state.reference_signal_2d[i, j]
            wfs.state.signal_2d[i, j] = corrected
            wfs.state.slopes[idx] = corrected
        else
            wfs.state.signal_2d[i, j] = zero(eltype(wfs.state.signal_2d))
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
        end
    end
    return wfs.state.slopes
end

function curvature_signal_from_planes!(style::AcceleratorStyle, wfs::CurvatureWFS)
    n_sub = wfs.params.pupil_samples
    epsval = eps(eltype(wfs.state.signal_2d))
    launch_kernel!(style, curvature_signal_from_frame_kernel!, wfs.state.signal_2d, wfs.state.slopes,
        wfs.state.reduced_plus, wfs.state.reduced_minus, wfs.state.reference_signal_2d, wfs.state.valid_mask, epsval, n_sub;
        ndrange=size(wfs.state.signal_2d))
    return wfs.state.slopes
end

function curvature_signal_from_current_frames!(wfs::CurvatureWFS)
    return curvature_signal_from_current_frames!(execution_style(wfs.state.camera_frame), wfs)
end

function curvature_signal_from_current_frames!(::ScalarCPUStyle, wfs::CurvatureWFS)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(wfs)
end

function curvature_signal_from_current_frames!(style::AcceleratorStyle, wfs::CurvatureWFS)
    reduce_curvature_frame_signal!(wfs)
    return curvature_signal_from_planes!(style, wfs)
end

function curvature_signal!(::ScalarCPUStyle, ::CurvatureCountingReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(wfs.state.camera_frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled channel readout"))
    n_sub = wfs.params.pupil_samples
    epsval = eps(eltype(wfs.state.signal_2d))
    @inbounds for i in 1:n_sub, j in 1:n_sub
        idx = (i - 1) * n_sub + j
        if wfs.state.valid_mask[i, j]
            plus = frame[1, idx]
            minus = frame[2, idx]
            signal = (plus - minus) / (plus + minus + epsval)
            corrected = signal - wfs.state.reference_signal_2d[i, j]
            wfs.state.signal_2d[i, j] = corrected
            wfs.state.slopes[idx] = corrected
        else
            wfs.state.signal_2d[i, j] = zero(eltype(wfs.state.signal_2d))
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
        end
    end
    return wfs.state.slopes
end

function curvature_signal!(style::AcceleratorStyle, ::CurvatureCountingReadout, wfs::CurvatureWFS,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(frame) == size(wfs.state.camera_frame) ||
        throw(DimensionMismatchError("CurvatureWFS frame size must match the sampled channel readout"))
    n_sub = wfs.params.pupil_samples
    epsval = eps(eltype(wfs.state.signal_2d))
    launch_kernel!(style, curvature_signal_from_channels_kernel!, wfs.state.signal_2d, wfs.state.slopes,
        frame, wfs.state.reference_signal_2d, wfs.state.valid_mask, epsval, n_sub;
        ndrange=size(wfs.state.signal_2d))
    return wfs.state.slopes
end

function ensure_curvature_calibration!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    λ = calibration_wavelength(src, typeof(src.params.wavelength))
    if calibration_matches(wfs.state.calibrated, wfs.state.calibration_wavelength, λ)
        return wfs
    end
    opd_saved = copy(tel.state.opd)
    reset_opd!(tel)
    update_valid_mask!(wfs, tel)
    curvature_intensity!(wfs, tel, src)
    curvature_signal!(wfs, wfs.state.camera_frame)
    copyto!(wfs.state.reference_signal_2d, wfs.state.signal_2d)
    copyto!(tel.state.opd, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
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
    curvature_intensity!(wfs, tel, src)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource, atm::AbstractAtmosphere;
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(wfs.state.frame_plus)))
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src, atm; propagation=propagation, model=model)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource, atm::AbstractAtmosphere, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng(),
    propagation::Union{Nothing,AtmosphericFieldPropagation}=nothing,
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(wfs.state.frame_plus)))
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src, atm; propagation=propagation, model=model)
    capture!(det, wfs.state.camera_frame; rng=rng)
    return curvature_signal!(wfs, output_frame(det))
end

function measure!(wfs::CurvatureWFS, tel::Telescope, ast::Asterism, atm::AbstractAtmosphere;
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(wfs.state.frame_plus)))
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    ensure_curvature_calibration!(wfs, tel, ast.sources[1])
    acc_plus = similar(wfs.state.frame_plus)
    acc_minus = similar(wfs.state.frame_minus)
    fill!(acc_plus, zero(eltype(acc_plus)))
    fill!(acc_minus, zero(eltype(acc_minus)))
    @inbounds for src in ast.sources
        curvature_intensity!(wfs, tel, src, atm; model=model)
        acc_plus .+= wfs.state.frame_plus
        acc_minus .+= wfs.state.frame_minus
    end
    copyto!(wfs.state.frame_plus, acc_plus)
    copyto!(wfs.state.frame_minus, acc_minus)
    pack_curvature_readout!(wfs)
    return curvature_signal_from_current_frames!(wfs)
end

function measure!(wfs::CurvatureWFS, tel::Telescope, ast::Asterism, atm::AbstractAtmosphere, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng(),
    model::AbstractAtmosphericFieldModel=LayeredFresnelAtmosphericPropagation(T=eltype(wfs.state.frame_plus)))
    measure!(wfs, tel, ast, atm; model=model)
    capture!(det, wfs.state.camera_frame; rng=rng)
    return curvature_signal!(wfs, output_frame(det))
end

function measure!(::Diffractive, wfs::CurvatureWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    return measure_detector_coupled!(wfs.params.readout_model, wfs, tel, src, det; rng=rng)
end

function measure_detector_coupled!(::CurvatureCountingReadout, wfs::CurvatureWFS, tel::Telescope,
    src::AbstractSource, det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    throw(InvalidConfiguration("CurvatureWFS detector coupling requires CurvatureFrameReadout"))
end

function measure_detector_coupled!(::CurvatureCountingReadout, wfs::CurvatureWFS, tel::Telescope,
    src::AbstractSource, det::AbstractCountingDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src)
    capture!(det, wfs.state.camera_frame; rng=rng)
    size(output_frame(det)) == size(wfs.state.camera_frame) ||
        throw(InvalidConfiguration("CurvatureWFS counting-detector output size must match the sampled channel readout"))
    return curvature_signal!(wfs, output_frame(det))
end

function measure_detector_coupled!(::CurvatureFrameReadout, wfs::CurvatureWFS, tel::Telescope,
    src::AbstractSource, det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    ensure_curvature_calibration!(wfs, tel, src)
    curvature_intensity!(wfs, tel, src)
    capture!(det, wfs.state.camera_frame; rng=rng)
    size(output_frame(det)) == size(wfs.state.camera_frame) ||
        throw(InvalidConfiguration("CurvatureWFS detector output size must match the sampled camera frame"))
    return curvature_signal!(wfs, output_frame(det))
end

@inline valid_subaperture_mask(wfs::CurvatureWFS) = wfs.state.valid_mask
@inline reference_signal(wfs::CurvatureWFS) = wfs.state.reference_signal_2d
@inline camera_frame(wfs::CurvatureWFS) = wfs.state.camera_frame

@inline wfs_output_frame(wfs::CurvatureWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame(wfs::CurvatureWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::CurvatureWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::CurvatureWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_metadata(wfs::CurvatureWFS) = wfs_output_metadata(wfs.params.readout_model, wfs)
@inline wfs_output_metadata(::CurvatureFrameReadout, wfs::CurvatureWFS) = nothing
@inline wfs_output_metadata(::CurvatureCountingReadout, wfs::CurvatureWFS) =
    CountingReadoutMetadata(:branch_by_channel, size(wfs.state.camera_frame), length(wfs.state.camera_frame))

@inline supports_prepared_runtime(::CurvatureWFS, ::AbstractSource) = true
@inline supports_detector_output(wfs::CurvatureWFS, det::AbstractDetector) = supports_detector_output(wfs.params.readout_model, det)
@inline supports_detector_output(::CurvatureFrameReadout, ::AbstractDetector) = true
@inline supports_detector_output(::CurvatureCountingReadout, ::AbstractDetector) = false
@inline supports_detector_output(::CurvatureCountingReadout, ::AbstractCountingDetector) = true

@inline function prepare_runtime_wfs!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepropagate_runtime_wfs!(wfs::CurvatureWFS, atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource, rng::AbstractRNG)
    reset_opd!(tel)
    advance!(atm, tel, rng)
    apply!(optic, tel, DMAdditive())
    return nothing
end

@inline function measure_runtime_wfs!(wfs::CurvatureWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource, rng::AbstractRNG)
    measure!(wfs, tel, src, atm)
    return nothing
end

@inline function measure_runtime_wfs!(wfs::CurvatureWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    measure!(wfs, tel, src, atm, det; rng=rng)
    return nothing
end

@inline function finish_runtime_wfs_sensing!(wfs::CurvatureWFS, atm::AbstractAtmosphere,
    tel::Telescope, src::AbstractSource)
    tel.state.opd .+= atm.state.opd .* tel.state.pupil
    return nothing
end
