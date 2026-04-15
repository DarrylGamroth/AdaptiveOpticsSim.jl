using Statistics

#
# Shack-Hartmann wavefront sensing
#
# Two sensing models live in this file:
#
# - `Geometric()`: local wavefront gradients are sampled directly.
# - `Diffractive()`: each subaperture field is propagated to the focal plane,
#   intensity is sampled/cropped, and centroids are converted to slopes.
#
# The diffractive path uses FFT-based Fraunhofer propagation on each lenslet.
# For LGS sensing, elongated spots are handled through focal-plane convolution.
# For asterisms and GPU execution, the implementation batches lenslet/source
# stacks so the algorithm stays mathematically the same while reducing launch
# and detector-processing overhead.
#
"""
    ShackHartmann

Shack-Hartmann wavefront sensor with geometric and diffractive sensing modes.

The diffractive model computes focal-plane lenslet spots and converts their
centroids into x/y slopes. The geometric model bypasses image formation and
samples the OPD gradient directly.
"""
struct ShackHartmannParams{T<:AbstractFloat,VP<:AbstractValidSubaperturePolicy}
    n_subap::Int
    threshold::T
    threshold_cog::T
    threshold_convolution::T
    half_pixel_shift::Bool
    diffraction_padding::Int
    pixel_scale::Union{T,Nothing}
    n_pix_subap::Union{Int,Nothing}
    shannon_sampling::Bool
    valid_subaperture_policy::VP
end

mutable struct ShackHartmannState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    CC<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RC3<:AbstractArray{T,3},
    RT3<:AbstractArray{T,3},
    RB<:AbstractMatrix{T},
    RS<:AbstractMatrix{T},
    RCS<:AbstractArray{T,3},
    RC<:AbstractArray{T,3},
    RCE<:AbstractArray{T,3},
    RA<:AbstractArray{T,3},
    RN<:AbstractArray{T,3},
    P,
    PS,
    Pi,
    PSi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3},
    L,
    CAL}
    valid_mask::A
    slopes::V
    field::C
    phasor::C
    fft_buffer::C
    fft_stack::CC
    intensity::R
    intensity_stack::RC3
    intensity_tmp_stack::RT3
    temp::R
    bin_buffer::RB
    spot::RS
    sampled_spot_cube::RCS
    spot_cube::RC
    exported_spot_cube::RCE
    spot_cube_accum::RA
    detector_noise_cube::RN
    fft_plan::P
    fft_stack_plan::PS
    ifft_plan::Pi
    ifft_stack_plan::PSi
    elongation_kernel::K
    lgs_kernel_fft::KF
    lgs_kernel_tag::UInt
    effective_padding::Int
    binning_pixel_scale::Int
    sampled_n_pix_subap::Int
    phasor_ratio::T
    reference_signal_2d::RB
    fft_asterism_stack::CC
    fft_asterism_plan::PS
    asterism_capacity::Int
    amp_scales::V
    amp_scales_host::Vector{T}
    opd_to_cycles::V
    opd_to_cycles_host::Vector{T}
    spot_stats::V
    spot_stats_accum::V
    valid_mask_host::Matrix{Bool}
    reference_signal_host::Vector{T}
    slopes_host::Vector{T}
    centroid_host::Matrix{T}
    slopes_units::T
    export_pixels_enabled::Bool
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    layout::L
    calibration::CAL
end

struct ShackHartmann{M<:SensingMode,P<:ShackHartmannParams,S<:ShackHartmannState,B<:AbstractArrayBackend} <: AbstractWFS
    params::P
    state::S
end

@inline backend(::ShackHartmann{<:Any,<:Any,<:Any,B}) where {B} = B()

"""
    ShackHartmann(tel; ...)

Construct a Shack-Hartmann WFS on the telescope pupil grid.

Important diffractive quantities:

- `diffraction_padding` controls the focal-plane FFT grid size
- `pixel_scale` and `shannon_sampling` determine detector-plane sampling
- `n_pix_subap` controls the cropped spot size used for centroiding

The constructor allocates the full lenslet workspace, including stacked FFT and
spot buffers used by the batched GPU/runtime paths.
"""
function ShackHartmann(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    threshold_cog::Real=0.01, threshold_convolution::Real=0.05, half_pixel_shift::Bool=false,
    diffraction_padding::Int=2, pixel_scale=nothing, n_pix_subap=nothing,
    shannon_sampling::Bool=true,
    valid_subaperture_policy::Union{Nothing,AbstractValidSubaperturePolicy}=nothing,
    mode::SensingMode=Geometric(),
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    pixel_scale_t = pixel_scale === nothing ? nothing : T(pixel_scale)
    raw_policy = isnothing(valid_subaperture_policy) ? GeometryValidSubapertures(threshold=threshold, T=T) : valid_subaperture_policy
    policy = convert_valid_subaperture_policy(raw_policy, T)
    params = ShackHartmannParams{T, typeof(policy)}(n_subap, T(threshold), T(threshold_cog), T(threshold_convolution), half_pixel_shift, diffraction_padding,
        pixel_scale_t, n_pix_subap, shannon_sampling, policy)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    sub = div(tel.params.resolution, n_subap)
    pad = max(sub, sub * diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    phasor = similar(field)
    fft_buffer = similar(field)
    fft_stack = backend{Complex{T}}(undef, pad, pad, n_subap * n_subap)
    intensity = backend{T}(undef, pad, pad)
    intensity_stack = backend{T}(undef, pad, pad, n_subap * n_subap)
    intensity_tmp_stack = backend{T}(undef, pad, pad, n_subap * n_subap)
    temp = similar(intensity)
    bin_buffer = backend{T}(undef, sub, sub)
    spot = similar(bin_buffer)
    sampled_spot_cube = backend{T}(undef, n_subap * n_subap, sub, sub)
    spot_cube = similar(sampled_spot_cube)
    exported_spot_cube = similar(sampled_spot_cube)
    spot_cube_accum = similar(sampled_spot_cube)
    detector_noise_cube = similar(sampled_spot_cube)
    fft_plan = plan_fft_backend!(fft_buffer)
    fft_stack_plan = plan_fft_backend!(fft_stack, (1, 2))
    ifft_plan = plan_ifft_backend!(fft_buffer)
    ifft_stack_plan = plan_ifft_backend!(fft_stack, (1, 2))
    fft_asterism_stack = similar(fft_stack)
    fft_asterism_plan = plan_fft_backend!(fft_asterism_stack, (1, 2))
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0, 0)
    amp_scales = backend{T}(undef, 1)
    amp_scales_host = Vector{T}(undef, 1)
    opd_to_cycles = backend{T}(undef, 1)
    opd_to_cycles_host = Vector{T}(undef, 1)
    spot_stats = backend{T}(undef, 3 * n_subap * n_subap)
    spot_stats_accum = backend{T}(undef, 3 * n_subap * n_subap)
    valid_mask_host = Matrix{Bool}(undef, n_subap, n_subap)
    reference_signal_host = Vector{T}(undef, 2 * n_subap * n_subap)
    slopes_host = Vector{T}(undef, 2 * n_subap * n_subap)
    centroid_host = Matrix{T}(undef, sub, sub)
    reference_signal_2d = backend{T}(undef, 2 * n_subap, n_subap)
    layout = SubapertureLayout(n_subap, tel.params.resolution, tel.params.diameter, threshold, valid_mask, valid_mask_host)
    calibration = SubapertureCalibration(reference_signal_2d, reference_signal_host, CenterOfGravityExtraction(T(threshold_cog); T=T))
    state = ShackHartmannState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(field),
        typeof(fft_stack),
        typeof(intensity),
        typeof(intensity_stack),
        typeof(intensity_tmp_stack),
        typeof(bin_buffer),
        typeof(spot),
        typeof(sampled_spot_cube),
        typeof(spot_cube),
        typeof(exported_spot_cube),
        typeof(spot_cube_accum),
        typeof(detector_noise_cube),
        typeof(fft_plan),
        typeof(fft_stack_plan),
        typeof(ifft_plan),
        typeof(ifft_stack_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
        typeof(layout),
        typeof(calibration),
    }(
        valid_mask,
        slopes,
        field,
        phasor,
        fft_buffer,
        fft_stack,
        intensity,
        intensity_stack,
        intensity_tmp_stack,
        temp,
        bin_buffer,
        spot,
        sampled_spot_cube,
        spot_cube,
        exported_spot_cube,
        spot_cube_accum,
        detector_noise_cube,
        fft_plan,
        fft_stack_plan,
        ifft_plan,
        ifft_stack_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
        diffraction_padding,
        1,
        sub,
        T(NaN),
        reference_signal_2d,
        fft_asterism_stack,
        fft_asterism_plan,
        1,
        amp_scales,
        amp_scales_host,
        opd_to_cycles,
        opd_to_cycles_host,
        spot_stats,
        spot_stats_accum,
        valid_mask_host,
        reference_signal_host,
        slopes_host,
        centroid_host,
        one(T),
        true,
        false,
        zero(T),
        UInt(0),
        layout,
        calibration,
    )
    wfs = ShackHartmann{typeof(mode), typeof(params), typeof(state), typeof(selector)}(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ShackHartmann{M}) where {M} = M()
@inline subaperture_layout(wfs::ShackHartmann) = wfs.state.layout
@inline subaperture_calibration(wfs::ShackHartmann) = wfs.state.calibration
@inline valid_subaperture_policy(wfs::ShackHartmann) = wfs.params.valid_subaperture_policy
@inline slope_extraction_model(wfs::ShackHartmann) = slope_extraction_model(subaperture_calibration(wfs))
@inline centroid_threshold(wfs::ShackHartmann) = slope_extraction_model(wfs).threshold

abstract type AbstractShackHartmannSensingPlan end
struct ShackHartmannScalarPlan <: AbstractShackHartmannSensingPlan end
struct ShackHartmannBatchedPlan <: AbstractShackHartmannSensingPlan end
struct ShackHartmannDeviceStatsPlan <: AbstractShackHartmannSensingPlan end
struct ShackHartmannRocmSafePlan <: AbstractShackHartmannSensingPlan end

@inline sh_sensing_execution_plan(style::ExecutionStyle, wfs::ShackHartmann) =
    sh_sensing_execution_plan(typeof(style), typeof(wfs))
@inline sh_sensing_execution_plan(::Type{<:ScalarCPUStyle}, ::Type{<:ShackHartmann}) = ShackHartmannScalarPlan()
@inline sh_sensing_execution_plan(::Type{<:AcceleratorStyle}, ::Type{<:ShackHartmann}) = ShackHartmannBatchedPlan()
@inline sh_sensing_execution_plan(wfs::ShackHartmann) =
    sh_sensing_execution_plan(execution_style(wfs.state.slopes), wfs)

@inline sh_uses_rocm_safe_sensing_plan(::AbstractShackHartmannSensingPlan) = false
@inline sh_uses_rocm_safe_sensing_plan(::ShackHartmannRocmSafePlan) = true
@inline sh_uses_rocm_safe_sensing_plan(wfs::ShackHartmann) =
    sh_uses_rocm_safe_sensing_plan(sh_sensing_execution_plan(wfs))
@inline sh_uses_batched_sensing_plan(::AbstractShackHartmannSensingPlan) = false
@inline sh_uses_batched_sensing_plan(::ShackHartmannBatchedPlan) = true
@inline sh_uses_batched_sensing_plan(::ShackHartmannDeviceStatsPlan) = true
@inline sh_uses_batched_sensing_plan(wfs::ShackHartmann) =
    sh_uses_batched_sensing_plan(sh_sensing_execution_plan(wfs))
@inline sh_uses_device_stats_sensing_plan(::AbstractShackHartmannSensingPlan) = false
@inline sh_uses_device_stats_sensing_plan(::ShackHartmannDeviceStatsPlan) = true
@inline sh_uses_device_stats_sensing_plan(wfs::ShackHartmann) =
    sh_uses_device_stats_sensing_plan(sh_sensing_execution_plan(wfs))

convert_valid_subaperture_policy(policy::GeometryValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    GeometryValidSubapertures(threshold=T(policy.threshold), T=T)

convert_valid_subaperture_policy(policy::FluxThresholdValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    FluxThresholdValidSubapertures(light_ratio=T(policy.light_ratio), T=T)

@inline function sh_safe_peak_value(A::AbstractArray{T}) where {T<:AbstractFloat}
    return backend_maximum_value(A)
end

@inline function sh_refresh_valid_mask_host!(wfs::ShackHartmann)
    copyto!(wfs.state.valid_mask_host, wfs.state.valid_mask)
    return wfs.state.valid_mask_host
end

function update_valid_mask!(wfs::ShackHartmann, tel::Telescope)
    policy = valid_subaperture_policy(wfs)
    if policy isa GeometryValidSubapertures
        update_subaperture_layout!(subaperture_layout(wfs), tel.state.pupil, policy)
    elseif policy isa FluxThresholdValidSubapertures
        update_subaperture_layout!(subaperture_layout(wfs), tel.state.pupil_reflectivity, policy)
    else
        throw(UnsupportedAlgorithm("unsupported valid subaperture policy $(typeof(policy))"))
    end
    return wfs
end

@inline function sh_grouped_stack_capacity(wfs::ShackHartmann)
    return wfs.params.n_subap * wfs.params.n_subap * wfs.state.asterism_capacity
end

function ensure_sh_buffers!(wfs::ShackHartmann, pad::Int)
    n_spots = wfs.params.n_subap * wfs.params.n_subap
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.fft_buffer = similar(wfs.state.fft_buffer, pad, pad)
        wfs.state.fft_stack = similar(wfs.state.fft_stack, Complex{eltype(wfs.state.field)}, pad, pad, n_spots)
        wfs.state.intensity = similar(wfs.state.intensity, pad, pad)
        wfs.state.intensity_stack = similar(wfs.state.intensity_stack, eltype(wfs.state.intensity), pad, pad, n_spots)
        total = sh_grouped_stack_capacity(wfs)
        wfs.state.intensity_tmp_stack = similar(wfs.state.intensity_tmp_stack, eltype(wfs.state.intensity), pad, pad, total)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.fft_buffer)
        wfs.state.fft_stack_plan = plan_fft_backend!(wfs.state.fft_stack, (1, 2))
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.fft_buffer)
        wfs.state.ifft_stack_plan = plan_ifft_backend!(wfs.state.fft_stack, (1, 2))
        wfs.state.fft_asterism_stack = similar(wfs.state.fft_asterism_stack, pad, pad, total)
        wfs.state.fft_asterism_plan = plan_fft_backend!(wfs.state.fft_asterism_stack, (1, 2))
        wfs.state.phasor_ratio = eltype(wfs.state.slopes)(NaN)
        wfs.state.calibrated = false
    end
    return wfs
end

function ensure_sh_asterism_buffers!(wfs::ShackHartmann, n_sources::Int)
    n_sources > 0 || throw(InvalidConfiguration("asterism must contain at least one source"))
    if n_sources > wfs.state.asterism_capacity
        pad = size(wfs.state.fft_stack, 1)
        n_spots = wfs.params.n_subap * wfs.params.n_subap
        total = n_spots * n_sources
        wfs.state.fft_asterism_stack = similar(wfs.state.fft_asterism_stack, pad, pad, total)
        wfs.state.intensity_tmp_stack = similar(wfs.state.intensity_tmp_stack, eltype(wfs.state.intensity_tmp_stack), pad, pad, total)
        wfs.state.fft_asterism_plan = plan_fft_backend!(wfs.state.fft_asterism_stack, (1, 2))
        wfs.state.asterism_capacity = n_sources
    end
    if length(wfs.state.amp_scales) != n_sources
        wfs.state.amp_scales = similar(wfs.state.amp_scales, n_sources)
    end
    if length(wfs.state.amp_scales_host) != n_sources
        wfs.state.amp_scales_host = Vector{eltype(wfs.state.amp_scales_host)}(undef, n_sources)
    end
    if length(wfs.state.opd_to_cycles) != n_sources
        wfs.state.opd_to_cycles = similar(wfs.state.opd_to_cycles, n_sources)
    end
    if length(wfs.state.opd_to_cycles_host) != n_sources
        wfs.state.opd_to_cycles_host = Vector{eltype(wfs.state.opd_to_cycles_host)}(undef, n_sources)
    end
    return wfs
end

function build_sh_phasor!(wfs::ShackHartmann, ratio::T) where {T<:AbstractFloat}
    if size(wfs.state.phasor, 1) == 0
        return wfs
    end
    if isequal(wfs.state.phasor_ratio, ratio)
        return wfs
    end
    n = size(wfs.state.phasor, 1)
    scale = -T(π) * (T(n) + one(T) + ratio) / T(n)
    host = Matrix{Complex{T}}(undef, n, n)
    @inbounds for i in 1:n, j in 1:n
        host[i, j] = cis(scale * (i + j - 2))
    end
    copyto!(wfs.state.phasor, host)
    wfs.state.phasor_ratio = ratio
    return wfs
end

"""
    prepare_sampling!(wfs, tel, src)

Resolve the diffractive lenslet sampling for a source/wavelength pair.

This chooses the effective FFT padding, detector binning, and cropped spot size
so that the propagated spot grid is consistent with the requested angular pixel
scale. The result is cached in the WFS state and reused across measurements.
"""
@inline function sh_pixel_scale_init(d_subap::Real, padding::Int, src::AbstractSource)
    return lgs_pixel_scale(d_subap, padding, wavelength(src))
end

function prepare_sampling!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    sub = div(tel.params.resolution, wfs.params.n_subap)
    padding = wfs.params.diffraction_padding
    pixel_scale_req = wfs.params.pixel_scale
    d_subap = tel.params.diameter / wfs.params.n_subap
    pixel_scale_init = sh_pixel_scale_init(d_subap, padding, src)

    if pixel_scale_req !== nothing
        while pixel_scale_req / pixel_scale_init < 0.95
            padding += 1
            pixel_scale_init = sh_pixel_scale_init(d_subap, padding, src)
        end
    end

    binning_pixel_scale = if pixel_scale_req === nothing
        wfs.params.shannon_sampling ? 1 : 2
    else
        factor = pixel_scale_req / pixel_scale_init
        lower = max(1, floor(Int, factor))
        upper = max(1, ceil(Int, factor))
        abs(lower * pixel_scale_init - pixel_scale_req) <= abs(upper * pixel_scale_init - pixel_scale_req) ? lower : upper
    end

    pad = sub * padding
    while pad % binning_pixel_scale != 0
        padding += 1
        pad = sub * padding
        pixel_scale_init = sh_pixel_scale_init(d_subap, padding, src)
        if pixel_scale_req !== nothing
            factor = pixel_scale_req / pixel_scale_init
            lower = max(1, floor(Int, factor))
            upper = max(1, ceil(Int, factor))
            binning_pixel_scale = abs(lower * pixel_scale_init - pixel_scale_req) <= abs(upper * pixel_scale_init - pixel_scale_req) ? lower : upper
        end
    end

    n_pix_subap = wfs.params.n_pix_subap === nothing ? sub : wfs.params.n_pix_subap
    if isodd(n_pix_subap)
        throw(InvalidConfiguration("n_pix_subap must be even"))
    end

    if padding != wfs.state.effective_padding || pad != size(wfs.state.field, 1)
        ensure_sh_buffers!(wfs, pad)
        wfs.state.lgs_kernel_fft = similar(wfs.state.fft_buffer, Complex{eltype(wfs.state.fft_buffer)}, 0, 0, 0)
        wfs.state.lgs_kernel_tag = UInt(0)
        wfs.state.effective_padding = padding
    end

    if n_pix_subap != wfs.state.sampled_n_pix_subap
        wfs.state.spot = similar(wfs.state.spot, n_pix_subap, n_pix_subap)
        wfs.state.sampled_spot_cube = similar(wfs.state.sampled_spot_cube, eltype(wfs.state.sampled_spot_cube),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.spot_cube = similar(wfs.state.spot_cube, eltype(wfs.state.spot_cube),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.exported_spot_cube = similar(wfs.state.exported_spot_cube, eltype(wfs.state.exported_spot_cube),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.spot_cube_accum = similar(wfs.state.spot_cube_accum, eltype(wfs.state.spot_cube_accum),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.detector_noise_cube = similar(wfs.state.detector_noise_cube, eltype(wfs.state.detector_noise_cube),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.sampled_n_pix_subap = n_pix_subap
    end

    wfs.state.binning_pixel_scale = binning_pixel_scale
    T = eltype(wfs.state.slopes)
    half_shift_ratio = wfs.params.half_pixel_shift ? T(binning_pixel_scale) : zero(T)
    build_sh_phasor!(wfs, half_shift_ratio)
    return wfs
end

@inline function sh_exported_spot_cube(wfs::ShackHartmann)
    return wfs.state.exported_spot_cube
end

@inline function sh_sampled_spot_cube(wfs::ShackHartmann)
    return wfs.state.sampled_spot_cube
end

@inline function sync_signal_spots_from_sampled!(wfs::ShackHartmann)
    copyto!(wfs.state.spot_cube, wfs.state.sampled_spot_cube)
    return wfs.state.spot_cube
end

@inline function capture_sampled_spot_stack!(wfs::ShackHartmann, det::AbstractDetector, rng::AbstractRNG)
    copyto!(wfs.state.spot_cube, wfs.state.sampled_spot_cube)
    capture_stack!(det, wfs.state.spot_cube, wfs.state.detector_noise_cube, rng)
    return wfs.state.spot_cube
end

@inline function sync_exported_spots!(wfs::ShackHartmann)
    wfs.state.export_pixels_enabled || return wfs.state.exported_spot_cube
    copyto!(wfs.state.exported_spot_cube, wfs.state.spot_cube)
    return wfs.state.exported_spot_cube
end

"""
    sample_spot!(wfs, intensity)

Convert an oversampled diffractive lenslet intensity into the cropped lenslet
spot used for centroiding.

The operation is:

1. optional detector-style binning on the FFT intensity grid
2. centered crop/resize onto the configured `n_pix_subap x n_pix_subap` spot
"""
