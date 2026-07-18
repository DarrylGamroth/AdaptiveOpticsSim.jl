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
Immutable model and numerical-sampling parameters for a square regular
microlens array.

`n_lenslets` is the number of lenslets along each pupil axis, and
`pixel_scale_arcsec` is the requested detector-plane angular sampling in
arcseconds per pixel. The associated `SubapertureLayout` supplies the physical
pitch. This model does not currently describe focal length, fill factor,
per-lenslet prescriptions, or manufacturing errors.
"""
struct MicrolensArrayParams{T<:AbstractFloat}
    n_lenslets::Int
    half_pixel_shift::Bool
    diffraction_padding::Int
    pixel_scale_arcsec::Union{T,Nothing}
    n_pix_subap::Union{Int,Nothing}
    shannon_sampling::Bool
end

"""
A regular angular-coordinate microlens-array model, independent of detector
acquisition and WFS estimation.
"""
struct MicrolensArray{P<:MicrolensArrayParams}
    params::P
end

function MicrolensArray(; n_lenslets::Integer,
    half_pixel_shift::Bool=false,
    diffraction_padding::Integer=2, pixel_scale_arcsec=nothing,
    n_pix_subap::Union{Nothing,Integer}=nothing,
    shannon_sampling::Bool=true,
    T::Type{<:AbstractFloat}=pixel_scale_arcsec === nothing ? Float64 :
        typeof(float(pixel_scale_arcsec)))
    lenslets = Int(n_lenslets)
    padding = Int(diffraction_padding)
    pixels = n_pix_subap === nothing ? nothing : Int(n_pix_subap)
    lenslets > 0 || throw(InvalidConfiguration(
        "n_lenslets must be positive"))
    padding > 0 || throw(InvalidConfiguration(
        "diffraction_padding must be positive"))
    if pixels !== nothing
        pixels > 0 || throw(InvalidConfiguration(
            "n_pix_subap must be positive"))
        iseven(pixels) || throw(InvalidConfiguration(
            "n_pix_subap must be even"))
    end
    scale_arcsec = pixel_scale_arcsec === nothing ? nothing :
        T(pixel_scale_arcsec)
    if scale_arcsec !== nothing
        isfinite(scale_arcsec) && scale_arcsec > zero(T) ||
            throw(InvalidConfiguration(
                "pixel_scale_arcsec must be finite and positive"))
    end
    params = MicrolensArrayParams{T}(lenslets, half_pixel_shift,
        padding, scale_arcsec, pixels, shannon_sampling)
    return MicrolensArray{typeof(params)}(params)
end

@inline microlens_numeric_type(
    ::MicrolensArray{<:MicrolensArrayParams{T}}) where {T} = T

"""Column-major linear index of lenslet `(i, j)` in an `n`-by-`n` array."""
@inline sh_lenslet_index(i::Integer, j::Integer, n::Integer) =
    i + (j - 1) * n

struct ShackHartmannWFSParams{T<:AbstractFloat,VP<:AbstractValidSubaperturePolicy}
    threshold_convolution::T
    valid_subaperture_policy::VP
end

"""
Backend-bound plans and preallocated scratch for microlens propagation.

This is prepared execution state, not part of the physical `MicrolensArray`
description. One optical front end owns one instance under a single-writer
contract.
"""
mutable struct PreparedMicrolensPropagation{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    CC<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RC3<:AbstractArray{T,3},
    RT3<:AbstractArray{T,3},
    RB<:AbstractMatrix{T},
    RS<:AbstractMatrix{T},
    RCS<:AbstractArray{T,3},
    RA<:AbstractArray{T,3},
    P,
    PS,
    Pi,
    PSi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3},
    V<:AbstractVector{T}}
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
    spot_cube_accum::RA
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
    fft_asterism_stack::CC
    fft_asterism_plan::PS
    asterism_capacity::Int
    amp_scales::V
    amp_scales_host::Vector{T}
    opd_to_cycles::V
    opd_to_cycles_host::Vector{T}
end

"""
An independently composable Shack-Hartmann optical front end.

Only microlens propagation and subaperture-layout state are retained here;
detector acquisition, calibration, and estimation remain separate stages.
"""
struct ShackHartmannOpticalFrontEnd{M,PR,L,T<:AbstractFloat,S}
    microlens_array::M
    propagation::PR
    layout::L
    threshold_convolution::T
    source::S
end

"""
The microlens geometry and subaperture layout used by the direct geometric
Shack-Hartmann measurement path. It intentionally owns no propagation or
detector-acquisition storage.
"""
struct ShackHartmannDirectFrontEnd{M,L}
    microlens_array::M
    layout::L
end

"""Mutable detector-facing observation storage for a Shack-Hartmann sensor."""
mutable struct ShackHartmannAcquisitionState{T<:AbstractFloat,
    RC<:AbstractArray{T,3},RCE<:AbstractArray{T,3},RN<:AbstractArray{T,3}}
    spot_cube::RC
    exported_spot_cube::RCE
    detector_noise_cube::RN
    export_pixels_enabled::Bool
end

"""Mutable slope-estimation storage, separate from optics and acquisition."""
mutable struct ShackHartmannEstimatorState{T<:AbstractFloat,
    V<:AbstractVector{T}}
    slopes::V
    spot_stats::V
    spot_stats_accum::V
    slopes_host::Vector{T}
    centroid_host::Matrix{T}
end

"""
    ShackHartmannWFS

Shack-Hartmann wavefront sensor composed from microlens optics, subaperture
layout and calibration, detector acquisition, and slope estimation. The
geometric mode uses the same layout but follows a direct-measurement path.
"""
struct ShackHartmannWFS{M<:SensingMode,P<:ShackHartmannWFSParams,
    F,AQ,E,C,B<:AbstractArrayBackend} <: AbstractWFS
    params::P
    front_end::F
    acquisition::AQ
    estimator::E
    calibration::C
end

@inline backend(::ShackHartmannWFS{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

function ShackHartmannOpticalFrontEnd(
    microlens_array::MicrolensArray,
    propagation::PreparedMicrolensPropagation,
    layout::SubapertureLayout,
    source=nothing;
    threshold_convolution::Real=0.05)
    layout.n_subap == microlens_array.params.n_lenslets ||
        throw(InvalidConfiguration(
            "microlens array and subaperture layout counts differ"))
    size(propagation.fft_stack, 3) == layout.n_subap^2 ||
        throw(InvalidConfiguration(
            "prepared microlens propagation does not match the layout"))
    T = eltype(propagation.intensity)
    threshold = T(threshold_convolution)
    isfinite(threshold) && zero(T) <= threshold <= one(T) ||
        throw(InvalidConfiguration(
            "threshold_convolution must lie in [0, 1]"))
    return ShackHartmannOpticalFrontEnd{
        typeof(microlens_array),typeof(propagation),typeof(layout),T,
        typeof(source),
    }(microlens_array, propagation, layout, threshold, source)
end

function ShackHartmannOpticalFrontEnd(
    front_end::ShackHartmannOpticalFrontEnd, source)
    return ShackHartmannOpticalFrontEnd(front_end.microlens_array,
        front_end.propagation, front_end.layout, source;
        threshold_convolution=front_end.threshold_convolution)
end

@inline backend(front_end::ShackHartmannOpticalFrontEnd) =
    backend(front_end.propagation.fft_stack)
@inline microlens_array(front_end::ShackHartmannOpticalFrontEnd) =
    front_end.microlens_array
@inline n_lenslets(front_end::ShackHartmannOpticalFrontEnd) =
    microlens_array(front_end).params.n_lenslets
@inline subaperture_layout(front_end::ShackHartmannOpticalFrontEnd) =
    front_end.layout
@inline backend(front_end::ShackHartmannDirectFrontEnd) =
    backend(front_end.layout.valid_mask)
@inline microlens_array(front_end::ShackHartmannDirectFrontEnd) =
    front_end.microlens_array
@inline n_lenslets(front_end::ShackHartmannDirectFrontEnd) =
    microlens_array(front_end).params.n_lenslets
@inline subaperture_layout(front_end::ShackHartmannDirectFrontEnd) =
    front_end.layout
@inline sh_threshold_convolution(front_end::ShackHartmannOpticalFrontEnd) =
    front_end.threshold_convolution

"""
    ShackHartmannWFS(tel; ...)

Construct a Shack-Hartmann WFS on the telescope pupil grid.

Important diffractive quantities:

- `diffraction_padding` controls the focal-plane FFT grid size
- `pixel_scale_arcsec` and `shannon_sampling` determine detector-plane sampling
- `n_pix_subap` controls the cropped spot size used for centroiding

Diffractive construction prepares the lenslet FFT and acquisition workspaces.
Geometric construction intentionally prepares neither: it retains only the
subaperture layout, calibration, and direct slope-estimation storage.
"""
function ShackHartmannWFS(tel::Telescope; n_lenslets::Int, threshold::Real=0.1,
    threshold_cog::Real=0.01, threshold_convolution::Real=0.05, half_pixel_shift::Bool=false,
    diffraction_padding::Int=2, pixel_scale_arcsec=nothing,
    n_pix_subap=nothing,
    shannon_sampling::Bool=true,
    slope_extraction::Union{Nothing,AbstractSlopeExtractionModel}=nothing,
    valid_subaperture_policy::Union{Nothing,AbstractValidSubaperturePolicy}=nothing,
    mode::SensingMode=Geometric(),
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    n_lenslets > 0 || throw(InvalidConfiguration(
        "n_lenslets must be positive"))
    if tel.params.resolution % n_lenslets != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_lenslets"))
    end
    raw_policy = isnothing(valid_subaperture_policy) ? GeometryValidSubapertures(threshold=threshold, T=T) : valid_subaperture_policy
    policy = convert_valid_subaperture_policy(raw_policy, T)
    convolution_threshold = T(threshold_convolution)
    isfinite(convolution_threshold) &&
        zero(T) <= convolution_threshold <= one(T) ||
        throw(InvalidConfiguration(
            "threshold_convolution must lie in [0, 1]"))
    params = ShackHartmannWFSParams{T,typeof(policy)}(
        convolution_threshold, policy)
    mla = MicrolensArray(; n_lenslets, half_pixel_shift,
        diffraction_padding, pixel_scale_arcsec, n_pix_subap,
        shannon_sampling, T)
    valid_mask = backend{Bool}(undef, n_lenslets, n_lenslets)
    fill!(valid_mask, false)
    slopes = backend{T}(undef, 2 * n_lenslets * n_lenslets)
    fill!(slopes, zero(T))
    sub = div(tel.params.resolution, n_lenslets)
    spot_stats = backend{T}(undef, 3 * n_lenslets * n_lenslets)
    spot_stats_accum = backend{T}(undef, 3 * n_lenslets * n_lenslets)
    valid_mask_host = Matrix{Bool}(undef, n_lenslets, n_lenslets)
    fill!(valid_mask_host, false)
    reference_signal_host = Vector{T}(undef, 2 * n_lenslets * n_lenslets)
    fill!(reference_signal_host, zero(T))
    slopes_host = Vector{T}(undef, 2 * n_lenslets * n_lenslets)
    centroid_host = Matrix{T}(undef, sub, sub)
    # Columns are the first and second pupil-axis centroid components. The
    # lenslet row uses Julia's column-major linear order for the n-by-n mask.
    reference_signal_2d = backend{T}(undef, n_lenslets * n_lenslets, 2)
    fill!(reference_signal_2d, zero(T))
    layout = SubapertureLayout(n_lenslets, tel.params.resolution, tel.params.diameter, threshold, valid_mask, valid_mask_host)
    extraction = prepare_sh_slope_extraction(slope_extraction,
        threshold_cog, T)
    calibration = SubapertureCalibration(reference_signal_2d,
        reference_signal_host, extraction)
    estimator = ShackHartmannEstimatorState(slopes, spot_stats,
        spot_stats_accum, slopes_host, centroid_host)
    front_end, acquisition = _prepare_sh_mode_storage(mode, backend,
        T, mla, layout, sub, convolution_threshold)
    wfs = ShackHartmannWFS{
        typeof(mode),typeof(params),typeof(front_end),typeof(acquisition),
        typeof(estimator),typeof(calibration),typeof(selector),
    }(params, front_end, acquisition, estimator, calibration)
    update_valid_mask!(wfs, tel)
    return wfs
end

@inline prepare_sh_slope_extraction(::Nothing, threshold::Real,
    ::Type{T}) where {T<:AbstractFloat} =
    CenterOfGravityExtraction(T(threshold); T=T)

function prepare_sh_slope_extraction(model::CenterOfGravityExtraction,
    ::Real, ::Type{T}) where {T<:AbstractFloat}
    window = model.window === nothing ? nothing : Matrix{T}(model.window)
    return CenterOfGravityExtraction(T(model.threshold); window, T=T)
end

@inline prepare_sh_slope_extraction(model::AbstractSlopeExtractionModel,
    ::Real, ::Type{T}) where {T<:AbstractFloat} = model

@inline _prepare_sh_mode_storage(::Geometric, backend,
    ::Type{T}, mla::MicrolensArray, layout::SubapertureLayout, sub::Int,
    threshold_convolution) where {T<:AbstractFloat} =
    (ShackHartmannDirectFrontEnd(mla, layout), nothing)

function _prepare_sh_mode_storage(::Diffractive, backend,
    ::Type{T}, mla::MicrolensArray, layout::SubapertureLayout, sub::Int,
    threshold_convolution) where {T<:AbstractFloat}
    propagation = _prepare_microlens_propagation(backend, T, mla, sub)
    spot_cube = similar(propagation.sampled_spot_cube)
    acquisition = ShackHartmannAcquisitionState(spot_cube,
        similar(spot_cube), similar(spot_cube), true)
    front_end = ShackHartmannOpticalFrontEnd(mla, propagation, layout,
        nothing; threshold_convolution)
    return front_end, acquisition
end

"""
    prepare_microlens_propagation(microlens_array, pupil_resolution;
        T=Float64, backend=CPUBackend())

Prepare backend-specific FFT plans and reusable scratch storage for an
independent `MicrolensArray`. This cold-path constructor performs no detector
acquisition or slope-estimator setup.
"""
function prepare_microlens_propagation(mla::MicrolensArray,
    pupil_resolution::Integer;
    T::Type{<:AbstractFloat}=microlens_numeric_type(mla),
    backend::AbstractArrayBackend=CPUBackend())
    resolution = Int(pupil_resolution)
    resolution > 0 || throw(InvalidConfiguration(
        "pupil_resolution must be positive"))
    resolution % mla.params.n_lenslets == 0 ||
        throw(InvalidConfiguration(
            "pupil_resolution must be divisible by n_lenslets"))
    selector = _resolve_backend_selector(backend)
    array_type = _resolve_array_backend(selector)
    return _prepare_microlens_propagation(array_type, T, mla,
        div(resolution, mla.params.n_lenslets))
end

function _prepare_microlens_propagation(backend, ::Type{T},
    mla::MicrolensArray, sub::Int) where {T<:AbstractFloat}
    n_lenslets = mla.params.n_lenslets
    pad = max(sub, sub * mla.params.diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    phasor = similar(field)
    fft_buffer = similar(field)
    fft_stack = backend{Complex{T}}(undef, pad, pad,
        n_lenslets * n_lenslets)
    intensity = backend{T}(undef, pad, pad)
    intensity_stack = backend{T}(undef, pad, pad,
        n_lenslets * n_lenslets)
    intensity_tmp_stack = similar(intensity_stack)
    temp = similar(intensity)
    bin_buffer = backend{T}(undef, sub, sub)
    spot = similar(bin_buffer)
    sampled_spot_cube = backend{T}(undef, n_lenslets * n_lenslets,
        sub, sub)
    spot_cube_accum = similar(sampled_spot_cube)
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
    propagation = PreparedMicrolensPropagation(field, phasor, fft_buffer,
        fft_stack, intensity, intensity_stack, intensity_tmp_stack, temp,
        bin_buffer, spot, sampled_spot_cube, spot_cube_accum, fft_plan,
        fft_stack_plan, ifft_plan, ifft_stack_plan, elongation_kernel,
        lgs_kernel_fft, UInt(0), mla.params.diffraction_padding, 1, sub,
        T(NaN), fft_asterism_stack, fft_asterism_plan, 1, amp_scales,
        amp_scales_host, opd_to_cycles, opd_to_cycles_host)
    return propagation
end

sensing_mode(::ShackHartmannWFS{M}) where {M} = M()
@inline n_lenslets(wfs::ShackHartmannWFS) =
    n_lenslets(wfs.front_end)
@inline subaperture_calibration(wfs::ShackHartmannWFS) = wfs.calibration
@inline valid_subaperture_policy(wfs::ShackHartmannWFS) = wfs.params.valid_subaperture_policy
@inline slope_extraction_model(wfs::ShackHartmannWFS) = slope_extraction_model(subaperture_calibration(wfs))
@inline centroid_threshold(wfs::ShackHartmannWFS) = slope_extraction_model(wfs).threshold

@inline function sh_common_spectral_grid_wavelength(
    wfs::ShackHartmannWFS, src::SpectralSource)
    T = eltype(wfs.estimator.slopes)
    samples = spectral_bundle(src).samples
    isempty(samples) && return (false, zero(T))
    wavelength_ref = T(first(samples).wavelength)
    isfinite(wavelength_ref) && wavelength_ref > zero(T) ||
        return (false, wavelength_ref)
    @inbounds for i in 2:length(samples)
        wavelength_i = T(samples[i].wavelength)
        if !isfinite(wavelength_i) || wavelength_i <= zero(T) ||
                wavelength_i != wavelength_ref
            return (false, wavelength_ref)
        end
    end
    return (true, wavelength_ref)
end

@inline function sh_has_common_spectral_grid(
    wfs::ShackHartmannWFS, src::SpectralSource)
    compatible, _ = sh_common_spectral_grid_wavelength(wfs, src)
    return compatible
end

function require_sh_common_spectral_grid(
    wfs::ShackHartmannWFS, src::SpectralSource)
    compatible, wavelength_ref = sh_common_spectral_grid_wavelength(wfs, src)
    compatible || throw(InvalidConfiguration(
        "diffractive ShackHartmannWFS spectral samples must share one " *
        "finite, positive wavelength on the WFS numerical grid; distinct " *
        "wavelengths require an explicit native-to-detector sampling map"))
    return wavelength_ref
end

abstract type AbstractShackHartmannWFSSensingPlan end
struct ShackHartmannWFSScalarPlan <: AbstractShackHartmannWFSSensingPlan end
struct ShackHartmannWFSBatchedPlan <: AbstractShackHartmannWFSSensingPlan end
struct ShackHartmannWFSDeviceStatsPlan <: AbstractShackHartmannWFSSensingPlan end
struct ShackHartmannWFSRocmSafePlan <: AbstractShackHartmannWFSSensingPlan end
struct ShackHartmannWFSRocmHostStatsPlan <: AbstractShackHartmannWFSSensingPlan end

@inline sh_sensing_execution_plan(style::ExecutionStyle, wfs::ShackHartmannWFS) =
    sh_sensing_execution_plan(typeof(style), typeof(wfs))
@inline sh_sensing_execution_plan(::Type{<:ScalarCPUStyle}, ::Type{<:ShackHartmannWFS}) = ShackHartmannWFSScalarPlan()
@inline sh_sensing_execution_plan(::Type{<:AcceleratorStyle}, ::Type{<:ShackHartmannWFS}) = ShackHartmannWFSBatchedPlan()
@inline sh_sensing_execution_plan(wfs::ShackHartmannWFS) =
    sh_sensing_execution_plan(execution_style(wfs.estimator.slopes), wfs)

@inline sh_uses_rocm_safe_sensing_plan(::AbstractShackHartmannWFSSensingPlan) = false
@inline sh_uses_rocm_safe_sensing_plan(::ShackHartmannWFSRocmSafePlan) = true
@inline sh_uses_rocm_safe_sensing_plan(wfs::ShackHartmannWFS) =
    sh_uses_rocm_safe_sensing_plan(sh_sensing_execution_plan(wfs))
@inline sh_uses_host_stats_sensing_plan(::AbstractShackHartmannWFSSensingPlan) = false
@inline sh_uses_host_stats_sensing_plan(::ShackHartmannWFSRocmSafePlan) = true
@inline sh_uses_host_stats_sensing_plan(::ShackHartmannWFSRocmHostStatsPlan) = true
@inline sh_uses_host_stats_sensing_plan(wfs::ShackHartmannWFS) =
    sh_uses_host_stats_sensing_plan(sh_sensing_execution_plan(wfs))
@inline sh_uses_batched_sensing_plan(::AbstractShackHartmannWFSSensingPlan) = false
@inline sh_uses_batched_sensing_plan(::ShackHartmannWFSBatchedPlan) = true
@inline sh_uses_batched_sensing_plan(::ShackHartmannWFSDeviceStatsPlan) = true
@inline sh_uses_batched_sensing_plan(wfs::ShackHartmannWFS) =
    sh_uses_batched_sensing_plan(sh_sensing_execution_plan(wfs))
@inline sh_uses_device_stats_sensing_plan(::AbstractShackHartmannWFSSensingPlan) = false
@inline sh_uses_device_stats_sensing_plan(::ShackHartmannWFSDeviceStatsPlan) = true
@inline sh_uses_device_stats_sensing_plan(wfs::ShackHartmannWFS) =
    sh_uses_device_stats_sensing_plan(sh_sensing_execution_plan(wfs))

convert_valid_subaperture_policy(policy::GeometryValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    GeometryValidSubapertures(threshold=T(policy.threshold), T=T)

convert_valid_subaperture_policy(policy::FluxThresholdValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    FluxThresholdValidSubapertures(light_ratio=T(policy.light_ratio), T=T)

@inline function sh_safe_peak_value(A::AbstractArray{T}) where {T<:AbstractFloat}
    return backend_maximum_value(A)
end

@inline function sh_refresh_valid_mask_host!(wfs::ShackHartmannWFS)
    layout = wfs.front_end.layout
    copyto!(layout.valid_mask_host, layout.valid_mask)
    return layout.valid_mask_host
end

function update_valid_mask!(wfs::ShackHartmannWFS, tel::Telescope)
    update_valid_mask!(wfs, tel, valid_subaperture_policy(wfs))
    return wfs
end

function update_valid_mask!(wfs::ShackHartmannWFS, tel::Telescope, policy::GeometryValidSubapertures)
    update_subaperture_layout!(wfs.front_end.layout, pupil_mask(tel), policy)
    invalidate_sh_calibration!(wfs)
    return wfs
end

function update_valid_mask!(wfs::ShackHartmannWFS, tel::Telescope, policy::FluxThresholdValidSubapertures)
    update_subaperture_layout!(wfs.front_end.layout, pupil_reflectivity(tel), policy)
    invalidate_sh_calibration!(wfs)
    return wfs
end

function update_valid_mask!(::ShackHartmannWFS, ::Telescope, policy)
    throw(UnsupportedAlgorithm("unsupported valid subaperture policy $(typeof(policy))"))
end

@inline function sh_grouped_stack_capacity(
    front_end::ShackHartmannOpticalFrontEnd)
    propagation = front_end.propagation
    return n_lenslets(front_end)^2 * propagation.asterism_capacity
end

function ensure_sh_buffers!(front_end::ShackHartmannOpticalFrontEnd,
    pad::Int)
    propagation = front_end.propagation
    n_spots = n_lenslets(front_end)^2
    if size(propagation.field) != (pad, pad)
        propagation.field = similar(propagation.field, pad, pad)
        propagation.phasor = similar(propagation.phasor, pad, pad)
        propagation.fft_buffer = similar(propagation.fft_buffer, pad, pad)
        propagation.fft_stack = similar(propagation.fft_stack,
            eltype(propagation.field), pad, pad, n_spots)
        propagation.intensity = similar(propagation.intensity, pad, pad)
        propagation.intensity_stack = similar(propagation.intensity_stack,
            eltype(propagation.intensity), pad, pad, n_spots)
        total = sh_grouped_stack_capacity(front_end)
        propagation.intensity_tmp_stack = similar(
            propagation.intensity_tmp_stack, eltype(propagation.intensity),
            pad, pad, total)
        propagation.temp = similar(propagation.temp, pad, pad)
        propagation.fft_plan = plan_fft_backend!(propagation.fft_buffer)
        propagation.fft_stack_plan = plan_fft_backend!(
            propagation.fft_stack, (1, 2))
        propagation.ifft_plan = plan_ifft_backend!(propagation.fft_buffer)
        propagation.ifft_stack_plan = plan_ifft_backend!(
            propagation.fft_stack, (1, 2))
        propagation.fft_asterism_stack = similar(
            propagation.fft_asterism_stack, pad, pad, total)
        propagation.fft_asterism_plan = plan_fft_backend!(
            propagation.fft_asterism_stack, (1, 2))
        propagation.phasor_ratio = eltype(propagation.intensity)(NaN)
    end
    return front_end
end

@inline function invalidate_sh_calibration!(wfs::ShackHartmannWFS)
    wfs.calibration.calibrated = false
    wfs.calibration.revision += UInt(1)
    return nothing
end

function ensure_sh_asterism_buffers!(
    front_end::ShackHartmannOpticalFrontEnd, n_sources::Int)
    propagation = front_end.propagation
    n_sources > 0 || throw(InvalidConfiguration("asterism must contain at least one source"))
    if n_sources > propagation.asterism_capacity
        pad = size(propagation.fft_stack, 1)
        n_spots = n_lenslets(front_end)^2
        total = n_spots * n_sources
        propagation.fft_asterism_stack = similar(
            propagation.fft_asterism_stack, pad, pad, total)
        propagation.intensity_tmp_stack = similar(
            propagation.intensity_tmp_stack,
            eltype(propagation.intensity_tmp_stack), pad, pad, total)
        propagation.fft_asterism_plan = plan_fft_backend!(
            propagation.fft_asterism_stack, (1, 2))
        propagation.asterism_capacity = n_sources
    end
    if length(propagation.amp_scales) != n_sources
        propagation.amp_scales = similar(propagation.amp_scales, n_sources)
    end
    if length(propagation.amp_scales_host) != n_sources
        propagation.amp_scales_host = Vector{
            eltype(propagation.amp_scales_host)}(undef, n_sources)
    end
    if length(propagation.opd_to_cycles) != n_sources
        propagation.opd_to_cycles = similar(
            propagation.opd_to_cycles, n_sources)
    end
    if length(propagation.opd_to_cycles_host) != n_sources
        propagation.opd_to_cycles_host = Vector{
            eltype(propagation.opd_to_cycles_host)}(undef, n_sources)
    end
    return front_end
end

function build_sh_phasor!(front_end::ShackHartmannOpticalFrontEnd,
    ratio::T) where {T<:AbstractFloat}
    propagation = front_end.propagation
    if size(propagation.phasor, 1) == 0
        return front_end
    end
    if isequal(propagation.phasor_ratio, ratio)
        return front_end
    end
    n = size(propagation.phasor, 1)
    scale = -T(π) * (T(n) + one(T) + ratio) / T(n)
    host = Matrix{Complex{T}}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        host[i, j] = cis(scale * (i + j - 2))
    end
    copyto!(propagation.phasor, host)
    propagation.phasor_ratio = ratio
    return front_end
end

"""
    prepare_sampling!(wfs, tel, src)

Resolve the diffractive lenslet sampling for a source/wavelength pair.

This chooses the effective FFT padding, detector binning, and cropped spot size
so that the propagated spot grid is consistent with the requested angular pixel
scale. The result is cached in prepared propagation state and reused across
measurements.
"""
@inline sh_pixel_scale_init(d_subap::Real, padding::Int,
    wavelength_m::Real) = lgs_pixel_scale(d_subap, padding, wavelength_m)

@inline sh_pixel_scale_init(d_subap::Real, padding::Int,
    src::AbstractSource) = sh_pixel_scale_init(d_subap, padding,
    wavelength(src))

function prepare_sampling!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    return prepare_sampling_wavelength!(wfs, tel, wavelength(src))
end

function prepare_sampling!(wfs::ShackHartmannWFS, tel::Telescope,
    src::SpectralSource)
    wavelength_ref = require_sh_common_spectral_grid(wfs, src)
    return prepare_sampling_wavelength!(wfs, tel, wavelength_ref)
end

function prepare_sampling_wavelength!(wfs::ShackHartmannWFS,
    tel::Telescope, wavelength_m::Real)
    return prepare_sampling_wavelength!(wfs, tel.params.resolution,
        tel.params.diameter, wavelength_m)
end

function prepare_sampling_wavelength!(wfs::ShackHartmannWFS,
    pupil_resolution::Int, pupil_diameter_m::Real, wavelength_m::Real)
    front_end = wfs.front_end
    propagation = front_end.propagation
    sampling_before = (size(propagation.field),
        propagation.effective_padding, propagation.binning_pixel_scale,
        propagation.sampled_n_pix_subap, propagation.phasor_ratio)
    _prepare_microlens_sampling_wavelength!(front_end, pupil_resolution,
        pupil_diameter_m, wavelength_m)
    sampling_after = (size(propagation.field),
        propagation.effective_padding, propagation.binning_pixel_scale,
        propagation.sampled_n_pix_subap, propagation.phasor_ratio)
    sampling_after == sampling_before || invalidate_sh_calibration!(wfs)
    ensure_sh_acquisition_buffers!(wfs,
        propagation.sampled_n_pix_subap)
    return wfs
end

function _prepare_microlens_sampling_wavelength!(
    front_end::ShackHartmannOpticalFrontEnd,
    pupil_resolution::Int, pupil_diameter_m::Real, wavelength_m::Real)
    propagation = front_end.propagation
    pupil_resolution % n_lenslets(front_end) == 0 ||
        throw(InvalidConfiguration(
            "pupil resolution must be divisible by n_lenslets"))
    sub = div(pupil_resolution, n_lenslets(front_end))
    padding = front_end.microlens_array.params.diffraction_padding
    pixel_scale_req = front_end.microlens_array.params.pixel_scale_arcsec
    d_subap = pupil_diameter_m / n_lenslets(front_end)
    pixel_scale_init = sh_pixel_scale_init(d_subap, padding, wavelength_m)

    if pixel_scale_req !== nothing
        while pixel_scale_req / pixel_scale_init < 0.95
            padding += 1
            pixel_scale_init = sh_pixel_scale_init(d_subap, padding,
                wavelength_m)
        end
    end

    binning_pixel_scale = if pixel_scale_req === nothing
        front_end.microlens_array.params.shannon_sampling ? 1 : 2
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
        pixel_scale_init = sh_pixel_scale_init(d_subap, padding,
            wavelength_m)
        if pixel_scale_req !== nothing
            factor = pixel_scale_req / pixel_scale_init
            lower = max(1, floor(Int, factor))
            upper = max(1, ceil(Int, factor))
            binning_pixel_scale = abs(lower * pixel_scale_init - pixel_scale_req) <= abs(upper * pixel_scale_init - pixel_scale_req) ? lower : upper
        end
    end

    n_pix_subap = front_end.microlens_array.params.n_pix_subap === nothing ?
        sub : front_end.microlens_array.params.n_pix_subap
    if isodd(n_pix_subap)
        throw(InvalidConfiguration("n_pix_subap must be even"))
    end

    if padding != propagation.effective_padding ||
            pad != size(propagation.field, 1)
        ensure_sh_buffers!(front_end, pad)
        propagation.lgs_kernel_fft = similar(propagation.fft_buffer,
            eltype(propagation.fft_buffer), 0, 0, 0)
        propagation.lgs_kernel_tag = UInt(0)
        propagation.effective_padding = padding
    end

    if n_pix_subap != propagation.sampled_n_pix_subap
        propagation.spot = similar(propagation.spot,
            n_pix_subap, n_pix_subap)
        propagation.sampled_spot_cube = similar(
            propagation.sampled_spot_cube,
            eltype(propagation.sampled_spot_cube),
            n_lenslets(front_end)^2, n_pix_subap, n_pix_subap)
        propagation.spot_cube_accum = similar(propagation.spot_cube_accum,
            eltype(propagation.spot_cube_accum),
            n_lenslets(front_end)^2, n_pix_subap, n_pix_subap)
        propagation.sampled_n_pix_subap = n_pix_subap
    end

    propagation.binning_pixel_scale = binning_pixel_scale
    T = eltype(propagation.intensity)
    half_shift_ratio = front_end.microlens_array.params.half_pixel_shift ?
        T(binning_pixel_scale) : zero(T)
    build_sh_phasor!(front_end, half_shift_ratio)
    return front_end
end


function ensure_sh_acquisition_buffers!(wfs::ShackHartmannWFS,
    n_pix_subap::Int)
    n_spots = n_lenslets(wfs) * n_lenslets(wfs)
    expected = (n_spots, n_pix_subap, n_pix_subap)
    size(wfs.acquisition.spot_cube) == expected && return wfs
    wfs.acquisition.spot_cube = similar(wfs.acquisition.spot_cube,
        eltype(wfs.acquisition.spot_cube), expected...)
    wfs.acquisition.exported_spot_cube = similar(
        wfs.acquisition.exported_spot_cube,
        eltype(wfs.acquisition.exported_spot_cube), expected...)
    wfs.acquisition.detector_noise_cube = similar(
        wfs.acquisition.detector_noise_cube,
        eltype(wfs.acquisition.detector_noise_cube), expected...)
    return wfs
end

@inline function sh_exported_spot_cube(wfs::ShackHartmannWFS)
    return wfs.acquisition.exported_spot_cube
end

@inline function sh_sampled_spot_cube(wfs::ShackHartmannWFS)
    return wfs.front_end.propagation.sampled_spot_cube
end

@inline function sync_signal_spots_from_sampled!(wfs::ShackHartmannWFS)
    copyto!(wfs.acquisition.spot_cube,
        wfs.front_end.propagation.sampled_spot_cube)
    return wfs.acquisition.spot_cube
end

@inline function capture_sampled_spot_stack!(wfs::ShackHartmannWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    copyto!(wfs.acquisition.spot_cube,
        wfs.front_end.propagation.sampled_spot_cube)
    capture_stack!(det, wfs.acquisition.spot_cube, wfs.acquisition.detector_noise_cube,
        src, rng)
    return wfs.acquisition.spot_cube
end

@inline function sync_exported_spots!(wfs::ShackHartmannWFS)
    wfs.acquisition.export_pixels_enabled || return wfs.acquisition.exported_spot_cube
    copyto!(wfs.acquisition.exported_spot_cube, wfs.acquisition.spot_cube)
    return wfs.acquisition.exported_spot_cube
end

"""
    sample_spot!(front_end, intensity)

Convert an oversampled diffractive lenslet intensity into the cropped lenslet
spot used for centroiding.

The operation is:

1. optional detector-style binning on the FFT intensity grid
2. centered crop/resize onto the configured `n_pix_subap x n_pix_subap` spot
"""
