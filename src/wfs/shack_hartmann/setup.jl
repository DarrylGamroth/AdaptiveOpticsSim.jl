#
# Shack-Hartmann wavefront sensing
#
# The physical path uses FFT-based Fraunhofer propagation on each lenslet.
# For LGS sensing, elongated spots are handled through focal-plane convolution.
# For asterisms and GPU execution, the implementation batches lenslet/source
# stacks so the algorithm stays mathematically the same while reducing launch
# and detector-processing overhead.
#
"""Column-major linear index of lenslet `(i, j)` in an `n`-by-`n` array."""
@inline sh_lenslet_index(i::Integer, j::Integer, n::Integer) =
    i + (j - 1) * n

struct ShackHartmannWFSParams{T<:AbstractFloat,VP<:AbstractValidSubaperturePolicy}
    threshold_convolution::T
    valid_subaperture_policy::VP
end

"""Physical and configuration definition for a diffractive SH front end."""
struct ShackHartmannOpticalFrontEnd{M,L,T<:AbstractFloat,S}
    microlens_array::M
    layout::L
    threshold_convolution::T
    source::S
end

"""
Execution composition of an SH front-end definition with replaceable
microlens execution resources. Prepared stage owners bind this model exactly.
"""
struct ShackHartmannOptics{F,PR}
    front_end::F
    propagation::PR
end

"""
    ShackHartmannWFS

Physical diffractive Shack-Hartmann front end. It owns microlens optics and
valid-subaperture geometry; detector acquisition and frame estimation are
explicitly composed by their respective owners.
"""
struct ShackHartmannWFS{P<:ShackHartmannWFSParams,F,O,
    B<:AbstractArrayBackend} <: AbstractWFS
    params::P
    front_end::F
    optics::O
end

@inline backend(::ShackHartmannWFS{P,F,O,B}) where {P,F,O,B} = B()

function ShackHartmannOpticalFrontEnd(
    microlens_array::MicrolensArray,
    layout::SubapertureLayout,
    source=nothing;
    threshold_convolution::Real=0.05)
    layout.n_subap == microlens_array.params.n_lenslets ||
        throw(InvalidConfiguration(
            "microlens array and subaperture layout counts differ"))
    T = microlens_numeric_type(microlens_array)
    threshold = T(threshold_convolution)
    isfinite(threshold) && zero(T) <= threshold <= one(T) ||
        throw(InvalidConfiguration(
            "threshold_convolution must lie in [0, 1]"))
    return ShackHartmannOpticalFrontEnd{
        typeof(microlens_array),typeof(layout),T,typeof(source),
    }(microlens_array, layout, threshold, source)
end

function ShackHartmannOpticalFrontEnd(
    front_end::ShackHartmannOpticalFrontEnd, source)
    return ShackHartmannOpticalFrontEnd(front_end.microlens_array,
        front_end.layout, source;
        threshold_convolution=front_end.threshold_convolution)
end

function ShackHartmannOptics(
    front_end::ShackHartmannOpticalFrontEnd,
    propagation::PreparedMicrolensPropagation,
)
    plan = microlens_propagation_plan(propagation)
    workspace = microlens_propagation_workspace(propagation)
    plan.microlens_array === front_end.microlens_array ||
        throw(InvalidConfiguration(
            "prepared microlens propagation has a different microlens definition"))
    plan.pupil_samples_per_lenslet == front_end.layout.subap_pixels ||
        throw(InvalidConfiguration(
            "prepared microlens propagation sampling does not match the layout"))
    size(workspace.fft_stack, 3) == front_end.layout.n_subap^2 ||
        throw(InvalidConfiguration(
            "prepared microlens propagation does not match the layout"))
    plan.numeric_type === eltype(workspace.intensity) &&
        eltype(workspace.intensity) ===
        microlens_numeric_type(front_end.microlens_array) ||
        throw(InvalidConfiguration(
            "prepared microlens propagation precision does not match the front end"))
    return ShackHartmannOptics{
        typeof(front_end),typeof(propagation)}(front_end, propagation)
end

function ShackHartmannOptics(
    model::ShackHartmannOptics, source)
    return ShackHartmannOptics(
        ShackHartmannOpticalFrontEnd(model.front_end, source),
        model.propagation)
end

@inline backend(model::ShackHartmannOptics) =
    backend(microlens_propagation_workspace(
        model.propagation).fft_stack)
@inline microlens_array(model::ShackHartmannOptics) =
    microlens_array(model.front_end)
@inline n_lenslets(model::ShackHartmannOptics) =
    n_lenslets(model.front_end)
@inline subaperture_layout(model::ShackHartmannOptics) =
    subaperture_layout(model.front_end)
@inline backend(front_end::ShackHartmannOpticalFrontEnd) =
    backend(front_end.layout.valid_mask)
@inline microlens_array(front_end::ShackHartmannOpticalFrontEnd) =
    front_end.microlens_array
@inline n_lenslets(front_end::ShackHartmannOpticalFrontEnd) =
    microlens_array(front_end).params.n_lenslets
@inline subaperture_layout(front_end::ShackHartmannOpticalFrontEnd) =
    front_end.layout
@inline sh_threshold_convolution(model::ShackHartmannOptics) =
    model.front_end.threshold_convolution

"""
    ShackHartmannWFS(tel; ...)

Construct a Shack-Hartmann WFS on the telescope pupil grid.

Important physical quantities:

- `diffraction_padding` controls the focal-plane FFT grid size
- `pixel_scale_arcsec` and `shannon_sampling` determine detector-plane sampling
- `n_pix_subap` controls the cropped focal-plane spot size
"""
function ShackHartmannWFS(tel::Telescope; n_lenslets::Int, threshold::Real=0.1,
    threshold_convolution::Real=0.05, half_pixel_shift::Bool=false,
    diffraction_padding::Int=2, pixel_scale_arcsec=nothing,
    n_pix_subap=nothing,
    shannon_sampling::Bool=true,
    valid_subaperture_policy::Union{Nothing,AbstractValidSubaperturePolicy}=nothing,
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
    sub = div(tel.params.resolution, n_lenslets)
    valid_mask_host = Matrix{Bool}(undef, n_lenslets, n_lenslets)
    fill!(valid_mask_host, false)
    layout = SubapertureLayout(n_lenslets, tel.params.resolution, tel.params.diameter, threshold, valid_mask, valid_mask_host)
    propagation = _prepare_microlens_propagation(backend, T, mla, sub)
    front_end = ShackHartmannOpticalFrontEnd(mla, layout, nothing;
        threshold_convolution=convolution_threshold)
    optics = ShackHartmannOptics(front_end, propagation)
    wfs = ShackHartmannWFS{
        typeof(params),typeof(front_end),typeof(optics),typeof(selector),
    }(params, front_end, optics)
    initialize_valid_mask!(wfs, tel, policy)
    return wfs
end

@inline shack_hartmann_optics(wfs::ShackHartmannWFS) = wfs.optics
@inline shack_hartmann_optics(wfs::ShackHartmannWFS, source) =
    ShackHartmannOptics(wfs.optics, source)
@inline n_lenslets(wfs::ShackHartmannWFS) =
    n_lenslets(wfs.front_end)
@inline valid_subaperture_policy(wfs::ShackHartmannWFS) = wfs.params.valid_subaperture_policy

@inline function sh_common_spectral_grid_wavelength(
    wfs::ShackHartmannWFS, src::SpectralSource)
    T = microlens_numeric_type(wfs.front_end.microlens_array)
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

convert_valid_subaperture_policy(policy::GeometryValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    GeometryValidSubapertures(threshold=T(policy.threshold), T=T)

convert_valid_subaperture_policy(policy::FluxThresholdValidSubapertures, ::Type{T}) where {T<:AbstractFloat} =
    FluxThresholdValidSubapertures(light_ratio=T(policy.light_ratio), T=T)

function initialize_valid_mask!(wfs::ShackHartmannWFS,
    tel::Telescope, policy::GeometryValidSubapertures)
    update_subaperture_layout!(wfs.front_end.layout, pupil_mask(tel),
        policy)
    return wfs
end

function initialize_valid_mask!(wfs::ShackHartmannWFS,
    tel::Telescope, policy::FluxThresholdValidSubapertures)
    update_subaperture_layout!(wfs.front_end.layout,
        pupil_reflectivity(tel), policy)
    return wfs
end

function update_valid_mask!(wfs::ShackHartmannWFS,
    pupil::PupilFunction)
    update_valid_mask!(wfs, pupil, valid_subaperture_policy(wfs))
    return wfs
end

function update_valid_mask!(wfs::ShackHartmannWFS,
    pupil::PupilFunction, policy::GeometryValidSubapertures)
    update_subaperture_layout!(wfs.front_end.layout, pupil.amplitude,
        policy)
    return wfs
end

function update_valid_mask!(wfs::ShackHartmannWFS,
    pupil::PupilFunction, policy::FluxThresholdValidSubapertures)
    update_subaperture_layout_from_amplitude!(wfs.front_end.layout,
        pupil.amplitude, policy)
    return wfs
end

function update_valid_mask!(::ShackHartmannWFS, ::PupilFunction, policy)
    throw(UnsupportedAlgorithm("unsupported valid subaperture policy $(typeof(policy))"))
end

@inline function sh_grouped_stack_capacity(
    optics::ShackHartmannOptics)
    propagation = microlens_propagation_workspace(optics.propagation)
    return n_lenslets(optics)^2 * propagation.asterism_capacity
end

function ensure_sh_buffers!(optics::ShackHartmannOptics,
    pad::Int)
    propagation = microlens_propagation_workspace(optics.propagation)
    n_spots = n_lenslets(optics)^2
    if size(propagation.field) != (pad, pad)
        propagation.field = similar(propagation.field, pad, pad)
        propagation.phasor = similar(propagation.phasor, pad, pad)
        propagation.fft_buffer = similar(propagation.fft_buffer, pad, pad)
        propagation.fft_stack = similar(propagation.fft_stack,
            eltype(propagation.field), pad, pad, n_spots)
        propagation.intensity = similar(propagation.intensity, pad, pad)
        propagation.intensity_stack = similar(propagation.intensity_stack,
            eltype(propagation.intensity), pad, pad, n_spots)
        total = sh_grouped_stack_capacity(optics)
        propagation.intensity_tmp_stack = similar(
            propagation.intensity_tmp_stack, eltype(propagation.intensity),
            pad, pad, total)
        propagation.temp = similar(propagation.temp, pad, pad)
        propagation.fft_plan = plan_repeated_fft_backend!(
            propagation.fft_buffer)
        propagation.fft_stack_plan = plan_repeated_fft_backend!(
            propagation.fft_stack, (1, 2))
        propagation.ifft_plan = plan_repeated_ifft_backend!(
            propagation.fft_buffer)
        propagation.ifft_stack_plan = plan_repeated_ifft_backend!(
            propagation.fft_stack, (1, 2))
        propagation.fft_asterism_stack = similar(
            propagation.fft_asterism_stack, pad, pad, total)
        propagation.fft_asterism_plan = plan_repeated_fft_backend!(
            propagation.fft_asterism_stack, (1, 2))
        propagation.phasor_ratio = eltype(propagation.intensity)(NaN)
    end
    return optics
end

"""
    set_valid_subapertures!(sensor, valid_subapertures)

Install an explicit Shack–Hartmann valid-subaperture mask. A prepared optical
owner that binds the preceding layout revision must be prepared again.
"""
function set_valid_subapertures!(
    wfs::ShackHartmannWFS,
    valid_subapertures::AbstractMatrix{Bool},
)
    set_valid_subapertures!(wfs.front_end.layout, valid_subapertures)
    return wfs
end

function ensure_sh_asterism_buffers!(
    optics::ShackHartmannOptics, n_sources::Int)
    propagation = microlens_propagation_workspace(optics.propagation)
    n_sources > 0 || throw(InvalidConfiguration("asterism must contain at least one source"))
    if n_sources > propagation.asterism_capacity
        pad = size(propagation.fft_stack, 1)
        n_spots = n_lenslets(optics)^2
        total = n_spots * n_sources
        propagation.fft_asterism_stack = similar(
            propagation.fft_asterism_stack, pad, pad, total)
        propagation.intensity_tmp_stack = similar(
            propagation.intensity_tmp_stack,
            eltype(propagation.intensity_tmp_stack), pad, pad, total)
        propagation.fft_asterism_plan = plan_repeated_fft_backend!(
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
    return optics
end

function build_sh_phasor!(optics::ShackHartmannOptics,
    ratio::T) where {T<:AbstractFloat}
    propagation = microlens_propagation_workspace(optics.propagation)
    if size(propagation.phasor, 1) == 0
        return optics
    end
    if isequal(propagation.phasor_ratio, ratio)
        return optics
    end
    n = size(propagation.phasor, 1)
    scale = -T(π) * (T(n) + one(T) + ratio) / T(n)
    host = Matrix{Complex{T}}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        host[i, j] = cis(scale * (i + j - 2))
    end
    copyto!(propagation.phasor, host)
    propagation.phasor_ratio = ratio
    return optics
end

@inline sh_pixel_scale_init(d_subap::Real, padding::Int,
    wavelength_m::Real) = lgs_pixel_scale(d_subap, padding, wavelength_m)

@inline sh_pixel_scale_init(d_subap::Real, padding::Int,
    src::AbstractSource) = sh_pixel_scale_init(d_subap, padding,
    wavelength(src))

function _prepare_microlens_sampling_wavelength!(
    optics::ShackHartmannOptics,
    pupil_resolution::Int, pupil_diameter_m::Real, wavelength_m::Real)
    sampling = _sh_microlens_sampling_configuration(optics,
        pupil_resolution, pupil_diameter_m, wavelength_m)
    propagation = microlens_propagation_workspace(optics.propagation)
    padding = sampling.padding
    pad = sampling.padded_subaperture_samples
    binning_pixel_scale = sampling.binning_pixel_scale
    n_pix_subap = sampling.spot_samples_per_axis

    if padding != propagation.effective_padding ||
            pad != size(propagation.field, 1)
        ensure_sh_buffers!(optics, pad)
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
            n_lenslets(optics)^2, n_pix_subap, n_pix_subap)
        propagation.spot_cube_accum = similar(propagation.spot_cube_accum,
            eltype(propagation.spot_cube_accum),
            n_lenslets(optics)^2, n_pix_subap, n_pix_subap)
        propagation.sampled_n_pix_subap = n_pix_subap
    end

    n_binned = div(pad, binning_pixel_scale)
    if size(propagation.bin_buffer) != (n_binned, n_binned)
        propagation.bin_buffer = similar(
            propagation.bin_buffer,
            n_binned,
            n_binned,
        )
    end

    propagation.binning_pixel_scale = binning_pixel_scale
    T = eltype(propagation.intensity)
    half_shift_ratio = microlens_array(optics).params.half_pixel_shift ?
        T(binning_pixel_scale) : zero(T)
    build_sh_phasor!(optics, half_shift_ratio)
    return optics
end

function _sh_microlens_sampling_configuration(
    optics::ShackHartmannOptics,
    pupil_resolution::Int, pupil_diameter_m::Real, wavelength_m::Real)
    pupil_resolution % n_lenslets(optics) == 0 ||
        throw(InvalidConfiguration(
            "pupil resolution must be divisible by n_lenslets"))
    sub = div(pupil_resolution, n_lenslets(optics))
    microlens = microlens_array(optics)
    padding = microlens.params.diffraction_padding
    pixel_scale_req =
        microlens.params.pixel_scale_arcsec
    d_subap = pupil_diameter_m / n_lenslets(optics)
    pixel_scale_init = sh_pixel_scale_init(d_subap, padding, wavelength_m)

    if pixel_scale_req !== nothing
        while pixel_scale_req / pixel_scale_init < 0.95
            padding += 1
            pixel_scale_init = sh_pixel_scale_init(d_subap, padding,
                wavelength_m)
        end
    end

    binning_pixel_scale = if pixel_scale_req === nothing
        microlens.params.shannon_sampling ? 1 : 2
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

    n_pix_subap = microlens.params.n_pix_subap === nothing ?
        sub : microlens.params.n_pix_subap
    if isodd(n_pix_subap)
        throw(InvalidConfiguration("n_pix_subap must be even"))
    end

    return (
        padding=padding,
        padded_subaperture_samples=pad,
        binning_pixel_scale=binning_pixel_scale,
        spot_samples_per_axis=n_pix_subap,
        pixel_scale_arcsec=pixel_scale_init * binning_pixel_scale,
    )
end
