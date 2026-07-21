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

"""Immutable physical definition of a pyramid focal-plane phase mask."""
struct PyramidPhaseMask{T<:AbstractFloat}
    old_mask::Bool
    rooftop::T
    theta_rotation::T
    mask_scale::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
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

"""Single-writer FFT and scratch state for pyramid-mask propagation."""
mutable struct PreparedPyramidPropagation{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}}}
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
end

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

"""Mutable detector-plane sampling storage owned by acquisition."""
mutable struct PyramidAcquisitionState{T<:AbstractFloat,R<:AbstractMatrix{T}}
    binned_intensity::R
    camera_frame::R
    nominal_detector_resolution::Int
end

struct PyramidDetectorAcquisition{S}
    binning::Int
    state::S
end

"""Mutable support, calibration, and output storage for slope estimation."""
mutable struct PyramidEstimatorState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},V<:AbstractVector{T},I<:AbstractVector{Int},
    R<:AbstractMatrix{T}}
    valid_mask::A
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
    shift_x::NTuple{4,Int}
    shift_y::NTuple{4,Int}
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
    calibration_revision::UInt
end

struct PyramidDifferentialEstimator{P<:PyramidEstimatorParams,S}
    params::P
    state::S
end

struct PyramidWFS{M<:SensingMode,F,A,E,B<:AbstractArrayBackend} <: AbstractWFS
    front_end::F
    acquisition::A
    estimator::E
end

@inline backend(::PyramidWFS{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

@inline pyramid_estimator_params(wfs::PyramidWFS) = wfs.estimator.params
@inline pyramid_estimator_state(wfs::PyramidWFS) = wfs.estimator.state
@inline pyramid_acquisition_state(wfs::PyramidWFS{<:Diffractive}) =
    wfs.acquisition.state
@inline pyramid_propagation(wfs::PyramidWFS{<:Diffractive}) =
    wfs.front_end.propagation
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
"""
function PyramidWFS(tel::Telescope; pupil_samples::Int, threshold::Real=0.1, modulation::Real=2.0,
    light_ratio::Real=0.0,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    modulation_points::Union{Int,Nothing}=nothing, extra_modulation_factor::Int=0,
    old_mask::Bool=false, rooftop::Real=0.0, theta_rotation::Real=0.0, delta_theta::Real=0.0,
    user_modulation_path=nothing, mask_scale::Real=1.0, diffraction_padding::Int=2,
    psf_centering::Bool=true, n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
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
        T(theta_rotation),
        typed_mask_scale,
        diffraction_padding, psf_centering, n_pix_separation, n_pix_edge)
    operating_policy = legacy_modulation_policy(T(modulation),
        modulation_points, extra_modulation_factor, T(delta_theta),
        user_modulation_path)
    calibration_policy = calibration_modulation_policy(operating_policy,
        T(calib_modulation), T(delta_theta))
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
    estimator_state = PyramidEstimatorState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(valid_signal_indices),
        typeof(signal_2d),
    }(
        valid_mask,
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
        (0, 0, 0, 0),
        (0, 0, 0, 0),
        false,
        zero(T),
        UInt(0),
        UInt(0),
    )
    estimator = PyramidDifferentialEstimator(estimator_params,
        estimator_state)
    front_end, acquisition = prepare_pyramid_mode(mode, backend, T, tel,
        phase_mask, operating_policy, calibration_policy, pupil_samples,
        binning)
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
    binning) where {T} =
    (nothing, nothing)

function prepare_pyramid_mode(::Diffractive, backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning) where {
    T<:AbstractFloat,
}
    return _prepare_pyramid_diffractive_storage(backend, T, tel, phase_mask,
        operating_policy, calibration_policy, pupil_samples, binning)
end

function _prepare_pyramid_diffractive_storage(backend, ::Type{T}, tel,
    phase_mask, operating_policy, calibration_policy, pupil_samples,
    binning) where {
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
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    propagation = PreparedPyramidPropagation(
        field, focal_field, pupil_field, mask, phasor, intensity, temp,
        scratch, asterism_stack, fft_plan, ifft_plan, elongation_kernel,
        lgs_kernel_fft, UInt(0), pad, 1, UInt(0))
    modulation = prepare_focal_plane_modulation(operating_policy,
        tel.params.resolution, field, T)
    calibration_modulation = prepare_focal_plane_modulation(
        calibration_policy, tel.params.resolution, field, T)
    front_end = PyramidOpticalFrontEnd(phase_mask, modulation,
        calibration_modulation, propagation, pupil_samples, binning, nothing)
    binned_intensity = similar(intensity)
    subaperture_pixels = div(tel.params.resolution, pupil_samples)
    nominal = div(pad, subaperture_pixels)
    camera_frame = backend{T}(undef, nominal, nominal)
    acquisition = PyramidDetectorAcquisition(binning,
        PyramidAcquisitionState(binned_intensity, camera_frame, nominal))
    return front_end, acquisition
end

@inline prepare_pyramid_front_end!(::Geometric, ::PyramidWFS,
    ::Telescope) = nothing

function prepare_pyramid_front_end!(::Diffractive, wfs::PyramidWFS,
    tel::Telescope)
    build_pyramid_phasor!(pyramid_propagation(wfs).phasor)
    build_pyramid_mask!(wfs, PupilFunction(tel))
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
    throw(WFSPreparationError(:optical_formation, :unsupported,
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
    if size(wfs.front_end.propagation.field) != (pad, pad)
        wfs.front_end.propagation.revision += UInt(1)
        wfs.front_end.propagation.field = similar(wfs.front_end.propagation.field, pad, pad)
        wfs.front_end.propagation.focal_field = similar(wfs.front_end.propagation.focal_field, pad, pad)
        wfs.front_end.propagation.pupil_field = similar(wfs.front_end.propagation.pupil_field, pad, pad)
        wfs.front_end.propagation.pyramid_mask = similar(wfs.front_end.propagation.pyramid_mask, pad, pad)
        wfs.front_end.propagation.phasor = similar(wfs.front_end.propagation.phasor, pad, pad)
        wfs.front_end.propagation.intensity = similar(wfs.front_end.propagation.intensity, pad, pad)
        wfs.front_end.propagation.temp = similar(wfs.front_end.propagation.temp, pad, pad)
        wfs.front_end.propagation.scratch = similar(wfs.front_end.propagation.scratch, pad, pad)
        wfs.acquisition.state.binned_intensity = similar(wfs.acquisition.state.binned_intensity, pad, pad)
        wfs.front_end.propagation.asterism_stack = similar(wfs.front_end.propagation.asterism_stack, pad, pad, wfs.front_end.propagation.asterism_capacity)
        wfs.front_end.propagation.fft_plan = plan_fft_backend!(wfs.front_end.propagation.focal_field)
        wfs.front_end.propagation.ifft_plan = plan_ifft_backend!(wfs.front_end.propagation.pupil_field)
        wfs.front_end.propagation.lgs_kernel_fft = similar(
            wfs.front_end.propagation.focal_field,
            eltype(wfs.front_end.propagation.focal_field), 0, 0)
        wfs.front_end.propagation.lgs_kernel_tag = UInt(0)
        wfs.front_end.propagation.effective_resolution = pad
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
        build_pyramid_phasor!(wfs.front_end.propagation.phasor)
        build_pyramid_mask!(wfs, pupil)
    end
    return wfs
end

function ensure_pyramid_asterism_stack!(wfs::PyramidWFS, n_src::Int)
    n_src >= 1 || throw(InvalidConfiguration("asterism source count must be >= 1"))
    pad = size(wfs.front_end.propagation.intensity, 1)
    if size(wfs.front_end.propagation.asterism_stack, 1) != pad || size(wfs.front_end.propagation.asterism_stack, 2) != pad ||
            size(wfs.front_end.propagation.asterism_stack, 3) < n_src
        capacity = max(n_src, wfs.front_end.propagation.asterism_capacity)
        wfs.front_end.propagation.asterism_stack = similar(wfs.front_end.propagation.asterism_stack, pad, pad, capacity)
        wfs.front_end.propagation.asterism_capacity = capacity
    end
    return wfs.front_end.propagation.asterism_stack
end

@inline grouped_staging_buffer(wfs::PyramidWFS, out::AbstractMatrix) = wfs.front_end.propagation.intensity

function accumulate_pyramid_asterism_intensity!(::ScalarCPUStyle, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, wfs.front_end.propagation.intensity, stack, ast.sources, pyramid_intensity!, wfs, pupil)
end

function accumulate_pyramid_asterism_intensity!(style::AcceleratorStyle, wfs::PyramidWFS, pupil::PupilFunction, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(style, wfs, wfs.front_end.propagation.intensity, stack, ast.sources, pyramid_intensity!, wfs, pupil)
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
    T = eltype(wfs.front_end.propagation.intensity)
    @inbounds for (sample_idx, sample) in pairs(src.bundle.samples)
        channel_qe = pyramid_spectral_component_qe(qe_model, sample, T)
        variant = source_with_wavelength_and_radiometric_value(src, sample.wavelength,
            T(total_irradiance * sample.weight * channel_qe))
        pyramid_intensity_core!(@view(stack[:, :, sample_idx]), wfs, pupil,
            variant, modulation)
    end
    return reduce_grouped_stack!(style, wfs.front_end.propagation.intensity, stack, count)
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
    out === wfs.front_end.propagation.intensity || copyto!(out, wfs.front_end.propagation.intensity)
    return out
end

@inline function pyramid_support_selection_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil, src,
        qe_model, pyramid_calibration_modulation(wfs))
    out === wfs.front_end.propagation.intensity || copyto!(out, wfs.front_end.propagation.intensity)
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
    out === wfs.front_end.propagation.intensity ||
        copyto!(out, wfs.front_end.propagation.intensity)
    return out
end

@inline function pyramid_calibration_intensity!(out::AbstractMatrix,
    wfs::PyramidWFS, pupil::PupilFunction, src::SpectralSource,
    qe_model::AbstractQuantumEfficiencyModel)
    accumulate_pyramid_spectral_intensity!(execution_style(out), wfs, pupil,
        src, qe_model, pyramid_operating_modulation(wfs))
    out === wfs.front_end.propagation.intensity ||
        copyto!(out, wfs.front_end.propagation.intensity)
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
    if pad % wfs.acquisition.binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide padded resolution"))
    end
    if _pupil_resolution(pupil) % wfs.acquisition.binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide telescope resolution"))
    end
    ensure_pyramid_buffers!(wfs, pad, pupil)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration("sample_pyramid_intensity! requires telescope context"))
end

function resize_pyramid_signal_buffers!(wfs::PyramidWFS, frame_size::Int)
    nominal = wfs.acquisition.state.nominal_detector_resolution
    frame_size >= 1 || throw(InvalidConfiguration(
        "pyramid camera frame size must be >= 1"))
    nominal > 0 || throw(InvalidConfiguration(
        "pyramid nominal detector resolution must be prepared before signal resizing"))
    nominal % wfs.acquisition.binning == 0 || throw(InvalidConfiguration(
        "pyramid binning must evenly divide the nominal detector resolution"))
    sampled_size = div(nominal, wfs.acquisition.binning)
    sampled_size % frame_size == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must evenly divide the pyramid camera frame"))
    reduction = div(sampled_size, frame_size)
    total_sampling = wfs.acquisition.binning * reduction
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
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        calibration_storage_changed = true
    elseif size(wfs.estimator.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.estimator.state.reference_signal_2d = similar(wfs.estimator.state.reference_signal_2d, 2 * n_pixels, n_pixels)
        fill!(wfs.estimator.state.reference_signal_2d,
            zero(eltype(wfs.estimator.state.reference_signal_2d)))
        calibration_storage_changed = true
    end
    if size(wfs.acquisition.state.camera_frame) != (frame_size, frame_size)
        wfs.acquisition.state.camera_frame = similar(wfs.acquisition.state.camera_frame, frame_size, frame_size)
    end
    update_pyramid_valid_signal!(wfs)
    if calibration_storage_changed
        wfs.estimator.state.calibrated = false
        wfs.estimator.state.calibration_revision += UInt(1)
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
    n_pixels = size(wfs.estimator.state.signal_2d, 2)
    nominal = wfs.acquisition.state.nominal_detector_resolution
    nominal > 0 || throw(InvalidConfiguration(
        "pyramid nominal detector resolution must be prepared before signal extraction"))
    nominal % wfs.acquisition.binning == 0 || throw(InvalidConfiguration(
        "pyramid binning must evenly divide the nominal detector resolution"))
    sampled_size = div(nominal, wfs.acquisition.binning)
    sampled_size % n_rows == 0 || throw(InvalidConfiguration(
        "detector sampling and binning must evenly divide the pyramid camera frame"))
    reduction = div(sampled_size, n_rows)
    total_sampling = wfs.acquisition.binning * reduction
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
    n_pixels = size(wfs.estimator.state.valid_i4q, 1)
    fill!(wfs.estimator.state.valid_signal, false)
    @views begin
        wfs.estimator.state.valid_signal[1:n_pixels, :] .= wfs.estimator.state.valid_i4q
        wfs.estimator.state.valid_signal[n_pixels+1:end, :] .= wfs.estimator.state.valid_i4q
    end
    return wfs
end

function update_pyramid_valid_signal_indices!(wfs::PyramidWFS)
    valid_host = wfs.estimator.state.valid_i4q_host
    if size(valid_host) != size(wfs.estimator.state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(wfs.estimator.state.valid_i4q)...)
        wfs.estimator.state.valid_i4q_host = valid_host
    end
    copyto!(valid_host, wfs.estimator.state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(wfs.estimator.state.valid_signal_indices) < n_valid
        wfs.estimator.state.valid_signal_indices = similar(wfs.estimator.state.valid_signal_indices, n_valid)
    end
    if length(wfs.estimator.state.valid_signal_indices_host) < n_valid
        wfs.estimator.state.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = wfs.estimator.state.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(wfs.estimator.state.valid_signal_indices, 1, host_indices, 1, n_valid)
    wfs.estimator.state.valid_signal_count = n_valid
    return n_valid
end

function resize_pyramid_slope_buffers!(wfs::PyramidWFS)
    n_valid = wfs.estimator.state.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("pyramid valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(wfs.estimator.state.slopes) != n_slopes
        wfs.estimator.state.slopes = similar(wfs.estimator.state.slopes, n_slopes)
    end
    if length(wfs.estimator.state.optical_gain) != n_slopes
        wfs.estimator.state.optical_gain = similar(wfs.estimator.state.optical_gain, n_slopes)
        fill!(wfs.estimator.state.optical_gain, one(eltype(wfs.estimator.state.optical_gain)))
    end
    return wfs
end

function pyramid_valid_flux_sum!(::ScalarCPUStyle, wfs::PyramidWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    return masked_sum2d(ScalarCPUStyle(), i4q, wfs.estimator.state.valid_i4q_host)
end

function pyramid_valid_flux_sum!(style::AcceleratorStyle, wfs::PyramidWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.estimator.state.valid_i4q,
        wfs.estimator.state.valid_i4q_host,
        wfs.estimator.state.valid_flux_sum_buffer,
        wfs.estimator.state.valid_flux_sum_host,
        wfs.estimator.state.valid_flux_i4q_host,
    )
    wfs.estimator.state.valid_flux_i4q_host = host_parent
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
    i4q = @view wfs.estimator.state.signal_2d[1:n_pixels, :]
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
    pyramid_support_selection_intensity!(wfs.front_end.propagation.temp,
        wfs, pupil, src, qe_model)
    sampled = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.temp)
    frame = detector_calibration_frame!(det, sampled,
        pyramid_detector_calibration_qe(src, det, eltype(det.state.frame)))
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
    n_pixels = div(wfs.estimator.params.pupil_samples, wfs.acquisition.binning)
    if size(wfs.estimator.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.valid_i4q = similar(wfs.estimator.state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(wfs.estimator.state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    pyramid_support_selection_intensity!(wfs.front_end.propagation.temp,
        wfs, pupil, src, qe_model)
    frame = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.temp)

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
        wfs.estimator.state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
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
    n_pixels = div(wfs.estimator.params.pupil_samples, wfs.acquisition.binning)
    if size(wfs.estimator.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.estimator.state.valid_i4q = similar(wfs.estimator.state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.estimator.params.light_ratio)
        fill!(wfs.estimator.state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    pyramid_support_selection_intensity!(wfs.front_end.propagation.temp,
        wfs, pupil, src, qe_model)
    frame = sample_pyramid_intensity!(wfs, pupil, wfs.front_end.propagation.temp)

    center, n_extra = require_pyramid_frame_geometry(wfs, frame)
    rows_lo = center - n_extra - n_pixels + 1:center - n_extra
    rows_hi = center + n_extra + 1:center + n_extra + n_pixels
    cols_lo = center - n_extra - n_pixels + 1:center - n_extra
    cols_hi = center + n_extra + 1:center + n_extra + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view wfs.estimator.state.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.estimator.params.light_ratio * maximum(i4q)
    @. wfs.estimator.state.valid_i4q = i4q >= cutoff
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, pupil::PupilFunction, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.acquisition.binning
    sub = div(_pupil_resolution(pupil), wfs.estimator.params.pupil_samples)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("pyramid intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    wfs.acquisition.state.nominal_detector_resolution = n_camera
    if size(wfs.acquisition.state.camera_frame) != (n_camera, n_camera)
        wfs.acquisition.state.camera_frame = similar(wfs.acquisition.state.camera_frame, n_camera, n_camera)
    end
    frame = wfs.acquisition.state.camera_frame
    if binning != 1
        if n_camera % binning != 0
            throw(InvalidConfiguration("pyramid binning must evenly divide detector resolution"))
        end
        n_binned = div(n_camera, binning)
        if size(wfs.acquisition.state.binned_intensity) != (n_binned, n_binned)
            wfs.acquisition.state.binned_intensity = similar(wfs.acquisition.state.binned_intensity, n_binned, n_binned)
        end
        bin2d!(wfs.acquisition.state.binned_intensity, intensity, sub * binning)
        frame = wfs.acquisition.state.binned_intensity
    else
        bin2d!(wfs.acquisition.state.camera_frame, intensity, sub)
    end
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    return frame
end
