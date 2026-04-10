using Statistics

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

@kernel function modulation_phases_kernel!(phases, mod_amp, center, n_pts::Int, n::Int)
    i, j, p = @index(Global, NTuple)
    if i <= n && j <= n && p <= n_pts
        x = (i - center) / n
        y = (j - center) / n
        angle = 2 * π * (p - 1) / n_pts
        s, c = sincos(angle)
        phase = 2 * π * mod_amp * (x * c + y * s)
        @inbounds phases[i, j, p] = cis(phase)
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

struct PyramidParams{T<:AbstractFloat,M,N<:WFSNormalization}
    n_subap::Int
    threshold::T
    light_ratio::T
    normalization::N
    calib_modulation::T
    modulation::T
    modulation_points::Int
    extra_modulation_factor::Int
    old_mask::Bool
    rooftop::T
    theta_rotation::T
    delta_theta::T
    user_modulation_path::M
    mask_scale::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
    binning::Int
end

mutable struct PyramidState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    I<:AbstractVector{Int},
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    RB<:AbstractMatrix{T},
    RS<:AbstractArray{T,3},
    C3<:AbstractArray{Complex{T},3},
    Pf,
    Pi,
    K<:AbstractVector{T},
    Kf<:AbstractMatrix{Complex{T}},
    Vg<:AbstractVector{T}}
    valid_mask::A
    slopes::V
    field::C
    focal_field::C
    pupil_field::C
    pyramid_mask::C
    phasor::C
    modulation_phases::C3
    intensity::R
    temp::R
    scratch::R
    binned_intensity::RB
    asterism_stack::RS
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::Kf
    lgs_kernel_tag::UInt
    optical_gain::Vg
    valid_i4q::A
    valid_i4q_host::Matrix{Bool}
    valid_signal::A
    valid_signal_indices::I
    valid_signal_indices_host::Vector{Int}
    valid_signal_count::Int
    valid_flux_sum_buffer::Vg
    valid_flux_sum_host::Vector{T}
    valid_flux_i4q_host::Matrix{T}
    signal_2d::R
    reference_signal_2d::R
    camera_frame::RB
    shift_x::NTuple{4,Int}
    shift_y::NTuple{4,Int}
    effective_resolution::Int
    nominal_detector_resolution::Int
    asterism_capacity::Int
    calibrated::Bool
    calibration_wavelength::T
    calibration_signature::UInt
end

struct PyramidWFS{M<:SensingMode,P<:PyramidParams,S<:PyramidState} <: AbstractWFS
    params::P
    state::S
end

"""
    PyramidWFS(tel; ...)

Construct a pyramid wavefront sensor.

The diffractive model forms four re-imaged pupil intensities through a
focal-plane pyramid mask. Slopes are obtained from left/right and top/bottom
intensity differences after optional modulation averaging and binning.
"""
function PyramidWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1, modulation::Real=2.0,
    light_ratio::Real=0.0,
    normalization::WFSNormalization=MeanValidFluxNormalization(),
    calib_modulation::Real=min(50.0, tel.params.resolution / 2 - 1),
    modulation_points::Union{Int,Nothing}=nothing, extra_modulation_factor::Int=0,
    old_mask::Bool=false, rooftop::Real=0.0, theta_rotation::Real=0.0, delta_theta::Real=0.0,
    user_modulation_path=nothing, mask_scale::Real=1.0, diffraction_padding::Int=2,
    psf_centering::Bool=true, n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())

    backend = _resolve_array_backend(backend)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_mod = resolve_modulation_points(T(modulation), modulation_points, extra_modulation_factor, user_modulation_path)
    params = PyramidParams{T,typeof(user_modulation_path),typeof(normalization)}(
        n_subap,
        T(threshold),
        T(light_ratio),
        normalization,
        T(calib_modulation),
        T(modulation),
        n_mod,
        extra_modulation_factor,
        old_mask,
        T(rooftop),
        T(theta_rotation),
        T(delta_theta),
        user_modulation_path,
        T(mask_scale),
        diffraction_padding, psf_centering, n_pix_separation, n_pix_edge, binning)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    pad = tel.params.resolution * diffraction_padding
    if n_pix_separation !== nothing
        edge = n_pix_edge === nothing ? div(n_pix_separation, 2) : n_pix_edge
        pad = Int(round((n_subap * 2 + n_pix_separation + 2 * edge) * tel.params.resolution / n_subap))
    end
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    pyramid_mask = similar(field)
    phasor = similar(field)
    modulation_phases = backend{Complex{T}}(undef, tel.params.resolution, tel.params.resolution, n_mod)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    scratch = similar(intensity)
    binned_intensity = similar(intensity)
    asterism_stack = backend{T}(undef, pad, pad, 1)
    n_pix_signal = cld(n_subap, binning)
    valid_i4q = backend{Bool}(undef, n_pix_signal, n_pix_signal)
    valid_i4q_host = Matrix{Bool}(undef, n_pix_signal, n_pix_signal)
    valid_signal = backend{Bool}(undef, 2 * n_pix_signal, n_pix_signal)
    valid_signal_indices = backend{Int}(undef, n_pix_signal * n_pix_signal)
    valid_signal_indices_host = Vector{Int}(undef, n_pix_signal * n_pix_signal)
    valid_flux_sum_buffer = backend{T}(undef, 1)
    valid_flux_sum_host = Vector{T}(undef, 1)
    valid_flux_i4q_host = Matrix{T}(undef, size(temp)...)
    signal_2d = backend{T}(undef, 2 * n_pix_signal, n_pix_signal)
    reference_signal_2d = similar(signal_2d)
    camera_frame = backend{T}(undef, max(1, n_subap), max(1, n_subap))
    fft_plan = plan_fft_backend!(focal_field)
    ifft_plan = plan_ifft_backend!(pupil_field)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0)
    optical_gain = similar(slopes)
    fill!(optical_gain, one(T))
    state = PyramidState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(valid_signal_indices),
        typeof(field),
        typeof(intensity),
        typeof(binned_intensity),
        typeof(asterism_stack),
        typeof(modulation_phases),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
        typeof(optical_gain),
    }(
        valid_mask,
        slopes,
        field,
        focal_field,
        pupil_field,
        pyramid_mask,
        phasor,
        modulation_phases,
        intensity,
        temp,
        scratch,
        binned_intensity,
        asterism_stack,
        fft_plan,
        ifft_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
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
        signal_2d,
        reference_signal_2d,
        camera_frame,
        (0, 0, 0, 0),
        (0, 0, 0, 0),
        pad,
        0,
        1,
        false,
        zero(T),
        UInt(0),
    )
    wfs = PyramidWFS{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    build_pyramid_phasor!(wfs.state.phasor)
    build_pyramid_mask!(wfs, tel)
    build_modulation_phases!(wfs, tel)
    return wfs
end

sensing_mode(::PyramidWFS{M}) where {M} = M()

function resolve_modulation_points(modulation::T, modulation_points::Union{Int,Nothing},
    extra_modulation_factor::Int, user_modulation_path) where {T<:AbstractFloat}
    if user_modulation_path !== nothing
        n = length(user_modulation_path)
        n >= 1 || throw(InvalidConfiguration("user_modulation_path must contain at least one point"))
        return n
    end
    if modulation == 0
        return 1
    end
    if modulation_points === nothing
        perimeter = T(2 * π) * modulation
        return max(1, 4 * Int(extra_modulation_factor + ceil(perimeter / 4)))
    end
    modulation_points >= 1 || throw(InvalidConfiguration("modulation_points must be >= 1"))
    return modulation_points
end

function update_valid_mask!(wfs::PyramidWFS, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function ensure_pyramid_buffers!(wfs::PyramidWFS, pad::Int, tel::Telescope)
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.focal_field = similar(wfs.state.focal_field, pad, pad)
        wfs.state.pupil_field = similar(wfs.state.pupil_field, pad, pad)
        wfs.state.pyramid_mask = similar(wfs.state.pyramid_mask, pad, pad)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.intensity = similar(wfs.state.intensity, pad, pad)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.scratch = similar(wfs.state.scratch, pad, pad)
        wfs.state.binned_intensity = similar(wfs.state.binned_intensity, pad, pad)
        wfs.state.asterism_stack = similar(wfs.state.asterism_stack, pad, pad, wfs.state.asterism_capacity)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.focal_field)
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.pupil_field)
        wfs.state.lgs_kernel_fft = similar(wfs.state.focal_field, Complex{eltype(wfs.state.focal_field)}, 0, 0)
        wfs.state.lgs_kernel_tag = UInt(0)
        wfs.state.effective_resolution = pad
        wfs.state.calibrated = false
        build_pyramid_phasor!(wfs.state.phasor)
        build_pyramid_mask!(wfs, tel)
    end
    return wfs
end

function ensure_pyramid_asterism_stack!(wfs::PyramidWFS, n_src::Int)
    n_src >= 1 || throw(InvalidConfiguration("asterism source count must be >= 1"))
    pad = size(wfs.state.intensity, 1)
    if size(wfs.state.asterism_stack, 1) != pad || size(wfs.state.asterism_stack, 2) != pad ||
            size(wfs.state.asterism_stack, 3) < n_src
        capacity = max(n_src, wfs.state.asterism_capacity)
        wfs.state.asterism_stack = similar(wfs.state.asterism_stack, pad, pad, capacity)
        wfs.state.asterism_capacity = capacity
    end
    return wfs.state.asterism_stack
end

@inline grouped_staging_buffer(wfs::PyramidWFS, out::AbstractMatrix) = wfs.state.intensity

function accumulate_pyramid_asterism_intensity!(::ScalarCPUStyle, wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, wfs.state.intensity, stack, ast.sources, pyramid_intensity!, wfs, tel)
end

function accumulate_pyramid_asterism_intensity!(style::AcceleratorStyle, wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(style, wfs, wfs.state.intensity, stack, ast.sources, pyramid_intensity!, wfs, tel)
end

function accumulate_pyramid_spectral_intensity!(::ScalarCPUStyle, wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    count = length(src.bundle.samples)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    total_flux = photon_flux(src)
    @inbounds for (sample_idx, sample) in pairs(src.bundle.samples)
        variant = source_with_wavelength_and_flux(src, sample.wavelength,
            eltype(wfs.state.intensity)(total_flux * sample.weight))
        pyramid_intensity!(@view(stack[:, :, sample_idx]), wfs, tel, variant)
    end
    return reduce_grouped_stack!(ScalarCPUStyle(), wfs.state.intensity, stack, count)
end

function accumulate_pyramid_spectral_intensity!(style::AcceleratorStyle, wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    count = length(src.bundle.samples)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    total_flux = photon_flux(src)
    @inbounds for (sample_idx, sample) in pairs(src.bundle.samples)
        variant = source_with_wavelength_and_flux(src, sample.wavelength,
            eltype(wfs.state.intensity)(total_flux * sample.weight))
        pyramid_intensity!(@view(stack[:, :, sample_idx]), wfs, tel, variant)
    end
    return reduce_grouped_stack!(style, wfs.state.intensity, stack, count)
end

function accumulate_pyramid_extended_intensity!(::ScalarCPUStyle, out::AbstractMatrix, wfs::PyramidWFS,
    tel::Telescope, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return pyramid_intensity!(out, wfs, tel, ast.sources[1])
    end
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(ScalarCPUStyle(), wfs, out, stack, ast.sources, pyramid_intensity!, wfs, tel)
end

function accumulate_pyramid_extended_intensity!(style::AcceleratorStyle, out::AbstractMatrix, wfs::PyramidWFS,
    tel::Telescope, src::ExtendedSource)
    ast = extended_source_asterism(src)
    if length(ast.sources) == 1
        return pyramid_intensity!(out, wfs, tel, ast.sources[1])
    end
    count = length(ast.sources)
    stack = grouped_stack_view(ensure_pyramid_asterism_stack!(wfs, count), count)
    return accumulate_grouped_sources!(style, wfs, out, stack, ast.sources, pyramid_intensity!, wfs, tel)
end

function prepare_pyramid_sampling!(wfs::PyramidWFS, tel::Telescope)
    n_sub = wfs.params.n_subap
    pad = tel.params.resolution * wfs.params.diffraction_padding
    if wfs.params.n_pix_separation !== nothing
        edge = wfs.params.n_pix_edge === nothing ? div(wfs.params.n_pix_separation, 2) : wfs.params.n_pix_edge
        pad = Int(round((n_sub * 2 + wfs.params.n_pix_separation + 2 * edge) * tel.params.resolution / n_sub))
    end
    if pad < tel.params.resolution
        throw(InvalidConfiguration("pyramid padding must be >= telescope resolution"))
    end
    if pad % wfs.params.binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide padded resolution"))
    end
    if tel.params.resolution % wfs.params.binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide telescope resolution"))
    end
    ensure_pyramid_buffers!(wfs, pad, tel)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration("sample_pyramid_intensity! requires telescope context"))
end

function resize_pyramid_signal_buffers!(wfs::PyramidWFS, frame_size::Int)
    nominal = wfs.state.nominal_detector_resolution
    n_pixels = cld(wfs.params.n_subap, wfs.params.binning)
    if size(wfs.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.state.valid_i4q = similar(wfs.state.valid_i4q, n_pixels, n_pixels)
        fill!(wfs.state.valid_i4q, false)
    end
    if size(wfs.state.valid_signal) != (2 * n_pixels, n_pixels)
        wfs.state.valid_signal = similar(wfs.state.valid_signal, 2 * n_pixels, n_pixels)
    end
    if size(wfs.state.signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.signal_2d = similar(wfs.state.signal_2d, 2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    elseif size(wfs.state.reference_signal_2d) != (2 * n_pixels, n_pixels)
        wfs.state.reference_signal_2d = similar(wfs.state.reference_signal_2d, 2 * n_pixels, n_pixels)
    end
    if size(wfs.state.camera_frame) != (frame_size, frame_size)
        wfs.state.camera_frame = similar(wfs.state.camera_frame, frame_size, frame_size)
    end
    update_pyramid_valid_signal!(wfs)
    return wfs
end

function update_pyramid_valid_signal!(wfs::PyramidWFS)
    n_pixels = size(wfs.state.valid_i4q, 1)
    fill!(wfs.state.valid_signal, false)
    @views begin
        wfs.state.valid_signal[1:n_pixels, :] .= wfs.state.valid_i4q
        wfs.state.valid_signal[n_pixels+1:end, :] .= wfs.state.valid_i4q
    end
    return wfs
end

function update_pyramid_valid_signal_indices!(wfs::PyramidWFS)
    valid_host = wfs.state.valid_i4q_host
    if size(valid_host) != size(wfs.state.valid_i4q)
        valid_host = Matrix{Bool}(undef, size(wfs.state.valid_i4q)...)
        wfs.state.valid_i4q_host = valid_host
    end
    copyto!(valid_host, wfs.state.valid_i4q)
    n_pixels = size(valid_host, 1)
    n_valid = count(valid_host)
    if length(wfs.state.valid_signal_indices) < n_valid
        wfs.state.valid_signal_indices = similar(wfs.state.valid_signal_indices, n_valid)
    end
    if length(wfs.state.valid_signal_indices_host) < n_valid
        wfs.state.valid_signal_indices_host = Vector{Int}(undef, n_valid)
    end
    host_indices = wfs.state.valid_signal_indices_host
    idx = 1
    @inbounds for i in 1:n_pixels, j in 1:n_pixels
        if valid_host[i, j]
            host_indices[idx] = i + (j - 1) * (2 * n_pixels)
            idx += 1
        end
    end
    copyto!(wfs.state.valid_signal_indices, 1, host_indices, 1, n_valid)
    wfs.state.valid_signal_count = n_valid
    return n_valid
end

function resize_pyramid_slope_buffers!(wfs::PyramidWFS)
    n_valid = wfs.state.valid_signal_count
    if n_valid == 0
        throw(InvalidConfiguration("pyramid valid pixel selection produced no valid signals"))
    end
    n_slopes = 2 * n_valid
    if length(wfs.state.slopes) != n_slopes
        wfs.state.slopes = similar(wfs.state.slopes, n_slopes)
    end
    if length(wfs.state.optical_gain) != n_slopes
        wfs.state.optical_gain = similar(wfs.state.optical_gain, n_slopes)
        fill!(wfs.state.optical_gain, one(eltype(wfs.state.optical_gain)))
    end
    return wfs
end

function pyramid_valid_flux_sum!(::ScalarCPUStyle, wfs::PyramidWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    return masked_sum2d(ScalarCPUStyle(), i4q, wfs.state.valid_i4q_host)
end

function pyramid_valid_flux_sum!(style::AcceleratorStyle, wfs::PyramidWFS, i4q::AbstractMatrix{T}) where {T<:AbstractFloat}
    summed, host_parent = masked_sum2d(
        style,
        i4q,
        wfs.state.valid_i4q,
        wfs.state.valid_i4q_host,
        wfs.state.valid_flux_sum_buffer,
        wfs.state.valid_flux_sum_host,
        wfs.state.valid_flux_i4q_host,
    )
    wfs.state.valid_flux_i4q_host = host_parent
    return summed
end

function select_pyramid_valid_i4q!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    return select_pyramid_valid_i4q!(execution_style(wfs.state.valid_i4q), wfs, tel, src)
end

function select_pyramid_valid_i4q!(::ScalarCPUStyle, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    n_pixels = cld(wfs.params.n_subap, wfs.params.binning)
    if size(wfs.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.state.valid_i4q = similar(wfs.state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.params.light_ratio)
        fill!(wfs.state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    copyto!(wfs.state.modulation_phases, host_modulation_phases(
        eltype(wfs.state.slopes),
        tel,
        wfs.params.calib_modulation,
        wfs.params.modulation_points,
        wfs.params.delta_theta,
        wfs.params.user_modulation_path,
    ))
    pyramid_intensity!(wfs.state.temp, wfs, tel, src)
    frame = sample_pyramid_intensity!(wfs, tel, wfs.state.temp)
    build_modulation_phases!(wfs, tel)

    n_extra = round(Int, (something(wfs.params.n_pix_separation, 0) / 2) / wfs.params.binning)
    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
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
    cutoff = wfs.params.light_ratio * max_i4q
    @inbounds for j in 1:n_pixels, i in 1:n_pixels
        q1 = frame[center - n_extra - n_pixels + i, center - n_extra - n_pixels + j]
        q2 = frame[center - n_extra - n_pixels + i, center + n_extra + j]
        q3 = frame[center + n_extra + i, center + n_extra + j]
        q4 = frame[center + n_extra + i, center - n_extra - n_pixels + j]
        wfs.state.valid_i4q[i, j] = (q1 + q2 + q3 + q4) >= cutoff
    end
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function select_pyramid_valid_i4q!(::AcceleratorStyle, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    n_pixels = cld(wfs.params.n_subap, wfs.params.binning)
    if size(wfs.state.valid_i4q) != (n_pixels, n_pixels)
        wfs.state.valid_i4q = similar(wfs.state.valid_i4q, n_pixels, n_pixels)
    end
    if iszero(wfs.params.light_ratio)
        fill!(wfs.state.valid_i4q, true)
        update_pyramid_valid_signal!(wfs)
        update_pyramid_valid_signal_indices!(wfs)
        resize_pyramid_slope_buffers!(wfs)
        return wfs
    end

    copyto!(wfs.state.modulation_phases, host_modulation_phases(
        eltype(wfs.state.slopes),
        tel,
        wfs.params.calib_modulation,
        wfs.params.modulation_points,
        wfs.params.delta_theta,
        wfs.params.user_modulation_path,
    ))
    pyramid_intensity!(wfs.state.temp, wfs, tel, src)
    frame = sample_pyramid_intensity!(wfs, tel, wfs.state.temp)
    build_modulation_phases!(wfs, tel)

    n_extra = round(Int, (something(wfs.params.n_pix_separation, 0) / 2) / wfs.params.binning)
    center = round(Int, wfs.state.nominal_detector_resolution / wfs.params.binning / 2)
    rows_lo = center - n_extra - n_pixels + 1:center - n_extra
    rows_hi = center + n_extra + 1:center + n_extra + n_pixels
    cols_lo = center - n_extra - n_pixels + 1:center - n_extra
    cols_hi = center + n_extra + 1:center + n_extra + n_pixels
    q1 = @view frame[rows_lo, cols_lo]
    q2 = @view frame[rows_lo, cols_hi]
    q3 = @view frame[rows_hi, cols_hi]
    q4 = @view frame[rows_hi, cols_lo]
    i4q = @view wfs.state.signal_2d[1:n_pixels, :]
    @. i4q = q1 + q2 + q3 + q4
    cutoff = wfs.params.light_ratio * maximum(i4q)
    @. wfs.state.valid_i4q = i4q >= cutoff
    update_pyramid_valid_signal!(wfs)
    update_pyramid_valid_signal_indices!(wfs)
    resize_pyramid_slope_buffers!(wfs)
    return wfs
end

function sample_pyramid_intensity!(wfs::PyramidWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.params.binning
    sub = div(tel.params.resolution, wfs.params.n_subap)
    if size(intensity, 1) % sub != 0
        throw(InvalidConfiguration("pyramid intensity size must be divisible by telescope pixels per subaperture"))
    end
    n_camera = div(size(intensity, 1), sub)
    wfs.state.nominal_detector_resolution = n_camera
    if size(wfs.state.camera_frame) != (n_camera, n_camera)
        wfs.state.camera_frame = similar(wfs.state.camera_frame, n_camera, n_camera)
    end
    frame = wfs.state.camera_frame
    if binning != 1
        if n_camera % binning != 0
            throw(InvalidConfiguration("pyramid binning must evenly divide detector resolution"))
        end
        n_binned = div(n_camera, binning)
        if size(wfs.state.binned_intensity) != (n_binned, n_binned)
            wfs.state.binned_intensity = similar(wfs.state.binned_intensity, n_binned, n_binned)
        end
        bin2d!(wfs.state.binned_intensity, intensity, sub * binning)
        frame = wfs.state.binned_intensity
    else
        bin2d!(wfs.state.camera_frame, intensity, sub)
    end
    resize_pyramid_signal_buffers!(wfs, size(frame, 1))
    return frame
end
