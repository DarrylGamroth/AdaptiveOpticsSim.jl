using Statistics

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
        phase = 2 * π * mod_amp * (x * cos(angle) + y * sin(angle))
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

struct PyramidParams{T<:AbstractFloat,M}
    n_subap::Int
    threshold::T
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
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    RB<:AbstractMatrix{T},
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
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::Kf
    lgs_kernel_tag::UInt
    optical_gain::Vg
    shift_x::NTuple{4,Int}
    shift_y::NTuple{4,Int}
    effective_resolution::Int
end

struct PyramidWFS{M<:SensingMode,P<:PyramidParams,S<:PyramidState} <: AbstractWFS
    params::P
    state::S
end

function PyramidWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1, modulation::Real=2.0,
    modulation_points::Union{Int,Nothing}=nothing, extra_modulation_factor::Int=0,
    old_mask::Bool=false, rooftop::Real=0.0, theta_rotation::Real=0.0, delta_theta::Real=0.0,
    user_modulation_path=nothing, mask_scale::Real=1.0, diffraction_padding::Int=2,
    psf_centering::Bool=true, n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    n_mod = resolve_modulation_points(T(modulation), modulation_points, extra_modulation_factor, user_modulation_path)
    params = PyramidParams{T,typeof(user_modulation_path)}(
        n_subap,
        T(threshold),
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
        typeof(field),
        typeof(intensity),
        typeof(binned_intensity),
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
        fft_plan,
        ifft_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
        optical_gain,
        (0, 0, 0, 0),
        (0, 0, 0, 0),
        pad,
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
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.focal_field)
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.pupil_field)
        wfs.state.lgs_kernel_fft = similar(wfs.state.focal_field, Complex{eltype(wfs.state.focal_field)}, 0, 0)
        wfs.state.lgs_kernel_tag = UInt(0)
        wfs.state.effective_resolution = pad
        build_pyramid_phasor!(wfs.state.phasor)
        build_pyramid_mask!(wfs, tel)
    end
    return wfs
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
    binning = wfs.params.binning
    if binning == 1
        return intensity
    end
    pad = size(intensity, 1)
    if pad % binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide padded resolution"))
    end
    n_binned = div(pad, binning)
    if size(wfs.state.binned_intensity) != (n_binned, n_binned)
        wfs.state.binned_intensity = similar(wfs.state.binned_intensity, n_binned, n_binned)
    end
    bin2d!(wfs.state.binned_intensity, intensity, binning)
    return wfs.state.binned_intensity
end

function measure!(mode::Geometric, wfs::PyramidWFS, tel::Telescope)
    geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask)
    gain = inv(1 + wfs.params.modulation)
    @. wfs.state.slopes = gain * wfs.state.slopes * wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.n_subap
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope)
    throw(InvalidConfiguration("Diffractive PyramidWFS requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::PyramidWFS, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::PyramidWFS, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    n_sub = wfs.params.n_subap
    pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, wfs.state.intensity)
    pyramid_slopes!(wfs, tel, intensity)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    n_sub = wfs.params.n_subap
    pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    pyramid_slopes!(wfs, tel, frame)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    fill!(wfs.state.intensity, zero(eltype(wfs.state.intensity)))
    for src in ast.sources
        pyramid_intensity!(wfs.state.temp, wfs, tel, src)
        wfs.state.intensity .+= wfs.state.temp
    end
    n_sub = wfs.params.n_subap
    intensity = sample_pyramid_intensity!(wfs, wfs.state.intensity)
    pyramid_slopes!(wfs, tel, intensity)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    fill!(wfs.state.intensity, zero(eltype(wfs.state.intensity)))
    for src in ast.sources
        pyramid_intensity!(wfs.state.temp, wfs, tel, src)
        wfs.state.intensity .+= wfs.state.temp
    end
    n_sub = wfs.params.n_subap
    intensity = sample_pyramid_intensity!(wfs, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    pyramid_slopes!(wfs, tel, frame)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function build_pyramid_phasor!(phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    _build_pyramid_phasor!(execution_style(phasor), phasor)
    return phasor
end

function _build_pyramid_phasor!(::ScalarCPUStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    @inbounds for i in 1:n, j in 1:n
        phase = scale * (i + j - 2)
        phasor[i, j] = cis(phase)
    end
    return phasor
end

function _build_pyramid_phasor!(style::AcceleratorStyle, phasor::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    n = size(phasor, 1)
    scale = -T(pi) * (n + 1) / n
    launch_kernel!(style, pyramid_phasor_kernel!, phasor, scale, n; ndrange=size(phasor))
    return phasor
end

function build_pyramid_mask!(wfs::PyramidWFS, tel::Telescope)
    mask = wfs.state.pyramid_mask
    copyto!(mask, host_pyramid_mask(wfs, tel))
    return mask
end

function host_pyramid_mask(wfs::PyramidWFS, tel::Telescope)
    n = size(wfs.state.pyramid_mask, 1)
    T = eltype(wfs.state.slopes)
    host = Matrix{Complex{T}}(undef, n, n)
    if wfs.params.old_mask
        build_pyramid_mask_old_host!(host, wfs, tel)
    else
        build_pyramid_mask_new_host!(host, wfs, tel)
    end
    return host
end

function build_pyramid_mask_new_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    rooftop_pixels = wfs.params.rooftop * wfs.params.diffraction_padding / sqrt(T(2))
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(π)
    if wfs.params.psf_centering
        xvals = range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        xvals = range(-lim, lim; length=n + 1)[1:n]
    end
    r = (T(n_sub) + T(sep)) / 2
    sx = wfs.state.shift_x
    sy = wfs.state.shift_y
    θ = wfs.params.theta_rotation
    cθ = cos(θ)
    sθ = sin(θ)
    @inbounds for i in 1:n, j in 1:n
        x_ = xvals[i]
        y_ = xvals[j]
        x = x_ * cθ - y_ * sθ
        y = y_ * cθ + x_ * sθ
        p1 = x * r + x_ * sx[1] + y * r - y_ * sy[1] + rooftop_pixels
        p2 = -x * r + x_ * sx[2] + y * r - y_ * sy[2]
        p3 = -x * r + x_ * sx[3] - y * r - y_ * sy[3] + rooftop_pixels
        p4 = x * r + x_ * sx[4] - y * r - y_ * sy[4]
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function build_pyramid_mask_old_host!(mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n_tot = size(mask, 1)
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    sx = wfs.state.shift_x
    sy = wfs.state.shift_y
    norma = (T(tel.params.resolution) / T(n_sub)) / 4
    fill!(mask, complex(zero(T), zero(T)))
    if wfs.params.psf_centering
        tip = centered_grid(T, n_tot ÷ 2, true)
        tilt = centered_grid(T, n_tot ÷ 2, true)
        Tip = repeat(tip', n_tot ÷ 2, 1)
        Tilt = repeat(tilt, 1, n_tot ÷ 2)
        Tip .-= mean(Tip)
        Tilt .-= mean(Tilt)
        q = T(n_sub + sep)
        @views begin
            mask[1:n_tot÷2, 1:n_tot÷2] .= cis.((Tip .* (q + sx[1]) .+ Tilt .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[2]) .+ Tilt .* (q - sy[2])) .* norma)
            mask[n_tot÷2+1:end, n_tot÷2+1:end] .= cis.((-Tip .* (q - sx[3]) .- Tilt .* (q + sy[3])) .* norma)
            mask[n_tot÷2+1:end, 1:n_tot÷2] .= cis.((Tip .* (q + sx[4]) .- Tilt .* (q + sy[4])) .* norma)
        end
    else
        d_pix = T(π) / T(n_tot)
        lim_p = T(π)
        lim_m = T(π) - 2 * d_pix
        tip1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip2 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tip3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tip4 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt1 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        tilt2 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt3 = axis_values(T, n_tot ÷ 2 - 1, -lim_m, lim_m; endpoint=false)
        tilt4 = axis_values(T, n_tot ÷ 2 + 1, -lim_p, lim_p; endpoint=true)
        Tip1 = repeat(tip1', length(tilt1), 1); Tilt1 = repeat(tilt1, 1, length(tip1))
        Tip2 = repeat(tip2', length(tilt2), 1); Tilt2 = repeat(tilt2, 1, length(tip2))
        Tip3 = repeat(tip3', length(tilt3), 1); Tilt3 = repeat(tilt3, 1, length(tip3))
        Tip4 = repeat(tip4', length(tilt4), 1); Tilt4 = repeat(tilt4, 1, length(tip4))
        Tip1 .-= mean(Tip1); Tilt1 .-= mean(Tilt1)
        Tip2 .-= mean(Tip2); Tilt2 .-= mean(Tilt2)
        Tip3 .-= mean(Tip3); Tilt3 .-= mean(Tilt3)
        Tip4 .-= mean(Tip4); Tilt4 .-= mean(Tilt4)
        q = T(n_sub + sep)
        @views begin
            mask[1:n_tot÷2+1, 1:n_tot÷2+1] .= cis.((Tip1 .* (q + sx[1]) .+ Tilt1 .* (q - sy[1])) .* norma)
            mask[1:n_tot÷2+1, n_tot÷2:end] .= cis.((-Tip4 .* (q - sx[2]) .+ Tilt4 .* (q - sy[2])) .* norma)
            mask[n_tot÷2:end, n_tot÷2:end] .= cis.((-Tip3 .* (q - sx[3]) .- Tilt3 .* (q + sy[3])) .* norma)
            mask[n_tot÷2:end, 1:n_tot÷2+1] .= cis.((Tip2 .* (q + sx[4]) .- Tilt2 .* (q + sy[4])) .* norma)
        end
    end
    return mask
end

centered_grid(::Type{T}, n::Int, endpoint::Bool) where {T<:AbstractFloat} =
    repeat(axis_values(T, n, -T(π), T(π); endpoint=endpoint), 1, 1)

function axis_values(::Type{T}, n::Int, lo::T, hi::T; endpoint::Bool) where {T<:AbstractFloat}
    if endpoint
        return reshape(collect(range(lo, hi; length=n)), n, 1)
    end
    vals = collect(range(lo, hi; length=n + 1))
    return reshape(vals[1:n], n, 1)
end

function _build_pyramid_mask!(::ScalarCPUStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.params.mask_scale / 2
    pix_per_subap = T(tel.params.resolution) / T(n_sub)
    norma = pix_per_subap
    lim = T(pi)
    x_vals = if wfs.params.psf_centering
        range(-lim * (one(T) - one(T) / T(n)), lim * (one(T) - one(T) / T(n)); length=n)
    else
        range(-lim, lim; length=n, endpoint=false)
    end
    @inbounds for i in 1:n, j in 1:n
        x = x_vals[i]
        y = x_vals[j]
        p1 = x * r + y * r
        p2 = -x * r + y * r
        p3 = -x * r - y * r
        p4 = x * r - y * r
        phase = -max(max(p1, p2), max(p3, p4)) * norma
        mask[i, j] = cis(phase)
    end
    return mask
end

function _build_pyramid_mask!(style::AcceleratorStyle, mask::AbstractMatrix{Complex{T}}, wfs::PyramidWFS, tel::Telescope) where {T<:AbstractFloat}
    n = size(mask, 1)
    n_sub = wfs.params.n_subap
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    r = (T(n_sub) + T(sep)) * wfs.params.mask_scale / 2
    norma = T(tel.params.resolution) / T(n_sub)
    lim = T(pi)
    start = wfs.params.psf_centering ? -lim * (one(T) - one(T) / T(n)) : -lim
    step = T(2) * lim / T(n)
    launch_kernel!(style, pyramid_mask_kernel!, mask, r, norma, start, step, n; ndrange=size(mask))
    return mask
end

function pyramid_intensity_core!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::AbstractSource) where {T<:AbstractFloat}
    prepare_pyramid_sampling!(wfs, tel)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    phase_scale = (2 * pi) / wavelength(src)

    fill!(out, zero(T))
    for p in 1:wfs.params.modulation_points
        fill!(wfs.state.field, zero(eltype(wfs.state.field)))
        @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = tel.state.pupil *
            wfs.state.modulation_phases[:, :, p] * cis(phase_scale * tel.state.opd)
        copyto!(wfs.state.focal_field, wfs.state.field)
        if wfs.params.psf_centering
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.phasor
            mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
            @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.pyramid_mask
            copyto!(wfs.state.pupil_field, wfs.state.focal_field)
            mul!(wfs.state.pupil_field, wfs.state.ifft_plan, wfs.state.pupil_field)
        else
            mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
            fftshift2d!(wfs.state.pupil_field, wfs.state.focal_field)
            @. wfs.state.pupil_field = wfs.state.pupil_field * wfs.state.pyramid_mask
            mul!(wfs.state.pupil_field, wfs.state.ifft_plan, wfs.state.pupil_field)
        end
        @. wfs.state.temp = abs2(wfs.state.pupil_field)
        out .+= wfs.state.temp
    end
    out ./= wfs.params.modulation_points

    return out
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::AbstractSource) where {T<:AbstractFloat}
    return pyramid_intensity_core!(out, wfs, tel, src)
end

function pyramid_intensity!(out::AbstractMatrix{T}, wfs::PyramidWFS, tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    pyramid_intensity_core!(out, wfs, tel, src)
    apply_lgs_elongation!(lgs_profile(src), out, wfs, tel, src)
    return out
end

function apply_lgs_elongation!(::LGSProfileNone, out::AbstractMatrix{T}, wfs::PyramidWFS,
    ::Telescope, src::LGSSource) where {T<:AbstractFloat}
    wfs.state.elongation_kernel = apply_elongation!(
        out,
        lgs_elongation_factor(src),
        wfs.state.scratch,
        wfs.state.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, out::AbstractMatrix{T}, wfs::PyramidWFS,
    tel::Telescope, src::LGSSource) where {T<:AbstractFloat}
    ensure_lgs_kernel!(wfs, tel, src)
    apply_lgs_convolution!(
        out,
        wfs.state.lgs_kernel_fft,
        wfs.state.focal_field,
        wfs.state.fft_plan,
        wfs.state.pupil_field,
        wfs.state.ifft_plan,
    )
    return wfs
end

function ensure_lgs_kernel!(wfs::PyramidWFS, tel::Telescope, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    pad = size(wfs.state.intensity, 1)
    tag = objectid(na_profile) ⊻ hash(src.params.laser_coordinates) ⊻ hash(src.params.fwhm_spot_up) ⊻
        hash(pad) ⊻ hash(wfs.params.n_subap)
    if size(wfs.state.lgs_kernel_fft, 1) == pad && wfs.state.lgs_kernel_tag == tag
        return wfs
    end
    padding = wfs.state.effective_resolution / tel.params.resolution
    pixel_scale = lgs_pixel_scale(tel.params.diameter, padding, wavelength(src))
    wfs.state.lgs_kernel_fft = lgs_average_kernel_fft(
        tel,
        src,
        pad,
        wfs.params.n_subap,
        pixel_scale,
        wfs.state.focal_field,
        wfs.state.fft_plan,
    )
    wfs.state.lgs_kernel_tag = tag
    return wfs
end

function build_modulation_phases!(wfs::PyramidWFS, tel::Telescope)
    phases = wfs.state.modulation_phases
    copyto!(phases, host_modulation_phases(eltype(wfs.state.slopes), tel, wfs.params.modulation,
        wfs.params.modulation_points, wfs.params.delta_theta, wfs.params.user_modulation_path))
    return phases
end

function host_modulation_phases(::Type{T}, tel::Telescope, modulation::T, n_pts::Int,
    delta_theta::T, user_modulation_path) where {T<:AbstractFloat}
    n = tel.params.resolution
    phases = Array{Complex{T}}(undef, n, n, n_pts)
    if n_pts == 1 && modulation == 0 && user_modulation_path === nothing
        fill!(phases, complex(one(T), zero(T)))
        return phases
    end
    tip_vals = collect(range(-T(π), T(π); length=n))
    Tilt = reshape(tip_vals, n, 1)
    Tip = reshape(tip_vals, 1, n)
    pupil = tel.state.pupil
    if user_modulation_path === nothing
        @inbounds for p in 1:n_pts
            θ = delta_theta + T(2 * π * (p - 1) / n_pts)
            mx = modulation * cos(θ)
            my = modulation * sin(θ)
            for i in 1:n, j in 1:n
                phase = pupil[i, j] ? (mx * Tip[j] + my * Tilt[i]) : zero(T)
                phases[i, j, p] = cis(phase)
            end
        end
    else
        @inbounds for p in 1:n_pts
            mx = T(user_modulation_path[p][1])
            my = T(user_modulation_path[p][2])
            for i in 1:n, j in 1:n
                phase = pupil[i, j] ? (mx * Tip[j] + my * Tilt[i]) : zero(T)
                phases[i, j, p] = cis(phase)
            end
        end
    end
    return phases
end

function _build_modulation_phases!(::ScalarCPUStyle, phases::AbstractArray{Complex{T},3}, mod, center, n_pts::Int, n::Int) where {T<:AbstractFloat}
    x_coords = Vector{T}(undef, n)
    for i in 1:n
        x_coords[i] = (i - center) / n
    end
    @inbounds for p in 1:n_pts
        angle = 2 * pi * (p - 1) / n_pts
        c = cos(angle)
        s = sin(angle)
        for i in 1:n, j in 1:n
            phase = 2 * pi * mod * (x_coords[i] * c + x_coords[j] * s)
            phases[i, j, p] = cis(phase)
        end
    end
    return phases
end

function _build_modulation_phases!(style::AcceleratorStyle, phases::AbstractArray{Complex{T},3}, mod, center, n_pts::Int, n::Int) where {T<:AbstractFloat}
    launch_kernel!(style, modulation_phases_kernel!, phases, T(mod), T(center), n_pts, n; ndrange=size(phases))
    return phases
end

function pyramid_slopes!(wfs::PyramidWFS, tel::Telescope)
    return pyramid_slopes!(wfs, tel, wfs.state.intensity)
end

function pyramid_slopes!(wfs::PyramidWFS, tel::Telescope, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    pad = size(intensity, 1)
    n_sub = wfs.params.n_subap
    if wfs.state.effective_resolution % pad != 0
        throw(InvalidConfiguration("pyramid intensity size must divide padded resolution"))
    end
    binning = div(wfs.state.effective_resolution, pad)
    if tel.params.resolution % binning != 0
        throw(InvalidConfiguration("pyramid binning must evenly divide telescope resolution"))
    end
    pupil_size = div(tel.params.resolution, binning)
    sep = wfs.params.n_pix_separation === nothing ? 0 : wfs.params.n_pix_separation
    pix_per_subap = tel.params.resolution / n_sub
    sep_pix = round(Int, sep * pix_per_subap / binning)
    edge_pix = div(pad - 2 * pupil_size - sep_pix, 2)
    if edge_pix < 0
        throw(InvalidConfiguration("pyramid padding is too small for pupil extraction"))
    end
    sub = div(pupil_size, n_sub)
    if sub < 1
        throw(InvalidConfiguration("pyramid pupil size is too small for n_subap"))
    end
    ox1 = edge_pix
    oy1 = edge_pix
    ox2 = edge_pix
    oy2 = edge_pix + pupil_size + sep_pix
    ox3 = edge_pix + pupil_size + sep_pix
    oy3 = edge_pix
    ox4 = edge_pix + pupil_size + sep_pix
    oy4 = edge_pix + pupil_size + sep_pix
    _pyramid_slopes!(execution_style(wfs.state.slopes), wfs.state.slopes, intensity, wfs.state.valid_mask, sub, n_sub, pad,
        n_sub * n_sub, ox1, oy1, ox2, oy2, ox3, oy3, ox4, oy4, wfs.state.shift_x, wfs.state.shift_y)
    return wfs.state.slopes
end

function _pyramid_slopes!(::ScalarCPUStyle, slopes::AbstractVector, intensity::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, ::Int, offset::Int, ox1::Int, oy1::Int, ox2::Int, oy2::Int, ox3::Int, oy3::Int, ox4::Int, oy4::Int,
    shift_x::NTuple{4,Int}, shift_y::NTuple{4,Int}) where {T<:AbstractFloat}
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if valid_mask[i, j]
            q1 = quadrant_sum(intensity, xs, ys, sub, ox1, oy1, shift_x[1], shift_y[1])
            q2 = quadrant_sum(intensity, xs, ys, sub, ox2, oy2, shift_x[2], shift_y[2])
            q3 = quadrant_sum(intensity, xs, ys, sub, ox3, oy3, shift_x[3], shift_y[3])
            q4 = quadrant_sum(intensity, xs, ys, sub, ox4, oy4, shift_x[4], shift_y[4])
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
        idx += 1
    end
    return slopes
end

function _pyramid_slopes!(style::AcceleratorStyle, slopes::AbstractVector, intensity::AbstractMatrix, valid_mask::AbstractMatrix{Bool},
    sub::Int, n_sub::Int, pad::Int, offset::Int, ox1::Int, oy1::Int, ox2::Int, oy2::Int, ox3::Int, oy3::Int, ox4::Int, oy4::Int,
    shift_x::NTuple{4,Int}, shift_y::NTuple{4,Int})
    launch_kernel!(style, pyramid_slopes_kernel!, slopes, intensity, valid_mask, sub, n_sub, pad, offset,
        ox1, oy1, ox2, oy2, ox3, oy3, ox4, oy4,
        shift_x[1], shift_y[1], shift_x[2], shift_y[2], shift_x[3], shift_y[3], shift_x[4], shift_y[4];
        ndrange=(n_sub, n_sub))
    return slopes
end

function quadrant_sum(intensity::AbstractMatrix{T}, xs::Int, ys::Int, sub::Int,
    ox::Int, oy::Int, sx::Int, sy::Int) where {T<:AbstractFloat}
    acc = zero(T)
    n = size(intensity, 1)
    @inbounds for i in 0:(sub - 1), j in 0:(sub - 1)
        xi = ox + xs + i - 1 + sx
        yj = oy + ys + j - 1 + sy
        if 1 <= xi <= n && 1 <= yj <= n
            acc += intensity[xi, yj]
        end
    end
    return acc
end

function apply_shift_wfs!(wfs::PyramidWFS; sx, sy)
    if sx isa Real && sy isa Real
        shift_x = (round(Int, sx), round(Int, sx), round(Int, sx), round(Int, sx))
        shift_y = (round(Int, sy), round(Int, sy), round(Int, sy), round(Int, sy))
    else
        if length(sx) != 4 || length(sy) != 4
            throw(InvalidConfiguration("pyramid shift must have 4 elements"))
        end
        shift_x = (round(Int, sx[1]), round(Int, sx[2]), round(Int, sx[3]), round(Int, sx[4]))
        shift_y = (round(Int, sy[1]), round(Int, sy[2]), round(Int, sy[3]), round(Int, sy[4]))
    end
    wfs.state.shift_x = shift_x
    wfs.state.shift_y = shift_y
    return wfs
end

function set_optical_gain!(wfs::PyramidWFS, gain::Real)
    fill!(wfs.state.optical_gain, gain)
    return wfs
end

function set_optical_gain!(wfs::PyramidWFS, gain::AbstractVector)
    if length(gain) != length(wfs.state.optical_gain)
        throw(InvalidConfiguration("optical_gain length must match slope vector"))
    end
    copyto!(wfs.state.optical_gain, gain)
    return wfs
end
