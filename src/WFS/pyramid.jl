using Statistics

struct PyramidParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    modulation::T
    modulation_points::Int
    mask_scale::T
    diffraction_padding::Int
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
    modulation_points::Int=1, mask_scale::Real=1.0, diffraction_padding::Int=2,
    n_pix_separation=nothing, n_pix_edge=nothing, binning::Int=1,
    mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    if modulation_points < 1
        throw(InvalidConfiguration("modulation_points must be >= 1"))
    end
    if binning < 1
        throw(InvalidConfiguration("binning must be >= 1"))
    end
    params = PyramidParams{T}(n_subap, T(threshold), T(modulation), modulation_points, T(mask_scale),
        diffraction_padding, n_pix_separation, n_pix_edge, binning)
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
    modulation_phases = backend{Complex{T}}(undef, tel.params.resolution, tel.params.resolution, modulation_points)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    scratch = similar(intensity)
    binned_intensity = similar(intensity)
    fft_plan = plan_fft!(focal_field)
    ifft_plan = plan_ifft!(pupil_field)
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
    build_pyramid_mask!(wfs.state.pyramid_mask, wfs.params.mask_scale)
    build_modulation_phases!(wfs, tel)
    return wfs
end

sensing_mode(::PyramidWFS{M}) where {M} = M()

function update_valid_mask!(wfs::PyramidWFS, tel::Telescope)
    Base.require_one_based_indexing(tel.state.pupil)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        subap_pupil = @view tel.state.pupil[xs:xe, ys:ye]
        wfs.state.valid_mask[i, j] = mean(subap_pupil) > wfs.params.threshold
    end
    return wfs
end

function ensure_pyramid_buffers!(wfs::PyramidWFS, pad::Int)
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.focal_field = similar(wfs.state.focal_field, pad, pad)
        wfs.state.pupil_field = similar(wfs.state.pupil_field, pad, pad)
        wfs.state.pyramid_mask = similar(wfs.state.pyramid_mask, pad, pad)
        wfs.state.intensity = similar(wfs.state.intensity, pad, pad)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.scratch = similar(wfs.state.scratch, pad, pad)
        wfs.state.binned_intensity = similar(wfs.state.binned_intensity, pad, pad)
        wfs.state.fft_plan = plan_fft!(wfs.state.focal_field)
        wfs.state.ifft_plan = plan_ifft!(wfs.state.pupil_field)
        wfs.state.lgs_kernel_fft = similar(wfs.state.focal_field, Complex{eltype(wfs.state.focal_field)}, 0, 0)
        wfs.state.lgs_kernel_tag = UInt(0)
        wfs.state.effective_resolution = pad
        build_pyramid_mask!(wfs.state.pyramid_mask, wfs.params.mask_scale)
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
    ensure_pyramid_buffers!(wfs, pad)
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
    Base.require_one_based_indexing(tel.state.opd)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            sx = 0.0
            sy = 0.0
            count_x = 0
            count_y = 0
            for x in xs:(xe - 1), y in ys:ye
                sx += tel.state.opd[x + 1, y] - tel.state.opd[x, y]
                count_x += 1
            end
            for x in xs:xe, y in ys:(ye - 1)
                sy += tel.state.opd[x, y + 1] - tel.state.opd[x, y]
                count_y += 1
            end
            gain = 1 / (1 + wfs.params.modulation)
            wfs.state.slopes[idx] = gain * sx / max(count_x, 1)
            wfs.state.slopes[idx + n_sub * n_sub] = gain * sy / max(count_y, 1)
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    @. wfs.state.slopes *= wfs.state.optical_gain
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
    pyramid_slopes!(wfs, n_sub, intensity)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    n_sub = wfs.params.n_subap
    pyramid_intensity!(wfs.state.intensity, wfs, tel, src)
    intensity = sample_pyramid_intensity!(wfs, wfs.state.intensity)
    frame = capture!(det, intensity; rng=rng)
    pyramid_slopes!(wfs, n_sub, frame)
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
    pyramid_slopes!(wfs, n_sub, intensity)
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
    pyramid_slopes!(wfs, n_sub, frame)
    @. wfs.state.slopes *= wfs.state.optical_gain
    return wfs.state.slopes
end

function build_pyramid_mask!(mask::AbstractMatrix{Complex{T}}, modulation::T) where {T<:AbstractFloat}
    n = size(mask, 1)
    center = (n + 1) / 2
    scale = modulation / center
    @inbounds for i in 1:n, j in 1:n
        x = (i - center) * scale
        y = (j - center) * scale
        phase = (sign(x) * x + sign(y) * y)
        mask[i, j] = cis(phase)
    end
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
        mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
        @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.pyramid_mask
        copyto!(wfs.state.pupil_field, wfs.state.focal_field)
        mul!(wfs.state.pupil_field, wfs.state.ifft_plan, wfs.state.pupil_field)
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
    if !(wfs.state.intensity isa Array)
        throw(InvalidConfiguration("LGS Na-profile kernels currently require Array backend"))
    end
    pad = size(wfs.state.intensity, 1)
    tag = objectid(na_profile) ⊻ hash(src.params.laser_coordinates) ⊻ hash(src.params.fwhm_spot_up) ⊻
        hash(pad) ⊻ hash(wfs.params.n_subap)
    if size(wfs.state.lgs_kernel_fft, 1) == pad && wfs.state.lgs_kernel_tag == tag
        return wfs
    end
    pixel_scale = lgs_pixel_scale(tel.params.diameter, wfs.params.diffraction_padding, wavelength(src))
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
    n = tel.params.resolution
    center = (n + 1) / 2
    mod = wfs.params.modulation
    n_pts = wfs.params.modulation_points
    phases = wfs.state.modulation_phases
    if n_pts == 1 || mod == 0
        fill!(phases, one(eltype(phases)))
        return phases
    end
    T = eltype(wfs.state.slopes)
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

function pyramid_slopes!(wfs::PyramidWFS, n_sub::Int)
    return pyramid_slopes!(wfs, n_sub, wfs.state.intensity)
end

function pyramid_slopes!(wfs::PyramidWFS, n_sub::Int, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    pad = size(intensity, 1)
    half = div(pad, 2)
    pupil_size = half
    sub = div(pupil_size, n_sub)
    if sub < 1
        throw(InvalidConfiguration("pyramid pupil size is too small for n_subap"))
    end
    shift_x = wfs.state.shift_x
    shift_y = wfs.state.shift_y
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        if wfs.state.valid_mask[i, j]
            q1 = quadrant_sum(intensity, xs, ys, sub, 1, 1, shift_x[1], shift_y[1])
            q2 = quadrant_sum(intensity, xs, ys, sub, 1, half + 1, shift_x[2], shift_y[2])
            q3 = quadrant_sum(intensity, xs, ys, sub, half + 1, 1, shift_x[3], shift_y[3])
            q4 = quadrant_sum(intensity, xs, ys, sub, half + 1, half + 1, shift_x[4], shift_y[4])
            left = q1 + q3
            right = q2 + q4
            bottom = q1 + q2
            top = q3 + q4
            total = left + right
            if total <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                wfs.state.slopes[idx] = (right - left) / total
                wfs.state.slopes[idx + n_sub * n_sub] = (top - bottom) / total
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    return wfs.state.slopes
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
