using Statistics

struct ShackHartmannParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    threshold_cog::T
    threshold_convolution::T
    half_pixel_shift::Bool
    diffraction_padding::Int
    pixel_scale::Union{T,Nothing}
    n_pix_subap::Union{Int,Nothing}
    shannon_sampling::Bool
end

mutable struct ShackHartmannState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    RB<:AbstractMatrix{T},
    RS<:AbstractMatrix{T},
    RC<:AbstractArray{T,3},
    P,
    Pi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3}}
    valid_mask::A
    slopes::V
    field::C
    phasor::C
    fft_buffer::C
    intensity::R
    temp::R
    bin_buffer::RB
    spot::RS
    spot_cube::RC
    fft_plan::P
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::KF
    lgs_kernel_tag::UInt
    effective_padding::Int
    binning_pixel_scale::Int
    sampled_n_pix_subap::Int
    phasor_ratio::T
    reference_signal_2d::RB
    slopes_units::T
    calibrated::Bool
    calibration_wavelength::T
end

struct ShackHartmann{M<:SensingMode,P<:ShackHartmannParams,S<:ShackHartmannState} <: AbstractWFS
    params::P
    state::S
end

function ShackHartmann(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    threshold_cog::Real=0.01, threshold_convolution::Real=0.05, half_pixel_shift::Bool=false,
    diffraction_padding::Int=2, pixel_scale=nothing, n_pix_subap=nothing,
    shannon_sampling::Bool=true, mode::SensingMode=Geometric(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    pixel_scale_t = pixel_scale === nothing ? nothing : T(pixel_scale)
    params = ShackHartmannParams{T}(n_subap, T(threshold), T(threshold_cog), T(threshold_convolution), half_pixel_shift, diffraction_padding,
        pixel_scale_t, n_pix_subap, shannon_sampling)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    sub = div(tel.params.resolution, n_subap)
    pad = max(sub, sub * diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    phasor = similar(field)
    fft_buffer = similar(field)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    bin_buffer = backend{T}(undef, sub, sub)
    spot = similar(bin_buffer)
    spot_cube = backend{T}(undef, n_subap * n_subap, sub, sub)
    fft_plan = plan_fft_backend!(fft_buffer)
    ifft_plan = plan_ifft_backend!(fft_buffer)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0, 0)
    state = ShackHartmannState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(field),
        typeof(intensity),
        typeof(bin_buffer),
        typeof(spot),
        typeof(spot_cube),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
    }(
        valid_mask,
        slopes,
        field,
        phasor,
        fft_buffer,
        intensity,
        temp,
        bin_buffer,
        spot,
        spot_cube,
        fft_plan,
        ifft_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
        diffraction_padding,
        1,
        sub,
        T(NaN),
        backend{T}(undef, 2 * n_subap, n_subap),
        one(T),
        false,
        zero(T),
    )
    wfs = ShackHartmann{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ShackHartmann{M}) where {M} = M()

function update_valid_mask!(wfs::ShackHartmann, tel::Telescope)
    set_valid_subapertures!(wfs.state.valid_mask, tel.state.pupil, wfs.params.threshold)
    return wfs
end

function ensure_sh_buffers!(wfs::ShackHartmann, pad::Int)
    if size(wfs.state.field) != (pad, pad)
        wfs.state.field = similar(wfs.state.field, pad, pad)
        wfs.state.phasor = similar(wfs.state.phasor, pad, pad)
        wfs.state.fft_buffer = similar(wfs.state.fft_buffer, pad, pad)
        wfs.state.intensity = similar(wfs.state.intensity, pad, pad)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.fft_buffer)
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.fft_buffer)
        wfs.state.phasor_ratio = eltype(wfs.state.slopes)(NaN)
        wfs.state.calibrated = false
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
        wfs.state.spot_cube = similar(wfs.state.spot_cube, eltype(wfs.state.spot_cube),
            wfs.params.n_subap * wfs.params.n_subap, n_pix_subap, n_pix_subap)
        wfs.state.sampled_n_pix_subap = n_pix_subap
    end

    wfs.state.binning_pixel_scale = binning_pixel_scale
    T = eltype(wfs.state.slopes)
    half_shift_ratio = wfs.params.half_pixel_shift ? T(binning_pixel_scale) : zero(T)
    build_sh_phasor!(wfs, half_shift_ratio)
    return wfs
end

function sample_spot!(wfs::ShackHartmann, intensity::AbstractMatrix{T}) where {T<:AbstractFloat}
    binning = wfs.state.binning_pixel_scale
    spot_in = intensity
    if binning > 1
        pad = size(intensity, 1)
        if pad % binning != 0
            throw(InvalidConfiguration("lenslet sampling is not divisible by binning_pixel_scale"))
        end
        n_binned = div(pad, binning)
        if size(wfs.state.bin_buffer) != (n_binned, n_binned)
            wfs.state.bin_buffer = similar(wfs.state.bin_buffer, n_binned, n_binned)
        end
        bin2d!(wfs.state.bin_buffer, intensity, binning)
        spot_in = wfs.state.bin_buffer
    end
    center_resize2d!(wfs.state.spot, spot_in)
    return wfs.state.spot
end

function measure!(mode::Geometric, wfs::ShackHartmann, tel::Telescope)
    geometric_slopes!(wfs.state.slopes, tel.state.opd, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure!(::Geometric, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    return measure!(Geometric(), wfs, tel)
end

function measure!(::Geometric, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    slopes = measure!(Geometric(), wfs, tel)
    n_sub = wfs.params.n_subap
    factor = lgs_elongation_factor(src)
    @views slopes[n_sub * n_sub + 1:end] .*= factor
    return slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope)
    throw(InvalidConfiguration("Diffractive ShackHartmann requires a source; call measure!(wfs, tel, src)."))
end

function measure!(wfs::ShackHartmann, tel::Telescope)
    return measure!(sensing_mode(wfs), wfs, tel)
end

function measure!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    return measure!(sensing_mode(wfs), wfs, tel, src)
end

function measure!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, src, det; rng=rng)
end

function measure!(wfs::ShackHartmann, tel::Telescope, ast::Asterism, det::AbstractDetector;
    rng::AbstractRNG=Random.default_rng())
    return measure!(sensing_mode(wfs), wfs, tel, ast, det; rng=rng)
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    Base.require_one_based_indexing(tel.state.opd)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src)
    sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
    subtract_reference_and_scale!(wfs)
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    peak = sampled_spots_peak!(wfs, tel, src, det, rng)
    sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
    subtract_reference_and_scale!(wfs)
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            total_sum = zero(eltype(wfs.state.slopes))
            sx_sum = zero(eltype(wfs.state.slopes))
            sy_sum = zero(eltype(wfs.state.slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
            end
            if total_sum <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                wfs.state.slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_2d[idx]) / wfs.state.slopes_units
                wfs.state.slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_2d[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, ast::Asterism,
    det::AbstractDetector; rng::AbstractRNG=Random.default_rng())
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            total_sum = zero(eltype(wfs.state.slopes))
            sx_sum = zero(eltype(wfs.state.slopes))
            sy_sum = zero(eltype(wfs.state.slopes))
            for src in ast.sources
                total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx, det, rng)
                total_sum += total
                sx_sum += sx
                sy_sum += sy
            end
            if total_sum <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                wfs.state.slopes[idx] = ((sy_sum / nsrc) - wfs.state.reference_signal_2d[idx]) / wfs.state.slopes_units
                wfs.state.slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - wfs.state.reference_signal_2d[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    return wfs.state.slopes
end

function compute_intensity!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int)
    opd_to_cycles = eltype(wfs.state.intensity)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.intensity)(photon_flux(src) * tel.params.sampling_time *
        (tel.params.diameter / tel.params.resolution)^2))
    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+sub, oy+1:oy+sub] =
        amp_scale * tel.state.pupil[xs:xe, ys:ye] * cispi(opd_to_cycles * tel.state.opd[xs:xe, ys:ye])
    @. wfs.state.field *= wfs.state.phasor
    copyto!(wfs.state.fft_buffer, wfs.state.field)
    mul!(wfs.state.fft_buffer, wfs.state.fft_plan, wfs.state.fft_buffer)
    @. wfs.state.intensity = abs2(wfs.state.fft_buffer)
    return wfs.state.intensity
end

@inline function centroid_from_intensity!(intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    return centroid_from_intensity!(execution_style(intensity), intensity, threshold)
end

@inline function centroid_from_intensity!(::ScalarCPUStyle, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    peak = maximum(intensity)
    if peak <= 0
        return zero(T), zero(T), zero(T)
    end
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), intensity, threshold * peak)
end

@inline function centroid_from_intensity_cutoff!(intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    return centroid_from_intensity_cutoff!(execution_style(intensity), intensity, cutoff)
end

@inline function centroid_from_intensity_cutoff!(::ScalarCPUStyle, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    total = zero(T)
    sx = zero(T)
    sy = zero(T)
    n1 = size(intensity, 1)
    n2 = size(intensity, 2)
    @inbounds for x in 1:n1, y in 1:n2
        val = intensity[x, y]
        if val < cutoff
            intensity[x, y] = zero(T)
        else
            total += val
            sx += T(x - 1) * val
            sy += T(y - 1) * val
        end
    end
    if total <= 0
        return zero(T), zero(T), zero(T)
    end
    return total, sx / total, sy / total
end

@inline function centroid_from_intensity!(::AcceleratorStyle, intensity::AbstractMatrix{T}, threshold::T) where {T<:AbstractFloat}
    host_intensity = Array(intensity)
    return centroid_from_intensity!(ScalarCPUStyle(), host_intensity, threshold)
end

@inline function centroid_from_intensity_cutoff!(::AcceleratorStyle, intensity::AbstractMatrix{T}, cutoff::T) where {T<:AbstractFloat}
    host_intensity = Array(intensity)
    return centroid_from_intensity_cutoff!(ScalarCPUStyle(), host_intensity, cutoff)
end

@inline function sh_spot_view(wfs::ShackHartmann, idx::Int)
    return @view wfs.state.spot_cube[idx, :, :]
end

host_mask_view(::ScalarCPUStyle, valid_mask::AbstractMatrix{Bool}) = valid_mask
host_mask_view(::AcceleratorStyle, valid_mask::AbstractMatrix{Bool}) = Array(valid_mask)

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    valid_mask = host_mask_view(execution_style(wfs.state.valid_mask), wfs.state.valid_mask)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            sample_spot!(wfs, wfs.state.intensity)
            copyto!(spot_view, wfs.state.spot)
            peak = max(peak, maximum(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    valid_mask = host_mask_view(execution_style(wfs.state.valid_mask), wfs.state.valid_mask)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
            sample_spot!(wfs, wfs.state.intensity)
            copyto!(spot_view, wfs.state.spot)
            peak = max(peak, maximum(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    valid_mask = host_mask_view(execution_style(wfs.state.valid_mask), wfs.state.valid_mask)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            sample_spot!(wfs, wfs.state.intensity)
            frame = capture!(det, wfs.state.spot; rng=rng)
            copyto!(spot_view, frame)
            peak = max(peak, maximum(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    valid_mask = host_mask_view(execution_style(wfs.state.valid_mask), wfs.state.valid_mask)
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if valid_mask[i, j]
            xs = (i - 1) * sub + 1
            ys = (j - 1) * sub + 1
            xe = min(i * sub, n)
            ye = min(j * sub, n)
            compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
            apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
            sample_spot!(wfs, wfs.state.intensity)
            frame = capture!(det, wfs.state.spot; rng=rng)
            copyto!(spot_view, frame)
            peak = max(peak, maximum(spot_view))
        else
            fill!(spot_view, zero(eltype(spot_view)))
        end
        idx += 1
    end
    return peak
end

function sh_signal_from_spots!(wfs::ShackHartmann, cutoff::T) where {T<:AbstractFloat}
    return sh_signal_from_spots!(execution_style(wfs.state.slopes), wfs, cutoff)
end

function sh_signal_from_spots!(::ScalarCPUStyle, wfs::ShackHartmann, cutoff::T) where {T<:AbstractFloat}
    n_sub = wfs.params.n_subap
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if wfs.state.valid_mask[i, j]
            total, sx, sy = centroid_from_intensity_cutoff!(sh_spot_view(wfs, idx), cutoff)
            if total <= 0
                wfs.state.slopes[idx] = zero(T)
                wfs.state.slopes[idx + n_sub * n_sub] = zero(T)
            else
                wfs.state.slopes[idx] = sy
                wfs.state.slopes[idx + n_sub * n_sub] = sx
            end
        else
            wfs.state.slopes[idx] = zero(T)
            wfs.state.slopes[idx + n_sub * n_sub] = zero(T)
        end
        idx += 1
    end
    return wfs.state.slopes
end

function sh_signal_from_spots!(::AcceleratorStyle, wfs::ShackHartmann, cutoff::T) where {T<:AbstractFloat}
    n_sub = wfs.params.n_subap
    offset = n_sub * n_sub
    host_valid_mask = Array(wfs.state.valid_mask)
    host_spots = Array(wfs.state.spot_cube)
    host_slopes = Array(wfs.state.slopes)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if host_valid_mask[i, j]
            spot = @view host_spots[idx, :, :]
            total, sx, sy = centroid_from_intensity_cutoff!(ScalarCPUStyle(), spot, cutoff)
            if total <= 0
                host_slopes[idx] = zero(T)
                host_slopes[idx + offset] = zero(T)
            else
                host_slopes[idx] = sy
                host_slopes[idx + offset] = sx
            end
        else
            host_slopes[idx] = zero(T)
            host_slopes[idx + offset] = zero(T)
        end
        idx += 1
    end
    copyto!(wfs.state.slopes, host_slopes)
    return wfs.state.slopes
end

function sh_signal_from_spots!(wfs::ShackHartmann, peak::T, threshold::T) where {T<:AbstractFloat}
    cutoff = peak <= 0 ? zero(T) : threshold * peak
    return sh_signal_from_spots!(wfs, cutoff)
end

function mean_valid_signal(signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return mean_valid_signal(execution_style(signal), signal, valid_mask)
end

function mean_valid_signal(::ScalarCPUStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    acc = zero(T)
    count = 0
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if valid_mask[i, j]
            acc += signal[idx]
            acc += signal[idx + offset]
            count += 2
        end
        idx += 1
    end
    return count == 0 ? one(T) : acc / count
end

function mean_valid_signal(::AcceleratorStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return mean_valid_signal(ScalarCPUStyle(), Array(signal), Array(valid_mask))
end

function subtract_reference_and_scale!(wfs::ShackHartmann)
    inv_units = inv(wfs.state.slopes_units)
    @. wfs.state.slopes = (wfs.state.slopes - wfs.state.reference_signal_2d) * inv_units
    return wfs.state.slopes
end

function subtract_reference!(wfs::ShackHartmann)
    @. wfs.state.slopes = wfs.state.slopes - wfs.state.reference_signal_2d
    return wfs.state.slopes
end

@kernel function calibration_ramp_kernel!(opd, scale, step, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        xi = -π + (i - 1) * step
        xj = -π + (j - 1) * step
        @inbounds opd[i, j] = (xj + xi) * scale
    end
end

function fill_calibration_ramp!(::ScalarCPUStyle, opd::AbstractMatrix{T}, scale::T, n::Int) where {T<:AbstractFloat}
    step = T(2π) / n
    xvals = range(-T(π); step=step, length=n)
    @inbounds for j in 1:n, i in 1:n
        opd[i, j] = (xvals[j] + xvals[i]) * scale
    end
    return opd
end

function fill_calibration_ramp!(style::AcceleratorStyle, opd::AbstractMatrix{T}, scale::T, n::Int) where {T<:AbstractFloat}
    step = T(2π) / n
    launch_kernel!(style, calibration_ramp_kernel!, opd, scale, step, n; ndrange=size(opd))
    return opd
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, ::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs, wfs.state.intensity)
    return centroid_from_intensity!(spot, wfs.params.threshold_cog)
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
    spot = sample_spot!(wfs, wfs.state.intensity)
    return centroid_from_intensity!(spot, wfs.params.threshold_cog)
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, ::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    spot = sample_spot!(wfs, wfs.state.intensity)
    frame = capture!(det, spot; rng=rng)
    return centroid_from_intensity!(frame, wfs.params.threshold_cog)
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, ::Int, idx::Int,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
    spot = sample_spot!(wfs, wfs.state.intensity)
    frame = capture!(det, spot; rng=rng)
    return centroid_from_intensity!(frame, wfs.params.threshold_cog)
end

function sh_reference_signal!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    peak = sampled_spots_peak!(wfs, tel, src)
    sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
    copyto!(wfs.state.reference_signal_2d, wfs.state.slopes)
    return wfs
end

function ensure_sh_calibration!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    λ = eltype(wfs.state.slopes)(wavelength(src))
    if wfs.state.calibrated && wfs.state.calibration_wavelength == λ
        return wfs
    end
    opd_saved = copy(tel.state.opd)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    sh_reference_signal!(wfs, tel, src)
    copyto!(tel.state.opd, opd_saved)

    T = eltype(wfs.state.slopes)
    n = tel.params.resolution
    pixel_scale_init = sh_pixel_scale_init(tel.params.diameter / wfs.params.n_subap, wfs.state.effective_padding, src)
    pixel_scale = T(wfs.state.binning_pixel_scale) * pixel_scale_init
    rad2arcsec = T(180 * 3600 / π)
    scale = T(T(tel.params.diameter) * pixel_scale / (T(2π) * rad2arcsec))
    fill_calibration_ramp!(execution_style(tel.state.opd), tel.state.opd, scale, n)
    peak = sampled_spots_peak!(wfs, tel, src)
    sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
    subtract_reference!(wfs)
    wfs.state.slopes_units = mean_valid_signal(wfs.state.slopes, wfs.state.valid_mask)
    if !isfinite(wfs.state.slopes_units) || wfs.state.slopes_units == zero(T)
        wfs.state.slopes_units = one(T)
    end
    copyto!(tel.state.opd, opd_saved)
    wfs.state.calibrated = true
    wfs.state.calibration_wavelength = λ
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNone, wfs::ShackHartmann, ::Telescope, src::LGSSource, ::Int)
    wfs.state.elongation_kernel = apply_elongation!(
        wfs.state.intensity,
        lgs_elongation_factor(src),
        wfs.state.temp,
        wfs.state.elongation_kernel,
    )
    return wfs
end

function apply_lgs_elongation!(::LGSProfileNaProfile, wfs::ShackHartmann, tel::Telescope, src::LGSSource, idx::Int)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution!(
        wfs.state.intensity,
        wfs.state.lgs_kernel_fft,
        wfs.state.fft_buffer,
        wfs.state.fft_plan,
        wfs.state.ifft_plan,
        idx,
    )
    return wfs
end

function ensure_lgs_kernels!(wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    na_profile = src.params.na_profile
    if na_profile === nothing
        return wfs
    end
    pad = size(wfs.state.intensity, 1)
    n_sub = wfs.params.n_subap
    tag = objectid(na_profile) ⊻ hash(src.params.laser_coordinates) ⊻ hash(src.params.fwhm_spot_up) ⊻ hash(pad)
    if size(wfs.state.lgs_kernel_fft, 1) == pad &&
        size(wfs.state.lgs_kernel_fft, 3) == n_sub * n_sub &&
        wfs.state.lgs_kernel_tag == tag
        return wfs
    end
    wfs.state.lgs_kernel_fft = lgs_spot_kernels_fft(tel, wfs, src, pad)
    wfs.state.lgs_kernel_tag = tag
    return wfs
end

function apply_lgs_convolution!(intensity::AbstractMatrix{T}, kernels_fft::AbstractArray{Complex{T},3},
    fft_buffer::AbstractMatrix{Complex{T}}, fft_plan, ifft_plan, idx::Int) where {T<:AbstractFloat}
    if size(kernels_fft, 3) < idx
        return intensity
    end
    kernel = @view kernels_fft[:, :, idx]
    apply_lgs_convolution!(intensity, kernel, fft_buffer, fft_plan, ifft_plan)
    return intensity
end

function lgs_spot_kernels_fft(tel::Telescope, wfs::ShackHartmann, src::LGSSource, pad::Int)
    T = eltype(wfs.state.intensity)
    n_sub = wfs.params.n_subap
    na_profile = src.params.na_profile
    altitudes = na_profile[1, :]
    weights = na_profile[2, :]
    if length(altitudes) == 0
        return similar(wfs.state.fft_buffer, Complex{T}, 0, 0, 0)
    end

    pixel_scale = lgs_pixel_scale(
        tel.params.diameter / wfs.params.n_subap,
        wfs.state.effective_padding,
        wavelength(src),
    )
    fwhm_px = src.params.fwhm_spot_up / pixel_scale
    sigma_px = fwhm_px / (2 * sqrt(2 * log(T(2))))
    center = (pad + 1) / 2

    x_subap = range(-tel.params.diameter / 2, tel.params.diameter / 2; length=n_sub)
    y_subap = range(-tel.params.diameter / 2, tel.params.diameter / 2; length=n_sub)
    kernels_fft = similar(wfs.state.fft_buffer, Complex{T}, pad, pad, n_sub * n_sub)

    x0 = src.params.laser_coordinates[2]
    y0 = -src.params.laser_coordinates[1]
    ref_idx = Int(cld(length(altitudes), 2))
    ref_vec = lgs_reference_vector(tel, x0, y0, altitudes[ref_idx])

    idx = 1
    kernel = Matrix{T}(undef, pad, pad)
    for iy in 1:n_sub, ix in 1:n_sub
        lgs_spot_kernel!(kernel, tel, src, altitudes, weights, pixel_scale, sigma_px, center,
            x_subap[ix], y_subap[iy], ref_vec)
        peak = maximum(kernel)
        if peak > 0
            cutoff = wfs.params.threshold_convolution * peak
            @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
                if kernel[i, j] < cutoff
                    kernel[i, j] = zero(T)
                end
            end
            total = sum(kernel)
            if total > 0
                kernel ./= total
            end
        end
        fft_buffer = wfs.state.fft_buffer
        @. fft_buffer = complex(kernel, zero(T))
        mul!(fft_buffer, wfs.state.fft_plan, fft_buffer)
        @views kernels_fft[:, :, idx] .= fft_buffer
        idx += 1
    end

    return kernels_fft
end
