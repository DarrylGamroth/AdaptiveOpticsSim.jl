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
    CC<:AbstractArray{Complex{T},3},
    R<:AbstractMatrix{T},
    RC3<:AbstractArray{T,3},
    RT3<:AbstractArray{T,3},
    RB<:AbstractMatrix{T},
    RS<:AbstractMatrix{T},
    RC<:AbstractArray{T,3},
    RN<:AbstractArray{T,3},
    P,
    PS,
    Pi,
    PSi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3}}
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
    spot_cube::RC
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
    intensity_asterism_stack::RC3
    fft_asterism_plan::PS
    asterism_capacity::Int
    valid_mask_host::Matrix{Bool}
    reference_signal_host::Vector{T}
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
    fft_stack = backend{Complex{T}}(undef, pad, pad, n_subap * n_subap)
    intensity = backend{T}(undef, pad, pad)
    intensity_stack = backend{T}(undef, pad, pad, n_subap * n_subap)
    intensity_tmp_stack = similar(intensity_stack)
    temp = similar(intensity)
    bin_buffer = backend{T}(undef, sub, sub)
    spot = similar(bin_buffer)
    spot_cube = backend{T}(undef, n_subap * n_subap, sub, sub)
    detector_noise_cube = similar(spot_cube)
    fft_plan = plan_fft_backend!(fft_buffer)
    fft_stack_plan = plan_fft_backend!(fft_stack, (1, 2))
    ifft_plan = plan_ifft_backend!(fft_buffer)
    ifft_stack_plan = plan_ifft_backend!(fft_stack, (1, 2))
    fft_asterism_stack = similar(fft_stack)
    intensity_asterism_stack = similar(intensity_stack)
    fft_asterism_plan = plan_fft_backend!(fft_asterism_stack, (1, 2))
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0, 0)
    valid_mask_host = Matrix{Bool}(undef, n_subap, n_subap)
    reference_signal_host = Vector{T}(undef, 2 * n_subap * n_subap)
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
        typeof(spot_cube),
        typeof(detector_noise_cube),
        typeof(fft_plan),
        typeof(fft_stack_plan),
        typeof(ifft_plan),
        typeof(ifft_stack_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
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
        spot_cube,
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
        backend{T}(undef, 2 * n_subap, n_subap),
        fft_asterism_stack,
        intensity_asterism_stack,
        fft_asterism_plan,
        1,
        valid_mask_host,
        reference_signal_host,
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
    copyto!(wfs.state.valid_mask_host, Array(wfs.state.valid_mask))
    return wfs
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
        wfs.state.intensity_tmp_stack = similar(wfs.state.intensity_tmp_stack, eltype(wfs.state.intensity), pad, pad, n_spots)
        wfs.state.temp = similar(wfs.state.temp, pad, pad)
        wfs.state.fft_plan = plan_fft_backend!(wfs.state.fft_buffer)
        wfs.state.fft_stack_plan = plan_fft_backend!(wfs.state.fft_stack, (1, 2))
        wfs.state.ifft_plan = plan_ifft_backend!(wfs.state.fft_buffer)
        wfs.state.ifft_stack_plan = plan_ifft_backend!(wfs.state.fft_stack, (1, 2))
        total = n_spots * wfs.state.asterism_capacity
        wfs.state.fft_asterism_stack = similar(wfs.state.fft_asterism_stack, Complex{eltype(wfs.state.field)}, pad, pad, total)
        wfs.state.intensity_asterism_stack = similar(wfs.state.intensity_asterism_stack, eltype(wfs.state.intensity), pad, pad, total)
        wfs.state.fft_asterism_plan = plan_fft_backend!(wfs.state.fft_asterism_stack, (1, 2))
        wfs.state.phasor_ratio = eltype(wfs.state.slopes)(NaN)
        wfs.state.calibrated = false
    end
    return wfs
end

function ensure_sh_asterism_buffers!(wfs::ShackHartmann, n_sources::Int)
    n_sources > 0 || throw(InvalidConfiguration("asterism must contain at least one source"))
    if n_sources <= wfs.state.asterism_capacity
        return wfs
    end
    pad = size(wfs.state.fft_stack, 1)
    n_spots = wfs.params.n_subap * wfs.params.n_subap
    total = n_spots * n_sources
    wfs.state.fft_asterism_stack = similar(wfs.state.fft_asterism_stack, Complex{eltype(wfs.state.field)}, pad, pad, total)
    wfs.state.intensity_asterism_stack = similar(wfs.state.intensity_asterism_stack, eltype(wfs.state.intensity), pad, pad, total)
    wfs.state.fft_asterism_plan = plan_fft_backend!(wfs.state.fft_asterism_stack, (1, 2))
    wfs.state.asterism_capacity = n_sources
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
    if sh_stacked_asterism_compatible(ast)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.state.slopes), wfs, tel, ast)
        sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
        subtract_reference_and_scale!(wfs)
        return wfs.state.slopes
    end
    return measure_sh_asterism_diffractive!(execution_style(wfs.state.slopes), wfs, tel, ast, n, n_sub, sub, pad, ox, oy)
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
    if sh_stacked_asterism_compatible(ast)
        peak = sampled_spots_peak_asterism_stacked!(execution_style(wfs.state.slopes), wfs, tel, ast, det, rng)
        sh_signal_from_spots!(wfs, peak, wfs.params.threshold_cog)
        subtract_reference_and_scale!(wfs)
        return wfs.state.slopes
    end
    return measure_sh_asterism_diffractive!(execution_style(wfs.state.slopes), wfs, tel, ast, det, rng, n, n_sub, sub, pad, ox, oy)
end

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
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
        idx += 1
    end
    zero_invalid_sh_slopes!(ScalarCPUStyle(), wfs.state.slopes, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope,
    ast::Asterism, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    host_valid_mask = wfs.state.valid_mask_host
    host_ref = wfs.state.reference_signal_host
    host_slopes = Array(wfs.state.slopes)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if host_valid_mask[i, j]
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
                host_slopes[idx] = zero(eltype(wfs.state.slopes))
                host_slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                host_slopes[idx] = ((sy_sum / nsrc) - host_ref[idx]) / wfs.state.slopes_units
                host_slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - host_ref[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            host_slopes[idx] = zero(eltype(wfs.state.slopes))
            host_slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.slopes, host_slopes)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
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
        idx += 1
    end
    zero_invalid_sh_slopes!(ScalarCPUStyle(), wfs.state.slopes, wfs.state.valid_mask)
    return wfs.state.slopes
end

function measure_sh_asterism_diffractive!(::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG, n::Int, n_sub::Int, sub::Int, pad::Int, ox::Int, oy::Int)
    host_valid_mask = wfs.state.valid_mask_host
    host_ref = wfs.state.reference_signal_host
    host_slopes = Array(wfs.state.slopes)
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if host_valid_mask[i, j]
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
                host_slopes[idx] = zero(eltype(wfs.state.slopes))
                host_slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                nsrc = length(ast.sources)
                host_slopes[idx] = ((sy_sum / nsrc) - host_ref[idx]) / wfs.state.slopes_units
                host_slopes[idx + n_sub * n_sub] = ((sx_sum / nsrc) - host_ref[idx + n_sub * n_sub]) / wfs.state.slopes_units
            end
        else
            host_slopes[idx] = zero(eltype(wfs.state.slopes))
            host_slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    copyto!(wfs.state.slopes, host_slopes)
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
    execute_fft_plan!(wfs.state.fft_buffer, wfs.state.fft_plan)
    @. wfs.state.intensity = abs2(wfs.state.fft_buffer)
    return wfs.state.intensity
end

@kernel function sh_field_stack_kernel!(fft_stack, valid_mask, pupil, opd, phasor,
    amp_scale, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_sub * n_sub
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(phasor)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            if ox + 1 <= x <= ox + sub && oy + 1 <= y <= oy + sub
                px = (i - 1) * sub + (x - ox)
                py = (j - 1) * sub + (y - oy)
                if px <= n && py <= n
                    amp = amp_scale * pupil[px, py]
                    val = amp * cispi(opd_to_cycles * opd[px, py]) * phasor[x, y]
                end
            end
        end
        @inbounds fft_stack[x, y, idx] = val
    end
end

@kernel function complex_abs2_stack_kernel!(intensity_stack, fft_stack, pad::Int, n_spots::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_spots
        @inbounds intensity_stack[x, y, idx] = abs2(fft_stack[x, y, idx])
    end
end

@kernel function sh_field_asterism_stack_kernel!(fft_stack, valid_mask, pupil, opd, phasor,
    amp_scales, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int, n_spots::Int)
    x, y, idx_total = @index(Global, NTuple)
    if x <= pad && y <= pad && idx_total <= size(fft_stack, 3)
        src_idx = (idx_total - 1) ÷ n_spots + 1
        idx = idx_total - (src_idx - 1) * n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(phasor)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            if ox + 1 <= x <= ox + sub && oy + 1 <= y <= oy + sub
                px = (i - 1) * sub + (x - ox)
                py = (j - 1) * sub + (y - oy)
                if px <= n && py <= n
                    amp = (@inbounds amp_scales[src_idx]) * pupil[px, py]
                    val = amp * cispi(opd_to_cycles * opd[px, py]) * phasor[x, y]
                end
            end
        end
        @inbounds fft_stack[x, y, idx_total] = val
    end
end

@kernel function reduce_asterism_intensity_stack_kernel!(intensity_stack, source_stack, n_spots::Int, n_src::Int, pad::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_spots
        acc = zero(eltype(intensity_stack))
        @inbounds for src_idx in 1:n_src
            acc += source_stack[x, y, idx + (src_idx - 1) * n_spots]
        end
        @inbounds intensity_stack[x, y, idx] = acc
    end
end

@kernel function sh_sample_spot_stack_kernel!(spot_cube, intensity_stack, valid_mask,
    binning::Int, n_sub::Int, n_binned::Int, n_out::Int)
    idx, u, v = @index(Global, NTuple)
    n_spots = n_sub * n_sub
    if idx <= n_spots && u <= n_out && v <= n_out
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(spot_cube)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            ox = div(n_out - n_binned, 2)
            oy = div(n_out - n_binned, 2)
            bu = u - ox
            bv = v - oy
            if 1 <= bu <= n_binned && 1 <= bv <= n_binned
                if binning == 1
                    @inbounds val = intensity_stack[bu, bv, idx]
                else
                    acc = zero(T)
                    @inbounds for ii in 1:binning, jj in 1:binning
                        acc += intensity_stack[(bu - 1) * binning + ii, (bv - 1) * binning + jj, idx]
                    end
                    val = acc
                end
            end
        end
        @inbounds spot_cube[idx, u, v] = val
    end
end

function compute_intensity_stack!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    pad = size(wfs.state.fft_stack, 1)
    n_sub = wfs.params.n_subap
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.state.intensity)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time *
        (tel.params.diameter / tel.params.resolution)^2))
    launch_kernel_async!(style, sh_field_stack_kernel!, wfs.state.fft_stack, wfs.state.valid_mask,
        tel.state.pupil, tel.state.opd, wfs.state.phasor, amp_scale, opd_to_cycles, n_sub, sub, ox, oy,
        tel.params.resolution, pad; ndrange=size(wfs.state.fft_stack))
    synchronize_backend!(style)
    execute_fft_plan!(wfs.state.fft_stack, wfs.state.fft_stack_plan)
    launch_kernel!(style, complex_abs2_stack_kernel!, wfs.state.intensity_stack, wfs.state.fft_stack,
        pad, n_sub * n_sub; ndrange=size(wfs.state.intensity_stack))
    return wfs.state.intensity_stack
end

@inline sh_stacked_asterism_compatible(ast::Asterism) =
    !isempty(ast.sources) && all(src -> !(src isa LGSSource), ast.sources)

function compute_intensity_asterism_stack!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    n_src = length(ast.sources)
    ensure_sh_asterism_buffers!(wfs, n_src)
    pad = size(wfs.state.fft_stack, 1)
    n_sub = wfs.params.n_subap
    n_spots = n_sub * n_sub
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.state.intensity)
    amp_scales = similar(wfs.state.slopes, n_src)
    @inbounds for i in eachindex(ast.sources)
        src = ast.sources[i]
        amp_scales[i] = sqrt(T(photon_flux(src) * tel.params.sampling_time *
            (tel.params.diameter / tel.params.resolution)^2))
    end
    total = n_spots * n_src
    fft_view = @view wfs.state.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view wfs.state.intensity_asterism_stack[:, :, 1:total]
    opd_to_cycles = T(2) / wavelength(ast)
    launch_kernel_async!(style, sh_field_asterism_stack_kernel!, fft_view, wfs.state.valid_mask,
        tel.state.pupil, tel.state.opd, wfs.state.phasor, amp_scales, opd_to_cycles,
        n_sub, sub, ox, oy, tel.params.resolution, pad, n_spots; ndrange=size(fft_view))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, wfs.state.fft_asterism_plan)
    launch_kernel_async!(style, complex_abs2_stack_kernel!, intensity_view, fft_view, pad, total; ndrange=size(intensity_view))
    launch_kernel!(style, reduce_asterism_intensity_stack_kernel!, wfs.state.intensity_stack, intensity_view,
        n_spots, n_src, pad; ndrange=size(wfs.state.intensity_stack))
    return wfs.state.intensity_stack
end

function sample_spot_stack!(::ScalarCPUStyle, wfs::ShackHartmann)
    n_spots = size(wfs.state.spot_cube, 1)
    @inbounds for idx in 1:n_spots
        sample_spot!(wfs, @view(wfs.state.intensity_stack[:, :, idx]))
        copyto!(sh_spot_view(wfs, idx), wfs.state.spot)
    end
    return wfs.state.spot_cube
end

function sample_spot_stack!(style::AcceleratorStyle, wfs::ShackHartmann)
    pad = size(wfs.state.intensity_stack, 1)
    binning = wfs.state.binning_pixel_scale
    @assert pad % binning == 0
    n_binned = div(pad, binning)
    n_out = size(wfs.state.spot_cube, 2)
    launch_kernel!(style, sh_sample_spot_stack_kernel!, wfs.state.spot_cube, wfs.state.intensity_stack,
        wfs.state.valid_mask, binning, wfs.params.n_subap, n_binned, n_out; ndrange=size(wfs.state.spot_cube))
    return wfs.state.spot_cube
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    fill!(wfs.state.detector_noise_cube, zero(eltype(wfs.state.detector_noise_cube)))
    for src in ast.sources
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, src)
        wfs.state.detector_noise_cube .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.detector_noise_cube)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    fill!(wfs.state.detector_noise_cube, zero(eltype(wfs.state.detector_noise_cube)))
    for src in ast.sources
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, src, det, rng)
        wfs.state.detector_noise_cube .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.detector_noise_cube)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    compute_intensity_asterism_stack!(style, wfs, tel, ast)
    sample_spot_stack!(style, wfs)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_asterism_stack!(style, wfs, tel, ast)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_subap
    capture_stack!(det, wfs.state.spot_cube, wfs.state.detector_noise_cube; rng=rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return maximum(wfs.state.spot_cube)
end

@kernel function sh_spot_centroid_kernel!(slopes, spot_cube, valid_mask, cutoff, n_sub::Int, offset::Int, n1::Int, n2::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(slopes)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                if val < cutoff
                    spot_cube[idx, x, y] = zero(T)
                else
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
            end
            if total <= 0
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            else
                slopes[idx] = sy / total
                slopes[idx + offset] = sx / total
            end
        else
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end

@kernel function zero_invalid_spots_kernel!(spot_cube, valid_mask, n_sub::Int, n1::Int, n2::Int)
    idx, x, y = @index(Global, NTuple)
    n_spots = n_sub * n_sub
    if idx <= n_spots && x <= n1 && y <= n2
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        if !(@inbounds valid_mask[i, j])
            @inbounds spot_cube[idx, x, y] = zero(eltype(spot_cube))
        end
    end
end

@kernel function zero_invalid_sh_slopes_kernel!(slopes, valid_mask, n_sub::Int, offset::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        if !(@inbounds valid_mask[i, j])
            T = eltype(slopes)
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
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

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
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

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    compute_intensity_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
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

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, tel, src)
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
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

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(style, wfs, tel, src)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_subap
    capture_stack!(det, wfs.state.spot_cube, wfs.state.detector_noise_cube; rng=rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak!(wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak!(execution_style(wfs.state.valid_mask), wfs, tel, src, det, rng)
end

function sampled_spots_peak!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    peak = zero(eltype(wfs.state.slopes))
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        spot_view = sh_spot_view(wfs, idx)
        if wfs.state.valid_mask[i, j]
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

function sampled_spots_peak!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak_lgs!(lgs_profile(src), style, wfs, tel, src, det, rng)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    compute_intensity_stack!(style, wfs, tel, src)
    apply_elongation_stack!(wfs.state.intensity_stack, lgs_elongation_factor(src),
        wfs.state.intensity_tmp_stack, wfs.state.elongation_kernel)
    sample_spot_stack!(style, wfs)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNone, style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_stack!(style, wfs, tel, src)
    apply_elongation_stack!(wfs.state.intensity_stack, lgs_elongation_factor(src),
        wfs.state.intensity_tmp_stack, wfs.state.elongation_kernel)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_subap
    capture_stack!(det, wfs.state.spot_cube, wfs.state.detector_noise_cube; rng=rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    n_sub = wfs.params.n_subap
    compute_intensity_stack!(style, wfs, tel, src)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution_stack!(wfs.state.intensity_stack, wfs.state.lgs_kernel_fft,
        wfs.state.fft_stack, wfs.state.fft_stack_plan, wfs.state.ifft_stack_plan)
    sample_spot_stack!(style, wfs)
    return maximum(wfs.state.spot_cube)
end

function sampled_spots_peak_lgs!(::LGSProfileNaProfile, style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    det::AbstractDetector, rng::AbstractRNG)
    n_sub = wfs.params.n_subap
    compute_intensity_stack!(style, wfs, tel, src)
    ensure_lgs_kernels!(wfs, tel, src)
    apply_lgs_convolution_stack!(wfs.state.intensity_stack, wfs.state.lgs_kernel_fft,
        wfs.state.fft_stack, wfs.state.fft_stack_plan, wfs.state.ifft_stack_plan)
    sample_spot_stack!(style, wfs)
    capture_stack!(det, wfs.state.spot_cube, wfs.state.detector_noise_cube; rng=rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return maximum(wfs.state.spot_cube)
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
    style = execution_style(wfs.state.slopes)
    launch_kernel!(style, sh_spot_centroid_kernel!, wfs.state.slopes, wfs.state.spot_cube,
        wfs.state.valid_mask, cutoff, n_sub, offset, size(wfs.state.spot_cube, 2),
        size(wfs.state.spot_cube, 3); ndrange=offset)
    return wfs.state.slopes
end

@inline function zero_invalid_sh_slopes!(::ScalarCPUStyle, slopes::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if !valid_mask[i, j]
            slopes[idx] = zero(T)
            slopes[idx + offset] = zero(T)
        end
        idx += 1
    end
    return slopes
end

@inline function zero_invalid_sh_slopes!(style::AcceleratorStyle, slopes::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    launch_kernel!(style, zero_invalid_sh_slopes_kernel!, slopes, valid_mask, n_sub, offset; ndrange=offset)
    return slopes
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
    count = sum(valid_mask)
    return count == 0 ? one(T) : sum(signal) / (T(2) * T(count))
end

function subtract_reference_and_scale!(wfs::ShackHartmann)
    inv_units = inv(wfs.state.slopes_units)
    ref = vec(wfs.state.reference_signal_2d)
    @. wfs.state.slopes = (wfs.state.slopes - ref) * inv_units
    return wfs.state.slopes
end

function subtract_reference!(wfs::ShackHartmann)
    ref = vec(wfs.state.reference_signal_2d)
    @. wfs.state.slopes = wfs.state.slopes - ref
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
    copyto!(wfs.state.reference_signal_host, vec(Array(wfs.state.reference_signal_2d)))
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
        _copy_real_kernel_to_complex!(execution_style(fft_buffer), fft_buffer, kernel)
        execute_fft_plan!(fft_buffer, wfs.state.fft_plan)
        @views kernels_fft[:, :, idx] .= fft_buffer
        idx += 1
    end

    return kernels_fft
end

function _copy_real_kernel_to_complex!(::ScalarCPUStyle, dest::AbstractMatrix{Complex{T}}, kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    @. dest = Complex{T}(kernel, zero(T))
    return dest
end

function _copy_real_kernel_to_complex!(::AcceleratorStyle, dest::AbstractMatrix{Complex{T}}, kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    host_complex = Matrix{Complex{T}}(undef, size(kernel)...)
    @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
        host_complex[i, j] = Complex{T}(kernel[i, j], zero(T))
    end
    copyto!(dest, host_complex)
    return dest
end
