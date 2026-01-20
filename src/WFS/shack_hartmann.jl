using FFTW
using Statistics

struct ShackHartmannParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    diffraction_padding::Int
end

mutable struct ShackHartmannState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    P,
    Pi,
    K<:AbstractVector{T},
    KF<:AbstractArray{Complex{T},3}}
    valid_mask::A
    slopes::V
    field::C
    fft_buffer::C
    intensity::R
    temp::R
    fft_plan::P
    ifft_plan::Pi
    elongation_kernel::K
    lgs_kernel_fft::KF
    lgs_kernel_tag::UInt
end

struct ShackHartmann{M<:SensingMode,P<:ShackHartmannParams,S<:ShackHartmannState} <: AbstractWFS
    params::P
    state::S
end

function ShackHartmann(tel::Telescope; n_subap::Int, threshold::Real=0.1,
    diffraction_padding::Int=2, mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)
    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    params = ShackHartmannParams{T}(n_subap, T(threshold), diffraction_padding)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    sub = div(tel.params.resolution, n_subap)
    pad = max(sub, sub * diffraction_padding)
    field = backend{Complex{T}}(undef, pad, pad)
    fft_buffer = similar(field)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    fft_plan = FFTW.plan_fft!(fft_buffer)
    ifft_plan = FFTW.plan_ifft!(fft_buffer)
    elongation_kernel = backend{T}(undef, 1)
    lgs_kernel_fft = backend{Complex{T}}(undef, 0, 0, 0)
    state = ShackHartmannState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(field),
        typeof(intensity),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
        typeof(lgs_kernel_fft),
    }(
        valid_mask,
        slopes,
        field,
        fft_buffer,
        intensity,
        temp,
        fft_plan,
        ifft_plan,
        elongation_kernel,
        lgs_kernel_fft,
        UInt(0),
    )
    wfs = ShackHartmann{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    return wfs
end

sensing_mode(::ShackHartmann{M}) where {M} = M()

function update_valid_mask!(wfs::ShackHartmann, tel::Telescope)
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

function measure!(mode::Geometric, wfs::ShackHartmann, tel::Telescope)
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
            wfs.state.slopes[idx] = sx / max(count_x, 1)
            wfs.state.slopes[idx + n_sub * n_sub] = sy / max(count_y, 1)
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end

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

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    Base.require_one_based_indexing(tel.state.opd)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    center = (pad + 1) / 2
    idx = 1

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1
        ys = (j - 1) * sub + 1
        xe = min(i * sub, n)
        ye = min(j * sub, n)
        if wfs.state.valid_mask[i, j]
            total, sx, sy = centroid_sums!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub, pad, idx)
            if total <= 0
                wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
                wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
            else
                wfs.state.slopes[idx] = sx / (total * center)
                wfs.state.slopes[idx + n_sub * n_sub] = sy / (total * center)
            end
        else
            wfs.state.slopes[idx] = zero(eltype(wfs.state.slopes))
            wfs.state.slopes[idx + n_sub * n_sub] = zero(eltype(wfs.state.slopes))
        end
        idx += 1
    end
    return wfs.state.slopes
end

function measure!(::Diffractive, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    Base.require_one_based_indexing(tel.state.opd)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    wavelength(ast)
    n = tel.params.resolution
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    pad = size(wfs.state.field, 1)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    center = (pad + 1) / 2
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
                wfs.state.slopes[idx] = sx_sum / (total_sum * center)
                wfs.state.slopes[idx + n_sub * n_sub] = sy_sum / (total_sum * center)
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
    phase_scale = (2 * pi) / wavelength(src)
    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+sub, oy+1:oy+sub] =
        tel.state.pupil[xs:xe, ys:ye] * cis(phase_scale * tel.state.opd[xs:xe, ys:ye])
    copyto!(wfs.state.fft_buffer, wfs.state.field)
    mul!(wfs.state.fft_buffer, wfs.state.fft_plan, wfs.state.fft_buffer)
    @. wfs.state.intensity = abs2(wfs.state.fft_buffer)
    return wfs.state.intensity
end

@inline function centroid_from_intensity(intensity::AbstractMatrix{T}, pad::Int) where {T<:AbstractFloat}
    total = sum(intensity)
    if total <= 0
        return zero(T), zero(T), zero(T)
    end
    sx = zero(T)
    sy = zero(T)
    center = (T(pad) + one(T)) / 2
    @inbounds for x in 1:pad, y in 1:pad
        val = intensity[x, y]
        sx += (x - center) * val
        sy += (y - center) * val
    end
    return total, sx, sy
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, ::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    return centroid_from_intensity(wfs.state.intensity, pad)
end

function centroid_sums!(wfs::ShackHartmann, tel::Telescope, src::LGSSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int, pad::Int, idx::Int)
    compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
    apply_lgs_elongation!(lgs_profile(src), wfs, tel, src, idx)
    return centroid_from_intensity(wfs.state.intensity, pad)
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
    if !(wfs.state.intensity isa Array)
        throw(InvalidConfiguration("LGS Na-profile kernels currently require Array backend"))
    end
    altitudes = na_profile[1, :]
    weights = na_profile[2, :]
    if length(altitudes) == 0
        return similar(wfs.state.fft_buffer, Complex{T}, 0, 0, 0)
    end

    pixel_scale = lgs_pixel_scale(
        tel.params.diameter / wfs.params.n_subap,
        wfs.params.diffraction_padding,
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
    kernel = similar(wfs.state.intensity)
    for iy in 1:n_sub, ix in 1:n_sub
        lgs_spot_kernel!(kernel, tel, src, altitudes, weights, pixel_scale, sigma_px, center,
            x_subap[ix], y_subap[iy], ref_vec)
        fft_buffer = wfs.state.fft_buffer
        @inbounds for i in 1:pad, j in 1:pad
            fft_buffer[i, j] = complex(kernel[i, j], zero(T))
        end
        mul!(fft_buffer, wfs.state.fft_plan, fft_buffer)
        @views kernels_fft[:, :, idx] .= fft_buffer
        idx += 1
    end

    return kernels_fft
end
