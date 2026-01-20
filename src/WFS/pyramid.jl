using FFTW
using Statistics

struct PyramidParams{T<:AbstractFloat}
    n_subap::Int
    threshold::T
    modulation::T
    diffraction_padding::Int
end

mutable struct PyramidState{T<:AbstractFloat,
    A<:AbstractMatrix{Bool},
    V<:AbstractVector{T},
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi,
    K<:AbstractVector{T}}
    valid_mask::A
    slopes::V
    field::C
    focal_field::C
    pupil_field::C
    pyramid_mask::C
    intensity::R
    temp::R
    fft_plan::Pf
    ifft_plan::Pi
    elongation_kernel::K
end

struct PyramidWFS{M<:SensingMode,P<:PyramidParams,S<:PyramidState} <: AbstractWFS
    params::P
    state::S
end

function PyramidWFS(tel::Telescope; n_subap::Int, threshold::Real=0.1, modulation::Real=2.0,
    diffraction_padding::Int=2, mode::SensingMode=Geometric(), T::Type{<:AbstractFloat}=Float64, backend=Array)

    if tel.params.resolution % n_subap != 0
        throw(InvalidConfiguration("telescope resolution must be divisible by n_subap"))
    end
    params = PyramidParams{T}(n_subap, T(threshold), T(modulation), diffraction_padding)
    valid_mask = backend{Bool}(undef, n_subap, n_subap)
    slopes = backend{T}(undef, 2 * n_subap * n_subap)
    fill!(slopes, zero(T))
    pad = tel.params.resolution * diffraction_padding
    field = backend{Complex{T}}(undef, pad, pad)
    focal_field = similar(field)
    pupil_field = similar(field)
    pyramid_mask = similar(field)
    intensity = backend{T}(undef, pad, pad)
    temp = similar(intensity)
    fft_plan = FFTW.plan_fft!(focal_field)
    ifft_plan = FFTW.plan_ifft!(pupil_field)
    elongation_kernel = backend{T}(undef, 1)
    state = PyramidState{
        T,
        typeof(valid_mask),
        typeof(slopes),
        typeof(field),
        typeof(intensity),
        typeof(fft_plan),
        typeof(ifft_plan),
        typeof(elongation_kernel),
    }(
        valid_mask,
        slopes,
        field,
        focal_field,
        pupil_field,
        pyramid_mask,
        intensity,
        temp,
        fft_plan,
        ifft_plan,
        elongation_kernel,
    )
    wfs = PyramidWFS{typeof(mode), typeof(params), typeof(state)}(params, state)
    update_valid_mask!(wfs, tel)
    build_pyramid_mask!(wfs.state.pyramid_mask, wfs.params.modulation)
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

function measure!(::Diffractive, wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    Base.require_one_based_indexing(tel.state.opd)
    n = tel.params.resolution
    pad = size(wfs.state.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    phase_scale = (2 * pi) / wavelength(src)
    center = (pad + 1) / 2
    n_sub = wfs.params.n_subap
    sub = div(n, n_sub)
    idx = 1

    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+n, oy+1:oy+n] = tel.state.pupil * cis(phase_scale * tel.state.opd)
    copyto!(wfs.state.focal_field, wfs.state.field)
    mul!(wfs.state.focal_field, wfs.state.fft_plan, wfs.state.focal_field)
    @. wfs.state.focal_field = wfs.state.focal_field * wfs.state.pyramid_mask
    copyto!(wfs.state.pupil_field, wfs.state.focal_field)
    mul!(wfs.state.pupil_field, wfs.state.ifft_plan, wfs.state.pupil_field)
    @. wfs.state.intensity = abs2(wfs.state.pupil_field)

    if src isa LGSSource
        wfs.state.elongation_kernel = apply_elongation!(
            wfs.state.intensity,
            lgs_elongation_factor(src),
            wfs.state.temp,
            wfs.state.elongation_kernel,
        )
    end

    @inbounds for i in 1:n_sub, j in 1:n_sub
        xs = (i - 1) * sub + 1 + ox
        ys = (j - 1) * sub + 1 + oy
        xe = xs + sub - 1
        ye = ys + sub - 1
        if wfs.state.valid_mask[i, j]
            total = zero(eltype(wfs.state.slopes))
            left = zero(eltype(wfs.state.slopes))
            right = zero(eltype(wfs.state.slopes))
            bottom = zero(eltype(wfs.state.slopes))
            top = zero(eltype(wfs.state.slopes))
            midx = div(xs + xe, 2)
            midy = div(ys + ye, 2)
            for x in xs:xe, y in ys:ye
                val = wfs.state.intensity[x, y]
                total += val
                if y <= midy
                    left += val
                else
                    right += val
                end
                if x <= midx
                    bottom += val
                else
                    top += val
                end
            end
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
