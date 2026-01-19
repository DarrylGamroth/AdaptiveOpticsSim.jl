using FFTW
import Base: filter!

struct SpatialFilterParams{T<:AbstractFloat}
    shape::Symbol
    diameter::T
    zero_padding::Int
    resolution::Int
end

mutable struct SpatialFilterState{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T}}
    mask::C
    mask_shifted::C
    field::C
    filtered_field::C
    phase::R
    amplitude::R
end

struct SpatialFilter{P<:SpatialFilterParams,S<:SpatialFilterState} <: AbstractOpticalElement
    params::P
    state::S
end

function SpatialFilter(tel::Telescope; shape::Symbol=:circular, diameter::Real=tel.params.resolution / 2,
    zero_padding::Int=2, T::Type{<:AbstractFloat}=Float64, backend=Array)
    n_pad = tel.params.resolution * zero_padding
    params = SpatialFilterParams{T}(shape, T(diameter), zero_padding, n_pad)
    mask = backend{Complex{T}}(undef, n_pad, n_pad)
    mask_shifted = similar(mask)
    field = similar(mask)
    filtered_field = similar(mask)
    phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    amplitude = similar(phase)
    state = SpatialFilterState{T, typeof(mask), typeof(phase)}(
        mask,
        mask_shifted,
        field,
        filtered_field,
        phase,
        amplitude,
    )
    sf = SpatialFilter(params, state)
    set_spatial_filter!(sf)
    return sf
end

function set_spatial_filter!(sf::SpatialFilter)
    n = sf.params.resolution
    diameter_padded = sf.params.diameter * sf.params.zero_padding
    shape = sf.params.shape
    center = (n + 1) / 2

    if shape === :circular
        @inbounds for i in 1:n, j in 1:n
            x = i - center
            y = j - center
            r2 = x^2 + y^2
            value = r2 <= diameter_padded^2
            sf.state.mask[i, j] = (value + im * value) / sqrt(2)
        end
    elseif shape === :square
        fill!(sf.state.mask, zero(eltype(sf.state.mask)))
        half = Int(round(diameter_padded / 2))
        cx = Int(round(center))
        xs = max(1, cx - half)
        xe = min(n, cx + half)
        @views @. sf.state.mask[xs:xe, xs:xe] = (1 + im) / sqrt(2)
    elseif shape === :foucault
        fill!(sf.state.mask, zero(eltype(sf.state.mask)))
        @inbounds for i in 1:floor(Int, n / 2), j in 1:n
            sf.state.mask[i, j] = (1 + im) / sqrt(2)
        end
    else
        throw(InvalidConfiguration("shape must be :circular, :square, or :foucault"))
    end

    @views sf.state.mask[1:end-1, 1:end-1] .= sf.state.mask[2:end, 2:end]
    FFTW.fftshift!(sf.state.mask_shifted, sf.state.mask)
    return sf
end

function ensure_spatial_filter_buffers!(sf::SpatialFilter, n::Int, n_pad::Int)
    if size(sf.state.field) != (n_pad, n_pad)
        sf.state.field = similar(sf.state.field, n_pad, n_pad)
        sf.state.filtered_field = similar(sf.state.filtered_field, n_pad, n_pad)
        sf.state.mask = similar(sf.state.mask, n_pad, n_pad)
        sf.state.mask_shifted = similar(sf.state.mask_shifted, n_pad, n_pad)
    end
    if size(sf.state.phase) != (n, n)
        sf.state.phase = similar(sf.state.phase, n, n)
        sf.state.amplitude = similar(sf.state.amplitude, n, n)
    end
    return sf
end

function filter!(sf::SpatialFilter, tel::Telescope, src::AbstractSource)
    n = tel.params.resolution
    n_pad = sf.params.resolution
    ensure_spatial_filter_buffers!(sf, n, n_pad)

    fill!(sf.state.field, zero(eltype(sf.state.field)))
    phase_scale = (2 * pi) / wavelength(src)
    ox = div(n_pad - n, 2)
    oy = div(n_pad - n, 2)
    @views @. sf.state.field[ox+1:ox+n, oy+1:oy+n] = tel.state.pupil * cis(phase_scale * tel.state.opd)

    sf.state.filtered_field .= ifft(fft(sf.state.field) .* sf.state.mask_shifted)
    @views begin
        region = sf.state.filtered_field[ox+1:ox+n, oy+1:oy+n]
        @. sf.state.phase = angle(region) * tel.state.pupil
        @. sf.state.amplitude = abs(region)
    end
    return sf.state.phase, sf.state.amplitude
end
