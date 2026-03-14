import Base: filter!

@kernel function circular_filter_mask_kernel!(mask, threshold2, center, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        x = i - center
        y = j - center
        r2 = x^2 + y^2
        value = r2 <= threshold2
        @inbounds mask[i, j] = (value + im * value) / sqrt(2)
    end
end

abstract type SpatialFilterShape end
struct CircularFilter <: SpatialFilterShape end
struct SquareFilter <: SpatialFilterShape end
struct FoucaultFilter <: SpatialFilterShape end

struct SpatialFilterParams{T<:AbstractFloat}
    diameter::T
    zero_padding::Int
    resolution::Int
end

mutable struct SpatialFilterState{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    Pf,
    Pi}
    mask::C
    mask_shifted::C
    field::C
    fft_buffer::C
    filtered_field::C
    fft_plan::Pf
    ifft_plan::Pi
    phase::R
    amplitude::R
end

struct SpatialFilter{S<:SpatialFilterShape,P<:SpatialFilterParams,Sf<:SpatialFilterState} <: AbstractOpticalElement
    params::P
    state::Sf
end

function SpatialFilter(tel::Telescope; shape::SpatialFilterShape=CircularFilter(),
    diameter::Real=tel.params.resolution / 2, zero_padding::Int=2,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    n_pad = tel.params.resolution * zero_padding
    params = SpatialFilterParams{T}(T(diameter), zero_padding, n_pad)
    mask = backend{Complex{T}}(undef, n_pad, n_pad)
    mask_shifted = similar(mask)
    field = similar(mask)
    fft_buffer = similar(mask)
    filtered_field = similar(mask)
    fft_plan = plan_fft_backend!(fft_buffer)
    ifft_plan = plan_ifft_backend!(filtered_field)
    phase = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    amplitude = similar(phase)
    state = SpatialFilterState{T, typeof(mask), typeof(phase), typeof(fft_plan), typeof(ifft_plan)}(
        mask,
        mask_shifted,
        field,
        fft_buffer,
        filtered_field,
        fft_plan,
        ifft_plan,
        phase,
        amplitude,
    )
    sf = SpatialFilter{typeof(shape), typeof(params), typeof(state)}(params, state)
    set_spatial_filter!(sf)
    return sf
end

function set_spatial_filter!(sf::SpatialFilter{CircularFilter})
    _set_spatial_filter!(execution_style(sf.state.mask), sf, CircularFilter())
    return sf
end

function _set_spatial_filter!(::ScalarCPUStyle, sf::SpatialFilter, ::CircularFilter)
    n = sf.params.resolution
    diameter_padded = sf.params.diameter * sf.params.zero_padding
    center = (n + 1) / 2

    @inbounds for i in 1:n, j in 1:n
        x = i - center
        y = j - center
        r2 = x^2 + y^2
        value = r2 <= diameter_padded^2
        sf.state.mask[i, j] = (value + im * value) / sqrt(2)
    end

    finalize_spatial_filter_mask!(sf)
    return sf
end

function _set_spatial_filter!(style::AcceleratorStyle, sf::SpatialFilter, ::CircularFilter)
    n = sf.params.resolution
    diameter_padded = sf.params.diameter * sf.params.zero_padding
    center = (n + 1) / 2
    launch_kernel!(style, circular_filter_mask_kernel!, sf.state.mask, diameter_padded^2, center, n; ndrange=size(sf.state.mask))
    finalize_spatial_filter_mask!(sf)
    return sf
end

function set_spatial_filter!(sf::SpatialFilter{SquareFilter})
    n = sf.params.resolution
    diameter_padded = sf.params.diameter * sf.params.zero_padding
    center = (n + 1) / 2
    fill!(sf.state.mask, zero(eltype(sf.state.mask)))
    half = Int(round(diameter_padded / 2))
    cx = Int(round(center))
    xs = max(1, cx - half)
    xe = min(n, cx + half)
    @views @. sf.state.mask[xs:xe, xs:xe] = (1 + im) / sqrt(2)

    finalize_spatial_filter_mask!(sf)
    return sf
end

function set_spatial_filter!(sf::SpatialFilter{FoucaultFilter})
    n = sf.params.resolution
    fill!(sf.state.mask, zero(eltype(sf.state.mask)))
    @views @. sf.state.mask[1:floor(Int, n / 2), :] = (1 + im) / sqrt(2)
    finalize_spatial_filter_mask!(sf)
    return sf
end

function finalize_spatial_filter_mask!(sf::SpatialFilter)
    tmp = similar(sf.state.mask)
    fill!(tmp, zero(eltype(tmp)))
    @views tmp[1:end-1, 1:end-1] .= sf.state.mask[2:end, 2:end]
    copyto!(sf.state.mask, tmp)
    fftshift2d!(sf.state.mask_shifted, sf.state.mask)
    return sf
end

function ensure_spatial_filter_buffers!(sf::SpatialFilter, n::Int, n_pad::Int)
    resized = false
    if size(sf.state.field) != (n_pad, n_pad)
        sf.state.field = similar(sf.state.field, n_pad, n_pad)
        sf.state.fft_buffer = similar(sf.state.fft_buffer, n_pad, n_pad)
        sf.state.filtered_field = similar(sf.state.filtered_field, n_pad, n_pad)
        sf.state.mask = similar(sf.state.mask, n_pad, n_pad)
        sf.state.mask_shifted = similar(sf.state.mask_shifted, n_pad, n_pad)
        sf.state.fft_plan = plan_fft_backend!(sf.state.fft_buffer)
        sf.state.ifft_plan = plan_ifft_backend!(sf.state.filtered_field)
        resized = true
    end
    if size(sf.state.phase) != (n, n)
        sf.state.phase = similar(sf.state.phase, n, n)
        sf.state.amplitude = similar(sf.state.amplitude, n, n)
    end
    if resized
        set_spatial_filter!(sf)
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

    copyto!(sf.state.fft_buffer, sf.state.field)
    mul!(sf.state.fft_buffer, sf.state.fft_plan, sf.state.fft_buffer)
    @. sf.state.filtered_field = sf.state.fft_buffer * sf.state.mask_shifted
    mul!(sf.state.filtered_field, sf.state.ifft_plan, sf.state.filtered_field)
    @views begin
        region = sf.state.filtered_field[ox+1:ox+n, oy+1:oy+n]
        @. sf.state.phase = angle(region) * tel.state.pupil
        @. sf.state.amplitude = abs(region)
    end
    return sf.state.phase, sf.state.amplitude
end
