struct ElectricFieldParams{T<:AbstractFloat}
    resolution::Int
    padded_resolution::Int
    zero_padding::Int
    wavelength::T
end

mutable struct ElectricFieldState{T<:AbstractFloat,
    C<:AbstractMatrix{Complex{T}},
    R<:AbstractMatrix{T},
    P}
    field::C
    fft_buffer::C
    intensity::R
    fft_plan::P
end

struct ElectricField{P<:ElectricFieldParams,S<:ElectricFieldState} <: AbstractOpticalElement
    params::P
    state::S
end

@inline function field_embedding_offsets(resolution::Int, padded_resolution::Int)
    return div(padded_resolution - resolution, 2), div(padded_resolution - resolution, 2)
end

@inline function field_active_axes(params::ElectricFieldParams)
    ox, oy = field_embedding_offsets(params.resolution, params.padded_resolution)
    return (ox + 1:ox + params.resolution, oy + 1:oy + params.resolution)
end

function ElectricField(tel::Telescope, src::AbstractSource;
    zero_padding::Int=1,
    T::Type{<:AbstractFloat}=eltype(tel.state.opd),
    backend=nothing)
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n = tel.params.resolution
    n_pad = n * zero_padding
    field = if backend === nothing
        similar(tel.state.opd, Complex{T}, n_pad, n_pad)
    else
        backend{Complex{T}}(undef, n_pad, n_pad)
    end
    fft_buffer = similar(field)
    intensity = if backend === nothing
        similar(tel.state.opd, T, n_pad, n_pad)
    else
        backend{T}(undef, n_pad, n_pad)
    end
    fft_plan = plan_fft_backend!(fft_buffer)
    params = ElectricFieldParams{T}(n, n_pad, zero_padding, T(wavelength(src)))
    state = ElectricFieldState{T, typeof(field), typeof(intensity), typeof(fft_plan)}(
        field,
        fft_buffer,
        intensity,
        fft_plan,
    )
    ef = ElectricField(params, state)
    fill_from_telescope!(ef, tel, src)
    return ef
end

function ensure_field_buffers!(ef::ElectricField)
    n_pad = ef.params.padded_resolution
    if size(ef.state.field) != (n_pad, n_pad)
        ef.state.field = similar(ef.state.field, n_pad, n_pad)
        ef.state.fft_buffer = similar(ef.state.fft_buffer, n_pad, n_pad)
        ef.state.intensity = similar(ef.state.intensity, n_pad, n_pad)
        ef.state.fft_plan = plan_fft_backend!(ef.state.fft_buffer)
    end
    return ef
end

function fill_telescope_field!(out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource;
    zero_padding::Int=1,
    center_even_grid::Bool=true) where {T<:AbstractFloat}
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n = tel.params.resolution
    n_pad = n * zero_padding
    size(out) == (n_pad, n_pad) ||
        throw(DimensionMismatchError("electric field size must match telescope resolution * zero_padding"))

    fill!(out, zero(eltype(out)))
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    @views @. out[ox+1:ox+n, oy+1:oy+n] = amp_scale * sqrt(tel.state.pupil_reflectivity) * cispi(opd_to_cycles * tel.state.opd)
    if center_even_grid && iseven(n_pad)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(execution_style(out), out, phase_shift)
    end
    return out
end

function fill_from_telescope!(ef::ElectricField, tel::Telescope, src::AbstractSource)
    wavelength(src) == ef.params.wavelength ||
        throw(InvalidConfiguration("source wavelength must match ElectricField wavelength"))
    tel.params.resolution == ef.params.resolution ||
        throw(DimensionMismatchError("ElectricField resolution must match telescope resolution"))
    ensure_field_buffers!(ef)
    fill_telescope_field!(ef.state.field, tel, src; zero_padding=ef.params.zero_padding)
    return ef
end

function field_target_view(field::ElectricField, input::AbstractMatrix)
    if size(input) == size(field.state.field)
        return field.state.field
    end
    active_axes = field_active_axes(field.params)
    if size(input) == (field.params.resolution, field.params.resolution)
        return @view field.state.field[active_axes...]
    end
    throw(DimensionMismatchError("field map size must match the active pupil or full padded field"))
end

function apply_phase!(field::ElectricField, phase_or_opd::AbstractMatrix; units::Symbol=:opd)
    target = field_target_view(field, phase_or_opd)
    if units === :opd
        T = eltype(field.state.intensity)
        opd_to_cycles = (T(2) / field.params.wavelength)
        @. target *= cispi(opd_to_cycles * phase_or_opd)
        return field
    elseif units === :phase
        @. target *= cis(phase_or_opd)
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function apply_amplitude!(field::ElectricField, amplitude::AbstractMatrix)
    target = field_target_view(field, amplitude)
    @. target *= amplitude
    return field
end

function intensity!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.state.field) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    @. out = abs2(field.state.field)
    return out
end

function intensity!(field::ElectricField)
    intensity!(field.state.intensity, field)
    return field.state.intensity
end
