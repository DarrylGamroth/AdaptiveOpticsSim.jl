@kernel function fill_telescope_field_kernel!(out, pupil_reflectivity, opd, phase_shift, amp_scale,
    opd_to_cycles, ox::Int, oy::Int, n::Int, n_pad::Int, center_even_grid::Bool)
    i, j = @index(Global, NTuple)
    if i <= n_pad && j <= n_pad
        xi = i - ox
        yj = j - oy
        val = zero(eltype(out))
        if 1 <= xi <= n && 1 <= yj <= n
            @inbounds val = amp_scale * sqrt(pupil_reflectivity[xi, yj]) * cispi(opd_to_cycles * opd[xi, yj])
        end
        if center_even_grid
            val *= cis(phase_shift * (i + j - 2))
        end
        @inbounds out[i, j] = val
    end
end

@kernel function apply_phase_opd_kernel!(field, phase_or_opd, opd_to_cycles, ox::Int, oy::Int, n::Int,
    full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(phase_or_opd, 1) && j <= size(phase_or_opd, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= cispi(opd_to_cycles * phase_or_opd[i, j])
    end
end

@kernel function apply_phase_rad_kernel!(field, phase_or_opd, ox::Int, oy::Int, n::Int, full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(phase_or_opd, 1) && j <= size(phase_or_opd, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= cis(phase_or_opd[i, j])
    end
end

@kernel function apply_amplitude_kernel!(field, amplitude, ox::Int, oy::Int, n::Int, full_field::Bool)
    i, j = @index(Global, NTuple)
    if i <= size(amplitude, 1) && j <= size(amplitude, 2)
        fi = full_field ? i : i + ox
        fj = full_field ? j : j + oy
        @inbounds field[fi, fj] *= amplitude[i, j]
    end
end

@kernel function intensity_kernel!(out, field, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds out[i, j] = abs2(field[i, j])
    end
end

@kernel function accumulate_abs2_kernel!(out, field, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds out[i, j] += abs2(field[i, j])
    end
end

struct ElectricFieldParams{T<:AbstractFloat}
    resolution::Int
    padded_resolution::Int
    zero_padding::Int
    wavelength::T
    sampling_m::T
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
    sampling_m = T(tel.params.diameter / tel.params.resolution)
    params = ElectricFieldParams{T}(n, n_pad, zero_padding, T(wavelength(src)), sampling_m)
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

    _fill_telescope_field!(execution_style(out), out, tel, src, zero_padding, center_even_grid)
    return out
end

function fill_telescope_field_async!(out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource;
    zero_padding::Int=1,
    center_even_grid::Bool=true) where {T<:AbstractFloat}
    zero_padding >= 1 || throw(InvalidConfiguration("zero_padding must be >= 1"))
    n = tel.params.resolution
    n_pad = n * zero_padding
    size(out) == (n_pad, n_pad) ||
        throw(DimensionMismatchError("electric field size must match telescope resolution * zero_padding"))
    _fill_telescope_field_async!(execution_style(out), out, tel, src, zero_padding, center_even_grid)
    return out
end

function _fill_telescope_field!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    fill!(out, zero(eltype(out)))
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    @views @. out[ox+1:ox+n, oy+1:oy+n] = amp_scale * sqrt(tel.state.pupil_reflectivity) * cispi(opd_to_cycles * tel.state.opd)
    if center_even_grid && iseven(n_pad)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(ScalarCPUStyle(), out, phase_shift)
    end
    return out
end

function _fill_telescope_field!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    phase_shift = center_even_grid && iseven(n_pad) ? -T(pi) * (T(n_pad) + one(T)) / T(n_pad) : zero(T)
    launch_kernel!(style, fill_telescope_field_kernel!, out, tel.state.pupil_reflectivity, tel.state.opd, phase_shift,
        amp_scale, opd_to_cycles, ox, oy, n, n_pad, center_even_grid && iseven(n_pad); ndrange=size(out))
    return out
end

function _fill_telescope_field_async!(::ScalarCPUStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    return _fill_telescope_field!(ScalarCPUStyle(), out, tel, src, zero_padding, center_even_grid)
end

function _fill_telescope_field_async!(style::AcceleratorStyle, out::AbstractMatrix{Complex{T}}, tel::Telescope, src::AbstractSource,
    zero_padding::Int, center_even_grid::Bool) where {T<:AbstractFloat}
    n = tel.params.resolution
    n_pad = n * zero_padding
    opd_to_cycles = T(2) / T(wavelength(src))
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2))
    ox, oy = field_embedding_offsets(n, n_pad)
    phase_shift = center_even_grid && iseven(n_pad) ? -T(pi) * (T(n_pad) + one(T)) / T(n_pad) : zero(T)
    launch_kernel_async!(style, fill_telescope_field_kernel!, out, tel.state.pupil_reflectivity, tel.state.opd, phase_shift,
        amp_scale, opd_to_cycles, ox, oy, n, n_pad, center_even_grid && iseven(n_pad); ndrange=size(out))
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

function fill_from_telescope_async!(ef::ElectricField, tel::Telescope, src::AbstractSource)
    wavelength(src) == ef.params.wavelength ||
        throw(InvalidConfiguration("source wavelength must match ElectricField wavelength"))
    tel.params.resolution == ef.params.resolution ||
        throw(DimensionMismatchError("ElectricField resolution must match telescope resolution"))
    ensure_field_buffers!(ef)
    fill_telescope_field_async!(ef.state.field, tel, src; zero_padding=ef.params.zero_padding)
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
    _apply_phase!(execution_style(field.state.field), field, phase_or_opd, units)
    return field
end

function apply_phase_async!(field::ElectricField, phase_or_opd::AbstractMatrix; units::Symbol=:opd)
    _apply_phase_async!(execution_style(field.state.field), field, phase_or_opd, units)
    return field
end

function _phase_target_layout(field::ElectricField, input::AbstractMatrix)
    if size(input) == size(field.state.field)
        return true, 0, 0
    elseif size(input) == (field.params.resolution, field.params.resolution)
        ox, oy = field_embedding_offsets(field.params.resolution, field.params.padded_resolution)
        return false, ox, oy
    end
    throw(DimensionMismatchError("field map size must match the active pupil or full padded field"))
end

function _apply_phase!(::ScalarCPUStyle, field::ElectricField, phase_or_opd::AbstractMatrix, units::Symbol)
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

function _apply_phase!(style::AcceleratorStyle, field::ElectricField, phase_or_opd::AbstractMatrix, units::Symbol)
    full_field, ox, oy = _phase_target_layout(field, phase_or_opd)
    if units === :opd
        T = eltype(field.state.intensity)
        opd_to_cycles = T(2) / field.params.wavelength
        launch_kernel!(style, apply_phase_opd_kernel!, field.state.field, phase_or_opd, opd_to_cycles,
            ox, oy, field.params.resolution, full_field; ndrange=size(phase_or_opd))
        return field
    elseif units === :phase
        launch_kernel!(style, apply_phase_rad_kernel!, field.state.field, phase_or_opd,
            ox, oy, field.params.resolution, full_field; ndrange=size(phase_or_opd))
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function _apply_phase_async!(::ScalarCPUStyle, field::ElectricField, phase_or_opd::AbstractMatrix, units::Symbol)
    return _apply_phase!(ScalarCPUStyle(), field, phase_or_opd, units)
end

function _apply_phase_async!(style::AcceleratorStyle, field::ElectricField, phase_or_opd::AbstractMatrix, units::Symbol)
    full_field, ox, oy = _phase_target_layout(field, phase_or_opd)
    if units === :opd
        T = eltype(field.state.intensity)
        opd_to_cycles = T(2) / field.params.wavelength
        launch_kernel_async!(style, apply_phase_opd_kernel!, field.state.field, phase_or_opd, opd_to_cycles,
            ox, oy, field.params.resolution, full_field; ndrange=size(phase_or_opd))
        return field
    elseif units === :phase
        launch_kernel_async!(style, apply_phase_rad_kernel!, field.state.field, phase_or_opd,
            ox, oy, field.params.resolution, full_field; ndrange=size(phase_or_opd))
        return field
    end
    throw(InvalidConfiguration("units must be :opd or :phase"))
end

function apply_amplitude!(field::ElectricField, amplitude::AbstractMatrix)
    _apply_amplitude!(execution_style(field.state.field), field, amplitude)
    return field
end

function _apply_amplitude!(::ScalarCPUStyle, field::ElectricField, amplitude::AbstractMatrix)
    target = field_target_view(field, amplitude)
    @. target *= amplitude
    return field
end

function _apply_amplitude!(style::AcceleratorStyle, field::ElectricField, amplitude::AbstractMatrix)
    full_field, ox, oy = _phase_target_layout(field, amplitude)
    launch_kernel!(style, apply_amplitude_kernel!, field.state.field, amplitude,
        ox, oy, field.params.resolution, full_field; ndrange=size(amplitude))
    return field
end

function intensity!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.state.field) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    _intensity!(execution_style(out), out, field.state.field)
    return out
end

function _intensity!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out = abs2(field_values)
    return out
end

function _intensity!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel!(style, intensity_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function accumulate_intensity!(out::AbstractMatrix{T}, field::ElectricField) where {T<:AbstractFloat}
    size(out) == size(field.state.field) ||
        throw(DimensionMismatchError("intensity output must match ElectricField size"))
    _accumulate_intensity!(execution_style(out), out, field.state.field)
    return out
end

function _accumulate_intensity!(::ScalarCPUStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. out += abs2(field_values)
    return out
end

function _accumulate_intensity!(style::AcceleratorStyle, out::AbstractMatrix{T}, field_values::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    launch_kernel!(style, accumulate_abs2_kernel!, out, field_values, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function intensity!(field::ElectricField)
    intensity!(field.state.intensity, field)
    return field.state.intensity
end
