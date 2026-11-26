@inline wrap_upper_index(i::Int, n::Int) = ifelse(i > n, i - n, i)

# Native pupil and screen arrays use [x,y], matching their plane orientation.
# Sampling is center + cone_scale * pupil_position - wind_offset + source_shift.
@kernel function moving_layer_extract_kernel!(out, screen, start_x, start_y, footprint_scale, output_scale, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = start_y + footprint_scale * T(j - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = start_x + footprint_scale * T(i - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[ix0, iy0]
            v01 = screen[ix1, iy0]
            v10 = screen[ix0, iy1]
            v11 = screen[ix1, iy1]
            out[i, j] = output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

@kernel function moving_layer_accumulate_kernel!(out, screen, start_x, start_y, footprint_scale, output_scale, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = start_y + footprint_scale * T(j - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = start_x + footprint_scale * T(i - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[ix0, iy0]
            v01 = screen[ix1, iy0]
            v10 = screen[ix0, iy1]
            v11 = screen[ix1, iy1]
            out[i, j] += output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

@kernel function moving_layer_replay_kernel!(
    out,
    screen,
    start_x,
    start_y,
    screen_index::Int,
    footprint_scale,
    output_scale,
    accumulate::Bool,
    n::Int,
    m::Int,
)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = @inbounds(start_y[screen_index]) + footprint_scale * T(j - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = @inbounds(start_x[screen_index]) + footprint_scale * T(i - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[ix0, iy0]
            v01 = screen[ix1, iy0]
            v10 = screen[ix0, iy1]
            v11 = screen[ix1, iy1]
            value = output_scale * (
                wy0 * (wx0 * v00 + fx * v01) +
                fy * (wx0 * v10 + fx * v11)
            )
            if accumulate
                out[i, j] += value
            else
                out[i, j] = value
            end
        end
    end
end

@inline function _normalize_replay_start(start::T, period::T) where {
    T<:AbstractFloat,
}
    return start - floor((start - one(T)) / period) * period
end

@kernel function advance_moving_screen_replay_kernel!(
    offset_x,
    offset_y,
    start_x,
    start_y,
    offset_increment_x,
    offset_increment_y,
    base_start_x,
    base_start_y,
    screen_period,
    sampling_ratio,
    screen_count::Int,
)
    i = @index(Global, Linear)
    if i <= screen_count
        @inbounds begin
            next_offset_x = offset_x[i] + offset_increment_x[i]
            next_offset_y = offset_y[i] + offset_increment_y[i]
            offset_x[i] = next_offset_x
            offset_y[i] = next_offset_y
            start_x[i] = _normalize_replay_start(
                base_start_x[i] - next_offset_x * sampling_ratio[i],
                screen_period[i],
            )
            start_y[i] = _normalize_replay_start(
                base_start_y[i] - next_offset_y * sampling_ratio[i],
                screen_period[i],
            )
        end
    end
end

@kernel function apply_atmosphere_pupil_kernel!(opd, pupil, count::Int)
    i = @index(Global, Linear)
    if i <= count
        @inbounds opd[i] *= pupil[i]
    end
end

struct MultiLayerParams{T<:AbstractFloat,
    V1<:AbstractVector{T},
    V2<:AbstractVector{T},
    V3<:AbstractVector{T},
    V4<:AbstractVector{T},
    V5<:AbstractVector{T},
    V6<:AbstractVector{T},
    I<:Tuple}
    cn2_fractions::V1
    wind_speed::V2
    wind_direction_deg::V3
    altitude::V4
    wind_velocity_x::V5
    wind_velocity_y::V6
    layer_ids::I
    r0::T
    reference_wavelength_m::T
    L0::T
end

struct MovingLayerParams{T<:AbstractFloat}
    cn2_amplitude_scale::T
    wind_velocity_x::T
    wind_velocity_y::T
    altitude::T
    sampling_m::T
end

mutable struct MovingLayerState{T<:AbstractFloat}
    offset_x::T
    offset_y::T
    initialized::Bool
    low_frequency_initialized::Bool
end

struct MovingAtmosphereLayer{
    P<:MovingLayerParams,
    S<:MovingLayerState,
    A<:KolmogorovAtmosphere,
    TT<:Telescope,
    F<:Union{Nothing,KolmogorovAtmosphere},
} <: AbstractAtmosphereLayer
    params::P
    generator::A
    generator_telescope::TT
    low_frequency::F
    state::S
end

mutable struct MultiLayerState{T<:AbstractFloat}
    timeline::AtmosphereTimelineState{T}
end

struct MultiLayerAtmosphere{
    P<:MultiLayerParams,
    S<:MultiLayerState,
    L<:AbstractVector{<:MovingAtmosphereLayer},
    I<:AtmosphereIdentity,
    B<:AbstractArrayBackend,
} <: AbstractTimedAtmosphere
    params::P
    layers::L
    state::S
    identity::I
end

"""Run-immutable moving-screen coefficients used by device-graph replay."""
struct _MovingScreenReplayPlan{T<:AbstractFloat,A<:AbstractVector{T}}
    offset_increment_x::A
    offset_increment_y::A
    base_start_x::A
    base_start_y::A
    screen_period::A
    sampling_ratio::A
end

"""Persistent moving-screen offsets in native pixels for one device-graph owner."""
mutable struct _MovingScreenReplayState{T<:AbstractFloat,A<:AbstractVector{T}}
    offset_x::A
    offset_y::A
end

"""Replaceable device storage derived from moving-screen replay state."""
struct _MovingScreenReplayWorkspace{T<:AbstractFloat,A<:AbstractVector{T}}
    start_x::A
    start_y::A
end

"""Exact prepared owner for one renderer's device-resident screen motion."""
struct _PreparedMovingScreenReplay{P,S,W}
    plan::P
    state::S
    workspace::W
end

@inline backend(::MultiLayerAtmosphere{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()
@inline atmosphere_numeric_type(atm::MultiLayerAtmosphere) =
    typeof(atm.params.r0)

@inline moving_layer_screen_resolution(n::Int) = 3 * n

function _prepare_low_frequency_screen(generator::KolmogorovAtmosphere)
    # Replace the central native Fourier cells with a 27× finer frequency grid.
    # Retain at least seven samples per wavelength of the highest refined mode
    # to limit attenuation from bilinear pupil extraction.
    period_ratio = 27
    native_n = size(generator.state.phase_rad, 1)
    native_max_mode = min(6, (native_n - 1) ÷ 2)
    samples_per_native_period = min(48, 8 * (native_max_mode + 1))
    n = period_ratio * samples_per_native_period
    params = generator.params
    native_period_m = native_n * params.sampling_m
    sampling_m = native_period_m / samples_per_native_period
    low_frequency = _prepare_kolmogorov_atmosphere(n, sampling_m, params.r0,
        params.reference_wavelength_m, params.L0, backend(generator))

    # FFT order is 0, 1, …, -2, -1. Zero the central square on the native grid
    # and keep the coarse indices inside its outer cell boundaries, located at
    # ±(native_max_mode + ½) in native frequency bins. DC stays zero on both.
    ensure_psd!(generator, params.sampling_m)
    native_psd = generator.state.psd
    native_band = (1:(native_max_mode + 1), (native_n - native_max_mode + 1):native_n)
    for rows in native_band, columns in native_band
        fill!(view(native_psd, rows, columns), zero(eltype(native_psd)))
    end
    ensure_psd!(low_frequency, sampling_m)
    max_mode = period_ratio * native_max_mode + period_ratio ÷ 2
    excluded = (max_mode + 2):(n - max_mode)
    psd = low_frequency.state.psd
    fill!(view(psd, excluded, :), zero(eltype(psd)))
    fill!(view(psd, :, excluded), zero(eltype(psd)))
    return low_frequency
end

@inline _moving_layer_screens(layer::MovingAtmosphereLayer) =
    _moving_layer_screens(layer.generator, layer.low_frequency)
@inline _moving_layer_screens(generator, ::Nothing) = (generator,)
@inline _moving_layer_screens(generator, low_frequency) = (generator, low_frequency)

@inline _moving_screen_sampling_ratio(layer::MovingAtmosphereLayer,
    generator::KolmogorovAtmosphere) =
    layer.generator.params.sampling_m / generator.params.sampling_m

function moving_layer_telescope(
    tel::Telescope;
    resolution::Int,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel),
)
    delta = tel.params.diameter / tel.params.resolution
    return Telescope(
        resolution=resolution,
        diameter=delta * resolution,
        central_obstruction=0.0,
        fov_arcsec=tel.params.fov_arcsec,
        T=T,
        backend=backend,
    )
end

function MovingAtmosphereLayer(
    tel::Telescope;
    r0::Real,
    reference_wavelength_m::Real,
    L0::Real,
    cn2_fraction::Real,
    wind_velocity_x::Real,
    wind_velocity_y::Real,
    altitude::Real,
    low_frequency_completion::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
)
    selector = require_same_backend(tel, backend)
    screen_resolution = moving_layer_screen_resolution(tel.params.resolution)
    screen_telescope = moving_layer_telescope(tel; resolution=screen_resolution, T=T, backend=selector)
    generator = KolmogorovAtmosphere(screen_telescope;
        r0,
        reference_wavelength_m,
        L0,
        T,
        backend=selector,
    )
    params = MovingLayerParams(T(sqrt(cn2_fraction)), T(wind_velocity_x),
        T(wind_velocity_y), T(altitude), T(tel.aperture.sampling_m[1]))
    low_frequency = low_frequency_completion ? _prepare_low_frequency_screen(generator) : nothing
    state = MovingLayerState(zero(T), zero(T), false, false)
    return MovingAtmosphereLayer(params, generator, screen_telescope, low_frequency, state)
end

function ensure_initialized!(layer::MovingAtmosphereLayer, rng::AbstractRNG)
    if !layer.state.initialized
        regenerate_phase_screen!(layer.generator, rng)
        layer.state.initialized = true
    end
    return layer
end

@inline wrap_index(i::Int, n::Int) = mod1(i, n)
@inline normalize_start_coordinate(start::T, m::Int) where {T<:AbstractFloat} = mod(start - one(T), T(m)) + one(T)

function extract_shifted_screen!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, output_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration("footprint_scale must be positive"))

    # This sampling/extraction helper is shared infrastructure for moving-screen
    # atmosphere models. The current finite backend uses periodic wraparound in
    # `_extract_shifted_screen!`; the planned infinite backend will reuse the
    # same extraction interface after updating its persistent buffer via
    # boundary injection instead of wraparound.
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    _extract_shifted_screen!(execution_style(out), out, screen, start_x, start_y, footprint_scale, output_scale, n, m)
    return out
end

function extract_shifted_screen_async!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, output_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration("footprint_scale must be positive"))
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    _extract_shifted_screen_async!(execution_style(out), out, screen, start_x, start_y, footprint_scale, output_scale, n, m)
    return out
end

function _extract_shifted_screen!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    # Finite moving-screen backend: move a periodic canvas under the pupil with
    # bilinear subpixel interpolation.
    @inbounds for j in 1:n
        y = start_y + footprint_scale * T(j - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        for i in 1:n
            x = start_x + footprint_scale * T(i - 1)
            x0 = floor(Int, x)
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_index(x0, m)
            ix1 = wrap_index(x0 + 1, m)

            v00 = screen[ix0, iy0]
            v01 = screen[ix1, iy0]
            v10 = screen[ix0, iy1]
            v11 = screen[ix1, iy1]
            out[i, j] = output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
    return out
end

function _extract_shifted_screen!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel!(style, moving_layer_extract_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, output_scale, n, m; ndrange=size(out))
    return out
end

function accumulate_shifted_screen!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, output_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration("footprint_scale must be positive"))
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    _accumulate_shifted_screen!(execution_style(out), out, screen, start_x, start_y, footprint_scale, output_scale, n, m)
    return out
end

@inline function _extract_shifted_screen_async!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    return _extract_shifted_screen!(ScalarCPUStyle(), out, screen, start_x, start_y, footprint_scale, output_scale, n, m)
end

function _extract_shifted_screen_async!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel_async!(style, moving_layer_extract_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, output_scale, n, m; ndrange=size(out))
    return out
end

function _accumulate_shifted_screen!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    @inbounds for j in 1:n
        y = start_y + footprint_scale * T(j - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        for i in 1:n
            x = start_x + footprint_scale * T(i - 1)
            x0 = floor(Int, x)
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_index(x0, m)
            ix1 = wrap_index(x0 + 1, m)

            v00 = screen[ix0, iy0]
            v01 = screen[ix1, iy0]
            v10 = screen[ix0, iy1]
            v11 = screen[ix1, iy1]
            out[i, j] += output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
    return out
end

function _accumulate_shifted_screen!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, output_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel!(style, moving_layer_accumulate_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, output_scale, n, m; ndrange=size(out))
    return out
end

@inline function _accumulate_shifted_screen_async!(::ScalarCPUStyle,
    out::AbstractMatrix{T}, screen::AbstractMatrix{T}, start_x::T,
    start_y::T, footprint_scale::T, output_scale::T, n::Int,
    m::Int) where {T<:AbstractFloat}
    return _accumulate_shifted_screen!(ScalarCPUStyle(), out, screen,
        start_x, start_y, footprint_scale, output_scale, n, m)
end

function _accumulate_shifted_screen_async!(style::AcceleratorStyle,
    out::AbstractMatrix{T}, screen::AbstractMatrix{T}, start_x::T,
    start_y::T, footprint_scale::T, output_scale::T, n::Int,
    m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel_async!(style, moving_layer_accumulate_kernel!, out, screen,
        start_x_wrapped, start_y_wrapped, footprint_scale, output_scale, n, m;
        ndrange=size(out))
    return out
end

function _accumulate_shifted_screen_async!(out::AbstractMatrix{T},
    screen::AbstractMatrix{T}, offset_x::T, offset_y::T, output_scale::T,
    footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError(
        "screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration(
        "footprint_scale must be positive"))
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    return _accumulate_shifted_screen_async!(execution_style(out), out,
        screen, start_x, start_y, footprint_scale, output_scale, n, m)
end

@inline _accumulate_low_frequency!(accumulate!, out, ::Nothing,
    offset_x, offset_y, footprint_scale, output_scale, native_sampling) = out

function _accumulate_low_frequency!(accumulate!, out::AbstractMatrix{T},
    screen::KolmogorovAtmosphere, offset_x::T, offset_y::T, footprint_scale::T,
    output_scale::T, native_sampling::T) where {T<:AbstractFloat}
    phase_rad = screen.state.phase_rad
    n = size(out, 1)
    m = size(phase_rad, 1)
    sampling_ratio = native_sampling / screen.params.sampling_m
    scale = footprint_scale * sampling_ratio
    center = T(m + 1) / 2 - scale * T(n - 1) / 2
    return accumulate!(execution_style(out), out, phase_rad,
        center - offset_x * sampling_ratio, center - offset_y * sampling_ratio,
        scale, output_scale, n, m)
end

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    offset_x = layer.state.offset_x - shift_x
    offset_y = layer.state.offset_y - shift_y
    extract_shifted_screen_async!(out, layer.generator.state.phase_rad,
        offset_x, offset_y, opd_scale, footprint_scale)
    return _accumulate_low_frequency!(_accumulate_shifted_screen_async!, out,
        layer.low_frequency, offset_x, offset_y, footprint_scale, opd_scale,
        layer.generator.params.sampling_m)
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer_accumulate!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    offset_x = layer.state.offset_x - shift_x
    offset_y = layer.state.offset_y - shift_y
    accumulate_shifted_screen!(out, layer.generator.state.phase_rad,
        offset_x, offset_y, opd_scale, footprint_scale)
    return _accumulate_low_frequency!(_accumulate_shifted_screen!, out,
        layer.low_frequency, offset_x, offset_y, footprint_scale, opd_scale,
        layer.generator.params.sampling_m)
end

function _render_layer_accumulate_async!(out::AbstractMatrix{T},
    layer::MovingAtmosphereLayer, shift_x::T, shift_y::T,
    footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    offset_x = layer.state.offset_x - shift_x
    offset_y = layer.state.offset_y - shift_y
    _accumulate_shifted_screen_async!(out, layer.generator.state.phase_rad,
        offset_x, offset_y, opd_scale, footprint_scale)
    return _accumulate_low_frequency!(_accumulate_shifted_screen_async!, out,
        layer.low_frequency, offset_x, offset_y, footprint_scale, opd_scale,
        layer.generator.params.sampling_m)
end

function evolve_layer!(layer::MovingAtmosphereLayer, duration::Real)
    T = typeof(layer.state.offset_x)
    dt = T(duration)
    delta = layer.params.sampling_m
    layer.state.offset_x += layer.params.wind_velocity_x * dt / delta
    layer.state.offset_y += layer.params.wind_velocity_y * dt / delta
    return layer
end

"""
    MultiLayerAtmosphere(telescope; low_frequency_completion=false, ...)

Prepare independent finite von Kármán phase screens advected by the specified
layer winds. With `low_frequency_completion=true`, each layer samples its
central square of Fourier cells at 27 times finer frequency spacing, excluding
the zero frequency. The two bands undergo the same physical advection and
source footprint geometry. Retained native Fourier coefficients continue to
use the original pupil sampling and seeded draws.
"""
function MultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    reference_wavelength_m::Real,
    L0::Real=25.0,
    fractional_cn2::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction_deg::AbstractVector,
    altitude::AbstractVector,
    layer_ids=nothing,
    low_frequency_completion::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))

    selector = require_same_backend(tel, backend)
    n_layers = length(fractional_cn2)
    n_layers > 0 || throw(InvalidConfiguration("fractional_cn2 cannot be empty"))
    if length(wind_speed) != n_layers || length(wind_direction_deg) != n_layers || length(altitude) != n_layers
        throw(InvalidConfiguration("layer parameter lengths must match fractional_cn2"))
    end
    all(>=(0), fractional_cn2) || throw(InvalidConfiguration("fractional_cn2 must be non-negative"))
    isapprox(sum(fractional_cn2), 1; atol=1e-6, rtol=1e-6) ||
        throw(InvalidConfiguration("fractional_cn2 must sum to 1"))
    prepared_layer_ids = _prepare_atmosphere_layer_ids(layer_ids, n_layers)
    r0_t = _converted_positive_finite(r0, T,
        "atmosphere Fried parameter r0")
    reference_wavelength_t = _converted_positive_finite(
        reference_wavelength_m,
        T,
        "atmosphere reference wavelength in metres",
    )
    L0_t = _converted_positive_finite(L0, T,
        "atmosphere outer scale L0")

    params = MultiLayerParams(
        T.(fractional_cn2),
        T.(wind_speed),
        T.(wind_direction_deg),
        T.(altitude),
        T[T(wind_speed[i]) * cosd(T(wind_direction_deg[i])) for i in 1:n_layers],
        T[T(wind_speed[i]) * sind(T(wind_direction_deg[i])) for i in 1:n_layers],
        prepared_layer_ids,
        r0_t,
        reference_wavelength_t,
        L0_t,
    )

    layers = [
        MovingAtmosphereLayer(
            tel;
            r0=params.r0,
            reference_wavelength_m=params.reference_wavelength_m,
            L0=params.L0,
            cn2_fraction=params.cn2_fractions[i],
            wind_velocity_x=params.wind_velocity_x[i],
            wind_velocity_y=params.wind_velocity_y[i],
            altitude=params.altitude[i],
            low_frequency_completion,
            T=T,
            backend=selector,
        ) for i in 1:n_layers
    ]

    state = MultiLayerState{T}(new_atmosphere_timeline(T))
    identity = AtmosphereIdentity()
    return MultiLayerAtmosphere{typeof(params), typeof(state), typeof(layers),
        typeof(identity), typeof(selector)}(params, layers, state, identity)
end

function initialize_atmosphere!(atm::MultiLayerAtmosphere, rng::AbstractRNG)
    @inbounds for layer in atm.layers
        ensure_initialized!(layer, rng)
    end
    # Preserve every native layer's RNG draws when completion is enabled.
    @inbounds for layer in atm.layers
        _ensure_low_frequency_initialized!(layer.state, layer.low_frequency, rng)
    end
    return atm
end

@inline _ensure_low_frequency_initialized!(::MovingLayerState, ::Nothing, ::AbstractRNG) = nothing

function _ensure_low_frequency_initialized!(state::MovingLayerState,
    screen::KolmogorovAtmosphere, rng::AbstractRNG)
    if !state.low_frequency_initialized
        regenerate_phase_screen!(screen, rng)
        state.low_frequency_initialized = true
    end
    return nothing
end

function evolve_atmosphere!(atm::MultiLayerAtmosphere, duration::Real,
    ::AbstractRNG)
    @inbounds for layer in atm.layers
        evolve_layer!(layer, duration)
    end
    return atm
end

function _copy_moving_screen_replay_vector(
    target::AbstractComputeDevice,
    values::Vector{T},
) where {T<:AbstractFloat}
    storage = allocate_device_array(target, T, length(values))
    copyto!(storage, values)
    return storage
end

function _prepare_moving_screen_replay(
    atm::MultiLayerAtmosphere,
    renderer::AtmosphereDirectionRenderer,
    duration::Real,
)
    _validate_atmosphere_renderer_binding(renderer, atm)
    T = atmosphere_numeric_type(atm)
    dt = _explicit_atmosphere_duration(duration, T)
    dt > zero(T) || throw(AtmosphereTimeError(
        "moving-screen replay duration must be positive",
    ))
    layers = atm.layers
    screen_count = sum(layer -> length(_moving_layer_screens(layer)), layers)
    n = size(renderer.pupil, 1)
    target = compute_device(renderer.pupil)
    host_offset_increment_x = Vector{T}(undef, screen_count)
    host_offset_increment_y = Vector{T}(undef, screen_count)
    host_base_start_x = Vector{T}(undef, screen_count)
    host_base_start_y = Vector{T}(undef, screen_count)
    host_screen_period = Vector{T}(undef, screen_count)
    host_sampling_ratio = Vector{T}(undef, screen_count)

    screen_index = 0
    @inbounds for i in eachindex(layers)
        layer = layers[i]
        for generator in _moving_layer_screens(layer)
            screen_index += 1
            screen = generator.state.phase_rad
            m = size(screen, 1)
            size(screen, 2) == m || throw(DimensionMismatchError(
                "moving atmosphere layer screen must be square",
            ))
            sampling_ratio = _moving_screen_sampling_ratio(layer, generator)
            footprint_scale = renderer.footprint_scale[i] * sampling_ratio
            footprint_scale * T(n - 1) <= T(m - 1) || throw(DimensionMismatchError(
                "moving atmosphere layer screen is smaller than the pupil footprint",
            ))
            screen_period = T(m)
            # Keep the same accumulated native-pixel displacement for both
            # bands. Convert to coarse pixels only when sampling that screen.
            host_offset_increment_x[screen_index] =
                layer.params.wind_velocity_x * dt / layer.params.sampling_m
            host_offset_increment_y[screen_index] =
                layer.params.wind_velocity_y * dt / layer.params.sampling_m
            center = T(m + 1) / T(2) - footprint_scale * T(n - 1) / T(2)
            host_base_start_x[screen_index] = _normalize_replay_start(
                center + renderer.shift_x[i] * sampling_ratio, screen_period,
            )
            host_base_start_y[screen_index] = _normalize_replay_start(
                center + renderer.shift_y[i] * sampling_ratio, screen_period,
            )
            host_screen_period[screen_index] = screen_period
            host_sampling_ratio[screen_index] = sampling_ratio
        end
    end

    offset_increment_x = _copy_moving_screen_replay_vector(
        target,
        host_offset_increment_x,
    )
    offset_increment_y = _copy_moving_screen_replay_vector(
        target,
        host_offset_increment_y,
    )
    base_start_x = _copy_moving_screen_replay_vector(
        target,
        host_base_start_x,
    )
    base_start_y = _copy_moving_screen_replay_vector(
        target,
        host_base_start_y,
    )
    screen_period = _copy_moving_screen_replay_vector(
        target,
        host_screen_period,
    )
    sampling_ratio = _copy_moving_screen_replay_vector(target, host_sampling_ratio)
    offset_x = allocate_device_array(target, T, screen_count)
    offset_y = allocate_device_array(target, T, screen_count)
    start_x = allocate_device_array(target, T, screen_count)
    start_y = allocate_device_array(target, T, screen_count)
    fill!(offset_x, zero(T))
    fill!(offset_y, zero(T))
    copyto!(start_x, base_start_x)
    copyto!(start_y, base_start_y)

    plan = _MovingScreenReplayPlan(
        offset_increment_x,
        offset_increment_y,
        base_start_x,
        base_start_y,
        screen_period,
        sampling_ratio,
    )
    state = _MovingScreenReplayState(offset_x, offset_y)
    workspace = _MovingScreenReplayWorkspace(start_x, start_y)
    return _PreparedMovingScreenReplay(plan, state, workspace)
end

function _enqueue_moving_screen_replay!(
    opd::AbstractMatrix{T},
    renderer::AtmosphereDirectionRenderer,
    atm::MultiLayerAtmosphere,
    replay::_PreparedMovingScreenReplay,
) where {T<:AbstractFloat}
    style = execution_style(opd)
    layers = atm.layers
    plan = replay.plan
    screen_count = length(plan.screen_period)
    state = replay.state
    workspace = replay.workspace
    launch_kernel_async!(
        style,
        advance_moving_screen_replay_kernel!,
        state.offset_x,
        state.offset_y,
        workspace.start_x,
        workspace.start_y,
        plan.offset_increment_x,
        plan.offset_increment_y,
        plan.base_start_x,
        plan.base_start_y,
        plan.screen_period,
        plan.sampling_ratio,
        screen_count;
        ndrange=screen_count,
    )

    n = size(opd, 1)
    screen_index = 0
    @inbounds for layer_index in eachindex(layers)
        layer = layers[layer_index]
        output_scale = T(layer.params.cn2_amplitude_scale) *
            T(layer.generator.params.opd_per_radian)
        for generator in _moving_layer_screens(layer)
            screen_index += 1
            screen = generator.state.phase_rad
            footprint_scale = T(renderer.footprint_scale[layer_index]) *
                _moving_screen_sampling_ratio(layer, generator)
            launch_kernel_async!(
                style,
                moving_layer_replay_kernel!,
                opd,
                screen,
                workspace.start_x,
                workspace.start_y,
                screen_index,
                footprint_scale,
                output_scale,
                screen_index != 1,
                n,
                size(screen, 1);
                ndrange=size(opd),
            )
        end
    end
    launch_kernel_async!(
        style,
        apply_atmosphere_pupil_kernel!,
        opd,
        renderer.pupil,
        length(opd);
        ndrange=length(opd),
    )
    return opd
end

function _reset_moving_screen_replay!(replay::_PreparedMovingScreenReplay)
    fill!(replay.state.offset_x, zero(eltype(replay.state.offset_x)))
    fill!(replay.state.offset_y, zero(eltype(replay.state.offset_y)))
    copyto!(replay.workspace.start_x, replay.plan.base_start_x)
    copyto!(replay.workspace.start_y, replay.plan.base_start_y)
    return replay
end

function _reset_multilayer_atmosphere!(
    atm::MultiLayerAtmosphere,
    rng::AbstractRNG,
)
    timeline = atm.state.timeline
    timeline.model_time = zero(timeline.model_time)
    timeline.sequence = UInt64(0)
    timeline.initialized = false
    @inbounds for layer in atm.layers
        layer.state.offset_x = zero(layer.state.offset_x)
        layer.state.offset_y = zero(layer.state.offset_y)
        layer.state.initialized = false
        layer.state.low_frequency_initialized = false
    end
    advance_by!(atm, zero(timeline.model_time), rng)
    return atm
end

function _preflight_atmosphere_replay_step(
    atm::MultiLayerAtmosphere,
    duration::Real,
)
    timeline = atm.state.timeline
    timeline.initialized || throw(AtmosphereTimeError(
        "moving-screen replay requires an initialized atmosphere",
    ))
    T = typeof(timeline.model_time)
    dt = _explicit_atmosphere_duration(duration, T)
    next_time = timeline.model_time + dt
    isfinite(next_time) || throw(AtmosphereTimeError(
        "atmosphere model time overflowed its numeric representation",
    ))
    timeline.sequence == typemax(UInt64) && throw(AtmosphereTimeError(
        "atmosphere epoch sequence is exhausted",
    ))
    return nothing
end
