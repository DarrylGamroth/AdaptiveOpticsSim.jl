@inline wrap_upper_index(i::Int, n::Int) = ifelse(i > n, i - n, i)

@kernel function moving_layer_extract_kernel!(out, screen, start_x, start_y, footprint_scale, output_scale, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = start_y + footprint_scale * T(i - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = start_x + footprint_scale * T(j - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
            out[i, j] = output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

@kernel function moving_layer_accumulate_kernel!(out, screen, start_x, start_y, footprint_scale, output_scale, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = start_y + footprint_scale * T(i - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = start_x + footprint_scale * T(j - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
            out[i, j] += output_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

@kernel function moving_layer_replay_kernel!(
    out,
    screen,
    start_x,
    start_y,
    layer_index::Int,
    footprint_scale,
    output_scale,
    accumulate::Bool,
    n::Int,
    m::Int,
)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = @inbounds(start_y[layer_index]) + footprint_scale * T(i - 1)
        y0 = unsafe_trunc(Int, floor(y))
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_upper_index(y0, m)
        iy1 = wrap_upper_index(y0 + 1, m)

        x = @inbounds(start_x[layer_index]) + footprint_scale * T(j - 1)
        x0 = unsafe_trunc(Int, floor(x))
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_upper_index(x0, m)
        ix1 = wrap_upper_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
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
    layer_count::Int,
)
    i = @index(Global, Linear)
    if i <= layer_count
        @inbounds begin
            next_offset_x = offset_x[i] + offset_increment_x[i]
            next_offset_y = offset_y[i] + offset_increment_y[i]
            offset_x[i] = next_offset_x
            offset_y[i] = next_offset_y
            start_x[i] = _normalize_replay_start(
                base_start_x[i] - next_offset_x,
                screen_period[i],
            )
            start_y[i] = _normalize_replay_start(
                base_start_y[i] - next_offset_y,
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
end

struct MovingAtmosphereLayer{
    P<:MovingLayerParams,
    S<:MovingLayerState,
    A<:KolmogorovAtmosphere,
    TT<:Telescope,
} <: AbstractAtmosphereLayer
    params::P
    generator::A
    generator_telescope::TT
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
end

"""Persistent device-resident moving-screen offsets for one graph owner."""
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
    state = MovingLayerState(zero(T), zero(T), false)
    return MovingAtmosphereLayer(params, generator, screen_telescope, state)
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
    @inbounds for i in 1:n
        y = start_y + footprint_scale * T(i - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        for j in 1:n
            x = start_x + footprint_scale * T(j - 1)
            x0 = floor(Int, x)
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_index(x0, m)
            ix1 = wrap_index(x0 + 1, m)

            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
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
    @inbounds for i in 1:n
        y = start_y + footprint_scale * T(i - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        for j in 1:n
            x = start_x + footprint_scale * T(j - 1)
            x0 = floor(Int, x)
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_index(x0, m)
            ix1 = wrap_index(x0 + 1, m)

            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
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

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    extract_shifted_screen_async!(out, layer.generator.state.phase_rad,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        opd_scale,
        footprint_scale)
    return out
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer_accumulate!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    accumulate_shifted_screen!(out, layer.generator.state.phase_rad,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        opd_scale,
        footprint_scale)
    return out
end

function _render_layer_accumulate_async!(out::AbstractMatrix{T},
    layer::MovingAtmosphereLayer, shift_x::T, shift_y::T,
    footprint_scale::T) where {T<:AbstractFloat}
    opd_scale = T(layer.params.cn2_amplitude_scale) *
        T(layer.generator.params.opd_per_radian)
    _accumulate_shifted_screen_async!(out, layer.generator.state.phase_rad,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        opd_scale,
        footprint_scale)
    return out
end

function evolve_layer!(layer::MovingAtmosphereLayer, duration::Real)
    T = typeof(layer.state.offset_x)
    dt = T(duration)
    delta = layer.params.sampling_m
    layer.state.offset_x += layer.params.wind_velocity_x * dt / delta
    layer.state.offset_y += layer.params.wind_velocity_y * dt / delta
    return layer
end

function MultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    reference_wavelength_m::Real,
    L0::Real=25.0,
    fractional_cn2::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction_deg::AbstractVector,
    altitude::AbstractVector,
    layer_ids=nothing,
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
    return atm
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
    layer_count = length(layers)
    n = size(renderer.pupil, 1)
    target = compute_device(renderer.pupil)
    host_offset_increment_x = Vector{T}(undef, layer_count)
    host_offset_increment_y = Vector{T}(undef, layer_count)
    host_base_start_x = Vector{T}(undef, layer_count)
    host_base_start_y = Vector{T}(undef, layer_count)
    host_screen_period = Vector{T}(undef, layer_count)

    @inbounds for i in eachindex(layers)
        layer = layers[i]
        screen = layer.generator.state.phase_rad
        m = size(screen, 1)
        size(screen, 2) == m || throw(DimensionMismatchError(
            "moving atmosphere layer screen must be square",
        ))
        m >= n || throw(DimensionMismatchError(
            "moving atmosphere layer screen is smaller than the pupil",
        ))
        screen_period = T(m)
        footprint_scale = renderer.footprint_scale[i]
        host_offset_increment_x[i] =
            layer.params.wind_velocity_x * dt / layer.params.sampling_m
        host_offset_increment_y[i] =
            layer.params.wind_velocity_y * dt / layer.params.sampling_m
        host_base_start_x[i] = _normalize_replay_start(
            T(m + 1) / T(2) - footprint_scale * T(n - 1) / T(2) +
            renderer.shift_x[i],
            screen_period,
        )
        host_base_start_y[i] = _normalize_replay_start(
            T(m + 1) / T(2) - footprint_scale * T(n - 1) / T(2) +
            renderer.shift_y[i],
            screen_period,
        )
        host_screen_period[i] = screen_period
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
    offset_x = allocate_device_array(target, T, layer_count)
    offset_y = allocate_device_array(target, T, layer_count)
    start_x = allocate_device_array(target, T, layer_count)
    start_y = allocate_device_array(target, T, layer_count)
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
    layer_count = length(layers)
    plan = replay.plan
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
        layer_count;
        ndrange=layer_count,
    )

    n = size(opd, 1)
    @inbounds for layer_index in eachindex(layers)
        layer = layers[layer_index]
        screen = layer.generator.state.phase_rad
        m = size(screen, 1)
        output_scale = T(layer.params.cn2_amplitude_scale) *
            T(layer.generator.params.opd_per_radian)
        launch_kernel_async!(
            style,
            moving_layer_replay_kernel!,
            opd,
            screen,
            workspace.start_x,
            workspace.start_y,
            layer_index,
            T(renderer.footprint_scale[layer_index]),
            output_scale,
            layer_index != firstindex(layers),
            n,
            m;
            ndrange=size(opd),
        )
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
