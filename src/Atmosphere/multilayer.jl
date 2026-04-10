@inline wrap_upper_index(i::Int, n::Int) = ifelse(i > n, i - n, i)

@kernel function moving_layer_extract_kernel!(out, screen, start_x, start_y, footprint_scale, amplitude_scale, n::Int, m::Int)
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
            out[i, j] = amplitude_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

@kernel function moving_layer_accumulate_kernel!(out, screen, start_x, start_y, footprint_scale, amplitude_scale, n::Int, m::Int)
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
            out[i, j] += amplitude_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
end

struct MultiLayerParams{T<:AbstractFloat,
    V1<:AbstractVector{T},
    V2<:AbstractVector{T},
    V3<:AbstractVector{T},
    V4<:AbstractVector{T},
    V5<:AbstractVector{T},
    V6<:AbstractVector{T}}
    cn2_fractions::V1
    wind_speed::V2
    wind_direction::V3
    altitude::V4
    wind_velocity_x::V5
    wind_velocity_y::V6
    r0::T
    L0::T
end

struct MovingLayerParams{T<:AbstractFloat}
    amplitude_scale::T
    wind_velocity_x::T
    wind_velocity_y::T
    altitude::T
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

mutable struct MultiLayerState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
    source_geometry::AtmosphereSourceGeometryCache{T,Vector{T}}
end

struct MultiLayerAtmosphere{
    P<:MultiLayerParams,
    S<:MultiLayerState,
    L<:AbstractVector{<:MovingAtmosphereLayer},
} <: AbstractAtmosphere
    params::P
    layers::L
    state::S
end

@inline moving_layer_screen_resolution(n::Int) = 3 * n

function moving_layer_telescope(
    tel::Telescope;
    resolution::Int,
    T::Type{<:AbstractFloat}=Float64,
    backend=CPUBackend(),
)
    delta = tel.params.diameter / tel.params.resolution
    return Telescope(
        resolution=resolution,
        diameter=delta * resolution,
        sampling_time=tel.params.sampling_time,
        central_obstruction=0.0,
        fov_arcsec=tel.params.fov_arcsec,
        T=T,
        backend=backend,
    )
end

function MovingAtmosphereLayer(
    tel::Telescope;
    r0::Real,
    L0::Real,
    cn2_fraction::Real,
    wind_velocity_x::Real,
    wind_velocity_y::Real,
    altitude::Real,
    T::Type{<:AbstractFloat}=Float64,
    backend=CPUBackend(),
)
    screen_resolution = moving_layer_screen_resolution(tel.params.resolution)
    screen_telescope = moving_layer_telescope(tel; resolution=screen_resolution, T=T, backend=backend)
    generator = KolmogorovAtmosphere(screen_telescope; r0=r0, L0=L0, T=T, backend=backend)
    params = MovingLayerParams(T(sqrt(cn2_fraction)), T(wind_velocity_x), T(wind_velocity_y), T(altitude))
    state = MovingLayerState(zero(T), zero(T), false)
    return MovingAtmosphereLayer(params, generator, screen_telescope, state)
end

function ensure_initialized!(layer::MovingAtmosphereLayer, rng::AbstractRNG)
    if !layer.state.initialized
        advance!(layer.generator, layer.generator_telescope, rng)
        layer.state.initialized = true
    end
    return layer
end

@inline wrap_index(i::Int, n::Int) = mod1(i, n)
@inline normalize_start_coordinate(start::T, m::Int) where {T<:AbstractFloat} = mod(start - one(T), T(m)) + one(T)

function extract_shifted_screen!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, amplitude_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
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
    _extract_shifted_screen!(execution_style(out), out, screen, start_x, start_y, footprint_scale, amplitude_scale, n, m)
    return out
end

function extract_shifted_screen_async!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, amplitude_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration("footprint_scale must be positive"))
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    _extract_shifted_screen_async!(execution_style(out), out, screen, start_x, start_y, footprint_scale, amplitude_scale, n, m)
    return out
end

function _extract_shifted_screen!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
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
            out[i, j] = amplitude_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
    return out
end

function _extract_shifted_screen!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel!(style, moving_layer_extract_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, amplitude_scale, n, m; ndrange=size(out))
    return out
end

function accumulate_shifted_screen!(out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    offset_x::T, offset_y::T, amplitude_scale::T, footprint_scale::T=one(T)) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))
    footprint_scale > zero(T) || throw(InvalidConfiguration("footprint_scale must be positive"))
    start_x = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_x
    start_y = T(m + 1) / 2 - footprint_scale * T(n - 1) / 2 - offset_y
    _accumulate_shifted_screen!(execution_style(out), out, screen, start_x, start_y, footprint_scale, amplitude_scale, n, m)
    return out
end

@inline function _extract_shifted_screen_async!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    return _extract_shifted_screen!(ScalarCPUStyle(), out, screen, start_x, start_y, footprint_scale, amplitude_scale, n, m)
end

function _extract_shifted_screen_async!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel_async!(style, moving_layer_extract_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, amplitude_scale, n, m; ndrange=size(out))
    return out
end

function _accumulate_shifted_screen!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
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
            out[i, j] += amplitude_scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
    return out
end

function _accumulate_shifted_screen!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, footprint_scale::T, amplitude_scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    start_x_wrapped = normalize_start_coordinate(start_x, m)
    start_y_wrapped = normalize_start_coordinate(start_y, m)
    launch_kernel!(style, moving_layer_accumulate_kernel!, out, screen, start_x_wrapped, start_y_wrapped, footprint_scale, amplitude_scale, n, m; ndrange=size(out))
    return out
end

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    amplitude_scale = T(layer.params.amplitude_scale)
    extract_shifted_screen_async!(out, layer.generator.state.opd,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        amplitude_scale,
        footprint_scale)
    return out
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope,
    src::Union{AbstractSource,Nothing}=nothing) where {T<:AbstractFloat}
    return render_layer_accumulate!(out, layer, layer_render_context(src, layer, tel, T))
end

function render_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer,
    shift_x::T, shift_y::T, footprint_scale::T) where {T<:AbstractFloat}
    amplitude_scale = T(layer.params.amplitude_scale)
    accumulate_shifted_screen!(out, layer.generator.state.opd,
        layer.state.offset_x - shift_x,
        layer.state.offset_y - shift_y,
        amplitude_scale,
        footprint_scale)
    return out
end

function sample_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope, rng::AbstractRNG) where {T<:AbstractFloat}
    ensure_initialized!(layer, rng)
    delta = T(tel.params.diameter / tel.params.resolution)
    dt = T(tel.params.sampling_time)
    layer.state.offset_x += T(layer.params.wind_velocity_x) * dt / delta
    layer.state.offset_y += T(layer.params.wind_velocity_y) * dt / delta
    return render_layer!(out, layer, tel)
end

function sample_layer_accumulate!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope, rng::AbstractRNG) where {T<:AbstractFloat}
    ensure_initialized!(layer, rng)
    delta = T(tel.params.diameter / tel.params.resolution)
    dt = T(tel.params.sampling_time)
    layer.state.offset_x += T(layer.params.wind_velocity_x) * dt / delta
    layer.state.offset_y += T(layer.params.wind_velocity_y) * dt / delta
    return render_layer_accumulate!(out, layer, tel)
end

function MultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    L0::Real=25.0,
    fractional_cn2::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction::AbstractVector,
    altitude::AbstractVector,
    T::Type{<:AbstractFloat}=Float64,
    backend=CPUBackend())

    backend = _resolve_array_backend(backend)

    n_layers = length(fractional_cn2)
    n_layers > 0 || throw(InvalidConfiguration("fractional_cn2 cannot be empty"))
    if length(wind_speed) != n_layers || length(wind_direction) != n_layers || length(altitude) != n_layers
        throw(InvalidConfiguration("layer parameter lengths must match fractional_cn2"))
    end
    all(>=(0), fractional_cn2) || throw(InvalidConfiguration("fractional_cn2 must be non-negative"))
    isapprox(sum(fractional_cn2), 1; atol=1e-6, rtol=1e-6) ||
        throw(InvalidConfiguration("fractional_cn2 must sum to 1"))

    params = MultiLayerParams(
        T.(fractional_cn2),
        T.(wind_speed),
        T.(wind_direction),
        T.(altitude),
        T[T(wind_speed[i]) * cosd(T(wind_direction[i])) for i in 1:n_layers],
        T[T(wind_speed[i]) * sind(T(wind_direction[i])) for i in 1:n_layers],
        T(r0),
        T(L0),
    )

    layers = [
        MovingAtmosphereLayer(
            tel;
            r0=r0,
            L0=L0,
            cn2_fraction=params.cn2_fractions[i],
            wind_velocity_x=params.wind_velocity_x[i],
            wind_velocity_y=params.wind_velocity_y[i],
            altitude=params.altitude[i],
            T=T,
            backend=backend,
        ) for i in 1:n_layers
    ]

    opd = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    fill!(opd, zero(T))
    state = MultiLayerState{T, typeof(opd)}(opd, AtmosphereSourceGeometryCache(n_layers, T))

    return MultiLayerAtmosphere(params, layers, state)
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope, rng::AbstractRNG)
    accumulate_sampled_layers!(atm.state.opd, atm.layers, tel, rng)
    return atm
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    return advance!(atm, tel, rng)
end

function propagate!(atm::MultiLayerAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end

function propagate!(atm::MultiLayerAtmosphere, tel::Telescope, src::AbstractSource)
    return propagate_source_aware!(atm, tel, src)
end
