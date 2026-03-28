@kernel function moving_layer_extract_kernel!(out, screen, start_x, start_y, scale, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        T = eltype(out)
        y = start_y + T(i - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        x = start_x + T(j - 1)
        x0 = floor(Int, x)
        fx = x - T(x0)
        wx0 = one(T) - fx
        ix0 = wrap_index(x0, m)
        ix1 = wrap_index(x0 + 1, m)

        @inbounds begin
            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
            out[i, j] = scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
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
}
    params::P
    generator::A
    generator_telescope::TT
    state::S
end

mutable struct MultiLayerState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
    layer_buffer::A
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
    backend=Array,
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
    T::Type{<:AbstractFloat}=Float64,
    backend=Array,
)
    screen_resolution = moving_layer_screen_resolution(tel.params.resolution)
    screen_telescope = moving_layer_telescope(tel; resolution=screen_resolution, T=T, backend=backend)
    generator = KolmogorovAtmosphere(screen_telescope; r0=r0, L0=L0, T=T, backend=backend)
    params = MovingLayerParams(T(sqrt(cn2_fraction)), T(wind_velocity_x), T(wind_velocity_y))
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

function extract_shifted_screen!(out::AbstractMatrix{T}, screen::AbstractMatrix{T}, offset_x::T, offset_y::T, scale::T) where {T<:AbstractFloat}
    n = size(out, 1)
    size(out, 2) == n || throw(DimensionMismatchError("output must be square"))
    m = size(screen, 1)
    size(screen, 2) == m || throw(DimensionMismatchError("screen must be square"))
    m >= n || throw(DimensionMismatchError("screen resolution must be at least as large as the pupil resolution"))

    start_x = T((m - n) / 2 + 1) - offset_x
    start_y = T((m - n) / 2 + 1) - offset_y
    _extract_shifted_screen!(execution_style(out), out, screen, start_x, start_y, scale, n, m)
    return out
end

function _extract_shifted_screen!(::ScalarCPUStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    # Move a finite periodic canvas under the pupil with bilinear subpixel interpolation.
    @inbounds for i in 1:n
        y = start_y + T(i - 1)
        y0 = floor(Int, y)
        fy = y - T(y0)
        wy0 = one(T) - fy
        iy0 = wrap_index(y0, m)
        iy1 = wrap_index(y0 + 1, m)

        for j in 1:n
            x = start_x + T(j - 1)
            x0 = floor(Int, x)
            fx = x - T(x0)
            wx0 = one(T) - fx
            ix0 = wrap_index(x0, m)
            ix1 = wrap_index(x0 + 1, m)

            v00 = screen[iy0, ix0]
            v01 = screen[iy0, ix1]
            v10 = screen[iy1, ix0]
            v11 = screen[iy1, ix1]
            out[i, j] = scale * (wy0 * (wx0 * v00 + fx * v01) + fy * (wx0 * v10 + fx * v11))
        end
    end
    return out
end

function _extract_shifted_screen!(style::AcceleratorStyle, out::AbstractMatrix{T}, screen::AbstractMatrix{T},
    start_x::T, start_y::T, scale::T, n::Int, m::Int) where {T<:AbstractFloat}
    launch_kernel!(style, moving_layer_extract_kernel!, out, screen, start_x, start_y, scale, n, m; ndrange=size(out))
    return out
end

function sample_layer!(out::AbstractMatrix{T}, layer::MovingAtmosphereLayer, tel::Telescope, rng::AbstractRNG) where {T<:AbstractFloat}
    ensure_initialized!(layer, rng)
    delta = T(tel.params.diameter / tel.params.resolution)
    dt = T(tel.params.sampling_time)
    layer.state.offset_x += T(layer.params.wind_velocity_x) * dt / delta
    layer.state.offset_y += T(layer.params.wind_velocity_y) * dt / delta
    scale = T(layer.params.amplitude_scale)
    extract_shifted_screen!(out, layer.generator.state.opd, layer.state.offset_x, layer.state.offset_y, scale)
    return out
end

function MultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    L0::Real=25.0,
    fractional_cn2::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction::AbstractVector,
    altitude::AbstractVector,
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)

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
            T=T,
            backend=backend,
        ) for i in 1:n_layers
    ]

    opd = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    layer_buffer = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    fill!(opd, zero(T))
    fill!(layer_buffer, zero(T))
    state = MultiLayerState{T, typeof(opd)}(opd, layer_buffer)

    return MultiLayerAtmosphere(params, layers, state)
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope, rng::AbstractRNG)
    fill!(atm.state.opd, zero(eltype(atm.state.opd)))

    for layer in atm.layers
        sample_layer!(atm.state.layer_buffer, layer, tel, rng)
        atm.state.opd .+= atm.state.layer_buffer
    end

    return atm
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    return advance!(atm, tel, rng)
end

function propagate!(atm::MultiLayerAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
