struct MultiLayerParams{T<:AbstractFloat,
    V1<:AbstractVector{T},
    V2<:AbstractVector{T},
    V3<:AbstractVector{T},
    V4<:AbstractVector{T}}
    r0_fractions::V1
    wind_speed::V2
    wind_direction::V3
    altitude::V4
    L0::T
end

mutable struct MultiLayerState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
    shift_buffer::A
end

struct MultiLayerAtmosphere{P<:MultiLayerParams,S<:MultiLayerState} <: AbstractAtmosphere
    params::P
    layers::Vector{KolmogorovAtmosphere}
    state::S
end

function MultiLayerAtmosphere(tel::Telescope;
    r0::Real,
    L0::Real=25.0,
    fractional_r0::AbstractVector,
    wind_speed::AbstractVector,
    wind_direction::AbstractVector,
    altitude::AbstractVector,
    T::Type{<:AbstractFloat}=Float64,
    backend=Array)

    n_layers = length(fractional_r0)
    if length(wind_speed) != n_layers || length(wind_direction) != n_layers || length(altitude) != n_layers
        throw(InvalidConfiguration("layer parameter lengths must match fractional_r0"))
    end

    params = MultiLayerParams(
        T.(fractional_r0),
        T.(wind_speed),
        T.(wind_direction),
        T.(altitude),
        T(L0),
    )

    layers = Vector{KolmogorovAtmosphere}(undef, n_layers)
    for i in 1:n_layers
        r0_layer = T(r0) / params.r0_fractions[i]
        layers[i] = KolmogorovAtmosphere(tel; r0=r0_layer, L0=L0, T=T, backend=backend)
    end

    opd = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    shift_buffer = backend{T}(undef, tel.params.resolution, tel.params.resolution)
    fill!(opd, zero(T))
    fill!(shift_buffer, zero(T))
    state = MultiLayerState{T, typeof(opd)}(opd, shift_buffer)

    return MultiLayerAtmosphere(params, layers, state)
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope, rng::AbstractRNG)
    n = tel.params.resolution
    delta = tel.params.diameter / n
    dt = tel.params.sampling_time

    fill!(atm.state.opd, zero(eltype(atm.state.opd)))

    for (i, layer) in enumerate(atm.layers)
        advance!(layer, tel, rng)
        dx = atm.params.wind_speed[i] * cosd(atm.params.wind_direction[i]) * dt / delta
        dy = atm.params.wind_speed[i] * sind(atm.params.wind_direction[i]) * dt / delta
        circshift!(atm.state.shift_buffer, layer.state.opd, (round(Int, dy), round(Int, dx)))
        atm.state.opd .+= atm.state.shift_buffer
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
