struct MultiLayerParams{T<:AbstractFloat}
    r0_fractions::Vector{T}
    wind_speed::Vector{T}
    wind_direction::Vector{T}
    altitude::Vector{T}
    L0::T
end

mutable struct MultiLayerState{T<:AbstractFloat,A<:AbstractMatrix{T}}
    opd::A
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

    params = MultiLayerParams{T}(
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
    fill!(opd, zero(T))
    state = MultiLayerState{T, typeof(opd)}(opd)

    return MultiLayerAtmosphere(params, layers, state)
end

function advance!(atm::MultiLayerAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng())
    n = tel.params.resolution
    delta = tel.params.diameter / n
    dt = tel.params.sampling_time

    fill!(atm.state.opd, zero(eltype(atm.state.opd)))

    for (i, layer) in enumerate(atm.layers)
        advance!(layer, tel; rng=rng)
        dx = atm.params.wind_speed[i] * cosd(atm.params.wind_direction[i]) * dt / delta
        dy = atm.params.wind_speed[i] * sind(atm.params.wind_direction[i]) * dt / delta
        shifted = circshift(layer.state.opd, (round(Int, dy), round(Int, dx)))
        atm.state.opd .+= shifted
    end

    return atm
end

function propagate!(atm::MultiLayerAtmosphere, tel::Telescope)
    tel.state.opd .= atm.state.opd .* tel.state.pupil
    return tel
end
