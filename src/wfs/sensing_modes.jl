abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric <: SensingMode end

sensing_mode(::AbstractWFS) = Diffractive()

@inline deterministic_frame_readout_gain(::CCDSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::CMOSSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::InGaAsSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::EMCCDSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = one(T)
@inline deterministic_frame_readout_gain(::AbstractHgCdTeSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(
    sensor::HgCdTeAvalancheArraySensor, gain,
    ::Type{T}) where {T<:AbstractFloat} =
    T(sensor.avalanche_gain) * T(gain)
