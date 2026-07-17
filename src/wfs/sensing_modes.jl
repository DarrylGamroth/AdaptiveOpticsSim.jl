abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric <: SensingMode end

abstract type WFSNormalization end

"""Normalize a WFS signal by the mean measured signal over valid samples."""
struct MeanValidFluxNormalization <: WFSNormalization end

"""
Normalize a WFS signal by the expected incident photon signal per pupil sample.

Detector-coupled measurements convert the denominator to the detector's
deterministic linear photon-detection units. Frame detectors include exposure,
QE, deterministic detector readout gain, and HgCdTe avalanche gain when
applicable. Counting detectors include exposure, QE, gate, fill factor, source
throughput, and post-detection gain.
Detector noise is never applied to a calibration reference. A zero incident or
detectable signal follows the WFS no-signal convention and returns finite zero
output.
"""
struct IncidenceFluxNormalization <: WFSNormalization end

sensing_mode(::AbstractWFS) = Diffractive()

@inline wfs_incident_photon_irradiance(src::AbstractSource, ::Type{T}) where {T<:AbstractFloat} =
    T(photon_irradiance(src))

function wfs_incident_photon_irradiance(ast::Asterism, ::Type{T}) where {T<:AbstractFloat}
    total = zero(T)
    @inbounds for src in ast.sources
        total += T(photon_irradiance(src))
    end
    return total
end

@inline wfs_detector_incidence_scale(::Nothing, ::AbstractSource,
    ::Type{T}) where {T<:AbstractFloat} = one(T)

@inline deterministic_frame_readout_gain(::CCDSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::CMOSSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::InGaAsSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = T(gain)
@inline deterministic_frame_readout_gain(::EMCCDSensor, gain,
    ::Type{T}) where {T<:AbstractFloat} = one(T)
@inline deterministic_frame_readout_gain(
    sensor::HgCdTeAvalancheArraySensor, gain,
    ::Type{T}) where {T<:AbstractFloat} =
    T(sensor.avalanche_gain) * T(gain)

@inline function wfs_detector_incidence_scale(det::Detector,
    src::AbstractSource, ::Type{T}) where {T<:AbstractFloat}
    gain = deterministic_frame_readout_gain(det.params.sensor,
        det.params.gain, T)
    return T(det.params.integration_time) * effective_qe(det, src, T) * gain
end

@inline function wfs_detector_incidence_scale(det::AbstractCountingDetector,
    src::AbstractSource, ::Type{T}) where {T<:AbstractFloat}
    return T(counting_exposure_time(det)) * counting_qe(det, T) *
        counting_fill_factor(det, T) * counting_source_throughput(det, src, T) *
        counting_post_gain(det, T)
end

@inline usable_wfs_normalization(value::T) where {T<:AbstractFloat} =
    isfinite(value) && value > eps(T)
