"""Base exception type for AdaptiveOpticsSim.jl."""
abstract type AdaptiveOpticsSimError <: Exception end

struct InvalidConfiguration <: AdaptiveOpticsSimError
    msg::String
end

struct DimensionMismatchError <: AdaptiveOpticsSimError
    msg::String
end

struct UnsupportedAlgorithm <: AdaptiveOpticsSimError
    msg::String
end

struct NumericalConditionError <: AdaptiveOpticsSimError
    msg::String
end

"""Invalid or non-monotonic explicit atmosphere model time."""
struct AtmosphereTimeError <: AdaptiveOpticsSimError
    msg::String
end

"""Missing, stale, or incompatible atmosphere epoch identity."""
struct AtmosphereEpochError <: AdaptiveOpticsSimError
    msg::String
end

"""Invalid canonical plant timestamp, duration, conversion, or arithmetic."""
struct PlantTimeError <: AdaptiveOpticsSimError
    operation::Symbol
    reason::Symbol
    msg::String
end

"""Invalid periodic schedule, event generator, or deterministic scheduler operation."""
struct PlantScheduleError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Prepared-contract violation at one semantic wavefront-sensor stage."""
struct WFSPreparationError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

"""Invalid immutable plant, optical-path, or acquisition declaration."""
struct PlantDefinitionError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Invalid preparation or prepared binding for a plant component."""
struct PlantPreparationError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

"""Invalid prepared detector-acquisition event or lifecycle transition."""
struct DetectorAcquisitionError <: AdaptiveOpticsSimError
    component::Symbol
    reason::Symbol
    msg::String
end

Base.showerror(io::IO, e::AdaptiveOpticsSimError) = print(io, e.msg)
