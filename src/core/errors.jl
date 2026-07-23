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

"""Prepared-contract violation at one semantic wavefront-sensor stage."""
struct WFSPreparationError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

Base.showerror(io::IO, e::AdaptiveOpticsSimError) = print(io, e.msg)
