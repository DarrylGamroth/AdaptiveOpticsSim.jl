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

Base.showerror(io::IO, e::AdaptiveOpticsSimError) = print(io, e.msg)
