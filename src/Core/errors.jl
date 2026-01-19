"""Base exception type for AdaptiveOptics.jl."""
abstract type AdaptiveOpticsError <: Exception end

struct InvalidConfiguration <: AdaptiveOpticsError
    msg::String
end

struct DimensionMismatchError <: AdaptiveOpticsError
    msg::String
end

Base.showerror(io::IO, e::AdaptiveOpticsError) = print(io, e.msg)
