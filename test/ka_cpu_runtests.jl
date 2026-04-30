using Test
using AdaptiveOpticsSim
using FFTW
using KernelAbstractions
using LinearAlgebra
using Random
using SpecialFunctions
using TOML

# The package exports only the user-facing API. This standalone coverage target
# shares reference helpers with the full test suite, so make internal test-only
# names available without expanding the public export list.
for name in names(AdaptiveOpticsSim; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

include("reference_harness.jl")
include("ka_cpu_matrix.jl")
