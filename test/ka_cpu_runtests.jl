using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim: Plant
using AdaptiveOpticsSim.Plant
using FFTW
using LinearAlgebra
using Random
using SpecialFunctions
using TOML

BLAS.set_num_threads(1)
AdaptiveOpticsSim.set_fft_provider_threads!(1)

# The package exports only the user-facing API. This standalone coverage target
# shares reference helpers with the full test suite, so make internal test-only
# names available without expanding the public export list.
for name in names(AdaptiveOpticsSim; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end


for name in names(Plant; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Plant, $(QuoteNode(name)))
    end
end

include("reference_harness.jl")
include("ka_cpu_style_fixture.jl")
include("ka_cpu_matrix.jl")
