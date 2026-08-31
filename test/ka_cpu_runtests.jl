using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
import AdaptiveOpticsSim.Optics: filter!
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Tomography
using AdaptiveOpticsSim.Ensembles
using FFTW
using LinearAlgebra
using Random
using SpecialFunctions
using TOML

BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

# The package exports only the user-facing API. This standalone coverage target
# shares reference helpers with the full test suite, so make internal test-only
# names available without expanding the public export list.
for name in names(AdaptiveOpticsSim; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

for name in names(Backends; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Backends, $(QuoteNode(name)))
    end
end

for name in names(Optics; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Optics, $(QuoteNode(name)))
    end
end

for name in names(Detectors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Detectors, $(QuoteNode(name)))
    end
end

for name in names(Atmospheres; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Atmospheres, $(QuoteNode(name)))
    end
end

for name in names(WavefrontSensors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) =
            getfield(WavefrontSensors, $(QuoteNode(name)))
    end
end

for name in names(Tomography; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Tomography, $(QuoteNode(name)))
    end
end

for name in names(Ensembles; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Ensembles, $(QuoteNode(name)))
    end
end

include("reference_harness.jl")
include("ka_cpu_style_fixture.jl")
include("ka_cpu_matrix.jl")
