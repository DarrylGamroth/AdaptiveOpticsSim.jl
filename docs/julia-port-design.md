# Julia Design Sketch for AdaptiveOptics.jl

This document sketches an idiomatic Julia design for porting OOPAO into
AdaptiveOptics.jl. It favors
multiple dispatch, explicit state, and preallocated workspaces over OO-style
inheritance and hidden mutation.

## Goals
- Preserve OOPAO capabilities with a Julia-centric API.
- Make hot paths allocation-free and type-stable.
- Enable GPU backends without rewriting core algorithms.
- Keep configuration explicit and reproducible.

## Non-goals
- 1:1 translation of Python classes and operator overloading.
- Implicit global state or hidden mutability.

## Package layout (proposed)
- `AdaptiveOptics.jl` root module
- `Core/` interfaces, shared types, traits, workspaces
- `Optics/` telescope, DM, Zernike, propagation utilities
- `Atmosphere/` phase screens, turbulence models
- `WFS/` Shack-Hartmann, Pyramid, shared WFS utilities
- `Calibration/` interaction matrices, reconstructors
- `Sim/` optical chains and loop orchestration
- `IO/` optional format helpers, configuration, serialization
- `Viz/` thin plotting helpers (optional)

## Core interfaces
Use abstract types to declare interfaces, then implement functionality through
multiple dispatch and small trait functions.

```julia
abstract type AbstractOpticalElement end
abstract type AbstractSource <: AbstractOpticalElement end
abstract type AbstractAtmosphere <: AbstractOpticalElement end
abstract type AbstractWFS <: AbstractOpticalElement end
abstract type AbstractDetector <: AbstractOpticalElement end
```

Concrete types (e.g., `Telescope`, `KolmogorovAtmosphere`, `ShackHartmann`) live
in domain modules and implement these interfaces via methods.

## Traits
Traits avoid flag-heavy structs and let dispatch pick the right algorithm.

```julia
abstract type SensingMode end
struct Diffractive <: SensingMode end
struct Geometric  <: SensingMode end

sensing_mode(::AbstractWFS) = Diffractive()

# WFS types carry the sensing mode as a type parameter (e.g.,
# `ShackHartmann{Geometric}`), not as a field.

abstract type NoiseModel end
struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end
struct NoiseReadout <: NoiseModel end
struct NoisePhotonReadout <: NoiseModel end

abstract type DMApplyMode end
struct DMAdditive <: DMApplyMode end
struct DMReplace <: DMApplyMode end
```

Example: `sense!(::Geometric, wfs, tel, ws)` uses gradients only; diffractive
uses a propagated intensity model.

Noise and apply-mode logic follow the same pattern:
- `Detector{NoiseNone}`/`Detector{NoisePhoton}`/... are set by the `noise=` argument
  (accepts either a single `NoiseModel` or a tuple like `(NoisePhoton(), NoiseReadout(0.5))`).
- `apply!(dm, tel, DMAdditive())` and `apply!(dm, tel, DMReplace())` dispatch
  without runtime conditionals.

## Dataflow and propagation
Pipeline modeled as explicit functions on an `OpticalChain` or direct calls.

```julia
propagate!(ws, src::AbstractSource, tel::Telescope) = # ...
propagate!(ws, atm::AbstractAtmosphere, tel::Telescope) = # ...
propagate!(ws, wfs::AbstractWFS, tel::Telescope) = sense!(ws, sensing_mode(wfs), wfs, tel)

struct OpticalChain
    elements::Vector{AbstractOpticalElement}
end
```

## State and workspace
- Parameter structs are immutable; evolving data lives in `mutable struct` state.
- A `Workspace` holds scratch buffers sized from `TelescopeParams`.
- All hot-path operations use `!`-mutating functions.

## Parametric types and backends
Parametrize by numeric type and array backend to support `Float32`/`Float64`
and CPU/GPU storage without abstract fields.

Guidelines:
- Keep parameters minimal (numeric type and concrete array types).
- Prefer concrete array fields over `AbstractArray` fields in structs.
- Accept `AbstractArray` in method signatures for flexibility.

Sketch:
```julia
struct TelescopeParams{T<:AbstractFloat}
    resolution::Int
    diameter::T
    sampling_time::T
    central_obstruction::T
    fov_arcsec::T
end

mutable struct TelescopeState{Aopd<:AbstractMatrix{T}, Apsf<:AbstractMatrix{T}, Amask<:AbstractMatrix{Bool}, T}
    pupil::Amask
    opd::Aopd
    psf::Apsf
end

struct Telescope{P,S}
    params::P
    state::S
end

function Telescope(; T=Float64, backend=Array, kwargs...)
    params = TelescopeParams{T}(; kwargs...)
    n = params.resolution
    pupil = backend{Bool}(undef, n, n)
    opd   = backend{T}(undef, n, n)
    psf   = backend{T}(undef, n, n)
    state = TelescopeState{typeof(opd), typeof(psf), typeof(pupil), T}(pupil, opd, psf)
    return Telescope(params, state)
end
```

This keeps the API consistent while allowing `backend=CuArray` for GPU runs.

## Performance guidelines
- Preallocate arrays in `Workspace` and reuse them in every loop iteration.
- Avoid captured variables in closures on hot paths.
- Use `@views` and `mul!` where possible; avoid temporary arrays.

## Parallelism strategy
Use a centralized parallel configuration and focus on coarse-grained parallelism
first (independent sources, time steps, Monte Carlo runs, parameter sweeps).
Avoid nested parallelism and oversubscription.

Guidelines:
- Prefer outer-loop threading; keep inner FFT/WFS kernels single-threaded unless
  profiling proves otherwise.
- Coordinate CPU thread counts with FFTW/BLAS threads.
- For deterministic mode, set all thread counts to 1.

Sketch:
```julia
struct ParallelConfig
    threads::Int
    fftw_threads::Int
    blas_threads::Int
end

function with_parallel_config(cfg::ParallelConfig, f)
    Base.Threads.nthreads() == cfg.threads || @warn "start Julia with JULIA_NUM_THREADS"
    FFTW.set_num_threads(cfg.fftw_threads)
    BLAS.set_num_threads(cfg.blas_threads)
    return f()
end

with_parallel_config(ParallelConfig(8, 2, 2)) do
    Threads.@threads :static for i in eachindex(tasks)
        run_task!(tasks[i])
    end
end
```

## Benchmarking strategy
Maintain a small benchmark suite with standardized scenarios to track runtime,
allocations, and GPU/CPU parity.

Suggested benchmarks:
- PSF generation (single source, fixed zero-padding).
- WFS slope computation (fixed subaperture count).
- Reconstructor application (matrix-vector and modal reconstructor).
- Closed-loop step (single iteration, deterministic inputs).

Metrics:
- Wall time, allocations, and memory footprint.
- CPU vs GPU speedup and accuracy deltas.

Tools:
- `BenchmarkTools.jl` for CPU.
- `CUDA.@benchmark` for GPU kernels.

## Tabular outputs (Tables.jl)
Tables.jl is useful for exporting per-iteration metrics without coupling the
core to DataFrames. Keep pixel data as arrays; expose summary telemetry via a
Tables-compatible interface.

Sketch:
```julia
struct TelemetryRow
    iter::Int
    t::Float64
    wfe_rms::Float64
    strehl::Float64
    loop_gain::Float64
end

struct Telemetry
    rows::Vector{TelemetryRow}
end

Tables.istable(::Type{Telemetry}) = true
Tables.rowaccess(::Type{Telemetry}) = true
Tables.rows(t::Telemetry) = t.rows
```

This keeps the core dependency light while enabling `DataFrames(tel)` or
CSV/Arrow output when desired.

## Logging and diagnostics
Use `Logging.jl` for messages instead of print statements. Keep logging out of
hot loops to avoid overhead.

Guidelines:
- `@debug` for detailed internal state; `@info` for high-level progress; `@warn`
  for recoverable issues; `@error` for failures.
- Allow user-configured loggers and levels; do not set global logger by default.
- For progress bars, prefer `ProgressMeter.jl` in optional, user-facing code.

## Error handling
Prefer structured exceptions over print-and-return. Define a small hierarchy of
domain errors for common failures.

Sketch:
```julia
abstract type AdaptiveOpticsError <: Exception end
struct InvalidConfiguration <: AdaptiveOpticsError
    msg::String
end
struct DimensionMismatchError <: AdaptiveOpticsError
    msg::String
end

Base.showerror(io::IO, e::AdaptiveOpticsError) = print(io, e.msg)
```

Use `throw(InvalidConfiguration("..."))` for invalid inputs and `@assert` only
for internal invariants.

## Additional considerations
- API stability: define a small, stable public API early; keep internals flexible.
- Units and metadata: use explicit naming for meters/arcsec; optionally integrate `Unitful.jl`.
- Config and provenance: provide config structs, TOML export, and run metadata capture.
- Reference tolerances: decide numeric tolerances per component (PSF energy, slope RMS).
- Precision strategy: define expected `Float32` vs `Float64` usage and impacts.
- Extensibility: keep WFS, DM, and atmosphere models open via traits and minimal interfaces.
 - Deterministic mode: first-class seed control and fixed thread counts for reproducibility.
 - Workspace discipline: explicit preallocation to keep hot paths allocation-free.
 - Parallelism discipline: avoid nested parallelism and oversubscription by design.
 - GPU readiness: backend-parameterized arrays with minimal host-device transfers.
 - Diagnostics: structured logging and typed exceptions, not print statements.

## I/O strategy (no baked-in FITS)
Keep the core independent of any specific file format. Provide optional helpers
or extension modules for FITS/HDF5/CSV as needed.

Guidelines:
- Core types operate on in-memory arrays and structured configs.
- Optional helpers can live under `IO/` or an extension package.
- Prefer explicit `read_*` / `write_*` utilities that return plain arrays and
  metadata, so downstream code stays format-agnostic.

## Extension modules
Optional features (I/O formats, visualization, ModelingToolkit control) should
live in extension modules or separate packages to keep core dependencies light.
Use Julia's package extension mechanism or a thin `AdaptiveOpticsIO.jl` / `AdaptiveOpticsViz.jl`
package if preferred.

## GPU strategy
Provide GPU-backed workspaces and dispatch based on a trait, e.g.
`supports_gpu(::WFS)` and `array_backend(::Workspace)`. Keep core algorithms
backend-agnostic with `AbstractArray`.

Guidelines:
- Make array backends explicit (`backend=Array` or `backend=CuArray`).
- Avoid scalar indexing on GPU; call `CUDA.allowscalar(false)` in tests.
- Preallocate GPU workspaces; minimize host-device transfers.
- Use `AbstractFFTs.jl` so `FFTW` and `CUFFT` share the same API.
- Use `KernelAbstractions.jl` for kernels that should run on CPU and GPU.
- Keep noise generation backend-aware (CPU RNG vs GPU RNG).

## ModelingToolkit integration (control layer)
ModelingToolkit.jl is a good fit for the AO control loop (filters, delays, DM
dynamics), but not for the FFT-heavy optics pipeline. The intended split:
- Optics, WFS, phase screens: direct numerical kernels using multiple dispatch.
- Control and loop dynamics: optional ModelingToolkit-based models that ingest WFS slopes
  and output DM commands.

Example sketch:
```julia
using ModelingToolkit

@parameters t K τ s
@variables i(t) dm(t)
D = Differential(t)

# Simple integrator + first-order DM dynamics:
# i' = K * slopes
# dm' = (i - dm) / τ
eqs = [
    D(i)  ~ K * s,
    D(dm) ~ (i - dm) / τ
]
sys = ODESystem(eqs)
```

Hook this into the loop via a `Controller` interface:
```julia
abstract type AbstractController end
update!(ctrl::AbstractController, slopes, dt) = # returns DM command
```

One possible integration pattern:
```julia
using DifferentialEquations

params = (K=0.3, τ=0.02, s=0.0)
u0 = [0.0, 0.0]  # i(0), dm(0)
prob = ODEProblem(sys, u0, (0.0, 0.0), params)
integrator = init(prob, Tsit5(); save_everystep=false)

function update!(ctrl, slopes, dt)
    # Treat slopes as piecewise-constant input over dt
    ctrl.integrator.p[ctrl.idx_s] = slopes
    step!(ctrl.integrator, dt)
    return ctrl.integrator.u[2] # dm state
end
```

For vector-valued slopes, make `s` a vector parameter and store `idx_s` as a
range into the parameter vector. For performance, use `SVector` or
`ComponentArrays.jl` for structured parameters.

### Control module stub (sketch)
```julia
module Control

export AbstractController, ModelingToolkitController, DiscreteIntegratorController, update!

using DifferentialEquations

abstract type AbstractController end

mutable struct ModelingToolkitController <: AbstractController
    integrator
    idx_s::UnitRange{Int}
end

function update!(ctrl::ModelingToolkitController, slopes, dt)
    ctrl.integrator.p[ctrl.idx_s] .= slopes
    step!(ctrl.integrator, dt)
    return ctrl.integrator.u[2]
end

mutable struct DiscreteIntegratorController <: AbstractController
    gain::Float64
    τ::Float64
    i::Vector{Float64}
    dm::Vector{Float64}
end

function update!(ctrl::DiscreteIntegratorController, slopes, dt)
    ctrl.i .+= ctrl.gain .* slopes .* dt
    ctrl.dm .+= (ctrl.i .- ctrl.dm) .* (dt / ctrl.τ)
    return ctrl.dm
end

end
```

### Multi-channel ModelingToolkit sketch
```julia
using ModelingToolkit

n = 64
@parameters t K τ
@parameters s[1:n]
@variables i(t)[1:n] dm(t)[1:n]
D = Differential(t)

eqs = [
    D.(i)  .~ K .* s,
    D.(dm) .~ (i .- dm) ./ τ
]
sys = ODESystem(eqs)
```

The ModelingToolkit model can live in a `Control/` submodule and be swapped with a
lightweight numeric controller for performance-critical runs.

## Testing and validation
- Unit tests for optics primitives (FFT normalization, PSF energy, Zernike orthogonality).
- Integration tests for a minimal AO loop with fixed seeds.
- A tutorial validation suite that compares OOPAO reference outputs within tolerance.
- Store reference outputs from OOPAO (PSF, slopes, reconstructor outputs) with fixed seeds.
See `docs/oopao-reference-datasets.md` for a proposed reference dataset plan.

## Julia ecosystem building blocks
These packages cover most of the OOPAO dependency surface area:
- FFTs: `FFTW.jl`, `AbstractFFTs.jl`
- GPU: `CUDA.jl`, `KernelAbstractions.jl`, `Adapt.jl`
- Linear algebra: `LinearAlgebra` (stdlib), `SparseArrays`
- Interpolation/resampling: `Interpolations.jl`, `ImageTransformations.jl`
- Image processing: `ImageFiltering.jl`, `Images.jl`
- FITS IO: `FITSIO.jl` (optional)
- Performance: `StaticArrays.jl`, `StructArrays.jl`, `LoopVectorization.jl`
- Config/serialization: `TOML` (stdlib), `JSON3.jl`
- Visualization: `Plots.jl` as default; `Makie.jl` optional for advanced use

If you have a specific Julia AO package in mind, I can check whether it provides
phase screens, WFS models, or reconstructor utilities we can reuse.
