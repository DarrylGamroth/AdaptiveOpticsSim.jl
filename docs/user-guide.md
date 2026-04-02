# User Guide

Status: active

Start with:

- [`documentation-map.md`](./documentation-map.md) for doc navigation
- [`platform-architecture.md`](./platform-architecture.md) for package-level
  synthesis
- [`platform-workflows.md`](./platform-workflows.md) for platform-scale usage
  and evidence workflows
- [`scenario-builder-style.md`](./scenario-builder-style.md) for scenario
  authoring conventions
- [`api-reference.md`](./api-reference.md) for maintained public APIs
- [`model-validity-matrix.md`](./model-validity-matrix.md) for model evidence
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md) for performance
  surfaces

AdaptiveOpticsSim.jl is an idiomatic Julia adaptive-optics toolkit. The package
keeps the OOPAO feature set recognizable, but shifts the design toward
multiple dispatch, explicit state, and deterministic execution.

## Public API tiers

The package now distinguishes between:

- stable exported workflow APIs
- advanced but maintained APIs that may require qualification as
  `AdaptiveOpticsSim.<name>`
- developer/backend support APIs used mainly by benchmark and extension code

For ordinary usage, start with the exported workflow surface shown in the quick
start and the main API sections in
[api-reference.md](./api-reference.md).

If you want the package-level picture before diving into workflows, read
[platform-architecture.md](./platform-architecture.md) first.

If you want the maintained script-first workflow for realistic simulations,
validation, or benchmark runs, then read
[platform-workflows.md](./platform-workflows.md) next.

If you are writing or reviewing a realistic scenario script, then use
[scenario-builder-style.md](./scenario-builder-style.md) as the maintained
style guide.

For advanced utilities such as telemetry/config helpers, scenario builders,
build-backend policy helpers, and misregistration-identification tools, use
namespaced access. Examples:

```julia
ws = AdaptiveOpticsSim.Workspace(tel; rng=MersenneTwister(0))
sim = AdaptiveOpticsSim.initialize_ao_shwfs(...)
sprint = AdaptiveOpticsSim.SPRINT(tel, dm, wfs, basis)
cpu_backend = AdaptiveOpticsSim.CPUBuildBackend()
```

## Core model

- Parameters live in immutable structs such as `TelescopeParams` and
  `ShackHartmannParams`.
- Mutable simulation data lives in state objects such as `TelescopeState`,
  `DetectorState`, and `PyramidState`.
- Hot-path methods use the `!` convention: `compute_psf!`, `advance!`,
  `measure!`, `reconstruct!`, `apply!`.
- Algorithm choice is carried by dispatch and traits:
  `Geometric()` vs `Diffractive()`, `NoisePhoton()` vs `NoiseReadout()`,
  `DMAdditive()` vs `DMReplace()`.

## Quick start

### PSF generation

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
psf = compute_psf!(tel, src; zero_padding=2)
```

### Atmosphere and sensing

```julia
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    L0=25.0,
    fractional_cn2=(0.6, 0.4),
    wind_speed=(8.0, 12.0),
    wind_direction=(0.0, 90.0),
    altitude=(0.0, 5000.0),
)

wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6)
src = Source(band=:I, magnitude=8.0)

advance!(atm, tel)
propagate!(atm, tel)
slopes = measure!(wfs, tel, src)
```

### Closed loop

```julia
dm = DeformableMirror(tel; n_act=3, influence_width=0.35)
imat = interaction_matrix(dm, wfs, tel, src; amplitude=1e-9)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(dm.state.coefs)

reconstruct!(cmd, recon, slopes)
dm.state.coefs .= -cmd
apply!(dm, tel, DMAdditive())
```

## Recommended learning path

Use this progression instead of starting from subsystem plans:

1. image formation and PSF generation
2. atmosphere plus one WFS
3. detector-backed sensing
4. closed-loop runtime
5. field propagation, polychromatic sensing, and extended-source cases

The examples under `examples/tutorials/` already follow this shape well enough
to be the primary learning surface.

## Common workflows

### Workflow 1: optics-only PSF studies

Start with:

- `Telescope`
- `Source`
- `compute_psf!`
- optional `NCPA`, `OPDMap`, `apply_opd!`

Use this when you care about:

- pupil construction
- PSF normalization
- image formation
- simple aberration studies

### Workflow 2: atmosphere and WFS studies

Start with:

- `MultiLayerAtmosphere` or `InfiniteMultiLayerAtmosphere`
- `advance!`
- `propagate!`
- one of `ShackHartmann`, `PyramidWFS`, `BioEdgeWFS`, `CurvatureWFS`,
  `ZernikeWFS`
- `measure!`

Use this when you care about:

- sensor behavior
- atmosphere evolution
- detector-coupled readout
- optical gain / calibration studies

### Workflow 3: closed-loop simulation

Start with:

- `DeformableMirror`
- calibration helpers such as `interaction_matrix`
- a reconstructor
- `ClosedLoopRuntime`
- `prepare!`
- `update!`, `simulation_readout`, `simulation_command`, `simulation_slopes`

Use this when you care about:

- loop staging
- latency
- exported runtime products
- HIL-style or realistic detector-backed sensing paths

### Workflow 4: field propagation and diffractive optics

Start with:

- `ElectricField`
- `FraunhoferPropagation`, `FresnelPropagation`
- `AtmosphericFieldPropagation`
- `SpectralSource`
- `ExtendedSource`

Use this when you care about:

- field-level propagation
- polychromatic sensing
- extended sources
- curvature or atmosphere-aware field propagation

## Determinism

- Use a fixed `MersenneTwister` and pass it into `advance!` or detector calls.
- Keep detector noise disabled when comparing against reference datasets unless
  the test is explicitly about noise.
- Run single-threaded when strict reproducibility matters. See
  [`deterministic-simulation.md`](./deterministic-simulation.md).

## Logging and errors

- Use `Logging.jl` macros such as `@info` and `@warn` in user-facing code.
- The core throws structured exceptions such as `InvalidConfiguration` and
  `DimensionMismatchError`.
- Avoid logging inside inner loops; return state and let the caller decide how
  much telemetry to emit.

## Tutorials and examples

Runnable example ports live under `examples/tutorials/`. Start with:

- `examples/tutorials/image_formation.jl`
- `examples/tutorials/detector.jl`
- `examples/tutorials/transfer_function.jl`
- `examples/tutorials/gain_sensing_camera.jl`
- `examples/tutorials/closed_loop_shack_hartmann.jl`
- `examples/tutorials/closed_loop_pyramid.jl`
- `examples/tutorials/closed_loop_bioedge.jl`

See `docs/julia-tutorial-mappings.md` for the full mapping back to OOPAO.

## Validation workflow

- `Pkg.test()` runs the tutorial smoke tests plus the reference harness.
- `test/reference_data/` stores a small deterministic OOPAO bundle for
  cross-validation.
- `ENV["ADAPTIVEOPTICS_REFERENCE_ROOT"]` can point to an alternate bundle.
- The committed bundle currently covers PSF, geometric/diffractive
  Shack-Hartmann, Pyramid, BioEdge, GSC optical gains, and the analytical
  transfer-function workflow. LiFT, closed-loop traces, and tomography remain
  follow-on parity work.

See [`oopao-reference-datasets.md`](./oopao-reference-datasets.md) for the
bundle contract.

For the maintained synthesis of model evidence, use
[`model-validity-matrix.md`](./model-validity-matrix.md).

## Benchmark and comparison workflow

Use the benchmark docs intentionally:

- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - maintained local runtime surfaces
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
  - archived cross-package OOPAO/SPECULA/REVOLT-like evidence

These are engineering evidence surfaces, not normal tutorial entry points and
not part of `Pkg.test()`.

## Plotting

The core package does not require plotting. For user-facing visualization, use
`Plots.jl` in examples, notebooks, or downstream tooling rather than wiring it
into hot paths.

## Backends and performance

Current code already parameterizes many core arrays by numeric type and array
backend. The intended execution model is:

- `AbstractFFTs.jl` for FFT portability across CPU and GPU backends.
- trait-driven algorithm selection for geometric vs diffractive sensing and for
  backend-aware kernels.
- `KernelAbstractions.jl` for non-FFT kernels that should run on CPU and GPU
  from a single implementation.

The package is not fully on `KernelAbstractions.jl` yet, but the design target
is explicit: stencil-like loops, binning, focal-plane mask application,
weighting updates, and per-subaperture reductions should move behind
trait-selected kernels instead of hard-coded CPU loops.
