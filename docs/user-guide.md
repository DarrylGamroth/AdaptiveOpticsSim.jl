# User Guide

Status: active

This is the main user-facing guide.

If you only need the normal package entry points, read these and stop there:

- [../README.md](../README.md)
- [model-cookbook.md](./model-cookbook.md)
- [api-reference.md](./api-reference.md)
- `examples/tutorials/`

You do not need the platform, benchmark, audit, or production-hardening docs to
build a normal AO model.

Use those deeper docs only when you are:

- validating backend parity
- working on production support claims
- benchmarking against OOPAO/SPECULA/REVOLT-like surfaces
- refactoring or maintaining the package internals

## What To Read Next

Choose one path:

- optics and PSF work:
  - `examples/tutorials/image_formation.jl`
  - `examples/tutorials/ncpa.jl`
- detector and WFS work:
  - `examples/tutorials/detector.jl`
  - `examples/tutorials/shack_hartmann_subapertures.jl`
  - `examples/tutorials/extended_source_sensing.jl`
- closed-loop AO work:
  - `examples/closed_loop_demo.jl`
  - `examples/tutorials/closed_loop_shack_hartmann.jl`
  - `examples/tutorials/closed_loop_pyramid.jl`
- calibration and identification work:
  - `examples/tutorials/gain_sensing_camera.jl`
  - `examples/tutorials/lift.jl`
  - `examples/tutorials/transfer_function.jl`

## Mental Model

The package is organized around a small set of modeling objects:

- `Telescope` and `Source`
- atmosphere objects such as `KolmogorovAtmosphere` and `MultiLayerAtmosphere`
- sensing objects such as `ShackHartmann`, `PyramidWFS`, and `BioEdgeWFS`
- `Detector` when the sensing path needs explicit detector physics
- controllable optics such as `DeformableMirror`, `ModalControllableOptic`, `TipTiltMirror`, `FocusStage`, and `CompositeControllableOptic`
- modal-optic basis specs such as `CartesianTiltBasis`, `ZernikeOpticBasis`, and `MatrixModalBasis` when you want to choose controlled modes explicitly
- `RuntimeScenario` when you want the maintained step-wise AO or HIL simulation surface

## Three Execution Layers

The package exposes three layers on purpose.

### 1. Primitive physics layer

Use this layer when you are studying or composing one physical subsystem at a
time.

Canonical verbs:

- `advance!`
- `propagate!`
- `apply!`
- `measure!`
- `capture!`

Typical use:

- PSF formation
- direct WFS studies
- calibration internals
- custom research scripts

### 2. Runtime execution layer

Use this layer when you want repeated AO execution with a stable command/input
and readout/output boundary.

Canonical verbs and accessors:

- `prepare!`
- `sense!`
- `step!`
- `set_command!`
- `update_command!`
- `readout`
- `command`, `slopes`, `wfs_frame`, `science_frame`

Semantics:

- `prepare!(...)` performs one-time runtime/WFS precomputation
- `sense!(...)` runs only the plant/sensor side
- `step!(...)` runs the full closed-loop update

### 3. Orchestration layer

Use this as the default public runtime assembly surface.

Canonical types and builder:

- `RuntimeBranch`
- `SingleRuntimeConfig`
- `GroupedRuntimeConfig`
- `RuntimeScenario`
- `build_runtime_scenario(...)`

This layer is the recommended path for:

- normal closed-loop examples
- HIL / RTC-boundary simulations
- grouped or multi-branch runtime composition

The lower-level `ClosedLoopRuntime` + `AdaptiveOpticsSim.simulation_interface(...)` path remains
available, but it is an advanced single-runtime surface rather than the default
user entry point.

For a compact recipe-first version of this guide, use [model-cookbook.md](./model-cookbook.md).

## Build A Model

### Workflow 1: Optics-only PSF

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
psf = compute_psf!(tel, src; zero_padding=2)
```

Use this when you care about:

- pupil construction
- PSF normalization
- image formation
- simple aberration studies

### Workflow 2: Atmosphere plus one WFS

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

advance!(atm, tel)
propagate!(atm, tel)
slopes = measure!(wfs, tel, src)
```

Use this when you care about:

- sensor behavior
- atmosphere evolution
- detector-coupled readout
- optical gain and calibration studies

### Workflow 3: Closed-loop AO simulation

```julia
using Random

rng = MersenneTwister(0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = RuntimeBranch(:main, sim, recon; rng=rng)

cfg = SingleRuntimeConfig(
    name=:closed_loop_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)

for _ in 1:5
    step!(scenario)
end

rt = readout(scenario)
cmd = command(rt)
slopes_vec = slopes(rt)
frame = wfs_frame(rt)
```

`DeformableMirror` keeps the concise public Gaussian inputs:

- `influence_width=...`
- `mechanical_coupling=...`

It also now accepts an explicit advanced model object through
`influence_model=...`, for example:

- `GaussianInfluenceWidth(0.3)`
- `GaussianMechanicalCoupling(0.08)`
- `DenseInfluenceMatrix(modes)`

Use the scalar keywords for normal work. Use `influence_model=...` when you
need an explicit DM influence representation.

Use this when you care about:

- loop staging
- latency
- exported runtime products
- HIL-style or detector-backed sensing paths

`prepare!(...)` performs any WFS/runtime precomputation once before repeated
`step!(...)` or `sense!(...)` calls. `step!(...)` runs the full closed-loop
update, while `sense!(...)` runs only the plant/sensor side and is the right
entry point when commands come from an external controller.

### Workflow 4: Field propagation and diffractive optics

Start with:

- `ElectricField`
- `FraunhoferPropagation`, `FresnelPropagation`
- `AtmosphericFieldPropagation`
- `ExtendedSource`

Use this when you care about:

- field-level propagation
- polychromatic sensing
- extended sources
- curvature or atmosphere-aware field propagation

## Choosing Components

### Atmosphere

- `KolmogorovAtmosphere`
  - compact single-screen studies
- `MultiLayerAtmosphere`
  - standard finite multilayer turbulence
- `InfiniteMultiLayerAtmosphere`
  - longer-running translated-screen studies

### Wavefront sensor

- `ShackHartmann`
  - general SH studies and HIL-style RTC surfaces
- `PyramidWFS`
  - pyramid sensing and modulation studies
- `BioEdgeWFS`
  - BioEdge variants
- `CurvatureWFS`
  - curvature sensing
- `ZernikeWFS`
  - Zernike WFS studies

### Detector

Use an explicit `Detector(...)` when the sensing path needs detector physics,
readout behavior, windowing, or exported frame products.

For counting-imager or counting-channel paths, use a maintained counting
detector family instead of the generic frame-detector surface:

- `APDDetector(...)` for channel-style counting readout
- `SPADArrayDetector(...)` for accumulated-count imaging arrays

Leave detector effects simple or disabled when the goal is deterministic model
comparison rather than detector realism.

## Public API Tiers

The package distinguishes between:

- stable exported workflow APIs
- advanced but maintained APIs that may require qualification as `AdaptiveOpticsSim.<name>`
- developer/backend support APIs used mainly by benchmark and extension code

For ordinary usage, start with the exported workflow surface shown above and in
[api-reference.md](./api-reference.md).

For advanced utilities such as telemetry/config helpers, scenario builders, and
some backend policy helpers, use namespaced access. Examples:

```julia
ws = AdaptiveOpticsSim.Workspace(tel; rng=MersenneTwister(0))
sim = AdaptiveOpticsSim.initialize_ao_shwfs(...)
sprint = AdaptiveOpticsSim.SPRINT(tel, dm, wfs, basis)
```

## Determinism

- Use a fixed `MersenneTwister` and pass it into `advance!` or detector calls.
- Keep detector noise disabled when comparing against reference datasets unless
  the test is explicitly about noise.
- Run single-threaded when strict reproducibility matters.

See [deterministic-simulation.md](./deterministic-simulation.md).

## Tutorials and Examples

Runnable example ports live under `examples/tutorials/`. Good starting points:

- `examples/tutorials/image_formation.jl`
- `examples/tutorials/detector.jl`
- `examples/tutorials/closed_loop_shack_hartmann.jl`
- `examples/tutorials/closed_loop_pyramid.jl`
- `examples/tutorials/closed_loop_bioedge.jl`
- `examples/tutorials/closed_loop_zernike.jl`

See [julia-tutorial-mappings.md](./julia-tutorial-mappings.md) for the mapping
back to OOPAO tutorials.

## When You Need More Than The User Guide

Only go deeper if your task actually needs it:

- public API details:
  - [api-reference.md](./api-reference.md)
- maintained validation status:
  - [model-validity-matrix.md](./model-validity-matrix.md)
- supported production scope:
  - [supported-production-surfaces.md](./supported-production-surfaces.md)
- benchmark and cross-package evidence:
  - [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)
  - [benchmark-matrix-plan.md](./benchmark-matrix-plan.md)
- maintainer/developer navigation:
  - [documentation-map.md](./documentation-map.md)
