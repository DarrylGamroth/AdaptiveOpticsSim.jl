# User Guide

Status: active

This is the main user-facing guide.

If you only need the normal package entry points, read these and stop there:

- [../README.md](../README.md)
- [model-cookbook.md](model-cookbook.md)
- [api-reference.md](api-reference.md)
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
- sensing objects such as `ShackHartmannWFS`, `PyramidWFS`, and `BioEdgeWFS`
- `Detector` when the sensing path needs explicit detector physics
- controllable optics such as `DeformableMirror`, `ModalControllableOptic`, `TipTiltMirror`, `FocusStage`, and `CompositeControllableOptic`
- modal-optic basis specs such as `CartesianTiltBasis`, `ZernikeOpticBasis`, and `MatrixModalBasis` when you want to choose controlled modes explicitly
- `ControlLoopScenario` when you want the maintained step-wise AO or HIL simulation surface

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

- `ControlLoopBranch`
- `SingleControlLoopConfig`
- `GroupedControlLoopConfig`
- `ControlLoopScenario`
- `build_control_loop_scenario(...)`

This layer is the recommended path for:

- normal closed-loop examples
- HIL / RTC-boundary simulations
- grouped or multi-branch runtime composition

The lower-level `ClosedLoopRuntime` + `AdaptiveOpticsSim.simulation_interface(...)` path remains
available, but it is an advanced single-runtime surface rather than the default
user entry point.

For a compact recipe-first version of this guide, use [model-cookbook.md](model-cookbook.md).

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

wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6)

advance!(atm, tel)
propagate!(atm, tel)
slopes = measure!(wfs, tel, src)
```

For HIL or RTC export, attach a detector and request the detector image after
measurement. `bits` defines the quantization depth, `full_well` defines the
analog-to-digital scaling, and `output_type` defines the Julia array element
type used for the exported frame:

```julia
det = Detector(noise=NoiseNone(), full_well=30_000.0, bits=12, output_type=UInt16)
measure!(wfs, tel, src, det; rng=rng)
adu = wfs_detector_image(wfs, det)
```

Here `adu` is a `UInt16` image containing 12-bit ADU values. For
Shack-Hartmann sensors this is the lenslet spot mosaic; for frame-style WFSs it
is the maintained detector/readout frame.

For a quantitative CMOS camera similar to the Hamamatsu ORCA-Quest family, use
the qCMOS convenience constructors. The camera preset supplies the family-level
readout noise, QE, full well, dark current, frame rate, and pixel metadata;
ORCA-Quest presets default to rolling-shutter timing. `bits` and `output_type`
still define the exported RTC/HIL frame:

```julia
det = ORCAQuest2Detector(
    scan_mode=QCMOSUltraQuietScan(),
    integration_time=1e-3,
    bits=16,
    output_type=UInt16,
)

measure!(wfs, tel, src, det; rng=rng)
adu = wfs_detector_image(wfs, det)
```

The qCMOS model intentionally keeps camera-specific defaults separate from
measured calibration data. Add `PixelResponseNonuniformity`,
`DarkSignalNonuniformity`, or `BadPixelMask` through the detector defect model
when you have a calibrated ORCA-Quest unit.

Rolling-shutter detectors can also capture a time-varying scene. Use
`InPlaceFrameSource` when the source can write into a preallocated frame, or
`FunctionFrameSource` when a function returns a frame for each sample time:

```julia
det = Detector(
    noise=NoiseNone(),
    sensor=CMOSSensor(timing_model=RollingShutter(25e-6)),
    response_model=NullFrameResponse(),
)

pulse = InPlaceFrameSource((out, t) -> begin
    fill!(out, 0.0)
    t >= 50e-6 && fill!(out, 1.0)
    return out
end, (64, 64))

frame = capture!(det, pulse; rng=rng)
```

This path samples each rolling-shutter row group at its own readout time, so it
can show transient illumination and rolling-shutter artifacts. Static
`capture!(det, image)` remains the preferred path when the scene does not vary
during the exposure.

Use this when you care about:

- sensor behavior
- atmosphere evolution
- detector-coupled readout
- optical gain and calibration studies

### Workflow 3: Closed-loop AO simulation

```julia
using Random

rng = runtime_rng(0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = ControlLoopBranch(:main, sim, recon; rng=rng)

cfg = SingleControlLoopConfig(
    name=:closed_loop_demo,
    branch_label=:main,
    outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true),
)

scenario = build_control_loop_scenario(cfg, branch)
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

It also accepts explicit advanced composition through:

- `topology=...`
- `influence_model=...`
- `actuator_model=...`

For example:

- `GaussianInfluenceWidth(0.3)`
- `GaussianMechanicalCoupling(0.08)`
- `DenseInfluenceMatrix(modes)`
- `MeasuredInfluenceFunctions(modes; metadata=...)`
- `ActuatorGridTopology(16)`
- `SampledActuatorTopology(coords; valid_actuators=mask, metadata=...)`
- `ClippedActuators(-0.2, 0.2)`
- `ActuatorHealthMap(gains)`
- `CompositeDMActuatorModel(...)`

Use the scalar keywords for normal work. Use `influence_model=...` when you
need an explicit DM influence representation, `topology=...` when the actuator
layout is not the default full square grid, and `actuator_model=...` when
command preprocessing should model clipping or actuator health without changing
the sampled influence basis.

Use this when you care about:

- loop staging
- latency
- exported runtime outputs
- HIL-style or detector-backed sensing paths

`prepare!(...)` performs any wfs/runtime precomputation once before repeated
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

- `ShackHartmannWFS`
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
readout behavior, windowing, or exported frame outputs.

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

Use [api-reference.md](api-reference.md) for the exported surface. If a name is
not exported, prefer qualifying it rather than adding it to the public namespace
unless it is part of an ordinary workflow or documented extension seam.

For ordinary usage, start with the exported workflow surface shown above and in
[api-reference.md](api-reference.md).

For advanced utilities such as telemetry/config helpers and some backend policy
helpers, use namespaced access. Examples:

```julia
ws = AdaptiveOpticsSim.Workspace(tel; rng=deterministic_reference_rng(0))
sprint = AdaptiveOpticsSim.SPRINT(tel, dm, wfs, basis)
```

## Determinism

- Use a fixed RNG and pass it into `advance!`, detector calls, or runtime
  constructors instead of relying on `Random.default_rng()`.
- Use `deterministic_reference_rng(seed)` for reference datasets and regression
  fixtures. This preserves the long-standing `MersenneTwister` stream.
- Use `runtime_rng(seed)` for new RTC/HIL-style simulations and benchmarks.
  This uses `Xoshiro`, which is a better default for throughput-oriented hot
  paths while still being repeatable for a fixed software stack.
- Keep detector noise disabled when comparing against reference datasets unless
  the test is explicitly about noise.
- Run single-threaded when strict reproducibility matters.

See [deterministic-simulation.md](deterministic-simulation.md).

## Tutorials and Examples

Runnable example ports live under `examples/tutorials/`. Good starting points:

- `examples/tutorials/image_formation.jl`
- `examples/tutorials/detector.jl`
- `examples/tutorials/closed_loop_shack_hartmann.jl`
- `examples/tutorials/closed_loop_pyramid.jl`
- `examples/tutorials/closed_loop_bioedge.jl`
- `examples/tutorials/closed_loop_zernike.jl`

See [julia-tutorial-mappings.md](julia-tutorial-mappings.md) for the mapping
back to OOPAO tutorials.

## When You Need More Than The User Guide

Only go deeper if your task actually needs it:

- public API details:
  - [api-reference.md](api-reference.md)
- maintained validation status:
  - [model-validity-matrix.md](model-validity-matrix.md)
- supported production scope:
  - [supported-production-surfaces.md](supported-production-surfaces.md)
- benchmark and cross-package evidence:
  - benchmark artifacts under `benchmarks/results/`
  - benchmark artifacts under `benchmarks/results/`
- maintainer/developer navigation:
  - [documentation-map.md](documentation-map.md)
