# Model Cookbook

Status: active

This cookbook is the shortest path from "I want to model X" to a working
AdaptiveOpticsSim.jl script.

Use this together with:

- [../README.md](../README.md)
- [user-guide.md](./user-guide.md)
- [api-reference.md](./api-reference.md)

Each recipe is intentionally small. Start here, then grow the script toward your
real instrument or benchmark surface.

## Recipe 1: PSF From A Telescope And Source

Use this when you want a minimal optics-only model.

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
psf = compute_psf!(tel, src; zero_padding=2)
```

Key objects:

- `Telescope`
- `Source`
- `compute_psf!`

## Recipe 2: Atmosphere Plus Shack-Hartmann Sensing

Use this when you want a single-sensor atmospheric simulation.

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
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

Swap `ShackHartmann(...)` for `PyramidWFS(...)`, `BioEdgeWFS(...)`,
`CurvatureWFS(...)`, or `ZernikeWFS(...)` when the sensing family changes.

## Recipe 3: Detector-Backed Sensing

Use this when the sensor path should include detector physics and exported
frames.

```julia
using AdaptiveOpticsSim
using Random

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6)
det = Detector(
    noise=NoiseReadout(1.0),
    integration_time=1.0,
    qe=1.0,
    binning=1,
)

advance!(atm, tel)
propagate!(atm, tel)
measure!(wfs, tel, src; detector=det, rng=MersenneTwister(0))
frame = output_frame(det)
```

Use this pattern when you care about:

- readout behavior
- exported pixels
- detector-backed WFS comparisons

## Recipe 4: Closed-Loop AO Model

Use this when you want the maintained public AO runtime surface rather than
isolated subsystem calls.

```julia
using AdaptiveOpticsSim
using Random

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive())
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = RuntimeBranch(:main, sim, recon; rng=MersenneTwister(0))

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
```

Use this pattern when you need:

- loop stepping
- runtime latency products
- command/slopes/frame export from one maintained runtime boundary

## Recipe 5: HIL / RTC Boundary Simulation

Use this when the simulation should behave like a hardware-in-the-loop or RTC
boundary surface: commands in, WFS/science products out.

For a single-branch RTC-style runtime with an external command source:

```julia
using AdaptiveOpticsSim
using Random

# Build the physical plant.
tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive())
sim = AOSimulation(tel, src, atm, dm, wfs)

# Use NullReconstructor() when the controller lives outside the package and
# commands are injected explicitly through set_command!(...).
recon = NullReconstructor()

# Attach detectors when the external interface should export pixel products.
wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
science_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

branch = RuntimeBranch(
    :main,
    sim,
    recon;
    wfs_detector=wfs_det,
    science_detector=science_det,
    rng=MersenneTwister(1),
)

cfg = SingleRuntimeConfig(
    name=:single_runtime_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)

# `prepare!` performs one-time runtime/WFS precomputation and enables the
# requested export surfaces before repeated `sense!` or `step!` calls.
prepare!(scenario)

# Inspect the external command contract before wiring in a controller.
# For this single-DM plant the boundary is one flat command vector.
n_command = length(command(scenario))
external_command = zeros(eltype(command(scenario)), n_command)

# Typical HIL / plant-facing flow per control step:
#   1. get a command from the external controller
#   2. inject it into the simulated plant
#   3. run only the sensing side of the package
#   4. read the exported products that go back to the controller
for k in 1:3
    fill!(external_command, 0)
    external_command[1] = 0.05 * k

    set_command!(scenario, external_command)
    sense!(scenario)

    command_vec = command(scenario)
    slopes_vec = slopes(scenario)
    wfs_img = wfs_frame(scenario)
    science_img = science_frame(scenario)
    wfs_meta = wfs_metadata(scenario)
    science_meta = science_metadata(scenario)

    @show k length(command_vec) length(slopes_vec)
    @show size(wfs_img) size(science_img)
    @show wfs_meta.output_size science_meta.output_size
end

# If the plant itself contains several controllable surfaces, make them explicit
# in the optic model and still drive them through one packed RTC boundary.
tiptilt = TipTiltMirror(tel; scale=0.1, label=:tiptilt)
high_order_dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
composite_optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => high_order_dm)
composite_sim = AOSimulation(tel, src, atm, composite_optic, wfs)

composite_branch = RuntimeBranch(
    :main,
    composite_sim,
    recon;
    wfs_detector=wfs_det,
    science_detector=science_det,
    rng=MersenneTwister(2),
)
composite_scenario = build_runtime_scenario(cfg, composite_branch)
prepare!(composite_scenario)
set_command!(composite_scenario, (
    tiptilt=fill(0.01, 2),
    dm=fill(0.02, 16),
))
sense!(composite_scenario)

# For a configurable modal surface, prefer an explicit basis spec. This example
# groups piston, tip, tilt, and focus-like Zernike commands into one segment.
pttf = ModalControllableOptic(
    tel,
    ZernikeOpticBasis([1, 2, 3, 4]; scale=0.02);
    labels=:pttf,
)

# A steering-style two-axis low-order optic can be spelled directly too.
steering = ModalControllableOptic(tel, CartesianTiltBasis(scale=0.1); labels=:steering)

# Partial updates do not require manual slice math either. Here only the
# tip/tilt segment changes; the DM segment remains staged in-place.
update_command!(composite_scenario, (
    tiptilt=fill(0.015, 2),
))
sense!(composite_scenario)

# If a runtime boundary only needs a logical split over one existing optic, you
# can still provide an explicit command layout without changing the plant.
split_branch = RuntimeBranch(
    :main,
    sim,
    recon;
    science_detector=science_det,
    rng=MersenneTwister(3),
    command_layout=RuntimeCommandLayout(:woofer => 8, :tweeter => 8),
)
split_scenario = build_runtime_scenario(cfg, split_branch)
prepare!(split_scenario)
set_command!(split_scenario, (
    woofer=fill(0.01, 8),
    tweeter=fill(0.02, 8),
))
sense!(split_scenario)

# Export semantics for the single-branch `scenario` above:
# - `command(scenario)` is one packed command vector
# - `slopes(scenario)` is one packed slopes vector
# - `wfs_frame(scenario)` is one WFS frame, or `nothing`
# - `science_frame(scenario)` is one science frame, or `nothing`
# - detector metadata describes the pixel contract seen by the controller
```

Use `step!(scenario)` only when you want the package to perform its own
reconstruction and command update. For an external RTC or HIL controller,
prefer the explicit `set_command!(scenario, cmd)` plus `sense!(scenario)` flow
above.

The exported science frame is the configured science detector's `output_frame`,
so its shape and element type are detector-defined. Inspect
`science_metadata(scenario)` to learn the exact export contract for:

- `output_size`
- `output_precision`
- `binning`
- `window_rows` / `window_cols`
- `bits` / `full_well`

For grouped RTC boundaries, keep the outer structure branch-oriented and nest
structured commands per branch when needed:

```julia
set_command!(grouped_scenario, (
    high_order=(; woofer=woofer_cmd, tweeter=tweeter_cmd),
    low_order=(; steering=tt_cmd, dm=ho_cmd),
))
```

Grouped output semantics differ from the single-branch case:

- `command(grouped_scenario)` and `slopes(grouped_scenario)` are still one
  packed vector each, aggregated across branches.
- `wfs_frame(grouped_scenario)` and `science_frame(grouped_scenario)` become
  per-branch frame collections when those grouped frame products are enabled.
- `grouped_wfs_stack(grouped_scenario)` and
  `grouped_science_stack(grouped_scenario)` provide stacked array forms when
  the grouped contract requests them.

That lets an external HIL client validate the runtime boundary before wiring the
frame into another controller, transport layer, or logging path.

For grouped or multi-branch RTC-style composition, start from:

- [platform_grouped_runtime.jl](../examples/closed_loop/platform_grouped_runtime.jl)
- [platform_single_runtime.jl](../examples/closed_loop/platform_single_runtime.jl)

Use this pattern when you need:

- a clear input/output simulation boundary
- exported WFS or science frames per step
- grouped runtime products
- HIL-style benchmarking or controller integration

## Recipe 6: GPU HIL Runtime

Use this when the plant should run on CUDA or AMDGPU and the external
controller still talks to the same HIL runtime boundary.

```julia
using AdaptiveOpticsSim
using CUDA
using Random

const GPU = CUDABackend()

# Choose backend and precision once near the top of the script.
tel = Telescope(
    resolution=16,
    diameter=8.0f0,
    sampling_time=1f-3,
    central_obstruction=0.0f0,
    T=Float32,
    backend=GPU,
)
src = Source(band=:I, magnitude=0.0f0, T=Float32)
atm = KolmogorovAtmosphere(tel; r0=0.2f0, L0=25.0f0, T=Float32, backend=GPU)

# Multi-surface plants work on GPU too.
tiptilt = TipTiltMirror(tel; scale=0.1f0, label=:tiptilt, T=Float32, backend=GPU)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3f0, T=Float32, backend=GPU)
optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)

wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=Float32, backend=GPU)
wfs_det = Detector(noise=NoiseNone(), integration_time=1.0f0, qe=1.0f0, binning=1; T=Float32, backend=GPU)
science_det = Detector(noise=NoiseNone(), integration_time=1.0f0, qe=1.0f0, binning=1; T=Float32, backend=GPU)

sim = AOSimulation(tel, src, atm, optic, wfs)
branch = RuntimeBranch(
    :main,
    sim,
    NullReconstructor();
    wfs_detector=wfs_det,
    science_detector=science_det,
    rng=MersenneTwister(1),
)
cfg = SingleRuntimeConfig(
    name=:hil_gpu,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)

# `prepare!` performs one-time runtime/WFS precomputation and enables the
# requested export surfaces before repeated `sense!` or `step!` calls.
prepare!(scenario)

# The external controller can still send structured commands by segment.
set_command!(scenario, (
    tiptilt=Float32[0.01, -0.01],
    dm=fill(0.02f0, 16),
))
sense!(scenario)

# These products stay backend-native until you explicitly copy them to CPU.
slopes_vec = slopes(scenario)
wfs_img = wfs_frame(scenario)
science_img = science_frame(scenario)

@show typeof(slopes_vec)
@show typeof(wfs_img)
@show typeof(science_img)

# Only copy to CPU when another process or transport layer really needs host
# memory.
wfs_host = Array(wfs_img)
science_host = Array(science_img)
```

Notes:

- Replace `CUDABackend()` with `AMDGPUBackend()` for the AMDGPU path.
- Keep `backend=GPU` on the long-lived plant objects so runtime buffers stay on
  device.
- If you need a GPU-built internal reconstructor as well, build that
  calibration surface intentionally with `AdaptiveOpticsSim.GPUArrayBuildBackend(...)`.
- `set_command!` accepts ordinary CPU vectors or tuples; the runtime stages the
  command into the backend-native optic state before `sense!(...)`.
- `CPUBackend()` is the default public constructor surface. Use `CUDABackend()` or
  `AMDGPUBackend()` to move the plant onto a GPU.

Use this pattern when you need:

- external-control / HIL semantics on GPU
- backend-native WFS and science exports
- realistic CUDA or AMDGPU runtime profiling
- plant models with several controllable surfaces on one RTC boundary

## Recipe 7: Advanced Direct Runtime Interface

Use this only when you are manually assembling or testing a single runtime and
do not need the full scenario/config layer.

```julia
runtime = ClosedLoopRuntime(sim, recon; rng=MersenneTwister(0))
interface = AdaptiveOpticsSim.simulation_interface(runtime)
prepare!(interface)
step!(interface)
rt = readout(interface)
```

This is a maintained advanced surface, but it is not the default public runtime
assembly path. Prefer `build_runtime_scenario(...)` for normal closed-loop and
HIL work.

## How To Choose The Right Entry Surface

Use:

- subsystem functions such as `compute_psf!`, `advance!`, `propagate!`, and
  `measure!`
  - when you are studying one physical layer
- `SingleRuntimeConfig` or `GroupedRuntimeConfig`
  - when you want the maintained public runtime/orchestration surface
- `ClosedLoopRuntime` plus `AdaptiveOpticsSim.simulation_interface(...)`
  - only when you are manually assembling or testing one low-level runtime

## Next Step

After the first working script:

1. move constants into a small `build_scenario(...)` helper
2. make RNG setup explicit when reproducibility matters
3. add detector realism only when the study needs it
4. move to [user-guide.md](./user-guide.md) for broader workflow guidance
5. move to [api-reference.md](./api-reference.md) when you need the full public
   surface
