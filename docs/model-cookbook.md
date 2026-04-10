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

Use this when you want an actual AO control step surface rather than isolated
subsystem calls.

```julia
using AdaptiveOpticsSim
using Random

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive())
sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
runtime = ClosedLoopRuntime(sim, recon; rng=MersenneTwister(0))
interface = simulation_interface(runtime)

for _ in 1:5
    step!(interface)
end

readout = simulation_readout(interface)
cmd = simulation_command(readout)
slopes = simulation_slopes(readout)
```

Use this pattern when you need:

- loop stepping
- runtime latency products
- command/slopes/frame export from one simulation object

## Recipe 5: HIL / RTC Boundary Simulation

Use this when the simulation should behave like a hardware-in-the-loop or RTC
boundary surface: commands in, WFS/science products out.

For a single-branch RTC-style runtime with an external DM command source:

```julia
using AdaptiveOpticsSim
using Random

# Build the physical plant.
tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive())
sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)

# Use NullReconstructor() when the control law lives outside the package and
# commands are injected through set_command!(...).
recon = NullReconstructor()

# Attach a science detector if the external interface should export a science
# image per sensing step.
science_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

branch = RuntimeBranch(
    :main,
    sim,
    recon;
    science_detector=science_det,
    rng=MersenneTwister(1),
)

cfg = SingleRuntimeConfig(
    name=:single_runtime_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)

# The external command vector must match the DM command length exported by the
# runtime boundary.
n_command = length(command(scenario))
external_command = zeros(eltype(command(scenario)), n_command)

# Example: set one actuator command from an external controller.
external_command[1] = 0.05

# HIL / plant-facing flow:
#   1. inject the external DM command
#   2. run the sensing side of the simulation
#   3. read the exported WFS/science products
set_command!(scenario, external_command)
sense!(scenario)

# If a runtime boundary represents several logical control surfaces packed into
# one RTC vector, provide an explicit command layout and then inject a
# structured NamedTuple instead of doing manual slice math. This is useful for
# tip/tilt mirrors, steering mirrors, or packed woofer/tweeter commands.
split_branch = RuntimeBranch(
    :main,
    sim,
    recon;
    science_detector=science_det,
    rng=MersenneTwister(2),
    command_layout=RuntimeCommandLayout(:woofer => 8, :tweeter => 8),
)
split_scenario = build_runtime_scenario(cfg, split_branch)
prepare!(split_scenario)
set_command!(split_scenario, (
    woofer=fill(0.01, 8),
    tweeter=fill(0.02, 8),
))
sense!(split_scenario)

# Exported boundary products after the sensing step.
command_vec = command(scenario)
slopes_vec = slopes(scenario)
wfs_img = wfs_frame(scenario)
science_img = science_frame(scenario)
wfs_meta = wfs_metadata(scenario)
science_meta = science_metadata(scenario)

# Typical controller-side checks:
@show length(command_vec)
@show length(slopes_vec)
@show size(wfs_img)
@show size(science_img)
@show science_meta.output_size
@show science_meta.output_precision
```

Use `step!(scenario)` only when you want the package to perform its own
reconstruction and DM update. For an external RTC or HIL controller, prefer the
explicit `set_command!(scenario, cmd)` plus `sense!(scenario)` flow above.

The exported science frame is the configured science detector's `output_frame`,
so its shape and element type are detector-defined. Inspect
`science_metadata(scenario)` to learn the exact export contract for:

For grouped RTC boundaries, keep the outer structure branch-oriented and nest
structured commands per branch when needed:

```julia
set_command!(grouped_scenario, (
    high_order=(; woofer=woofer_cmd, tweeter=tweeter_cmd),
    low_order=(; steering=tt_cmd, dm=ho_cmd),
))
```

- `output_size`
- `output_precision`
- `binning`
- `window_rows` / `window_cols`
- `bits` / `full_well`

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

## How To Choose The Right Entry Surface

Use:

- subsystem functions such as `compute_psf!`, `advance!`, `propagate!`, and
  `measure!`
  - when you are studying one physical layer
- `ClosedLoopRuntime`
  - when you want an AO control loop step surface
- `SinglePlatformConfig` or `GroupedPlatformConfig`
  - when you need explicit exported products and branch composition at the
    runtime boundary

## Next Step

After the first working script:

1. move constants into a small `build_scenario(...)` helper
2. make RNG setup explicit when reproducibility matters
3. add detector realism only when the study needs it
4. move to [user-guide.md](./user-guide.md) for broader workflow guidance
5. move to [api-reference.md](./api-reference.md) when you need the full public
   surface
