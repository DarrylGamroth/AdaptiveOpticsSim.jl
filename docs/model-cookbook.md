# Model Cookbook

Status: active

This cookbook is the shortest path from “I want to model X” to a working
AdaptiveOpticsSim.jl script. Start with one recipe, then move to the
[`user-guide.md`](user-guide.md) or
[`api-reference.md`](api-reference.md) when the model needs more detail.

## Recipe 1: Direct Image From A Telescope And Source

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
pupil = PupilFunction(tel)
imaging = prepare_direct_imaging(pupil, src; zero_padding=2)

form_direct_image!(imaging)
photon_rate_image = direct_imaging_output(imaging)
```

The output is a source-scaled, cell-integrated photon-arrival-rate
`IntensityMap` on focal-plane angular coordinates before detector exposure. It
is not an implicitly normalized PSF. Keep `pupil`, the prepared output, and
workspace owned by the path that writes them.

Update the path-owned pupil with `apply_opd!`, `render_atmosphere!`, or
`apply_surface!` before the next `form_direct_image!`.

## Recipe 2: Atmosphere Plus Shack–Hartmann Sensing

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction=[0.0, 90.0],
    altitude=[0.0, 5000.0],
)
wfs = ShackHartmannWFS(
    tel;
    n_lenslets=4,
    mode=Diffractive(),
    pixel_scale_arcsec=0.1,
    n_pix_subap=6,
)

rng = runtime_rng(0)
renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)
epoch = advance_by!(atm, 1e-3; rng)
render_atmosphere!(pupil, renderer, atm, epoch)
signal = measure!(wfs, pupil, src)
```

Use `PyramidWFS`, `BioEdgeWFS`, `CurvatureWFS`, or `ZernikeWFS` when the
sensing physics changes. A `ShackHartmannWFS` composes a `MicrolensArray`;
optical formation, detector acquisition, and estimation remain distinct
stages.

## Recipe 3: Detector-Backed Sensing

```julia
using AdaptiveOpticsSim

detector = Detector(
    noise=NoiseReadout(1.0),
    integration_time=1.0,
    qe=1.0,
    binning=1,
)

measure!(wfs, pupil, src, detector; rng)
frame = output_frame(detector)
```

Use the generic frame detector for CCD, EMCCD, CMOS, sCMOS, configured
quantitative low-noise CMOS, HgCdTe avalanche-array, and Skipper-CCD models.
Choose the sensor and sampling/readout model independently of the optical
front end. `detector_mtf` reports the realized discrete presampling response's
normalized interior MTF.

Counting and channel detectors remain explicit:

```julia
counting_wfs = CurvatureWFS(
    tel;
    pupil_samples=8,
    readout_model=CurvatureCountingReadout(),
)
spad = SPADArrayDetector(
    integration_time=1.0,
    noise=NoiseNone(),
    sensor=SPADArraySensor(
        pde=0.5,
        fill_factor=0.8,
        dark_count_rate=0.0,
    ),
)
signal = measure!(counting_wfs, pupil, src, spad; rng)
counts = output_frame(spad)
```

Use `LinearAPDDetector` for analog single-element or fixed-bank APD channels,
`APDDetector` for Geiger-mode channel counting, `SPADArrayDetector` for
accumulated-count images, and `MKIDArrayDetector` for accumulated count images
with MKID energy-resolution and timing metadata.

For event-driven detector timing, prepare the appropriate global-shutter,
rolling-shutter, or frame-transfer acquisition definition. Nondestructive reads
and up-the-ramp sampling are scheduled detector events; frame transfer changes
acquisition timing, not optical performance or MTF.

## Recipe 4: Explicit Closed-Loop Composition

Use a direct loop when you want the reusable numerical and control primitives
without HIL scheduling:

```julia
using AdaptiveOpticsSim

dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
interaction = interaction_matrix(
    dm,
    wfs,
    PupilFunction(tel),
    src;
    amplitude=0.1,
)
reconstructor = ModalReconstructor(interaction; gain=0.5)
command = similar(dm.state.coefs)

for _ in 1:100
    epoch = advance_by!(atm, 1e-3; rng)
    render_atmosphere!(pupil, renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(command, reconstructor, slopes(wfs))
    @. command = -command
    set_command!(dm, command)
end
```

Model packages may wrap a fixed composition in model-specific `prepare!`,
`step!`, and `readout` functions. The Subaru AO188/AO3k example modules use
that pattern.

Use `VectorDelayLine`, `DiscreteIntegratorController`, or a custom controller
between reconstruction and `set_command!` when the numerical experiment needs
latency or control dynamics.

## Recipe 5: Independent Controllable Optics

Represent a woofer/tweeter pair as independent optics, even when both are
conjugated to the same altitude:

```julia
woofer = DeformableMirror(tel; n_act=4, influence_width=0.5)
tweeter = DeformableMirror(tel; n_act=8, influence_width=0.2)

set_command!(woofer, woofer_command)
set_command!(tweeter, tweeter_command)

update_surface!(woofer)
update_surface!(tweeter)
apply_surface!(pupil, woofer, DMAdditive())
apply_surface!(pupil, tweeter, DMAdditive())
```

Each optic retains its own command basis, state, response, cadence, and
event-loop endpoint. Optical addition does not imply synchronized command
application. Apply the same rule to low-order steering stages, focus stages,
several MCAO DMs, or path-specific MOAO DMs.

## Recipe 6: Prepared HIL Command Routing

The HIL-neutral runtime lives in `AdaptiveOpticsSim.Plant`. A model or companion
package declares a `PlantDefinition`, prepares it with a run seed and one
`CommandEndpointConfiguration` per endpoint, then prepares the event loop.

If an RTC computes one flat vector, expose semantic views and route them to
independent endpoints:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant

rtc_output = zeros(Float32, 5)
products = (
    woofer=@view(rtc_output[1:2]),
    tweeter=@view(rtc_output[3:5]),
)

routing = prepare_controller_output_routing(
    plant,
    products,
    ControllerOutputRoute(:woofer, :woofer_command),
    ControllerOutputRoute(:tweeter, :tweeter_command),
)

route = Plant.controller_output_route(routing, Val(:woofer_command))
command = PlantCommand(
    Plant.controller_output_schema(route),
    42,
    PlantTimestamp(1_000_000),
    Plant.controller_output_payload(route),
)
```

Preparation validates exact type, shape, backend, and physical device. The
route is zero-copy; successful command admission performs the bounded endpoint
copy. Construct a separate `PlantCommand` for each due endpoint so its sequence
and requested effective timestamp remain independent.

Use `PlantCommandTransaction` only when all-or-none application is a real
physical requirement. Flat-buffer adjacency, equal timestamps, and common
conjugation do not create transaction semantics.

Transport is application-owned. TCP, UDP, Aeron, iceoryx2, ZeroMQ, shared
memory, and server/client roles do not change the core command or acquisition
contracts.

## Recipe 7: Mixed Fidelity And Acquisition Rates

One prepared plant may declare several paths and acquisitions:

- NGS and LGS WFS paths in different directions
- science cameras
- native or PROPER-backed coronagraph paths
- calibration-illumination paths
- full-optical, reduced-order, synthetic, or bounded-replay providers

Each acquisition owns its schedule or trigger relationship, detector lifecycle,
product contract, and provider. A fast synthetic WFS can therefore exercise RTC
latency while a slower science arm retains full optical propagation.

Fidelity is fixed during preparation. Re-prepare instead of changing provider
shape or semantics during a run.

## Recipe 8: GPU-Resident Work

Choose the backend when constructing every owner:

```julia
import CUDA

backend = CUDABackend()
tel_gpu = Telescope(
    resolution=128,
    diameter=8.0f0,
    central_obstruction=0.1f0,
    T=Float32,
    backend=backend,
)
pupil_gpu = PupilFunction(tel_gpu; T=Float32, backend=backend)
```

Use `AMDGPUBackend()` on ROCm hardware. Keep the atmosphere, optics, WFS,
detectors, reconstruction storage, controller-output products, and command
endpoints on the same backend and physical device. Copy to host only at a
deliberate transport or inspection boundary.

For a device-resident offline simulation, synchronize only when a dependency or
measurement requires it. For CPU-paced HIL, measure host-ready and device-ready
latency separately.

## Recipe 9: External Coronagraph Model With Proper.jl

Detailed coronagraph propagation belongs at the prepared external-optics seam,
not in a second coronagraph implementation in core. A typical composition:

1. form the corrected, path-specific pupil field in AdaptiveOpticsSim.jl
2. pass a typed caller-owned field and immutable prescription parameters to a
   prepared Proper.jl model
3. return a declared photon-arrival-rate `IntensityMap`
4. apply the ordinary detector acquisition path

See
[`examples/support/proper_hil_coronagraph_common.jl`](../examples/support/proper_hil_coronagraph_common.jl)
for the maintained integration example. NCPA and path-specific static surfaces
are applied on the branch where they physically occur.

## Choosing The Entry Surface

Use:

- subsystem functions for optical or detector experiments
- an explicit model-specific loop for one fixed offline AO model
- `AdaptiveOpticsSim.Plant` for independent virtual-time commands,
  acquisitions, triggers, detector lifecycles, and HIL-neutral execution
- `SimulationEnsemble` for coarse independent sweeps, not a HIL deadline path

See [`runtime-dataflow.md`](runtime-dataflow.md) for ownership and execution
order, and [`hil-package-boundary.md`](hil-package-boundary.md) for the
core-versus-integration boundary.
