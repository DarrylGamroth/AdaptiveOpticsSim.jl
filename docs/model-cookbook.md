# Model Cookbook

Status: active

This cookbook is the shortest path from “I want to model X” to a working
AdaptiveOpticsSim.jl script. Start with one recipe, then move to the
[`user-guide.md`](user-guide.md) or
[`api-reference.md`](api-reference.md) when the model needs more detail.

## Recipe 1: Direct Image From A Telescope And Source

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics

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
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    reference_wavelength_m=500e-9,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction_deg=[0.0, 90.0],
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

Use `PyramidWFS`, `BiOEdgeWFS`, `CurvatureWFS`, or `ZernikeWFS` when the
sensing physics changes. A `ShackHartmannWFS` composes a `MicrolensArray`;
WFS optics, detector acquisition, and estimation remain distinct
stages.

## Recipe 3: Detector-Backed Sensing

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.WavefrontSensors

detector = Detector(
    noise=NoiseReadout(1.0),
    exposure_duration=1.0,
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

For a fast deterministic RTC throughput load, retain the same prepared
frame-detector boundary but omit optional physical effects:

```julia
load_detector = Detector(
    exposure_duration=1e-3,
    qe=1.0,
    noise=NoiseNone(),
    response_model=NullFrameResponse(),
)
load_acquisition = prepare_detector_acquisition(load_detector, detector_rate)
frame = capture!(load_acquisition, rng)
```

Preparation returns the exact detector/input/plan/state/workspace/product
owner. Reuse that value while the input storage is updated in place; prepare a
new owner after replacing any bound storage. The run-immutable plan and
ownership accessors are qualified-public under `Detectors`, not root exports.

This path still applies declared radiometry and exposure. It is not a
camera-specific profile and does not imply reduced optical fidelity upstream.

Counting and channel detectors remain explicit:

```julia
counting_wfs = CurvatureWFS(
    tel;
    pupil_samples=8,
    readout_model=CurvatureChannelReadout(),
)
spad = SPADArrayDetector((2, 64);
    exposure_duration=1.0,
    noise=NoiseNone(),
    sensor=SPADArraySensor(
        active_area_detection_efficiency=0.5,
        fill_factor=0.8,
        dark_count_rate=0.0,
    ),
)
signal = measure!(counting_wfs, pupil, src, spad; rng)
counts = output_frame(spad)
```

Use `LinearAPDDetector` for analog single-element or fixed-bank APD channels,
`SPADArrayDetector` for Geiger-mode accumulated-count images, and
`MKIDArrayDetector` for accumulated-count images. Optional
`MKIDArrayCharacteristics` describe physical resolving power,
photon-arrival-time resolution, and a wavelength passband; they do not add
per-photon energy or timestamp products.

For detailed detector timing, compose the appropriate global-shutter,
rolling-shutter, or frame-transfer acquisition explicitly. Nondestructive reads
and up-the-ramp sampling are detector lifecycle state; frame transfer changes
acquisition timing, not optical performance or MTF.

## Recipe 4: Explicit Closed-Loop Composition

Use a direct loop when you want the reusable numerical and control primitives
without HIL scheduling:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control

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

Use `AdaptiveOpticsSim.Control.VectorDelayLine`,
`AdaptiveOpticsSim.Control.DiscreteIntegratorController`, or a custom controller
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

Each optic retains its own command basis, state, response, and cadence. Optical
addition does not imply synchronized command application. Apply the same rule
to low-order steering stages, focus stages, several MCAO DMs, or path-specific
MOAO DMs.

## Recipe 6: External RTC Lockstep

Prepare a static sensor graph and bind its complete command input and detector
frame output:

~~~julia
using AdaptiveOpticsSim.AlgorithmGraphs

definition = load_algorithm_graph(
    "examples/graphs/revolt_classic_hil_coordinate_gaussian.toml";
    bindings=(;
        pdm_command=zeros(Float32, 277),
        pdm_actuator_coordinates,
    ),
)
graph = prepare_algorithm_graph(definition)
boundary = prepare_graph_hil_boundary(
    graph;
    command_input=:pdm_command,
    frame_output=:shwfs_frame,
)

sequence = step_hil_frame!(boundary)
send_frame(sequence, hil_frame_buffer(boundary))
receive_command!(hil_command_buffer(boundary), sequence)
adopt_hil_command!(boundary, sequence)
~~~

The boundary validates finite complete commands, enforces sequence matching, and
keeps host exchange buffers distinct from graph-owned arrays. PipeWireAO,
shared memory, or another application owns transport and pacing.

A graph with several DMs may expose one structured or packed command input if
that is the declared RTC contract. Preserve the semantic command map and units;
flat-buffer adjacency does not imply that the physical optics are one device.

## Recipe 7: Mixed Fidelity And Rates

Use direct Julia composition when optical paths or acquisitions run at different
rates. The application can combine:

- NGS and LGS paths in different directions
- WFS and science-camera cadences
- native and Proper-backed optical paths
- full optical and explicitly named reduced-order providers
- row-time rolling-shutter samples

Keep one writer for each path, detector, controller, and RNG state. Publish only
complete products. Fidelity and buffer shape should be fixed for a prepared run;
re-prepare to change them.

A version 1 `AlgorithmGraphs` graph is single-rate. Several separately prepared
graphs may be stepped by application code at explicit model times, but the
application then owns their synchronization and shared-state rules.

## Recipe 8: GPU-Resident Work

Choose a backend and exact device when constructing the scientific storage:

~~~julia
import CUDA
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Optics

backend = CUDABackend()
telescope_gpu = Telescope(
    resolution=128,
    diameter=8.0f0,
    central_obstruction=0.1f0,
    T=Float32,
    backend=backend,
)
pupil_gpu = PupilFunction(telescope_gpu; T=Float32, backend=backend)
target = compute_device(opd_map(pupil_gpu))
~~~

Prepare covered graph nodes with `prepare_algorithm_graph(definition; target)`.
Every graph input and parameter must already use native storage on that exact
target. There is no implicit CPU fallback or mixed-device copy.

For a CPU RTC, `PreparedGraphHILBoundary` performs explicit complete-frame and
complete-command copies between the device graph and host exchange arrays.

## Recipe 9: External Coronagraph Model With Proper.jl

Detailed coronagraph propagation belongs at an optional prepared external-
optics seam. A typical composition:

1. form the corrected, path-specific pupil field in AdaptiveOpticsSim.jl
2. pass a typed caller-owned field and immutable prescription parameters to a
   prepared Proper.jl context
3. return a declared photon-arrival-rate `IntensityMap`
4. apply the ordinary detector acquisition path

See
[`examples/support/proper_hil_coronagraph_common.jl`](../examples/support/proper_hil_coronagraph_common.jl)
and [`proper-integration-guide.md`](proper-integration-guide.md). NCPA and
path-specific static surfaces are applied on the branch where they physically
occur.

## Choosing The Entry Surface

Use:

- subsystem functions for optical, detector, calibration, or controller studies
- an explicit Julia loop for generated, multi-rate, conditional, or sub-frame
  models
- `AlgorithmGraphs` for a static single-rate complete-frame graph
- `PreparedGraphHILBoundary` for external-RTC frame/command lockstep
- `SimulationEnsemble` for coarse independent sweeps

See [`runtime-dataflow.md`](runtime-dataflow.md) for graph ownership and
execution order.
