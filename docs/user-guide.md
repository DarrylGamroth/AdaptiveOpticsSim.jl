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

- direct-imaging and optical-product work:
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

- `advance_by!` or `advance_to!`
- `render_atmosphere!`
- `propagate!`
- `update_surface!` and `apply_surface!`
- `measure!`
- `capture!`

Typical use:

- direct photon-arrival-rate image formation
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

### Workflow 1: Optics-only direct image

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
pupil = PupilFunction(tel)
imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
form_direct_image!(imaging)
photon_rate_image = intensity_values(direct_imaging_output(imaging))
```

Use this when you care about:

- pupil construction
- source-scaled focal-plane photon-arrival rates
- image formation
- simple aberration studies

`direct_imaging_output(imaging)` is a caller-owned `IntensityMap` on focal-plane
angular coordinates. Its values are source-scaled, cell-integrated photon
arrival rates before detector exposure, not an inherently unit-normalized PSF.
Preparation binds the pupil, work field, output storage, FFT workspace,
numeric type, backend, and physical device. The `PupilFunction` owns the
mutable path OPD and amplitude; update that same product with `apply_opd!`,
`render_atmosphere!`, or `apply_surface!` before calling
`form_direct_image!` again. `Telescope` owns aperture geometry, not a mutable
optical path.

### Workflow 2: Atmosphere plus one WFS

```julia
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction=[0.0, 90.0],
    altitude=[0.0, 5000.0],
)

wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale_arcsec=0.1, n_pix_subap=6)

rng = runtime_rng(0)
renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)
epoch = advance_by!(atm, 1e-3; rng=rng)
render_atmosphere!(pupil, renderer, atm, epoch)
slopes = measure!(wfs, pupil, src)
```

Atmosphere time is explicit and belongs to the caller. A prepared renderer is
bound to one atmosphere and one frozen source direction; it can render the
current epoch repeatedly without advancing time or consuming RNG. Use
`prepare_atmosphere_renderers` for an `Asterism` or `ExtendedSource`.

For HIL or RTC export, attach a detector and request the detector image after
measurement. `bits` defines the quantization depth, `full_well` defines the
analog-to-digital scaling, and `output_type` defines the Julia array element
type used for the exported frame:

```julia
det = Detector(noise=NoiseNone(), full_well=30_000.0, bits=12, output_type=UInt16)
rng = runtime_rng(0)
measure!(wfs, pupil, src, det; rng=rng)
adu = wfs_detector_image(wfs, det)
```

Here `adu` is a `UInt16` image containing 12-bit ADU values. For
Shack-Hartmann sensors this is the lenslet spot mosaic; for frame-style WFSs it
is the maintained detector/readout frame. Quantized detectors require a fixed
positive `full_well`; the capture path never rescales each frame by its own
peak.

Detector QE may be a scalar or a sampled QE curve. Scalar QE is the simplest
and fastest path. Use a sampled curve when the source spectrum matters:

```julia
qe = AdaptiveOpticsSim.SampledQuantumEfficiency(
    [0.50e-6, 0.60e-6, 0.70e-6],
    [0.35, 0.85, 0.60],
)

det = Detector(noise=NoiseNone(), qe=qe)
image = ones(32, 32) # cell-integrated photon-arrival rate per represented cell
rng = runtime_rng(1)
frame = capture!(det, image, src; rng=rng)
```

For `Source`, detector capture evaluates the curve at `wavelength(src)`. The
generic `capture!(det, image, spectral_source)` boundary uses the
flux-weighted effective QE over the spectral bundle. Pyramid acquisition with
a frame `Detector` instead folds the sampled QE into each wavelength's
optical-rate contribution before the incoherent sum. Prepared diffractive
Shack–Hartmann retains distinct wavelength-rate products in a bundle so each
can use its channel-specific acquisition mapping; the legacy single-product
path remains common-grid only. Matrix-only
capture without a source uses the detector's scalar reference QE, which is the
peak value of a sampled curve.
The bare-matrix path treats its values as cell-integrated photon-arrival rates;
use `IntensityMap` plus `DetectorAcquisitionPlan` when spatial-density versus
cell-integrated semantics must be checked explicitly.

For CMOS, sCMOS, and quantitative low-noise CMOS sensors, compose the generic
`CMOSSensor` from measured properties. Core does not contain camera names or
vendor presets; those belong in a companion profiles package. Uniform
independent read noise uses `NoiseReadout` or `NoisePhotonReadout`, while
row/column components and a heteroscedastic per-pixel map are sensor
properties:

```julia
sigma_map = fill(0.30, size(image))
sensor = CMOSSensor(
    row_readout_sigma=0.05,
    column_readout_sigma=0.08,
    readout_noise_model=CMOSReadNoiseMap(sigma_map),
    timing_model=RollingShutter(10e-6; row_group_size=2),
)
det = Detector(
    sensor=sensor,
    noise=NoisePhoton(),
    qe=0.85,
    full_well=7_000.0,
    bits=16,
    output_type=UInt16,
)
```

`CMOSReadNoiseMap` is an absolute per-pixel read-noise sigma map and is added
to any uniform `NoiseReadout` component. Use `PixelResponseNonuniformity`,
`DarkSignalNonuniformity`, `BadPixelMask`, and `StaticCMOSOutputPattern` for
measured gain, dark, bad-pixel, and output-amplifier structure. CMOS has no
implicit blur: select `RectangularPixelAperture`, another frame response, or
`InterpixelCapacitance` only when the detector sampling or calibration supports
it.

Skipper CCD readout is a CCD sampling mode rather than a photon-counting
detector. It averages nondestructive samples online and retains only the mean,
so memory remains proportional to frame size instead of sample count:

```julia
skipper = Detector(
    sensor=CCDSensor(
        sampling_mode=SkipperSampling(64),
        read_time=20e-6,
    ),
    noise=NoisePhotonReadout(3.0),
)
frame = capture!(skipper, image; rng=runtime_rng(3))
```

The reported effective read-noise sigma scales as `1/sqrt(n_samples)`, and
sampling wall-clock metadata includes all reads. The core model assumes
independent read samples; calibrated correlated-noise and adaptive-read
policies remain future extensions.

For EMCCD cameras, the core package models the generic sensor physics rather
than vendor camera presets. Use `EMOutput()` for the electron-multiplication
register path and `ConventionalOutput()` for a conventional output channel.
Linear EM operation is the default; use `excess_noise_factor=sqrt(2)` when you
want the common high-gain linear-mode excess-noise approximation:

```julia
det = Detector(
    noise=NoisePhotonReadout(30.0),
    gain=300.0,
    sensor=EMCCDSensor(
        output_path=EMOutput(),
        operating_mode=LinearEMMode(),
        excess_noise_factor=sqrt(2.0),
        clock_induced_charge_per_frame=0.01,
        em_gain_range=(1.0, 5000.0),
    ),
)
```

Use `PhotonCountingEMMode` when you want a thresholded photon-counting
approximation for low-flux operation. The threshold is applied after EM gain and
readout noise, so it is expressed in post-EM frame units. Detection efficiency
is a Bernoulli probability for each threshold crossing; accepted events have
unit amplitude rather than using efficiency as an output scale:

```julia
det = Detector(
    noise=NoiseReadout(30.0),
    gain=1000.0,
    sensor=EMCCDSensor(
        operating_mode=PhotonCountingEMMode(threshold=300.0),
        output_path=EMOutput(),
    ),
)
```

`emccd_snr(...)` provides a lightweight analytic check for linear EM,
conventional output, and photon-counting operating modes. Treat it as a design
and validation helper, not a replacement for a calibrated camera model.
`clock_induced_charge_per_frame` is explicitly per frame and is not scaled by
integration time. The default excess-noise model is the fast moment
approximation. `AdaptiveOpticsSim.StochasticMultiplicationRegister` uses a
conditional Gamma model on CPU and a nonnegative moment approximation on
accelerators. Camera-specific parameter packs belong in a companion profiles
package.

Frame transfer is an acquisition-timing policy, not an optical response. Set a
pixel readout rate and transfer time when frame latency or sustained cadence
matters:

```julia
frame_transfer_emccd = Detector(
    integration_time=1e-3,
    sensor=EMCCDSensor(
        readout_rate_hz=10e6,
        acquisition_mode=FrameTransferAcquisition(transfer_time=20e-6),
    ),
)
```

After the first capture, `detector_export_metadata(frame_transfer_emccd)`
reports one-frame output latency as `sampling_wallclock_time` and the overlapped
cadence as `steady_state_frame_period`. `SequentialAcquisition()` instead adds
integration and readout durations. Both modes run the same optical, charge, EM
gain, and noise pipeline.

HgCdTe avalanche arrays likewise have no implicit optical blur or interpixel
coupling. Configure presampling detector response and post-collection IPC as
separate effects. `detector_mtf` reports the normalized discrete-space transfer
magnitude of the realized response kernel on its shift-invariant interior.
Finite frames use zero extension, so edge response is boundary-dependent and
can lose signal outside detector support. The diagnostic does not substitute
for a continuous subpixel-aperture model on an oversampled optical grid:

```julia
det = Detector(
    sensor=HgCdTeAvalancheArraySensor(avalanche_gain=20.0),
    response_model=RectangularPixelAperture(fill_factor_x=0.9,
        fill_factor_y=0.9),
    charge_coupling_model=InterpixelCapacitance(
        [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]),
)
```

`RectangularPixelAperture` records pitch and fill-factor configuration and
applies the resulting discrete detector-grid kernel. It deliberately reports
no subpixel-geometry capability: at unit detector-grid sampling, distinct
continuous apertures can collapse to the same discrete kernel. Prepare an
explicitly oversampled optical mapping when those differences must affect the
image or MTF.

Conventional gain-one HgCdTe arrays and avalanche/SAPHIRA-style arrays support
up-the-ramp fitting:

```julia
ramp_detector = Detector(
    integration_time=1.0,
    noise=NoisePhotonReadout(8.0),
    sensor=HgCdTeAvalancheArraySensor(
        avalanche_gain=1.0,
        read_time=20e-3,
        sampling_mode=UpTheRampSampling(16),
    ),
)

integrated = capture!(ramp_detector, image; rng=runtime_rng(5))
slope = detector_ramp_slope(ramp_detector)
intercept = detector_ramp_intercept(ramp_detector)
read_cube = detector_ramp_cube(ramp_detector)
read_times = detector_ramp_times(ramp_detector)
```

Reads are evenly spaced from zero through the integration time. The returned
frame is `slope * integration_time`, while the fitted slope, intercept, read
cube, and timestamps remain in detector-owned reusable products. The read time
must not exceed the spacing between ramp samples. This linear estimator does
not yet perform cosmic-ray segmentation, saturation-aware fitting, or
correlated-noise estimation.

For a linear-mode single-element APD, use `LinearAPDDetector`. Its channel
storage is a vector rather than a fake 1×1 image. `SingleElementAPD()` accepts
either a scalar photon flux or a one-element vector; `APDChannelBank(n)` uses a
fixed-size vector suitable for preallocated channel readout:

```julia
apd = LinearAPDDetector(
    topology=SingleElementAPD(),
    integration_time=100e-6,
    qe=0.75,
    avalanche_gain=30.0,
    excess_noise_factor=1.2,
    noise=NoisePhotonReadout(2.0),
)
value = only(capture!(apd, 2.0e5; rng=runtime_rng(4)))
```

The existing `APDDetector` remains the Geiger/counting channel model; SPAD and
MKID detectors remain accumulated counting-array models.

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

rng = runtime_rng(3)
frame = capture!(det, pulse; rng=rng)
```

This path samples each rolling-shutter row group at its own readout time, so it
can show transient illumination and rolling-shutter artifacts. Static
`capture!(det, image)` remains the preferred path when the scene does not vary
during the exposure.

For transient sources where the flux rate depends on the full exposure
interval, use `InPlaceExposureFrameSource` or `FunctionExposureFrameSource`.
These receive `start_time` and `exposure_time`, which is important for
global-reset rolling readout:

```julia
pulse = FunctionExposureFrameSource((start_time, exposure_time) -> begin
    active = start_time <= 50e-6 < start_time + exposure_time
    return fill(active ? 1.0 : 0.0, 64, 64)
end)
```

For cameras that use global reset with rolling readout, set
`exposure_mode=GlobalResetExposure()`. This starts all row groups together and
then increases the effective exposure time for later row groups as the rolling
readout reaches them:

```julia
det = Detector(
    noise=NoiseNone(),
    sensor=CMOSSensor(
        timing_model=RollingShutter(25e-6; exposure_mode=GlobalResetExposure()),
    ),
    response_model=NullFrameResponse(),
)
```

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

calibration_pupil = PupilFunction(tel)
imat = interaction_matrix(dm, wfs, calibration_pupil, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = ControlLoopBranch(:main, sim, recon; rng=rng)

cfg = SingleControlLoopConfig(
    atmosphere_step=1e-3,
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

Analytic Gaussian DMs keep a lazy actuator-to-OPD operator instead of storing a
`resolution^2 × n_actuators` sampled matrix. A complete regular grid applies
through the allocation-free factored `X * C * Y'` path; masked, sampled, or
non-separable misregistered Gaussian layouts use a fused matrix-free path.
Dense storage is retained only when explicitly supplied through
`DenseInfluenceMatrix` or `MeasuredInfluenceFunctions`, and calibration code
may materialize the lazy operator during setup when a dense matrix is genuinely
required.
Regular-grid command vectors follow Julia column-major `C[x, y]` ordering, so
the x actuator coordinate is the first and fastest-varying axis.

Actuator print-through is represented only when it is already present in a
sampled influence basis supplied through `DenseInfluenceMatrix` or
`MeasuredInfluenceFunctions`; the built-in Gaussian DM path does not add a
separate print-through model.

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

`SpectralSource`, `ExtendedSource`, and `Asterism` are alternative top-level
source expansions in the maintained API. An `Asterism` is a flat,
common-wavelength directional list and rejects `SpectralSource`,
`ExtendedSource`, and nested `Asterism` children. `with_spectrum` accepts a
`Source` or `LGSSource` leaf and rejects an existing spectral, extended, or
directional expansion. Each `SpectralBundle` weight is the normalized
photon-number fraction of `photon_irradiance(source)` assigned to that sample;
ordinary constructors normalize nonnegative proportional weights to sum to
one. The weights are not radiant-energy fractions. When a model needs a
spectral-by-spatial-by-directional Cartesian quadrature, prepare the components
explicitly and accumulate only metadata-compatible intensity products; there is
not yet a nested convenience API for that product space.

Prepared diffractive Shack–Hartmann formation keeps distinct wavelength grids
as separate native-sampling products in an `OpticalProductBundle`; it never
index-adds or implicitly resamples them. A single output therefore accepts only
a source compatible with its declared wavelength and sampling. The legacy
single-product `measure!` convenience path remains restricted to a common
wavelength grid. Model bundled channels with independent acquisition mappings
unless an application prepares an explicit flux-conserving resampler.

Prepared Pyramid and BioEdge formation follows the same photon-arrival-rate
boundary while retaining distinct physical masks. Their zero, circular, and
sampled focal-plane modulation policies are optical cycle averages in λ/D;
they contain no exposure time or trigger semantics. A spectral source produces
one four-pupil `IntensityMap` per wavelength in an `OpticalProductBundle`.
Directional `Asterism` and `ExtendedSource` inputs require a matching tuple or
vector of path-rendered pupil functions and remain separate bundle products so
direction-dependent atmosphere states are not silently combined. Apply an
explicit compatible incoherent sum or detector mapping only when the intended
instrument path requires it.

The acquired Pyramid/BioEdge estimator accepts a real, square
`:four_pupil_mosaic`, including integer ADU/count frames. It converts samples to
the estimator's floating-point precision before differential arithmetic, and
preparation rejects incompatible frame geometry, backend, or physical device.
Reprepare an optical plan after changing the front-end propagation sampling.

A single diffractive Shack–Hartmann, Pyramid, BioEdge, or atmosphere-aware
Curvature acquisition also requires every asterism leaf to share one optical
calibration signature. In particular, mixed NGS/LGS lists and LGS leaves with
different elongation or sodium-profile geometry belong on independently
calibrated WFS paths.

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
  - composes an independent `MicrolensArray`, prepared optical workspace,
    layout/calibration, detector acquisition, and estimator state
  - use `ShackHartmannOpticalFrontEnd` and `shack_hartmann_rate_map` with the
    prepared WFS stage API when optical formation, acquisition, and estimation
    must be scheduled independently
- `PyramidWFS`
  - pyramid sensing and modulation studies
  - use `PyramidOpticalFrontEnd`, `pyramid_rate_map`, and
    `set_pyramid_calibration!` when optical formation, acquisition, and
    differential estimation must be scheduled independently
- `BioEdgeWFS`
  - BioEdge variants
  - use `BioEdgeOpticalFrontEnd`, `bioedge_rate_map`, and
    `set_bioedge_calibration!` for the corresponding staged path
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
- `MKIDArrayDetector(...)` for accumulated-count imaging with MKID
  energy/timing metadata and an optional meter-valued source passband

MKID source filtering is applied by source-aware `capture!` and WFS/runtime
paths that carry a source through detector capture. Matrix-only capture assumes
the input was already spectrally filtered.

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
pupil = PupilFunction(tel)
ws = AdaptiveOpticsSim.Workspace(pupil.opd, size(pupil.opd, 1);
    rng=deterministic_reference_rng(0))
sprint = AdaptiveOpticsSim.SPRINT(tel, dm, wfs, basis)
```

On CPU, SPRINT uses ForwardDiff-backed DM misregistration sensitivity by
default for grid-backed Gaussian mirrors. Use `sensitivity=:finite_difference`
for validation runs, supported WFS-misregistration finite differences, or
accelerator-backed arrays.

`compute_meta_sensitivity_matrix` returns a structured `MetaSensitivity`, and
`snapshot_config` returns an ordinary string-keyed dictionary. Core neither
reads nor writes a sensitivity cache and does not choose a configuration file
format. If a workflow needs persistence, serialize these returned values
explicitly in application code or through an optional format extension.

## Determinism

- Use a fixed RNG and pass it into `advance_by!`/`advance_to!`, detector calls, or runtime
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
- maintainer/developer navigation:
  - [documentation-map.md](documentation-map.md)
