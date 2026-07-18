# HIL Plant And Optical Paths

Status: active

## Purpose

This specification defines the shared physical plant, reusable optical paths,
independently scheduled or triggered acquisitions, optical modeling boundary,
and placement of independently controllable optics. Start with
[`../hil-package-boundary.md`](../hil-package-boundary.md) for the product
boundary and document map.

## Plant, Path, And Acquisition Model

A plant declaration contains one telescope and atmosphere, independently
placed controllable optics, reusable optical-path definitions, and
independently scheduled or triggered acquisition endpoints. Separating a path
from its acquisitions prevents a second camera or readout cadence from forcing
duplicate propagation. The names below illustrate the intended decomposition;
they are not committed public API names.

```julia
ground_plane = AtmosphericConjugate(0.0)

plant = HILPlant(
    telescope=tel,
    atmosphere=atm,
    optics=(
        low_order=PlacedControllableOptic(dm_low, ground_plane),
        high_order=PlacedControllableOptic(dm_high, ground_plane),
        high_altitude=PlacedControllableOptic(dm8k,
            AtmosphericConjugate(8_000.0)),
        tiptilt=PlacedControllableOptic(ttm, PupilConjugate()),
    ),
    paths=(
        lgs1=OpticalPath(lgs1_source; train=lgs_train),
        ngs=OpticalPath(ngs_source; train=ngs_train),
        science=OpticalPath(science_source; train=science_train),
    ),
    acquisitions=(
        lgs1_wfs=AcquisitionEndpoint(:lgs1, lgs1_wfs;
            detector=lgs1_detector,
            schedule=PeriodicSchedule(period_ns=2_000_000)),
        ngs_wfs=AcquisitionEndpoint(:ngs, ngs_wfs;
            detector=ngs_detector,
            schedule=PeriodicSchedule(period_ns=1_000_000)),
        science=AcquisitionEndpoint(:science, science_detector;
            schedule=PeriodicSchedule(period_ns=10_000_000)),
    ),
)
```

`OpticalPath`, `AcquisitionEndpoint`, and `PeriodicSchedule` are illustrative
configuration concepts. An acquisition endpoint does not itself prescribe a
cross-owner handoff. The HIL data-plane boundary uses ports backed by the
sequenced SPSC rings specified in [`rtc-ports.md`](rtc-ports.md).

An optical-path declaration owns immutable configuration:

- source direction and spectral definition
- optical-train selection and path visibility
- propagation model, output plane, and sampling geometry
- keys needed to prepare compatible shared results

A prepared path executor owns its mutable OPD, field, detector-facing optical
product, FFT plan, and other scratch state. That workspace has one execution owner and is not stored in the
immutable path declaration.

An acquisition endpoint owns:

- its WFS and/or detector state
- exposure, optical-sample, nondestructive-read, readout, and publication
  schedule or trigger binding
- product sequence, readiness, and acquisition timing metadata
- its output port and buffer-pool policy
- a deterministic RNG stream derived from the central run seed and stable
  endpoint identity

Its declared product geometry may be an area frame, a channel vector, or a
scalar. The schedule and completion-port contracts therefore do not assume an
area sensor; only large products require a leased payload buffer.

Several acquisitions may consume one prepared optical result only when their
source geometry, spectral sampling, optical train, optical sampling instant,
propagation model, and output plane are compatible. Sharing is based on this
prepared path key, not detector type, object identity, or a common nominal
frame rate.

## Optical Plane And Workspace Ownership

The pre-HIL refactor establishes explicit ownership before paths are scheduled
or parallelized:

```text
prepared telescope definition + explicitly advanced atmosphere epoch + optic surfaces
  -> prepared direction/path renderer
  -> caller-owned pupil function or electric field
  -> prepared WFS or science optical front end
  -> caller-owned detector-plane photon-arrival-rate intensity map or bundle
  -> detector-owned temporal integration and acquisition
```

Concrete product names follow the maintained
[`glossary`](../glossary.md); the ownership boundary is normative. A prepared telescope definition owns aperture,
reflectivity, spatial sampling, and geometry. It does not own temporal cadence,
exposure duration, any path's current OPD or electric field, a focal-plane
photon-arrival-rate result, or FFT/intensity scratch. Transitional telescope
OPD remains only for WFS families awaiting decomposition. Scenario, atmosphere, acquisition, and control
owners receive explicit model times or durations instead of reading time from
the telescope. Optical data products likewise do not own propagation plans or
temporary arrays. Prepared workspaces own those plans and scratch arrays, have
one execution writer, and remain fixed in shape, numeric type, backend, and
physical device for a prepared run.

Every prepared optical product carries run-immutable compatibility metadata.
The metadata declares plane kind, coordinate domain, dimensions, sampling,
origin and centering, axis orientation, wavelength or spectral channel, units or
normalization, whether samples represent a spatial density or a cell-integrated
quantity, coherence and combination policy, numeric type, backend, and physical
device where applicable. Preparation rejects incompatible handoffs. Metadata
is not constructed, interpreted through dynamic dispatch, or wrapped around a
newly allocated product on every execution.

Static aberrations and controllable optics form or expose optical surfaces that
apply to an explicit caller-owned path product. They do not mutate a shared
telescope scratch plane. This boundary permits independent paths to share one
telescope definition and permits co-conjugated optics to remain independent;
the later HIL gates still define placement, visibility, command timing, and
effective-state semantics.

Reusable propagating optics follow the same split. In particular, a spatial
filter owns its prepared mask definition separately from its single-writer FFT
workspace and caller-owned filtered-field result. It consumes an explicit field
and does not resize buffers, rebuild masks, or create FFT plans during repeated
execution.

Atmosphere evolution also remains separate from direction rendering. One
writer advances shared layers to an explicit model time and publishes a stable
epoch identity. Advancement never infers a duration from telescope sampling,
detector cadence, or wall clock. Prepared path-local renderers own mutable
single-writer propagation workspace where required, plus frozen source
descriptions and precomputed direction shifts and finite-height footprint
scales. Direction renderers write caller-owned pupil products; field renderers
own their reusable field scratch. Several NGS, LGS, spectral, or extended-source
renderers may consume the same frozen epoch in any order without mutating shared
output, RNG, or geometry-cache state. Source definitions are run-immutable
execution descriptions; mutable user profile/image inputs are frozen during
preparation, and extended-source quadrature expansion belongs to the prepared
renderer or path.

A native direct-science front end consumes an explicit pupil function or field
and writes a caller-owned focal-plane photon-arrival-rate `IntensityMap`. One
compatible rate product may then feed several detectors with independent state
and exposure durations. Each applies its presampling response on the
declared optical grid, integrates over physical pixels, applies QE and elapsed
time exactly once, and then applies its coupling, stochastic response, and
readout. The transitional frame-step shared runtime draws detector noise
sequentially from one runtime RNG, so tuple-order-independent stochastic streams
are not yet claimed. It preflights every detector's exact prepared binding and
idle state before advancing the atmosphere, but unforeseen execution failures
are fail-stop rather than transactional rollback. Changing a
detector exposure cannot change or recompute the optical result. A prepared
external-optics result, including one produced through `Proper.jl`, enters at
the same arrival-rate/acquisition boundary after declaring either photon
irradiance or cell-integrated photon-rate units, or an explicit prepared
conversion from its documented normalization; this does not make `Proper.jl` a
required core dependency.

The maintained native contract is prepared with `prepare_direct_imaging` and
executed with `form_direct_image!`. Same-wavelength `Asterism` components form
one compatible incoherent sum and remain individually accessible through
`direct_imaging_components`; a `SpectralSource` retains wavelength-dependent
grids in an `OpticalProductBundle`. An `ExtendedSource` is expanded explicitly
with `extended_source_asterism` before preparation. The current off-axis model
resolves a finite integer displacement during preparation, honors the focal
grid's declared `:x`/`:y` axis order and signs, and applies a periodic shift;
it does not claim subpixel interpolation or finite-field loss.

Physical complex fields declare a normalization whose squared magnitude is a
photon-arrival-rate product on the represented grid and state whether each
sample is photon irradiance or is integrated over its represented cell.
Dimensionless or otherwise normalized fields and intensity maps remain valid for calibration, estimation, or
synthetic tests, but cannot enter a physical detector acquisition without an
explicit prepared scaling stage. Source photon irradiance also has declared units;
elapsed time is never hidden in telescope or optical-front-end normalization.

Spectral or extended-source components may be accumulated elementwise only
when their physical coordinate grids are compatible and their declared policy
is incoherent intensity addition. Otherwise the optical front end retains a
typed plane bundle or uses an explicitly prepared mapping. Gate 0 does not add
a general image-resampling or coherent-source framework.

These concrete products, prepared functions, and workspaces are sufficient for
the maintained native pipelines. This design does not introduce a universal
optical graph, a general relay language, or a replacement for detailed PROPER
instrument and coronagraph modeling.

## WFS Optical, Acquisition, And Estimation Boundary

Before implementation of the proposed general HIL runtime begins, every
maintained wavefront-sensor family is decomposed across the same semantic
stages without forcing the families into one concrete representation:

```text
caller-owned pupil function or electric field
  -> WFS optical front end
  -> detector-plane photon-arrival-rate intensity map, static tuple, or bundle
  -> detector-owned temporal integration and acquisition
  -> WFS observation or concrete tuple of observations
  -> estimator and calibration
  -> WFS measurement
```

The implemented generic protocol names this boundary directly:

- `prepare_wfs_optical_formation` / `form_wfs_optical_products!`
- `prepare_wfs_acquisition` / `acquire_wfs_observation!`
- `prepare_wfs_estimation` / `estimate_wfs_measurement!`

`WFSObservation` and `WFSMeasurement` bind caller-owned storage to explicit
units and concrete layout/kind, backend, and device metadata. Multi-plane,
spectral, and branch optical sets remain concrete tuples or
`OpticalProductBundle` values. Detector observation sets use one
`WFSObservation` or a concrete tuple of them; there is no runtime-selected
universal graph.
An estimator reports `AcquiredObservationPath()` or
`DirectMeasurementPath()`. The latter intentionally omits optical/acquisition
intermediates rather than representing absent stages with placeholders.

The optical front end owns the sensor-specific propagation and physical optics.
It may produce one detector plane, several simultaneous planes, or a prepared
bundle of channels. A physical detector-facing product is a photon-arrival-rate
product, not an exposure-integrated frame. The front end does not apply
elapsed-time integration, presampling detector response, pixel integration,
charge coupling, stochastic detector response, readout packing, or signal
estimation.

Detector acquisition consumes the rate product and applies each selected model
in its declared spatial, temporal, or charge domain. The baseline prepared
frame path applies presampling response before physical-pixel integration, then
applies QE and its explicit whole-exposure or incremental optical-sample
duration exactly once. A different ordering requires a separately prepared
model whose contract and validation establish the physical meaning. The
contract
supports one detector, several independent detectors, or several regions or
channels of one detector without assuming that every WFS has one two-dimensional
camera. Presampling response, physical-pixel integration, QE and elapsed-time
integration therefore remain downstream of optical spot or pupil-image
formation and upstream of the estimator.

The estimator consumes acquired observations and owns reference subtraction,
normalization, calibration, masking, centroiding, correlation, differential
combination, or phase retrieval as appropriate. Its output is a typed WFS
measurement; the common contract does not claim that every measurement is a
geometric slope vector even when a compatibility accessor uses the name
`slopes`.

The intended family decomposition is:

| WFS family | Optical front end | Estimator |
|---|---|---|
| Shack-Hartmann | independent microlens array and lenslet propagation | centroid, correlation, matched-filter, or other spot estimator |
| Pyramid | modulation policy, pyramid focal-plane optic, and pupil relay | calibrated normalized pupil-image differences |
| BioEdge | modulation policy, complementary edge filters, and pupil relay | calibrated edge-image differences |
| Zernike | phase-shifting focal-plane spot and pupil relay | referenced and normalized pupil-intensity estimator |
| Curvature | prepared intra-/extra-focal propagation branches | calibrated branch-difference estimator |
| LiFT | focal-plane forward model supplied independently of acquisition | iterative phase-retrieval estimator over an acquired observation |

Pyramid and BioEdge may share prepared modulation and focal-plane-filter
machinery where their semantics are genuinely identical, but their physical
optics remain distinct. Curvature acquisition must permit separate branch
detectors as well as packed regions or channels on one detector. LiFT is not
forced into the ordinary `AbstractWFS` hierarchy; its observation acquisition,
forward model, and inverse estimator are independently composable.

The maintained Curvature front end publishes a fixed tuple ordered as positive
then negative defocus. Separate detector acquisitions may have different
presampling response or MTF, physical-pixel integration, QE, exposure,
stochastic response, and readout. A packed acquisition is an explicit mapping,
not an implicit concatenation: both planes must agree in geometry, radiometry,
numeric type, backend, and device, and both branches share the one detector's
exposure. A packed frame uses spatial regions; a packed counting detector uses
channels and requires cell-integrated photon rates. When branch exposures
differ, the estimator receives their prepared inverse rate scales and the
simulation uses separate detectors.

The legacy Curvature branch throughput/background parameters remain an optical
relay rate adjustment for compatibility with the frozen oracle. They can
approximate unequal relay transmission or fixed stray-light offsets, but they
are not a detector response or detector MTF. All detector effects remain after
the branch rate planes and before estimation.

Implementations use dispatch and traits with statically typed prepared
composition. Immutable parameters remain separate from single-writer mutable
workspace and calibration state. Repeated execution writes to caller-owned
storage through explicit mutating functions, performs no scalar GPU indexing,
and does not allocate a new product wrapper or bundle on every sample.

Geometric and reduced-order policies do not allocate or initialize unused
diffractive FFT workspaces. An approximate raw-pixel policy still forms its
declared approximate photon-arrival-rate product and passes through the selected detector
acquisition. A direct slope or modal provider may intentionally bypass optical
and detector stages, but it declares that validity boundary rather than
constructing fictitious intermediate arrays.

This is a narrow WFS composition contract, not a universal optical graph.
Detailed relays and coronagraphs remain eligible for the prepared `Proper.jl`
boundary, and new reusable optical primitives are introduced only when a
maintained WFS or independent optical use justifies them.

## Prepared Product Fidelity

RTC latency and throughput work must not require the full optical plant when
the physical result is not under test. Fidelity is therefore a prepared
per-acquisition product-provider policy, not one global `fast` flag. The names
below describe semantic tiers rather than committed public API types:

| Provider tier | Prepared work | Claims it can support |
|---|---|---|
| Synthetic boundary source | Publish preallocated deterministic patterns, constant products, or recorded products without evaluating the optical path | Canonical-boundary latency, capacity, overload, transport integration, and RTC product-processing throughput |
| Reduced-order plant | Evaluate a command-responsive geometric, linear, modal, low-rank, or reduced-resolution surrogate | RTC reconstruction, tomography, control, loop stability, command-response latency, and surrogate behavior inside its validated envelope |
| Full optical plant | Evaluate the selected atmospheric, diffractive/Fresnel, WFS, detector, and optional PROPER models | Maintained physical and optical claims supported by the model-validity evidence |

Every tier uses the same acquisition schedule or trigger relationship, product
shape and element type, sequence domain, timestamps, completion descriptor,
buffer lease, memory-domain declaration, and full/drop/failure semantics. A
user RTC adapter therefore does not change when the provider changes. Product
providers expose prepared mutating execution through dispatch and traits;
immutable parameters remain separate from single-writer mutable state. The
runtime must not select providers through an `isa` ladder or allocate a new
product object on every event.

A synthetic source declares how it treats its payload: reuse unchanged
contents, touch or regenerate selected elements, copy a prepared product, or
replay a bounded sequence. Its prepared corpus size, rotation, and reuse
distance make cache residency explicit. Metadata-only publication is valid
only when the claimed RTC boundary never reads or transfers the payload. A
raw-pixel RTC throughput claim requires production-shaped buffers in the
correct memory domain and representative payload reads, writes, copies, and
transport; merely rotating descriptors through a ring does not establish pixel
throughput.

A reduced-order provider remains responsive to effective optic commands. It
may use existing geometric WFS measurements for slope products, a calibrated
interaction operator for modal or slope products, or a base image plus
precomputed response images for approximate pixels. Selected detector timing,
noise, presampling response, coupling, and readout stages may still follow that
approximation.
Its validated input domain, approximation error, and excluded physical effects
are recorded. A constant or replayed product is not a closed-loop plant and
must not be used to claim controller stability or optical performance.

### Meaningful Reduced-Order AO Contract

The reduced-order tier is an AO plant emulator, not an arbitrary changing
frame generator. It preserves the causal relationship among disturbance,
effective controllable-optic state, exposure, and sensor product so the
external RTC must perform genuine reconstruction and control. A useful linear
form for path `i` is:

```text
residual_i(t) = projection_i * disturbance(t)
                + sum(command_response_i,j * effective_command_j(t))
product_i(t)  = sensor_i * residual_i(t) + noise_i(t)
```

Each command-response operator includes the calibrated sign convention; the
formula does not assume whether positive RTC coefficients produce positive or
negative optical correction.

The exact implementation may instead use geometric sensing, a nonlinear local
surrogate, or a reduced-resolution model. In every case a prepared reduced-
order scenario declares:

- a deterministic or seeded, time-correlated disturbance state advanced from
  explicit plant time; prerecorded modal disturbances are also valid
- path-specific projections for every active sensing direction, including the
  declared NGS/LGS geometry and tomography relationships
- the same admitted/effective command timeline, hold behavior, device response,
  clipping, misregistration, and path visibility used by the full plant at the
  selected fidelity
- independently timed contributions from every declared controllable optic,
  including co-conjugated, woofer/tweeter, or multi-conjugate devices
- a calibrated sensor operator and conventions matching the product presented
  to the RTC
- scheduled exposure integration, triggering, detector/readout timing, noise,
  and fault behavior selected for the acquisition

Replaying a disturbance sequence is still command-responsive because the
effective correction is applied afterward. Replaying completed sensor products
is not command-responsive and belongs to the synthetic boundary tier.

For a slope boundary, the product can be computed directly with geometric or
calibrated response operators. For a raw-pixel boundary, a geometric spot model
or a base image plus modal-response images can produce approximate pixels
before the selected detector pipeline. The RTC then still performs its real
pixel calibration, slope extraction, reconstruction, tomography, control law,
command splitting or offload, and fault handling.

An accepted command affects only sensor samples whose physical integration
includes the resulting effective optic state. A known stabilizing reference
controller, and later an external RTC, should reduce a declared residual metric
under a matched reduced-order calibration; wrong command sign, excessive delay,
dropped frames, stale commands, or an intentional calibration mismatch must
remain observable rather than being normalized away by the provider.

This tier can therefore support meaningful loop-closure, stability, latency
sensitivity, multi-rate control, tomography, and recovery claims inside its
validated envelope. It does not establish diffraction, aliasing, pyramid
optical-gain, coronagraph, detector-physics, or science-performance behavior
that the selected surrogate omits.

Different acquisitions may deliberately use different tiers in one plant. For
example, full-sized synthetic LGS frames can load an RTC while a reduced-order
NGS pyramid path remains command-responsive and a low-rate science path uses
full PROPER propagation. Provider choice is resolved during preparation,
recorded in the manifest, and immutable while running. A lower-fidelity
overload run is prepared separately; the initial HIL design has no runtime
provider downgrade whose changed physical claim could be mistaken for the
original experiment.

Each controllable-optic declaration owns one physical device, its placement and
path visibility, its response state, and normally one command endpoint. A
device with independently latched command segments may expose more than one
endpoint. Two declarations may reference the same modeled plane without
sharing a command cadence. During preparation the executor may group
compatible declarations by plane and path to avoid repeated propagation; that
derived execution group is not a user-visible command or clock group.

## AO Modes And Optical Path Types

One atmosphere epoch may feed several different path families:

```text
shared atmosphere + telescope epoch
├── NGS or LGS WFS path -> WFS model -> optional detector -> RTC frame
├── another WFS path    -> different WFS/direction/rate -> RTC frame
├── direct science path -> focal-plane photon-arrival-rate map -> science camera
└── coronagraph path    -> residual pupil -> PROPER model -> science camera
```

WFS paths may differ in guide-star direction, NGS/LGS source type, wavelength,
WFS implementation, detector, sampling cadence, and optical train. The current
source and WFS implementations already provide NGS and finite-height LGS
surfaces, including maintained LGS elongation and tomography primitives. The
multi-rate HIL executor must allow mixed NGS and LGS paths without evolving a
separate atmosphere for each guide star.

The paths need not use the same WFS family: a prepared plant may, for example,
combine Shack-Hartmann and pyramid channels when those model/backend
combinations are covered by the support matrix. Path composition does not by
itself promote every WFS/source/backend combination to production support.

Science products are also paths rather than one distinguished camera. A direct
science path can feed one or more compatible maintained detector models from a
shared photon-arrival-rate product. A coronagraph path hands the appropriate
residual pupil field or OPD to a prepared `Proper.jl` model and then applies
the science-camera acquisition semantics to a declared photon-arrival-rate result or
an explicitly prepared normalization conversion. Slow science or coronagraph
work must not delay a higher-rate WFS publication.

Controllable optics have a visibility set in addition to their physical plane:

- SCAO and co-conjugated DM surfaces are visible to every selected path through
  the common pupil train. Each DM retains its own state and command timing even
  when several devices share one modeled plane.
- MCAO uses multiple common DMs conjugated to different altitudes; each source
  sees the appropriate shifted or scaled footprint of every common DM.
- MOAO combines common tomographic WFS paths with a path-specific open-loop DM
  for each corrected science target. A command for one MOAO channel must not
  alter another target's optical path.
- An actuated plane may contain several independently commanded optics inside
  SCAO, MCAO, or an MOAO target path. Sharing a plane is an optical-evaluation
  optimization, not a command group or clock domain.

The conventional terms `woofer` and `tweeter` may be used as configuration
labels for two co-conjugated DMs with different stroke, spatial bandwidth, or
response, but they do not define distinct core types or force synchronized
operation.

Co-conjugated additive pupil behavior, NGS, LGS, multiple WFS arms, direct
science detectors, and the PROPER coronagraph handoff are current foundations.
True altitude-conjugated MCAO and path-specific MOAO are target capabilities
and must not be described as production-supported until their geometry,
command isolation, and hardware paths have maintained evidence.

## Calibration Illumination Without Instrument Assumptions

Calibration illumination is a scenario role, not a privileged core source type
or detector mode. A normal source model or user-provided prepared evaluator
feeds a declared entry boundary of an optical path. Detector-only calibration
may instead use the existing static or temporal detector-input frame seams.
The resulting acquisition retains its ordinary exposure, readout, detector,
trigger, product, and port semantics.

The generic contract records only what the simulator must execute:

- the typed entry boundary and payload it accepts, such as a native electric
  field, a photon-arrival-rate or explicitly normalized intensity map, detector
  input, or prepared external-optics result
- the downstream optical segment and explicit path visibility
- the source evaluator's immutable parameters, single-writer mutable state,
  explicit plant time, deterministic RNG identity, workspace, and backend
  requirements
- source-state changes or setpoints declared by the scenario, without assuming
  whether an RTC, supervisor, script, or free-running device owns them
- explicit combination semantics when more than one illumination contribution
  is active at the same boundary

There is no global `is_calibration` propagation branch and no collection of
`bypass_atmosphere`, `bypass_telescope`, or instrument-selection flags. If a
calibration unit enters downstream of those elements, its path simply begins
at that declared boundary. A descriptive calibration role may be recorded in
the manifest, but it does not change physics through conditionals.

Core does not assume a lamp, laser, fiber, pinhole, integrating sphere,
coherence model, pupil or focal geometry, spectrum, source count, insertion
mechanism, or control protocol. Reusable native models may be added when they
have independent value and evidence. A user model supplies unusual source
physics through dispatch and writes into caller-owned prepared storage. A
detailed relay or coherent calibration unit may be a prepared `Proper.jl` or
other external optical executor.

The core also does not infer how simultaneous contributions combine. A model
that needs coherent field addition, incoherent intensity addition, switching,
or a more specialized rule declares that rule or produces the already-combined
entry payload. Unsupported combinations fail during preparation rather than
being assigned a convenient default.

This is intentionally an incremental extension seam, not a universal optical
graph introduced solely for calibration. New entry-plane payloads are added
only for demonstrated source/path models. Instrument configurations and bench
procedures remain in user code or companion profile packages.

## Pyramid-WFS Modulation

A modulated pyramid WFS contains a fast path-local steering optic, but that
optic is not normally commanded point by point by the external RTC. A local
waveform generator repeatedly drives the modulation path. RTC or supervisory
integration may change modulation radius, frequency, phase, dither, mode, or
start/stop state through a bounded setpoint command endpoint.

The modulation waveform and detector acquisition share an explicit trigger
relationship. The modulator may generate the detector trigger, a detector
trigger may reset modulation phase, both may follow a common trigger source, or
a free-running waveform may expose the phase reference used to place detector
triggers. Per-branch distribution delay, skew, jitter, and dropped or duplicate
edges are modeled by the trigger topology in
[`time-and-scheduling.md`](time-and-scheduling.md), not by streaming waveform
points through an RTC port.

The existing `PyramidWFS` diffractive implementation averages a prepared set of
modulation positions against one unchanged telescope OPD. That remains the
fast baseline cycle-averaged fidelity policy. Its normalized quadrature weights
form an optical average and do not supply an exposure-time factor; detector
acquisition still integrates the resulting rate. A later, profile-driven
time-resolved capability additionally evaluates the modulation phase,
atmosphere epoch, and effective optic state at explicit integration sample
times. It is required only when partial cycles, dither, mirror response, or
trigger-synchronization faults cannot be represented by the cycle-averaged
contract; it does not block the first HIL vertical slice.

Cycle-averaged and time-resolved policies share the same optical calibration
conventions and detector pipeline. Their fidelity, quadrature or sample count,
cycles per exposure, phase reference, and trigger relationship are prepared
configuration and recorded in the run manifest.

## Optical Modeling And NCPA Boundary

`AdaptiveOpticsSim.jl` is an AO-plant optical simulator, not only a timing
shell. Its native scope includes telescope pupil, reflectivity and OPD;
source-aware atmospheric propagation; controllable-surface formation;
electric-field, Fraunhofer, and Fresnel propagation; WFS optics; direct
photon-arrival-rate image formation; spatial filtering; and detector response. Presampling
detector response and charge-coupling effects remain detector models applied
after the incident optical product has been formed. Its reported MTF is the
interior, infinite-grid transfer magnitude of the realized discrete acquisition
kernel. Finite frames use non-amplifying zero extension and therefore have
boundary-dependent response; a continuous subpixel-aperture MTF belongs to an
explicitly oversampled optical-grid model.

The core should own aberrations that are naturally expressed on its sampled
optical planes:

- common or path-specific OPD/phase maps
- modal static NCPA and slowly varying NCPA drift
- sampled amplitude or throughput maps where the native field model is
  sufficient
- explicit path visibility, so a science-only aberration cannot accidentally
  contaminate a WFS path

The current `NCPA` model already synthesizes a static pupil OPD from KL,
Zernike, or external modal-to-command bases and applies it additively or by
replacement. It is a useful optical primitive, but applying it to the shared
telescope state does not by itself express a non-common path. The multi-path
runtime must attach that aberration to the selected science or WFS branch.
Static sampled maps can also be represented directly with `OPDMap`.

Detailed relay and instrument prescriptions remain outside the core. Use
`Proper.jl` when the required result depends on a sequence of physical planes,
powered optics, coronagraph masks, Lyot stops, instrument-specific Fresnel
distances, or an existing PROPER prescription. The preferred seam passes a
typed pupil and residual-OPD or electric-field payload without host copies when
both packages use the same backend.

These choices compose rather than compete. A detailed `Proper.jl` model may
generate an instrument NCPA map that is cached as a fast native HIL surrogate,
or it may propagate the complete science path at that path's lower cadence.
Neither detailed science propagation nor slow NCPA drift may block a high-rate
WFS deadline.

## Conjugated And General Controllable Optics

The current composite-optic behavior sums child OPD surfaces at the telescope
pupil. That numerical behavior remains appropriate for multiple co-conjugated
DMs and low-order pupil-plane optics, but the type also reflects the current
`AOSimulation` constraint that exactly one `AbstractControllableOptic` is
stored. It conflates physical collection, optical application, and command
payload packing, and it cannot represent multi-conjugate placement or
independent timing cleanly.

`CompositeControllableOptic` is therefore removed from the target design rather
than retained as a compatibility wrapper. Its responsibilities separate into:

- a named tuple or registry of individual controllable optics
- placement and path-visibility traits for each optic
- one bounded command endpoint per independently timed optic or segment
- prepared, internal plane groups derived solely for optical execution
- explicit multi-optic transactions when hardware really requires atomic latch
- a user-integration command schema when an external protocol uses packed
  vectors

No new aggregate type should combine all of those responsibilities again.
`RuntimeCommandLayout` should also leave the physical optic API. Its legitimate
responsibilities split into prepared control-output routing in core and a
canonical command descriptor in the HIL companion, with transport packing in
user integration code. None of those responsibilities implies physical
grouping, atomic application, an endpoint, or timing. During implementation the
current composite scenarios remain short-lived numerical references. The
replacement stage deletes `CompositeControllableOptic` and replaces
`RuntimeCommandLayout` rather than retaining either for source compatibility.

Every DM is an independent controllable optic. Its actuator pitch, influence
functions, stroke, spatial bandwidth, misregistration, dynamics, command
endpoint, and update timing remain device-specific. Multiple DMs may reference
the same modeled plane and conjugate transform; their held OPD maps then add.
A prepared field-propagation path should transform to that plane once, apply
every currently effective co-conjugated surface, and continue. This grouping
must not synchronize commands: a faster device may update while a slower one
continues to hold its previous surface.

An optic may also expose one or more sampled feedback acquisitions: encoder or
capacitive-sensor position, measured surface, local-controller state, scalar
health, or another prepared device observation. Feedback has its own product
schema, schedule or trigger, latency, sequence, criticality, and optional noise
model. It is not inferred from command completion and need not share the
command cadence. The same acquisition contract can therefore represent M4/M5,
DM, tip/tilt, or calibration-device feedback without adding a transport- or
instrument-specific interface to core.

`woofer` and `tweeter` are therefore optional roles or configuration labels,
not optical types. Existing segmented-command scenarios define reference
outputs during the replacement, but the target implementation uses individual
devices and endpoints. Any command split, offload, or cross-over policy between
those devices belongs to the external RTC; the plant models the independent
physical responses it is commanded to produce.

A controllable optical plane needs two independent definitions:

- physical optic and actuator response
- geometric or optical placement in each propagated path

For geometric AO propagation, an atmospheric-conjugate mapping applies the
surface using the source footprint at the conjugate altitude. It must account
for source angle, finite-height LGS cone scaling, magnification, rotation,
parity, offsets, and existing misregistration terms.

For diffractive propagation, conjugate altitude alone is insufficient. The
relay and transformations between optical planes determine the field. The core
may provide reusable field and controllable-plane seams, while detailed relay,
coronagraph, and instrument propagation should be delegated to `Proper.jl` or
another optical companion.

MCAO planes are common to every optical path selected by the train. MOAO planes
are path-specific: each science target can have an independently commanded
open-loop DM while sharing the atmosphere and tomographic WFS measurements.
Plane visibility and conjugate geometry are therefore separate traits. Any
common MCAO plane or target-specific MOAO plane may contain more than one
independently commanded device.

Command transport remains independent of plane geometry. The same generic
time-stamped endpoint can drive a pupil DM, atmospheric-conjugate DM, tip/tilt
mirror, focus stage, spatial light modulator, or another supported
`AbstractControllableOptic`; dispatch selects how that optic affects a path.
