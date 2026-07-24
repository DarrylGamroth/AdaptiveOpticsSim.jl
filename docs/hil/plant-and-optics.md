# HIL Plant And Optical Paths

Status: active

## Purpose

This specification defines the shared physical plant, reusable optical paths,
independently scheduled or triggered acquisitions, optical modeling boundary,
and placement of independently controllable optics. Start with
[`../hil-package-boundary.md`](../hil-package-boundary.md) for the product
boundary and document map.

## Plant, Path, And Acquisition Model

A plant ultimately contains one telescope and atmosphere, independently placed
controllable optics, reusable optical paths, and independently scheduled or
triggered acquisitions. The implemented Gate 2 declaration boundary commits
the telescope, atmosphere, path, and acquisition topology. The first four
Gate 4 slices additionally declare every physical controllable optic, define
one immutable versioned semantic schema for each independently timed or
latched command endpoint, and provide standalone bounded admission and
effective-command state owners. The fifth slice prepares every declared optic
and endpoint, composes physical command publication into the serial event
loop, and provides explicit all-or-none transactions across distinct optics.
The current physical composition is deliberately limited to one common
co-conjugated group applied to every path. Device-specific dynamics, explicit
conjugate placement, and path visibility remain their assigned later slices;
none is hidden in an identity or schema.
Separating a path from its acquisitions prevents a second camera or readout
cadence from forcing duplicate propagation.

Gate 2 is validated by the focused topology, preparation, provider, RNG, and
illumination testsets plus the clean [serial plant CPU
artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml). The
artifact composes science, off-axis NGS Shack-Hartmann, and finite-height LGS
pyramid directions with four acquisitions and one shared science result. Its
zero-allocation and HdrHistogram results characterize warmed, self-paced,
in-process service time only; event scheduling, external RTC response, and
fixed-rate load remain later gates.

Optical and acquisition model declarations fail closed until their concrete
types explicitly assert the cold configuration contract:

```julia
const Plant = AdaptiveOpticsSim.Plant

Plant.plant_model_definition_style(
    ::Type{InstrumentOpticalModelDefinition},
) = Plant.ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{InstrumentAcquisitionModelDefinition},
) = Plant.ColdPlantModelDefinition()
```

The assertion excludes prepared plans and workspaces, mutable simulation or
detector state, schedules, RNG streams, queues, transport, and HIL descriptors.
It is intentionally opt-in instead of a shallow mutability heuristic.

```julia
function dm_command_schema(id, endpoint, actuator_count)
    return Plant.PlantCommandSchema(
        Float32,
        (actuator_count,);
        id,
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=Plant.CommandBasis(:actuator, id),
        basis_revision=1,
        semantics=Plant.AbsoluteCommand,
        bounds=Plant.UniformCommandBounds(-5f-6, 5f-6),
        value_policy=Plant.CommandValuePolicy(
            Plant.RejectInvalidCommand,
            Plant.RejectInvalidCommand,
            Plant.ValidateOnPresentation,
        ),
        sequence_policy=Plant.CommandSequencePolicy(
            Plant.RejectSequence,
            Plant.RejectSequence,
            Plant.RejectSequence,
            Plant.AcceptSequence,
        ),
        effective_time_policy=Plant.CommandEffectiveTimePolicy(
            Plant.AllowFutureCommand,
            Plant.RejectLateCommand,
            Plant.PreservePendingCommands,
        ),
        silence_policy=Plant.CommandSilencePolicy(
            Plant.HoldLastCommand,
            Plant.AgeFromApplication,
        ),
    )
end

plant = Plant.PlantDefinition(
    telescope=tel,
    atmosphere=atm,
    controllable_optics=(
        woofer=Plant.ControllableOpticDefinition(
            :woofer,
            woofer_model;
            command_schemas=(
                woofer_command=dm_command_schema(
                    :woofer_v1,
                    :woofer_command,
                    97,
                ),
            ),
        ),
        tweeter=Plant.ControllableOpticDefinition(
            :tweeter,
            tweeter_model;
            command_schemas=(
                tweeter_command=dm_command_schema(
                    :tweeter_v1,
                    :tweeter_command,
                    2_041,
                ),
            ),
        ),
    ),
    paths=(
        lgs1=Plant.OpticalPathDefinition(:lgs1, lgs1_source, lgs_train),
        ngs=Plant.OpticalPathDefinition(:ngs, ngs_source, ngs_train),
        science=Plant.OpticalPathDefinition(
            :science,
            science_source,
            science_train,
        ),
    ),
    acquisitions=(
        lgs1_wfs=Plant.AcquisitionDefinition(
            :lgs1_wfs, :lgs1, lgs1_acquisition),
        ngs_wfs=Plant.AcquisitionDefinition(
            :ngs_wfs, :ngs, ngs_acquisition),
        fast_science=Plant.AcquisitionDefinition(
            :fast_science,
            :science,
            fast_science_acquisition,
        ),
        slow_science=Plant.AcquisitionDefinition(
            :slow_science,
            :science,
            slow_science_acquisition,
        ),
    ),
)
```

These declarations are owned by `AdaptiveOpticsSim.Plant`; the root package
provides no forwarding aliases. `PlantDefinition`,
`ControllableOpticDefinition`, `PlantCommandSchema`,
`OpticalPathDefinition`, `AcquisitionDefinition`, `ControllableOpticID`,
`CommandEndpointID`, `PlantCommandSchemaID`, `OpticalPathID`,
`AcquisitionID`, `RNGOwnerIdentity`, `RNGDerivationVersion`,
`plant_model_definition_style`, and `ColdPlantModelDefinition` are committed
public declaration names. The atmosphere-owned `AtmosphereLayerID` remains in
the root scientific API.
`PreparedPlant`, `PreparedPathExecutor`, `PreparedAcquisitionOwner`,
`PreparedAcquisitionProvider`, `AcquisitionProducts`,
`AcquisitionProductContract`, `PathResultKey`, `prepare_plant`,
`rng_replay_metadata`, `execute_path!`, and `execute_acquisition!` are the
corresponding schedule-free prepared boundary.
A symbol passed as an identity is normalized to the corresponding typed ID. A
tuple or named tuple is only declaration organization: every definition carries
its own explicit identity, a named-tuple key must agree with it, and reordering
cannot change a reference. A named command-schema key must match its target
endpoint. Every command-endpoint identity and every active schema identity is
unique in one plant, and each endpoint has exactly one physical optic owner.
One optic may carry several independently latched endpoint schemas, but shared
optic ownership implies neither packed layout nor atomic application. Multiple
acquisitions may reference the same path, as the two
science acquisitions do above; sampled device feedback likewise remains an
ordinary acquisition rather than a command outcome.

These Julia structs are immutable topology records: their field bindings
cannot be reassigned. The model-definition trait makes controllable-optic,
path, and acquisition model ownership an explicit extension contract;
unrecognized types, including live controllable optics, detectors, and mutable
wrappers, are rejected. The telescope and atmosphere remain separately owned
scientific models with their documented state semantics; in particular, the
atmosphere has one evolution writer. Preparation freezes compatible
configuration and constructs separately owned plans, single-writer workspaces,
and state for the components supported by the current gate.

The definitions contain immutable command payload schemas but no mutable
command state, admission calendar, trigger binding, RNG stream, propagation
workspace, queue, transport, or HIL descriptor. Those mutable and boundary
concepts are attached by their assigned preparation, scheduling, and
HIL-boundary layers. Core names do not carry an `HIL` prefix
because the same plant can run in deterministic virtual time, offline, or
behind the HIL companion. An acquisition definition does not itself prescribe
a cross-owner handoff. The HIL data-plane boundary uses ports backed by the
sequenced SPSC rings specified in [`rtc-ports.md`](rtc-ports.md).

The third and fourth Gate 4 slices prepare standalone endpoint admission and
effective-command owners with fixed payload slots, accepted-sequence history,
future calendar, terminal-disposition storage, and replayable silence policy.
The fifth slice attaches those owners to prepared physical optics and the
virtual-time event loop. `prepare_plant` now requires one
`CommandEndpointConfiguration` for every declared endpoint and fails rather
than silently omitting a physical device or schema.

Preparation turns immutable definitions into backend-, device-, shape-, and
capacity-bound plans plus explicitly owned mutable state, workspaces, and
stateful RNG streams derived from a required central run seed and version.
Repeated execution calls mutating operations over those prepared owners and
caller-owned products. Owner construction validates every concrete execution
plan against the exact input, result, detector state, observation, and optional
estimator storage it retains; warmed calls reject a stale or substituted
binding before mutation. The implemented preparation canonicalizes every
optic and endpoint, freezes every path source, builds a concrete tuple of path
executors with exact atmosphere-input materialization bindings, resolves
acquisitions by stable ID, and builds a concrete tuple of separate acquisition
owners. A cold acquisition selection
canonicalizes caller-selected IDs to stable-ID order and deduplicates their
referenced paths. Its serial executor preflights every path, acquisition,
revision, exact RNG binding, and current epoch before output mutation,
materializes all unique path inputs, forms each path result once, and only then
executes each acquisition. Model construction is multiple dispatch on the opted-in cold
type; there is no
central registry, universal optical graph, abstract executor vector, or stored
closure. The target API does not reintroduce the removed generic frame-step
runtime as a universal object, retain the OOPAO class hierarchy, or hide
scheduling inside optical elements. Frozen numerical oracles remain in the
test harness rather than in the production API.

An optical-path declaration owns immutable configuration:

- source direction and spectral definition
- optical-train selection and path visibility
- propagation model, output plane, and sampling geometry
- keys needed to prepare compatible shared results

A prepared path executor owns its mutable OPD, field, detector-facing optical
product, atmosphere-materialization binding, FFT plan, and other scratch
state. That workspace has one execution owner and is not stored in the
immutable path declaration.

A prepared core acquisition in the current schedule-free slice owns:

- one run-immutable full-optical, command-responsive reduced-order, or
  nonresponsive synthetic/replay provider
- that provider's WFS/detector, reduced-order, payload, or replay state
- caller-owned observation and optional measurement products with required
  logical metadata and one prepared compatibility contract
- an exact read-only binding to one prepared path result and compatibility key

The low-level acquisition receives an explicit RNG when
`execute_acquisition!` is called. Selected serial execution instead supplies
the exact acquisition/detector stream prepared from its stable identity; path
providers and extension-declared device roles have independent streams. The
selection references those plant-owned streams and does not create another RNG
writer.
Exposure sampling events, nondestructive-read events, readiness, publication,
and trigger bindings are Gate 3 or later concerns and are not hidden in this
owner.

The HIL companion binds that core acquisition to a boundary endpoint identity,
completion stream sequence, output port, payload pool/lease policy, external
timing metadata, and resource-specific full behavior. Those boundary owners do
not become fields of the core optical or detector model.

Its declared product geometry may be an area frame, a channel vector, or a
scalar. The schedule and completion-port contracts therefore do not assume an
area sensor; only large products require a leased payload buffer.

Several acquisitions may consume one prepared optical result only when their
source geometry, spectral sampling, optical train, optical sampling instant,
propagation model, and output plane are compatible. Sharing is based on this
prepared path key, not detector type, object identity, or a common nominal
frame rate. `PathResultKey` also binds radiometry, telescope/model revisions,
numeric/output geometry, backend, and physical device. The default
`InstantaneousOpticalSample` describes a photon-arrival-rate result at one plant
instant; it is neither a detector exposure nor a cadence. An acquisition calls
`require_path_result` before constructing destinations when it has stricter
requirements. The prepared owner then retains the exact key and result storage,
so two unequal-exposure detectors can share one formed optical rate without
sharing detector state.

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
[`glossary`](../glossary.md); the ownership boundary is normative. A prepared
telescope definition owns aperture, reflectivity, spatial sampling, and
geometry. It does not own temporal cadence, exposure duration, any path's
current OPD or electric field, a focal-plane photon-arrival-rate result, or
FFT/intensity scratch. Every maintained WFS, science, calibration, atmosphere,
and controllable-optic path now consumes an explicit path product; there is no
transitional telescope-owned OPD. Scenario, atmosphere, acquisition, and
control owners receive explicit model times or durations instead of reading
time from the telescope. Optical data products likewise do not own propagation
plans or temporary arrays. Prepared workspaces own those plans and scratch
arrays, have one execution writer, and remain fixed in shape, numeric type,
backend, and physical device for a prepared run.

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
Gate 4 defines independent endpoint command timing and effective-state
semantics. Later HIL gates still define physical placement, path visibility,
device dynamics, and boundary pacing.

Reusable propagating optics follow the same split. In particular, a spatial
filter owns its prepared mask definition separately from its single-writer FFT
workspace and caller-owned filtered-field result. It consumes an explicit field
and does not resize buffers, rebuild masks, or create FFT plans during repeated
execution.

Atmosphere evolution also remains separate from direction rendering. One
writer advances shared layers to an explicit model time and publishes an
`AtmosphereEpoch` token for the current state. Advancement never infers a
duration from telescope sampling, detector cadence, or wall clock. The token
is identity and time metadata, not retained layer storage: the current
renderers reject it after the same atmosphere advances.

Prepared path-local renderers own mutable single-writer propagation workspace
where required, plus frozen source descriptions and precomputed direction
shifts and finite-height footprint scales. Direction renderers write
caller-owned pupil products; field renderers own their reusable field scratch.
Several NGS, LGS, spectral, or extended-source renderers may consume the same
current epoch in any order without mutating shared output, RNG, or geometry-
cache state. The serial scheduled executor therefore materializes every due
path's atmosphere-dependent OPD, field, or model-specific input before the
atmosphere writer advances to a later timestamp. Downstream optical and
detector work may outlive that epoch once it reads only the materialized
caller-owned product.

The implemented schedule-free serial seam prepares an explicit acquisition
selection and executes it from either a caller-supplied current
`AtmosphereEpoch` or an explicit absolute atmosphere model time. It resolves
each referenced path once and performs a complete mutation-free preflight over
the selection before writing the first input. It does not decide which
acquisitions are due, infer time, retain old atmosphere layers, or introduce a
scheduler. The later event scheduler supplies its due set to this boundary.

A path that must render from older mutable layer state requires a separately
prepared, bounded retained atmosphere-state snapshot whose representation is
chosen by the atmosphere model. A parallel executor may instead hold the
writer at one epoch while same-epoch readers finish. It MUST NOT treat an
`AtmosphereEpoch` token as an immutable copy, and it MUST NOT advance the
writer while an unmaterialized consumer still reads the current layers.
Source definitions are run-immutable execution descriptions; mutable user
profile/image inputs are frozen during preparation, and extended-source
quadrature expansion belongs to the prepared renderer or path.

A native direct-science front end consumes an explicit pupil function or field
and writes a caller-owned focal-plane photon-arrival-rate `IntensityMap`. One
compatible rate product may then feed several detectors with independent state
and exposure durations. Each applies its presampling response on the
declared optical grid, integrates over physical pixels, applies QE and elapsed
time exactly once, and then applies its coupling, stochastic response, and
readout. The schedule-free prepared plant requires a run seed and derivation
version, binds stable identities to atmosphere layers, providers, detectors,
and declared device roles, and supplies their independent streams as specified
by the [`determinism policy`](../deterministic-simulation.md). Multilayer plant
atmospheres require explicit layer identities so reordering layers does not
reassign their streams. The selected serial execution performs the binding and
epoch preflights described above before product mutation; unforeseen execution
failures remain fail-stop rather than transactional rollback. Changing a
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

An `OpticalProductBundle` fixes its product membership when constructed while
leaving each product leaf's numerical array mutable for its single prepared
writer. A prepared owner therefore cannot silently substitute, reorder, append,
or remove acquisition-facing branches after its compatibility contract is
built.

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

Before implementation of the general HIL runtime began, every maintained
wavefront-sensor family was decomposed across the same semantic stages without
forcing the families into one concrete representation:

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

The current `capture!(...; integration_duration=...)` and
`capture_incremental!` methods are frame-step convenience APIs. Their duration
is expressed in seconds; it is not an absolute sample timestamp. They
accumulate until the configured exposure duration is reached and then finalize
automatically using floating-point tolerance. The virtual-time event layer
does not use that tolerance as its completion authority. Its initial global-
shutter implementation prepares explicit
begin-exposure, accumulate-interval, nondestructive-read, close-exposure,
readout-complete, and publish operations, with integer plant timestamps owned
by the scheduler and exact interval durations supplied to detector physics.
The prepared definition and separately owned state reject busy retriggers,
noncontiguous intervals, exposure-boundary crossings, missed ramp reads, and
out-of-order close/readout/readiness transitions before detector mutation.

The existing `UpTheRampSampling` whole-exposure path synthesizes its read cube
after final integration by scaling the completed frame and applying per-read
noise. It is retained as a declared post-exposure, lower-fidelity convenience
and frozen oracle. A scheduled up-the-ramp acquisition now snapshots the
evolving accumulated charge at each nondestructive-read event without ending
the integration. Stochastic avalanche response is applied to each newly
generated charge interval and retained in later reads rather than redrawn for
previously observed charge. Atmosphere, source, and effective-command changes
inside the exposure therefore affect the applicable reads. Only that path supports
time-resolved nondestructive-read claims.

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
Its prepared forward model publishes a cell-integrated focal-plane photon rate
and binds deterministic spatial preprocessing without exposure or QE. A
caller-owned `LiFTObservation` declares whether acquired values are photon
rates, expected detected counts, or normalized intensity, including the
explicit conversion back to the estimator's canonical photon-rate domain.
Modal selection is prepared with the estimator and is not supplied on each
reconstruction call.

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
per-acquisition product-provider policy, not one global `fast` flag. The
implemented trait values are `FullOpticalProviderStyle`,
`CommandResponsiveReducedOrderProviderStyle`, and
`SyntheticReplayProviderStyle`:

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

The core implementation binds one `PreparedAcquisitionProvider` into each
immutable `PreparedAcquisitionOwner`. `AcquisitionProductContract` snapshots
the observation and measurement shape, numeric type, backend/device, typed
units and metadata, and required acquisition-level metadata. Providers mutate
and return the exact caller-owned `AcquisitionProducts` value; incompatible
shape, numeric type, memory domain, units, or metadata fails during
preparation. A selected reduced-order or synthetic/replay acquisition does not
cause its referenced full-optical path to execute.

A synthetic source declares how it treats its payload: reuse unchanged
contents, touch or regenerate selected elements, copy a prepared product, or
replay a bounded sequence. Its prepared corpus size, rotation, and reuse
distance make cache residency explicit. Metadata-only publication is valid
only when the claimed RTC boundary never reads or transfers the payload. A
raw-pixel RTC throughput claim requires production-shaped buffers in the
correct memory domain and representative payload reads, writes, copies, and
transport; merely rotating descriptors through a ring does not establish pixel
throughput.

Core maintains unchanged-reuse, owned-snapshot copy, and fixed-size cyclic
completed-product replay implementations. The replay corpus has fixed capacity
and a separate mutable single-writer cursor. Provider payload-work metadata
records `:reuse_unchanged`, `:copy_prepared_product`, or
`:bounded_cyclic_replay`; these declarations are not RTC load evidence.

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

Core implements this linear form for intentional direct measurements.
`HarmonicDisturbanceModel` supplies exactly replayable explicit-time
disturbance coordinates, `ReducedOrderCommandResponse` binds every currently
visible endpoint to an operator and exact command-schema semantics, and
`LinearReducedOrderAcquisitionModel` supplies the path projection and sensor
operator. `DirectMeasurementAcquisitionDefinition` exposure-averages the held
instantaneous `WFSMeasurement` and retains ordinary periodic or delivered
trigger start, readout-completion, readiness, product-sequence, and
right-continuous effective-command timing. Reduced-order-only samples bypass
unused full-optical path and atmosphere execution.

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
calibrated response operators. For a future raw-pixel boundary, a geometric
spot model or a base image plus modal-response images could produce
approximate pixels before the selected detector pipeline. The maintained
built-in linear provider currently supports direct measurement vectors only;
it does not claim approximate raw-pixel or detector behavior. When such a
model is added, the RTC can still perform its real pixel calibration, slope
extraction, reconstruction, tomography, control law, command splitting or
offload, and fault handling.

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

The maintained acceptance fixture currently establishes deterministic modal
loop closure, independent endpoint timing, clipping through the ordinary
command layer, and expected degradation under wrong sign, excess delay, stale
feedback, dropped updates, and a response-calibration mismatch. It does not
yet establish tomography, raw-pixel processing, external-RTC behavior,
instrument-scale capacity, integrated GPU event-loop execution, or comparison
to a full-optical reference.

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

That schedule-free contract is now implemented by
`prepare_illumination_entry` and the typed
`PupilFunctionIlluminationEntry`, `ElectricFieldIlluminationEntry`,
`IntensityMapIlluminationEntry`, `ExternalOpticsResultIlluminationEntry`, and
`DetectorInputIlluminationEntry` tags. A `PreparedIlluminationEntry` is an
ordinary path materializer: it receives the selected atmosphere epoch's plant
time and a stable path-owned `:illumination` RNG, writes one exact caller-owned
product, and hands that product to the path's normal prepared execution.
`SingleIllumination`, `ExclusiveIlluminationSelection`,
`CoherentFieldCombination`, and `IncoherentIntensityAddition` are explicit
declarations rather than inferred behavior. The small native
`UniformIntensityIllumination` model has CPU and maintained accelerator-path
evidence through ordinary detector acquisition; a user-defined stateful pupil
model separately covers explicit time/RNG use and normal downstream optical
formation on CPU.

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

This implementation does not yet add HIL scheduling, triggers, calibration
setpoint commands, sequence/timestamp descriptors, transport, or a source-state
mode manager. Those remain in their owning later gates; the evaluator's
explicit plant time and state are sufficient for deterministic schedule-free
execution without assigning control authority to core.

## Pyramid WFS Modulation

A modulated Pyramid WFS contains a fast path-local steering optic, but that
optic is not normally commanded point by point by the external RTC. A local
waveform generator repeatedly drives the modulation path. RTC or supervisory
integration changes declared bounded setpoints rather than sending the
waveform samples.

The native contract is `CircularPyramidModulator`. It requires four distinct
absolute scalar endpoint schemas in the
`CommandBasis(:waveform_setpoint, :pyramid_modulation)` basis:

- nonnegative radius in λ/D
- nonnegative frequency in hertz
- counterclockwise phase offset in radians
- binary disabled/enabled state

Each schema has finite bounds and its own endpoint state, future calendar,
command-silence policy, optional safe value, and effective timestamp. Enabled
state is the native start/stop setpoint. Dither, waveform mode, or another
device-specific setpoint requires a model that explicitly declares its own
endpoint rather than a universal heterogeneous command vector.

`AutonomousPeriodicOpticDefinition` binds the prepared modulator to exactly one
full-optical Pyramid path, one immutable fidelity policy, and one phase
relationship:

- `FreeRunningPhaseReference` uses a configured plant-time origin
- `TriggerSourcePhaseReference` follows realization of one trigger source,
  representing a local waveform reference that can also feed a fan-out
- `TriggerResetPhaseReference` follows physical delivery to one consumer, so
  its branch independently observes link delay, skew, jitter, phase steps,
  drops, and duplicates

Preparation rejects an unknown source or consumer, a trigger consumer assigned
to more than one event owner, a path-local autonomous optic without a binding,
and two native modulators that would write the same prepared Pyramid
quadrature. Timestamp-label faults remain reporting effects and do not move a
physical reset. These relationships reuse the trigger topology in
[`time-and-scheduling.md`](time-and-scheduling.md); no RTC port streams
waveform points.

The implemented `CycleAveragedModulationFidelity` regenerates the already
allocated circular quadrature only when radius or enabled state changes. It
preserves point count, backend, numerical quadrature origin, and normalized
weights. Disabled state centers every quadrature point. Frequency, phase, and
the trigger-relative analytic phase remain observable device state but cannot
alter a complete-cycle average. A command effective at an optical-sample
timestamp changes that sample because trigger updates precede command
application and both precede optical execution.

The existing diffractive `PyramidWFS` then averages that quadrature against one
unchanged path-owned pupil OPD. The result remains a detector-plane
photon-arrival-rate product: quadrature weights do not supply exposure time,
and detector acquisition still integrates its half-open exposure exactly once.
Changing the operating setpoint does not silently rebuild a WFS calibration
reference or change the immutable fidelity policy.

Time-resolved modulation remains profile-driven. A profile needs that separate
fidelity when it claims partial-cycle effects, within-exposure atmosphere or
optic evolution, dither trajectory, mirror servo response, or phase-dependent
trigger faults that a complete-cycle average cannot represent. No
time-resolved or physical steering-servo claim is implied by the native
baseline.

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
replacement. It is a path-agnostic optical primitive; applying it to a
caller-owned `PupilFunction` selects the affected path. A multi-path runtime
must therefore attach that aberration to the selected science or WFS branch.
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

Co-conjugated DM and low-order OPD surfaces still add at the applicable pupil
plane, but physical collection, optical application, and command payload
mapping are separate concerns. The former aggregate optic and generic
single-optic runtime conflated those responsibilities and could not represent
multi-conjugate placement or independent timing cleanly.

The aggregate optic has been removed rather than retained as a compatibility
wrapper. Its responsibilities are now separated into:

- a named tuple or registry of individual controllable optics
- placement and path-visibility traits for each optic
- one immutable semantic command schema per independently timed optic or
  segment, followed by one bounded prepared endpoint owner
- prepared, internal plane groups derived solely for optical execution
- explicit multi-optic transactions when hardware really requires atomic latch
- a user-integration command schema when an external protocol uses packed
  vectors

No new aggregate type should combine all of those responsibilities again. The
former packed command-layout API has also left the physical-optic surface. Its
legitimate responsibilities are split into implemented prepared controller-
output routing and plant command schemas in core, a later command-submission
descriptor schema in the HIL companion, and transport packing in user
integration code. None implies physical grouping, atomic application, an
endpoint, or timing. Frozen composite cases remain test-only numerical
references; no compatibility adapter is retained.

`ControllerOutputRoute` and `prepare_controller_output_routing` now bind every
named caller-owned controller product to one distinct prepared endpoint.
Products may be views into a larger controller output and are borrowed without
packing or copying. Preparation validates exact numeric type, shape, backend,
and physical device. Each integration layer then constructs an independent
`PlantCommand` with its endpoint-local sequence and requested effective
timestamp; successful admission is the bounded payload-copy boundary.

The implemented `PlantCommandSchema` defines semantic interpretation:
stable schema/version and endpoint identities, exact scalar or backend-neutral
array element type and dimensions, physical-unit and sign-convention
identities, basis and revision, absolute or incremental meaning, uniform or
unbounded values, and explicit value, sequence, effective-time, supersession,
and plant-time-silence policy vocabulary. A complete fixed-shape replacement
uses absolute semantics; it is not a third semantic mode. Stroke or slew checks
that depend on the currently effective state are application-stage device
rules. `validate_plant_command_payload` checks presentation compatibility
without mutation, clipping, admission, sequencing, scheduling, or application.

Core now represents one mapped presentation with `PlantCommand` and prepares a
standalone `PreparedCommandEndpoint`. Its separately owned state uses fixed
endpoint-owned payload slots, a fixed accepted-sequence window, and a flat
time-ordered future calendar. Successful admission copies the caller payload;
only that success enters sequence history. A full or time-rejected command can
therefore retry the same endpoint sequence as a new presentation identity.
Presentation-time clipping occurs only in transactional staging storage. An
array copy or clip failure cannot corrupt an existing pending payload.

Duplicate, stale, reordered, skipped, future, late, full-capacity, and absolute-
command supersession decisions are explicit. Pending incremental deltas may not
be discarded by the supersession policy. Late apply-now retains the requested
timestamp, schedules at the current plant timestamp, and records lateness; no
transition backdates plant state. Calendar order is scheduled plant time,
stable prepared endpoint ordinal, then endpoint-local command sequence. Every
admission rejection and displaced pending command writes one terminal
disposition into caller-owned fixed storage.

The single writer may claim at most one opaque application-ready command,
borrow its endpoint-owned payload read-only, and then report applied or failed
exactly once. Claim completion revalidates its identity, order key, and
timestamps against endpoint-owned state, and terminal records are constructed
from that state rather than caller-supplied claim fields. A bounded drain fails
all still-pending commands in calendar order.

The qualified `CommandApplicationState` now separately owns the endpoint's
explicit initial/effective command, preallocated application staging, optional
copied safe command, last application time, and command-silence latch.
`apply_claimed_plant_command!` transactionally replaces an absolute value or
adds an incremental delta, rejects any nonfinite result, enforces the declared
application-stage bound action, requires application at the command's
immutable scheduled timestamp, and publishes exactly one disposition. A
missed timestamp fails the presentation instead of backdating state. Rejected
or failed candidates do not change the effective command. Array
storage remains on the endpoint's prepared backend and commits by exchanging
the effective and staging buffers.

Before the first qualifying command event, the endpoint's configured initial
timestamp is the common silence-age origin. Thereafter a successful admission
rebases `AgeFromAdmission`, while successful effective application rebases
`AgeFromApplication`. Indefinite hold schedules no transition. Safe and fail
policies schedule one exact transition for each new selected origin. A due
command at the same timestamp must be applied or failed first. The safe action
copies the prepared safe value into effective state without manufacturing a
command disposition; the fail action drains all pending presentations and
marks the endpoint terminally failed. The returned
`PlantCommandSilenceTransition` is the replayable state/fault record.

Direct use of the standalone application layer does not itself mutate a
physical optic or model stroke, slew, settling, hysteresis, or other dynamics
beyond the declared effective-value policy. The composed event loop now binds
each endpoint to one prepared physical optic: it stages effective-command and
physical state before bounded publication and applies visible surfaces before
path execution. `PlantCommandTransaction` explicitly provides all-or-none
admission and publication for endpoints on distinct optics at one common
timestamp. Atomicity is not inferred from placement, equal time, or packed
storage.

The HIL command-submission descriptor wraps a mapped plant command with
run/session correlation, source timestamp-domain and mapping metadata, payload-
lease ownership, and outcome credit. Its paired command outcome wraps the core
disposition with boundary timing and returns that credit. Enqueue, payload
validation, semantic admission, physical application, and boundary completion
therefore remain separate observable operations, and core never imports a HIL
descriptor or lease type.

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
