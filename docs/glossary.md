# Project Glossary

Status: active

This is the canonical terminology guide for public APIs, implementation names,
documentation, requirements, and tests. Plane names describe optical location;
product names describe the represented quantity. New scientific terms should be
defined here before they become public API.

## Naming Rules

- Do not combine a plane location and an assumed physical quantity in one type
  name. Use, for example, an `IntensityMap` whose metadata says `FocalPlane()`.
- State coordinate domain and units explicitly. Do not infer metres or radians
  from plane kind alone.
- Use *rate*, *exposure*, and *count* according to whether time has not yet been
  integrated, has been integrated per area, or has been integrated per sample.
- Use *MTF* only for a frequency-domain modulation transfer function. A kernel
  applied to an image is a spatial-domain response model.
- Use *WFS signal* for a sensor-independent typed measurement. Use *slope* only
  when the values are wavefront gradients or calibrated centroid-equivalent
  slopes.

## Optical Planes, Coordinates, And Products

| Term | Project meaning |
|---|---|
| Pupil plane | A transverse plane conjugate to the aperture stop. Represented by `PupilPlane()`. |
| Focal plane | A physical focus or focal-mask plane. Represented by `FocalPlane()`. Its coordinates may be angular when that is the declared propagation output. |
| Image plane | An image-forming plane. Use this term in prose when the distinction from a focal-mask plane matters; introduce a concrete marker only with an implemented contract. |
| Detector plane | The optical sampling plane incident on a detector, represented by `DetectorPlane()`. It is not the detector readout frame itself. |
| Intermediate plane | A known transverse optical plane that is not more specifically a pupil, focal, image, or detector plane. Represented by `IntermediatePlane()`; do not use *propagation plane*. |
| Metric coordinates | Plane `sampling` and `origin` are in metres. Represented by `MetricCoordinates()`. |
| Angular coordinates | Plane `sampling` and `origin` are in radians. Represented by `AngularCoordinates()`. |
| Normalized pupil coordinates | Dimensionless coordinates normalized to the entrance-pupil diameter, represented by `NormalizedPupilCoordinates()`. For a re-imaged-pupil detector mosaic they describe the optical relay sampling, not a manufactured detector's physical pixel pitch. |
| Sample-centered | A geometric origin lies on a stored sample. |
| Inter-sample-centered | A geometric origin lies between stored samples. This is not necessarily an inter-pixel convention. |
| Telescope aperture | The revisioned support, intensity reflectivity, diameter, and spatial sampling prepared by `TelescopeAperture`. It is geometry shared while paths are prepared; it is not a mutable wavefront or OPD buffer. |
| Pupil function | The complex pupil transmission in Fourier optics. `PupilFunction` snapshots one telescope-aperture revision, stores its wavelength-independent support and field amplitude, and owns the mutable OPD for one path; field formation realizes the wavelength-specific complex pupil. |
| Pupil reflectivity | The finite dimensionless intensity-throughput fraction in `[0, 1]` stored by `TelescopeAperture`. Optical field formation applies its square root as an amplitude transmission. It is not an amplitude coefficient itself. |
| Electric field | In this package, a scalar complex optical field envelope/phasor sampled on a declared plane. `ElectricField` is not an SI vector electromagnetic field in V/m unless its metadata explicitly establishes such a normalization. |
| Intensity | The computational quantity `abs2(E)`. Its physical units and normalization follow the field and metadata; the word alone does not assert radiometric irradiance. |
| Intensity map | A real sampled array of intensity values plus `OpticalPlaneMetadata`, represented by `IntensityMap`. The metadata, not the type name, states plane location and physical semantics. |
| Photon-arrival-rate intensity map | An `IntensityMap` with `PhotonRateNormalization()` and an explicit spatial measure. It represents a rate before detector exposure: either a spatial density or a rate integrated over each represented cell. |
| Direct imaging | Native image formation from a declared pupil function or preformed pupil-plane electric field to a caller-owned focal-plane photon-arrival-rate `IntensityMap`. `prepare_direct_imaging` fixes storage, grid, backend, device, and workspace; `form_direct_image!` executes without applying detector exposure time. The current off-axis model resolves a finite integer displacement during preparation, honors the focal grid's declared `:x`/`:y` axis order and signs, and applies a periodic shift. |
| Optical path difference (OPD) | Scalar optical path-length difference in metres. Positive-sign convention must be validated at integration boundaries. |
| Point-spread function (PSF) | The image-plane response to a point source under an explicitly declared normalization. A source-brightness- and exposure-scaled image is not intrinsically a PSF. |
| Non-common-path aberration (NCPA) | Aberration visible to one optical branch but not to its reference or sensing branch. Model it in the branch where it physically occurs. |

`UnspecifiedSpectralCoordinate()` means that no spectral coordinate has yet been
declared; it does not mean wavelength independence.
`AchromaticSpectralCoordinate()` explicitly declares wavelength-independent
geometry or data. `MonochromaticChannel(λ)` declares one wavelength in metres.
`IntegratedSpectralChannel(id)` declares that values have already been integrated
over an application-defined passband identified by `id`; equal identifiers are
the compatibility contract, while the passband definition remains application
owned.

## Radiometry And Photon Quantities

| Term | Unit and use |
|---|---|
| Radiant irradiance | Incident radiant flux density, W·m⁻². Use *irradiance* without a photon qualifier only when this radiometric meaning is intended. |
| Photon irradiance | Incident photon-rate density, photons·s⁻¹·m⁻². `photon_irradiance(source)` uses this quantity. |
| Spectral sample weight | The normalized photon-number fraction of `photon_irradiance(source)` assigned to one `SpectralSample` in a `SpectralBundle`. It is not a radiant-energy fraction. Ordinary bundle constructors normalize nonnegative proportional input weights to sum to one. |
| Photon flux | Photon rate, photons·s⁻¹, integrated over a stated receiving area or channel. Do not use it for a per-area source value. |
| Photon-arrival-rate product | Detector-facing umbrella term when metadata may declare either a spatial density or a cell-integrated photon rate. |
| Photon exposure | Time-integrated photon irradiance, photons·m⁻². |
| Expected photon count | Expected photons integrated over both time and represented cell area. It is the mean supplied to photon-counting statistics. |
| Relative or normalized intensity | Dimensionless intensity under a declared normalization, such as unit total power, peak, or contrast. It cannot enter a physical detector acquisition without an explicit scaling contract. |
| Point-sampled measure | `PointSampledMeasure()`: each value is evaluated at a coordinate and includes no represented cell area. It cannot enter detector acquisition without an explicit spatial mapping or integration. |
| Spatial-density measure | `SpatialDensityMeasure()`: each value is per unit area of the declared coordinate domain. Detector acquisition multiplies it by the represented cell area before exposure integration. A photon-rate density on metric coordinates is photon irradiance. |
| Cell-integrated measure | `CellIntegratedMeasure()`: each value has already been integrated over its represented spatial cell. A photon-rate-normalized value therefore has units photons·s⁻¹ per cell. |
| Coherent field combination | `CoherentFieldCombination()`: compatible complex fields may be combined before intensity formation, so relative phase is retained. |
| Incoherent intensity addition | `IncoherentIntensityAddition()`: compatible intensity maps may be added elementwise after their grid, spectral, radiometric, backend, and device contracts have been checked. |
| Non-combinable product | `NonCombinableProduct()`: the product must remain separate unless an explicit physically justified mapping is prepared. |
| Optical product bundle | `OpticalProductBundle`: a typed, fixed-membership collection that preserves products on incompatible grids or channels without implicit resampling, radiometric conversion, or accumulation. Product membership cannot change after construction; each leaf's caller-owned numerical array remains mutable. |
| Normalized test source | A source carrying `NormalizedTestSource()` radiometry and dimensionless relative power. It is useful for calibration and deterministic tests, but does not claim physical photon irradiance and requires an explicit scale before physical detector acquisition. |

The detector boundary distinguishes photon-rate density from cell-integrated
photon rate. Optical formation applies no elapsed-time factor. A prepared
detector acquisition validates the map metadata, applies the presampling
response before physical-pixel integration, and integrates its explicit
exposure duration exactly once. Incompatible maps remain separate in an
`OpticalProductBundle` unless an explicit conversion is prepared.

The definitions of [irradiance](https://cie.co.at/eilvterm/17-21-053),
[photon irradiance](https://cie.co.at/eilvterm/17-21-058),
[photon flux](https://cie.co.at/eilvterm/17-21-040), and
[photon exposure](https://cie.co.at/eilvterm/17-21-073) follow the CIE
International Lighting Vocabulary.

## Detector Response

| Term | Project meaning |
|---|---|
| Presampling detector response | A spatial-domain response applied before discrete detector sampling or binning. `GaussianPixelResponse`, `SampledFrameResponse`, and `RectangularPixelAperture` are response models. |
| Detector acquisition plan | `DetectorAcquisitionPlan`: a cold-path contract prepared for one frame detector and one immutable `IntensityMap` description. It validates plane, radiometry, current finite nonnegative samples, spatial measure, declared-channel QE, storage backend, and physical device, and sizes work buffers before repeated capture. Repeated capture trusts later producer writes to retain the sample-value contract. |
| Integration duration | The positive duration in seconds over which one detector-facing photon-arrival-rate sample contributes to an exposure. The current incremental convenience API names this `integration_duration`; it is neither an absolute sample timestamp nor the sample period between consecutive samples. A virtual-time scheduler owns the interval boundaries. |
| Nondestructive read | A detector read that observes accumulated charge without ending or resetting the active integration. A scheduled up-the-ramp model records the charge state at each declared read event; a post-exposure cube synthesized from one final frame is a lower-fidelity convenience, not a time-resolved nondestructive-read simulation. |
| Pixel aperture | The photosensitive spatial support and fill factor of one pixel. `RectangularPixelAperture` is separable in the implemented rectangular case, but its detector-grid kernel does not claim continuous subpixel geometry; that requires an explicit oversampled optical-grid mapping. |
| Optical transfer function (OTF) | The complex spatial-frequency response associated with an imaging response. |
| Modulation transfer function (MTF) | The magnitude of the OTF. `detector_mtf` is the normalized interior, infinite-grid transfer magnitude of the realized discrete response kernel, in cycles per detector pixel. A finite detector uses zero extension and therefore has boundary-dependent response rather than one global MTF. The diagnostic does not claim the continuous subpixel-aperture MTF; that requires an explicit oversampled optical-grid mapping. No spatial kernel type is itself named an MTF. |
| Quantum efficiency (QE) | Expected collected signal carriers per incident photon under the declared spectral model. |
| PRNU / DSNU | Pixel-response nonuniformity / dark-signal nonuniformity. |
| Interpixel capacitance (IPC) | Post-collection capacitive charge coupling between detector nodes, not a presampling optical blur. |

CCD, EMCCD, CMOS, sCMOS, quantitative low-noise CMOS, HgCdTe avalanche arrays,
Skipper CCD, SPAD, MKID, and APD describe detector technology or readout
families. Configured CMOS variants share the CMOS architecture in core; named
commercial camera profiles belong in a companion package. A single-pixel APD is
a channel detector and is not forced into an area-frame API.

## WFS And HIL Terms

| Term | Project meaning |
|---|---|
| Stable declared identity | An explicit typed component identity that remains attached to the same physical or logical owner when declarations are reordered. `OpticalPathID` and `AcquisitionID` wrap symbols for the implemented plant topology. Tuple position, task/thread identity, ring position, and device ordinal are not stable declared identities. |
| Plant definition | `PlantDefinition`: the structurally immutable binding of one telescope and atmosphere to reusable optical-path definitions and independent acquisition definitions. It is a topology declaration, not prepared execution, and owns no schedule, RNG stream, queue, transport, or HIL descriptor. The separately owned telescope and atmosphere retain their documented state semantics. |
| Cold plant-model definition | A configuration-only value whose type explicitly returns `ColdPlantModelDefinition()` from `plant_model_definition_style`. The opt-in contract excludes prepared plans and workspaces, mutable simulation or detector state, schedules, RNG streams, queues, transport, and HIL descriptors. Unrecognized types fail closed. |
| Preparation | The cold phase that validates configuration and binds it to concrete storage, array backend, physical device, plans, workspaces, and exact ownership relationships before repeated execution. A `prepare_*` function may allocate and plan; its corresponding warmed `!` operation reuses the prepared state. Preparation is distinct from simulation initialization that evolves physical state. |
| Optical-path definition | `OpticalPathDefinition`: an explicitly identified source and opted-in cold optical-model declaration that can be prepared once and reused by compatible acquisitions. It owns no propagation workspace or acquisition state. |
| Acquisition definition | `AcquisitionDefinition`: an explicitly identified, opted-in cold acquisition-model declaration plus a reference to one `OpticalPathID`. It is independent of the path definition and carries no schedule or prepared mutable acquisition state. |
| Prepared plant | `PreparedPlant`: the schedule-free result of preparing a `PlantDefinition`. It contains concrete tuples of path executors and acquisition owners plus exact stateful RNG streams derived from one run seed, derivation version, and stable owner identities. It does not advance the atmosphere, decide due work, own a clock, or publish through a HIL port. |
| Prepared path executor | `PreparedPathExecutor`: the single-writer owner of one frozen source, exact plant-atmosphere binding, path-local pupil-function, declared-plane field, or intensity input, prepared input-materialization operation, acquisition-facing photon-rate result, and concrete prepared optical execution state. Its result may be borrowed read-only by several acquisitions. |
| Calibration illumination | Illumination used in a calibration scenario. It is an application role, not a privileged source technology, optical propagation mode, detector mode, or control authority. It enters an ordinary prepared path at an explicitly declared boundary. |
| Illumination entry boundary | The typed location/product contract at which a prepared illumination evaluator writes a caller-owned `PupilFunction`, `ElectricField`, or `IntensityMap`. The boundary excludes upstream optics by construction; it does not request that core bypass named components. |
| Prepared illumination entry | `PreparedIlluminationEntry`: an internal schedule-free path materializer binding one entry tag, exact caller-owned payload, immutable evaluator wrapper, separate single-writer state/workspace, defensive visibility and payload contracts, explicit combination semantics, plant time, and a path-owned `:illumination` RNG stream. It is not an acquisition, scheduler, transport, or universal optical graph. |
| Downstream visibility | An application-owned description of which declared optical branches may consume a product or see an optic. The Gate 2 illumination seam snapshots this value but does not interpret topology; later placed-optics preparation validates physical visibility. |
| Acquisition product contract | `AcquisitionProductContract`: the run-immutable compatibility snapshot for one caller-owned `AcquisitionProducts` destination. It covers observation and measurement shape, numeric type, backend/device, typed units and metadata, and required acquisition-level metadata. That metadata carries any geometry, radiometry, layout, or semantics not already present in the typed products. A provider mutates and returns the exact destination; the contract is not a HIL descriptor, lease, or transport schema. |
| Acquisition product provider | A prepared per-acquisition implementation selected as full optical, command-responsive reduced order, or nonresponsive synthetic/replay through `acquisition_provider_style`. Every provider writes the same logical acquisition product contract. Provider selection is run-immutable; a different fidelity requires another preparation. Its `acquisition_provider_payload_work` declaration describes payload work but is not latency, throughput, cache-residency, or optical evidence. |
| Prepared acquisition owner | `PreparedAcquisitionOwner`: an independent single-writer acquisition endpoint bound to one exact prepared path result and one run-immutable `PreparedAcquisitionProvider`. It owns caller-visible `AcquisitionProducts` plus any provider, detector/readout, or WFS-estimator state and receives exact plant-owned RNG streams for execution. It does not own the path workspace, schedule, trigger, queue, or transport. |
| Prepared acquisition selection | `PreparedAcquisitionSelection`: a cold, schedule-free selection of exact acquisition owners, only the deduplicated full-optical paths their provider styles require, and references to the corresponding plant-owned RNG groups, all canonicalized by stable declared identity so declaration and caller-selection order cannot alter direct serial order. It owns no independent RNG state, cadence, trigger state, queue, transport, or retained atmosphere snapshot. |
| Path-result key | `PathResultKey`: the preparation-time value contract for compatible optical-result reuse, covering source geometry, spectral sampling, radiometry, optical/propagation model keys, sampling semantics, output plane, revisions, backend, and physical device. Its descriptive values are defensively snapshotted. Value equality and Julia's `hash` operation support cold keyed lookup; hashing is not performed by warmed path execution. `InstantaneousOpticalSample` denotes a rate sample at one plant instant, not an exposure or cadence. |
| WFS optical front end | Physical propagation and optics that consume an explicit pupil function or electric field and produce one or more detector-facing photon-arrival-rate products. Prepared implementations use `prepare_wfs_optical_formation` and `form_wfs_optical_products!`. |
| WFS observation | Acquired detector scalar, vector, frame, stack, packed product, or concrete tuple of products before estimation. `WFSObservation` binds caller-owned storage to explicit units, layout, backend, and physical-device metadata. |
| WFS measurement | Typed estimator output represented by `WFSMeasurement`, with declared semantic kind and units independently of storage rank. |
| WFS signal | Generic WFS measurement supplied to calibration or reconstruction. |
| Direct WFS measurement | An intentionally declared geometric or reduced-order path whose estimator consumes an explicit pupil function or electric field and produces a measurement without fictitious optical-rate or detector-observation products. Its prepared estimator returns `DirectMeasurementPath()`. |
| LiFT forward model | A prepared monochromatic focal-plane model that freezes the pupil transmission, dimensionless modal OPD shapes, diversity OPD in metres, source photon irradiance, optional object kernel, observation-grid mapping, backend, and device. LiFT coefficients are OPD metres, so a coefficient times its dimensionless mode shape contributes an OPD in metres. `PreparedLiFTForwardModel` owns reusable single-writer propagation state and predicts a cell-integrated photon-arrival rate; it owns no detector or acquisition cadence. |
| LiFT observation | Independently acquired caller-owned values bound by `LiFTObservation` to the prepared forward model's geometry, wavelength, spatial measure, deterministic preprocessing signature, backend, and physical device. The observation may declare photon-rate, expected-count, or normalized-intensity values; it is not produced or mutated by estimation. |
| LiFT observation domain | The explicit conversion between native observation values and the estimator's canonical incident-photon-rate domain. `LiFTPhotonRate` leaves values as rates, `LiFTExpectedCounts` declares exposure and QE, and `LiFTNormalizedIntensity` declares the photon rate represented by one normalized unit. Exposure and QE enter this conversion and noise weighting, not the focal-plane optical model. |
| LiFT frame mapping | `LiFTFrameMapping`: deterministic spatial preprocessing shared by forward preparation and external acquisition, ordered as presampling response followed by cell-summing sampling and binning. It excludes QE, exposure, gain, noise, readout windowing, and cadence. |
| LiFT estimator | `LiFT`: an iterative analytic- or numerical-Jacobian phase-retrieval estimator over one prepared forward model and independently supplied compatible observation. The estimated modal subset is fixed during estimator preparation. LiFT is not an `AbstractWFS` and does not trigger acquisition. |
| WFS slopes | Gradient or centroid-equivalent slope signals where the sensor and calibration establish that meaning. The existing generic `slopes` accessor is transitional and must not be extended to new non-slope sensors. |
| WFS flux normalization | Historical AO shorthand for normalization by summed detector signal or incident expected photon count. `MeanValidFluxNormalization` uses the mean measured valid-sample signal. `IncidenceFluxNormalization` uses expected incident signal per pupil sample; a detector-coupled measurement converts that denominator into the detector's deterministic signal units using exposure and detection efficiency, without applying stochastic detector effects to the calibration reference. Zero incident or detectable light returns a finite zero WFS signal. Existing `*FluxNormalization` names are transitional; they do not establish SI photon-flux units and should not be copied into new interfaces. |
| Microlens array | An independent optical element that partitions and focuses a field. The current `MicrolensArray` is a regular square angular-coordinate model: the subaperture layout supplies the entrance-pupil-referred subaperture pitch, while its parameters select the lenslet count and numerical focal-plane sampling. That pitch is not a claim about the manufactured lenslet pitch without an explicit relay magnification. The model does not yet include focal length, fill factor, per-lenslet prescription, or manufacturing errors. A Shack–Hartmann WFS composes it with acquisition and spot estimation. |
| Shack–Hartmann front end | The explicit component stored in `ShackHartmannWFS.front_end`. `ShackHartmannDirectFrontEnd` contains microlens geometry and subaperture layout for the intentional direct-measurement path and owns no propagation workspace. `ShackHartmannOpticalFrontEnd` additionally owns one prepared microlens propagation workspace for diffractive formation. Neither term denotes detector acquisition or centroid estimation. |
| Prepared microlens propagation | Backend- and sampling-bound FFT plans plus reusable optical scratch for executing a microlens array. It is execution state, not the microlens optic, detector acquisition, calibration, or estimator. |
| Focal-plane modulation | A prepared optical quadrature over focal-plane tip/tilt offsets in units of λ/D. `NoModulation`, `CircularModulation`, and `SampledModulation` describe the path; normalized weights average instantaneous intensities over a modulation cycle and do not represent exposure duration, trigger phase, or an RTC waveform. |
| Pyramid phase mask | The physical focal-plane phase-ramp optic of a Pyramid WFS, represented by `PyramidPhaseMask`. It is distinct from modulation and from the re-imaged-pupil differential estimator. |
| BioEdge amplitude mask | The physical family of complementary focal-plane amplitude filters of a BioEdge WFS, represented by `BioEdgeAmplitudeMask`. It is not a Pyramid phase mask even when both families reuse prepared modulation and propagation helpers. |
| Four-pupil mosaic | A detector-plane arrangement of four re-imaged pupil intensities used by Pyramid and BioEdge differential estimation. Its photon-arrival-rate values remain optical products until a detector acquisition applies response, QE, and exposure. |
| Curvature branch rate planes | The ordered positive- then negative-defocus detector-facing photon-arrival-rate products formed by a Curvature optical front end. They remain separate optical products until independent detectors or one explicit packed mapping acquire them. |
| Packed Curvature observation | One detector observation containing both compatible Curvature branches as spatial regions (`:curvature_branch_regions`) or counting channels (`:curvature_branch_channels`). Both branches share that detector's exposure; unequal exposures require separate detector acquisitions. |
| Subaperture layout revision | A monotonically advanced cold-configuration revision for synchronized valid-mask storage. Prepared Shack–Hartmann optical and estimator plans bind it and must be prepared again after a maintained layout update. Direct mask or host-mirror mutation is outside that contract. |
| Shack–Hartmann measurement order | A slope or centroid vector is `[axis 1; axis 2]`. Each square lenslet block uses Julia column-major order and therefore reshapes directly to the `(i, j)` indexing of the subaperture mask. This is distinct from OOPAO's axis-block and NumPy row-major convention. |
| Centroid response | The finite nonzero calibration factor converting a reference-subtracted detector centroid into one reported Shack–Hartmann measurement unit. It does not by itself assert angular units. |
| Atmosphere model time | Explicit simulation time, in seconds, owned by the scenario or scheduler and passed to atmosphere evolution. It is not inferred from telescope sampling, detector cadence, or wall time. |
| Atmosphere step | The positive model-time duration in seconds advanced for one runtime sensing update, configured explicitly as `atmosphere_step`. It is independent of detector exposure and is not a telescope property. |
| Atmosphere epoch token | `AtmosphereEpoch`: an immutable identity for the currently published atmosphere state, including atmosphere identity, model time, and a monotonic publication sequence. The current implementation validates it against the atmosphere writer's current state; it is not retained layer storage and does not by itself permit rendering after that writer advances. |
| Atmosphere direction renderer | A prepared, path-local mapping from one frozen source direction and telescope sampling geometry to caller-owned atmospheric output for a compatible epoch. It does not advance atmosphere time or consume RNG. |
| Atmosphere layer identity | `AtmosphereLayerID`: the stable declared name of one stochastic layer in a multilayer atmosphere prepared for plant execution. It is distinct from layer tuple position, altitude, and physical-device placement. |
| Materialized atmospheric path product | Caller-owned OPD, electric field, or model-specific prepared input containing the atmosphere-dependent data needed by one optical path at one epoch. Once materialized, downstream path work does not read mutable atmosphere layers for that sample. |
| Retained atmosphere-state snapshot | Optional bounded, versioned storage sufficient for a compatible renderer to evaluate an older atmosphere epoch after the atmosphere writer advances. It is distinct from an epoch token and need not copy every layer when a model provides a smaller correct representation. |
| Run seed | The explicit non-negative integer seed shared by one prepared plant run. It is combined with a derivation version and stable RNG owner identity; it is not consumed sequentially according to declaration order. |
| RNG derivation version | `RNGDerivationVersion`: the recorded positive version domain used when deriving prepared owner streams. Changing it intentionally selects different streams even when the run seed and topology are unchanged. |
| RNG owner identity | `RNGOwnerIdentity`: a stable `(category, component, role)` identity for one stochastic state owner, unique within a prepared plant and unchanged across equivalent prepared runs. It derives from declared component identity rather than tuple position, task, thread, completion order, or device placement. Together with the run seed and derivation version it selects that owner's reproducible stream or future addressable random domain. |
| Trigger | An acquisition-start event delivered to a detector endpoint. It is distinct from an internal camera oscillator or the simulation execution clock. |
| Cadence | Intended recurrence of an operation or acquisition. Cadence does not imply common trigger phase. |
| Sample period | The elapsed time between corresponding instants of consecutive samples; for uniform sampling it is the reciprocal of sample rate. It is distinct from the integration duration represented by one detector-facing optical sample. |
| Clock source | The timebase used to timestamp or schedule events. Detector trigger skew is modeled separately from clock-source error. |
| HIL execution time | Host-side time at which simulation work runs. It is distinct from model time and reported detector timestamps. |
| Plant command | A core virtual-time request expressed in the prepared plant command schema after any external timestamp has been mapped to plant time. Core owns its validation, admission, application, held state, and model disposition. |
| Command submission descriptor | A HIL-boundary descriptor that correlates one transferred payload lease and external timing metadata with a prepared command endpoint. Enqueueing the descriptor transfers ownership but does not imply plant-command admission or application. |
| Command outcome | A HIL-boundary completion that returns the submission credit and reports the correlated core model disposition plus boundary timing and failure metadata. It is not sampled device feedback. |

Use NGS for a natural guide star, LGS for a laser guide star, MCAO for
multi-conjugate adaptive optics, and MOAO for multi-object adaptive optics.
Woofer/tweeter behavior is represented by independent controllable optics, often
conjugated to the same altitude, rather than by a forced composite device.
