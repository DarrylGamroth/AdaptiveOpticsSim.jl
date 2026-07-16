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
- Use *WFS signal* for a sensor-independent measurement vector. Use *slope* only
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
| Sample-centered | A geometric origin lies on a stored sample. |
| Inter-sample-centered | A geometric origin lies between stored samples. This is not necessarily an inter-pixel convention. |
| Pupil function | The complex pupil transmission in Fourier optics. `PupilFunction` stores its wavelength-independent amplitude and OPD ingredients; field formation realizes the wavelength-specific complex pupil. |
| Pupil reflectivity | The finite dimensionless intensity-throughput fraction in `[0, 1]` stored by `TelescopeAperture`. Optical field formation applies its square root as an amplitude transmission. It is not an amplitude coefficient itself. |
| Electric field | In this package, a scalar complex optical field envelope/phasor sampled on a declared plane. `ElectricField` is not an SI vector electromagnetic field in V/m unless its metadata explicitly establishes such a normalization. |
| Intensity | The computational quantity `abs2(E)`. Its physical units and normalization follow the field and metadata; the word alone does not assert radiometric irradiance. |
| Intensity map | A real sampled array of intensity values plus `OpticalPlaneMetadata`, represented by `IntensityMap`. The metadata, not the type name, states plane location and physical semantics. |
| Photon-arrival-rate intensity map | An `IntensityMap` with `PhotonRateNormalization()` and an explicit spatial measure. It represents a rate before detector exposure: either a spatial density or a rate integrated over each represented cell. |
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
| Optical product bundle | `OpticalProductBundle`: a typed collection that preserves products on incompatible grids or channels without implicit resampling, radiometric conversion, or accumulation. |
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
| WFS optical front end | Physical propagation and optics that produce one or more detector-facing optical products. |
| WFS observation | Acquired detector frame, frame bundle, or direct reduced-order observation before estimation. |
| WFS signal | Generic estimator output supplied to calibration or reconstruction. |
| WFS slopes | Gradient or centroid-equivalent slope signals where the sensor and calibration establish that meaning. The existing generic `slopes` accessor is transitional and must not be extended to new non-slope sensors. |
| WFS flux normalization | Historical AO shorthand for normalization by summed detector signal or incident expected photon count. `MeanValidFluxNormalization` uses the mean measured valid-sample signal. `IncidenceFluxNormalization` uses expected incident signal per pupil sample; a detector-coupled measurement converts that denominator into the detector's deterministic signal units using exposure and detection efficiency, without applying stochastic detector effects to the calibration reference. Zero incident or detectable light returns a finite zero WFS signal. Existing `*FluxNormalization` names are transitional; they do not establish SI photon-flux units and should not be copied into new interfaces. |
| Microlens array | An independent optical element that partitions and focuses a field. A Shack–Hartmann WFS composes it with acquisition and spot estimation. |
| Atmosphere model time | Explicit simulation time, in seconds, owned by the scenario or scheduler and passed to atmosphere evolution. It is not inferred from telescope sampling, detector cadence, or wall time. |
| Atmosphere step | The positive model-time duration in seconds advanced for one runtime sensing update, configured explicitly as `atmosphere_step`. It is independent of detector exposure and is not a telescope property. |
| Atmosphere epoch | An immutable identity for one published atmosphere state, including atmosphere identity, model time, and a monotonic publication sequence. It is not an independently retained copy of every layer array. |
| Atmosphere direction renderer | A prepared, path-local mapping from one frozen source direction and telescope sampling geometry to caller-owned atmospheric output for a compatible epoch. It does not advance atmosphere time or consume RNG. |
| Trigger | An acquisition-start event delivered to a detector endpoint. It is distinct from an internal camera oscillator or the simulation execution clock. |
| Cadence | Intended recurrence of an operation or acquisition. Cadence does not imply common trigger phase. |
| Clock source | The timebase used to timestamp or schedule events. Detector trigger skew is modeled separately from clock-source error. |
| HIL execution time | Host-side time at which simulation work runs. It is distinct from model time and reported detector timestamps. |

Use NGS for a natural guide star, LGS for a laser guide star, MCAO for
multi-conjugate adaptive optics, and MOAO for multi-object adaptive optics.
Woofer/tweeter behavior is represented by independent controllable optics, often
conjugated to the same altitude, rather than by a forced composite device.
