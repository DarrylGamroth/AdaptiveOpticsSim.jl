# Core Optics Expansion Roadmap

This roadmap defines the next ranked core-optics work after the atmosphere and
runtime repairs. It is intentionally focused on the simulation core rather than
`Proper.jl` integration or richer controller families.

The ranking is based on:

1. scientific leverage on the AO plant and sensing stack,
2. how much new work it unlocks downstream,
3. how directly it closes gaps relative to SPECULA and the stronger parts of
   OOPAO,
4. how naturally it fits the current Julia module structure.

## Current Baseline

AdaptiveOpticsSim.jl already has:

- maintained telescope/source/asterism/DM/NCPA primitives,
- maintained finite and infinite atmosphere backends on CPU and GPU,
- maintained SH, Pyramid, BioEdge, Zernike, and Curvature WFS families,
- maintained tomography, LiFT, detector, and HIL-style runtime surfaces.

What it does not yet have in a first-class way is:

- an explicit electric-field representation,
- Fresnel/chromatic field propagation as a maintained core surface,
- polychromatic sensor propagation,
- generic extended-source sensor optics,
- a unified aperture/mask primitive layer,
- a maintained subaperture-calibration surface for SH beyond the current
  in-sensor logic.

That makes the next useful roadmap much closer to SPECULA's field/propagation
surface than to additional OOPAO parity work.

## Ground Rules

- Keep core optics in-core; keep coronagraphs and science-prescription logic out
  of core.
- Start each milestone with a CPU reference implementation and backend-generic
  data layout.
- Port hot paths to GPU only after the numerical contract is fixed and tested.
- Keep params/state separation explicit.
- Use preallocated workspaces and `!` mutating methods on hot paths.
- Preserve deterministic replay and stable RNG ownership.
- Reuse the existing interface-contract style in `docs/api-reference.md`.

## Ranked Execution Order

1. CO-1: Electric-field core
2. CO-2: Fresnel and chromatic propagation
3. CO-3: Polychromatic wavefront-sensor execution
4. CO-4: Extended-source sensing optics
5. CO-5: Aperture and mask primitives
6. CO-6: Shack-Hartmann subaperture calibration and slope extraction

## CO-1: Electric-Field Core

Goal: add a maintained complex-field representation so propagation and sensing
are not forced to pass only through OPD maps.

Why first:

- This is the largest remaining core-optics gap relative to SPECULA.
- Fresnel propagation, polychromatic sensing, extended-source sensing, and any
  later science-path bridge all depend on it.
- It can be introduced without destabilizing the existing telescope/OPD API if
  the conversion boundary is kept explicit.

Deliverables:

- Add a maintained electric-field type family:
  - `ElectricFieldParams`
  - `ElectricFieldState`
  - `ElectricField`
- Represent:
  - wavelength,
  - pupil sampling,
  - complex amplitude,
  - optional cached intensity/phase views,
  - reusable FFT work buffers/plans where needed.
- Add explicit builders and mutators:
  - `ElectricField(tel, src; backend=Array, T=Float64, pad=1)`
  - `fill_from_telescope!(field, tel, src)`
  - `apply_phase!(field, opd_or_phase; units=:opd)`
  - `apply_amplitude!(field, amp)`
  - `intensity!(out, field)`
- Keep telescope OPD as the current maintained plant state; the field is a
  derived propagation object, not a replacement for `TelescopeState`.

Primary file targets:

- New: `src/Optics/electric_field.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Optics/telescope.jl`
- Update: `src/Optics/source.jl`
- Update: `src/Optics/psf.jl`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Introduce the new type family and export surface in
   `src/Optics/electric_field.jl`.
2. Factor the existing PSF input preparation into reusable field-construction
   helpers instead of duplicating telescope-to-field code.
3. Keep units explicit at the boundary:
   - OPD in meters in telescope state,
   - phase in radians only inside field-application helpers when requested.
4. Add CPU reference tests for:
   - zero-OPD field build,
   - OPD-to-phase application,
   - amplitude masking,
   - Fraunhofer-limit PSF equivalence with the current `compute_psf!` path.
5. Add backend-smoke coverage once the CPU contract is stable.

Acceptance gate:

- `ElectricField` can reproduce the existing monochromatic PSF path within
  tolerance.
- Field construction and simple phase/amplitude application are deterministic.
- The new path does not force allocations into the current telescope hot path.

## CO-2: Fresnel And Chromatic Propagation

Goal: add maintained field propagation models beyond the current Fraunhofer-like
PSF path.

Why second:

- It is the first major consumer of the electric-field core.
- It improves the realism of curvature sensing, finite-height propagation, and
  future science handoff without pulling in `Proper.jl`.

Deliverables:

- Add a propagation model interface:
  - `AbstractPropagationModel`
  - `FraunhoferPropagation`
  - `FresnelPropagation`
- Add reusable propagation entry points:
  - `propagate_field!(out, field, model)`
  - `propagate_field!(field, model)`
- Support:
  - monochromatic pupil-to-focal propagation,
  - Fresnel propagation over explicit distance,
  - wavelength-aware propagation metadata.
- Reuse the same field/workspace surface on CPU and GPU.

Primary file targets:

- New: `src/Optics/propagation.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Optics/electric_field.jl`
- Update: `src/Optics/psf.jl`
- Update: `src/WFS/curvature.jl`
- Update: `docs/api-reference.md`
- Update: `docs/benchmark-matrix-plan.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Extract the current pupil-to-PSF FFT path in `src/Optics/psf.jl` into a
   reusable Fraunhofer propagation helper.
2. Add a Fresnel transfer-function path with precomputed phase factors and
   explicit propagation distance.
3. Introduce a small workspace layer so propagation plans are reused rather than
   rebuilt per call.
4. Wire the curvature WFS path to opt into the maintained propagation API
   instead of owning bespoke diffraction-only code where that improves clarity.
5. Add benchmark and regression coverage for:
   - Fraunhofer equivalence to the current PSF path,
   - Fresnel round-trip sanity cases,
   - CPU/GPU parity on a small maintained case.

Acceptance gate:

- Fraunhofer propagation matches the existing PSF behavior within tolerance.
- Fresnel propagation is deterministic and backend-safe.
- Curvature-path coupling can use the maintained propagation API without a
  performance regression that breaks the current runtime targets.

## CO-3: Polychromatic Wavefront-Sensor Execution

Goal: support maintained multi-wavelength sensing paths instead of forcing WFS
algorithms to be effectively monochromatic.

Why third:

- This extends the field and propagation work directly into the sensing stack.
- SPECULA has a stronger surface here than the current project.
- It materially improves SH and Pyramid realism without dragging in science-only
  optics.

Deliverables:

- Add a compact spectral bundle type, for example:
  - `SpectralSample`
  - `SpectralBundle`
- Add source-side helpers for weighted wavelength samples.
- Add maintained polychromatic execution for:
  - `ShackHartmann`
  - `PyramidWFS`
- Keep the single-wavelength case allocation-free and on the current fast path.

Primary file targets:

- New: `src/Optics/spectrum.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Optics/source.jl`
- Update: `src/Optics/electric_field.jl`
- Update: `src/Optics/propagation.jl`
- Update: `src/WFS/shack_hartmann.jl`
- Update: `src/WFS/pyramid.jl`
- Update: `src/Control/runtime.jl`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Add source-side spectral sampling types that reduce to the existing
   monochromatic `Source` behavior when only one wavelength is present.
2. Build a grouped wavelength-accumulation path for SH and Pyramid that reuses
   prepared workspaces rather than allocating one field object per sample.
3. Preserve a strict fast path:
   - one wavelength keeps the current execution shape,
   - many wavelengths opt into grouped propagation explicitly.
4. Add regression tests for:
   - one-sample equivalence with current monochromatic behavior,
   - weighted broadband accumulation,
   - deterministic grouped-runtime behavior.
5. Add a maintained profile surface for broad-band SH and Pyramid cases.

Acceptance gate:

- The one-wavelength path reproduces current results exactly or within current
  tolerance.
- Broad-band SH and Pyramid runs are regression-backed and benchmarked.

## CO-4: Extended-Source Sensing Optics

Goal: support maintained source-distribution models for SH and Pyramid sensing.

Why fourth:

- Extended-source sensing becomes much cleaner once polychromatic field
  propagation exists.
- It addresses a practical realism gap in both natural and laser-guide-star
  workflows.

Deliverables:

- Add source-distribution models, for example:
  - `PointCloudSourceModel`
  - `GaussianDiskSourceModel`
  - `SampledImageSourceModel`
- Add maintained extended-source execution for:
  - `ShackHartmann` first,
  - `PyramidWFS` second.
- Reuse the same grouped accumulation machinery used for polychromatic
  propagation where possible.

Primary file targets:

- New: `src/Optics/extended_source.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Optics/source.jl`
- Update: `src/WFS/shack_hartmann.jl`
- Update: `src/WFS/pyramid.jl`
- Update: `src/Control/runtime.jl`
- Update: `examples/tutorials/`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Introduce source-distribution models separately from `Source` so point-source
   behavior remains simple and fast.
2. Start with SH because its spot-formation and slope extraction path is the
   easiest maintained target.
3. Implement grouped source-sample accumulation using the same runtime-prepared
   pattern as the polychromatic path.
4. Add Pyramid support after the SH path is numerically pinned down.
5. Cover:
   - point-source limit,
   - symmetric extended-source broadening,
   - deterministic replay for grouped source execution.

Acceptance gate:

- Point-source behavior is unchanged when no extended-source model is attached.
- SH extended-source optics are regression-backed.
- Pyramid extended-source support lands only after the SH surface is stable.

## CO-5: Aperture And Mask Primitives

Goal: unify mask and support-map construction across telescope, WFS, spatial
filter, and tomography surfaces.

Why fifth:

- The project already has many local mask builders, but they are spread across
  telescope, WFS, and calibration code.
- A maintained primitive layer will simplify later field propagation, source
  models, and calibration builders.

Deliverables:

- Add maintained mask/aperture builders, for example:
  - `CircularAperture`
  - `AnnularAperture`
  - `SpiderMask`
  - `RectangularROI`
  - `SubapertureGridMask`
- Add backend-safe builders for:
  - boolean masks,
  - weighted amplitude masks,
  - support/ROI maps.
- Replace ad hoc duplicated geometry assembly where it is clearly shared.

Primary file targets:

- New: `src/Optics/aperture_masks.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Optics/telescope.jl`
- Update: `src/Optics/spatial_filter.jl`
- Update: `src/WFS/shack_hartmann.jl`
- Update: `src/WFS/pyramid.jl`
- Update: `src/Tomography/reconstructors.jl`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Add a small mask-builder API in one place rather than continuing to expand
   local geometry helpers.
2. Move shared circular/annular/spider/support-map logic into the new module.
3. Convert the easiest consumers first:
   - telescope pupil/spiders,
   - spatial-filter focal masks where appropriate,
   - WFS valid/support masks.
4. Only then touch tomography mask assembly if the new primitives actually
   reduce duplication.
5. Add conformance tests for CPU/GPU mask builds and shape semantics.

Acceptance gate:

- Maintained mask builders replace clear duplication without slowing existing hot
  paths.
- Telescope and WFS geometry continue to match the current regression suite.

## CO-6: Shack-Hartmann Subaperture Calibration And Slope Extraction

Goal: make the SH optical/calibration surface more explicit and maintainable.

Why sixth:

- The current SH implementation is already strong, so this is a hardening and
  realism milestone rather than a missing primitive.
- SPECULA has a stronger dedicated surface for subaperture handling and
  calibration than the current project.

Deliverables:

- Add a maintained SH subaperture description layer, for example:
  - `SubapertureLayout`
  - `SubapertureCalibration`
- Separate:
  - lenslet/support geometry,
  - calibration/reference spot data,
  - slope extraction weighting/windowing,
  - detector-coupled output reshaping.
- Make the calibration surface reusable across geometric and diffractive SH
  modes where possible.

Primary file targets:

- New: `src/WFS/subapertures.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/WFS/shack_hartmann.jl`
- Update: `src/Calibration/interaction_matrix.jl`
- Update: `src/Control/runtime.jl`
- Update: `docs/api-reference.md`
- Update: `examples/tutorials/`
- Update: `test/runtests.jl`

Implementation plan:

1. Extract subaperture geometry and valid-support bookkeeping from the monolithic
   SH path into a maintained layout object.
2. Add explicit calibration/reference structures so geometric and diffractive
   SH share the same outer contract.
3. Move slope-extraction weighting/windowing into named maintained components
   rather than leaving it as purely internal numeric glue.
4. Add regression cases covering:
   - geometric SH,
   - diffractive SH,
   - detector-coupled SH,
   - grouped-runtime SH if affected.

Acceptance gate:

- SH calibration and slope extraction are easier to extend without duplicating
  mode-specific logic.
- Existing SH regression coverage stays green through the refactor.

## Recommended Delivery Plan

Stage A: foundation

1. CO-1
2. CO-2

Stage B: sensor realism

1. CO-3
2. CO-4

Stage C: geometry and calibration hardening

1. CO-5
2. CO-6

Recommended stop points:

- After CO-1: decide whether the field API is minimal and stable enough.
- After CO-2: re-benchmark curvature and PSF paths before expanding into WFS.
- After CO-4: compare the new sensing surface against SPECULA again.
- After CO-6: decide whether the remaining core-optics gaps justify `Proper.jl`
  or whether control/runtime work becomes the new bottleneck.

## Explicit Deferrals

This roadmap does not include:

- coronagraph models in core,
- broad controller-family work,
- transport/HIL integration layers,
- distributed scheduling,
- a full `Proper.jl` bridge.

Those should stay deferred until the ranked core-optics milestones above have
either landed or been explicitly rejected.
