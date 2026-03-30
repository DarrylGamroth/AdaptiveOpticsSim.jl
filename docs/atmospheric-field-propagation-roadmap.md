# Atmospheric Field Propagation Roadmap

This roadmap defines the next algorithmic target after the completed
core-optics expansion work. It is based on the current re-baseline against
`../OOPAO` and `../SPECULA`.

The conclusion of that re-baseline is:

- OOPAO is no longer the limiting reference for core optics.
- AdaptiveOpticsSim now has the main standalone field/propagation primitives:
  `ElectricField`, Fraunhofer/Fresnel propagation, polychromatic SH/Pyramid,
  extended-source sensing, mask primitives, and SH subaperture calibration.
- The largest remaining algorithmic gap relative to SPECULA is not another
  standalone optic, but a maintained atmosphere-aware field propagation layer
  like
  [specula/processing_objects/atmo_propagation.py](/home/dgamroth/workspaces/codex/SPECULA/specula/processing_objects/atmo_propagation.py).

That makes the next useful milestone sequence a new atmospheric
field-propagation surface rather than immediate `Proper.jl` integration or
broader controller work.

## Goal

Add a maintained atmosphere-aware electric-field propagation layer that combines
the existing atmosphere backends, source geometry, and field/propagation stack
into one reusable core-optics surface.

This surface should support:

- geometric phase-only propagation through atmospheric layers,
- physical Fresnel propagation between layers,
- monochromatic and wavelength-aware execution,
- both finite and infinite atmosphere backends,
- CPU and GPU targets as maintained execution paths, not optional follow-on
  ports,
- reuse by curvature sensing and future science-path handoff.

## Non-Goals

- Do not replace `TelescopeState.opd` as the primary AO plant state.
- Do not implement coronagraphs or `Proper.jl` integration here.
- Do not widen controller scope here.
- Do not duplicate the existing SH/Pyramid grouped-spectral machinery unless the
  atmospheric propagation layer needs a shared helper.

## Reference Strategy

- Use OOPAO as a behavior sanity reference where it has an equivalent path.
- Use SPECULA as the primary algorithm/decomposition reference for this work.
- Keep the Julia implementation idiomatic:
  - params/state split,
  - mutating hot paths,
  - backend-generic arrays,
  - deterministic replay and central RNG ownership.
- CPU may be the first correctness reference for a slice, but AFP-2 onward is
  not complete until the maintained AMDGPU and CUDA paths exist and pass the
  relevant smoke/regression coverage.

## Ranked Execution Order

1. AFP-1: Shared atmosphere-field propagation contract
2. AFP-2: Geometric layered field propagation
3. AFP-3: Fresnel layered field propagation
4. AFP-4: Chromatic atmospheric propagation
5. AFP-5: Runtime and sensing integration
6. AFP-6: Benchmarks, regressions, and backend hardening

## AFP-1: Shared Atmosphere-Field Propagation Contract

Goal: define the maintained public and internal surface that connects the
current atmosphere backends to `ElectricField`.

Why first:

- The project now has atmosphere evolution and field propagation, but no stable
  contract between them.
- This is the right point to normalize layer ordering, source-footprint
  geometry, and the boundary between OPD state and derived field state.
- It keeps later Fresnel/chromatic work from growing ad hoc paths in curvature
  or tutorials first.

Deliverables:

- Add a maintained propagation type family, for example:
  - `AtmosphericFieldPropagationParams`
  - `AtmosphericFieldPropagationState`
  - `AtmosphericFieldPropagation`
- Add explicit entry points, for example:
  - `AtmosphericFieldPropagation(atm, tel, src; model=:geometric, zero_padding=1, ...)`
  - `propagate_atmosphere_field!(field, prop, atm, tel, src)`
  - `propagate_atmosphere_field!(prop, atm, tel, src)`
- Define the layer-side internal contract needed by both finite and infinite
  atmosphere backends:
  - layer altitude,
  - source footprint geometry,
  - phase sampling into a provided buffer,
  - layer ordering for up/down propagation.
- Keep the field object derived and reusable:
  - build once,
  - refill/repropagate many times.

Primary file targets:

- New: `src/Optics/atmospheric_field_propagation.jl`
- Update: `src/AdaptiveOpticsSim.jl`
- Update: `src/Atmosphere/source_geometry.jl`
- Update: `src/Atmosphere/multilayer.jl`
- Update: `src/Atmosphere/infinite_screen.jl`
- Update: `src/Optics/electric_field.jl`
- Update: `docs/api-reference.md`
- Update: `docs/atmosphere-runtime-spec.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Introduce the new type family and export surface in
   `src/Optics/atmospheric_field_propagation.jl`.
2. Factor any missing shared layer helpers out of the finite/infinite backends
   so the propagation layer does not special-case one backend.
3. Normalize layer ordering and direction semantics:
   - source to telescope,
   - telescope to source if later needed,
   - no implicit reversal hidden in call sites.
4. Add CPU reference tests for:
   - vacuum/no-atmosphere no-op,
   - one-layer geometric phase application,
   - equivalence to the current OPD accumulation path in the purely geometric
     limit.

Acceptance gate:

- A field propagated through atmospheric layers in the geometric limit matches
  the current OPD-based path within tolerance.
- Finite and infinite atmosphere backends conform to one maintained propagation
  contract.
- The surface is deterministic and backend-ready before Fresnel work starts.

## AFP-2: Geometric Layered Field Propagation

Goal: implement the fast phase-only layered atmospheric field path on top of the
new contract.

Why second:

- This is the simplest physically meaningful coupled path.
- It provides immediate value for curvature and future field-based sensing
  without waiting for full Fresnel propagation.
- It gives a strong regression anchor before more expensive propagation is added.

Deliverables:

- Apply per-layer phase to the field in correct geometric order.
- Support:
  - finite and infinite atmospheres,
  - on-axis and off-axis sources,
  - finite-height sources and cone footprints,
  - CPU and GPU execution as first-class maintained paths.
- Reuse shared field buffers rather than building a new `ElectricField` per
  layer or per step.

Primary file targets:

- New: `src/Optics/atmospheric_field_propagation.jl`
- Update: `src/Atmosphere/source_geometry.jl`
- Update: `src/Atmosphere/multilayer.jl`
- Update: `src/Atmosphere/infinite_screen.jl`
- Update: `src/Optics/electric_field.jl`
- Update: `docs/api-reference.md`
- Update: `scripts/gpu_smoke_contract.jl`
- Update: `test/runtests.jl`

Implementation plan:

1. Add a reusable phase-sampling buffer sized to the active field aperture.
2. Sample each atmospheric layer into that buffer using the existing
   source-aware geometry helpers.
3. Apply the sampled phase to the field with a maintained mutating helper.
4. Add regressions for:
   - finite vs infinite atmosphere agreement over short maintained runs,
   - off-axis geometric propagation,
   - AMDGPU and CUDA parity on maintained contract cases.

Acceptance gate:

- Geometric atmospheric field propagation reproduces the current OPD-based
  optical state in the no-Fresnel limit.
- The one-source hot path stays allocation-controlled.
- AMDGPU and CUDA smoke coverage both pass, and GPU execution is not treated as
  a deferred follow-on.

## AFP-3: Fresnel Layered Field Propagation

Goal: add maintained physical propagation between atmospheric layers, not only
phase application at a single plane.

Why third:

- This is the main algorithmic gap relative to SPECULA.
- It lets the atmosphere participate directly in wave propagation rather than
  acting only as a collapsed OPD screen at the pupil.
- Curvature and future physical sensor paths benefit immediately.

Deliverables:

- Add a layered Fresnel execution mode that:
  - propagates the field between layer heights,
  - supports downwards and later upwards propagation,
  - reuses precomputed transfer functions/workspaces.
- Add band-limited angular-spectrum support where needed to avoid obvious
  aliasing/pathology at larger propagation distances.
- Keep geometric propagation as the explicit fast path.
- Keep GPU execution first-class here as well; layered Fresnel propagation is
  not considered complete on CPU alone.

Primary file targets:

- Update: `src/Optics/propagation.jl`
- Update: `src/Optics/atmospheric_field_propagation.jl`
- Update: `src/Atmosphere/multilayer.jl`
- Update: `src/Atmosphere/infinite_screen.jl`
- Update: `src/WFS/curvature.jl`
- Update: `docs/api-reference.md`
- Update: `docs/benchmark-matrix-plan.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Build per-layer propagation distances from the atmosphere/source/telescope
   geometry.
2. Reuse the current `FresnelPropagation` machinery where it is sufficient, and
   extend it only where atmosphere-coupled execution needs additional controls.
3. Add band-limit/padding controls at the propagation-model level instead of
   burying them in a one-off sensor path.
4. Add regressions for:
   - zero-distance equivalence to the geometric path,
   - deterministic multi-layer Fresnel replay,
   - maintained CPU/AMDGPU/CUDA parity cases.

Acceptance gate:

- Layered Fresnel propagation is deterministic and regression-backed.
- Curvature can use the maintained layer-aware field path.
- The geometric path remains available and faster for HIL-oriented cases.
- AMDGPU and CUDA both pass maintained Fresnel smoke coverage.

## AFP-4: Chromatic Atmospheric Propagation

Goal: extend the layered field path to wavelength-aware atmospheric propagation,
including chromatic anisoplanatic effects where justified.

Why fourth:

- The project already has `SpectralSource`; the missing piece is atmospheric
  wavelength-aware propagation rather than standalone spectral accumulation.
- This closes a remaining SPECULA-style capability gap without requiring an
  external science package.

Deliverables:

- Support grouped spectral atmospheric propagation on top of `SpectralBundle`.
- Keep the one-wavelength path on the current fast path.
- Add optional chromatic layer-shift/refraction support with explicit reference
  wavelength semantics.
- Preserve backend-generic grouped execution so spectral atmosphere propagation
  remains maintained on GPU rather than CPU-only glue around GPU optics.

Primary file targets:

- Update: `src/Optics/spectrum.jl`
- Update: `src/Optics/atmospheric_field_propagation.jl`
- Update: `src/Optics/propagation.jl`
- Update: `src/Atmosphere/source_geometry.jl`
- Update: `src/Control/runtime.jl`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Extend the propagation workspace so grouped spectral execution reuses field
   and transfer-function buffers.
2. Add explicit reference-wavelength metadata where chromatic layer shifts are
   modeled.
3. Reuse the existing spectral grouping model from SH/Pyramid where it fits, but
   keep the atmospheric propagation layer as the owner of wavelength-dependent
   atmospheric geometry.
4. Add regressions for:
   - one-sample equivalence to monochromatic propagation,
   - weighted broadband accumulation through atmosphere,
   - deterministic grouped-runtime behavior.

Acceptance gate:

- One-sample spectral runs match monochromatic runs within tolerance.
- Broad-band atmosphere propagation is regression-backed and benchmarked.
- Chromatic effects are explicit and optional, not hidden in the default path.
- The maintained spectral atmosphere path executes on CPU, AMDGPU, and CUDA.

## AFP-5: Runtime And Sensing Integration

Goal: connect the new atmospheric field propagation layer to the parts of the
sensing/runtime surface that actually need it.

Why fifth:

- The core algorithm should be proven before it is wired broadly into runtime.
- Curvature is the first obvious consumer, but not the only one.

Deliverables:

- Integrate prepared atmospheric field propagation into runtime preparation.
- Use it in `CurvatureWFS` as the maintained physical path.
- Add grouped-source support where the new layer makes missing diffractive
  execution practical, especially the current curvature asterism gap.

Primary file targets:

- Update: `src/Control/runtime.jl`
- Update: `src/WFS/curvature.jl`
- Update: `src/Optics/atmospheric_field_propagation.jl`
- Update: `examples/tutorials/`
- Update: `docs/api-reference.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Add runtime-preparation hooks for atmospheric field propagation workspaces.
2. Refactor `CurvatureWFS` to consume the maintained path rather than a bespoke
   propagation flow.
3. Add curvature asterism support if the new shared grouped-field path makes it
   straightforward.
4. Add tutorial and regression coverage for:
   - curvature with atmosphere-aware field propagation,
   - prepared-runtime replay,
   - grouped-source curvature if implemented in this milestone.

Acceptance gate:

- Curvature uses the maintained atmosphere-aware field path.
- Prepared runtime can reuse the new propagation workspace without per-step
  rebuilds.
- The current curvature asterism gap is either closed or explicitly documented
  as the remaining blocker.

## AFP-6: Benchmarks, Regressions, And Backend Hardening

Goal: make the new coupled propagation layer a maintained surface rather than a
fragile feature.

Why last:

- The atmospheric field path will cross several of the project's historically
  risky seams: atmosphere, FFT-based optics, runtime preparation, and GPU
  backends.

Deliverables:

- Add maintained CPU, AMDGPU, and CUDA smoke coverage.
- Add representative benchmark surfaces for:
  - geometric atmospheric field propagation,
  - layered Fresnel atmospheric propagation,
  - curvature with the maintained coupled path.
- Add direct regression coverage for:
  - finite vs infinite backend agreement in maintained cases,
  - deterministic replay,
  - no-atmosphere and zero-distance degeneracies.

Primary file targets:

- Update: `scripts/gpu_smoke_contract.jl`
- Update: `scripts/profile_atmosphere_runtime.jl`
- Update: `docs/benchmark-matrix-plan.md`
- Update: `docs/atmosphere-runtime-spec.md`
- Update: `test/runtests.jl`

Implementation plan:

1. Extend the current atmosphere and GPU contract suite to include the new field
   path.
2. Add representative profile scripts or cases so geometric and Fresnel coupled
   propagation can be tracked over time.
3. Record the maintained benchmark surfaces and acceptable tolerances in docs.

Acceptance gate:

- The new propagation layer is covered by smoke tests, regressions, and
  maintained benchmark surfaces on CPU, AMDGPU, and CUDA.

## Immediate Sprint Recommendation

Start with AFP-1 and AFP-2 together:

1. define the maintained contract and workspace types,
2. implement the geometric layered field path first,
3. prove equivalence to the current OPD-based geometric limit,
4. only then add layered Fresnel propagation.

That keeps the next slice scientifically meaningful without forcing the Fresnel
path onto an unproven interface.
