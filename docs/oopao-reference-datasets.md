# OOPAO Reference Datasets

Status: active

This document describes the deterministic OOPAO bundles used to cross-validate
AdaptiveOpticsSim.jl.

For the separate SPECULA-targeted frozen bundle, see
[specula-reference-datasets.md](./specula-reference-datasets.md).

## Principles
- Use fixed seeds and deterministic settings.
- Store small, representative datasets.
- Compare with tolerances, not bitwise equality.
- Validate user-visible workflows and scientifically relevant outputs; do not
  require internal implementation details to match OOPAO exactly.

## Current committed bundle

The repository now commits a deterministic OOPAO bundle under
`test/reference_data/` covering:

1. `psf_baseline`
2. `shack_hartmann_geometric_ramp_xy`
3. `shack_hartmann_geometric_ramp_y`
4. `shack_hartmann_diffractive_ramp`
5. `shack_hartmann_diffractive_tip_mode`
6. `shack_hartmann_diffractive_tilt_mode`
7. `shack_hartmann_diffractive_tiptilt_dm`
8. `pyramid_diffractive_ramp`
9. `pyramid_diffractive_tip_mode`
10. `pyramid_diffractive_tiptilt_dm`
11. `bioedge_diffractive_ramp`
12. `bioedge_diffractive_tip_mode`
13. `bioedge_diffractive_tiptilt_dm`
14. `gain_sensing_camera_optical_gains`
15. `transfer_function_rejection`
16. `lift_interaction_matrix`
17. `closed_loop_shack_hartmann_trace`
18. `closed_loop_pyramid_trace`
19. `closed_loop_bioedge_trace`
20. `gsc_closed_loop_trace`
21. `gsc_atmosphere_replay_trace_bounded`
22. `gsc_branch_step_modulation_frame`
23. `gsc_branch_step_optical_gains`
24. `gsc_branch_step_signal`

These cases are stable enough to keep in CI and now cover the main image
formation, diffractive WFS, narrow modal-optic tip/tilt responses and one
maintained `tiptilt + dm` composite plant on Shack-Hartmann, Pyramid, and
BioEdge, LiFT Jacobians, compact closed-loop traces, GSC optical-gain,
bounded atmosphere-replay, gain-updated closed-loop behavior, the first
nonlinear GSC branch step, and analytical transfer-function paths.

The controllable-optic claim is intentionally narrow:
- committed OOPAO parity now exists for `TipTiltMirror`-equivalent Cartesian
  tip/tilt modes on diffractive `ShackHartmann`, `Pyramid`, and `BioEdge`
- committed OOPAO parity also now exists for one representative composite plant:
  `tiptilt + dm` on diffractive `ShackHartmann`, `Pyramid`, and `BioEdge`
- broader composite families such as `focus + dm` or richer grouped runtime
  surfaces are still validated through internal artifacts and backend parity
  rather than external OOPAO equivalence

The same `test/reference_data/` bundle also commits deterministic pyTomoAO
tomography references for:

1. `tomography_model_gamma`
2. `tomography_model_cxx`
3. `tomography_model_cox`
4. `tomography_model_cnz`
5. `tomography_model_reconstructor`
6. `tomography_model_wavefront`
7. `tomography_im_reconstructor`
8. `tomography_im_wavefront`
9. `tomography_kapa_model_wavefront`
10. `tomography_kapa_model_dm_commands`
11. `tomography_kapa_model_dm_commands_masked`

Diagnostic-only, not committed as a parity gate:
- Full long-horizon atmosphere-replay closed-loop tutorial traces. Local audits
  now show that the delayed-control path and residual-RMS telemetry track OOPAO
  closely, but Strehl and later slope norms still drift once both codes enter a
  huge-OPD nonlinear regime. Row-4 focal-plane/GSC diagnostics also show that
  very small control-vector differences can push the gain-sensing frame onto a
  different nonlinear branch. The committed `gsc_branch_step_*` cases isolate
  that first branch point and regression-test the residual OPD, Pyramid signal,
  modulation frame, and optical gains at fixed inputs.

## Reference configuration
- Seed: fixed RNG seed (documented per dataset).
- Noise: disabled unless explicitly validating noise models.
- Threads: single-threaded for determinism.

## Storage format
- Prefer a lightweight, format-agnostic bundle:
  - Arrays as `.npy` or `.npz` (easy to generate from Python).
  - Metadata as JSON or TOML (units, seeds, tolerances).
- Optionally mirror in FITS for astronomy workflows, but keep optional.

## Harness contract
AdaptiveOpticsSim.jl now includes a test-side reference harness in
`test/reference_harness.jl`.

Current expectations:
- The harness looks for a bundle root in `ENV["ADAPTIVEOPTICS_REFERENCE_ROOT"]`.
- If the environment variable is unset, it falls back to `test/reference_data/`.
- The bundle must contain a `manifest.toml`.
- Each case in the manifest records:
  - `kind` (`psf`, `shack_hartmann_slopes`, `pyramid_slopes`, `bioedge_slopes`)
  - `data` path
  - `shape`
  - optional `storage_order` (`F` / `julia_column_major` by default,
    `C` / `numpy_row_major` for NumPy row-major arrays)
  - `atol` / `rtol`
  - nested config tables (`telescope`, `source`, `opd`, `wfs`, `compute`)
  - optional `controllable_optic` table for cases that should be reconstructed
    through the maintained modal/control-surface API instead of direct OPD
    injection
  - optional `compare` rules for known convention adapters

This lets the Julia test suite reconstruct the scenario, compute the local
result, and compare it to the OOPAO-generated reference array.

Example manifest sketch:
```toml
version = 1

[cases.psf_baseline]
kind = "psf"
data = "psf_baseline.txt"
shape = [64, 64]
storage_order = "C"
atol = 1e-8
rtol = 1e-8

[cases.psf_baseline.telescope]
resolution = 32
diameter = 8.0
sampling_time = 1e-3
central_obstruction = 0.2

[cases.psf_baseline.source]
kind = "ngs"
band = "I"
magnitude = 0.0

[cases.psf_baseline.compute]
zero_padding = 2
```

For OOPAO geometric Shack-Hartmann data, the committed bundle also uses:

```toml
[cases.shack_hartmann_geometric_ramp_xy.compare]
convention = "oopao_geometric_sh_signal_2d"
```

That adapter is intentional. OOPAO’s exported `signal_2D` ordering and
calibrated slope units do not match the raw OPD-gradient convention currently
returned by AdaptiveOpticsSim.jl. The reference harness now treats this as an
explicit named convention rather than as anonymous `swap_halves` / `scale`
flags in the manifest.

For now the harness uses plain text arrays (`DelimitedFiles`) because that keeps
the validation path dependency-light. The loader now treats image-storage
ordering as an explicit convention (`julia_column_major` or
`numpy_row_major`), while still accepting the shorter legacy `F` / `C`
spellings. If `.npz` becomes more convenient, we can add a test-only reader
later without changing the manifest structure.

LiFT reference cases also treat detector-coupled image sizing explicitly: if a
case does not provide `compute.img_resolution`, the harness uses
`detector.psf_sampling * telescope.resolution`, which matches the practical
detector-pixel convention used by the current OOPAO and Julia LiFT surfaces.

Modal reference cases can also carry both:
- an `opd`/`basis` description for the OOPAO-side provenance
- a `controllable_optic` description for the Julia-side maintained API path

That lets the same frozen case compare OOPAO’s direct injected optical phase to
the Julia package’s modal controllable-optic execution path without pretending
the two codes assemble the surface internally in the same way.

## Generation workflow

The OOPAO generator lives in `scripts/generate_oopao_reference_bundle.py`.
The pyTomoAO generator lives in `scripts/generate_pytomoao_reference_bundle.py`.

Recommended approach:
1. Build a reproducible Python 3.8 image with the OOPAO dependency stack.
2. Run the generator with a pinned OOPAO upstream repo/ref. The script will
   clone OOPAO into a writable temporary directory before import. This matters
   because OOPAO writes `precision_oopao.npy` into its package root at import
   time.
3. Point the generator at a bundle directory.

The committed geometric SH bundle was generated this way and then copied into
`test/reference_data/`.

Recommended CI/local command:

```bash
python3 scripts/generate_oopao_reference_bundle.py /tmp/oopao-bundle \
  --oopao-repo https://github.com/cheritier/OOPAO.git \
  --oopao-ref 085d5e50ace0d20fe13cc2da20129d5400166973
```

If you already have a local writable OOPAO checkout and want to use that
instead, pass `--oopao-path /path/to/OOPAO`.

The generated `manifest.toml` now records:
- `metadata.oopao.repo_url`
- `metadata.oopao.requested_ref`
- `metadata.oopao.resolved_commit`

That provenance is intended to make CI artifacts and reference-bundle refreshes
auditable across machines.

## Tolerances
- Record per-dataset tolerances (relative and absolute).
- Note expected deviations for CPU vs GPU and `Float32` vs `Float64`.

For the current committed SH cases, tolerances are effectively exact once the
documented convention adapter is applied.
