# OOPAO Reference Dataset Plan (Draft)

This document outlines which reference outputs to capture from OOPAO to validate
AdaptiveOptics.jl with deterministic comparisons.

## Principles
- Use fixed seeds and deterministic settings.
- Store small, representative datasets.
- Compare with tolerances, not bitwise equality.

## Candidate datasets
1) PSF baseline
   - Tutorial: `tutorials/image_formation.py`
   - Outputs: PSF array, cropped PSF, pixel scale, Strehl (if computed)
   - Notes: Use a single source and fixed zero-padding.

2) Zernike modes
   - Tutorial: `tutorials/image_formation.py`
   - Outputs: First N Zernike modes, cross-product matrix stats

3) Shack-Hartmann slopes
   - Tutorial: `tutorials/AO_closed_loop_ShackHartmann_WFS.py`
   - Outputs: slopes vector, valid subaperture mask, subaperture spot cube

4) Asterism PSF grid
   - Tutorial: `tutorials/how_to_asterism.py`
   - Outputs: PSF grid (per-source), combined PSF

5) Atmospheric phase screen
   - Tutorial: `tutorials/image_formation.py` or minimal script
   - Outputs: OPD screen at fixed time step(s)

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
AdaptiveOptics.jl now includes a test-side reference harness in
`test/reference_harness.jl`.

Current expectations:
- The harness looks for a bundle root in `ENV["ADAPTIVEOPTICS_REFERENCE_ROOT"]`.
- If the environment variable is unset, it falls back to `test/reference_data/`.
- The bundle must contain a `manifest.toml`.
- Each case in the manifest records:
  - `kind` (`psf`, `shack_hartmann_slopes`, `pyramid_slopes`, `bioedge_slopes`)
  - `data` path
  - `shape`
  - `atol` / `rtol`
  - nested config tables (`telescope`, `source`, `opd`, `wfs`, `compute`)

This lets the Julia test suite reconstruct the scenario, compute the local
result, and compare it to the OOPAO-generated reference array.

Example manifest sketch:
```toml
version = 1

[cases.psf_baseline]
kind = "psf"
data = "psf_baseline.txt"
shape = [64, 64]
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

For now the harness uses plain text arrays (`DelimitedFiles`) because that keeps
the validation path dependency-light. If `.npz` becomes more convenient, we can
add a test-only reader later without changing the manifest structure.

## Podman container (recommended)
Use the repo's `Containerfile` + `environment.yml` to build a reproducible
environment for generating reference outputs.

Note: `aotools` is pinned via pip in `environment.yml` to ensure availability.

Build:
```
podman build -t oopao-ref -f Containerfile .
```

Run (mount local repo, install editable, generate outputs):
```
podman run --rm -v "$PWD":/workspace:Z oopao-ref \
  python -m pip install -e .
```

## Tolerances
- Record per-dataset tolerances (relative and absolute).
- Note expected deviations for CPU vs GPU and `Float32` vs `Float64`.
