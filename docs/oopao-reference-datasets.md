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
