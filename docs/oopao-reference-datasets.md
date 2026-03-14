# OOPAO Reference Datasets

This document describes the deterministic OOPAO bundles used to cross-validate
AdaptiveOptics.jl.

## Principles
- Use fixed seeds and deterministic settings.
- Store small, representative datasets.
- Compare with tolerances, not bitwise equality.

## Current committed bundle

The repository now commits a small bundle under `test/reference_data/` with two
OOPAO-generated geometric Shack-Hartmann cases:

1. `shack_hartmann_geometric_ramp_xy`
2. `shack_hartmann_geometric_ramp_y`

These cases are small, deterministic, and currently stable enough to keep in
CI. They deliberately target the part of the port where we already understand
the OOPAO convention mapping.

Not yet committed:
- PSF parity, which still differs in normalization and shape details.
- Diffractive Shack-Hartmann, Pyramid, and BioEdge parity, which still differ in
  units and/or signal construction.

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
swap_halves = true
scale = 7.5949367088607559e6
```

That adapter is intentional. OOPAO’s exported `signal_2D` ordering and
calibrated slope units do not match the raw OPD-gradient convention currently
returned by AdaptiveOptics.jl.

For now the harness uses plain text arrays (`DelimitedFiles`) because that keeps
the validation path dependency-light. If `.npz` becomes more convenient, we can
add a test-only reader later without changing the manifest structure.

## Generation workflow

The generator lives in `scripts/generate_oopao_reference_bundle.py`.

Recommended approach:
1. Build a reproducible Python 3.8 image with the OOPAO dependency stack.
2. Mount both the OOPAO checkout and the AdaptiveOptics.jl checkout.
3. Copy OOPAO to a writable temporary directory inside the container before
   import. OOPAO writes `precision_oopao.npy` into its package root at import
   time.
4. Run the generator and point it at a bundle directory.

The committed geometric SH bundle was generated this way and then copied into
`test/reference_data/`.

## Tolerances
- Record per-dataset tolerances (relative and absolute).
- Note expected deviations for CPU vs GPU and `Float32` vs `Float64`.

For the current committed SH cases, tolerances are effectively exact once the
documented convention adapter is applied.
