Status: active

Plan traceability:

- [`PVP-07`](./post-review-platform-plan.md)
- [`CP-02`](./cross-package-benchmark-inventory.md)

# REVOLT-Like SH Benchmark Contract

## Purpose

This note defines the maintained comparability contract for the `CP-02`
REVOLT-like Shack-Hartmann runtime benchmark.

The goal is not to claim full scenario identity between:

- `AdaptiveOpticsSim.jl` `main`
- `../AdaptiveOpticsSim.jl-revolt-real`
- `../REVOLT/Python/specula`

The goal is to make the current Julia-to-Julia comparison stronger and more
honest by distinguishing:

- fields that are normalized and expected to match
- fields that are intentionally different today

## Normalized Fields

The maintained `CP-02` contract currently expects the following Julia-side
fields to match:

- backend
- response mode
- number of Shack-Hartmann subapertures
- WFS detector sensor family

The `main` side now runs the REVOLT-like HIL benchmark with the WFS sensor set
to `cmos` so that the detector family matches the REVOLT-aligned branch.

## Known Differences

The following differences are currently accepted and recorded rather than
treated as parity failures:

- pupil resolution
- exact frame-response model
- detector output shape and payload layout
- deformable-mirror command layout
- science-detector payload

These differences mean `CP-02` is still a normalized medium-class runtime
comparison, not a full numerical-equivalence or deployment-parity benchmark.

## Interpretation Rule

When reading `CP-02` results:

- matching normalized fields means the benchmark is comparing like-for-like at
  the coarse scenario-contract level
- timing and allocation differences should still be interpreted in the presence
  of the known structural differences above
- a future stronger `CP-02` baseline should reduce the known-difference list,
  not silently ignore it

## Reproduction

Run the maintained harness:

```bash
julia --project=. --startup-file=no scripts/run_cross_package_benchmarks.jl
```

The archived `CP-02` result is written under:

- [benchmarks/results/cross_package](../benchmarks/results/cross_package)

The executable contract source is:

- [cross_package.toml](../benchmarks/contracts/cross_package.toml)
