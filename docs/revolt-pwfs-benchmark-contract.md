Status: active

Plan traceability:

- [`PVP-08`](./post-review-platform-plan.md)
- [`CP-03`](./cross-package-benchmark-inventory.md)

# REVOLT-Like PWFS Benchmark Contract

## Purpose

This note defines the maintained comparability contract for the `CP-03`
REVOLT-like representative Pyramid WFS benchmark.

The current maintained surface compares:

- `AdaptiveOpticsSim.jl` `main`
- `../AdaptiveOpticsSim.jl-revolt-real`

on the unmodulated representative PWFS contract.

## Normalized Fields

The maintained `CP-03` contract currently expects the following fields to match
between the two Julia implementations:

- backend
- response mode
- model label
- WFS family
- number of subapertures
- WFS detector sensor family
- science detector sensor family
- pupil resolution
- WFS frame-response family
- detector output shape

This makes `CP-03` stronger than the earlier contract-only state: the main
package now has a first-class REVOLT-like PWFS runner with the same coarse
scenario identity as the REVOLT-aligned branch.

## Known Differences

The following differences are still accepted and documented:

- gain-detector handling is metadata-only on `main`, not a maintained runtime
  path in this benchmark
- reconstructor build methods are similar in intent but not normalized as a
  cross-package parity surface
- timing breakdown internals are implementation-specific even when the total
  runtime surface is comparable

So `CP-03` is now a maintained representative runtime comparison, but it is
not yet a claim of full cross-platform numerical equivalence.

## External SPECULA Status

`CP-03` is first-class in `main`, but the external SPECULA/REVOLT PWFS runner
is still not a maintained execution participant in this contract.

Reason:

- the available REVOLT SPECULA tree does not currently provide a maintained
  runnable PWFS benchmark entry equivalent to the Julia contract surface
- the visible PWFS config asset exists, but the surrounding executable path is
  not yet normalized into the harness

That external alignment remains future benchmark work rather than part of the
current acceptance rule for `PVP-08`.

## Reproduction

Run the maintained harness:

```bash
julia --project=. --startup-file=no scripts/run_cross_package_benchmarks.jl
```

The executable contract source is:

- [cross_package.toml](../benchmarks/contracts/cross_package.toml)

Archived results are stored under:

- [benchmarks/results/cross_package](../benchmarks/results/cross_package)
