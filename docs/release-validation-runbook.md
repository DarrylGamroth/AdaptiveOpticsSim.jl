# Release Validation Runbook

Status: active

## Purpose

Use this runbook for a release candidate or another cross-cutting gate. During
ordinary development, run the smallest registered selectors that cover the
change.

## Before Validation

1. Start from the exact candidate commit with a clean worktree.
2. Record `git rev-parse HEAD`, `git status --short`, Julia version, platform,
   thread settings, and package environment.
3. Inspect the diff for unintended manifest, artifact, or generated-file
   changes.
4. Set deterministic CPU thread counts to one unless the test explicitly studies
   parallelism.
5. Confirm that optional hardware is healthy before starting a billed or remote
   run.

## CPU Gate

~~~bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
~~~

Then run the core examples when their surface changed:

~~~bash
./scripts/run_core_examples.sh
~~~

The full suite covers the registered composition in `test/test_selection.jl`.
Passing it does not qualify an accelerator that was skipped.

## AMDGPU Gate

Run on the local ROCm host when accelerator-relevant source, graph nodes,
detectors, WFS paths, or backend infrastructure changed:

~~~bash
julia --project=test/amdgpu --startup-file=no   -e 'using Pkg; Pkg.instantiate()'
julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu.jl
~~~

Record AMDGPU/ROCm/driver/device versions and confirm scalar indexing is
disabled.

## CUDA Gate

Run on the WSL CUDA host for the same accelerator-relevant scope:

~~~bash
ssh wsl 'cd /home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl &&
  julia --project=test/cuda --startup-file=no -e "using Pkg; Pkg.instantiate()" &&
  julia --project=test/cuda --startup-file=no test/runtests_cuda.jl'
~~~

Record CUDA.jl, toolkit, driver, and device versions. CUDA remains a manual
validation target rather than a continuously available release gate.

## AArch64 Gate

Use the Raspberry Pi only for generic CPU, packaging, or architecture-sensitive
changes. Record the exact remote revision and environment. Do not substitute a
different commit or dirty tree for the candidate.

## Algorithm Graph And HIL Checks

For graph-runtime changes, retain focused evidence in addition to the full CPU
suite:

~~~bash
julia --project=. test/ci/run_selected_tests.jl algorithm-graphs
julia --project=. benchmarks/benchmark_pre_hil_backend_latency.jl
~~~

Record warmed service-time distributions and graph configuration. The result
is a self-paced service-cost measurement, not fixed-arrival PipeWire/RTC
latency. Run instrument-specific graph and pyRTC gates from each downstream
instrument package when that package is part of the release candidate.

The lockstep HIL tests must cover complete-frame publication, same-sequence
command adoption, nonfinite-command rejection, failure stop, reset, exact
host/device copies, and warmed allocation for the claimed path.

## Optional Proper.jl Check

When Proper integration changes, use the explicit environment containing
Proper.jl:

~~~bash
julia --project=. examples/integrations/proper_hil_coronagraph.jl
julia --project=. scripts/profile_proper_hil_coronagraph.jl
~~~

Record the Proper revision and prescription inputs. Passing integration tests
does not qualify SPIDERS optical fidelity.

## External Reference Checks

Run only when the relevant sibling repositories and datasets are present:

~~~bash
ADAPTIVEOPTICS_VALIDATE_TRUTH=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1 ./scripts/run_release_validation.sh
~~~

Reference regeneration is not a routine release step. Review provenance and
scientific tolerances before replacing a frozen artifact.

## Unified Entry Point

~~~bash
./scripts/run_release_validation.sh
~~~

Optional flags:

- `ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1`
- `ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1`
- `ADAPTIVEOPTICS_VALIDATE_AMDGPU=1`
- `ADAPTIVEOPTICS_VALIDATE_CUDA=1`
- `ADAPTIVEOPTICS_VALIDATE_TRUTH=1`
- `ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1`

Do not skip the CPU gate for a release candidate unless an equivalent exact-head
result is already recorded and no CPU-relevant content changed.

## PR Evidence

Record:

- exact commit
- local/remote host identity and platform
- Julia and backend versions
- commands and selectors
- pass/fail/skip counts
- benchmark configuration and summary when applicable
- known untested surfaces
- whether GitHub Actions was intentionally not run

Use at most one justified exact-head GitHub validation run after local review.
Do not rerun an unchanged failure as a debugging method.

## Acceptance

A release candidate is acceptable only when required local gates pass, skips are
consistent with the documented support boundary, no stale current claims point
to removed implementation, and the candidate is the exact revision described by
the evidence.
