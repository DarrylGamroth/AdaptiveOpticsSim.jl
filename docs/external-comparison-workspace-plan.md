# External Comparison Workspace Plan

Date: 2026-04-03

Status: active

## Purpose

This plan defines how to replace the long-lived
[`../AdaptiveOpticsSim.jl-revolt-real`](../AdaptiveOpticsSim.jl-revolt-real)
fork with a dedicated external comparison workspace.

The goal is to:

- keep `AdaptiveOpticsSim.jl` focused on core simulation/runtime capability
- move cross-package comparison concerns out of the core package
- preserve the useful REVOLT/OOPAO/SPECULA comparison assets
- stop maintaining a drifting fork as if it were part of the core product

## Decision Summary

The default decision is:

- do **not** merge `../AdaptiveOpticsSim.jl-revolt-real` into `main`
- create a separate comparison workspace in `../`
- upstream only code that is clearly generic core capability

At the time of this plan, the default upstream set is intentionally empty.

That means:

- assume nothing from `revolt-real` belongs in `main` by default
- require an explicit justification before upstreaming any piece

## Proposed Workspace

Recommended name:

- `../AdaptiveOpticsComparisons`

Recommended role:

- engineering/comparison workspace, not core package replacement

The workspace should own:

- OOPAO comparison runners
- SPECULA / REVOLT comparison runners
- Julia wrapper scripts that call into `AdaptiveOpticsSim.jl`
- scenario normalization contracts
- environment setup for Python-side assets
- archived comparison results
- comparison-specific docs and reproducibility instructions

The workspace should not become:

- a second core AO simulator
- a long-lived fork of `AdaptiveOpticsSim.jl`
- the place where core runtime improvements are implemented first

## Recommended Layout

Suggested top-level layout:

```text
AdaptiveOpticsComparisons/
  Project.toml
  Manifest.toml
  README.md
  docs/
    workspace-overview.md
    scenario-contracts.md
    env-setup.md
    result-policy.md
  julia/
    runners/
    support/
    adapters/
  python/
    env/
    wrappers/
    revolt/
    oopao/
    specula/
  contracts/
    revolt.toml
    oopao.toml
    specula.toml
  results/
    archived/
    manifests/
  assets/
    revolt_like/
    normalization/
```

Design rules:

- keep Julia-side comparison logic separate from Python-side wrappers
- keep environment/bootstrap logic separate from scenario contracts
- keep archived results versioned and discoverable
- keep all comparison-specific assets out of `AdaptiveOpticsSim.jl` unless they
  are required by core package tests or docs

## Interface Between The Core Package And The Comparison Workspace

The comparison workspace should consume `AdaptiveOpticsSim.jl` through:

- exported APIs
- stable scripts
- stable runtime/readout surfaces
- typed platform/config builders

It should not depend on:

- internal scratch buffers
- non-exported internal helper functions
- branch-local implementation details that are not part of the maintained
  public or scripted workflow surface

If the comparison workspace discovers a missing stable hook, the default action
should be:

1. identify the missing core interface
2. upstream that interface to `AdaptiveOpticsSim.jl`
3. keep the comparison-specific logic external

## Upstream Decision Rule

Use this rule for every artifact currently in `revolt-real`.

### Upstream to `main` only if all are true

- it improves core simulation/runtime capability
- it is generic and not REVOLT-specific
- it is useful without OOPAO/SPECULA comparison context
- it fits the core package boundary
- it can be validated within the normal package validation model

### Keep external if any are true

- it exists mainly to compare against another package
- it normalizes external scenario conventions
- it depends on Python env setup or neighboring trees
- it is REVOLT/OOPAO/SPECULA specific scaffolding
- it is benchmark/archive/result management rather than simulator capability

## Current Default Classification

### Keep external

These should move to the new comparison workspace, not to `main`:

- `scripts/revolt/`
- `scripts/revolt/parameter_files/`
- REVOLT-like benchmark asset generation and normalization logic
- cross-package environment wrappers
- Python/REVOLT/SPECULA invocation scripts
- archived comparison result handling
- comparison-only docs and normalization notes
- fork-local benchmark harness code whose purpose is package-to-package
  comparison rather than core performance evidence

### Candidate upstream only if generalized

These may contain pieces worth upstreaming later, but only after explicit
generalization:

- benchmark runners that reveal missing stable core hooks
- reusable typed scenario builders that are not REVOLT-specific
- generic platform/runtime entry points discovered during comparison work
- core bug fixes or performance fixes first found while running comparison
  scenarios

### Do not upstream as-is

These should not be copied into `main` unchanged:

- REVOLT-specific TOML scenario parameter files
- fork-specific scenario builders tightly coupled to REVOLT naming/assets
- comparison-specific benchmark result schemas
- external environment bootstrap logic
- comparison-side copies of package docs/plans

## Migration Phases

### `ECW-1`: Workspace bootstrap

Outputs:

- new `../AdaptiveOpticsComparisons` workspace
- top-level README
- clear ownership note: comparison workspace, not simulator fork

Tasks:

- create the new workspace
- add Julia `Project.toml`
- add environment/bootstrap docs
- record supported external neighbors:
  - `AdaptiveOpticsSim.jl`
  - `REVOLT`
  - OOPAO
  - SPECULA

Acceptance:

- the workspace can be checked out and understood without reading the old fork

### `ECW-2`: Migrate REVOLT comparison assets

Outputs:

- migrated REVOLT-oriented scenario builders and assets
- no comparison-only files left unique to `revolt-real` without a disposition

Tasks:

- move `scripts/revolt/` content into the new workspace
- move REVOLT-like parameter files into comparison-side contracts/assets
- move comparison-side asset generation scripts there
- inventory everything left in `revolt-real`

Acceptance:

- the new workspace can run the current REVOLT-like Julia scenarios without
  requiring the fork as the execution home

### `ECW-3`: Re-anchor the maintained comparison harness

Outputs:

- comparison-side harness entry points
- explicit contract files
- archived result policy

Tasks:

- decide which harnesses remain in `AdaptiveOpticsSim.jl` and which move out
- keep only narrow core-facing contracts in `AdaptiveOpticsSim.jl`
- move heavier cross-package orchestration into the comparison workspace
- make the comparison workspace the default home for:
  - external env setup
  - Python runner invocation
  - multi-package result archival

Acceptance:

- cross-package comparisons no longer depend on `revolt-real` as a forked code
  host

### `ECW-4`: Upstream interface-only gaps

Outputs:

- a small list of actual upstream patches, if needed

Tasks:

- audit which comparison workflows still require unstable or awkward core hooks
- upstream only:
  - missing stable APIs
  - generic builders
  - generic readout/benchmark surfaces
- do not upstream comparison logic itself

Acceptance:

- the comparison workspace uses stable core interfaces without copying
  comparison logic into `main`

### `ECW-5`: Retire or freeze `revolt-real`

Outputs:

- explicit disposition of `../AdaptiveOpticsSim.jl-revolt-real`

Recommended end state:

- archive/freeze it as a transitional historical branch, or
- reduce it to a minimal reference with a README pointing to the new
  comparison workspace

Acceptance:

- `revolt-real` is no longer treated as an active long-term integration fork

## Immediate Upstream Decision

Current recommendation:

- upstream nothing from `revolt-real` right now

Rationale:

- the core package already has the main platform/runtime/orchestration surfaces
  we want to preserve
- `revolt-real` is now primarily comparison scaffolding and scenario
  normalization
- moving that scaffolding into `main` would worsen package boundaries

If a later audit finds a genuinely generic builder or API gap, handle it as a
small targeted upstream patch, not as a fork merge.

## First Concrete Tasks

Recommended near-term order:

1. create `../AdaptiveOpticsComparisons`
2. migrate `scripts/revolt/` and related REVOLT assets there
3. move cross-package environment/bootstrap logic there
4. keep `AdaptiveOpticsSim.jl` with only the narrow maintained comparison
   surfaces that still belong in core
5. freeze `revolt-real`

## Acceptance For The Overall Migration

The migration is complete when:

- cross-package benchmarking no longer depends on `revolt-real` as an active
  fork
- the new workspace can run the maintained REVOLT/OOPAO/SPECULA comparison
  surfaces
- `AdaptiveOpticsSim.jl` remains focused on core simulation/runtime concerns
- any upstreamed items are small, generic, and explicitly justified
