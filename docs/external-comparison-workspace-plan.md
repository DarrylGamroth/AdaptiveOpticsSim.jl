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

Status:

- completed 2026-04-03

Objective:

- create a clean external home for comparison work so future REVOLT/OOPAO/
  SPECULA integration no longer depends on a simulator fork

Outputs:

- new `../AdaptiveOpticsComparisons` workspace
- top-level README
- clear ownership note: comparison workspace, not simulator fork

Recommended execution order:

1. create the new directory and top-level Julia project
2. add baseline docs and ownership rules
3. define the minimal workspace structure
4. record external neighbor assumptions explicitly

Primary paths in scope:

- new external workspace root:
  `../AdaptiveOpticsComparisons`
- expected initial files:
  - `../AdaptiveOpticsComparisons/Project.toml`
  - `../AdaptiveOpticsComparisons/README.md`
  - `../AdaptiveOpticsComparisons/docs/workspace-overview.md`
  - `../AdaptiveOpticsComparisons/docs/env-setup.md`
  - `../AdaptiveOpticsComparisons/docs/result-policy.md`
  - `../AdaptiveOpticsComparisons/julia/`
  - `../AdaptiveOpticsComparisons/python/`
  - `../AdaptiveOpticsComparisons/contracts/`
  - `../AdaptiveOpticsComparisons/results/`
  - `../AdaptiveOpticsComparisons/assets/`

Tasks:

- create the new workspace
- add Julia `Project.toml`
- add environment/bootstrap docs
- record supported external neighbors:
  - `AdaptiveOpticsSim.jl`
  - `REVOLT`
  - OOPAO
  - SPECULA

Detailed work:

- create `../AdaptiveOpticsComparisons` as a new standalone workspace rather
  than as a fork of `AdaptiveOpticsSim.jl`
- add a minimal Julia `Project.toml` whose first responsibility is orchestration
  and reproducibility, not package development
- add a README that states:
  - this workspace is for comparison engineering
  - it is not the core simulator
  - it consumes `AdaptiveOpticsSim.jl` as a dependency/input
- create the initial subdirectories:
  - `docs/`
  - `julia/runners/`
  - `julia/support/`
  - `python/wrappers/`
  - `contracts/`
  - `results/archived/`
  - `results/manifests/`
  - `assets/`
- write `docs/workspace-overview.md` covering:
  - workspace role
  - supported comparison participants
  - ownership rules
  - boundary with `AdaptiveOpticsSim.jl`
- write `docs/env-setup.md` covering the expected external neighbors and how
  they are discovered:
  - `../AdaptiveOpticsSim.jl`
  - `../REVOLT`
  - OOPAO checkout or environment
  - SPECULA/REVOLT Python env
- write `docs/result-policy.md` covering:
  - what gets archived
  - which results are authoritative
  - how reruns supersede older artifacts

Deliverables:

- bootstrapped comparison workspace root
- initial Julia project file
- baseline ownership and environment docs
- empty-but-meaningful directory skeleton for later migration phases

Minimum README schema:

- purpose
- non-goals
- relationship to `AdaptiveOpticsSim.jl`
- supported external neighbors
- quick-start commands
- where migrated REVOLT assets will live

Minimum `Project.toml` expectations:

- package/environment name for the workspace
- dependency on `AdaptiveOpticsSim`
- only the minimum direct Julia dependencies needed for orchestration and
  result writing

Acceptance:

- the workspace can be checked out and understood without reading the old fork

Acceptance criteria:

- a new maintainer can identify:
  - what the workspace is for
  - what stays in `AdaptiveOpticsSim.jl`
  - where REVOLT/OOPAO/SPECULA assets are expected
- the workspace layout is sufficient to begin `ECW-2` without redesign
- no copied simulator source tree is introduced just to bootstrap the workspace

Evidence to record:

- exact workspace root path
- initial `Project.toml`
- README path
- baseline docs added
- any intentionally deferred bootstrap items

Recorded evidence:

- workspace root:
  `../AdaptiveOpticsComparisons`
- initial project:
  `../AdaptiveOpticsComparisons/Project.toml`
- bootstrap README:
  `../AdaptiveOpticsComparisons/README.md`
- baseline docs:
  - `../AdaptiveOpticsComparisons/docs/workspace-overview.md`
  - `../AdaptiveOpticsComparisons/docs/env-setup.md`
  - `../AdaptiveOpticsComparisons/docs/result-policy.md`
- intentionally deferred bootstrap items:
  - Python environment locking
  - containerization
  - CI automation
  - REVOLT asset migration

Notes:

- do not migrate assets during `ECW-1`
- do not introduce Python env management policy beyond the minimal discovery/
  setup contract needed for later phases

### `ECW-2`: Migrate REVOLT comparison assets

Status:

- completed 2026-04-03

Objective:

- move REVOLT-specific comparison scaffolding out of the fork and into the new
  comparison workspace without changing core package boundaries

Outputs:

- migrated REVOLT-oriented scenario builders and assets
- no comparison-only files left unique to `revolt-real` without a disposition

Recommended execution order:

1. inventory the REVOLT-specific comparison surface in `revolt-real`
2. move the low-risk comparison-only assets first
3. migrate the Julia-side REVOLT builders and parameter files
4. add compatibility shims or path updates in the new workspace
5. write a residual-in-fork inventory

Primary source paths in scope:

- `../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/`
- `../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/parameter_files/`
- `../AdaptiveOpticsSim.jl-revolt-real/benchmarks/assets/revolt_like/`
- comparison-side runner scripts that are only useful in REVOLT/OOPAO/SPECULA
  context

Primary destination paths in scope:

- `../AdaptiveOpticsComparisons/julia/runners/revolt/`
- `../AdaptiveOpticsComparisons/julia/support/revolt/`
- `../AdaptiveOpticsComparisons/contracts/revolt/`
- `../AdaptiveOpticsComparisons/assets/revolt_like/`
- `../AdaptiveOpticsComparisons/docs/`

Tasks:

- move `scripts/revolt/` content into the new workspace
- move REVOLT-like parameter files into comparison-side contracts/assets
- move comparison-side asset generation scripts there
- inventory everything left in `revolt-real`

Detailed work:

- create a migration inventory table for every file under:
  - `scripts/revolt/`
  - `scripts/revolt/parameter_files/`
  - `benchmarks/assets/revolt_like/`
- classify each item as:
  - move to comparison workspace
  - keep temporarily in `revolt-real`
  - candidate for later upstream review
- migrate the REVOLT Julia scenario builders into the new workspace with the
  same behavior but comparison-workspace-relative paths
- migrate the REVOLT TOML parameter files into `contracts/revolt/` or
  `assets/revolt_like/` depending on whether they are:
  - scenario definitions
  - camera/config normalization assets
- migrate any comparison-only helper scripts that generate or normalize
  REVOLT-like assets
- update path handling so the new workspace becomes the execution home for
  these scenarios
- write a `docs/revolt-migration-inventory.md` in the new workspace recording:
  - original path
  - new path
  - status
  - notes
- leave behind, in `revolt-real`, only:
  - files needed temporarily for transition, or
  - files explicitly classified for later upstream review

Deliverables:

- migrated REVOLT scenario-builder subtree in the comparison workspace
- migrated parameter/config asset subtree
- migrated REVOLT-like asset subtree
- residual inventory showing what remains in `revolt-real` and why

Minimum migration inventory schema:

- original path
- destination path
- classification
- migration status
- rationale

Residual fork inventory categories:

- transitional keep
- likely removable after `ECW-3`
- explicit upstream-review candidate
- intentionally historical

Acceptance:

- the new workspace can run the current REVOLT-like Julia scenarios without
  requiring the fork as the execution home

Acceptance criteria:

- the comparison workspace owns the maintained REVOLT Julia scenario builders
  and their parameter assets
- path handling no longer assumes `revolt-real` as the scenario home
- there is an explicit written record of everything still left in
  `revolt-real`
- no file is left behind in the fork simply because it was not classified

Evidence to record:

- migrated path list
- residual path list
- first successful runner command executed from the new workspace
- any blockers that prevent complete REVOLT-asset migration

Recorded evidence:

- migrated inventory:
  `../AdaptiveOpticsComparisons/docs/revolt-migration-inventory.md`
- maintained REVOLT support subtree:
  - `../AdaptiveOpticsComparisons/julia/support/revolt/common.jl`
  - `../AdaptiveOpticsComparisons/julia/runners/revolt/pwfs.jl`
  - `../AdaptiveOpticsComparisons/julia/runners/revolt/pwfs_unmod.jl`
  - `../AdaptiveOpticsComparisons/julia/runners/revolt/shwfs.jl`
  - `../AdaptiveOpticsComparisons/julia/runners/revolt/generate_revolt_like_benchmark_assets.jl`
- migrated parameter/config assets:
  `../AdaptiveOpticsComparisons/contracts/revolt/`
- migrated REVOLT-like HIL assets:
  `../AdaptiveOpticsComparisons/assets/revolt_like/`
- residual fork inventory:
  recorded in `../AdaptiveOpticsComparisons/docs/revolt-migration-inventory.md`
- first successful new-workspace runner command:
  `julia --project=. --startup-file=no -e 'include("julia/runners/revolt/shwfs.jl"); setup = revolt_setup(); println(setup.label)'`
- transition compatibility check:
  `julia --project=. --startup-file=no -e 'include("scripts/revolt/shwfs.jl"); setup = revolt_setup(); println(setup.label)'`
- blockers:
  none for the scoped REVOLT asset migration

Notes:

- do not refactor scenario semantics during `ECW-2`
- preserve behavior first, improve structure later
- keep upstream review separate from asset migration; `ECW-2` is not the phase
  where core code moves into `main`

### `ECW-3`: Re-anchor the maintained comparison harness

Status:

- completed 2026-04-03

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

Recorded evidence:

- comparison-workspace harness runner:
  `../AdaptiveOpticsComparisons/julia/runners/run_cross_package_benchmarks.jl`
- comparison-workspace contract:
  `../AdaptiveOpticsComparisons/contracts/cross_package.toml`
- comparison-workspace REVOLT runtime runner:
  `../AdaptiveOpticsComparisons/julia/runners/profile_revolt_runtime.jl`
- comparison-workspace harness guide:
  `../AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md`
- first archived re-anchored result:
  `../AdaptiveOpticsComparisons/results/archived/2026-04-03-ecw3-baseline.toml`
- manifest for preferred results:
  `../AdaptiveOpticsComparisons/results/manifests/cross_package.toml`
- acceptance command:
  `julia --project=. --startup-file=no julia/runners/run_cross_package_benchmarks.jl`
  from `../AdaptiveOpticsComparisons`

Notes:

- the core package still retains package-local profile scripts and historical
  archived cross-package evidence
- the comparison workspace is now the default home for new multi-package
  contracts and archived comparison runs
- `CP-06` is explicitly deferred in the comparison workspace until a distinct
  grouped-platform external runner is normalized

### `ECW-4`: Upstream interface-only gaps

Status:

- completed 2026-04-03

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

Recorded evidence:

- interface audit:
  `../AdaptiveOpticsComparisons/docs/upstream-interface-audit.md`
- localized comparison-side reference harness:
  `../AdaptiveOpticsComparisons/julia/support/reference_harness.jl`
- updated comparison harness runner:
  `../AdaptiveOpticsComparisons/julia/runners/run_cross_package_benchmarks.jl`
- updated REVOLT runtime runner:
  `../AdaptiveOpticsComparisons/julia/runners/profile_revolt_runtime.jl`
- post-audit archived result:
  `../AdaptiveOpticsComparisons/results/archived/2026-04-03-ecw4-baseline.toml`

Decision:

- no immediate upstream patch is required

Notes:

- exported detector/runtime output surfaces are sufficient for the current
  maintained comparison workflows
- the comparison workspace no longer reaches directly into
  `AdaptiveOpticsSim.jl/test/reference_harness.jl`
- no generic builder or benchmark API gap is currently blocking the migrated
  comparison flows

### `ECW-5`: Retire or freeze `revolt-real`

Status:

- completed 2026-04-03

Outputs:

- explicit disposition of `../AdaptiveOpticsSim.jl-revolt-real`

Recommended end state:

- archive/freeze it as a transitional historical branch, or
- reduce it to a minimal reference with a README pointing to the new
  comparison workspace

Acceptance:

- `revolt-real` is no longer treated as an active long-term integration fork

Recorded evidence:

- frozen-fork README redirect:
  `../AdaptiveOpticsSim.jl-revolt-real/README.md`
- comparison-workspace migration plan and inventories:
  - `../AdaptiveOpticsComparisons/docs/revolt-migration-inventory.md`
  - `../AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md`
  - `../AdaptiveOpticsComparisons/docs/upstream-interface-audit.md`

Disposition:

- `revolt-real` remains as a transitional historical fork with compatibility
  shims, not as the default home for new comparison work

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
