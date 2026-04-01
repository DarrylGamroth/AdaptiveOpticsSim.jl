# Future Platform Direction

Status: active

Plan traceability:

- [`PLAN-38`](./package-review-action-plan.md)
- [`PLAN-39`](./package-review-action-plan.md)
- [`PLAN-40`](./package-review-action-plan.md)
- review IDs: `PR-35`, `PR-36`, `PR-37`, `PR-38`

## Purpose

This document captures the package-level direction after the review cleanup
phases.

Its job is to prevent future work from drifting back toward:

- parity chasing for its own sake
- unplanned control-framework expansion before the plant/runtime is ready
- science-path integration pressure leaking into the core package boundary

The phased execution plan for the next post-review work lives in
[post-review-platform-plan.md](./post-review-platform-plan.md).

## Direction Summary

The package should continue to evolve as:

- an idiomatic Julia adaptive-optics simulation platform
- a scientifically defensible and benchmarked runtime core
- a system that uses OOPAO and SPECULA as references, not as architectural
  ceilings

## OOPAO: Reference, Not Ceiling

OOPAO remains important for:

- behavior checks
- tutorial parity
- frozen regression bundles
- cross-package sanity comparisons

OOPAO should not be treated as:

- the final architectural authority
- a reason to preserve Python-era decomposition in Julia
- a reason to block improvements that exceed OOPAO in runtime structure,
  backend strategy, or validation discipline

The correct rule is:

- preserve or explain behavior where it matters
- exceed OOPAO where Julia-native design is better justified

## SPECULA: Stronger Reference For Platform Breadth

After the current cleanup work, the next broader reference should be SPECULA
more than OOPAO.

SPECULA is the stronger external reference for:

- atmosphere-aware field propagation
- broader orchestration/process layering
- richer sensing/control platform breadth
- future multi-object or platform-scale expansion ideas

This does not mean “port SPECULA.” It means:

- use SPECULA where it is the better comparison baseline
- prefer Julia-native architecture when implementing the equivalent capability

## What To Prioritize Next

If the package expands again after the cleanup phases, the preferred order is:

1. strengthen validation and benchmark evidence further where needed
2. implement targeted platform-breadth areas where SPECULA is the stronger
   reference
3. revisit broader controller/process families only after the plant/runtime
   remains stable under those new demands

Examples of better next breadth areas than more OOPAO parity:

- additional atmospheric field propagation scenarios
- richer grouped sensing/runtime orchestration
- carefully selected controller/process families with clear benchmark and
  validation goals

The currently selected next breadth milestone is:

- grouped sensing/runtime orchestration
  - plan: [grouped-runtime-plan.md](./grouped-runtime-plan.md)
  - decision memo: [next-milestone-decision-2026-04.md](./next-milestone-decision-2026-04.md)

## What To Avoid

Avoid these patterns unless a later explicit plan justifies them:

- implementing additional OOPAO features only because OOPAO has them
- pulling optional science-path tooling into core
- adding a broad control framework before the runtime/product/validation
  surfaces demand it
- using parity language as a substitute for benchmark or validation evidence

## Controller And Process Breadth

Controller/process expansion remains a future platform-level topic, not the
default next step.

That work should come only after:

- runtime ownership is clear
- benchmark surfaces exist for the target use case
- validation expectations are explicit

The package should prefer a small number of well-benchmarked, well-validated
controller/process additions over a large surface of partially integrated
families.

## Science-Path Integration Boundary

Science-path integrations should remain optional unless a later boundary review
changes that decision explicitly.

Examples:

- `Proper.jl` handoff layers
- coronagraph or downstream science-camera pipelines
- heavy ecosystem adapters used mainly for science post-processing

The preferred form is:

- extension module, or
- separate adapter package

Core should expose a narrow, stable handoff boundary rather than embedding the
science tool directly.

## How To Use This Guide

Use this document when deciding whether a new idea belongs in:

- the core roadmap
- a subsystem plan
- an optional extension
- a separate adapter package

If a proposed feature conflicts with this guide, the default answer should be
to write an explicit plan first rather than allowing the work to enter the core
incrementally.
