# Platform Manifest Defer

Date: 2026-04-03

Status: active

Plan traceability:

- [`PSP-20`](./platform-strengthening-plan.md)
- direction IDs: `PSR-01`, `PSR-02`

## Purpose

This note records the explicit defer decision for scenario manifests and
config-file-style orchestration after the platform-strengthening plan.

## Decision

Scenario manifests remain deferred.

They are not the next main-platform priority, and they are not the primary
maintained interface.

The maintained interface remains:

- Julia scripts
- typed constructors
- typed platform/scenario objects
- maintained benchmark and validation runners

## Why The Defer Is Correct

The package’s strengths are:

- Julia-native composition
- type-driven specialization
- script-first reproducibility
- benchmark and validation workflows that already compose from code

Switching the package toward config-file-first orchestration now would:

- flatten the type-driven design
- create a second scenario model before the platform model truly needs it
- risk reintroducing ambiguity around ownership and defaults

## Allowed Future Form

If manifests are revisited later, the allowed form is:

- optional
- secondary
- layered on top of the typed platform model

That means:

- manifests should lower into `SinglePlatformConfig`,
  `GroupedPlatformConfig`, and related typed builders
- manifests should not replace the script-first workflow
- manifests should not become the authoritative runtime ownership model

## Revisit Conditions

This defer should be revisited only if at least one of these becomes true:

1. a repeated external orchestration use case cannot be kept clean with typed
   Julia builders alone
2. benchmark/validation scenario exchange would materially improve from a
   lightweight declarative layer
3. an optional adapter package needs a stable serialization format on top of
   the typed platform model

If revisited, it should start as:

- a separate plan, and likely
- an optional layer or adapter, not a core workflow replacement

## Current Rule

Until a later plan changes it explicitly:

- no config-file-first orchestration
- no manifest-driven primary user workflow
- no parallel scenario model that bypasses typed Julia builders
