# AdaptiveOpticsSim.jl-revolt-real

This repository is now a transitional historical fork.

Active simulator development belongs in:

- [../AdaptiveOpticsSim.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl)

Active cross-package comparison work belongs in:

- [../AdaptiveOpticsComparisons](/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons)

## Current Role

This fork is retained only so older comparison workflows and historical
benchmarks continue to have a concrete reference tree during the migration away
from a fork-based comparison model.

It should not be treated as the default execution home for new comparison work.

## What Moved

Maintained REVOLT comparison assets and harness ownership now live in
`AdaptiveOpticsComparisons`, including:

- REVOLT scenario builders
- REVOLT parameter contracts
- REVOLT-like HIL assets
- the cross-package benchmark harness
- new archived cross-package results

## Expected End State

The intended end state is that this repository remains only as:

- a transitional compatibility reference, or
- a frozen historical branch

If you are starting new work, use `AdaptiveOpticsComparisons` instead of this
fork.
