# LLM Agent Instructions for AdaptiveOpticsSim.jl
#
# Project name: AdaptiveOpticsSim.jl (Julia port of OOPAO).
#
# Core design principles:
# - Use multiple dispatch + traits; avoid OO-style inheritance patterns.
# - Prefer dispatch/traits or small helper functions over `isa` checks in package code.
# - Separate params (immutable structs) from state (mutable structs).
# - Use explicit `!`-mutating functions for hot paths.
# - Preallocate workspaces; avoid allocations in inner loops.
# - Keep core free of hard-coded file formats (no baked-in FITS).
# - Prefer structured errors (custom exception types) over print-and-return.
# - Use Logging.jl for diagnostics; avoid logging inside hot loops.
# - Prefer Plots.jl for visualization; keep plotting optional.
# - Favor idiomatic Julia patterns in API design and implementation.
# - Follow SciML style conventions: https://docs.sciml.ai/SciMLStyle/stable/
# - Use lower-case directory names; Julia type names may use CamelCase, but
#   source tree directories should stay lower-case for portability and clarity.
#
# Parallelism:
# - Focus on coarse-grained parallelism (sources, time steps, sweeps).
# - Avoid nested parallelism and thread oversubscription.
# - Deterministic mode should pin thread counts to 1.
#
# GPU readiness:
# - Parametrize structs by numeric type and array backend.
# - Avoid scalar indexing on GPU; minimize host-device transfers.
# - Use AbstractFFTs for FFT portability and KernelAbstractions for kernels.
#
# Determinism and validation:
# - Centralize RNG in workspace; fixed seeds for reproducibility.
# - Compare outputs against OOPAO reference datasets within tolerance.
#
# Dependencies:
# - Keep core dependencies minimal.
# - Optional features (I/O helpers, ModelingToolkit control, plotting extras)
#   should be in extension modules or separate packages.
#
# Documentation:
# - Keep docs/ intentionally small and navigable.
# - Prefer updating an existing maintained guide over adding a new document.
# - Do not add one-off plan, audit, triage, or inventory documents under docs/.
# - Use PR descriptions, issues, or git history for temporary planning records.
#
# Tutorials:
# - Port most OOPAO tutorials to Julia scripts/notebooks.
# - Use deterministic inputs for reference comparisons where possible.

# Reference docs:
# - docs/documentation-map.md
# - docs/extension-guide.md
# - docs/julia-port-design.md
# - docs/julia-tutorial-mappings.md
# - docs/roadmap.md
# - docs/model-validity-matrix.md
# - docs/supported-production-surfaces.md
# - docs/deterministic-simulation.md
