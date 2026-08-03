# LLM Agent Instructions for AdaptiveOpticsSim.jl
#
# Project name: AdaptiveOpticsSim.jl (Julia port of OOPAO).
#
# Core design principles:
# - Use multiple dispatch + traits; avoid OO-style inheritance patterns.
# - Prefer dispatch/traits or small helper functions over `isa` checks in package code.
# - Separate params (immutable structs) from state (mutable structs).
# - Use `Plan` only for a run-immutable numerical or physical execution
#   contract. Keep persistent scientific state, replaceable scratch, caller-
#   visible products, and exact input/output/context bindings in separate
#   owners.
# - Use explicit `!`-mutating functions for hot paths.
# - Preallocate workspaces; avoid allocations in inner loops.
# - Do not use `Memory` directly in package or extension implementation.
#   Use FixedSizeArrays.jl's `FixedSizeArray` with a concrete element type for
#   homogeneous armed host storage; a `Vector` may be a cold builder but must be
#   sealed before execution. Use tuples, concrete unions, family-grouped
#   owners, or purpose-built owners for bounded heterogeneity.
# - Do not store `Any` in a struct field or collection element type. Do not
#   return `Any` from preparation into execution ownership. Unconstrained
#   `::Any`, `<:Any`, or `Vararg{Any}` method signatures are dispatch wildcards,
#   not stored representations, and should be judged separately.
# - Do not store abstract element types or uninstantiated parametric families in
#   hot-path registries. A container change does not repair type erasure; use a
#   concrete owner or an explicit function barrier.
# - Keep core free of hard-coded file formats (no baked-in FITS).
# - Prefer structured errors (custom exception types) over print-and-return.
# - Use Logging.jl for diagnostics; avoid logging inside hot loops.
# - Prefer Plots.jl for visualization; keep plotting optional.
# - Favor idiomatic Julia patterns in API design and implementation.
# - Follow SciML style conventions: https://docs.sciml.ai/SciMLStyle/stable/
# - Use lower-case directory names; Julia type names may use CamelCase, but
#   source tree directories should stay lower-case for portability and clarity.
# - This package is undergoing a breaking refactor. Remove superseded APIs and
#   update callers directly; do not add deprecated aliases, synthetic property
#   forwarding, state views, or compatibility adapters unless the user
#   explicitly requests a compatibility surface.
# - Keep the root export surface for routine cross-domain workflows and
#   canonical domain modules. Dense domain APIs belong to their real owner
#   module: export routine vocabulary there, mark stable advanced seams
#   `public`, and leave implementation details unmarked.
#
# Agent reasoning budget:
# - Reserve `max` reasoning for architecture, planning, cross-cutting design
#   decisions, and difficult gate or final reviews. Do not spend `max`
#   reasoning on routine implementation, test execution, CI monitoring, or
#   mechanical documentation changes.
# - When parallel agents are authorized, the coordinating agent may remain at
#   `max` while it delegates concrete, bounded implementation work at `high` by
#   default. Use `xhigh` only for unusually subtle numerical, concurrency,
#   inference, or zero-allocation work.
# - Use a smaller context fork for bounded workers when they do not need the
#   full conversation. Keep planning and acceptance decisions with the
#   coordinating agent.
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
# - Centralize RNG state in an explicit single-writer state or prepared owner;
#   fixed seeds for reproducibility. Do not classify an evolving RNG as
#   replaceable workspace scratch.
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
# Terminology:
# - Use the canonical scientific and API terms in docs/glossary.md.
# - Plane names describe optical location; product names describe quantities.
# - Define new public scientific terms in the glossary before adding APIs.
#
# Tutorials:
# - Port most OOPAO tutorials to Julia scripts/notebooks.
# - Use deterministic inputs for reference comparisons where possible.

# Reference docs:
# - docs/documentation-map.md
# - docs/glossary.md
# - docs/extension-guide.md
# - docs/julia-tutorial-mappings.md
# - docs/roadmap.md
# - docs/model-validity-matrix.md
# - docs/supported-production-surfaces.md
# - docs/deterministic-simulation.md
