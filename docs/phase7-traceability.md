# Phase 7 Traceability

This matrix records the concrete artifacts that close phase 7.

| ID | Requirement | Implementation | Verification | Status |
| --- | --- | --- | --- | --- |
| P7-1 | Port most high-value tutorials as runnable Julia scripts | `examples/tutorials/*.jl`, `docs/julia-tutorial-mappings.md` | `test/runtests.jl` tutorial smoke testset | Covered |
| P7-2 | Add deterministic regression against OOPAO outputs | `test/reference_harness.jl`, `test/reference_data/` | `Pkg.test()` reference-harness testsets | Covered |
| P7-3 | Maintain a small committed OOPAO reference bundle | `test/reference_data/manifest.toml` plus array files | `load_reference_bundle`, `validate_reference_case` | Covered |
| P7-4 | Publish a user guide | `docs/user-guide.md` | Manual review plus README links | Covered |
| P7-5 | Publish an API reference | `docs/api-reference.md` | Manual review plus README links | Covered |
| P7-6 | Document backend-selection strategy around traits and `KernelAbstractions.jl` | `docs/user-guide.md`, `docs/julia-port-design.md` | Manual review | Covered |
| P7-7 | Keep tutorial ports tied to deterministic, testable workflows | `examples/tutorials/common.jl`, `test/runtests.jl` | Tutorial smoke testset | Covered |

Current OOPAO bundle scope:
- PSF, geometric Shack-Hartmann, diffractive Shack-Hartmann, Pyramid, BioEdge,
  GSC optical gains, and transfer-function rejection.
- Geometric Shack-Hartmann still uses an explicit convention adapter where OOPAO
  and AdaptiveOptics.jl intentionally expose different public slope ordering and
  units.
- LiFT, closed-loop traces, and tomography remain on the forward roadmap.
