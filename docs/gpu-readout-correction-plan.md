## GPU Readout Correction Plan

Goal: move the maintained batched frame readout-correction models onto stack-aware GPU paths without changing user-facing detector semantics.

Scope:
- `ReferencePixelCommonModeCorrection`
- `ReferenceRowCommonModeCorrection`
- `ReferenceColumnCommonModeCorrection`
- `ReferenceOutputCommonModeCorrection`
- `CompositeFrameReadoutCorrection`

Rules:
- keep `Detector(...)`, `capture!`, and `capture_stack!` APIs unchanged
- prefer stack-aware kernels over per-frame or per-slice correction loops
- use existing batched detector scratch buffers; do not allocate in hot paths
- require CPU semantic parity and backend parity before keeping a GPU path

Execution order:
1. implement GPU batched kernels for row, column, and output common-mode correction
2. make composite correction dispatch through stage-local batched GPU paths
3. add CPU regression tests comparing batched corrected stacks against per-frame capture on distinct inputs
4. add CUDA and AMDGPU parity checks for all maintained correction models
5. benchmark synthetic corrected stacks on both GPU backends to confirm the stack-aware path is a real win

Acceptance:
- corrected batched stacks match per-frame CPU capture for each maintained model
- CUDA and AMDGPU corrected stacks match CPU within maintained tolerances
- no regressions in existing optional backend checks
- synthetic GPU correction benchmarks show that the GPU path avoids the old generic 3D correction overhead
