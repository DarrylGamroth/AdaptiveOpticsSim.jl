## GPU Detector Follow-On Plan

Goal: continue replacing coarse detector host-mirror and generic per-frame GPU execution with maintained device-native paths, without changing user-facing detector APIs.

Terminology:
- "batched stack" means multiple detector frames are processed together in one stack-shaped array, typically `batch x rows x cols`.
- On the maintained Shack-Hartmann path, that batch dimension is the set of lenslet spot frames for a single simulation step, not multiple time steps of the same output frame.
- In other detector surfaces, the same stack shape can also represent multiple frames/products processed together.

Current state:
- maintained readout-correction models now have device-native GPU paths for:
  - direct single-frame capture
  - fixed batched stack capture
- maintained direct detector finalization now uses device-native GPU saturation/quantization on maintained frame surfaces
- public detector APIs stay unchanged:
  - `Detector(...)`
  - `capture!`
  - `capture_stack!`
- generic host-mirror and generic per-slice correction code remains as fallback-only behavior for unsupported surfaces
- generalized batched detector fast-path prototypes for windowed/integer output were evaluated and not kept:
  - they preserved CPU parity
  - they did not improve warmed runtime on CUDA
  - they regressed warmed runtime on AMDGPU relative to the kept per-frame loop
- direct `write_output!` profiling on the maintained windowed/integer HgCdTe surface did not justify more code:
  - AMDGPU: `capture ≈ 0.488 ms`, `write_output! ≈ 0.0245 ms`
  - CUDA: `capture ≈ 4.91 ms`, `write_output! ≈ 0.0119 ms`
  - conclusion: window extraction / integer packing is already device-native and not the dominant cost on the maintained direct-frame path
- a windowed HgCdTe readout-product fast path was prototyped and rejected:
  - it preserved deterministic CPU/GPU parity on AMDGPU and CUDA
  - but it regressed warmed runtime against the old full-cube branch on both backends
  - measured microbenchmark:
    - AMDGPU: `new ≈ 62.31 ms`, `old ≈ 26.29 ms`
    - CUDA: `new ≈ 34.56 ms`, `old ≈ 1.08 ms`
  - conclusion: repeated windowed slice copies are a worse GPU shape than the existing full-cube path on this maintained surface

Highest-value next targets:
1. Detector finalization on maintained frame surfaces
   - move `apply_saturation!` and `apply_quantization!` off coarse host-mirror behavior where semantics already match on GPU
   - acceptance:
     - CPU vs AMDGPU parity
     - CPU vs CUDA parity
     - no regression on maintained HEART surfaces

2. HgCdTe readout-product surfaces
   - only revisit if a maintained workload shows the current full-cube path is a real bottleneck
   - any replacement should avoid per-read windowed slice copies and prove a warmed runtime win on both AMDGPU and CUDA

Design rules:
- preserve detector semantics and user-facing APIs
- add internal execution seams, not public backend knobs
- keep GPU scratch/device buffers reused and preallocated
- benchmark after warmup; do not trust first-call timings
- treat old host-mirror/generic logic as compatibility fallback, not the target maintained GPU path

Execution order:
1. profile maintained detector finalization on AMDGPU and CUDA after warmup
2. implement device-native saturation/quantization where the current plan still host-mirrors
3. add focused parity tests for direct frame finalization on GPU
4. profile `_capture_stack_generalized!` and identify whether generalized batching is worth a maintained GPU path
   - current result: not worth keeping on maintained AMDGPU/CUDA surfaces yet
5. profile direct `write_output!` on maintained GPU detector surfaces
   - current result: already cheap; no further action justified
6. only revisit HgCdTe readout products if they show up as a maintained workload bottleneck

Stop/go criteria:
- keep a new GPU path only if:
  - it is numerically equivalent to CPU within maintained tolerances
  - it improves or at least does not materially regress warmed runtime
  - it does not broaden backend-specific special cases unnecessarily
