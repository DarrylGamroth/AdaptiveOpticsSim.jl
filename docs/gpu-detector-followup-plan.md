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
- public detector APIs stay unchanged:
  - `Detector(...)`
  - `capture!`
  - `capture_stack!`
- generic host-mirror and generic per-slice correction code remains as fallback-only behavior for unsupported surfaces

Highest-value next targets:
1. Detector finalization on maintained frame surfaces
   - move `apply_saturation!` and `apply_quantization!` off coarse host-mirror behavior where semantics already match on GPU
   - acceptance:
     - CPU vs AMDGPU parity
     - CPU vs CUDA parity
     - no regression on maintained HEART surfaces

2. Generalized batched detector capture
   - reduce or replace `_capture_stack_generalized!` batch-by-batch frame loops for maintained surfaces involving:
     - `psf_sampling > 1`
     - `binning > 1`
     - `readout_window`
     - `output_precision`
   - acceptance:
     - stack result matches per-frame CPU semantics
     - backend parity on maintained GPU hardware

3. Window extraction and packed output writeout
   - add stack-aware/device-aware paths for:
     - `readout_window`
     - integer output packing in `write_output!`
   - acceptance:
     - output product matches current CPU behavior
     - no hidden host round-trips on maintained GPU surfaces

4. HgCdTe readout-product surfaces
   - profile whether readout-product construction is a maintained bottleneck
   - only optimize if it shows up in maintained workloads

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
5. only then add stack-aware window/output packing if the benchmarks justify it

Stop/go criteria:
- keep a new GPU path only if:
  - it is numerically equivalent to CPU within maintained tolerances
  - it improves or at least does not materially regress warmed runtime
  - it does not broaden backend-specific special cases unnecessarily
