# AOS plant with an FGA RTC

This integration fixture is the S1 package-boundary reference. AOS forms a
complete physical Shack–Hartmann detector frame, AdaptiveOpticsCalibration
prepares the reconstructor, and the existing FilterGraphAlgorithms chain
performs centroiding, reconstruction, control, controller-to-VDM projection,
VDM-to-PDM projection, and PDM command conditioning.

The packages remain independent. This nested environment is only the
composing and verification layer; it does not provide another algorithm API.

The contract is deliberately narrow:

- all plant, calibration, and RTC arrays are `Float32`;
- AOS detector storage uses `(x, y)` axes, so the composing layer transposes
  into a preallocated FGA `(row=y, column=x)` frame;
- frame sequence and timestamp come from the AOS-side lockstep owner;
- demanded-command output metadata must retain the input sequence and time;
- the cold calibration result is bound to an immutable identity containing
  detector and estimator axes, ordered subapertures and slope pairs, physical
  actuator order, units, numeric type, and deterministic plant and estimator
  signatures;
- command `n` is adopted only after frame `n` completes and becomes active for
  frame `n + 1`;
- invalid frames block plant advancement without changing the adopted command;
  there is no implicit hold-last-command policy; and
- the warmed composed CPU step is inferred and allocates zero Julia heap
  bytes.

The composing bridge is intentionally CPU-resident and uses preallocated host
storage. AOS plant execution and FGA RTC execution are qualified independently
on CUDA and AMDGPU. This fixture does not claim shared device storage, a shared
accelerator stream, zero-copy exchange, or one CUDA Graph or HIP Graph spanning
both packages.

For local migration work, develop the three package worktrees into an
environment, then run:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --startup-file=no --project=/path/to/integration-env \
  examples/integrations/filter_graph_algorithms/runtests.jl
```

The nested `Project.toml` records package identities without embedding paths
to sibling worktrees. A checked-out FGA/JFG development tree must therefore be
developed locally until those packages are registered or a reviewed Git
source revision is selected.
