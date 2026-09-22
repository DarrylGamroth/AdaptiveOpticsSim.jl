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
  signatures; the plant signature fingerprints the explicit optical, WFS
  sampling, observation, and detector configuration, while the subaperture
  order is taken from the prepared FGA measurement plan;
- command `n` is adopted only after frame `n` completes and becomes active for
  frame `n + 1`;
- invalid frames block plant advancement without changing the adopted command;
  there is no implicit hold-last-command policy; and
- the warmed composed CPU step is inferred and allocates zero Julia heap
  bytes.

The integration tests also retain the pre-split controller-command routing,
constraint-feedback, fixed-frame-delay, and deterministic S1 boundary
trajectory as numerical migration oracles. They run the public FGA/JFG and
AOC interfaces only; AOS remains the plant owner.

S4 is a separate complete-frame Pyramid fixture. AOS forms a diffractive
four-pupil photon-rate frame and performs one explicit noiseless detector
acquisition. The composing layer transposes the acquired `(x, y)` detector
storage into a preallocated `(row=y, column=x)` FGA image, then uses FGA 0.3.0's
`PyramidImageF32`, `PyramidReconstructorF32`, and leaky integrator to produce
one adopted command. Its flat acquired frame defines the explicit reference
I4Q signal; the support, unity optical gain, pupil order `q1, q2, q3, q4`,
component order `(X, Y)`, units, type, model timestamps, and FGA release/tree
claim are held in one immutable identity. An independent direct I4Q equation
checks the FGA result.

S4 is intentionally a package boundary, not a shared estimator. It does not
call an AOS Pyramid estimator, measurement API, slope product, calibration
state, or compatibility API. It has no shared GPU buffers, streams, CUDA/HIP
graph, or zero-copy claim. Wrong sequence, model timestamp, identity,
calibration or order signature, schema, shape, numeric type, non-finite data,
discontinuity, corruption, or an FGA failure blocks the exchange before the
adopted command changes. `reset_s4_pyramid!` discards that blocked exchange and
restores FGA control state; it does not silently hold the previous command.

The composing bridge is intentionally CPU-resident and uses preallocated host
storage. AOS plant execution and FGA RTC execution are qualified independently
on CUDA and AMDGPU. This fixture does not claim shared device storage, a shared
accelerator stream, zero-copy exchange, or one CUDA Graph or HIP Graph spanning
both packages.

S4 also contains a Bi-O-edge detector-frame fixture. AOS owns only the
diffractive four-pupil photon-rate formation and noiseless detector acquisition;
FGA 0.3.0 owns the calibrated Bi-O-edge differential signal. The boundary uses
the same explicit `(x, y)` AOS to `(row=y, column=x)` FGA transpose. Its frozen
asymmetric oracle binds q1/q2/q3/q4 as top-left, bottom-left, bottom-right,
top-right, respectively; origins `((0,0), (2,0), (2,2), (0,2))`; and the
non-rectangular support `[true true; false true]`. It checks both FGA
normalization policies with nonzero reference signal and unequal optical gain.

The nested environment resolves AdaptiveOpticsCalibration,
FilterGraphAlgorithms, and JuliaFilterGraph from the configured registry and
uses this AOS checkout as its path source. Run:

```sh
JULIA_NUM_THREADS=1 julia --startup-file=no \
  --project=examples/integrations/filter_graph_algorithms \
  -e 'using Pkg; Pkg.instantiate()'

JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --startup-file=no \
    --project=examples/integrations/filter_graph_algorithms \
    examples/integrations/filter_graph_algorithms/runtests.jl
```

The CPU profile exercises both warmed S1 and S4 composed steps and verifies
zero Julia heap allocation before sampling:

```sh
AOS_FGA_PROFILE_STEPS=500000 JULIA_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --startup-file=no \
    --project=examples/integrations/filter_graph_algorithms \
    examples/integrations/filter_graph_algorithms/profile_cpu.jl
```

Develop a sibling package into this environment only when testing an
unpublished change to that package; do not commit the generated manifest.
