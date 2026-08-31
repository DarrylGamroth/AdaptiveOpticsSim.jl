# User Guide

Status: active

AdaptiveOpticsSim.jl provides adaptive-optics physics, prepared numerical
operations, and a small complete-frame graph runtime. It supports direct Julia
workflows on CPU and covered CUDA/AMDGPU paths.

## Choose A Composition Surface

Use direct Julia when the model has generated topology, multiple rates,
conditional execution, row-time optical sampling, or research-specific control
flow.

Use `AlgorithmGraphs` when the model is a static, single-rate sequence that
publishes complete frames. Graphs can be declared in Julia or in a versioned
TOML file.

Use `PreparedGraphHILBoundary` when an external RTC consumes each completed
simulated frame and returns one complete command for the next frame. PipeWireAO
or another application owns transport and wall-clock pacing.

## Direct Optical Model

~~~julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.WavefrontSensors

telescope = Telescope(
    resolution=32,
    diameter=8.0,
    central_obstruction=0.1,
)
source = Source(band=:I, magnitude=8.0)
pupil = PupilFunction(telescope)

atmosphere = MultiLayerAtmosphere(
    telescope;
    r0=0.15,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction=[0.0, 90.0],
    altitude=[0.0, 5000.0],
)

renderer = prepare_atmosphere_renderer(
    atmosphere,
    telescope,
    source,
)
rng = runtime_rng(1)
epoch = advance_by!(atmosphere, 1e-3; rng)
render_atmosphere!(pupil, renderer, atmosphere, epoch)

wfs = ShackHartmannWFS(
    telescope;
    n_lenslets=4,
    mode=Diffractive(),
    pixel_scale_arcsec=0.1,
    n_pix_subap=6,
)
measurement = measure!(wfs, pupil, source)
~~~

The atmosphere owns evolving turbulence state. The renderer owns a prepared
direction mapping and writes the caller's pupil product. Rendering does not
advance time or consume RNG.

## Explicit Closed Loop

A direct loop keeps the science visible and is often the best starting point:

~~~julia
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control

dm = DeformableMirror(telescope; n_act=4, influence_width=0.3)
interaction = interaction_matrix(
    dm,
    wfs,
    PupilFunction(telescope),
    source;
    amplitude=0.1,
)
reconstructor = ModalReconstructor(interaction; gain=0.5)
command = similar(dm.state.coefs)

for _ in 1:100
    epoch = advance_by!(atmosphere, 1e-3; rng)
    render_atmosphere!(pupil, renderer, atmosphere, epoch)

    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())

    measure!(wfs, pupil, source)
    reconstruct!(command, reconstructor, slopes(wfs))
    @. command = -command
    set_command!(dm, command)
end
~~~

Controller delay, detector acquisition, and command timing should be explicit
when they matter to the experiment. Use caller-owned state and separate RNGs
for independently replayable stochastic components.

## Detector Acquisition

~~~julia
using AdaptiveOpticsSim.Detectors

detector = Detector(
    noise=NoisePhotonReadout(2.0),
    full_well=30_000.0,
    bits=12,
    output_type=UInt16,
)

detector_rng = runtime_rng(2)
measure!(wfs, pupil, source, detector; rng=detector_rng)
frame = wfs_detector_image(wfs, detector)
~~~

Detector quantization (`bits` and `full_well`) is separate from the Julia output
element type. Use `output_type=nothing` for floating-point internal readout.

For staged WFS models, keep optical formation, detector acquisition, and
measurement estimation separate:

1. `form_wfs_optical_products!`
2. `acquire_wfs_observation!`
3. `estimate_wfs_measurement!`

## Complete-Frame Algorithm Graphs

A graph definition contains a concrete tuple of nodes, graph inputs/outputs,
same-frame links, explicit delayed links, and startup ndarray parameters.
Preparation validates every port and binds exact storage.

~~~julia
using AdaptiveOpticsSim.AlgorithmGraphs

definition = algorithm_graph(
    (
        discrete_integrator_node(
            :controller;
            extent=2,
            sample_period_s=1e-3,
            input_schema="example.residual.f32/1",
            output_schema="example.command.f32/1",
            gain=0.2f0,
            tau_s=0.02f0,
        ),
    );
    name=:controller,
    inputs=(
        graph_input(
            :residual,
            :controller => :input,
            zeros(Float32, 2),
        ),
    ),
    outputs=(
        graph_output(:command, :controller => :output),
    ),
)

graph = prepare_algorithm_graph(definition)
step_graph!(graph)
command = graph_output(graph, :command)
~~~

Use a delayed link for feedback. Direct links must follow node order and cannot
form a cycle.

## TOML Graphs

TOML is the convenient static authoring surface:

~~~julia
bindings = (
    pdm_command=zeros(Float32, 277),
    pdm_actuator_coordinates=actuator_coordinates,
)

definition = load_algorithm_graph(
    "examples/graphs/revolt_classic_hil.toml";
    bindings,
)
graph = prepare_algorithm_graph(definition)
~~~

Graph files embed scalar configuration but not large calibration arrays or
arbitrary filenames. `bindings` supplies caller-owned inputs, delayed initial
values, and large ndarray parameters by stable name.

Maintained examples include:

- `examples/graphs/revolt_classic_hil.toml`
- `examples/graphs/revolt_classic_hil_fast_dm.toml`
- `examples/graphs/revolt_copper_hil.toml`
- `examples/graphs/revolt_copper_hil_fast_dm.toml`

The fast profile changes the provisional Gaussian DM evaluation to a separable
regular-grid implementation without reducing detector or pupil resolution.

## External RTC Lockstep

~~~julia
boundary = prepare_graph_hil_boundary(
    graph;
    command_input=:pdm_command,
    frame_output=:shwfs_frame,
)

while running
    sequence = step_hil_frame!(boundary)
    send_frame(sequence, hil_frame_buffer(boundary))

    receive_command!(
        hil_command_buffer(boundary),
        sequence,
    )
    adopt_hil_command!(boundary, sequence)
end
~~~

The first frame uses the command present at preparation. Each subsequent frame
is blocked until a finite same-sequence command is adopted. The boundary owns
no socket, PipeWire object, retry policy, or wall-clock schedule.

For GPU graphs, graph arrays remain device-resident and the HIL boundary makes
the completed frame/device-command copies explicit through ordinary host
`Array` buffers.

## Use AOS With pyRTC

The maintained pyRTC integration is a Linux, CPU-shared-memory validation path.
AdaptiveOpticsSim runs the optical plant in Julia. A separate Python process
runs pyRTC's actual `SlopesProcess` and `Loop`, consumes complete detector
frames, and returns complete deformable-mirror commands. No camera or
deformable-mirror hardware is required for the reference systems.

The validated examples are controlled, noiseless, 25-command systems:

- a 64-by-64 Shack–Hartmann detector frame with 104 pyRTC signals
- a 36-by-36 Pyramid detector frame with 344 pyRTC signals

These reference systems establish transport, calibration, and closed-loop
behavior without depending on unavailable instrument calibration. They are not
REVOLT instrument models.

### Obtain And Install The Software

Install Git, Python, and [Julia 1.12 or newer](https://julialang.org/downloads/).
pyRTC currently supports Python 3.9 through 3.13; use a dedicated virtual
environment. Clone both repositories into one working directory:

~~~sh
mkdir aos-pyrtc-work
cd aos-pyrtc-work
export AO_HIL_ROOT="$PWD"

git clone https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl.git
git clone https://github.com/jacotay7/pyRTC.git

python3 -m venv "$AO_HIL_ROOT/.venv"
"$AO_HIL_ROOT/.venv/bin/python" -m pip install --upgrade pip
"$AO_HIL_ROOT/.venv/bin/python" -m pip install -e "$AO_HIL_ROOT/pyRTC"
~~~

Enter the AOS checkout and select the same Python interpreter for the native
worker and the optional PythonCall oracle:

~~~sh
cd "$AO_HIL_ROOT/AdaptiveOpticsSim.jl"
export PYRTC_ROOT="$AO_HIL_ROOT/pyRTC"
export PYRTC_PYTHON="$AO_HIL_ROOT/.venv/bin/python"
export JULIA_PYTHONCALL_EXE="$PYRTC_PYTHON"

julia --project=. --startup-file=no \
  -e 'using Pkg; Pkg.instantiate(); using AdaptiveOpticsSim'
julia --project=examples/integrations/pyrtc --startup-file=no \
  -e 'using Pkg; Pkg.instantiate()'
~~~

Re-export `PYRTC_ROOT`, `PYRTC_PYTHON`, and `JULIA_PYTHONCALL_EXE` after
opening a new shell.

### Verify Shared-Memory Compatibility

Run the short interoperability suite first:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/runtests.jl
~~~

This checks Julia-to-pyRTC and pyRTC-to-Julia vectors and matrices, NumPy C
order, metadata admission, segment ownership, cleanup, and zero warmed Julia
heap allocation for immediate publication and reads. A passing command exits
with status zero and prints passing test summaries.

The integration uses pyRTC's conventional `wfs`, `signal`, `signal2D`, and
`wfc` shared-memory names. Stop another pyRTC session before running it. The
launcher refuses to reuse an occupied name.

### Run Both Closed-Loop Reference Systems

Run the Shack–Hartmann and Pyramid systems against a separate pyRTC process:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/run_process_reference_matrix.jl
~~~

Each run performs push-pull calibration, verifies that the interaction matrix
has rank 25, installs the matrix in pyRTC, injects a known disturbance, and
requires the returned command to close the AOS optical loop. The command prints
the interaction condition, initial and final residuals, convergence ratio, and
command-recovery error for both sensors.

Run the same systems as assertions in the optional integration-test matrix:

~~~sh
AOS_PYRTC_PROCESS_TESTS=1 \
  julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/runtests.jl
~~~

For debugging, `run_reference_matrix.jl` provides an in-process PythonCall
oracle. The separate-process runner above is the representative native Julia
shared-memory path.

### Adapt A Scientific Graph

A graph used with pyRTC needs one complete detector-frame output and one
complete command input bound through `PreparedGraphHILBoundary`. Keep the
scientific algorithms and graph independent of pyRTC. The application adapter
then performs one ordered exchange:

1. Step the graph and obtain its frame sequence.
2. Publish the complete host frame to `wfs`.
3. Let pyRTC publish `signal` and the next command to `wfc`.
4. Read the complete command into `hil_command_buffer(boundary)`.
5. Adopt it with the same AOS frame sequence before stepping again.

Use
[`pyrtc_process_hil.jl`](../examples/integrations/pyrtc/pyrtc_process_hil.jl)
as the executable reference for stream creation, process lifecycle,
interaction-matrix calibration, and command adoption. Use
[`pyrtc_shared_memory.jl`](../examples/integrations/pyrtc/pyrtc_shared_memory.jl)
directly only when implementing another application adapter.

The current shared-memory layout has one overwriteable slot. It is valid for
the demonstrated one-producer, one-outstanding-frame lockstep exchange. It is
not a free-running asynchronous transport, a deadline guarantee, or a dropped-
frame policy. AOS GPU graphs can copy their completed frames through the host
HIL boundary, but the GPU-to-CPU pyRTC composition has not yet been qualified.
See [`runtime-dataflow.md`](runtime-dataflow.md) for the protocol and ownership
details.

## Model Time

`ModelTimestamp` and `ModelDuration` use exact integer nanoseconds.
`PeriodicSchedule` and the model-time drivers provide deterministic frame
coordinates:

~~~julia
driver = FixedStepModelTimeDriver(
    PeriodicSchedule(ModelDuration(1_000_000)),
)
timestamp = step_graph_at!(graph, driver)
~~~

Model time is not host execution time. An application may pace these steps
against a clock without changing the graph's scientific timestamps.

## Rolling Shutters And Sub-Frame Optics

A rolling-shutter detector still publishes one complete atomic frame. To model
moving turbulence, direct Julia code or one atomic custom node evaluates the
optical path at each row or row-group integration time, accumulates private
detector state, and publishes only after every row is complete.

Do not expose partially written frames. Mid-frame DM changes and multiple
cadences are possible in direct Julia, but they are outside the version 1 static
graph scheduler until an instrument requirement justifies a bounded extension.

## GPU Execution

Prepare the graph for one exact target:

~~~julia
using AdaptiveOpticsSim.Backends

target = compute_device(device_array)
graph = prepare_algorithm_graph(definition; target)
~~~

Every bound graph array must already use native storage on that target. There is
no hidden transfer or CPU fallback. Only nodes whose algorithms support the
selected backend can be prepared.

## Proper.jl

Proper.jl remains optional. Use it directly between explicit optical products or
wrap one prepared prescription in an application-specific graph node. See
[`proper-integration-guide.md`](proper-integration-guide.md).

## Next Reading

- [`model-cookbook.md`](model-cookbook.md) for task-oriented recipes
- [`runtime-dataflow.md`](runtime-dataflow.md) for graph and ownership details
- [`api-reference.md`](api-reference.md) for module navigation
- [`supported-production-surfaces.md`](supported-production-surfaces.md) for
  qualified scope
- `examples/tutorials/` for runnable examples
