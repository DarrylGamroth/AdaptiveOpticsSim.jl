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
    reference_wavelength_m=500e-9,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction_deg=[0.0, 90.0],
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
advance time or consume RNG. `r0` is defined at `reference_wavelength_m`;
rendering converts the generated phase to wavelength-independent OPD in metres.

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

For a GPU-resident application that has independent CPU work before it needs
the frame, split submission from completion:

~~~julia
ticket = step_graph_async!(graph)
do_independent_cpu_work()
wait_graph_step!(ticket)
command = graph_output(graph, :command)
~~~

One graph permits only one pending frame. Its inputs, outputs, node state,
workspaces, delayed values, and RNG state remain owned by that frame until
`wait_graph_step!` returns. Use separate prepared graphs for genuinely
independent concurrent work. The lockstep HIL API already waits at its required
host-visible frame boundary.

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

The short validated references are controlled, noiseless, 25-command systems:

- a 64-by-64 Shack–Hartmann detector frame with 104 pyRTC signals
- a 36-by-36 Pyramid detector frame with 344 pyRTC signals

These reference systems establish transport, calibration, and closed-loop
behavior without depending on unavailable instrument calibration. They are not
REVOLT instrument models.

The separate REVOLT Classic model-validation path uses the maintained
352-by-352 C-BLUE One IMX425 global-shutter CMOS graph, 188 valid
Shack–Hartmann subapertures, 376
pair-interleaved slope components, and a complete 277-element HSDM command. It
self-calibrates against the AOS model and therefore validates the complete
model/pyRTC loop without claiming compatibility with an operational REVOLT
interaction or command matrix. The normalized HSDM placement and Gaussian
influence function remain provisional.

The detector identity and operating geometry follow maintained HEART sources:
the C-BLUE One uses the global-shutter Sony IMX425, and REVOLT selects a
352-by-352 ROI at 500 Hz with 1896 µs exposure and 24 dB analog gain. The
underlying sensor has 1608-by-1104 9 µm pixels; the graph's
`pixel_scale_arcsec` is the projected optical sampling, not that physical pixel
pitch.

The maintained HEART snapshot does not state the operational pixel format or
unit calibration. The graph therefore selects the published 12-bit
high-sensitivity characterization: 2.33 e⁻ RMS read noise, 0.96 e⁻/pixel/s
dark current, 94 ke⁻ zero-gain image full well, and an approximate 0.27 QE at
the modeled 800 nm channel. It converts 24 dB to a linear post-readout voltage
gain of 15.8489. The internal unsigned ADC result is copied as integer-valued
`Float32` samples for the graph and pyRTC boundary. Replace these typical
values with the camera's measured pixel format, QE, bias/flat, noise,
saturation, and ROI calibration before claiming instrument parity.

### Obtain And Install The Software

Install Git, Python, and [Julia 1.12 or newer](https://julialang.org/downloads/).
pyRTC currently supports Python 3.9 through 3.13; use a dedicated virtual
environment. Clone AOS into one working directory:

~~~sh
mkdir aos-pyrtc-work
cd aos-pyrtc-work
export AO_HIL_ROOT="$PWD"

git clone https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl.git

python3 -m venv "$AO_HIL_ROOT/.venv"
"$AO_HIL_ROOT/.venv/bin/python" -m pip install --upgrade pip
"$AO_HIL_ROOT/.venv/bin/python" -m pip install -r \
  "$AO_HIL_ROOT/AdaptiveOpticsSim.jl/examples/integrations/pyrtc/requirements.txt"
~~~

The pyRTC
[`requirements.txt`](../examples/integrations/pyrtc/requirements.txt) file
installs the tested revision directly from the official
[pyRTC GitHub repository](https://github.com/jacotay7/pyRTC). A separate local
pyRTC checkout is not required.

Enter the AOS checkout and select the same Python interpreter for the native
worker and the optional PythonCall oracle:

~~~sh
cd "$AO_HIL_ROOT/AdaptiveOpticsSim.jl"
export PYRTC_PYTHON="$AO_HIL_ROOT/.venv/bin/python"
export JULIA_PYTHONCALL_EXE="$PYRTC_PYTHON"

julia --project=. --startup-file=no \
  -e 'using Pkg; Pkg.instantiate(); using AdaptiveOpticsSim'
julia --project=examples/integrations/pyrtc --startup-file=no \
  -e 'using Pkg; Pkg.instantiate()'
~~~

Re-export `PYRTC_PYTHON` and `JULIA_PYTHONCALL_EXE` after opening a new shell.

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

### Run The Closed-Loop Systems

Run the Shack–Hartmann and Pyramid systems against a separate pyRTC process:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/run_process_reference_matrix.jl
~~~

Each run performs push-pull calibration, verifies that the interaction matrix
has rank 25, installs the matrix in pyRTC, injects a known static disturbance,
and requires the returned command to close the AOS optical loop. It then uses
the same measured matrix with a deterministic, evolving four-layer atmosphere
and requires the corrected 750 nm on-axis Strehl ratio to exceed both the
open-loop result and an absolute acceptance floor. The command prints the
interaction condition, static-disturbance convergence and command error, and
the mean open-loop and corrected atmospheric Strehl ratios for both sensors.

Run the larger REVOLT Classic model-validation system separately:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/run_revolt_classic_hil.jl
~~~

This performs push-pull calibration for all 277 physical command elements,
uses a relative singular-value cutoff of `0.02`, and closes a 500-frame
(one modeled second) loop against the deterministic maintained five-layer
atmosphere. The cutoff rejects poorly conditioned directions that amplify the
modeled 12-bit detector and centroid noise while retaining at least 80% of the
277-dimensional command space. The test requires a finite retained interaction
subspace, corrected mean 750 nm on-axis Strehl above 0.35, at least tenfold
improvement over the open-loop mean, and more than a factor-of-two reduction
in mean pupil OPD RMS. These are simulation acceptance checks, not REVOLT
instrument qualification.

Run the same systems as assertions in the optional integration-test matrix:

~~~sh
AOS_PYRTC_PROCESS_TESTS=1 \
  julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/runtests.jl
~~~

The longer REVOLT Classic assertion is independently opt-in:

~~~sh
AOS_PYRTC_REVOLT_CLASSIC_TESTS=1 \
  julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/runtests.jl
~~~

For debugging, `run_reference_matrix.jl` provides an in-process PythonCall
oracle. The separate-process runner above is the representative native Julia
shared-memory path.

### View The Live Optical Loop

The integration environment installs pyRTC's optional Qt viewer. From a local
graphical session or an X11/Wayland-forwarded shell, run a paced Pyramid
demonstration and open a live eight-panel mosaic with:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/run_process_viewer_demo.jl pyramid 60 10
~~~

Use `shack_hartmann` in place of `pyramid` for the Shack–Hartmann reference
system. Use `revolt_classic` to run the maintained REVOLT Classic model:

~~~sh
julia --project=examples/integrations/pyrtc --startup-file=no \
  examples/integrations/pyrtc/run_process_viewer_demo.jl \
  revolt_classic 60 10
~~~

The final two arguments are the run duration in seconds and requested frame
rate. REVOLT Classic first calibrates 277 push-pull command columns, so its
viewer takes longer to start than either 25-command reference. The demo opens
the official `pyrtc-view` application and displays:

- the complete AOS detector frame in `wfs`
- pyRTC's `signal2D` slope product
- the current pyRTC command in `wfc2D` (5-by-5 for a reference system or the
  exact 19-by-19 HSDM physical grid for REVOLT Classic)
- AOS uncompensated, deformable-mirror, and residual pupil OPD maps
- diffraction-limited-normalized open-loop and corrected 750 nm science PSFs

The console reports pupil-supported, piston-removed open-loop and residual OPD
RMS, the open-loop and corrected on-axis Strehl ratios, and command RMS. The
PSF normalization uses the exact unsampled diffraction-limited peak; the
Strehl calculation uses the coherent on-axis pupil sum and therefore includes
image motion from residual tip and tilt.

The viewer polls at half the requested wall-clock demonstration rate, so each
normal refresh interval spans two graph periods. This avoids pyRTC's transient
`PAUSED` label between frames without altering graph ordering or RTC command
adoption. The compact reference graphs advance atmosphere model time by 1 ms
per frame; REVOLT Classic advances it by 2 ms to match its configured 500 Hz
C-BLUE One cadence. The lower wall-clock rate intentionally slows the modeled
loop for inspection.

The reference graphs evolve the maintained four-layer reference atmosphere;
REVOLT Classic evolves its maintained five-layer atmosphere. None substitutes
a deformable-mirror-shaped synthetic disturbance. Close the viewer window to
stop the demonstration before its requested duration.

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
