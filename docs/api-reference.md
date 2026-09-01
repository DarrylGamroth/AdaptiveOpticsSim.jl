# API Reference

Status: active

This is a navigation guide, not a duplicated symbol inventory. Routine
cross-domain vocabulary is exported by `AdaptiveOpticsSim`. Dense APIs are
exported from their canonical modules. Stable advanced seams are marked
`public` and should be imported explicitly from the owner module.

For exact signatures and keyword defaults, use Julia help mode:

~~~julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics

?prepare_direct_imaging
?AdaptiveOpticsSim.AlgorithmGraphs.prepare_graph_node
~~~

## Root Package

~~~julia
using AdaptiveOpticsSim
~~~

The root exports:

- canonical modules: `Backends`, `Optics`, `Atmospheres`, `Detectors`,
  `WavefrontSensors`, `Calibration`, `Control`, `Tomography`, `Ensembles`, and
  `AlgorithmGraphs`
- shared structured errors
- `ScientificProfile`, `FastProfile`, and `default_fidelity_profile`
- `SplitMix64RNG`, `runtime_rng`, and `deterministic_reference_rng`
- the routine graph-definition, preparation, execution, model-time, and
  lockstep HIL vocabulary described below

The root does not forward each dense domain API.

## Backends

~~~julia
using AdaptiveOpticsSim.Backends
~~~

Routine entry points include `CPUBackend`, `CUDABackend`,
`AMDGPUBackend`, `MetalBackend`, `backend`, and `compute_device`. Advanced
public seams describe exact devices, availability, backend identity, and
device-array allocation.

Backend selection identifies array and execution mechanics. It does not select
AO physics, graph topology, or an implicit host/device transfer.

## Optics

~~~julia
using AdaptiveOpticsSim.Optics
~~~

Major families:

- `Telescope`, `Source`, `LGSSource`, `Asterism`, and spectral/extended sources
- explicit products: `PupilFunction`, `ElectricField`, `IntensityMap`
- plane, coordinate, spectral, normalization, and spatial-measure metadata
- aperture masks, NCPA, OPD maps, spatial filters, and Zernike bases
- Fraunhofer/Fresnel propagation and direct imaging
- `DeformableMirror`, modal controllable optics, tip/tilt and focus stages
- reusable WFS optics such as microlens arrays and focal-plane masks

Prepared paths use names such as `prepare_direct_imaging`,
`form_direct_image!`, `prepare_pupil_field`, `fill_electric_field!`,
`propagate_field!`, `set_command!`, `update_surface!`, and `apply_surface!`.

## Atmospheres

~~~julia
using AdaptiveOpticsSim.Atmospheres
~~~

Use `KolmogorovAtmosphere`, `MultiLayerAtmosphere`, or
`InfiniteMultiLayerAtmosphere` for model construction. The explicit time flow
is:

~~~julia
renderer = prepare_atmosphere_renderer(atmosphere, telescope, source)
epoch = advance_by!(atmosphere, 1e-3; rng=rng)
render_atmosphere!(pupil, renderer, atmosphere, epoch)
~~~

Direction batches are prepared with `prepare_atmosphere_direction_batch` and
executed with `render_atmosphere_directions!`.

## Detectors

~~~julia
using AdaptiveOpticsSim.Detectors
~~~

The module owns conventional frame sensors and counting/channel detectors.
Major sensor families include `CCDSensor`, `EMCCDSensor`, `CMOSSensor`,
`HgCdTeSensor`, `HgCdTeAvalancheArraySensor`, `InGaAsSensor`, SPAD, MKID, and
linear APD models.

Use `Detector` for the configured acquisition owner, `capture!` for direct
capture, and `prepare_detector_acquisition` for a prepared acquisition path.
Shutter, sampling, response, defects, persistence, thermal, and readout models
are explicit values.

## Wavefront Sensors

~~~julia
using AdaptiveOpticsSim.WavefrontSensors
~~~

Implemented WFS families include Shack-Hartmann, Pyramid, Bi-O-edge, Zernike,
Curvature, and LiFT paths. High-level `measure!` remains convenient. The
maintained staged contract is:

1. `prepare_wfs_optics` / `form_wfs_optical_products!`
2. `prepare_wfs_acquisition` / `acquire_wfs_observation!`
3. `prepare_wfs_estimation` / `estimate_wfs_measurement!`

`WFSObservation` and `WFSMeasurement` are typed products with explicit
metadata and units.

## Calibration And Control

~~~julia
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control
~~~

`Calibration` owns interaction/control matrices, modal bases, inverse policies,
fitting, optical-gain calibration, and identification workflows. Common entry
points include `interaction_matrix`, `control_matrix`, `ao_calibration`,
`kl_modal_basis`, and `compute_optical_gains!`.

`Control` owns reconstructors, delay lines, discrete integration, closed-loop
correction, and controller/DM mappings. Common entry points include
`reconstruct!`, `shift_delay!`, `apply_closed_loop_correction!`, and the
projection/range operations.

## Tomography

~~~julia
using AdaptiveOpticsSim.Tomography
~~~

Tomography owns guide-star and layer geometry, model-based or
interaction-matrix reconstruction, fitting, and DM command projection. Use
`build_reconstructor`, `assemble_reconstructor_and_fitting`,
`reconstruct_wavefront_map`, and `dm_commands`.

## Ensembles

~~~julia
using AdaptiveOpticsSim.Ensembles
~~~

`SimulationEnsemble` applies an explicit coarse-grained execution policy to
independent model members. Policies include sequential, deterministic, threaded,
backend-stream, AcceleratedKernels, and Dagger execution. This facility is for
offline/coarse parallelism, not RTC deadline scheduling.

## Algorithm Graphs

~~~julia
using AdaptiveOpticsSim.AlgorithmGraphs
~~~

### Definitions

- `AlgorithmNodeDefinition`
- `AlgorithmGraphDefinition`
- `GraphInputDefinition` and `GraphOutputDefinition`
- `AlgorithmLink` and `DelayedAlgorithmLink`
- `algorithm_node`, `algorithm_graph`, `graph_input`, `graph_output`,
  `link`, and `delayed_link`

Use `load_algorithm_graph(path; bindings=...)` to compile a versioned TOML graph
or construct the same definition directly in Julia.

### Native Nodes

Factory functions include:

- `multilayer_atmosphere_opd_node`
- `deformable_mirror_surface_node`
- `gaussian_deformable_mirror_surface_node`
- `grid_gaussian_deformable_mirror_surface_node`
- `pupil_opd_composition_node`
- `shack_hartmann_rate_node` and `pyramid_rate_node`
- `ccd_detector_acquisition_node`, `cmos_detector_acquisition_node`, and
  `emccd_detector_acquisition_node`
- `shack_hartmann_centroid_node` and
  `shack_hartmann_slope_selection_node`
- `control_matrix_reconstruction_node`
- `discrete_integrator_node` and `closed_loop_correction_node`
- `modal_opd_expansion_node`

`builtin_graph_node_types()` reports the file-loader type map.

### Preparation And Execution

~~~julia
definition = load_algorithm_graph(path; bindings=bindings)
graph = prepare_algorithm_graph(definition)
step_graph!(graph)
reset_graph!(graph)
~~~

`StreamGraphExecution()` is the default. On CUDA or AMDGPU, an application may
request native command-graph replay when every node in the graph is explicitly
qualified:

~~~julia
graph = prepare_algorithm_graph(
    definition;
    target,
    execution=CapturedGraphExecution(),
)

captured_graph_node_count(graph)
~~~

Capture is strict: preparation records the complete node sequence and
delayed-link commits as one native graph, and fails when any node is unqualified
or when the selected backend cannot capture the operation. The current built-in
qualified operations are the finite multilayer-atmosphere node, coordinate and
regular-grid Gaussian DM evaluation, pupil-OPD composition, diffractive
Shack-Hartmann rate formation, and Pyramid rate formation with either
maintained modulation strategy. Unit-binning full-frame CCD, global-shutter
CMOS, and linear sequential EMCCD acquisition nodes also qualify for their
built-in simple response models. The coordinate-Gaussian node currently
qualifies only with the built-in linear-static actuator response. Preparation
rejects richer or unsupported detector and actuator models instead of mixing
stream execution into the captured graph. Consequently, both maintained
REVOLT Gaussian profiles can be captured in full.

`step_graph!` completes the frame before returning. Advanced GPU applications
may use the bounded capacity-one split boundary:

~~~julia
ticket = step_graph_async!(graph)
wait_graph_step!(ticket)
~~~

`graph_step_pending(graph)` reports ticket ownership. The committed
`graph_step_sequence(graph)` changes only after the wait succeeds. Do not read
outputs or mutate inputs while a ticket is pending.

Use `prepared_graph_node` for qualified inspection. Node adapters implement the
public `graph_node_ports`, `prepare_graph_node`, `step_graph_node!`, and
`reset_graph_node!` protocol. An adapter may additionally specialize
`enqueue_graph_node!` to defer a proven same-context completion boundary; the
default calls `step_graph_node!`. `graph_node_capture_capability` remains
conservative unless an adapter proves fixed addresses, device-resident evolving
state, and replay without host mutation, allocation, result queries, or
synchronization.

### Model Time

`ModelTimestamp` and `ModelDuration` are exact integer-nanosecond values.
`PeriodicSchedule` maps frame sequence to time. Drivers include
`FixedStepModelTimeDriver`, `PreparedBoundaryModelTimeDriver`, and captured
timestamp drivers. Use `step_graph_at!` when a graph step must receive an
explicit model timestamp.

### Lockstep HIL

~~~julia
boundary = prepare_graph_hil_boundary(
    graph;
    command_input=:dm_command,
    frame_output=:detector_frame,
)

sequence = step_hil_frame!(boundary)
copyto!(hil_command_buffer(boundary), rtc_command)
adopt_hil_command!(boundary, sequence)
frame = hil_frame_buffer(boundary)
~~~

The boundary exchanges complete host arrays and enforces one command response
per completed frame. It does not implement transport or wall-clock pacing.

## Optional Integrations

Optional GPU packages extend `Backends` and covered scientific operations.
Proper.jl is application-owned; see
[`proper-integration-guide.md`](proper-integration-guide.md). Plotting belongs
in the companion `AdaptiveOpticsSimPlots.jl` package.
