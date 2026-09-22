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
  `WavefrontSensors`, `Calibration`, `Tomography`, `Ensembles`, and
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
Curvature, and LiFT paths. The general staged contract is:

1. `prepare_wfs_optics` / `form_wfs_optical_products!`
2. `prepare_wfs_acquisition` / `acquire_wfs_observation!`
3. `prepare_wfs_estimation` / `estimate_wfs_measurement!` when the WFS owner
   provides an estimation stage

`WFSObservation` and `WFSMeasurement` are typed products with explicit
metadata and units.

The AOS Shack–Hartmann, Pyramid, Bi-O-edge, and Zernike surfaces implement
stages 1 and 2 only.
Their complete detector observations are package boundaries; maintained
operational estimation is provided by FilterGraphAlgorithms. Registered
FilterGraphAlgorithms v0.5.0 provides the Zernike complete-frame normalized
pupil-signal estimator and the Curvature paired-image and paired-channel
estimators; JuliaFilterGraph v0.2.3 provides their graph runtime. AOS Curvature
ends at complete detector acquisition. LiFT retains
its AOS physical forward and inverse APIs until its approved
AdaptiveOpticsCalibration inverse target is released and adopted; LiFT is a
phase-retrieval workflow, not a slope or centroid estimator. AOS has no
weighted-centroid estimator or graph surface. The explicitly named
`geometric_wavefront_slopes!` operation remains available for plant truth and
reference calculations, not as a detector-derived RTC measurement.

A prepared four-pupil LGS optics plan owns the sodium-profile kernel derived
during preparation. To change sodium-layer altitudes or relative weights,
construct the changed `LGSSource` and prepare a new optics plan; profile changes
are not a repeated-path cache refresh.

## Calibration And RTC Integration

~~~julia
using AdaptiveOpticsSim.Calibration
~~~

`Calibration` owns simulated interaction-response acquisition, physical
calibration observables, modal bases, fitting, plant-side runtime
materialization, optical-gain calibration, and identification workflows.
Reusable inverse methods and compact-SVD products come from
`AdaptiveOpticsCalibration.Reconstructors`. Common entry points include
`interaction_matrix`, `control_matrix`, `ao_calibration`, `modal_basis`, and
`compute_optical_gains!`.

`KarhunenLoeveBasis` is owned by
`AdaptiveOpticsCalibration.ModalBases`. AOS supplies the sampled physical DM,
pupil measure, and atmospheric OPD covariance, then materializes the accepted
calibration product on the plant's runtime backend.

The maintained cross-package reference is
[`examples/integrations/filter_graph_algorithms/`](../examples/integrations/filter_graph_algorithms/),
where AOS supplies the plant and FilterGraphAlgorithms/JuliaFilterGraph supplies
the RTC chain. The composing fixture is CPU-resident and retains concrete
floating-point element types and caller-owned arrays for allocation-free
execution. AOS and FGA qualify their accelerator paths independently; the
fixture does not claim shared device storage or one captured graph across
packages.

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
- `modal_opd_expansion_node`

The built-in graph catalog stops at complete detector frames. For a maintained
Shack–Hartmann estimator and RTC composition, use the
[`filter_graph_algorithms` fixture](../examples/integrations/filter_graph_algorithms/),
where FilterGraphAlgorithms owns the RTC operations.

`builtin_graph_node_types()` reports the file-loader type map.

### Preparation And Execution

~~~julia
definition = load_algorithm_graph(path; bindings=bindings)
graph = prepare_algorithm_graph(definition)
step_graph!(graph)
reset_graph!(graph)
~~~

`StreamGraphExecution()` is the default. A graph with independent branches may
declare bounded coarse execution groups. On CUDA or AMDGPU, nodes in each group
use distinct retained streams; the host preserves the same schedule serially:

~~~julia
graph = prepare_algorithm_graph(
    definition;
    target,
    execution=GroupedStreamGraphExecution(
        (:guide_star_1, :guide_star_2),
        (:wfs_1, :wfs_2),
        (:tomography,),
    ),
)
~~~

Flattened group names must exactly match declaration order, and direct links
must cross group boundaries. Grouped execution uses fixed reusable device
events; it does not infer groups or create CPU tasks.

On CUDA or AMDGPU, an application may instead request native command-graph
replay when every node in the graph is explicitly qualified:

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
stream execution into the captured graph. A downstream package must qualify its
complete instrument graph separately.

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
