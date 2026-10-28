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
  `WavefrontSensors`, `Calibration`, `Tomography`, and
  `AlgorithmGraphs`
- shared structured errors
- `ScientificProfile`, `FastProfile`, and `default_fidelity_profile`
- `SplitMix64RNG`, `runtime_rng`, and `deterministic_reference_rng`
- the routine graph-definition, preparation, execution, model-time, and
  lockstep HIL vocabulary described below

The root does not forward each dense domain API.
`runtime_timing` is an internal benchmark helper, not a supported package
entry point.

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
Prepared owners and immutable plan accessors are supported advanced seams;
replaceable propagation/direct-imaging scratch and batch compatibility
signatures are implementation details.

### Spectral synthesis

`prepare_spectral_synthesis(envelopes, basis, coefficients, visibility,
channels; response, output, summation)` prepares channel synthesis from
calibrated path-envelope matrices already sampled on one common physical
detector grid. Envelopes are ordered path-major, then anchor-minor. The real
interpolation basis, complex path coefficients, and one unit-diagonal Hermitian
positive-semidefinite coherence matrix per channel are explicit inputs;
`channels` supplies normalized photon-number weights. Optional `response`
multiplies those weights without renormalization. `:coherent_modes` is the
default nonnegative summation; `:path_pairs` is the reference arithmetic and
can produce small negative values from finite-precision cancellation. This
operation performs no propagation, resampling, or inferred photometry.

The qualified advanced API includes `SpectralSynthesisPlan`,
`SpectralSynthesisWorkspace`, `SpectralSynthesisProducts`, and
`PreparedSpectralSynthesis`. `synthesize_spectrum!` executes on the caller's
current stream; the caller owns completion before reading output or changing
bound inputs.

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
The prepared batch and its count, capacity, metadata, and output accessors are
supported; its built-in extraction capability and scratch owners are internal.

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

`emccd_mean(signal_electrons, sensor; dark_electrons=0,
cic_electrons=sensor.clock_induced_charge_per_frame, gain=1,
poisson_cutoff=64)` returns `mean_electrons` and
`truncation_bound_electrons` for independent Poisson charge and the configured
linear EMCCD multiplication model. It includes the clipped-Gaussian model's
positive clipping bias; Gamma and conventional output means are exact. This
scalar statistic excludes saturation, read noise and ADC. Its Poisson
truncation bound is separate from floating-point and RNG approximation error.

`pixel_mtf_transfer(shape; box_width_pixels, diffusion_variance_pixels2,
ipc_per_neighbor, T)` builds an explicit corner-origin, unit-DC frequency
response. `prepare_pixel_mtf(input, transfer; padding, output)` binds that
calibration and exact input/output storage; `apply_pixel_mtf!` applies it with
prepared FFT workspace. Apply the deterministic response before acquisition and
avoid repeating response effects already represented by the calibration. It
does not change any sensor default or define detector noise covariance. The
qualified advanced owner types are `PixelMTFPlan` and `PreparedPixelMTF`.

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
`RelativeIlluminationValidSubapertures(peak_fraction=0.5)` selects
Shack–Hartmann layout cells relative to peak illumination.

The AOS Shack–Hartmann, Pyramid, Bi-O-edge, and Zernike surfaces implement
stages 1 and 2 only.
Their complete detector observations are package boundaries; maintained
operational estimation is provided by FilterGraphAlgorithms. Registered
FilterGraphAlgorithms provides the Zernike complete-frame normalized
pupil-signal estimator and the Curvature paired-image and paired-channel
estimators; JuliaFilterGraph provides their graph runtime. AOS Curvature
ends at complete detector acquisition. AOS LiFT provides physical forward
prediction and an analytic Jacobian callback; registered
`AdaptiveOpticsCalibration.PhaseRetrieval` owns inverse policies, solves,
products, and diagnostics through its direct API. LiFT is a phase-retrieval
workflow, not a slope or centroid estimator. AOS has no
weighted-centroid estimator or graph surface. The explicitly named
`geometric_wavefront_slopes!` operation remains available for plant truth and
reference calculations, not as a detector-derived RTC measurement.

The AOS `LiFTObservation` retains its AOS metadata and observation contract.
`PhaseRetrieval.LiFTSpecification(forward, observation)` checks that the
observation belongs to the physical forward contract and transfers its domain
and read-noise metadata into the AOC specification. Callers provide
`PhaseRetrieval.LiFTInputs(observation.values)` to AdaptiveOpticsCalibration.
Direct construction from `LiFTForwardModel(forward)` remains available when
the caller owns observation compatibility validation. See
[`examples/tutorials/lift.jl`](../examples/tutorials/lift.jl) for the registered
CPU direct API workflow. Accelerator execution requires an explicit
`AdaptiveOpticsCalibration.KernelExecution` policy.

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
materialization, and identification workflows.
Its qualified-public `ModalOPDExpansionPlan` and `combine_basis!` form the
advanced prepared plant-graph executor; routine users construct modal OPD
products through the calibration workflow instead of importing these names.
Reusable inverse methods and compact-SVD products come from
`AdaptiveOpticsCalibration.Reconstructors`. Common AOS entry points include
`interaction_matrix` and `modal_basis`. Compose modal-basis construction,
plant response acquisition, and inverse construction explicitly so the AOC
numerical products remain separate from AOS plant behavior.
For `interaction_matrix` and `compute_meta_sensitivity_matrix`, pass a typed
caller-owned `WFSMeasurement` vector and a callable
`(measurement, pupil)` that writes its exact storage. The callable owns the
chosen optical, detector, source, and estimator composition.

`KarhunenLoeveBasis` is owned by
`AdaptiveOpticsCalibration.ModalBases`. AOS supplies the sampled physical DM,
pupil measure, and atmospheric OPD covariance, then materializes the accepted
calibration product on the plant's runtime backend.

For NCPA, acquire a physical sampled basis in AOS, pass that basis, explicit
modal coefficients, and pupil support to
`AdaptiveOpticsCalibration.ModalBases.ModalOPDExpansion`, then construct
`Optics.NCPA` from the returned OPD matrix. The numerical product does not
choose an optical branch or normalize a stochastic coefficient profile;
the composing workflow owns those choices. See
[`examples/tutorials/ncpa.jl`](../examples/tutorials/ncpa.jl).

The maintained cross-package reference is
[`examples/integrations/filter_graph_algorithms/`](../examples/integrations/filter_graph_algorithms/),
where AOS supplies the plant and FilterGraphAlgorithms/JuliaFilterGraph supplies
the RTC chain. The composing fixture is CPU-resident and retains concrete
floating-point element types and caller-owned arrays for allocation-free
execution. AOS and FGA qualify their accelerator paths independently; the
fixture does not claim shared device storage or one captured graph across
packages.

For Pyramid complete-image modal optical gains, AOS supplies the physical
modulation-cycle image with `pyramid_modulation_frame` or
`pyramid_modulation_frame!` and exposes its matching sampled complex
transmission through `pyramid_focal_mask`. The returned mask is borrowed; pass
it, an explicitly aligned modal basis, and a reference image to
`AdaptiveOpticsCalibration.OpticalGains.GainSensingSpecification`, then
prepare `GainSensing()` and process each complete current frame through
`GainSensingInputs`. AOC copies the calibration inputs during preparation. AOS
does not provide a gain-sensing camera, calibration state, or optical-gain
estimator.

## Tomography

~~~julia
using AdaptiveOpticsSim.Tomography
~~~

Tomography owns guide-star and layer geometry, model-based or
interaction-matrix reconstruction, fitting, and cold DM command-matrix
projection. Use `build_reconstructor` and
`assemble_reconstructor_and_fitting` to prepare matrices for an external RTC.
FilterGraphAlgorithms/JuliaFilterGraph applies these matrices to incoming
measurements and produces per-frame commands. CPU covariance reconstruction
uses AdaptiveOpticsCalibration's explicit-array numerical solve; AOS retains
the physical covariance construction, wavelength conversion, slope ordering,
sampled DM influence construction, and command-matrix composition. The cold
CPU fitting pseudoinverse also uses AdaptiveOpticsCalibration's compact-SVD
inverse; AOS retains the physical `TomographyFitting` adapter and its legacy
threshold convention. The GPU covariance build path retains its existing solve
pending separate accelerator qualification.

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

Use `load_algorithm_graph(path; bindings=...)` to compile a versioned YAML graph
(after `using YAML` activates the optional extension)
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
- `ingaas_detector_acquisition_node` for complete-frame acquisition on an
  explicit metric pixel grid and integrated spectral channel
- `pixel_mtf_node` for a caller-calibrated response on exact image buffers
- `modal_opd_expansion_node`
- `second_order_response_node` for an explicit unity-static-gain mechanical
  response before modal expansion or other physical coordinate use

`Optics.SecondOrderResponsePlan(; natural_frequency_hz, damping_ratio,
sample_period_s, T=Float32)` prepares exact zero-order-hold coefficients for
independent coordinates. `prepare_second_order_response(input, plan; output)`
binds distinct same-precision vectors and zero-initialized persistent position
and velocity. `advance_response!` emits the end-of-hold position;
`reset_response!` restores rest without replacing storage. The graph adapter's
`:input` and `:output` use one caller-declared `schema` and `count`. Applications
must match `sample_period_s` to their graph cadence and declare whether chosen
frequency/damping parameters are provisional or measured. No instrument values,
extra pure delay, or controller filter are supplied. Host stream and CUDA
captured execution are covered by focused response tests.

`ingaas_detector_acquisition_node` declares one complete-frame, global-shutter
InGaAs acquisition. The photon-rate input is cell-integrated photons per
second on a metric detector grid, and `spectral_channel` explicitly identifies
the already integrated passband. The node uses scalar effective-band QE; if an
upstream spectral response already includes QE, set detector QE to one. It
preserves the frame axes and applies its declared exposure, noise, gain, dark
current, glow, saturation, and optional quantization. Pixel MTF remains a
separate upstream operation. Captured execution is restricted to the admitted
full-frame null-persistence and null-nonlinearity path.

`hgcdte_avalanche_detector_acquisition_node` uses the same metric photon-rate
and frame contracts for `HgCdTeAvalancheArraySensor`, including caller-supplied
SAPHIRA configurations. It prepares the existing numerical detector with
single-read sampling. `avalanche_gain` multiplies charge before read noise;
`gain` is the separate post-readout conversion gain. `full_well_e` bounds
multiplied charge. `excess_noise_factor` uses F = E[M²]/E[M]², with F = 1
deterministic. `multiplication_model` accepts the CPU conditional-Gamma reference
or the clipped-Gaussian approximation described in the glossary. No empirical
camera defaults are supplied. CPU stream execution is covered; accelerator
capture qualification is unchanged.

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

Dense and measured influence-matrix deformable-mirror nodes also admit capture
with the built-in linear-static, clipping, and actuator-health stages or their
static compositions. Preparation retains the multiplication coefficients for
replay. Stateful or custom actuator stages need separate qualification.

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

For calibration runs that hold one probe over several exposures, use
`prepare_graph_calibration_boundary`, `adopt_hil_probe!`,
`step_hil_exposure!`, and `step_hil_exposure_at!`. The probe sequence is a
positive, strictly advancing `UInt64`, separate from graph exposure sequence.
Each accepted exposure advances the graph and optional model-time driver once;
writing the host command buffer does not adopt a new probe. The caller owns
transport completion, settling, exposure duration, and response aggregation.
`PreparedGraphCalibrationBoundary` is the supported advanced owner type.

## Optional Integrations

Optional GPU packages extend `Backends` and covered scientific operations.
Proper.jl is application-owned; see
[`proper-integration-guide.md`](proper-integration-guide.md). Plotting belongs
in the companion `AdaptiveOpticsSimPlots.jl` package.
