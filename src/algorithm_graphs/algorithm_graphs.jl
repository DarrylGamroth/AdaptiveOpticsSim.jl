"""
    AlgorithmGraphs

Portable, deterministic composition of AOS-native algorithm nodes. Scientific
modules retain their domain APIs; this module owns only the graph-node adapter
protocol, static topology, exact buffer bindings, delayed links, and model-time
execution policy.
"""
module AlgorithmGraphs

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault
using KernelAbstractions
using SHA: sha256
import ..Backends

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError
import ..Backends:
    AbstractComputeDevice,
    AcceleratorComputeDevice,
    AcceleratorStyle,
    HostComputeDevice,
    ScalarCPUStyle,
    _prepare_device_execution_context,
    _prepare_device_execution_event,
    _prepare_graph_rng,
    _prepare_mul_coefficients,
    _capture_prepared_device_graph,
    _launch_prepared_device_graph!,
    _record_prepared_device_execution_event!,
    _reset_graph_rng!,
    _synchronize_prepared_device_execution_context!,
    _synchronize_prepared_device_execution_context_blocking!,
    _wait_prepared_device_execution_event!,
    _with_prepared_device_execution_context,
    allocate_device_array,
    array_backend_type,
    compute_device,
    compute_device_backend,
    copyto_backend_async!,
    execution_style,
    launch_kernel!,
    launch_kernel_async!
import ..Calibration:
    ModalOPDExpansionPlan,
    combine_basis!
import ..Atmospheres:
    MultiLayerAtmosphereDefinition,
    SpatialAORejection,
    _prepare_spatial_ao_rejection,
    _apply_spatial_ao_rejection!,
    _enqueue_moving_screen_replay!,
    _preflight_atmosphere_replay_step,
    _prepare_moving_screen_replay,
    _render_atmosphere_async!,
    _reset_multilayer_atmosphere!,
    _reset_moving_screen_replay!,
    advance_by!,
    prepare_atmosphere_renderer,
    prepare_timed_atmosphere,
    render_atmosphere!
import ..Optics:
    ActuatorGridTopology,
    influence_model,
    MetricCoordinates,
    IntegratedSpectralChannel,
    DenseInfluenceMatrix,
    CompositeDMActuatorModel,
    ClippedActuators,
    ActuatorHealthMap,
    AbstractDMActuatorModel,
    AngularCoordinates,
    CellIntegratedMeasure,
    DeformableMirror,
    GaussianInfluenceOperator,
    DetectorPlane,
    GaussianInfluenceWidth,
    IncoherentIntensityAddition,
    IntensityMap,
    LinearStaticActuators,
    MeasuredInfluenceFunctions,
    MonochromaticChannel,
    NormalizedPupilCoordinates,
    OpticalPlaneMetadata,
    PhotonRateNormalization,
    PupilFunction,
    SampledActuatorTopology,
    Source,
    TelescopeDefinition,
    _update_dm_surface_async!,
    actuator_model,
    command_storage,
    photon_irradiance,
    prepare_telescope,
    set_command!,
    surface_opd,
    update_surface!
import ..Detectors:
    CCDSensor,
    CMOSSensor,
    Detector,
    EMCCDSensor,
    GlobalShutter,
    NoiseModel,
    NoiseNone,
    NoisePhoton,
    NoisePhotonReadout,
    NoiseReadout,
    NullFrameResponse,
    SaturatingFrameNonlinearity,
    NullPersistence,
    NullFrameNonlinearity,
    InGaAsSensor,
    HgCdTeAvalancheArraySensor,
    ConditionalGammaAvalancheMultiplication,
    ClippedGaussianAvalancheMultiplicationApproximation,
    ExponentialPersistence,
    _complete_enqueued_detector_acquisition!,
    _enqueue_detector_acquisition!,
    _preflight_enqueued_detector_acquisition!,
    _supports_enqueued_detector_acquisition,
    capture!,
    detector_acquisition_detector,
    detector_acquisition_state,
    output_frame,
    prepare_detector_acquisition,
    reset_integration!
import ..WavefrontSensors:
    _pyramid_face_shifts_pixels,
    Diffractive,
    AbstractPyramidModulationPropagationStrategy,
    NoPreparedFourPupilLGS,
    PyramidModulationBatchWorkspace,
    PyramidOpticalFrontEnd,
    PyramidPupilTiltStrategy,
    PyramidWFS,
    ShackHartmannWFS,
    enqueue_wfs_optical_products!,
    form_wfs_optical_products!,
    prepare_wfs_optics,
    pyramid_rate_map,
    shack_hartmann_optics,
    shack_hartmann_rate_map

include("definitions.jl")
include("model_time_values.jl")
include("preparation.jl")
include("native_nodes.jl")
include("second_order_response_node.jl")
include("ingaas_detector_node.jl")
include("hgcdte_avalanche_detector_node.jl")
include("pixel_mtf_node.jl")
include("captured_execution.jl")
include("proper_node.jl")
include("ideal_pupil_relay.jl")
include("graph_files.jl")
include("execution.jl")
include("model_time.jl")
include("hil_boundary.jl")
include("calibration_boundary.jl")
include("output_snapshot.jl")
include("external_pupil.jl")

include("api.jl")

end # module AlgorithmGraphs
