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
using TOML

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError
import ..Backends:
    AbstractComputeDevice,
    AcceleratorComputeDevice,
    AcceleratorStyle,
    HostComputeDevice,
    ScalarCPUStyle,
    _prepare_device_execution_context,
    _prepare_graph_rng,
    _capture_prepared_device_graph,
    _launch_prepared_device_graph!,
    _reset_graph_rng!,
    _synchronize_prepared_device_execution_context!,
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
import ..Control:
    ClosedLoopCorrectionPlan,
    ClosedLoopCorrectionState,
    ClosedLoopCorrectionWorkspace,
    ControlMatrixPlan,
    DiscreteIntegratorPlan,
    DiscreteIntegratorState,
    DiscreteIntegratorWorkspace,
    apply_closed_loop_correction!,
    reconstruct!,
    reset_closed_loop_correction!,
    reset_controller!,
    update!
import ..Optics:
    ActuatorGridTopology,
    AngularCoordinates,
    CellIntegratedMeasure,
    DeformableMirror,
    DetectorPlane,
    GaussianInfluenceWidth,
    IncoherentIntensityAddition,
    IntensityMap,
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
    capture!,
    detector_acquisition_detector,
    detector_acquisition_state,
    output_frame,
    prepare_detector_acquisition,
    reset_integration!
import ..WavefrontSensors:
    Diffractive,
    AbstractPyramidModulationPropagationStrategy,
    PyramidOpticalFrontEnd,
    PyramidPupilTiltStrategy,
    PyramidShiftedMaskStrategy,
    PyramidWFS,
    ShackHartmannSlopeSelectionPlan,
    ShackHartmannWFS,
    WFSMeasurement,
    WFSObservation,
    estimate_wfs_measurement!,
    form_wfs_optical_products!,
    prepare_wfs_estimation,
    prepare_wfs_optics,
    pyramid_rate_map,
    select_shack_hartmann_slopes!,
    set_subaperture_calibration!,
    set_valid_subapertures!,
    shack_hartmann_optics,
    shack_hartmann_rate_map,
    subaperture_calibration,
    selected_lenslet_count

include("definitions.jl")
include("model_time_values.jl")
include("preparation.jl")
include("native_nodes.jl")
include("captured_execution.jl")
include("graph_files.jl")
include("execution.jl")
include("model_time.jl")
include("hil_boundary.jl")
include("api.jl")

end # module AlgorithmGraphs
