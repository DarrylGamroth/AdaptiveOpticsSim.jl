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
using Random: seed!
using TOML

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError, runtime_rng
import ..Backends:
    AbstractComputeDevice,
    AcceleratorStyle,
    HostComputeDevice,
    ScalarCPUStyle,
    _prepare_device_execution_context,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    allocate_device_array,
    array_backend_type,
    compute_device,
    compute_device_backend,
    execution_style,
    launch_kernel!
import ..Plant:
    PeriodicSchedule,
    PlantDuration,
    PlantTimestamp,
    schedule_period,
    schedule_timestamp
import ..Calibration:
    ModalOPDExpansionPlan,
    combine_basis!
import ..Atmospheres:
    MultiLayerAtmosphereDefinition,
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
    command_storage,
    photon_irradiance,
    prepare_telescope,
    set_command!,
    surface_opd,
    update_surface!
import ..Detectors:
    CCDSensor,
    Detector,
    EMCCDSensor,
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
include("preparation.jl")
include("native_nodes.jl")
include("graph_files.jl")
include("execution.jl")
include("model_time.jl")
include("hil_boundary.jl")
include("api.jl")

end # module AlgorithmGraphs
