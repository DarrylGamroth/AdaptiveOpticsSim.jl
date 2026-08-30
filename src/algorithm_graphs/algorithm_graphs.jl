"""
    AlgorithmGraphs

Portable, deterministic composition of AOS-native algorithm nodes. Scientific
modules retain their domain APIs; this module owns only the graph-node adapter
protocol, static topology, exact buffer bindings, delayed links, and model-time
execution policy.
"""
module AlgorithmGraphs

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault
using Random: seed!
using TOML

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError, runtime_rng
import ..Backends:
    AbstractComputeDevice,
    HostComputeDevice,
    _prepare_device_execution_context,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    allocate_device_array,
    array_backend_type,
    compute_device,
    compute_device_backend
import ..Plant:
    PeriodicSchedule,
    PlantDuration,
    PlantTimestamp,
    schedule_period,
    schedule_timestamp
import ..Calibration:
    ModalOPDExpansionPlan,
    combine_basis!
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
    AngularCoordinates,
    CellIntegratedMeasure,
    DetectorPlane,
    IncoherentIntensityAddition,
    IntensityMap,
    MonochromaticChannel,
    NormalizedPupilCoordinates,
    OpticalPlaneMetadata,
    PhotonRateNormalization,
    PupilFunction,
    Source,
    TelescopeDefinition,
    photon_irradiance,
    prepare_telescope
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
    PyramidOpticalFrontEnd,
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
include("api.jl")

end # module AlgorithmGraphs
