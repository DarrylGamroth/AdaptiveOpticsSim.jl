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
    DiscreteIntegratorPlan,
    DiscreteIntegratorState,
    DiscreteIntegratorWorkspace,
    reset_controller!,
    update!
import ..Optics:
    AngularCoordinates,
    CellIntegratedMeasure,
    DetectorPlane,
    IncoherentIntensityAddition,
    IntensityMap,
    MonochromaticChannel,
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
    ShackHartmannWFS,
    form_wfs_optical_products!,
    prepare_wfs_optics,
    shack_hartmann_optics,
    shack_hartmann_rate_map

include("definitions.jl")
include("preparation.jl")
include("native_nodes.jl")
include("graph_files.jl")
include("execution.jl")
include("model_time.jl")
include("api.jl")

end # module AlgorithmGraphs
