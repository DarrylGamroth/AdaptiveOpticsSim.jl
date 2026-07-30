"""
    WavefrontSensors

Canonical owner of composed wavefront-sensor contracts, typed observations
and measurements, prepared stage protocols, and family-neutral execution
capabilities. Reusable physical WFS components remain owned by `Optics`, while
detector response and acquisition physics remain owned by `Detectors`.
"""
module WavefrontSensors

using KernelAbstractions
using Random

import ..AdaptiveOpticsSim:
    AdaptiveOpticsSimError,
    AbstractDetector,
    InvalidConfiguration

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    ExecutionStyle,
    HostComputeDevice,
    ScalarCPUStyle,
    backend,
    compute_device,
    execution_style,
    launch_kernel!

import ..Optics:
    AbstractCombinationPolicy,
    AbstractOpticalElement,
    AbstractOpticalNormalization,
    AbstractOpticalPlaneKind,
    AbstractOpticalProduct,
    AbstractSource,
    AbstractSpatialMeasure,
    AbstractSpectralCoordinate,
    Asterism,
    CellIntegratedMeasure,
    DetectorPlane,
    ElectricField,
    IncoherentIntensityAddition,
    IntegratedSpectralChannel,
    IntensityMap,
    MonochromaticChannel,
    OpticalPlaneMetadata,
    OpticalProductBundle,
    PhotonRateNormalization,
    PupilFunction,
    PupilPlane,
    SpatialDensityMeasure,
    photon_irradiance,
    require_leaf_source,
    validate_plane_storage,
    wavelength

import ..Detectors:
    AbstractCountingDetector,
    CCDSensor,
    CMOSSensor,
    Detector,
    EMCCDSensor,
    HgCdTeAvalancheArraySensor,
    InGaAsSensor,
    MKIDArrayDetector,
    capture!,
    counting_array,
    counting_exposure_time,
    counting_fill_factor,
    counting_post_gain,
    counting_qe,
    counting_source_throughput,
    effective_qe,
    ensure_buffers!,
    output_frame,
    prepare_detector_acquisition

include("abstract_wfs.jl")
include("sensing_modes.jl")
include("stage_contracts.jl")
include("grouped.jl")
include("interface.jl")
include("api.jl")

end # module WavefrontSensors
