module AOSFGALockstep

using AdaptiveOpticsCalibration
using AdaptiveOpticsCalibration.Reconstructors
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using FilterGraphAlgorithms
using JuliaFilterGraph
using LinearAlgebra
using Random

export DETECTOR_FRAME_SCHEMA
export S1DetectorFrame, S1FrameDisposition, FrameAccepted
export FrameRejectedBlocked, FrameRejectedNoOutstandingFrame
export FrameRejectedSequence, FrameRejectedTimestamp, FrameRejectedIncomplete
export FrameRejectedDiscontinuous, FrameRejectedCorrupted
export FrameRejectedSchema, FrameRejectedCalibrationSignature
export FrameRejectedNumericType, FrameRejectedShape, FrameRejectedNonFinite
export FrameRejectedRTCFailure, FrameRejectedRTCUnavailable
export FrameRejectedRTCMetadata, FrameRejectedCommand
export prepare_s1_lockstep, produce_detector_frame!, process_detector_frame!
export step_lockstep!, reset_s1_lockstep!

const DETECTOR_FRAME_SCHEMA =
    "org.adaptiveopticssim.integration.shack-hartmann-electron-counts-row-column.f32/1"
const _RNG_SEED = UInt64(0x5331)
const _FRAME_PERIOD_NANOSECONDS = Int64(1_000_000)
const _SIGNATURE_OFFSET = UInt64(0xcbf29ce484222325)
const _SIGNATURE_PRIME = UInt64(0x00000100000001b3)

const _TELESCOPE_RESOLUTION = 8
const _TELESCOPE_DIAMETER_M = 1.0f0
const _CENTRAL_OBSTRUCTION_RATIO = 0.0f0
const _TELESCOPE_FOV_ARCSEC = 0.0f0
const _TELESCOPE_REFLECTIVITY = 1.0f0
const _SOURCE_BAND = :custom
const _SOURCE_MAGNITUDE = 0.0f0
const _SOURCE_COORDINATES_ARCSEC_DEG = (0.0f0, 0.0f0)
const _SOURCE_WAVELENGTH_M = 750.0f-9
const _SOURCE_PHOTON_IRRADIANCE_M2_S = 2.0f8
const _SOURCE_RADIOMETRY = PhysicalPhotonIrradianceSource()
const _ACTUATOR_COORDINATES = ((0.0f0, 0.0f0),)
const _DM_INFLUENCE_WIDTH = 0.35f0
const _SHACK_HARTMANN_LENSLETS = 2
const _SHACK_HARTMANN_PIXELS_PER_SUBAPERTURE = 4
const _SHACK_HARTMANN_VALID_THRESHOLD = 0.1f0
const _SHACK_HARTMANN_COG_THRESHOLD = 0.01f0
const _SHACK_HARTMANN_CONVOLUTION_THRESHOLD = 0.05f0
const _SHACK_HARTMANN_HALF_PIXEL_SHIFT = false
const _SHACK_HARTMANN_DIFFRACTION_PADDING = 2
const _SHACK_HARTMANN_PIXEL_SCALE = nothing
const _SHACK_HARTMANN_SHANNON_SAMPLING = true
const _DETECTOR_EXPOSURE_S = 1.0f-3
const _DETECTOR_QUANTUM_EFFICIENCY = 1.0f0
const _WFS_FORMATION_MODEL = :diffractive_shack_hartmann
const _DETECTOR_NOISE_MODEL = :none
const _DETECTOR_RESPONSE_MODEL = :null_frame_response
const _CALIBRATION_POKE_M = 2.0f-8
const _DETECTOR_AXES = (:x, :y)
const _ESTIMATOR_FRAME_AXES = (:row, :column)
const _SLOPE_PAIR_ORDER = (:x, :y)
const _PDM_ACTUATOR_ORDER = (1,)
const _DETECTOR_UNITS = :electron_count
const _SLOPE_UNITS = :pixel
const _PDM_COMMAND_UNITS = :metre
const _OBSERVATION_LAYOUT = :lenslet_mosaic

@enum S1FrameDisposition::UInt8 begin
    FrameAccepted = 0
    FrameRejectedBlocked = 1
    FrameRejectedNoOutstandingFrame = 2
    FrameRejectedSequence = 3
    FrameRejectedTimestamp = 4
    FrameRejectedIncomplete = 5
    FrameRejectedDiscontinuous = 6
    FrameRejectedCorrupted = 7
    FrameRejectedSchema = 8
    FrameRejectedCalibrationSignature = 9
    FrameRejectedNumericType = 10
    FrameRejectedShape = 11
    FrameRejectedNonFinite = 12
    FrameRejectedRTCFailure = 13
    FrameRejectedRTCUnavailable = 14
    FrameRejectedRTCMetadata = 15
    FrameRejectedCommand = 16
end

"""Complete AOS detector frame plus the identity facts required by S1."""
struct S1DetectorFrame{A<:AbstractMatrix}
    values::A
    schema::String
    calibration_signature::UInt64
    sequence::UInt64
    timestamp_nanoseconds::Int64
    terminal::Bool
    discontinuity::Bool
    corrupted::Bool
end

"""Exact cross-package interpretation of the cold S1 calibration products."""
struct S1CalibrationIdentity{T<:AbstractFloat}
    detector_axes::NTuple{2,Symbol}
    estimator_frame_axes::NTuple{2,Symbol}
    subaperture_order::NTuple{4,NTuple{2,Int}}
    slope_pair_order::NTuple{2,Symbol}
    pdm_actuator_order::NTuple{1,Int}
    detector_units::Symbol
    slope_units::Symbol
    pdm_command_units::Symbol
    numeric_type::Type{T}
    plant_signature::UInt64
    estimator_signature::UInt64
    signature::UInt64
end

"""AOC numerical result bound to the exact S1 plant and RTC interpretation."""
struct S1CalibrationProduct{T<:AbstractFloat,I,R}
    identity::S1CalibrationIdentity{T}
    interaction_matrix::I
    reconstructor_product::R
end

mutable struct S1LockstepState
    sequence::UInt64
    outstanding::Bool
    blocked::Bool
end

struct PreparedS1Lockstep{
    Pupil,
    DM,
    Rate,
    Observation,
    OpticsPlan,
    AcquisitionPlan,
    RNG,
    Graph,
    Outputs,
    Inputs,
    Calibration,
}
    pupil::Pupil
    dm::DM
    rate::Rate
    observation::Observation
    optics_plan::OpticsPlan
    acquisition_plan::AcquisitionPlan
    rng::RNG
    graph::Graph
    outputs::Outputs
    inputs::Inputs
    calibration::Calibration
    disturbance_opd::Matrix{Float32}
    fga_frame::Matrix{Float32}
    adopted_command::Vector{Float32}
    applied_command::Vector{Float32}
    state::S1LockstepState
end

# This is a deterministic configuration fingerprint for invalidation and
# compatibility checks, not a cryptographic integrity mechanism.
@inline function _signature_byte(signature::UInt64, byte::UInt8)
    return xor(signature, UInt64(byte)) * _SIGNATURE_PRIME
end

function _signature_bytes(signature::UInt64, bytes)
    value = signature
    @inbounds for byte in bytes
        value = _signature_byte(value, UInt8(byte))
    end
    return value
end

function _signature_blob(signature::UInt64, bytes)
    result = _signature_uint64(signature, UInt64(length(bytes)))
    return _signature_bytes(result, bytes)
end

@inline function _signature_uint64(signature::UInt64, value::UInt64)
    result = signature
    @inbounds for shift in 0:8:56
        result = _signature_byte(result, UInt8((value >> shift) & 0xff))
    end
    return result
end

@inline _signature_integer(signature::UInt64, value::Integer) =
    _signature_uint64(signature, reinterpret(UInt64, Int64(value)))

@inline _signature_float32(signature::UInt64, value::Float32) =
    _signature_uint64(signature, UInt64(reinterpret(UInt32, value)))

@inline _signature_symbol(signature::UInt64, value::Symbol) =
    _signature_blob(signature, codeunits(String(value)))

@inline _signature_value(signature::UInt64, ::Nothing) =
    _signature_blob(signature, codeunits("Nothing"))

@inline function _signature_value(signature::UInt64, value::Bool)
    result = _signature_blob(signature, codeunits("Bool"))
    return _signature_byte(result, UInt8(value))
end

@inline function _signature_value(signature::UInt64, value::Integer)
    result = _signature_blob(signature, codeunits(string(typeof(value))))
    return _signature_blob(result, codeunits(string(value)))
end

@inline function _signature_value(signature::UInt64, value::AbstractFloat)
    result = _signature_blob(signature, codeunits(string(typeof(value))))
    return _signature_blob(result, codeunits(bitstring(value)))
end

@inline function _signature_value(signature::UInt64, value::Symbol)
    result = _signature_blob(signature, codeunits("Symbol"))
    return _signature_symbol(result, value)
end

@inline function _signature_value(signature::UInt64, value::Type)
    result = _signature_blob(signature, codeunits("Type"))
    return _signature_blob(result, codeunits(string(value)))
end

function _signature_value(signature::UInt64, values::Tuple)
    result = _signature_blob(signature, codeunits(string(typeof(values))))
    result = _signature_integer(result, length(values))
    for value in values
        result = _signature_value(result, value)
    end
    return result
end

function _signature_value(signature::UInt64, values::AbstractArray)
    result = _signature_blob(signature, codeunits(string(typeof(values))))
    result = _signature_integer(result, ndims(values))
    for extent in size(values)
        result = _signature_integer(result, extent)
    end
    @inbounds for value in values
        result = _signature_value(result, value)
    end
    return result
end

function _signature_value(signature::UInt64, value)
    value_type = typeof(value)
    isstructtype(value_type) || error(
        "unsupported S1 signature value type $value_type",
    )
    result = _signature_blob(signature, codeunits(string(value_type)))
    for name in fieldnames(value_type)
        result = _signature_symbol(result, name)
        result = _signature_value(result, getfield(value, name))
    end
    return result
end

function _signature_float32_array(signature::UInt64, values)
    result = _signature_integer(signature, ndims(values))
    for extent in size(values)
        result = _signature_integer(result, extent)
    end
    @inbounds for value in values
        result = _signature_float32(result, Float32(value))
    end
    return result
end

function _plant_signature(;
    telescope_resolution=_TELESCOPE_RESOLUTION,
    telescope_diameter_m=_TELESCOPE_DIAMETER_M,
    central_obstruction_ratio=_CENTRAL_OBSTRUCTION_RATIO,
    telescope_fov_arcsec=_TELESCOPE_FOV_ARCSEC,
    telescope_reflectivity=_TELESCOPE_REFLECTIVITY,
    source_band=_SOURCE_BAND,
    source_magnitude=_SOURCE_MAGNITUDE,
    source_coordinates_arcsec_deg=_SOURCE_COORDINATES_ARCSEC_DEG,
    source_wavelength_m=_SOURCE_WAVELENGTH_M,
    source_photon_irradiance_m2_s=_SOURCE_PHOTON_IRRADIANCE_M2_S,
    source_radiometry=_SOURCE_RADIOMETRY,
    actuator_coordinates=_ACTUATOR_COORDINATES,
    dm_influence_width=_DM_INFLUENCE_WIDTH,
    shack_hartmann_lenslets=_SHACK_HARTMANN_LENSLETS,
    shack_hartmann_pixels_per_subaperture=
        _SHACK_HARTMANN_PIXELS_PER_SUBAPERTURE,
    shack_hartmann_valid_threshold=_SHACK_HARTMANN_VALID_THRESHOLD,
    shack_hartmann_cog_threshold=_SHACK_HARTMANN_COG_THRESHOLD,
    shack_hartmann_convolution_threshold=
        _SHACK_HARTMANN_CONVOLUTION_THRESHOLD,
    shack_hartmann_half_pixel_shift=
        _SHACK_HARTMANN_HALF_PIXEL_SHIFT,
    shack_hartmann_diffraction_padding=
        _SHACK_HARTMANN_DIFFRACTION_PADDING,
    shack_hartmann_pixel_scale=_SHACK_HARTMANN_PIXEL_SCALE,
    shack_hartmann_shannon_sampling=
        _SHACK_HARTMANN_SHANNON_SAMPLING,
    detector_exposure_s=_DETECTOR_EXPOSURE_S,
    detector_quantum_efficiency=_DETECTOR_QUANTUM_EFFICIENCY,
    detector_units=_DETECTOR_UNITS,
    observation_layout=_OBSERVATION_LAYOUT,
    detector_metadata=nothing,
)
    signature =
        _signature_blob(_SIGNATURE_OFFSET, codeunits("AOS-S1-PLANT/2"))
    signature = _signature_integer(signature, telescope_resolution)
    signature = _signature_float32(signature, telescope_diameter_m)
    signature = _signature_float32(signature, central_obstruction_ratio)
    signature = _signature_float32(signature, telescope_fov_arcsec)
    signature = _signature_float32(signature, telescope_reflectivity)
    signature = _signature_symbol(signature, source_band)
    signature = _signature_float32(signature, source_magnitude)
    signature = _signature_value(signature, source_coordinates_arcsec_deg)
    signature = _signature_float32(signature, source_wavelength_m)
    signature = _signature_float32(signature, source_photon_irradiance_m2_s)
    signature = _signature_value(signature, source_radiometry)
    signature = _signature_value(signature, actuator_coordinates)
    signature = _signature_float32(signature, dm_influence_width)
    signature = _signature_integer(signature, shack_hartmann_lenslets)
    signature = _signature_integer(
        signature,
        shack_hartmann_pixels_per_subaperture,
    )
    signature = _signature_float32(signature, shack_hartmann_valid_threshold)
    signature = _signature_float32(signature, shack_hartmann_cog_threshold)
    signature = _signature_float32(
        signature,
        shack_hartmann_convolution_threshold,
    )
    signature = _signature_value(signature, shack_hartmann_half_pixel_shift)
    signature = _signature_integer(
        signature,
        shack_hartmann_diffraction_padding,
    )
    signature = _signature_value(signature, shack_hartmann_pixel_scale)
    signature = _signature_value(signature, shack_hartmann_shannon_sampling)
    signature = _signature_float32(signature, detector_exposure_s)
    signature = _signature_float32(signature, detector_quantum_efficiency)
    signature = _signature_symbol(signature, _WFS_FORMATION_MODEL)
    signature = _signature_symbol(signature, _DETECTOR_NOISE_MODEL)
    signature = _signature_symbol(signature, _DETECTOR_RESPONSE_MODEL)
    signature = _signature_uint64(signature, _RNG_SEED)
    signature = _signature_symbol(signature, detector_units)
    signature = _signature_symbol(signature, observation_layout)
    signature = _signature_value(signature, detector_metadata)
    return signature
end

function _estimator_signature(
    graph_path,
    reference_slopes,
    reconstructor_matrix,
    controller_to_vdm,
    active_to_full_vdm,
    vdm_to_pdm,
    subaperture_order,
)
    signature = _signature_blob(
        _SIGNATURE_OFFSET,
        codeunits("FGA-S1-ESTIMATOR/2"),
    )
    signature = _signature_blob(signature, read(graph_path))
    signature = _signature_float32_array(signature, reference_slopes)
    signature = _signature_float32_array(signature, reconstructor_matrix)
    signature = _signature_float32_array(signature, controller_to_vdm)
    signature = _signature_float32_array(signature, active_to_full_vdm)
    signature = _signature_float32_array(signature, vdm_to_pdm)
    signature = _signature_value(signature, subaperture_order)
    signature = _signature_float32(signature, _CALIBRATION_POKE_M)
    return signature
end

function _calibration_identity(
    plant_signature,
    estimator_signature;
    detector_axes=_DETECTOR_AXES,
    estimator_frame_axes=_ESTIMATOR_FRAME_AXES,
    subaperture_order,
    slope_pair_order=_SLOPE_PAIR_ORDER,
    pdm_actuator_order=_PDM_ACTUATOR_ORDER,
    detector_units=_DETECTOR_UNITS,
    slope_units=_SLOPE_UNITS,
    pdm_command_units=_PDM_COMMAND_UNITS,
    numeric_type=Float32,
)
    signature = _signature_blob(
        _SIGNATURE_OFFSET,
        codeunits("AOS-FGA-S1-CALIBRATION/1"),
    )
    for axis in detector_axes
        signature = _signature_symbol(signature, axis)
    end
    for axis in estimator_frame_axes
        signature = _signature_symbol(signature, axis)
    end
    for origin in subaperture_order, coordinate in origin
        signature = _signature_integer(signature, coordinate)
    end
    for component in slope_pair_order
        signature = _signature_symbol(signature, component)
    end
    for actuator in pdm_actuator_order
        signature = _signature_integer(signature, actuator)
    end
    signature = _signature_symbol(signature, detector_units)
    signature = _signature_symbol(signature, slope_units)
    signature = _signature_symbol(signature, pdm_command_units)
    signature = _signature_blob(signature, codeunits(string(numeric_type)))
    signature = _signature_uint64(signature, plant_signature)
    signature = _signature_uint64(signature, estimator_signature)
    return S1CalibrationIdentity(
        detector_axes,
        estimator_frame_axes,
        subaperture_order,
        slope_pair_order,
        pdm_actuator_order,
        detector_units,
        slope_units,
        pdm_command_units,
        numeric_type,
        plant_signature,
        estimator_signature,
        signature,
    )
end

function _prepared_subaperture_order(graph)
    measure = first(graph.nodes)
    measure.name === :measure || error(
        "the first S1 FGA Node must be the Shack-Hartmann measurement Node",
    )
    origins = measure.prepared.plan.regions.origins
    length(origins) == 4 || error(
        "the prepared S1 measurement plan must contain four subapertures",
    )
    return (origins[1], origins[2], origins[3], origins[4])
end

@inline function _all_finite(values)
    @inbounds for value in values
        isfinite(value) || return false
    end
    return true
end

@inline function _interleave_slopes!(ordered, slopes)
    @inbounds for subaperture in axes(slopes, 1)
        ordered[2 * subaperture - 1] = slopes[subaperture, 1]
        ordered[2 * subaperture] = slopes[subaperture, 2]
    end
    return ordered
end

function _prepare_plant()
    T = Float32
    telescope = Telescope(
        resolution=_TELESCOPE_RESOLUTION,
        diameter=_TELESCOPE_DIAMETER_M,
        central_obstruction=_CENTRAL_OBSTRUCTION_RATIO,
        fov_arcsec=_TELESCOPE_FOV_ARCSEC,
        pupil_reflectivity=_TELESCOPE_REFLECTIVITY,
        T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=_SOURCE_BAND,
        magnitude=_SOURCE_MAGNITUDE,
        coordinates=_SOURCE_COORDINATES_ARCSEC_DEG,
        wavelength=_SOURCE_WAVELENGTH_M,
        photon_irradiance=_SOURCE_PHOTON_IRRADIANCE_M2_S,
        radiometry=_SOURCE_RADIOMETRY,
        T=T,
    )
    topology = SampledActuatorTopology(
        reshape(T[first(_ACTUATOR_COORDINATES)...], 2, 1);
        T=T,
    )
    dm = DeformableMirror(
        telescope;
        topology,
        influence_width=_DM_INFLUENCE_WIDTH,
        T=T,
    )
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=_SHACK_HARTMANN_LENSLETS,
        threshold=_SHACK_HARTMANN_VALID_THRESHOLD,
        threshold_cog=_SHACK_HARTMANN_COG_THRESHOLD,
        threshold_convolution=_SHACK_HARTMANN_CONVOLUTION_THRESHOLD,
        half_pixel_shift=_SHACK_HARTMANN_HALF_PIXEL_SHIFT,
        diffraction_padding=_SHACK_HARTMANN_DIFFRACTION_PADDING,
        pixel_scale_arcsec=_SHACK_HARTMANN_PIXEL_SCALE,
        n_pix_subap=_SHACK_HARTMANN_PIXELS_PER_SUBAPERTURE,
        shannon_sampling=_SHACK_HARTMANN_SHANNON_SAMPLING,
        mode=Diffractive(),
        T=T,
    )
    rate = shack_hartmann_rate_map(sensor, pupil, source)
    optics = shack_hartmann_optics(sensor, source)
    optics_plan = prepare_wfs_optics(optics, pupil, rate)
    detector = Detector(
        noise=NoiseNone(),
        exposure_duration=_DETECTOR_EXPOSURE_S,
        qe=_DETECTOR_QUANTUM_EFFICIENCY,
        response_model=NullFrameResponse(),
        T=T,
    )
    observation = WFSObservation(
        similar(intensity_values(rate));
        units=_DETECTOR_UNITS,
        layout=_OBSERVATION_LAYOUT,
    )
    acquisition_plan = prepare_wfs_acquisition(
        detector,
        rate,
        observation;
        source,
    )
    return (;
        pupil,
        dm,
        rate,
        observation,
        optics_plan,
        acquisition_plan,
        rng=Xoshiro(_RNG_SEED),
        plant_signature=_plant_signature(
            detector_metadata=detector_export_metadata(detector),
        ),
    )
end

@inline function _form_plant_frame!(prepared)
    copyto!(opd_map(prepared.pupil), prepared.disturbance_opd)
    copyto!(prepared.applied_command, prepared.adopted_command)
    update_surface!(prepared.dm)
    apply_surface!(prepared.pupil, prepared.dm, DMAdditive())
    form_wfs_optical_products!(
        prepared.rate,
        prepared.pupil,
        prepared.optics_plan,
    )
    acquire_wfs_observation!(
        prepared.observation,
        prepared.rate,
        prepared.acquisition_plan,
        prepared.rng,
    )
    return observation_storage(prepared.observation)
end

function _calibrate_reconstructor!(
    prepared,
    graph,
    outputs,
    inputs,
    graph_path,
    controller_to_vdm,
    active_to_full_vdm,
    vdm_to_pdm,
    subaperture_order,
)
    fill!(prepared.disturbance_opd, 0.0f0)
    fill!(prepared.adopted_command, 0.0f0)
    set_command!(prepared.dm, prepared.adopted_command)

    _form_plant_frame!(prepared)
    permutedims!(prepared.fga_frame, observation_storage(prepared.observation), (2, 1))
    JuliaFilterGraph.process!(
        outputs,
        graph,
        inputs,
        (image=SampleMetadata(1; terminal=true),),
    )
    reference_slopes = copy(outputs.slopes)
    replace_parameters!(graph, Symbol("reference-slopes") => reference_slopes)

    probe_plan = AdaptiveOpticsCalibration.prepare(
        AdaptiveOpticsCalibration.ProbeBases.ZonalPushPull(
            Float32[_CALIBRATION_POKE_M],
        ),
        nothing,
    )
    probe_product = AdaptiveOpticsCalibration.allocate_result(probe_plan)
    AdaptiveOpticsCalibration.process!(
        probe_product,
        AdaptiveOpticsCalibration.allocate_workspace(probe_plan),
        probe_plan,
        nothing,
    )
    probe_commands = AdaptiveOpticsCalibration.ProbeBases.probe_commands(
        probe_product,
    )
    responses = zeros(Float32, size(probe_commands, 1), 8)
    for exposure in axes(probe_commands, 1)
        copyto!(prepared.adopted_command, @view probe_commands[exposure, :])
        set_command!(prepared.dm, prepared.adopted_command)
        _form_plant_frame!(prepared)
        permutedims!(
            prepared.fga_frame,
            observation_storage(prepared.observation),
            (2, 1),
        )
        JuliaFilterGraph.process!(
            outputs,
            graph,
            inputs,
            (image=SampleMetadata(exposure + 1; terminal=true),),
        )
        _interleave_slopes!(@view(responses[exposure, :]), outputs.slopes)
    end

    interaction_plan = AdaptiveOpticsCalibration.prepare(
        AdaptiveOpticsCalibration.InteractionMatrices.ZonalPushPull(
            Float32[_CALIBRATION_POKE_M],
        ),
        AdaptiveOpticsCalibration.InteractionMatrices.ResponseSpecification(8),
    )
    interaction_product = AdaptiveOpticsCalibration.allocate_result(interaction_plan)
    AdaptiveOpticsCalibration.process!(
        interaction_product,
        AdaptiveOpticsCalibration.allocate_workspace(interaction_plan),
        interaction_plan,
        responses,
    )
    interaction = AdaptiveOpticsCalibration.InteractionMatrices.interaction_matrix(
        interaction_product,
    )
    specification = ReconstructorSpecification(8, 1, Float32)
    method = StrokeWeightedTikhonov(1; tikhonov_scale=0.0f0)
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    product = AdaptiveOpticsCalibration.allocate_result(plan)
    calibration_inputs = StrokeWeightedReconstructorInputs(
        interaction,
        Matrix{Float32}(I, 8, 8),
        Matrix{Float32}(I, 1, 1),
    )
    AdaptiveOpticsCalibration.process!(
        product,
        nothing,
        plan,
        calibration_inputs,
    )
    replace_parameters!(graph, :reconstructor => reconstructor(product))
    identity = _calibration_identity(
        prepared.plant_signature,
        _estimator_signature(
            graph_path,
            reference_slopes,
            reconstructor(product),
            controller_to_vdm,
            active_to_full_vdm,
            vdm_to_pdm,
            subaperture_order,
        ),
        subaperture_order=subaperture_order,
    )
    return S1CalibrationProduct(identity, interaction, product)
end

"""
    prepare_s1_lockstep(; disturbance_command=3f-8)

Prepare the deterministic S1 Float32 composition. AOS owns the physical
SHWFS, one-actuator plant, and probe exposures. AdaptiveOpticsCalibration
constructs the zonal probe basis, estimates the interaction matrix from the
complete response cycle, and prepares the cold reconstructor. FGA/JFG own the
per-frame estimator and RTC chain.
"""
function prepare_s1_lockstep(; disturbance_command::Float32=3.0f-8)
    plant = _prepare_plant()
    graph_path = joinpath(@__DIR__, "s1-shwfs-f32.conf")
    graph = prepare_graph(
        graph_path;
        algorithms=algorithms(),
    )
    graph.input_formats.image.schema == DETECTOR_FRAME_SCHEMA || error(
        "S1 detector-frame schema does not match the prepared FGA Graph",
    )
    graph.output_formats.demanded.schema == DEMANDED_PDM_COMMAND_V1 || error(
        "S1 demanded-PDM-command schema does not match FGA",
    )
    controller_to_vdm = ones(Float32, 1, 1)
    active_to_full_vdm = ones(Float32, 1, 1)
    vdm_to_pdm = ones(Float32, 1, 1)
    replace_parameters!(
        graph,
        Symbol("controller-to-vdm") => controller_to_vdm,
        Symbol("active-to-full") => active_to_full_vdm,
        Symbol("vdm-to-pdm") => vdm_to_pdm,
    )
    subaperture_order = _prepared_subaperture_order(graph)
    fga_frame = zeros(Float32, 8, 8)
    outputs = (
        demanded=zeros(Float32, 1),
        slopes=zeros(Float32, 4, 2),
        flux=zeros(Float32, 4),
        validity=fill(false, 4),
    )
    inputs = (image=fga_frame,)
    adopted_command = zeros(Float32, 1)
    disturbance_opd = zeros(Float32, 8, 8)
    applied_command = similar(adopted_command)
    calibration_owner = (;
        pupil=plant.pupil,
        dm=plant.dm,
        rate=plant.rate,
        observation=plant.observation,
        optics_plan=plant.optics_plan,
        acquisition_plan=plant.acquisition_plan,
        rng=plant.rng,
        plant_signature=plant.plant_signature,
        disturbance_opd,
        fga_frame,
        adopted_command,
        applied_command,
    )
    calibration = _calibrate_reconstructor!(
        calibration_owner,
        graph,
        outputs,
        inputs,
        graph_path,
        controller_to_vdm,
        active_to_full_vdm,
        vdm_to_pdm,
        subaperture_order,
    )
    JuliaFilterGraph.reset!(graph)

    prepared = PreparedS1Lockstep(
        plant.pupil,
        plant.dm,
        plant.rate,
        plant.observation,
        plant.optics_plan,
        plant.acquisition_plan,
        plant.rng,
        graph,
        outputs,
        inputs,
        calibration,
        disturbance_opd,
        fga_frame,
        adopted_command,
        applied_command,
        S1LockstepState(0, false, false),
    )

    adopted_command[1] = disturbance_command
    set_command!(prepared.dm, adopted_command)
    update_surface!(prepared.dm)
    copyto!(prepared.disturbance_opd, surface_opd(prepared.dm))
    fill!(adopted_command, 0.0f0)
    set_command!(prepared.dm, adopted_command)
    fill!(prepared.applied_command, 0.0f0)
    Random.seed!(prepared.rng, _RNG_SEED)
    return prepared
end

"""Form the next complete AOS detector frame without yet invoking the RTC."""
function produce_detector_frame!(prepared::PreparedS1Lockstep)
    state = prepared.state
    state.blocked && throw(ArgumentError("S1 lockstep is blocked; reset is required"))
    state.outstanding && throw(ArgumentError(
        "the outstanding detector frame must be accepted or reset before advancing the plant",
    ))
    values = _form_plant_frame!(prepared)
    sequence = state.sequence + UInt64(1)
    timestamp = Int64(sequence) * _FRAME_PERIOD_NANOSECONDS
    state.sequence = sequence
    state.outstanding = true
    return S1DetectorFrame(
        values,
        DETECTOR_FRAME_SCHEMA,
        prepared.calibration.identity.signature,
        sequence,
        timestamp,
        true,
        false,
        false,
    )
end

@inline function _reject!(prepared, disposition::S1FrameDisposition)
    prepared.state.blocked = true
    return disposition
end

"""Validate, process, and adopt one outstanding detector frame."""
function process_detector_frame!(
    prepared::PreparedS1Lockstep,
    frame::S1DetectorFrame,
)
    state = prepared.state
    state.blocked && return FrameRejectedBlocked
    state.outstanding || return _reject!(prepared, FrameRejectedNoOutstandingFrame)
    frame.sequence == state.sequence || return _reject!(prepared, FrameRejectedSequence)
    expected_timestamp = Int64(frame.sequence) * _FRAME_PERIOD_NANOSECONDS
    frame.timestamp_nanoseconds == expected_timestamp ||
        return _reject!(prepared, FrameRejectedTimestamp)
    frame.terminal || return _reject!(prepared, FrameRejectedIncomplete)
    !frame.discontinuity || return _reject!(prepared, FrameRejectedDiscontinuous)
    !frame.corrupted || return _reject!(prepared, FrameRejectedCorrupted)
    frame.schema == DETECTOR_FRAME_SCHEMA ||
        return _reject!(prepared, FrameRejectedSchema)
    frame.calibration_signature == prepared.calibration.identity.signature ||
        return _reject!(prepared, FrameRejectedCalibrationSignature)
    eltype(frame.values) === Float32 ||
        return _reject!(prepared, FrameRejectedNumericType)
    size(frame.values) == (8, 8) || return _reject!(prepared, FrameRejectedShape)
    _all_finite(frame.values) || return _reject!(prepared, FrameRejectedNonFinite)

    permutedims!(prepared.fga_frame, frame.values, (2, 1))
    metadata = SampleMetadata(
        frame.sequence;
        timestamp_nanoseconds=frame.timestamp_nanoseconds,
        terminal=true,
    )
    result = JuliaFilterGraph.process!(
        prepared.outputs,
        prepared.graph,
        prepared.inputs,
        (image=metadata,),
    )
    result isa NamedTuple || return _reject!(prepared, FrameRejectedRTCFailure)
    result.demanded && result.slopes && result.flux && result.validity ||
        return _reject!(prepared, FrameRejectedRTCUnavailable)

    published = JuliaFilterGraph.output_metadata(prepared.graph).demanded
    published isa SampleMetadata ||
        return _reject!(prepared, FrameRejectedRTCMetadata)
    published.sequence == frame.sequence && published.has_timestamp &&
        published.timestamp_nanoseconds == frame.timestamp_nanoseconds &&
        published.terminal && !published.discontinuity && !published.corrupted ||
        return _reject!(prepared, FrameRejectedRTCMetadata)
    demanded = prepared.outputs.demanded
    length(demanded) == 1 && _all_finite(demanded) ||
        return _reject!(prepared, FrameRejectedCommand)

    copyto!(prepared.adopted_command, demanded)
    set_command!(prepared.dm, prepared.adopted_command)
    state.outstanding = false
    return FrameAccepted
end

"""Run one accepted frame-to-command lockstep cycle and return its sequence."""
function step_lockstep!(prepared::PreparedS1Lockstep)
    frame = produce_detector_frame!(prepared)
    disposition = process_detector_frame!(prepared, frame)
    disposition === FrameAccepted || error(
        "the internally formed S1 frame was rejected with $disposition",
    )
    return frame.sequence
end

"""Restore the RTC and exchange state while retaining the calibrated products."""
function reset_s1_lockstep!(prepared::PreparedS1Lockstep)
    JuliaFilterGraph.reset!(prepared.graph)
    fill!(prepared.adopted_command, 0.0f0)
    fill!(prepared.applied_command, 0.0f0)
    fill!(prepared.fga_frame, 0.0f0)
    fill!(prepared.outputs.demanded, 0.0f0)
    fill!(prepared.outputs.slopes, 0.0f0)
    fill!(prepared.outputs.flux, 0.0f0)
    fill!(prepared.outputs.validity, false)
    set_command!(prepared.dm, prepared.adopted_command)
    reset_opd!(prepared.pupil)
    Random.seed!(prepared.rng, _RNG_SEED)
    prepared.state.sequence = 0
    prepared.state.outstanding = false
    prepared.state.blocked = false
    return prepared
end

end # module AOSFGALockstep
