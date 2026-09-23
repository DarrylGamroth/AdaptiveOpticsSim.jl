module AOSFGAPyramid

using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using FilterGraphAlgorithms
using JuliaFilterGraph
using Random

import FilterGraphAlgorithms: commit!, process!, reset!

export S4_DETECTOR_FRAME_SCHEMA, S4DetectorFrame, S4FrameDisposition
export S4FrameAccepted, S4FrameRejectedBlocked, S4FrameRejectedNoOutstandingFrame
export S4FrameRejectedSequence, S4FrameRejectedTimestamp, S4FrameRejectedIncomplete
export S4FrameRejectedDiscontinuous, S4FrameRejectedCorrupted, S4FrameRejectedSchema
export S4FrameRejectedIdentity, S4FrameRejectedCalibrationSignature
export S4FrameRejectedOrder, S4FrameRejectedNumericType, S4FrameRejectedShape
export S4FrameRejectedNonFinite, S4FrameRejectedFGA, S4FrameRejectedCommand
export prepare_s4_pyramid, produce_s4_detector_frame!, process_s4_detector_frame!
export step_s4_pyramid!, reset_s4_pyramid!, direct_s4_signal!

const S4_DETECTOR_FRAME_SCHEMA =
    "org.adaptiveopticssim.integration.pyramid-electron-counts-xy.f32/1"
const _FRAME_PERIOD_NANOSECONDS = Int64(1_000_000)
const _MODEL_TIMESTAMP_ORIGIN_NANOSECONDS = Int64(0)
const _FGA_VERSION_CLAIM = v"0.5.2"
const _JFG_VERSION_CLAIM = v"0.2.3"
const _FGA_TREE_CLAIM =
    "922a65f3ac9486abeb7baa3e29f08794133b5bc6"
const _SIGNATURE_OFFSET = UInt64(0xcbf29ce484222325)
const _SIGNATURE_PRIME = UInt64(0x00000100000001b3)

@enum S4FrameDisposition::UInt8 begin
    S4FrameAccepted = 0
    S4FrameRejectedBlocked = 1
    S4FrameRejectedNoOutstandingFrame = 2
    S4FrameRejectedSequence = 3
    S4FrameRejectedTimestamp = 4
    S4FrameRejectedIncomplete = 5
    S4FrameRejectedDiscontinuous = 6
    S4FrameRejectedCorrupted = 7
    S4FrameRejectedSchema = 8
    S4FrameRejectedIdentity = 9
    S4FrameRejectedCalibrationSignature = 10
    S4FrameRejectedOrder = 11
    S4FrameRejectedNumericType = 12
    S4FrameRejectedShape = 13
    S4FrameRejectedNonFinite = 14
    S4FrameRejectedFGA = 15
    S4FrameRejectedCommand = 16
end

"""One complete acquired physical Pyramid detector frame at the AOS/FGA boundary."""
struct S4DetectorFrame{A<:AbstractMatrix}
    values::A
    schema::String
    identity_signature::UInt64
    calibration_signature::UInt64
    order_signature::UInt64
    sequence::UInt64
    timestamp_nanoseconds::Int64
    terminal::Bool
    discontinuity::Bool
    corrupted::Bool
end

"""Immutable interpretation and release identity for the S4 plant/RTC boundary."""
struct S4CalibrationIdentity{T<:AbstractFloat}
    detector_axes::NTuple{2,Symbol}
    estimator_frame_axes::NTuple{2,Symbol}
    pupil_origins::NTuple{4,NTuple{2,Int}}
    pupil_order::NTuple{4,Symbol}
    support_order::NTuple{4,NTuple{2,Int}}
    component_order::NTuple{2,Symbol}
    normalization::Symbol
    numeric_type::Type{T}
    detector_units::Symbol
    signal_units::Symbol
    command_units::Symbol
    model_timestamp_origin_nanoseconds::Int64
    model_timestamp_period_nanoseconds::Int64
    fga_version::VersionNumber
    fga_tree_claim::String
    order_signature::UInt64
    calibration_signature::UInt64
    signature::UInt64
end

mutable struct S4State
    sequence::UInt64
    outstanding::Bool
    blocked::Bool
end

struct PreparedS4Pyramid{
    Pupil,
    Rate,
    Observation,
    OpticsPlan,
    AcquisitionPlan,
    RNG,
    Image,
    Reconstructor,
    Integrator,
}
    pupil::Pupil
    rate::Rate
    observation::Observation
    optics_plan::OpticsPlan
    acquisition_plan::AcquisitionPlan
    rng::RNG
    image::Image
    reconstructor::Reconstructor
    integrator::Integrator
    identity::S4CalibrationIdentity{Float32}
    fga_frame::Matrix{Float32}
    signal::Matrix{Float32}
    divisor::Vector{Float32}
    residual::Vector{Float32}
    candidate_command::Vector{Float32}
    adopted_command::Vector{Float32}
    disturbance_opd::Matrix{Float32}
    state::S4State
end

@inline _signature_byte(signature::UInt64, byte::UInt8) =
    xor(signature, UInt64(byte)) * _SIGNATURE_PRIME

function _signature_bytes(signature::UInt64, bytes)
    @inbounds for byte in bytes
        signature = _signature_byte(signature, UInt8(byte))
    end
    return signature
end

function _signature_text(signature::UInt64, text::AbstractString)
    signature = _signature_bytes(signature, codeunits(text))
    return _signature_byte(signature, 0x00)
end

function _signature_integer(signature::UInt64, value::Integer)
    encoded = reinterpret(UInt64, Int64(value))
    @inbounds for shift in 0:8:56
        signature = _signature_byte(signature, UInt8((encoded >> shift) & 0xff))
    end
    return signature
end

function _signature_integer(signature::UInt64, value::Unsigned)
    encoded = UInt64(value)
    @inbounds for shift in 0:8:56
        signature = _signature_byte(signature, UInt8((encoded >> shift) & 0xff))
    end
    return signature
end

function _signature_float32(signature::UInt64, value::Float32)
    encoded = reinterpret(UInt32, value)
    @inbounds for shift in 0:8:24
        signature = _signature_byte(signature, UInt8((encoded >> shift) & 0xff))
    end
    return signature
end

function _signature_value(signature::UInt64, value::Symbol)
    return _signature_text(signature, String(value))
end

function _signature_value(signature::UInt64, values::Tuple)
    signature = _signature_integer(signature, length(values))
    for value in values
        signature = _signature_value(signature, value)
    end
    return signature
end

_signature_value(signature::UInt64, value::Integer) =
    _signature_integer(signature, value)

function _signature_matrix(signature::UInt64, values::AbstractMatrix{Float32})
    signature = _signature_integer(signature, size(values, 1))
    signature = _signature_integer(signature, size(values, 2))
    @inbounds for value in values
        signature = _signature_float32(signature, value)
    end
    return signature
end

@inline function _all_finite(values)
    @inbounds for value in values
        isfinite(value) || return false
    end
    return true
end

const _PUPIL_ORIGINS = ((0, 0), (2, 0), (2, 2), (0, 2))
const _PUPIL_ORDER = (:q1_top_left, :q2_bottom_left, :q3_bottom_right, :q4_top_right)
const _SUPPORT_ORDER = ((1, 1), (2, 1), (1, 2), (2, 2))

function _order_signature()
    signature = _signature_text(_SIGNATURE_OFFSET, "AOS-FGA-S4-PYRAMID-ORDER/1")
    signature = _signature_value(signature, (:x, :y))
    signature = _signature_value(signature, (:row, :column))
    signature = _signature_value(signature, _PUPIL_ORIGINS)
    signature = _signature_value(signature, _PUPIL_ORDER)
    signature = _signature_value(signature, _SUPPORT_ORDER)
    return _signature_value(signature, (:x, :y))
end

function _calibration_signature(reference::AbstractMatrix{Float32},
    gain::AbstractMatrix{Float32})
    signature = _signature_text(_SIGNATURE_OFFSET, "AOS-FGA-S4-PYRAMID-CALIBRATION/1")
    signature = _signature_matrix(signature, reference)
    signature = _signature_matrix(signature, gain)
    signature = _signature_text(signature, "mean-valid-flux")
    return signature
end

function _identity(reference::AbstractMatrix{Float32},
    gain::AbstractMatrix{Float32})
    Base.pkgversion(FilterGraphAlgorithms) == _FGA_VERSION_CLAIM || error(
        "S4 requires FilterGraphAlgorithms $_FGA_VERSION_CLAIM",
    )
    Base.pkgversion(JuliaFilterGraph) == _JFG_VERSION_CLAIM || error(
        "S4 requires JuliaFilterGraph $_JFG_VERSION_CLAIM",
    )
    order_signature = _order_signature()
    calibration_signature = _calibration_signature(reference, gain)
    signature = _signature_text(_SIGNATURE_OFFSET, "AOS-FGA-S4-PYRAMID-IDENTITY/1")
    signature = _signature_integer(signature, order_signature)
    signature = _signature_integer(signature, calibration_signature)
    signature = _signature_integer(signature, _MODEL_TIMESTAMP_ORIGIN_NANOSECONDS)
    signature = _signature_integer(signature, _FRAME_PERIOD_NANOSECONDS)
    signature = _signature_text(signature, string(_FGA_VERSION_CLAIM))
    signature = _signature_text(signature, _FGA_TREE_CLAIM)
    return S4CalibrationIdentity{Float32}(
        (:x, :y),
        (:row, :column),
        _PUPIL_ORIGINS,
        _PUPIL_ORDER,
        _SUPPORT_ORDER,
        (:x, :y),
        :mean_valid_flux,
        Float32,
        :electron_count,
        :normalized_pyramid_i4q,
        :metre,
        _MODEL_TIMESTAMP_ORIGIN_NANOSECONDS,
        _FRAME_PERIOD_NANOSECONDS,
        _FGA_VERSION_CLAIM,
        _FGA_TREE_CLAIM,
        order_signature,
        calibration_signature,
        signature,
    )
end

function _prepare_plant()
    T = Float32
    telescope = Telescope(
        resolution=4,
        diameter=one(T),
        central_obstruction=zero(T),
        fov_arcsec=zero(T),
        pupil_reflectivity=one(T),
        T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=:custom,
        magnitude=zero(T),
        separation_arcsec=zero(T), position_angle_deg=zero(T),
        wavelength=750.0f-9,
        photon_irradiance=2.0f8,
        radiometry=PhysicalPhotonIrradianceSource(),
        T=T,
    )
    # This constructs only the diffractive four-pupil optical front end.  The
    # FGA calibration owns differential estimation; no AOS estimator is used.
    sensor = PyramidWFS(
        telescope;
        pupil_samples=2,
        modulation=zero(T),
        diffraction_padding=2,
        T=T,
    )
    front_end = PyramidOpticalFrontEnd(sensor, source)
    rate = pyramid_rate_map(front_end, pupil)
    optics_plan = prepare_wfs_optics(front_end, pupil, rate)
    detector = Detector(
        noise=NoiseNone(),
        exposure_duration=1.0f-3,
        qe=one(T),
        response_model=NullFrameResponse(),
        T=T,
    )
    observation = WFSObservation(
        similar(intensity_values(rate));
        units=:electron_count,
        layout=:four_pupil_mosaic,
    )
    acquisition_plan = prepare_wfs_acquisition(
        detector,
        rate,
        observation;
        source,
    )
    return (; pupil, rate, observation, optics_plan, acquisition_plan,
        rng=Xoshiro(0x5334))
end

@inline function _form_s4_detector_frame!(prepared)
    copyto!(opd_map(prepared.pupil), prepared.disturbance_opd)
    form_wfs_optical_products!(prepared.rate, prepared.pupil, prepared.optics_plan)
    acquire_wfs_observation!(
        prepared.observation,
        prepared.rate,
        prepared.acquisition_plan,
        prepared.rng,
    )
    return observation_storage(prepared.observation)
end

@inline function _transpose_detector_frame!(
    destination::AbstractMatrix{Float32}, source)
    # AOS detector axes are (x, y); the FGA calibrated-image contract is
    # (row=y, column=x).  This writes the persistent exchange buffer.
    permutedims!(destination, source, (2, 1))
    return destination
end

"""Independent I4Q equation used as the S4 boundary oracle."""
function direct_s4_signal!(
    signal::AbstractMatrix{Float32},
    divisor::AbstractVector{Float32},
    frame::AbstractMatrix{Float32},
    reference::AbstractMatrix{Float32},
    gain::AbstractMatrix{Float32},
)
    total = 0.0f0
    @inbounds for (row, column) in _SUPPORT_ORDER
        total += frame[row, column] + frame[row + 2, column] +
            frame[row + 2, column + 2] + frame[row, column + 2]
    end
    divisor[1] = total / Float32(length(_SUPPORT_ORDER))
    index = 1
    @inbounds for (row, column) in _SUPPORT_ORDER
        q1 = frame[row, column]
        q2 = frame[row + 2, column]
        q3 = frame[row + 2, column + 2]
        q4 = frame[row, column + 2]
        signal[index, 1] = gain[index, 1] *
            ((q1 - q2 + q4 - q3) / divisor[1] - reference[index, 1])
        signal[index, 2] = gain[index, 2] *
            ((q1 - q4 + q2 - q3) / divisor[1] - reference[index, 2])
        index += 1
    end
    return signal
end

function _flat_reference!(prepared, reference::AbstractMatrix{Float32})
    fill!(prepared.disturbance_opd, 0.0f0)
    flat = _form_s4_detector_frame!(prepared)
    _transpose_detector_frame!(prepared.fga_frame, flat)
    raw = zeros(Float32, 4, 2)
    divisor = zeros(Float32, 1)
    direct_s4_signal!(raw, divisor, prepared.fga_frame, raw, ones(Float32, 4, 2))
    copyto!(reference, raw)
    return reference
end

function _prepare_reconstructor()
    # One explicit command coordinate reconstructs the mean X I4Q component.
    matrix = Float32[0.25, 0.25, 0.25, 0.25, 0, 0, 0, 0]
    return prepare_algorithm(
        PyramidReconstructorF32,
        PyramidReconstructorF32Configuration(
            reconstructed_count=1,
            selected_coordinate_count=4,
            initial_reconstructor=matrix,
            reconstructed_schema=CONTROLLER_RESIDUAL_ERROR_V1,
        ),
    )
end

"""
    prepare_s4_pyramid()

Prepare the independent S4 composition: physical diffractive Pyramid photon
rate formation, explicit detector acquisition, a fixed AOS-to-FGA transpose,
and FGA calibrated I4Q/reconstruction/integration.  The flat detector frame
defines the reference signal; all calibration arrays are explicit.
"""
function prepare_s4_pyramid()
    plant = _prepare_plant()
    owner = (; plant..., fga_frame=zeros(Float32, 4, 4),
        disturbance_opd=zeros(Float32, 4, 4))
    reference = zeros(Float32, 4, 2)
    _flat_reference!(owner, reference)
    gain = ones(Float32, 4, 2)
    image = prepare_algorithm(
        PyramidImageF32,
        PyramidImageF32Configuration(
            image_rows=4,
            image_columns=4,
            pupil_rows=2,
            pupil_columns=2,
            pupil_origins=[[0, 0], [2, 0], [2, 2], [0, 2]],
            i4q_support=Bool[true, true, true, true],
            initial_reference_signal=vec(reference),
            initial_optical_gain=vec(gain),
            normalization_policy="mean-valid-flux",
            incidence_flux=0.0f0,
        ),
    )
    reconstructor = _prepare_reconstructor()
    integrator = (
        plan=LeakyIntegratorPlan(-1.0f0, 0.0f0),
        workspace=LeakyIntegratorWorkspace(Float32, 1),
    )
    identity = _identity(reference, gain)
    prepared = PreparedS4Pyramid(
        plant.pupil,
        plant.rate,
        plant.observation,
        plant.optics_plan,
        plant.acquisition_plan,
        plant.rng,
        image,
        reconstructor,
        integrator,
        identity,
        owner.fga_frame,
        zeros(Float32, 4, 2),
        zeros(Float32, 1),
        zeros(Float32, 1),
        zeros(Float32, 1),
        zeros(Float32, 1),
        owner.disturbance_opd,
        S4State(0, false, false),
    )
    # A deterministic non-flat pupil OPD makes the qualified command path
    # observable while retaining the flat frame as calibration reference.
    @inbounds for column in axes(prepared.disturbance_opd, 2),
        row in axes(prepared.disturbance_opd, 1)
        prepared.disturbance_opd[row, column] =
            30.0f-9 * Float32(row - 1) / 3.0f0
    end
    Random.seed!(prepared.rng, 0x5334)
    return prepared
end

"""Form one complete physical detector frame; the FGA command is not run here."""
function produce_s4_detector_frame!(prepared::PreparedS4Pyramid)
    state = prepared.state
    state.blocked && throw(ArgumentError("S4 Pyramid exchange is blocked; reset is required"))
    state.outstanding && throw(ArgumentError(
        "the outstanding S4 detector frame must be accepted or reset before advancing the plant",
    ))
    values = _form_s4_detector_frame!(prepared)
    sequence = state.sequence + UInt64(1)
    timestamp = _MODEL_TIMESTAMP_ORIGIN_NANOSECONDS +
        Int64(sequence) * _FRAME_PERIOD_NANOSECONDS
    state.sequence = sequence
    state.outstanding = true
    identity = prepared.identity
    return S4DetectorFrame(
        values,
        S4_DETECTOR_FRAME_SCHEMA,
        identity.signature,
        identity.calibration_signature,
        identity.order_signature,
        sequence,
        timestamp,
        true,
        false,
        false,
    )
end

@inline function _reject!(prepared, disposition::S4FrameDisposition)
    prepared.state.blocked = true
    return disposition
end

"""Validate one boundary frame, execute FGA, and atomically adopt its command."""
function process_s4_detector_frame!(
    prepared::PreparedS4Pyramid,
    frame::S4DetectorFrame,
)
    state = prepared.state
    state.blocked && return S4FrameRejectedBlocked
    state.outstanding || return _reject!(prepared, S4FrameRejectedNoOutstandingFrame)
    frame.sequence == state.sequence || return _reject!(prepared, S4FrameRejectedSequence)
    expected_timestamp = _MODEL_TIMESTAMP_ORIGIN_NANOSECONDS +
        Int64(frame.sequence) * _FRAME_PERIOD_NANOSECONDS
    frame.timestamp_nanoseconds == expected_timestamp ||
        return _reject!(prepared, S4FrameRejectedTimestamp)
    frame.terminal || return _reject!(prepared, S4FrameRejectedIncomplete)
    !frame.discontinuity || return _reject!(prepared, S4FrameRejectedDiscontinuous)
    !frame.corrupted || return _reject!(prepared, S4FrameRejectedCorrupted)
    frame.schema == S4_DETECTOR_FRAME_SCHEMA || return _reject!(prepared, S4FrameRejectedSchema)
    identity = prepared.identity
    frame.identity_signature == identity.signature ||
        return _reject!(prepared, S4FrameRejectedIdentity)
    frame.calibration_signature == identity.calibration_signature ||
        return _reject!(prepared, S4FrameRejectedCalibrationSignature)
    frame.order_signature == identity.order_signature ||
        return _reject!(prepared, S4FrameRejectedOrder)
    eltype(frame.values) === Float32 || return _reject!(prepared, S4FrameRejectedNumericType)
    size(frame.values) == (4, 4) || return _reject!(prepared, S4FrameRejectedShape)
    _all_finite(frame.values) || return _reject!(prepared, S4FrameRejectedNonFinite)

    _transpose_detector_frame!(prepared.fga_frame, frame.values)
    process!(prepared.signal, prepared.divisor, prepared.image, prepared.fga_frame) === nothing ||
        return _reject!(prepared, S4FrameRejectedFGA)
    process!(prepared.residual, prepared.reconstructor, prepared.signal) === nothing ||
        return _reject!(prepared, S4FrameRejectedFGA)
    result = process!(
        prepared.candidate_command,
        prepared.integrator.workspace,
        prepared.integrator.plan,
        prepared.residual,
    )
    result isa AbstractVector || return _reject!(prepared, S4FrameRejectedFGA)
    _all_finite(prepared.candidate_command) || return _reject!(prepared, S4FrameRejectedCommand)

    # The command and committed controller state change only after every FGA
    # stage succeeds.  Rejections leave the prior adopted command intact.
    copyto!(prepared.adopted_command, prepared.candidate_command)
    commit!(prepared.integrator.workspace)
    state.outstanding = false
    return S4FrameAccepted
end

"""Run one accepted S4 plant-to-command frame and return its sequence."""
function step_s4_pyramid!(prepared::PreparedS4Pyramid)
    frame = produce_s4_detector_frame!(prepared)
    disposition = process_s4_detector_frame!(prepared, frame)
    disposition === S4FrameAccepted || error("the internally formed S4 frame was rejected with $disposition")
    return frame.sequence
end

"""Reset the exchange, FGA control state, and deterministic plant input."""
function reset_s4_pyramid!(prepared::PreparedS4Pyramid)
    reset!(prepared.image.workspace)
    reset!(prepared.reconstructor.workspace)
    reset!(prepared.integrator.workspace)
    fill!(prepared.fga_frame, 0.0f0)
    fill!(prepared.signal, 0.0f0)
    fill!(prepared.divisor, 0.0f0)
    fill!(prepared.residual, 0.0f0)
    fill!(prepared.candidate_command, 0.0f0)
    fill!(prepared.adopted_command, 0.0f0)
    Random.seed!(prepared.rng, 0x5334)
    prepared.state.sequence = 0
    prepared.state.outstanding = false
    prepared.state.blocked = false
    return prepared
end

end # module AOSFGAPyramid
