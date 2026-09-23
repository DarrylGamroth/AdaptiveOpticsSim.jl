module AOSFGACurvature

using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using FilterGraphAlgorithms
using JuliaFilterGraph
using Random

import FilterGraphAlgorithms: process!

export prepare_s4_curvature
export produce_s4_curvature_frame!, produce_s4_curvature_channels!
export process_s4_curvature_image!, process_s4_curvature_channels!
export step_s4_curvature_frame!, step_s4_curvature_channels!
export reset_s4_curvature_exchange!
export transpose_aos_curvature_frame!, frozen_aos_curvature_frame
export frozen_aos_curvature_channels
export S4CurvatureObservation, S4CurvatureDisposition
export S4CurvatureAccepted, S4CurvatureRejectedBlocked
export S4CurvatureRejectedNoOutstandingObservation
export S4CurvatureRejectedSequence, S4CurvatureRejectedModelTimestamp
export S4CurvatureRejectedExposure, S4CurvatureRejectedLayout
export S4CurvatureRejectedCalibrationSignature, S4CurvatureRejectedEstimator

const _FGA_VERSION_CLAIM = v"0.5.2"
const _JFG_VERSION_CLAIM = v"0.2.3"
const _FGA_TREE_CLAIM = "922a65f3ac9486abeb7baa3e29f08794133b5bc6"
const _JFG_TREE_CLAIM = "1893cea2d29f0cd89303939308e23ce7546bc240"
const _CALIBRATION_SIGNATURE = UInt64(0x43555256)
const _NANOSECONDS_PER_SECOND = Int64(1_000_000_000)
const _FRAME_LAYOUT = :curvature_branch_regions
const _CHANNEL_LAYOUT = :curvature_branch_channels

@enum S4CurvatureDisposition::UInt8 begin
    S4CurvatureAccepted = 0
    S4CurvatureRejectedBlocked = 1
    S4CurvatureRejectedNoOutstandingObservation = 2
    S4CurvatureRejectedSequence = 3
    S4CurvatureRejectedModelTimestamp = 4
    S4CurvatureRejectedExposure = 5
    S4CurvatureRejectedLayout = 6
    S4CurvatureRejectedCalibrationSignature = 7
    S4CurvatureRejectedEstimator = 8
end

"""One complete acquired Curvature observation at the AOS/FGA boundary."""
struct S4CurvatureObservation{A<:AbstractMatrix,T<:AbstractFloat}
    values::A
    layout::Symbol
    calibration_signature::UInt64
    sequence::UInt64
    exposure_duration::T
    model_timestamp_nanoseconds::Int64
end

"""Single-writer association state for one Curvature observation layout."""
mutable struct S4CurvatureExchangeState{T<:AbstractFloat}
    sequence::UInt64
    outstanding::Bool
    blocked::Bool
    published_sequence::UInt64
    published_exposure_duration::T
    published_model_timestamp_nanoseconds::Int64
end

S4CurvatureExchangeState(::Type{T}) where {T<:AbstractFloat} =
    S4CurvatureExchangeState(UInt64(0), false, false, UInt64(0), zero(T), Int64(0))

const _FGA_SUPPORT = Bool[
    true false
    true true
]

const _FGA_REFERENCE = Float32[
    0.1 0.0
    -0.1 0.2
]

"""Preallocated AOS `(x, y)` to FGA `(row=y, column=x)` transfer."""
@inline function transpose_aos_curvature_frame!(
    destination::AbstractMatrix{T}, source::AbstractMatrix{T},
) where {T}
    size(destination) == reverse(size(source)) || throw(DimensionMismatch(
        "AOS Curvature source and FGA destination axes are incompatible",
    ))
    permutedims!(destination, source, (2, 1))
    return destination
end

"""Asymmetric packed AOS frame: positive branch followed by negative branch."""
frozen_aos_curvature_frame() = Float32[
    10 30
    20 40
    2 9
    7 5
]

"""The same paired branches in AOS branch-by-channel order."""
frozen_aos_curvature_channels() = Float32[
    10 30 20 40
    2 9 7 5
]

function _image_owner(; support=_FGA_SUPPORT, reference=_FGA_REFERENCE,
    branch_scales=(1.25f0, 0.75f0))
    plan = CurvaturePairedImagePlan(
        (2, 4), (2, 2), (2, 2), ((0, 0), (0, 2)), support,
        reference, branch_scales, _CALIBRATION_SIGNATURE,
    )
    return (; plan, workspace=CurvaturePairedSignalWorkspace(plan),
        signal=zeros(Float32, 4), signature=UInt64[_CALIBRATION_SIGNATURE])
end

function _channel_owner(; support=_FGA_SUPPORT, reference=_FGA_REFERENCE,
    branch_scales=(1.25f0, 0.75f0))
    plan = CurvaturePairedChannelPlan(
        (2, 4), (2, 2), (0, 1), collect(0:3), support,
        reference, branch_scales, _CALIBRATION_SIGNATURE,
    )
    return (; plan, workspace=CurvaturePairedSignalWorkspace(plan),
        signal=zeros(Float32, 4), signature=UInt64[_CALIBRATION_SIGNATURE])
end

function _prepare_plant()
    T = Float32
    telescope = Telescope(
        resolution=8, diameter=T(4), central_obstruction=zero(T),
        fov_arcsec=zero(T), pupil_reflectivity=one(T), T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=:custom, magnitude=zero(T), separation_arcsec=zero(T), position_angle_deg=zero(T),
        wavelength=750.0f-9, photon_irradiance=12.0f0,
        radiometry=PhysicalPhotonIrradianceSource(), T=T,
    )
    sensor = CurvatureWFS(
        telescope; pupil_samples=2, diffraction_padding=2, T=T,
    )
    front_end = CurvatureOpticalFrontEnd(sensor, source)
    rates = curvature_rate_maps(front_end, pupil)
    optics_plan = prepare_wfs_optics(front_end, pupil, rates)

    frame_detector = Detector(
        noise=NoiseNone(), exposure_duration=one(T), qe=one(T),
        response_model=NullFrameResponse(), T=T,
    )
    frame_observation = WFSObservation(
        zeros(T, 4, 2); units=:electron_count,
        layout=:curvature_branch_regions,
    )
    frame_acquisition_plan = prepare_wfs_acquisition(
        CurvaturePackedAcquisition(frame_detector), rates,
        frame_observation,
    )

    channel_detector = LinearAPDDetector(
        topology=LinearAPDChannelBank(8), exposure_duration=one(T),
        qe=one(T), avalanche_gain=one(T), conversion_gain=one(T),
        dark_current=zero(T), noise=NoiseNone(), T=T,
    )
    channel_observation = WFSObservation(
        zeros(T, 2, 4); units=:detector_count,
        layout=:curvature_branch_channels,
    )
    channel_acquisition_plan = prepare_wfs_acquisition(
        CurvaturePackedAcquisition(channel_detector;
            readout_model=CurvatureChannelReadout()),
        rates, channel_observation,
    )
    return (; pupil, rates, frame_observation, channel_observation,
        optics_plan, frame_acquisition_plan, channel_acquisition_plan,
        rng=Xoshiro(0x43555256),
        frame_state=S4CurvatureExchangeState(T),
        channel_state=S4CurvatureExchangeState(T))
end

"""
Prepare one AOS Curvature plant and independent registered FGA paired-signal
estimators. AOS owns the two optical branches and complete detector acquisition;
FGA owns support, reference, branch-scale, and differential-signal calibration.
"""
function prepare_s4_curvature()
    Base.pkgversion(FilterGraphAlgorithms) == _FGA_VERSION_CLAIM || error(
        "S4 Curvature requires FilterGraphAlgorithms $_FGA_VERSION_CLAIM",
    )
    Base.pkgversion(JuliaFilterGraph) == _JFG_VERSION_CLAIM || error(
        "S4 Curvature requires JuliaFilterGraph $_JFG_VERSION_CLAIM",
    )
    plant = _prepare_plant()
    image = _image_owner()
    channels = _channel_owner()
    plant_image = _image_owner(; support=trues(2, 2),
        reference=zeros(Float32, 2, 2), branch_scales=(1.0f0, 1.0f0))
    plant_channels = _channel_owner(; support=trues(2, 2),
        reference=zeros(Float32, 2, 2), branch_scales=(1.0f0, 1.0f0))
    return (; plant..., frozen_fga_frame=zeros(Float32, 2, 4),
        plant_fga_frame=zeros(Float32, 2, 4), image, channels,
        plant_image, plant_channels, fga_tree_claim=_FGA_TREE_CLAIM,
        jfg_tree_claim=_JFG_TREE_CLAIM)
end

@inline function _model_timestamp_nanoseconds(sequence::UInt64,
    exposure_duration::AbstractFloat)
    period = round(Int64,
        Float64(exposure_duration) * Float64(_NANOSECONDS_PER_SECOND))
    return Int64(sequence) * period
end

@inline function _publish_observation!(state::S4CurvatureExchangeState,
    values::AbstractMatrix, layout::Symbol, exposure_duration::T,
) where {T<:AbstractFloat}
    _require_curvature_acquisition_ready(state)
    sequence = state.sequence + UInt64(1)
    state.sequence = sequence
    state.outstanding = true
    return S4CurvatureObservation(
        values,
        layout,
        _CALIBRATION_SIGNATURE,
        sequence,
        exposure_duration,
        _model_timestamp_nanoseconds(sequence, exposure_duration),
    )
end

@inline function _require_curvature_acquisition_ready(
    state::S4CurvatureExchangeState)
    state.blocked && throw(ArgumentError(
        "the Curvature exchange is blocked; reset is required"))
    state.outstanding && throw(ArgumentError(
        "the outstanding Curvature observation must be accepted or reset before acquisition"))
    return nothing
end

"""Form and acquire one complete packed Curvature detector image in AOS."""
@inline function produce_s4_curvature_frame!(prepared)
    _require_curvature_acquisition_ready(prepared.frame_state)
    form_wfs_optical_products!(prepared.rates, prepared.pupil,
        prepared.optics_plan)
    acquire_wfs_observation!(prepared.frame_observation, prepared.rates,
        prepared.frame_acquisition_plan, prepared.rng)
    return _publish_observation!(
        prepared.frame_state,
        observation_storage(prepared.frame_observation),
        _FRAME_LAYOUT,
        prepared.frame_acquisition_plan.detector_exposure_duration,
    )
end

"""Form and acquire one complete branch-by-channel Curvature readout in AOS."""
@inline function produce_s4_curvature_channels!(prepared)
    _require_curvature_acquisition_ready(prepared.channel_state)
    form_wfs_optical_products!(prepared.rates, prepared.pupil,
        prepared.optics_plan)
    acquire_wfs_observation!(prepared.channel_observation, prepared.rates,
        prepared.channel_acquisition_plan, prepared.rng)
    return _publish_observation!(
        prepared.channel_state,
        observation_storage(prepared.channel_observation),
        _CHANNEL_LAYOUT,
        prepared.channel_acquisition_plan.detector_exposure_duration,
    )
end

"""Transpose a complete AOS packed image and run the FGA paired-image plan."""
@inline function process_s4_curvature_image!(fga_frame, owner, aos_frame)
    transpose_aos_curvature_frame!(fga_frame, aos_frame)
    return process!(owner.signal, owner.workspace, owner.plan, fga_frame,
        owner.signature)
end

"""Run the FGA paired-channel plan on the complete AOS channel readout."""
@inline function process_s4_curvature_channels!(owner, aos_channels)
    return process!(owner.signal, owner.workspace, owner.plan, aos_channels,
        owner.signature)
end

@inline function _reject_curvature!(state::S4CurvatureExchangeState,
    disposition::S4CurvatureDisposition)
    state.blocked = true
    return disposition
end

@inline function _validate_curvature_observation!(
    state::S4CurvatureExchangeState,
    observation::S4CurvatureObservation,
    expected_layout::Symbol,
    expected_exposure_duration::AbstractFloat,
)
    state.blocked && return S4CurvatureRejectedBlocked
    state.outstanding || return _reject_curvature!(state,
        S4CurvatureRejectedNoOutstandingObservation)
    observation.sequence == state.sequence || return _reject_curvature!(state,
        S4CurvatureRejectedSequence)
    expected_timestamp = _model_timestamp_nanoseconds(
        observation.sequence, expected_exposure_duration)
    observation.model_timestamp_nanoseconds == expected_timestamp ||
        return _reject_curvature!(state,
            S4CurvatureRejectedModelTimestamp)
    isequal(observation.exposure_duration, expected_exposure_duration) ||
        return _reject_curvature!(state, S4CurvatureRejectedExposure)
    observation.layout === expected_layout || return _reject_curvature!(state,
        S4CurvatureRejectedLayout)
    observation.calibration_signature == _CALIBRATION_SIGNATURE ||
        return _reject_curvature!(state,
            S4CurvatureRejectedCalibrationSignature)
    return S4CurvatureAccepted
end

@inline function _accept_curvature!(state::S4CurvatureExchangeState,
    observation::S4CurvatureObservation)
    state.published_sequence = observation.sequence
    state.published_exposure_duration = observation.exposure_duration
    state.published_model_timestamp_nanoseconds =
        observation.model_timestamp_nanoseconds
    state.outstanding = false
    return S4CurvatureAccepted
end

"""Validate an acquired image's association and run the FGA image estimator."""
@inline function process_s4_curvature_image!(prepared,
    observation::S4CurvatureObservation)
    state = prepared.frame_state
    disposition = _validate_curvature_observation!(
        state,
        observation,
        _FRAME_LAYOUT,
        prepared.frame_acquisition_plan.detector_exposure_duration,
    )
    disposition === S4CurvatureAccepted || return disposition
    process_s4_curvature_image!(prepared.plant_fga_frame,
        prepared.plant_image, observation.values) === nothing ||
        return _reject_curvature!(state, S4CurvatureRejectedEstimator)
    return _accept_curvature!(state, observation)
end

"""Validate acquired channels' association and run the FGA channel estimator."""
@inline function process_s4_curvature_channels!(prepared,
    observation::S4CurvatureObservation)
    state = prepared.channel_state
    disposition = _validate_curvature_observation!(
        state,
        observation,
        _CHANNEL_LAYOUT,
        prepared.channel_acquisition_plan.detector_exposure_duration,
    )
    disposition === S4CurvatureAccepted || return disposition
    process_s4_curvature_channels!(prepared.plant_channels,
        observation.values) === nothing ||
        return _reject_curvature!(state, S4CurvatureRejectedEstimator)
    return _accept_curvature!(state, observation)
end

"""Run one associated AOS-image/FGA-estimator Curvature step."""
@inline function step_s4_curvature_frame!(prepared)
    observation = produce_s4_curvature_frame!(prepared)
    disposition = process_s4_curvature_image!(prepared, observation)
    disposition === S4CurvatureAccepted || error(
        "the internally formed Curvature image was rejected with $disposition")
    return observation.sequence
end

"""Run one associated AOS-channel/FGA-estimator Curvature step."""
@inline function step_s4_curvature_channels!(prepared)
    observation = produce_s4_curvature_channels!(prepared)
    disposition = process_s4_curvature_channels!(prepared, observation)
    disposition === S4CurvatureAccepted || error(
        "the internally formed Curvature channels were rejected with $disposition")
    return observation.sequence
end

@inline function _reset_curvature_state!(state::S4CurvatureExchangeState)
    state.sequence = UInt64(0)
    state.outstanding = false
    state.blocked = false
    state.published_sequence = UInt64(0)
    state.published_exposure_duration = zero(state.published_exposure_duration)
    state.published_model_timestamp_nanoseconds = Int64(0)
    return state
end

"""Reset Curvature frame/channel association state without replacing storage."""
function reset_s4_curvature_exchange!(prepared)
    _reset_curvature_state!(prepared.frame_state)
    _reset_curvature_state!(prepared.channel_state)
    return prepared
end

end # module AOSFGACurvature
