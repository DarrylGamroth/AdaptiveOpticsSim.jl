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

export DETECTOR_FRAME_SCHEMA, DETECTOR_FRAME_CALIBRATION_SIGNATURE
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
const DETECTOR_FRAME_CALIBRATION_SIGNATURE = UInt64(1)
const _RNG_SEED = UInt64(0x5331)
const _FRAME_PERIOD_NANOSECONDS = Int64(1_000_000)

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
    disturbance_opd::Matrix{Float32}
    fga_frame::Matrix{Float32}
    adopted_command::Vector{Float32}
    applied_command::Vector{Float32}
    state::S1LockstepState
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
        resolution=8,
        diameter=T(1),
        central_obstruction=zero(T),
        T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=:custom,
        wavelength=T(750e-9),
        photon_irradiance=T(2e8),
        T=T,
    )
    topology = SampledActuatorTopology(reshape(T[0, 0], 2, 1); T=T)
    dm = DeformableMirror(
        telescope;
        topology,
        influence_width=T(0.35),
        T=T,
    )
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=2,
        n_pix_subap=4,
        mode=Diffractive(),
        T=T,
    )
    rate = shack_hartmann_rate_map(sensor, pupil, source)
    optics = shack_hartmann_optics(sensor, source)
    optics_plan = prepare_wfs_optics(optics, pupil, rate)
    detector = Detector(
        noise=NoiseNone(),
        exposure_duration=T(1e-3),
        qe=one(T),
        response_model=NullFrameResponse(),
        T=T,
    )
    observation = WFSObservation(
        similar(intensity_values(rate));
        units=:electron_count,
        layout=:lenslet_mosaic,
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

function _calibrate_reconstructor!(prepared, graph, outputs, inputs)
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
    replace_parameters!(graph, Symbol("reference-slopes") => copy(outputs.slopes))

    poke = 2.0f-8
    positive = zeros(Float32, 8)
    negative = zeros(Float32, 8)
    prepared.adopted_command[1] = poke
    set_command!(prepared.dm, prepared.adopted_command)
    _form_plant_frame!(prepared)
    permutedims!(prepared.fga_frame, observation_storage(prepared.observation), (2, 1))
    JuliaFilterGraph.process!(
        outputs,
        graph,
        inputs,
        (image=SampleMetadata(2; terminal=true),),
    )
    _interleave_slopes!(positive, outputs.slopes)

    prepared.adopted_command[1] = -poke
    set_command!(prepared.dm, prepared.adopted_command)
    _form_plant_frame!(prepared)
    permutedims!(prepared.fga_frame, observation_storage(prepared.observation), (2, 1))
    JuliaFilterGraph.process!(
        outputs,
        graph,
        inputs,
        (image=SampleMetadata(3; terminal=true),),
    )
    _interleave_slopes!(negative, outputs.slopes)

    interaction = reshape((positive .- negative) ./ (2.0f0 * poke), 8, 1)
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
    return product
end

"""
    prepare_s1_lockstep(; disturbance_command=3f-8)

Prepare the deterministic S1 Float32 composition. AOS owns the physical
SHWFS and one-actuator plant, AdaptiveOpticsCalibration prepares the cold
reconstructor, and FGA/JFG own the per-frame estimator and RTC chain.
"""
function prepare_s1_lockstep(; disturbance_command::Float32=3.0f-8)
    plant = _prepare_plant()
    graph = prepare_graph(
        joinpath(@__DIR__, "s1-shwfs-f32.conf");
        algorithms=algorithms(),
    )
    graph.input_formats.image.schema == DETECTOR_FRAME_SCHEMA || error(
        "S1 detector-frame schema does not match the prepared FGA Graph",
    )
    graph.output_formats.demanded.schema == DEMANDED_PDM_COMMAND_V1 || error(
        "S1 demanded-PDM-command schema does not match FGA",
    )
    replace_parameters!(
        graph,
        Symbol("controller-to-vdm") => ones(Float32, 1, 1),
        Symbol("active-to-full") => ones(Float32, 1, 1),
        Symbol("vdm-to-pdm") => ones(Float32, 1, 1),
    )
    fga_frame = zeros(Float32, 8, 8)
    outputs = (
        demanded=zeros(Float32, 1),
        slopes=zeros(Float32, 4, 2),
        flux=zeros(Float32, 4),
        validity=fill(false, 4),
    )
    inputs = (image=fga_frame,)
    adopted_command = zeros(Float32, 1)
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
        zeros(Float32, 8, 8),
        fga_frame,
        adopted_command,
        similar(adopted_command),
        S1LockstepState(0, false, false),
    )
    _calibrate_reconstructor!(prepared, graph, outputs, inputs)
    JuliaFilterGraph.reset!(graph)

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
        DETECTOR_FRAME_CALIBRATION_SIGNATURE,
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
    frame.calibration_signature == DETECTOR_FRAME_CALIBRATION_SIGNATURE ||
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
