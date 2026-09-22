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
export transpose_aos_curvature_frame!, frozen_aos_curvature_frame
export frozen_aos_curvature_channels

const _FGA_VERSION_CLAIM = v"0.5.0"
const _JFG_VERSION_CLAIM = v"0.2.3"
const _FGA_TREE_CLAIM = "0ce4ceae4904446fba170ee25324d65b435c06e6"
const _JFG_TREE_CLAIM = "1893cea2d29f0cd89303939308e23ce7546bc240"
const _CALIBRATION_SIGNATURE = UInt64(0x43555256)

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
    2 6
    4 8
]

"""The same paired branches in AOS branch-by-channel order."""
frozen_aos_curvature_channels() = Float32[
    10 30 20 40
    2 6 4 8
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
        band=:custom, magnitude=zero(T), coordinates=(zero(T), zero(T)),
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
        rng=Xoshiro(0x43555256))
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

"""Form and acquire one complete packed Curvature detector image in AOS."""
@inline function produce_s4_curvature_frame!(prepared)
    form_wfs_optical_products!(prepared.rates, prepared.pupil,
        prepared.optics_plan)
    acquire_wfs_observation!(prepared.frame_observation, prepared.rates,
        prepared.frame_acquisition_plan, prepared.rng)
    return observation_storage(prepared.frame_observation)
end

"""Form and acquire one complete branch-by-channel Curvature readout in AOS."""
@inline function produce_s4_curvature_channels!(prepared)
    form_wfs_optical_products!(prepared.rates, prepared.pupil,
        prepared.optics_plan)
    acquire_wfs_observation!(prepared.channel_observation, prepared.rates,
        prepared.channel_acquisition_plan, prepared.rng)
    return observation_storage(prepared.channel_observation)
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

end # module AOSFGACurvature
