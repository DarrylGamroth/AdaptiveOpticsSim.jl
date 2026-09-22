module AOSFGAZernike

using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using FilterGraphAlgorithms
using JuliaFilterGraph
using Random

import FilterGraphAlgorithms: process!

export prepare_s4_zernike, produce_s4_zernike_frame!, process_s4_zernike!
export transpose_aos_zernike_frame!, frozen_aos_zernike_frame

const _FGA_VERSION_CLAIM = v"0.4.0"
const _JFG_VERSION_CLAIM = v"0.2.2"
const _FGA_TREE_CLAIM = "598a10402690250190ad1dd976b03577281e6388"
const _JFG_TREE_CLAIM = "edaa8b14e2a2b9266cff5674ecbe9c9f6c20c17c"
const _CALIBRATION_SIGNATURE = UInt64(781)

const _FGA_SUPPORT = Bool[
    true false true
    true true false
    false true true
]

const _FGA_REFERENCE = Float32[
    0.05 0.0 0.1
    -0.2 0.15 0.0
    0.0 -0.05 0.25
]

"""Preallocated AOS `(x, y)` to FGA `(row=y, column=x)` frame transfer."""
@inline function transpose_aos_zernike_frame!(
    destination::AbstractMatrix{T}, source::AbstractMatrix{T},
) where {T}
    size(destination) == reverse(size(source)) || throw(DimensionMismatch(
        "AOS Zernike source and FGA destination axes are incompatible",
    ))
    permutedims!(destination, source, (2, 1))
    return destination
end

"""Asymmetric AOS `(x, y)` frame whose transpose is the frozen FGA oracle."""
frozen_aos_zernike_frame() = Float32[
    2 3 17
    5 7 19
    11 13 23
]

function _signal_owner(normalization::ZernikePupilSignalNormalization;
    image_size=(3, 3), support=_FGA_SUPPORT,
    reference=_FGA_REFERENCE, incidence_flux=0.0f0,
)
    plan = ZernikePupilSignalPlan(
        image_size, support, reference, normalization,
        _CALIBRATION_SIGNATURE, incidence_flux,
    )
    return (; plan, workspace=ZernikePupilSignalWorkspace(plan),
        signal=zeros(Float32, count(support)), divisor=zeros(Float32, 1),
        signature=UInt64[_CALIBRATION_SIGNATURE])
end

function _prepare_plant()
    T = Float32
    telescope = Telescope(
        resolution=12, diameter=T(6), central_obstruction=zero(T),
        fov_arcsec=zero(T), pupil_reflectivity=one(T), T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=:custom, magnitude=zero(T), coordinates=(zero(T), zero(T)),
        wavelength=750.0f-9, photon_irradiance=12.0f0,
        radiometry=PhysicalPhotonIrradianceSource(), T=T,
    )
    sensor = ZernikeWFS(
        telescope; pupil_samples=3, diffraction_padding=2, T=T,
    )
    front_end = ZernikeOpticalFrontEnd(sensor, source)
    rate = zernike_rate_map(front_end, pupil)
    optics_plan = prepare_wfs_optics(front_end, pupil, rate)
    detector = Detector(
        noise=NoiseNone(), exposure_duration=0.375f0, qe=one(T),
        response_model=NullFrameResponse(), T=T,
    )
    observation = WFSObservation(
        similar(intensity_values(rate)); units=:electron_count,
        layout=:zernike_pupil_image,
    )
    acquisition_plan = prepare_wfs_acquisition(
        detector, rate, observation; source,
    )
    return (; pupil, rate, observation, optics_plan, acquisition_plan,
        rng=Xoshiro(0x5a45524e))
end

"""
Prepare a physical AOS Zernike plant and independent registered FGA pupil-signal
estimators. AOS owns phase-spot optics and detector acquisition; FGA owns
reference subtraction and normalization.
"""
function prepare_s4_zernike()
    Base.pkgversion(FilterGraphAlgorithms) == _FGA_VERSION_CLAIM || error(
        "S4 Zernike requires FilterGraphAlgorithms $_FGA_VERSION_CLAIM",
    )
    Base.pkgversion(JuliaFilterGraph) == _JFG_VERSION_CLAIM || error(
        "S4 Zernike requires JuliaFilterGraph $_JFG_VERSION_CLAIM",
    )
    plant = _prepare_plant()
    mean = _signal_owner(ZernikePupilSignalMeanValidFlux)
    incidence = _signal_owner(
        ZernikePupilSignalIncidenceFlux; incidence_flux=13.5f0,
    )
    plant_signal = _signal_owner(
        ZernikePupilSignalMeanValidFlux; support=trues(3, 3),
        reference=zeros(Float32, 3, 3),
    )
    return (; plant..., frozen_fga_frame=zeros(Float32, 3, 3),
        plant_fga_frame=zeros(Float32, 3, 3), mean, incidence,
        plant_signal, fga_tree_claim=_FGA_TREE_CLAIM,
        jfg_tree_claim=_JFG_TREE_CLAIM)
end

"""Form one actual AOS acquired Zernike pupil frame without estimation."""
@inline function produce_s4_zernike_frame!(prepared)
    form_wfs_optical_products!(prepared.rate, prepared.pupil,
        prepared.optics_plan)
    acquire_wfs_observation!(prepared.observation, prepared.rate,
        prepared.acquisition_plan, prepared.rng)
    return observation_storage(prepared.observation)
end

"""Transfer a complete AOS pupil frame into FGA and run one FGA signal plan."""
@inline function process_s4_zernike!(fga_frame, owner, aos_frame)
    transpose_aos_zernike_frame!(fga_frame, aos_frame)
    return process!(owner.signal, owner.divisor, owner.workspace, owner.plan,
        fga_frame, owner.signature)
end

end # module AOSFGAZernike
