module AOSFGABiOEdge

using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using FilterGraphAlgorithms
using JuliaFilterGraph
using Random

import FilterGraphAlgorithms: process!

export prepare_s4_bi_o_edge, produce_s4_bi_o_edge_frame!, process_s4_bi_o_edge!
export transpose_aos_bi_o_edge_frame!, frozen_aos_bi_o_edge_frame

const _FGA_VERSION_CLAIM = v"0.4.0"
const _JFG_VERSION_CLAIM = v"0.2.2"
const _ORIGINS = ((0, 0), (2, 0), (2, 2), (0, 2))
const _SUPPORT = Bool[true true; false true]
const _CALIBRATION_SIGNATURE = UInt64(7177611906121727310)

"""Preallocated AOS `(x, y)` to FGA `(row=y, column=x)` frame transfer."""
@inline function transpose_aos_bi_o_edge_frame!(
    destination::AbstractMatrix{T}, source::AbstractMatrix{T},
) where {T}
    size(destination) == reverse(size(source)) || throw(DimensionMismatch(
        "AOS Bi-O-edge source and FGA destination axes are incompatible",
    ))
    permutedims!(destination, source, (2, 1))
    return destination
end

"""Asymmetric AOS `(x, y)` frame whose transpose is the frozen FGA oracle."""
frozen_aos_bi_o_edge_frame() = Float32[
    10 20 1 2
    30 40 3 4
    5 6 9 10
    7 8 11 12
]

function _image_owner(
    normalization_policy::AbstractString;
    image_size=(4, 4),
    pupil_size=(2, 2),
    origins=_ORIGINS,
    support=_SUPPORT,
    reference=Float32[0.1 -0.1; 0.3 -0.3; 0.4 -0.4],
    gain=Float32[1.5 0.25; 0.5 1.25; 2.0 0.75],
    incidence_flux=0.0f0,
)
    plan = BiOEdgeImagePlan(
        image_size,
        pupil_size,
        origins,
        support,
        reference,
        gain,
        normalization_policy == "mean-valid-flux" ? BiOEdgeImageMeanValidFlux :
        BiOEdgeImageIncidenceFlux,
        _CALIBRATION_SIGNATURE,
        incidence_flux,
    )
    return (; plan, workspace=BiOEdgeImageWorkspace(plan),
        signal=zeros(Float32, size(reference)), divisor=zeros(Float32, 1),
        signature=UInt64[_CALIBRATION_SIGNATURE])
end

function _prepare_plant()
    T = Float32
    telescope = Telescope(
        resolution=4, diameter=one(T), central_obstruction=zero(T),
        fov_arcsec=zero(T), pupil_reflectivity=one(T), T=T,
    )
    pupil = PupilFunction(telescope; T=T)
    source = Source(
        band=:custom, magnitude=zero(T), coordinates=(zero(T), zero(T)),
        wavelength=750.0f-9, photon_irradiance=2.0f8,
        radiometry=PhysicalPhotonIrradianceSource(), T=T,
    )
    sensor = BiOEdgeWFS(
        telescope; pupil_samples=2, modulation=zero(T), diffraction_padding=2, T=T,
    )
    front_end = BiOEdgeOpticalFrontEnd(sensor, source)
    rate = bi_o_edge_rate_map(front_end, pupil)
    optics_plan = prepare_wfs_optics(front_end, pupil, rate)
    detector = Detector(
        noise=NoiseNone(), exposure_duration=1.0f-3, qe=one(T),
        response_model=NullFrameResponse(), T=T,
    )
    observation = WFSObservation(
        similar(intensity_values(rate)); units=:electron_count, layout=:four_pupil_mosaic,
    )
    acquisition_plan = prepare_wfs_acquisition(detector, rate, observation; source)
    return (; pupil, rate, observation, optics_plan, acquisition_plan, rng=Xoshiro(0xB10E))
end

"""
Prepare a physical AOS Bi-O-edge plant and the independent FGA 0.4.0 image
estimator. AOS owns detector formation/acquisition; FGA owns differential
signal calibration and normalization.
"""
function prepare_s4_bi_o_edge()
    Base.pkgversion(FilterGraphAlgorithms) == _FGA_VERSION_CLAIM || error(
        "S4 Bi-O-edge requires FilterGraphAlgorithms $_FGA_VERSION_CLAIM",
    )
    Base.pkgversion(JuliaFilterGraph) == _JFG_VERSION_CLAIM || error(
        "S4 Bi-O-edge requires JuliaFilterGraph $_JFG_VERSION_CLAIM",
    )
    plant = _prepare_plant()
    mean = _image_owner("mean-valid-flux")
    incidence = _image_owner("incidence-flux"; incidence_flux=50.0f0)
    plant_image = _image_owner(
        "mean-valid-flux";
        image_size=(8, 8), pupil_size=(4, 4),
        origins=((0, 0), (4, 0), (4, 4), (0, 4)), support=trues(4, 4),
        reference=zeros(Float32, 16, 2), gain=ones(Float32, 16, 2),
    )
    return (; plant..., frozen_fga_frame=zeros(Float32, 4, 4),
        plant_fga_frame=zeros(Float32, 8, 8), mean, incidence, plant_image)
end

"""Form one actual AOS acquired Bi-O-edge detector frame without estimation."""
@inline function produce_s4_bi_o_edge_frame!(prepared)
    form_wfs_optical_products!(prepared.rate, prepared.pupil, prepared.optics_plan)
    acquire_wfs_observation!(
        prepared.observation, prepared.rate, prepared.acquisition_plan, prepared.rng,
    )
    return observation_storage(prepared.observation)
end

"""Transfer a complete AOS frame into FGA and run one chosen FGA image plan."""
@inline function process_s4_bi_o_edge!(fga_frame, owner, aos_frame)
    transpose_aos_bi_o_edge_frame!(fga_frame, aos_frame)
    return process!(
        owner.signal, owner.divisor, owner.workspace, owner.plan,
        fga_frame, owner.signature,
    )
end

end # module AOSFGABiOEdge
