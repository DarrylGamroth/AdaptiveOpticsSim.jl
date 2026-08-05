include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(magnitude=6.0)
    model = GaussianDiskSourceModel(sigma_arcsec=0.4, n_side=5)
    ext = with_extended_source(src, model)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    pupil = PupilFunction(tel)
    @. pupil.opd = 5e-8 * zb.modes[:, :, 5]

    sh_point = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive())
    sh_ext = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive())
    point_peak = WavefrontSensors.sampled_spots_peak!(sh_point, pupil, src)
    ext_peak = WavefrontSensors.sampled_spots_peak!(sh_ext, pupil, ext)
    point_slopes = copy(measure!(sh_point, pupil, src))
    ext_slopes = copy(measure!(sh_ext, pupil, ext))

    pyr_point = PyramidWFS(tel; pupil_samples=6, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=6, mode=Diffractive(), modulation=1.0)
    pyr_point_detector = AdaptiveOpticsSim.Detectors.Detector(
        noise=AdaptiveOpticsSim.Detectors.NoiseNone(), qe=1.0,
        integration_time=1.0,
        response_model=AdaptiveOpticsSim.Detectors.NullFrameResponse())
    pyr_extended_detector = AdaptiveOpticsSim.Detectors.Detector(
        noise=AdaptiveOpticsSim.Detectors.NoiseNone(), qe=1.0,
        integration_time=1.0,
        response_model=AdaptiveOpticsSim.Detectors.NullFrameResponse())
    pyr_point_slopes = copy(measure!(pyr_point, pupil, src,
        pyr_point_detector))
    pyr_ext_slopes = copy(measure!(pyr_ext, pupil, ext,
        pyr_extended_detector))
    pyr_point_frame = copy(wfs_detector_image(pyr_point,
        pyr_point_detector))
    pyr_extended_frame = copy(wfs_detector_image(pyr_ext,
        pyr_extended_detector))

    point_rate = shack_hartmann_rate_map(sh_point, pupil, src)
    point_optics = prepare_wfs_optics(
        shack_hartmann_optics(sh_point, src), pupil, point_rate)
    form_wfs_optical_products!(point_rate, pupil, point_optics)
    extended_asterism = extended_source_asterism(ext)
    extended_rate = shack_hartmann_rate_map(
        sh_ext, pupil, extended_asterism)
    extended_optics = prepare_wfs_optics(
        shack_hartmann_optics(sh_ext, extended_asterism),
        pupil, extended_rate)
    form_wfs_optical_products!(
        extended_rate, pupil, extended_optics)
    point_spots = point_rate.values
    extended_spots = extended_rate.values
    sh_delta = copy(extended_spots .- point_spots)
    pyramid_frame_delta = pyr_extended_frame .- pyr_point_frame
    sh_relative_morphology = norm(sh_delta) / norm(point_spots)
    pyramid_relative_morphology = norm(pyramid_frame_delta) /
                                  norm(pyr_point_frame)
    @info(
        "Extended-source sensing tutorial complete",
        sh_rate_ratio=sum(extended_spots) / sum(point_spots),
        pyramid_count_ratio=sum(pyr_extended_frame) /
                            sum(pyr_point_frame),
        sh_relative_morphology=sh_relative_morphology,
        pyramid_relative_morphology=pyramid_relative_morphology,
    )
    return (
        sh_point_peak=point_peak,
        sh_extended_peak=ext_peak,
        sh_point_slopes=point_slopes,
        sh_extended_slopes=ext_slopes,
        sh_point_rate=sum(point_spots),
        sh_extended_rate=sum(extended_spots),
        sh_spot_delta=sh_delta,
        sh_relative_morphology=sh_relative_morphology,
        pyramid_point_slopes=pyr_point_slopes,
        pyramid_extended_slopes=pyr_ext_slopes,
        pyramid_point_counts=sum(pyr_point_frame),
        pyramid_extended_counts=sum(pyr_extended_frame),
        pyramid_frame_delta=pyramid_frame_delta,
        pyramid_relative_morphology=pyramid_relative_morphology,
        n_samples=length(extended_source_asterism(ext)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
