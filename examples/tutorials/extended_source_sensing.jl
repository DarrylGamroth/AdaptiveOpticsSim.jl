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

    sh_point = ShackHartmannWFS(tel; n_lenslets=6)
    sh_ext = ShackHartmannWFS(tel; n_lenslets=6)

    pyr_point = PyramidWFS(tel; pupil_samples=6, modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=6, modulation=1.0)

    pyramid_point_front_end = PyramidOpticalFrontEnd(pyr_point, src)
    pyramid_point_rate = pyramid_rate_map(pyramid_point_front_end, pupil)
    pyramid_point_optics = prepare_wfs_optics(
        pyramid_point_front_end, pupil, pyramid_point_rate)
    form_wfs_optical_products!(
        pyramid_point_rate, pupil, pyramid_point_optics)

    extended_paths = extended_source_asterism(ext)
    extended_pupils = ntuple(_ -> pupil, length(extended_paths))
    pyramid_extended_front_end = PyramidOpticalFrontEnd(pyr_ext, ext)
    pyramid_extended_rates = pyramid_rate_map(
        pyramid_extended_front_end, extended_pupils)
    pyramid_extended_optics = prepare_wfs_optics(
        pyramid_extended_front_end, extended_pupils,
        pyramid_extended_rates)
    form_wfs_optical_products!(pyramid_extended_rates, extended_pupils,
        pyramid_extended_optics)

    pyramid_point_frame = copy(intensity_values(pyramid_point_rate))
    pyramid_extended_frame = zeros(eltype(pyramid_point_frame),
        size(pyramid_point_frame))
    for rate in pyramid_extended_rates
        pyramid_extended_frame .+= intensity_values(rate)
    end

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
    point_peak = maximum(point_spots)
    ext_peak = maximum(extended_spots)
    sh_delta = copy(extended_spots .- point_spots)
    pyramid_frame_delta = pyramid_extended_frame .- pyramid_point_frame
    sh_relative_morphology = norm(sh_delta) / norm(point_spots)
    pyramid_relative_morphology = norm(pyramid_frame_delta) /
                                  norm(pyramid_point_frame)
    @info(
        "Extended-source sensing tutorial complete",
        sh_rate_ratio=sum(extended_spots) / sum(point_spots),
        pyramid_rate_ratio=sum(pyramid_extended_frame) /
                           sum(pyramid_point_frame),
        sh_relative_morphology=sh_relative_morphology,
        pyramid_relative_morphology=pyramid_relative_morphology,
    )
    return (
        sh_point_peak=point_peak,
        sh_extended_peak=ext_peak,
        sh_point_rate=sum(point_spots),
        sh_extended_rate=sum(extended_spots),
        sh_spot_delta=sh_delta,
        sh_relative_morphology=sh_relative_morphology,
        pyramid_point_rate=sum(pyramid_point_frame),
        pyramid_extended_rate=sum(pyramid_extended_frame),
        pyramid_frame_delta=pyramid_frame_delta,
        pyramid_relative_morphology=pyramid_relative_morphology,
        n_samples=length(extended_source_asterism(ext)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
