include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra

function source_pupil(pupil, source, reference)
    component = PupilFunction(pupil)
    direction = Optics.coordinates_xy_arcsec(source)
    pointing = Optics.coordinates_xy_arcsec(reference)
    # The sensor tracks the reference source. Supply the complete pupil OPD,
    # including W(r) = (θ_source - θ_reference) ⋅ r, with angles in radians.
    angle = (x=(direction[1] - pointing[1]) * pi / 648000,
        y=(direction[2] - pointing[2]) * pi / 648000)
    metadata = pupil.metadata
    slopes = ntuple(2) do d
        metadata.orientation.signs[d] *
            getproperty(angle, metadata.orientation.axes[d])
    end
    for j in axes(component.opd, 2), i in axes(component.opd, 1)
        r1 = metadata.origin[1] + (i - 1) * metadata.sampling[1]
        r2 = metadata.origin[2] + (j - 1) * metadata.sampling[2]
        component.opd[i, j] += slopes[1] * r1 + slopes[2] * r2
    end
    return component
end

function main(; resolution::Int=24,
    model=GaussianDiskSourceModel(sigma_arcsec=0.04, n_side=5))
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(magnitude=6.0)
    # At this grid spacing, the default ±0.08 arcsec source samples are below
    # the I-band pupil angular Nyquist limit (about 0.24 arcsec).
    ext = with_extended_source(src, model)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    pupil = PupilFunction(tel)
    @. pupil.opd = 5e-8 * zb.modes[:, :, 5]

    # Retain the full padded FFT footprint: a cropped spot window loses
    # photons when source tilt moves light toward its edge.
    spot_samples = 2 * div(resolution, 6)
    sh_point = ShackHartmannWFS(tel; n_lenslets=6, n_pix_subap=spot_samples,
        diffraction_padding=2)
    sh_ext = ShackHartmannWFS(tel; n_lenslets=6, n_pix_subap=spot_samples,
        diffraction_padding=2)

    pyr_point = PyramidWFS(tel; pupil_samples=6, modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=6, modulation=1.0)

    pyramid_point_front_end = PyramidOpticalFrontEnd(pyr_point, src)
    pyramid_point_rate = pyramid_rate_map(pyramid_point_front_end, pupil)
    pyramid_point_optics = prepare_wfs_optics(
        pyramid_point_front_end, pupil, pyramid_point_rate)
    form_wfs_optical_products!(
        pyramid_point_rate, pupil, pyramid_point_optics)

    extended_paths = extended_source_asterism(ext).sources
    extended_pupils = ntuple(length(extended_paths)) do i
        source_pupil(pupil, extended_paths[i], src)
    end
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
    point_spots = point_rate.values
    extended_spots = zeros(eltype(point_spots), size(point_spots))
    # Source samples are mutually incoherent: sum photon rates, not fields.
    # Each quadrature source already carries its normalized flux weight.
    for i in eachindex(extended_pupils)
        component_pupil = extended_pupils[i]
        component_source = extended_paths[i]
        component_rate = shack_hartmann_rate_map(
            sh_ext, component_pupil, component_source)
        component_optics = prepare_wfs_optics(
            shack_hartmann_optics(sh_ext, component_source),
            component_pupil, component_rate)
        form_wfs_optical_products!(
            component_rate, component_pupil, component_optics)
        extended_spots .+= component_rate.values
    end
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
