include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.1)
    src = base_source(magnitude=7.0)

    zb = ZernikeBasis(tel, 4)
    compute_zernike!(zb, tel)
    pupil = PupilFunction(tel)
    @. pupil.opd = 3e-8 * zb.modes[:, :, 4]

    sh = ShackHartmannWFS(tel; n_lenslets=6,
        pixel_scale_arcsec=0.06, n_pix_subap=8)
    rate = shack_hartmann_rate_map(sh, pupil, src)
    optics_plan = prepare_wfs_optics(
        shack_hartmann_optics(sh, src), pupil, rate)
    form_wfs_optical_products!(rate, pupil, optics_plan)

    layout = subaperture_layout(sh.front_end)
    metadata = AdaptiveOpticsSim.WavefrontSensors.wfs_output_metadata(sh)

    @info "Shack-Hartmann subaperture tutorial complete" n_valid=n_valid_subapertures(layout) pitch_m=layout.pitch_m total_photon_rate=sum(rate.values)
    return (
        n_valid=n_valid_subapertures(layout),
        pitch_m=layout.pitch_m,
        subap_pixels=layout.subap_pixels,
        valid_indices=copy(valid_subaperture_indices(layout)),
        metadata=metadata,
        photon_rate=copy(rate.values),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
