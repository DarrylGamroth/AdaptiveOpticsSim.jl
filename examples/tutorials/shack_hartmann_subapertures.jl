include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.1)
    src = base_source(magnitude=7.0)

    zb = ZernikeBasis(tel, 4)
    compute_zernike!(zb, tel)
    @. tel.state.opd = 3e-8 * zb.modes[:, :, 4]

    sh = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8)
    prepare_runtime_wfs!(sh, tel, src)
    slopes = copy(measure!(sh, tel, src))

    layout = subaperture_layout(sh.front_end)
    calibration = subaperture_calibration(sh)
    metadata = AdaptiveOpticsSim.wfs_output_metadata(sh)

    @info "Shack-Hartmann subaperture tutorial complete" n_valid=n_valid_subapertures(layout) pitch_m=layout.pitch_m centroid_response=calibration.centroid_response
    return (
        n_valid=n_valid_subapertures(layout),
        pitch_m=layout.pitch_m,
        subap_pixels=layout.subap_pixels,
        valid_indices=copy(valid_subaperture_indices(layout)),
        centroid_response=calibration.centroid_response,
        calibrated=calibration.calibrated,
        threshold=slope_extraction_model(sh).threshold,
        metadata=metadata,
        slopes=slopes,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
