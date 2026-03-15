include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, n_subap::Int=4)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(band=:R, magnitude=8.0)
    wfs = PyramidWFS(tel; n_subap=n_subap, mode=Diffractive(), modulation=3.0,
        modulation_points=8, diffraction_padding=2, n_pix_separation=2, n_pix_edge=1)

    zb = ZernikeBasis(tel, 8)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 2:5]
    gsc = GainSensingCamera(wfs.state.pyramid_mask, basis)

    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, src)
    calibrate!(gsc, calibration_frame)

    coeffs = [20e-9, -12e-9, 8e-9, 0.0]
    apply_opd!(tel, combine_modes(basis, coeffs))
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, tel, src)
    optical_gains = copy(compute_optical_gains!(gsc, frame))

    @info "Gain sensing camera tutorial complete" n_modes=length(optical_gains)
    return (
        coeffs=coeffs,
        calibration_frame=calibration_frame,
        frame=frame,
        optical_gains=optical_gains,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
