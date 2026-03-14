include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=32, zero_padding::Int=2, n_modes::Int=12,
    mode_id::Int=6, amplitude::Real=80e-9)
    tel = base_telescope(resolution=resolution)
    src = base_source()
    psf_nominal = copy(compute_psf!(tel, src; zero_padding=zero_padding))

    zb = ZernikeBasis(tel, n_modes)
    compute_zernike!(zb, tel)
    apply_opd!(tel, zb.modes[:, :, mode_id] .* amplitude)
    psf_aberrated = copy(compute_psf!(tel, src; zero_padding=zero_padding))
    pixel_scale = psf_pixel_scale_arcsec(tel, src, zero_padding)

    @info "Image-formation tutorial complete" mode_id amplitude pixel_scale
    return (
        psf_nominal=psf_nominal,
        psf_aberrated=psf_aberrated,
        pixel_scale_arcsec=pixel_scale,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
