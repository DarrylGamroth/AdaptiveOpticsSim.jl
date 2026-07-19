include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=32, zero_padding::Int=2, n_modes::Int=12,
    mode_id::Int=6, amplitude::Real=80e-9)
    tel = base_telescope(resolution=resolution)
    src = base_source()
    pupil = PupilFunction(tel)
    imaging = prepare_direct_imaging(pupil, src; zero_padding=zero_padding)
    nominal_image = copy(intensity_values(form_direct_image!(imaging)))

    zb = ZernikeBasis(tel, n_modes)
    compute_zernike!(zb, tel)
    apply_opd!(pupil, zb.modes[:, :, mode_id] .* amplitude)
    aberrated_image = copy(intensity_values(form_direct_image!(imaging)))
    pixel_scale = focal_plane_pixel_scale_arcsec(direct_imaging_output(imaging))

    @info "Image-formation tutorial complete" mode_id amplitude pixel_scale
    return (
        nominal_image=nominal_image,
        aberrated_image=aberrated_image,
        pixel_scale_arcsec=pixel_scale,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
