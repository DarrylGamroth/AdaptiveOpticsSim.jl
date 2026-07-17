include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    dm = DeformableMirror(tel; n_act=4, influence_width=0.35)
    atm = KolmogorovAtmosphere(tel; r0=0.18, L0=25.0)
    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 30e-9, -20e-9, 10e-9])
    apply!(ncpa, tel, DMReplace())
    pupil = PupilFunction(tel)
    apply_opd!(pupil, opd_map(tel))
    imaging = prepare_direct_imaging(tel, pupil, src; zero_padding=2)
    image = copy(intensity_values(form_direct_image!(imaging)))

    @info "NCPA tutorial complete" opd_rms=pupil_rms(tel.state.opd, pupil_mask(tel))
    return (
        ncpa_opd=copy(tel.state.opd),
        image=image,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
