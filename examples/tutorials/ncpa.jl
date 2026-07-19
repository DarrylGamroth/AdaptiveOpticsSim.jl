include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    dm = DeformableMirror(tel; n_act=4, influence_width=0.35)
    atm = KolmogorovAtmosphere(tel; r0=0.18, L0=25.0)
    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 30e-9, -20e-9, 10e-9])
    pupil = PupilFunction(tel)
    apply_surface!(pupil, ncpa, DMReplace())
    imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
    image = copy(intensity_values(form_direct_image!(imaging)))

    @info "NCPA tutorial complete" opd_rms=pupil_rms(pupil.opd, pupil_support(pupil))
    return (
        ncpa_opd=copy(pupil.opd),
        image=image,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
