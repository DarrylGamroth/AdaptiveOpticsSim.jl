include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    dm = DeformableMirror(tel; n_act=4, influence_width=0.35)
    atm = KolmogorovAtmosphere(tel; r0=0.18, L0=25.0)
    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 30e-9, -20e-9, 10e-9])
    apply!(ncpa, tel, DMReplace())
    psf = copy(compute_psf!(tel, src; zero_padding=2))

    @info "NCPA tutorial complete" opd_rms=pupil_rms(tel.state.opd, tel.state.pupil)
    return (
        ncpa_opd=copy(tel.state.opd),
        psf=psf,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
