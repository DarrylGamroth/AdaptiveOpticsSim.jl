using AdaptiveOpticsSim
using LinearAlgebra
using Random

rng = MersenneTwister(0)

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
src = Source(band=:I, magnitude=0.0)

imat = interaction_matrix(dm, wfs, PupilFunction(tel), src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
pupil = PupilFunction(tel)
renderer = prepare_atmosphere_renderer(atm, tel, src)
command = similar(dm.state.coefs)

n_iter = 5
for k in 1:n_iter
    epoch = advance_by!(atm, 1e-3; rng=rng)
    render_atmosphere!(pupil, renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(command, recon, slopes(wfs))
    @. command = -command
    set_command!(dm, command)
end

println("Closed-loop demo complete, command_norm=$(norm(command))")
