using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(1)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[0.7, 0.3],
    wind_speed=[5.0, 10.0], wind_direction=[0.0, 90.0], altitude=[0.0, 8000.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)
atmosphere_renderer = prepare_atmosphere_renderer(sim.atm, sim.tel, sim.src)
pupil = PupilFunction(sim.tel)

imat = interaction_matrix(sim.optic, sim.wfs, pupil; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(sim.optic.state.coefs)

n_half = div(length(cmd), 2)
for _ in 1:5
    epoch = advance_by!(sim.atm, atmosphere_step; rng=rng)
    render_atmosphere!(pupil, atmosphere_renderer, sim.atm, epoch)
    update_surface!(sim.optic)
    apply_surface!(pupil, sim.optic, DMAdditive())
    measure!(sim.wfs, pupil)
    reconstruct!(cmd, recon, slopes(sim.wfs))
    sim.optic.state.coefs .= 0
    sim.optic.state.coefs[1:n_half] .= -cmd[1:n_half]
end

@info "Closed-loop first-stage run complete"
