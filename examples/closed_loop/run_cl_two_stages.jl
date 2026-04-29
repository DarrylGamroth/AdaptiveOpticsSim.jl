using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(5)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[7.0], wind_direction=[10.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)

dm_coarse = DeformableMirror(sim.tel; n_act=2, influence_width=0.6)
dm_fine = sim.optic

imat_coarse = interaction_matrix(dm_coarse, sim.wfs, sim.tel; amplitude=0.1)
imat_fine = interaction_matrix(dm_fine, sim.wfs, sim.tel; amplitude=0.1)
recon_coarse = ModalReconstructor(imat_coarse; gain=0.3)
recon_fine = ModalReconstructor(imat_fine; gain=0.5)
cmd_coarse = similar(dm_coarse.state.coefs)
cmd_fine = similar(dm_fine.state.coefs)

for _ in 1:5
    advance!(sim.atm, sim.tel; rng=rng)
    propagate!(sim.atm, sim.tel)
    apply!(dm_coarse, sim.tel, DMAdditive())
    apply!(dm_fine, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd_coarse, recon_coarse, sim.wfs.state.slopes)
    reconstruct!(cmd_fine, recon_fine, sim.wfs.state.slopes)
    dm_coarse.state.coefs .= -cmd_coarse
    dm_fine.state.coefs .= -cmd_fine
end

@info "Closed-loop two stages complete"
