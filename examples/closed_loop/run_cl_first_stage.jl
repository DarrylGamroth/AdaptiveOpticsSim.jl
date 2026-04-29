using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(1)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[0.7, 0.3],
    wind_speed=[5.0, 10.0], wind_direction=[0.0, 90.0], altitude=[0.0, 8000.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(sim.optic, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(sim.optic.state.coefs)

n_half = div(length(cmd), 2)
for _ in 1:5
    advance!(sim.atm, sim.tel; rng=rng)
    propagate!(sim.atm, sim.tel)
    apply!(sim.optic, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, sim.wfs.state.slopes)
    sim.optic.state.coefs .= 0
    sim.optic.state.coefs[1:n_half] .= -cmd[1:n_half]
end

@info "Closed-loop first-stage run complete"
