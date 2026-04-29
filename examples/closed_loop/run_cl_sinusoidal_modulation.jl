using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(4)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[6.0], wind_direction=[30.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(sim.optic, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(sim.optic.state.coefs)

n_iter = 6
for k in 1:n_iter
    advance!(sim.atm, sim.tel; rng=rng)
    propagate!(sim.atm, sim.tel)
    apply!(sim.optic, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, sim.wfs.state.slopes)
    modulation = 1 + 0.2 * sin(2 * pi * k / n_iter)
    sim.optic.state.coefs .= -cmd .* modulation
end

@info "Closed-loop sinusoidal modulation complete"
