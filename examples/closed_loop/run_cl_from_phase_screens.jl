using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(2)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[8.0], wind_direction=[45.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)

screens = Vector{Matrix{Float64}}(undef, 5)
for k in 1:5
    advance!(sim.atm, sim.tel; rng=rng)
    screens[k] = copy(sim.atm.state.opd)
end

imat = interaction_matrix(sim.optic, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
cmd = similar(sim.optic.state.coefs)

for k in 1:length(screens)
    sim.tel.state.opd .= screens[k] .* sim.tel.state.pupil
    apply!(sim.optic, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, sim.wfs.state.slopes)
    sim.optic.state.coefs .= -cmd
end

@info "Closed-loop from phase screens complete"
