using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Logging

rng = MersenneTwister(0)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[5.0], wind_direction=[0.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(sim.optic, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = ControlLoopBranch(:main, sim, recon; rng=rng)
cfg = SingleControlLoopConfig(name=:run_cl_demo, branch_label=:main)
scenario = build_control_loop_scenario(cfg, branch)
prepare!(scenario)

for _ in 1:5
    step!(scenario)
end

@info "Closed-loop run_cl complete" residual_norm=norm(command(readout(scenario)))
