using AdaptiveOpticsSim
using LinearAlgebra
using Random

rng = MersenneTwister(0)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4)
src = Source(band=:I, magnitude=0.0)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
branch = RuntimeBranch(:main, sim, recon; rng=rng)
cfg = SingleRuntimeConfig(name=:closed_loop_demo, branch_label=:main)
scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)

n_iter = 5
for k in 1:n_iter
    step!(scenario)
end

println("Closed-loop demo complete, command_norm=$(norm(command(readout(scenario))))")
