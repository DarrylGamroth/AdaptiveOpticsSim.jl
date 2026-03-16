using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(3)

sim = initialize_ao_shwfs(
    resolution=32,
    diameter=8.0,
    sampling_time=1e-3,
    r0=0.2,
    L0=25.0,
    fractional_r0=[1.0],
    wind_speed=[5.0],
    wind_direction=[0.0],
    altitude=[0.0],
    n_act=4,
    n_subap=4,
)

imat1 = interaction_matrix(sim.dm, sim.wfs, sim.tel; amplitude=0.05)
imat2 = interaction_matrix(sim.dm, sim.wfs, sim.tel; amplitude=0.1)
mat = 0.5 .* (imat1.matrix .+ imat2.matrix)
recon = ModalReconstructor(InteractionMatrix(mat, 0.1); gain=0.4)
cmd = similar(sim.dm.state.coefs)

for _ in 1:5
    advance!(sim.atm, sim.tel; rng=rng)
    propagate!(sim.atm, sim.tel)
    apply!(sim.dm, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, sim.wfs.state.slopes)
    sim.dm.state.coefs .= -cmd
end

@info "Closed-loop long push-pull complete"
