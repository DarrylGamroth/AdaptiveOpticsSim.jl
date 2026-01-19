using AdaptiveOptics
using Random
using Logging

rng = MersenneTwister(1)

sim = initialize_ao_shwfs(
    resolution=32,
    diameter=8.0,
    sampling_time=1e-3,
    r0=0.2,
    L0=25.0,
    fractional_r0=[0.7, 0.3],
    wind_speed=[5.0, 10.0],
    wind_direction=[0.0, 90.0],
    altitude=[0.0, 8000.0],
    n_act=4,
    n_subap=4,
)

imat = interaction_matrix(sim.dm, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(sim.dm.state.coefs)

n_half = div(length(cmd), 2)
for _ in 1:5
    advance!(sim.atm, sim.tel; rng=rng)
    propagate!(sim.atm, sim.tel)
    apply!(sim.dm, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, sim.wfs.state.slopes)
    sim.dm.state.coefs .= 0
    sim.dm.state.coefs[1:n_half] .= -cmd[1:n_half]
end

@info "Closed-loop first-stage run complete"
