using AdaptiveOptics
using Random

rng = MersenneTwister(0)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4)

imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
cmd = similar(dm.state.coefs)

n_iter = 5
for k in 1:n_iter
    advance!(atm, tel; rng=rng)
    propagate!(atm, tel)
    apply!(dm, tel; additive=true)
    measure!(wfs, tel)

    reconstruct!(cmd, recon, wfs.state.slopes)
    dm.state.coefs .= -cmd
end

println("Closed-loop demo complete")
