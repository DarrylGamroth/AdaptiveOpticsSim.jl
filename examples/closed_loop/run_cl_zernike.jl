using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Logging

rng = MersenneTwister(4)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ZernikeWFS(tel; n_subap=8, diffraction_padding=2)
det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=1e-8)
recon = ModalReconstructor(imat; gain=0.5)
runtime = ClosedLoopRuntime(sim, recon; rng=rng, wfs_detector=det)
interface = simulation_interface(runtime)

prepare!(interface)
for _ in 1:5
    step!(interface)
end

@info "Closed-loop ZernikeWFS tutorial complete" residual_norm=norm(command(readout(interface))) slope_norm=norm(slopes(readout(interface)))
