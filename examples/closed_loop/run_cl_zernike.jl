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
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=1e-8)
recon = ModalReconstructor(imat; gain=0.5)
branch = RuntimeBranch(:main, sim, recon; rng=rng, wfs_detector=det)
cfg = SingleRuntimeConfig(
    name=:run_cl_zernike_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true),
)
scenario = build_runtime_scenario(cfg, branch)

prepare!(scenario)
for _ in 1:5
    step!(scenario)
end

rt = readout(scenario)
@info "Closed-loop ZernikeWFS tutorial complete" residual_norm=norm(command(rt)) slope_norm=norm(slopes(rt))
