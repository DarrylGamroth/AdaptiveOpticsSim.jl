using AdaptiveOpticsSim
using Random

tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ShackHartmann(tel; n_subap=4)
sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
science_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

branch = RuntimeBranch(
    :main,
    sim,
    recon;
    science_detector=science_det,
    rng=MersenneTwister(1),
)

cfg = SingleRuntimeConfig(
    name=:single_runtime_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=false, science_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)
step!(scenario)

println("platform_single_runtime")
println("  name: ", platform_name(scenario))
println("  branch_labels: ", platform_branch_labels(scenario))
println("  command_length: ", length(command(scenario)))
println("  slopes_length: ", length(slopes(scenario)))
println("  science_frame_shape: ", size(science_frame(scenario)))
