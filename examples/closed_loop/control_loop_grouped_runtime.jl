using AdaptiveOpticsSim
using Random

function build_branch(label::Symbol, seed::Integer)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    sim = AOSimulation(tel, src, atm, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    return ControlLoopBranch(label, sim, recon; wfs_detector=det, rng=MersenneTwister(seed))
end

cfg = GroupedControlLoopConfig(
    (:high, :low);
    name=:grouped_runtime_demo,
    outputs=GroupedRuntimeOutputRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
)

scenario = build_control_loop_scenario(cfg, build_branch(:high, 1), build_branch(:low, 2))
prepare!(scenario)
step!(scenario)

println("control_loop_grouped_runtime")
println("  name: ", control_loop_name(scenario))
println("  branch_labels: ", control_loop_branch_labels(scenario))
println("  command_length: ", length(command(scenario)))
println("  slopes_length: ", length(slopes(scenario)))
println("  grouped_wfs_stack_shape: ", size(grouped_wfs_stack(scenario)))
