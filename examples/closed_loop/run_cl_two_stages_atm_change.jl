using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(6)

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[7.0], wind_direction=[10.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
sim = AOSimulation(tel, src, atm, dm, wfs)

dm_coarse = DeformableMirror(sim.tel; n_act=2, influence_width=0.6)
dm_fine = sim.optic

imat_coarse = interaction_matrix(dm_coarse, sim.wfs, sim.tel; amplitude=0.1)
imat_fine = interaction_matrix(dm_fine, sim.wfs, sim.tel; amplitude=0.1)
recon_coarse = ModalReconstructor(imat_coarse; gain=0.3)
recon_fine = ModalReconstructor(imat_fine; gain=0.5)
cmd_coarse = similar(dm_coarse.state.coefs)
cmd_fine = similar(dm_fine.state.coefs)

function run_segment!(atm, renderer, atmosphere_output, sim, dm_coarse,
    dm_fine, recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, n_iter)
    for _ in 1:n_iter
        epoch = advance_by!(atm, sim.tel.params.sampling_time; rng=rng)
        render_atmosphere!(atmosphere_output, renderer, atm, epoch)
        copyto!(sim.tel.state.opd, atmosphere_output.opd)
        apply!(dm_coarse, sim.tel, DMAdditive())
        apply!(dm_fine, sim.tel, DMAdditive())
        measure!(sim.wfs, sim.tel)
        reconstruct!(cmd_coarse, recon_coarse, sim.wfs.state.slopes)
        reconstruct!(cmd_fine, recon_fine, sim.wfs.state.slopes)
        dm_coarse.state.coefs .= -cmd_coarse
        dm_fine.state.coefs .= -cmd_fine
    end
    return nothing
end

atmosphere_output = PupilFunction(sim.tel)
renderer = prepare_atmosphere_renderer(sim.atm, sim.tel, sim.src)
run_segment!(sim.atm, renderer, atmosphere_output, sim, dm_coarse, dm_fine,
    recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, 3)

changed_atm = KolmogorovAtmosphere(sim.tel; r0=0.1, L0=25.0)
changed_renderer = prepare_atmosphere_renderer(changed_atm, sim.tel, sim.src)
@info "Atmosphere r0 updated to 0.1"
run_segment!(changed_atm, changed_renderer, atmosphere_output, sim, dm_coarse,
    dm_fine, recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, 3)

@info "Closed-loop two stages with atmosphere change complete"
