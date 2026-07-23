using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(6)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[7.0], wind_direction=[10.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)

dm_coarse = DeformableMirror(tel; n_act=2, influence_width=0.6)
dm_fine = dm
pupil = PupilFunction(tel)

imat_coarse = interaction_matrix(dm_coarse, wfs, pupil, src; amplitude=0.1)
imat_fine = interaction_matrix(dm_fine, wfs, pupil, src; amplitude=0.1)
recon_coarse = ModalReconstructor(imat_coarse; gain=0.3)
recon_fine = ModalReconstructor(imat_fine; gain=0.5)
cmd_coarse = similar(dm_coarse.state.coefs)
cmd_fine = similar(dm_fine.state.coefs)

function run_segment!(atm, renderer, pupil, source, wfs, dm_coarse,
    dm_fine, recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, n_iter)
    for _ in 1:n_iter
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(pupil, renderer, atm, epoch)
        update_surface!(dm_coarse)
        update_surface!(dm_fine)
        apply_surface!(pupil, dm_coarse, DMAdditive())
        apply_surface!(pupil, dm_fine, DMAdditive())
        measure!(wfs, pupil, source)
        reconstruct!(cmd_coarse, recon_coarse, slopes(wfs))
        reconstruct!(cmd_fine, recon_fine, slopes(wfs))
        @. dm_coarse.state.coefs = -cmd_coarse
        @. dm_fine.state.coefs = -cmd_fine
    end
    return nothing
end

renderer = prepare_atmosphere_renderer(atm, tel, src)
run_segment!(atm, renderer, pupil, src, wfs, dm_coarse, dm_fine,
    recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, 3)

changed_atm = KolmogorovAtmosphere(tel; r0=0.1, L0=25.0)
changed_renderer = prepare_atmosphere_renderer(changed_atm, tel, src)
@info "Atmosphere r0 updated to 0.1"
run_segment!(changed_atm, changed_renderer, pupil, src, wfs, dm_coarse,
    dm_fine, recon_coarse, recon_fine, cmd_coarse, cmd_fine, rng, 3)

@info "Closed-loop two stages with atmosphere change complete"
