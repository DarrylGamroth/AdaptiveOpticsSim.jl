using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(5)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[7.0], wind_direction=[10.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)

dm_coarse = DeformableMirror(tel; n_act=2, influence_width=0.6)
dm_fine = dm

imat_coarse = interaction_matrix(dm_coarse, wfs, pupil, src; amplitude=0.1)
imat_fine = interaction_matrix(dm_fine, wfs, pupil, src; amplitude=0.1)
recon_coarse = ModalReconstructor(imat_coarse; gain=0.3)
recon_fine = ModalReconstructor(imat_fine; gain=0.5)
cmd_coarse = similar(dm_coarse.state.coefs)
cmd_fine = similar(dm_fine.state.coefs)

for _ in 1:5
    epoch = advance_by!(atm, atmosphere_step; rng=rng)
    render_atmosphere!(pupil, atmosphere_renderer, atm, epoch)
    update_surface!(dm_coarse)
    update_surface!(dm_fine)
    apply_surface!(pupil, dm_coarse, DMAdditive())
    apply_surface!(pupil, dm_fine, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(cmd_coarse, recon_coarse, slopes(wfs))
    reconstruct!(cmd_fine, recon_fine, slopes(wfs))
    @. dm_coarse.state.coefs = -cmd_coarse
    @. dm_fine.state.coefs = -cmd_fine
end

@info "Closed-loop two stages complete"
