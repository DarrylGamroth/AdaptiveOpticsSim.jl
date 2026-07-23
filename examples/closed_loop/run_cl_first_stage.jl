using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(1)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[0.7, 0.3],
    wind_speed=[5.0, 10.0], wind_direction=[0.0, 90.0], altitude=[0.0, 8000.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)

imat = interaction_matrix(dm, wfs, pupil, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(dm.state.coefs)

n_half = div(length(cmd), 2)
for _ in 1:5
    epoch = advance_by!(atm, atmosphere_step; rng=rng)
    render_atmosphere!(pupil, atmosphere_renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(cmd, recon, slopes(wfs))
    fill!(dm.state.coefs, 0)
    @views @. dm.state.coefs[1:n_half] = -cmd[1:n_half]
end

@info "Closed-loop first-stage run complete"
