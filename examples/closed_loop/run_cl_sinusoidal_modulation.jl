using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(4)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[6.0], wind_direction=[30.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0)
atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)

imat = interaction_matrix(dm, wfs, pupil, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(dm.state.coefs)

n_iter = 6
for k in 1:n_iter
    epoch = advance_by!(atm, atmosphere_step; rng=rng)
    render_atmosphere!(pupil, atmosphere_renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(cmd, recon, slopes(wfs))
    modulation = 1 + 0.2 * sin(2 * pi * k / n_iter)
    @. dm.state.coefs = -cmd * modulation
end

@info "Closed-loop sinusoidal modulation complete"
