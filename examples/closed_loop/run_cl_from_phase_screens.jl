using AdaptiveOpticsSim
using Random
using Logging

rng = MersenneTwister(2)
atmosphere_step = 1e-3

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[8.0], wind_direction=[45.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)
atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)

screens = Vector{Matrix{Float64}}(undef, 5)
for k in 1:5
    epoch = advance_by!(atm, atmosphere_step; rng=rng)
    render_atmosphere!(pupil, atmosphere_renderer, atm, epoch)
    screens[k] = copy(pupil.opd)
end

imat = interaction_matrix(dm, wfs, pupil, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
cmd = similar(dm.state.coefs)

for k in 1:length(screens)
    copyto!(pupil.opd, screens[k])
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(cmd, recon, slopes(wfs))
    @. dm.state.coefs = -cmd
end

@info "Closed-loop from phase screens complete"
