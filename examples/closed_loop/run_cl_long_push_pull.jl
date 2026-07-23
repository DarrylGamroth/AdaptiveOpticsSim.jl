using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Logging

rng = MersenneTwister(3)

tel = Telescope(resolution=32, diameter=8.0)
src = Source()
atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
    wind_speed=[5.0], wind_direction=[0.0], altitude=[0.0])
dm = DeformableMirror(tel; n_act=4, influence_width=0.2)
wfs = ShackHartmannWFS(tel; n_lenslets=4)

calibration_pupil = PupilFunction(tel)
imat1 = interaction_matrix(dm, wfs, calibration_pupil, src;
    amplitude=0.05)
imat2 = interaction_matrix(dm, wfs, calibration_pupil, src;
    amplitude=0.1)
mat = 0.5 .* (imat1.matrix .+ imat2.matrix)
recon = ModalReconstructor(InteractionMatrix(mat, 0.1); gain=0.4)
pupil = PupilFunction(tel)
renderer = prepare_atmosphere_renderer(atm, tel, src)
command = similar(dm.state.coefs)

for _ in 1:5
    epoch = advance_by!(atm, 1e-3; rng=rng)
    render_atmosphere!(pupil, renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src)
    reconstruct!(command, recon, slopes(wfs))
    @. command = -command
    set_command!(dm, command)
end

@info "Closed-loop long push-pull complete" command_norm=norm(command)
