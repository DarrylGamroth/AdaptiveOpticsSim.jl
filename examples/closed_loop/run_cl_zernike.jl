using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Logging

rng = MersenneTwister(4)

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
src = Source(band=:I, magnitude=0.0)
atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
wfs = ZernikeWFS(tel; pupil_samples=8, diffraction_padding=2)
det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

imat = interaction_matrix(dm, wfs, PupilFunction(tel), src; amplitude=1e-8)
recon = ModalReconstructor(imat; gain=0.5)
pupil = PupilFunction(tel)
renderer = prepare_atmosphere_renderer(atm, tel, src)
prepare_runtime_wfs!(wfs, pupil, src)
command = similar(dm.state.coefs)
for _ in 1:5
    epoch = advance_by!(atm, 1e-3; rng=rng)
    render_atmosphere!(pupil, renderer, atm, epoch)
    update_surface!(dm)
    apply_surface!(pupil, dm, DMAdditive())
    measure!(wfs, pupil, src, det; rng=rng)
    reconstruct!(command, recon, slopes(wfs))
    @. command = -command
    set_command!(dm, command)
end

@info "Closed-loop ZernikeWFS tutorial complete" command_norm=norm(command) slope_norm=norm(slopes(wfs))
