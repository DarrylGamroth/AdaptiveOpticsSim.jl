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
sim = AOSimulation(tel, src, atm, dm, wfs)
atmosphere_renderer = prepare_atmosphere_renderer(sim.atm, sim.tel, sim.src)
atmosphere_output = PupilFunction(sim.tel)

screens = Vector{Matrix{Float64}}(undef, 5)
for k in 1:5
    epoch = advance_by!(sim.atm, atmosphere_step; rng=rng)
    render_atmosphere!(atmosphere_output, atmosphere_renderer, sim.atm, epoch)
    screens[k] = copy(atmosphere_output.opd)
end

imat = interaction_matrix(sim.optic, sim.wfs, sim.tel; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
cmd = similar(sim.optic.state.coefs)

for k in 1:length(screens)
    copyto!(sim.tel.state.opd, screens[k])
    apply!(sim.optic, sim.tel, DMAdditive())
    measure!(sim.wfs, sim.tel)
    reconstruct!(cmd, recon, slopes(sim.wfs))
    sim.optic.state.coefs .= -cmd
end

@info "Closed-loop from phase screens complete"
