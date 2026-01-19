using AdaptiveOptics
using BenchmarkTools

function bench_psf()
    tel = Telescope(resolution=64, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    ws = Workspace(128; T=Float64)
    return @benchmark compute_psf!($tel, $src; zero_padding=2, ws=$ws)
end

function bench_wfs()
    tel = Telescope(resolution=64, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    wfs = ShackHartmann(tel; n_subap=8)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    return @benchmark measure!($wfs, $tel)
end

function bench_reconstructor()
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4)
    wfs = ShackHartmann(tel; n_subap=4)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = wfs.state.slopes
    return @benchmark reconstruct($recon, $slopes)
end

println("PSF benchmark:")
display(bench_psf())

println("WFS benchmark:")
display(bench_wfs())

println("Reconstructor benchmark:")
display(bench_reconstructor())
