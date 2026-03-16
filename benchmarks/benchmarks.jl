using AdaptiveOpticsSim
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

function bench_wfs_lgs()
    tel = Telescope(resolution=48, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    wfs = ShackHartmann(tel; n_subap=6, mode=Diffractive())
    lgs = LGSSource(elongation_factor=1.3)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i - j
    end
    return @benchmark measure!($wfs, $tel, $lgs)
end

function bench_pyramid()
    tel = Telescope(resolution=48, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    wfs = PyramidWFS(tel; n_subap=6, modulation=3.0, modulation_points=4, mode=Diffractive())
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i + j
    end
    return @benchmark measure!($wfs, $tel, $src)
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

function bench_reconstructor_inplace()
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4)
    wfs = ShackHartmann(tel; n_subap=4)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = wfs.state.slopes
    out = similar(slopes, size(recon.reconstructor, 1))
    return @benchmark reconstruct!($out, $recon, $slopes)
end

function bench_lift(numerical::Bool)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(16, 16, 6)
    diversity = zeros(16, 16)
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=1,
        img_resolution=16, numerical=numerical)
    coeffs = zeros(6)
    return @benchmark lift_interaction_matrix($lift, $coeffs, 1:3)
end

function bench_lift_inplace(numerical::Bool)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(16, 16, 6)
    diversity = zeros(16, 16)
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=1,
        img_resolution=16, numerical=numerical)
    coeffs = zeros(6)
    mode_ids = 1:3
    H = @view lift.state.H_buffer[:, 1:length(mode_ids)]
    return @benchmark lift_interaction_matrix!($H, $lift, $coeffs, $mode_ids)
end

function alloc_checks()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(16, 16, 6)
    diversity = zeros(16, 16)
    lift_a = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=1,
        img_resolution=16, numerical=false)
    lift_n = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=1,
        img_resolution=16, numerical=true)
    coeffs = zeros(6)
    mode_ids = 1:3
    H_a = @view lift_a.state.H_buffer[:, 1:length(mode_ids)]
    H_n = @view lift_n.state.H_buffer[:, 1:length(mode_ids)]
    lift_interaction_matrix!(H_a, lift_a, coeffs, mode_ids)
    lift_interaction_matrix!(H_n, lift_n, coeffs, mode_ids)
    alloc_lift_a = @allocated lift_interaction_matrix!(H_a, lift_a, coeffs, mode_ids)
    alloc_lift_n = @allocated lift_interaction_matrix!(H_n, lift_n, coeffs, mode_ids)

    tel_r = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel_r; n_act=4)
    wfs = ShackHartmann(tel_r; n_subap=4)
    imat = interaction_matrix(dm, wfs, tel_r; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = wfs.state.slopes
    out = similar(slopes, size(recon.reconstructor, 1))
    reconstruct!(out, recon, slopes)
    alloc_recon = @allocated reconstruct!(out, recon, slopes)

    println("Allocation checks:")
    println("  LiFT analytic (in-place): $(alloc_lift_a) bytes")
    println("  LiFT numerical (in-place): $(alloc_lift_n) bytes")
    println("  Reconstructor! (in-place): $(alloc_recon) bytes")
end

println("PSF benchmark:")
display(bench_psf())

println("WFS benchmark:")
display(bench_wfs())

println("WFS LGS benchmark:")
display(bench_wfs_lgs())

println("Pyramid benchmark:")
display(bench_pyramid())

println("Reconstructor benchmark:")
display(bench_reconstructor())

println("Reconstructor in-place benchmark:")
display(bench_reconstructor_inplace())

println("LiFT analytic benchmark:")
display(bench_lift(false))

println("LiFT numerical benchmark:")
display(bench_lift(true))

println("LiFT analytic in-place benchmark:")
display(bench_lift_inplace(false))

println("LiFT numerical in-place benchmark:")
display(bench_lift_inplace(true))

alloc_checks()
