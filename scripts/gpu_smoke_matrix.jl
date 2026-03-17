using AdaptiveOpticsSim
using Random

try
    using CUDA
catch err
    error("gpu_smoke_matrix.jl requires CUDA.jl: $(sprint(showerror, err))")
end

CUDA.allowscalar(false)

function record!(f::Function, failures::Vector{String}, name::AbstractString)
    try
        result = f()
        println(name, " ok :: ", typeof(result))
    catch err
        push!(failures, string(name, " ERROR :: ", sprint(showerror, err)))
        println(name, " ERROR :: ", sprint(showerror, err))
    end
    return nothing
end

function main()
    failures = String[]
    rng = MersenneTwister(1)
    T = Float32

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=CuArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    lgs = LGSSource(; magnitude=0.0, wavelength=589e-9, altitude=90_000.0,
        laser_coordinates=(0.0, 0.0), T=T)

    record!(failures, "psf_source") do
        psf = compute_psf!(tel, src; zero_padding=2)
        @assert psf isa CuArray
        return psf
    end

    record!(failures, "psf_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0)),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0)),
        ])
        psf = compute_psf!(tel, ast; zero_padding=2)
        @assert psf isa CuArray
        return psf
    end

    record!(failures, "detector_capture_none") do
        psf = compute_psf!(tel, src; zero_padding=2)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=2, T=T, backend=CuArray)
        frame = capture!(det, psf; rng=rng)
        @assert frame isa CuArray
        return frame
    end

    record!(failures, "detector_capture_noise") do
        psf = compute_psf!(tel, src; zero_padding=2)
        det = Detector(noise=(NoisePhoton(), NoiseReadout(T(1e-3))), integration_time=1.0, qe=1.0,
            binning=2, background_flux=T(0.5), dark_current=T(0.1), T=T, backend=CuArray)
        frame = capture!(det, psf; rng=rng)
        @assert frame isa CuArray
        return frame
    end

    record!(failures, "atmosphere_step") do
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=CuArray)
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        @assert tel.state.opd isa CuArray
        return tel.state.opd
    end

    record!(failures, "dm_apply") do
        dm = DeformableMirror(tel; n_act=4, influence_width=0.3, T=T, backend=CuArray)
        fill!(dm.state.coefs, T(0.05))
        apply!(dm, tel, DMReplace())
        @assert tel.state.opd isa CuArray
        return tel.state.opd
    end

    record!(failures, "measure_shack_geometric") do
        wfs = ShackHartmann(tel; n_subap=4, T=T, backend=CuArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_shack_diffractive") do
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_shack_lgs") do
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        slopes = measure!(wfs, tel, lgs)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_shack_diffractive_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        slopes = measure!(wfs, tel, ast)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_shack_diffractive_asterism_detector") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=CuArray)
        slopes = measure!(wfs, tel, ast, det; rng=rng)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_pyramid_geometric") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, T=T, backend=CuArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_pyramid_diffractive") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=CuArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_pyramid_detector") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=CuArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=CuArray)
        slopes = measure!(wfs, tel, src, det; rng=rng)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_bioedge_geometric") do
        wfs = BioEdgeWFS(tel; n_subap=4, modulation=0.0, T=T, backend=CuArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "measure_bioedge_diffractive") do
        wfs = BioEdgeWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=CuArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa CuArray
        return slopes
    end

    record!(failures, "closed_loop_step") do
        step_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=CuArray)
        atm = KolmogorovAtmosphere(step_tel; r0=0.2, L0=25.0, T=T, backend=CuArray)
        dm = DeformableMirror(step_tel; n_act=4, influence_width=0.3, T=T, backend=CuArray)
        wfs = ShackHartmann(step_tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=CuArray)
        advance!(atm, step_tel; rng=rng)
        propagate!(atm, step_tel)
        apply!(dm, step_tel, DMAdditive())
        slopes = measure!(wfs, step_tel, src, det; rng=rng)
        psf = compute_psf!(step_tel, src; zero_padding=2)
        frame = capture!(det, psf; rng=rng)
        @assert slopes isa CuArray
        @assert frame isa CuArray
        return frame
    end

    record!(failures, "closed_loop_runtime") do
        rt_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=CuArray)
        rt_src = Source(band=:I, magnitude=0.0, T=T)
        rt_atm = KolmogorovAtmosphere(rt_tel; r0=0.2, L0=25.0, T=T, backend=CuArray)
        rt_dm = DeformableMirror(rt_tel; n_act=4, influence_width=0.3, T=T, backend=CuArray)
        rt_wfs = ShackHartmann(rt_tel; n_subap=4, T=T, backend=CuArray)
        rt_sim = AOSimulation(rt_tel, rt_atm, rt_src, rt_dm, rt_wfs)
        rt_imat = interaction_matrix(rt_dm, rt_wfs, rt_tel; amplitude=T(0.05))
        rt_recon = ModalReconstructor(rt_imat; gain=T(0.5))
        runtime = ClosedLoopRuntime(rt_sim, rt_recon; rng=rng)
        step!(runtime)
        @assert runtime.command isa CuArray
        @assert rt_dm.state.coefs isa CuArray
        return runtime.command
    end

    record!(failures, "closed_loop_runtime_science") do
        rt_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=CuArray)
        rt_src = Source(band=:I, magnitude=0.0, T=T)
        rt_atm = KolmogorovAtmosphere(rt_tel; r0=0.2, L0=25.0, T=T, backend=CuArray)
        rt_dm = DeformableMirror(rt_tel; n_act=4, influence_width=0.3, T=T, backend=CuArray)
        rt_wfs = ShackHartmann(rt_tel; n_subap=4, T=T, backend=CuArray)
        rt_det = Detector(NoiseNone(); psf_sampling=2, T=T, backend=CuArray)
        rt_sim = AOSimulation(rt_tel, rt_atm, rt_src, rt_dm, rt_wfs)
        rt_imat = interaction_matrix(rt_dm, rt_wfs, rt_tel; amplitude=T(0.05))
        rt_recon = ModalReconstructor(rt_imat; gain=T(0.5))
        runtime = ClosedLoopRuntime(rt_sim, rt_recon; rng=rng, science_detector=rt_det)
        step!(runtime)
        frame = output_frame(rt_det)
        @assert frame isa CuArray
        return frame
    end

    record!(failures, "interaction_matrix_reconstructor") do
        cal_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=CuArray)
        dm = DeformableMirror(cal_tel; n_act=4, influence_width=0.3, T=T, backend=CuArray)
        wfs = ShackHartmann(cal_tel; n_subap=4, mode=Diffractive(), T=T, backend=CuArray)
        imat = interaction_matrix(dm, wfs, cal_tel, src; amplitude=T(0.05))
        recon = ModalReconstructor(imat; gain=one(T))
        measure!(wfs, cal_tel, src)
        cmds = reconstruct(recon, wfs.state.slopes)
        @assert imat.matrix isa CuArray
        @assert cmds isa CuArray
        return cmds
    end

    record!(failures, "gain_sensing_camera") do
        mask = CUDA.fill(one(T), 8, 8)
        basis = CUDA.rand(T, 8, 8, 3)
        frame = abs.(CUDA.randn(T, 8, 8))
        gsc = GainSensingCamera(mask, basis; T=T)
        calibrate!(gsc, frame)
        og = compute_optical_gains!(gsc, frame)
        @assert og isa CuArray
        return og
    end

    record!(failures, "lift") do
        lift_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=CuArray)
        lift_src = Source(band=:I, magnitude=8.0, T=T)
        basis = CUDA.rand(T, 16, 16, 3)
        diversity = CUDA.zeros(T, 16, 16)
        det = Detector(NoiseNone(); psf_sampling=2, T=T, backend=CuArray)
        lift = LiFT(lift_tel, lift_src, basis, det; diversity_opd=diversity,
            iterations=2, img_resolution=32, solve_mode=LiFTSolveAuto())
        psf = compute_psf!(lift_tel, lift_src; zero_padding=2)
        coeffs = reconstruct(lift, psf, [1, 2])
        @assert coeffs isa CuArray
        return coeffs
    end

    if !isempty(failures)
        error("GPU smoke matrix failed:\n" * join(failures, "\n"))
    end

    println("gpu_smoke_matrix complete")
    return nothing
end

main()
