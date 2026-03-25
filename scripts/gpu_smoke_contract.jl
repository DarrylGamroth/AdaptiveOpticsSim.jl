using AdaptiveOpticsSim
using Random

function record_gpu_smoke!(f::Function, failures::Vector{String}, name::AbstractString)
    try
        result = f()
        println(name, " ok :: ", typeof(result))
    catch err
        push!(failures, string(name, " ERROR :: ", sprint(showerror, err)))
        println(name, " ERROR :: ", sprint(showerror, err))
    end
    return nothing
end

function run_gpu_smoke_matrix(::Type{B}) where {B<:GPUBackendTag}
    disable_scalar_backend!(B)
    failures = String[]
    rng = MersenneTwister(1)
    T = Float32
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    lgs = LGSSource(; magnitude=0.0, wavelength=589e-9, altitude=90_000.0,
        laser_coordinates=(0.0, 0.0), T=T)

    record_gpu_smoke!(failures, "psf_source") do
        psf = compute_psf!(tel, src; zero_padding=2)
        @assert psf isa BackendArray
        return psf
    end

    record_gpu_smoke!(failures, "psf_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0)),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0)),
        ])
        psf = compute_psf!(tel, ast; zero_padding=2)
        @assert psf isa BackendArray
        return psf
    end

    record_gpu_smoke!(failures, "detector_capture_none") do
        psf = compute_psf!(tel, src; zero_padding=2)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=2, T=T, backend=BackendArray)
        frame = capture!(det, psf; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "detector_capture_noise") do
        psf = compute_psf!(tel, src; zero_padding=2)
        det = Detector(noise=(NoisePhoton(), NoiseReadout(T(1e-3))), integration_time=1.0, qe=1.0,
            binning=2, background_flux=T(0.5), dark_current=T(0.1), T=T, backend=BackendArray)
        frame = capture!(det, psf; rng=rng)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "atmosphere_step") do
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        @assert tel.state.opd isa BackendArray
        return tel.state.opd
    end

    record_gpu_smoke!(failures, "dm_apply") do
        dm = DeformableMirror(tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        fill!(dm.state.coefs, T(0.05))
        apply!(dm, tel, DMReplace())
        @assert tel.state.opd isa BackendArray
        return tel.state.opd
    end

    record_gpu_smoke!(failures, "measure_shack_geometric") do
        wfs = ShackHartmann(tel; n_subap=4, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive") do
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_lgs") do
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, lgs)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_asterism") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, ast)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_shack_diffractive_asterism_detector") do
        ast = Asterism([
            Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0), T=T),
            Source(band=:I, magnitude=0.0, coordinates=(1.0, 45.0), T=T),
        ])
        wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, ast, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_geometric") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_diffractive") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_pyramid_detector") do
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=BackendArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_bioedge_geometric") do
        wfs = BioEdgeWFS(tel; n_subap=4, modulation=0.0, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_bioedge_diffractive") do
        wfs = BioEdgeWFS(tel; n_subap=4, modulation=2.0, mode=Diffractive(), T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_zernike_diffractive") do
        wfs = ZernikeWFS(tel; n_subap=4, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "measure_zernike_detector") do
        wfs = ZernikeWFS(tel; n_subap=4, T=T, backend=BackendArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=BackendArray)
        slopes = measure!(wfs, tel, src, det; rng=rng)
        @assert slopes isa BackendArray
        return slopes
    end

    record_gpu_smoke!(failures, "closed_loop_step") do
        step_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        atm = KolmogorovAtmosphere(step_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        dm = DeformableMirror(step_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        wfs = ShackHartmann(step_tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1, T=T, backend=BackendArray)
        advance!(atm, step_tel; rng=rng)
        propagate!(atm, step_tel)
        apply!(dm, step_tel, DMAdditive())
        slopes = measure!(wfs, step_tel, src, det; rng=rng)
        psf = compute_psf!(step_tel, src; zero_padding=2)
        frame = capture!(det, psf; rng=rng)
        @assert slopes isa BackendArray
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "closed_loop_runtime") do
        rt_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        rt_src = Source(band=:I, magnitude=0.0, T=T)
        rt_atm = KolmogorovAtmosphere(rt_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        rt_dm = DeformableMirror(rt_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        rt_wfs = ShackHartmann(rt_tel; n_subap=4, T=T, backend=BackendArray)
        rt_sim = AOSimulation(rt_tel, rt_atm, rt_src, rt_dm, rt_wfs)
        rt_imat = interaction_matrix(rt_dm, rt_wfs, rt_tel; amplitude=T(0.05))
        rt_recon = ModalReconstructor(rt_imat; gain=T(0.5))
        runtime = ClosedLoopRuntime(rt_sim, rt_recon; rng=rng)
        step!(runtime)
        @assert runtime.command isa BackendArray
        @assert rt_dm.state.coefs isa BackendArray
        return runtime.command
    end

    record_gpu_smoke!(failures, "closed_loop_runtime_science") do
        rt_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        rt_src = Source(band=:I, magnitude=0.0, T=T)
        rt_atm = KolmogorovAtmosphere(rt_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        rt_dm = DeformableMirror(rt_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        rt_wfs = ShackHartmann(rt_tel; n_subap=4, T=T, backend=BackendArray)
        rt_det = Detector(NoiseNone(); psf_sampling=2, T=T, backend=BackendArray)
        rt_sim = AOSimulation(rt_tel, rt_atm, rt_src, rt_dm, rt_wfs)
        rt_imat = interaction_matrix(rt_dm, rt_wfs, rt_tel; amplitude=T(0.05))
        rt_recon = ModalReconstructor(rt_imat; gain=T(0.5))
        runtime = ClosedLoopRuntime(rt_sim, rt_recon; rng=rng, science_detector=rt_det)
        step!(runtime)
        frame = output_frame(rt_det)
        @assert frame isa BackendArray
        return frame
    end

    record_gpu_smoke!(failures, "closed_loop_multi_boundary") do
        rt1_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        rt1_src = Source(band=:I, magnitude=0.0, T=T)
        rt1_atm = KolmogorovAtmosphere(rt1_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        rt1_dm = DeformableMirror(rt1_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        rt1_wfs = ShackHartmann(rt1_tel; n_subap=4, T=T, backend=BackendArray)
        rt1_sim = AOSimulation(rt1_tel, rt1_atm, rt1_src, rt1_dm, rt1_wfs)
        rt1_recon = ModalReconstructor(interaction_matrix(rt1_dm, rt1_wfs, rt1_tel; amplitude=T(0.05)); gain=T(0.5))
        rt1 = ClosedLoopRuntime(rt1_sim, rt1_recon; rng=rng)

        rt2_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        rt2_src = Source(band=:I, magnitude=0.0, T=T)
        rt2_atm = KolmogorovAtmosphere(rt2_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        rt2_dm = DeformableMirror(rt2_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        rt2_wfs = PyramidWFS(rt2_tel; n_subap=4, mode=Geometric(), T=T, backend=BackendArray)
        rt2_sim = AOSimulation(rt2_tel, rt2_atm, rt2_src, rt2_dm, rt2_wfs)
        rt2_recon = ModalReconstructor(interaction_matrix(rt2_dm, rt2_wfs, rt2_tel; amplitude=T(0.05)); gain=T(0.5))
        rt2 = ClosedLoopRuntime(rt2_sim, rt2_recon; rng=rng)

        boundary = CompositeSimulationInterface(rt1, rt2)
        step!(boundary)
        @assert simulation_command(boundary) isa BackendArray
        @assert simulation_slopes(boundary) isa BackendArray
        return simulation_command(boundary)
    end

    record_gpu_smoke!(failures, "runtime_reconstructor_refresh") do
        rt_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        rt_src = Source(band=:I, magnitude=0.0, T=T)
        rt_atm = KolmogorovAtmosphere(rt_tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        rt_dm = DeformableMirror(rt_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        rt_wfs = ShackHartmann(rt_tel; n_subap=4, T=T, backend=BackendArray)
        rt_sim = AOSimulation(rt_tel, rt_atm, rt_src, rt_dm, rt_wfs)
        rt_imat = interaction_matrix(rt_dm, rt_wfs, rt_tel; amplitude=T(0.05))
        runtime = ClosedLoopRuntime(rt_sim, ModalReconstructor(rt_imat; gain=T(0.5)); rng=rng)
        refreshed = with_reconstructor(runtime, ModalReconstructor(rt_imat; gain=T(0.25)))
        step!(refreshed)
        @assert refreshed.command isa BackendArray
        return refreshed.command
    end

    record_gpu_smoke!(failures, "interaction_matrix_reconstructor") do
        cal_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        dm = DeformableMirror(cal_tel; n_act=4, influence_width=0.3, T=T, backend=BackendArray)
        wfs = ShackHartmann(cal_tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
        imat = interaction_matrix(dm, wfs, cal_tel, src; amplitude=T(0.05))
        recon = ModalReconstructor(imat; gain=one(T))
        measure!(wfs, cal_tel, src)
        cmds = reconstruct(recon, wfs.state.slopes)
        @assert imat.matrix isa BackendArray
        @assert cmds isa BackendArray
        return cmds
    end

    record_gpu_smoke!(failures, "gain_sensing_camera") do
        mask = backend_fill(B, one(T), 8, 8)
        basis = backend_rand(B, T, 8, 8, 3)
        frame = abs.(backend_randn(B, T, 8, 8))
        gsc = GainSensingCamera(mask, basis; T=T)
        calibrate!(gsc, frame)
        og = compute_optical_gains!(gsc, frame)
        @assert og isa BackendArray
        return og
    end

    record_gpu_smoke!(failures, "lift") do
        lift_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
            central_obstruction=0.0f0, T=T, backend=BackendArray)
        lift_src = Source(band=:I, magnitude=8.0, T=T)
        basis = backend_rand(B, T, 16, 16, 3)
        diversity = backend_zeros(B, T, 16, 16)
        det = Detector(NoiseNone(); psf_sampling=2, T=T, backend=BackendArray)
        lift = LiFT(lift_tel, lift_src, basis, det; diversity_opd=diversity,
            iterations=2, img_resolution=32, solve_mode=LiFTSolveAuto())
        psf = compute_psf!(lift_tel, lift_src; zero_padding=2)
        coeffs = reconstruct(lift, psf, [1, 2])
        @assert coeffs isa BackendArray
        return coeffs
    end

    if !isempty(failures)
        error("GPU smoke matrix failed:\n" * join(failures, "\n"))
    end

    println("gpu_smoke_matrix complete")
    return nothing
end
