using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

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

    record_gpu_smoke!(failures, "atmosphere_multilayer_step") do
        atm = MultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            T=T,
            backend=BackendArray,
        )
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        @assert atm.state.opd isa BackendArray
        @assert tel.state.opd isa BackendArray
        return atm.state.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_multilayer_step") do
        atm = InfiniteMultiLayerAtmosphere(tel;
            r0=T(0.2),
            L0=T(25.0),
            fractional_cn2=T[0.7, 0.3],
            wind_speed=T[8.0, 4.0],
            wind_direction=T[0.0, 90.0],
            altitude=T[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
            T=T,
            backend=BackendArray,
        )
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel, src)
        @assert atm.state.opd isa BackendArray
        @assert tel.state.opd isa BackendArray
        @assert atm.layers[1].screen.state.screen isa BackendArray
        @assert atm.layers[1].screen.state.screen_scratch isa BackendArray
        return atm.state.opd
    end

    record_gpu_smoke!(failures, "atmosphere_infinite_statistical_agreement") do
        function trajectory_stats(backend; steps::Int=10)
            local_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
                central_obstruction=0.0f0, T=T, backend=backend)
            local_src = Source(band=:I, magnitude=0.0, coordinates=(30.0, 20.0), T=T)
            local_atm = InfiniteMultiLayerAtmosphere(local_tel;
                r0=T(0.2),
                L0=T(25.0),
                fractional_cn2=T[0.7, 0.3],
                wind_speed=T[8.0, 4.0],
                wind_direction=T[0.0, 90.0],
                altitude=T[0.0, 5000.0],
                screen_resolution=33,
                stencil_size=35,
                T=T,
                backend=backend,
            )
            local_rng = MersenneTwister(21)
            stds = Float64[]
            corrs = Float64[]
            previous = nothing
            for _ in 1:steps
                advance!(local_atm, local_tel; rng=local_rng)
                propagate!(local_atm, local_tel, local_src)
                opd = Array(local_tel.state.opd)
                push!(stds, std(vec(opd)))
                if previous !== nothing
                    push!(corrs, dot(vec(previous), vec(opd)) / sqrt(dot(vec(previous), vec(previous)) * dot(vec(opd), vec(opd))))
                end
                previous = opd
            end
            return (; std_mean=mean(stds), corr_mean=mean(corrs))
        end

        cpu_stats = trajectory_stats(Array)
        gpu_stats = trajectory_stats(BackendArray)
        std_rel = abs(gpu_stats.std_mean - cpu_stats.std_mean) / max(cpu_stats.std_mean, eps(Float64))
        corr_abs = abs(gpu_stats.corr_mean - cpu_stats.corr_mean)
        @assert std_rel < 0.5
        @assert corr_abs < 0.05
        return gpu_stats
    end

    record_gpu_smoke!(failures, "atmosphere_phase_helpers") do
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, T=T, backend=BackendArray)
        ws = PhaseStatsWorkspace(tel.params.resolution; T=T, backend=BackendArray)
        screen, psd = ft_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, return_psd=true)
        sh_screen = ft_sh_phase_screen(atm, tel.params.resolution, tel.params.diameter / tel.params.resolution;
            rng=rng, ws=ws, subharmonics=true, n_levels=2, subharmonic_radius=1)
        rho = similar(atm.state.opd, T, 4, 4)
        copyto!(rho, reshape(T[0.0, 0.02, 0.05, 0.1, 0.15, 0.25, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0], 4, 4))
        cov = phase_covariance(rho, atm)
        freqs = similar(atm.state.freqs, T, 4)
        copyto!(freqs, T[0.1, 0.2, 0.3, 0.4])
        covmat = covariance_matrix(freqs, freqs, atm)
        spectrum = phase_spectrum(freqs, atm)
        freq_grid = similar(atm.state.opd, T, 2, 2)
        copyto!(freq_grid, reshape(T[0.1, 0.2, 0.3, 0.4], 2, 2))
        spectrum_grid = phase_spectrum(freq_grid, atm)
        @assert screen isa BackendArray
        @assert psd isa BackendArray
        @assert sh_screen isa BackendArray
        @assert cov isa BackendArray
        @assert covmat isa BackendArray
        @assert spectrum isa BackendArray
        @assert spectrum_grid isa BackendArray
        cov_ref = phase_covariance(Array(rho), atm)
        covmat_ref = covariance_matrix(Array(freqs), Array(freqs), atm)
        cov_rel = maximum(abs.(Array(cov) .- cov_ref)) / maximum(abs.(cov_ref))
        covmat_rel = maximum(abs.(Array(covmat) .- covmat_ref)) / maximum(abs.(covmat_ref))
        @assert cov_rel < 5e-4
        @assert covmat_rel < 5e-4
        return screen
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
