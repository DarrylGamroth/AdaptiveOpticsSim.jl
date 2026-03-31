@testset "Deformable mirror and WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    dm.state.coefs .= 0.1
    apply!(dm, tel, DMReplace())
    @test sum(abs.(tel.state.opd)) > 0

    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    wfs = ShackHartmann(tel; n_subap=4)
    slopes = measure!(wfs, tel)
    @test length(slopes) == 2 * 4 * 4
    @test maximum(slopes) > 0
end

@testset "Calibration and control" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    @test size(imat.matrix) == (length(wfs.state.slopes), length(dm.state.coefs))

    recon = ModalReconstructor(imat; gain=1.0)
    cmd = reconstruct(recon, wfs.state.slopes)
    @test length(cmd) == length(dm.state.coefs)

    mapped = MappedReconstructor(Matrix{Float64}(I, length(dm.state.coefs), length(dm.state.coefs)), imat; gain=0.5)
    mapped_cmd = reconstruct(mapped, wfs.state.slopes)
    @test length(mapped_cmd) == length(dm.state.coefs)
    @test mapped.n_control_modes == length(dm.state.coefs)

    delay = VectorDelayLine(dm.state.coefs, 1)
    delayed0 = shift_delay!(delay, fill(1.0, length(dm.state.coefs)))
    @test all(iszero, delayed0)
    delayed1 = shift_delay!(delay, fill(2.0, length(dm.state.coefs)))
    @test all(==(1.0), delayed1)

    ctrl = DiscreteIntegratorController(length(wfs.state.slopes); gain=0.1, tau=0.02)
    dm_cmd = update!(ctrl, wfs.state.slopes, 0.01)
    @test length(dm_cmd) == length(wfs.state.slopes)

    vault = CalibrationVault(imat.matrix; build_backend=CPUBuildBackend())
    @test vault.M isa Matrix
    recon_cpu = ModalReconstructor(imat; build_backend=CPUBuildBackend())
    @test recon_cpu.reconstructor isa Matrix
end

@testset "AO188/3k simulation" begin
    default_params = AO188SimulationParams()
    @test default_params.n_act == 64
    @test default_params.n_active_actuators == 3228
    @test default_params.n_control_modes == 188
    @test default_params.n_low_order_subap == 2
    @test default_params.low_order_resolution == 28
    @test default_params.n_low_order_modes == 4
    @test default_params.latency.high_measurement_delay_frames == 1
    @test default_params.high_detector.noise isa NoisePhotonReadout
    @test default_params.branch_execution isa SequentialExecution
    @test default_params.replay_mode isa DirectReplayMode

    params = AO188SimulationParams(
        T=Float32,
        resolution=48,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=24,
        control_grid_side=6,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    @test params.low_order_resolution == 12
    surrogate = subaru_ao188_simulation(; params=params, rng=MersenneTwister(1))
    @test count(surrogate.active_mask) == params.n_active_actuators
    @test length(surrogate.active_indices) == params.n_active_actuators
    @test !isnothing(surrogate.dm.state.separable_x)
    @test !isnothing(surrogate.dm.state.separable_y_t)
    @test size(surrogate.high_M2C) == (params.n_act^2, params.n_control_modes)
    @test size(surrogate.low_M2C) == (params.n_act^2, params.n_low_order_modes)
    @test size(surrogate.high_reconstructor.command_basis, 1) == params.n_act^2
    @test size(surrogate.high_reconstructor.command_basis, 2) == params.n_control_modes
    @test size(surrogate.high_reconstructor.reconstructor, 1) == params.n_control_modes
    @test size(surrogate.high_reconstructor.reconstructor, 2) == length(surrogate.high_wfs.state.slopes)
    @test size(surrogate.low_reconstructor.command_basis, 2) == params.n_low_order_modes
    @test size(surrogate.low_reconstructor.reconstructor, 2) == length(surrogate.low_wfs.state.slopes)
    @test surrogate.high_reconstructor.n_control_modes == params.n_control_modes
    @test surrogate.low_reconstructor.n_control_modes == params.n_low_order_modes

    shifted = DeformableMirror(Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, T=Float32);
        n_act=8, T=Float32, misregistration=Misregistration(shift_x=0.01, shift_y=-0.02, T=Float32))
    @test !isnothing(shifted.state.separable_x)

    rotated = DeformableMirror(Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, T=Float32);
        n_act=8, T=Float32, misregistration=Misregistration(rotation_deg=1.0, T=Float32))
    @test isnothing(rotated.state.separable_x)

    step!(surrogate)
    step!(surrogate)
    step!(surrogate)
    @test length(surrogate.command) == params.n_act^2
    @test maximum(abs, surrogate.command) > 0
    @test supports_prepared_runtime(typeof(surrogate))
    @test supports_detector_output(typeof(surrogate))
    @test supports_grouped_execution(typeof(surrogate))
    simif = simulation_interface(surrogate)
    @test simulation_command(simif) === surrogate.command
    @test length(simulation_slopes(simif)) == 2
    timing = runtime_timing(surrogate; warmup=1, samples=2, gc_before=false)
    @test timing.samples == 2
    phase = subaru_ao188_phase_timing(surrogate; warmup=1, samples=2, gc_before=false)
    @test phase.samples == 2
    @test phase.delay_mean_ns >= 0

    experimental = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
        branch_execution=ThreadedExecution(),
        replay_mode=PreparedReplayMode(),
    )
    @test experimental.low_order_resolution == 8
    surrogate_exp = subaru_ao188_simulation(; params=experimental, rng=MersenneTwister(2))
    @test surrogate_exp.replay_prepared
    step!(surrogate_exp)
    @test maximum(abs, surrogate_exp.command) >= 0

    stream_mode = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
        branch_execution=BackendStreamExecution(),
    )
    surrogate_stream = subaru_ao188_simulation(; params=stream_mode, rng=MersenneTwister(3))
    step!(surrogate_stream)
    @test maximum(abs, surrogate_stream.command) >= 0

    curvature_params = AO188CurvatureSimulationParams(
        T=Float32,
        resolution=48,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=24,
        control_grid_side=6,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    curvature_sim = subaru_ao188_curvature_simulation(; params=curvature_params, rng=MersenneTwister(4))
    @test curvature_sim.high_wfs isa CurvatureWFS
    @test curvature_sim.high_wfs.params.readout_model isa CurvatureCountingReadout
    @test curvature_sim.high_wfs.params.readout_crop_resolution == 32
    @test curvature_sim.high_wfs.params.readout_pixels_per_subap == 1
    @test curvature_params.high_detector isa AO188APDDetectorConfig
    @test curvature_sim.high_detector isa APDDetector
    step!(curvature_sim)
    @test length(curvature_sim.command) == curvature_params.n_act^2
    curvature_readout = simulation_interface(curvature_sim)
    @test size(simulation_wfs_frame(curvature_readout)[1]) == (2, curvature_params.n_subap^2)
    @test simulation_wfs_metadata(curvature_readout)[1] isa CountingDetectorExportMetadata

    ao3k_params = AO3kSimulationParams(
        T=Float32,
        resolution=64,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=32,
        control_grid_side=6,
        n_subap=8,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    ao3k_sim = subaru_ao3k_simulation(; params=ao3k_params, rng=MersenneTwister(5))
    @test ao3k_sim.high_wfs isa PyramidWFS
    @test ao3k_sim.high_detector isa Detector
    @test ao3k_sim.high_detector.params.sensor isa HgCdTeAvalancheArraySensor
    @test ao3k_params.high_detector.thermal_model isa FixedTemperature
    step!(ao3k_sim)
    @test length(ao3k_sim.command) == ao3k_params.n_act^2
    ao3k_iface = simulation_interface(ao3k_sim)
    ao3k_high_meta = simulation_wfs_metadata(ao3k_iface)[1]
    @test ao3k_high_meta.sensor == :hgcdte_avalanche_array
    @test ao3k_high_meta.thermal_model == :fixed_temperature
    @test ao3k_high_meta.detector_temperature_K == 80.0f0
    @test ao3k_high_meta.provides_combined_frame
    @test ao3k_high_meta.provides_reference_frame
    @test ao3k_high_meta.provides_signal_frame
end

function closed_loop_runtime_allocations()
    rng = MersenneTwister(0)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmann(tel; n_subap=4)
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    runtime = ClosedLoopRuntime(sim, recon; rng=rng)
    step!(runtime)
    step!(runtime)
    return @allocated step!(runtime)
end

@testset "Closed-loop runtime" begin
    rng = MersenneTwister(0)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmann(tel; n_subap=4)
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime = ClosedLoopRuntime(sim, recon; rng=rng, science_detector=det)
    @test runtime isa AbstractControlSimulation
    @test !supports_prepared_runtime(runtime)
    @test supports_detector_output(runtime)
    @test runtime_profile(runtime) isa ScientificRuntimeProfile
    @test runtime_latency(runtime).measurement_delay_frames == 0
    @test runtime.science_zero_padding == 2
    @test !runtime.prepared

    step!(runtime)
    @test length(runtime.command) == length(dm.state.coefs)
    @test size(output_frame(det)) == (32, 32)
    @test closed_loop_runtime_allocations() == 0
    @test simulation_command(runtime) === runtime.command
    @test simulation_science_frame(runtime) === output_frame(det)

    boundary = SimulationInterface(runtime)
    readout = simulation_readout(boundary)
    @test length(simulation_slopes(boundary)) == length(wfs.state.slopes)
    @test simulation_slopes(readout) === simulation_slopes(boundary)
    @test size(simulation_science_frame(boundary)) == size(output_frame(det))
    science_metadata = simulation_science_metadata(boundary)
    @test science_metadata isa DetectorExportMetadata
    @test science_metadata.output_size == size(output_frame(det))
    @test science_metadata.frame_size == size(det.state.frame)
    step!(boundary)
    @test simulation_command(boundary) == runtime.command

    rng2 = MersenneTwister(2)
    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src2 = Source(band=:I, magnitude=0.0)
    atm2 = KolmogorovAtmosphere(tel2; r0=0.2, L0=25.0)
    dm2 = DeformableMirror(tel2; n_act=4, influence_width=0.3)
    wfs2 = ShackHartmann(tel2; n_subap=4, mode=Diffractive())
    sim2 = AdaptiveOpticsSim.AOSimulation(tel2, atm2, src2, dm2, wfs2)
    imat2 = interaction_matrix(dm2, wfs2, tel2, src2; amplitude=0.1)
    recon2 = ModalReconstructor(imat2; gain=0.5)
    wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime2 = ClosedLoopRuntime(sim2, recon2; rng=rng2, wfs_detector=wfs_det)
    @test supports_prepared_runtime(runtime2)
    prepare!(runtime2)
    @test runtime2.wfs.state.calibrated
    @test supports_detector_output(runtime2)
    step!(runtime2)
    @test simulation_wfs_frame(runtime2) === wfs2.state.exported_spot_cube
    @test simulation_wfs_frame(runtime2) !== wfs2.state.spot_cube
    @test wfs2.state.spot_cube !== wfs2.state.sampled_spot_cube
    @test simulation_wfs_frame(runtime2) !== wfs2.state.sampled_spot_cube
    @test size(simulation_wfs_frame(runtime2)) == size(wfs2.state.spot_cube)
    @test all(simulation_wfs_frame(runtime2) .>= wfs2.state.spot_cube)
    boundary2 = SimulationInterface(runtime2)
    @test ndims(simulation_wfs_frame(boundary2)) == 3
    @test size(simulation_wfs_frame(boundary2), 1) == wfs2.params.n_subap^2

    composite = CompositeSimulationInterface(boundary, boundary2)
    @test supports_grouped_execution(composite)
    @test length(simulation_command(composite)) == length(simulation_command(boundary)) + length(simulation_command(boundary2))
    @test length(simulation_slopes(composite)) == length(simulation_slopes(boundary)) + length(simulation_slopes(boundary2))
    @test length(simulation_wfs_frame(composite)) == 2
    step!(composite)
    @test simulation_command(composite) == vcat(simulation_command(boundary), simulation_command(boundary2))

    runtime2_slopes_only = ClosedLoopRuntime(
        sim2,
        recon2;
        rng=MersenneTwister(22),
        wfs_detector=wfs_det,
        products=RuntimeProductRequirements(slopes=true, wfs_pixels=false, science_pixels=false),
    )
    prepare!(runtime2_slopes_only)
    @test runtime_products(runtime2_slopes_only).wfs_pixels == false
    @test !supports_detector_output(runtime2_slopes_only)
    @test !runtime2_slopes_only.wfs.state.export_pixels_enabled
    step!(runtime2_slopes_only)
    @test simulation_wfs_frame(runtime2_slopes_only) === nothing
    slopes_only_boundary = SimulationInterface(runtime2_slopes_only)
    @test simulation_wfs_frame(slopes_only_boundary) === nothing

    timing = runtime_timing(runtime; warmup=1, samples=5, gc_before=false)
    @test timing.samples == 5
    @test timing.min_ns >= 0
    @test timing.max_ns >= timing.min_ns
    phase_runtime = runtime_phase_timing(runtime; warmup=1, samples=3, gc_before=false)
    @test phase_runtime.delay_mean_ns >= 0

    rng3 = MersenneTwister(3)
    tel3 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src3 = Source(band=:I, magnitude=0.0)
    atm3 = KolmogorovAtmosphere(tel3; r0=0.2, L0=25.0)
    dm3 = DeformableMirror(tel3; n_act=4, influence_width=0.3)
    wfs3 = ZernikeWFS(tel3; n_subap=4, diffraction_padding=2)
    sim3 = AdaptiveOpticsSim.AOSimulation(tel3, atm3, src3, dm3, wfs3)
    imat3 = interaction_matrix(dm3, wfs3, tel3, src3; amplitude=1e-8)
    recon3 = ModalReconstructor(imat3; gain=0.5)
    wfs_det3 = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime3 = ClosedLoopRuntime(sim3, recon3; rng=rng3, wfs_detector=wfs_det3)
    @test supports_prepared_runtime(runtime3)
    @test supports_detector_output(runtime3)
    prepare!(runtime3)
    @test runtime3.prepared
    @test runtime3.wfs.state.calibrated
    step!(runtime3)
    @test length(runtime3.command) == length(dm3.state.coefs)
    @test simulation_wfs_frame(runtime3) === wfs3.state.camera_frame
    @test simulation_wfs_frame(runtime3) == output_frame(wfs_det3)
    @test size(simulation_wfs_frame(runtime3)) == size(wfs3.state.camera_frame)
    @test all(isfinite, simulation_slopes(runtime3))

    rng4 = MersenneTwister(4)
    tel4 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src4 = Source(band=:I, magnitude=0.0)
    atm4 = KolmogorovAtmosphere(tel4; r0=0.2, L0=25.0)
    dm4 = DeformableMirror(tel4; n_act=4, influence_width=0.3)
    wfs4 = ShackHartmann(tel4; n_subap=4)
    sim4 = AdaptiveOpticsSim.AOSimulation(tel4, atm4, src4, dm4, wfs4)
    imat4 = interaction_matrix(dm4, wfs4, tel4; amplitude=0.1)
    recon4 = ModalReconstructor(imat4; gain=0.5)
    runtime4 = ClosedLoopRuntime(sim4, recon4;
        rng=rng4,
        profile=HILRuntimeProfile(),
        latency=RuntimeLatencyModel(
            measurement_delay_frames=1,
            readout_delay_frames=1,
            reconstruction_delay_frames=1,
            dm_delay_frames=1,
        ),
    )
    @test runtime_profile(runtime4) isa HILRuntimeProfile
    @test runtime4.science_zero_padding == 0
    slope_norms = Float64[]
    command_norms = Float64[]
    dm_norms = Float64[]
    for _ in 1:5
        step!(runtime4)
        push!(slope_norms, norm(simulation_slopes(runtime4)))
        push!(command_norms, norm(simulation_command(runtime4)))
        push!(dm_norms, norm(runtime4.dm.state.coefs))
    end
    @test slope_norms[1] == 0
    @test slope_norms[2] == 0
    @test slope_norms[3] > 0
    @test command_norms[1] == 0
    @test command_norms[2] == 0
    @test command_norms[3] == 0
    @test command_norms[4] > 0
    @test dm_norms[1] == 0
    @test dm_norms[2] == 0
    @test dm_norms[3] == 0
    @test dm_norms[4] == 0
    @test dm_norms[5] > 0
end

