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

@testset "Semantic backend descriptors" begin
    @test array_backend_type(CPUBackend()) === Array

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0, backend=CPUBackend())
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0, backend=CPUBackend())
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3, backend=CPUBackend())
    tt = TipTiltMirror(tel; scale=0.1, backend=CPUBackend())
    focus = FocusStage(tel; scale=0.2, backend=CPUBackend())
    wfs = ShackHartmann(tel; n_subap=4, backend=CPUBackend())
    det = Detector(backend=CPUBackend())
    sf = SpatialFilter(tel; backend=CPUBackend())
    ef = ElectricField(tel, Source(band=:I, magnitude=0.0); backend=CPUBackend())

    @test tel.state.opd isa Matrix
    @test atm.state.opd isa Matrix
    @test dm.state.coefs isa Vector
    @test tt.state.coefs isa Vector
    @test focus.state.coefs isa Vector
    @test wfs.state.slopes isa Vector
    @test det.state.frame isa Matrix
    @test sf.state.phase isa Matrix
    @test ef.state.field isa Matrix{ComplexF64}
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

    vault = CalibrationVault(imat.matrix; build_backend=AdaptiveOpticsSim.CPUBuildBackend())
    @test vault.M isa Matrix
    recon_cpu = ModalReconstructor(imat; build_backend=AdaptiveOpticsSim.CPUBuildBackend())
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
    @test command(simif) === surrogate.command
    @test length(slopes(simif)) == 2
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
    @test size(wfs_frame(curvature_readout)[1]) == (2, curvature_params.n_subap^2)
    @test wfs_metadata(curvature_readout)[1] isa CountingDetectorExportMetadata

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
    ao3k_high_meta = wfs_metadata(ao3k_iface)[1]
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
    @test command(runtime) === runtime.command
    @test science_frame(runtime) === output_frame(det)
    @test command_segment_labels(command_layout(runtime)) == (:dm,)
    @test command(runtime) === command(runtime)
    @test science_frame(runtime) === science_frame(runtime)
    @test readout(runtime).science_frame === science_frame(runtime)

    tel_ext = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src_ext = Source(band=:I, magnitude=0.0)
    atm_ext = KolmogorovAtmosphere(tel_ext; r0=0.2, L0=25.0)
    dm_ext = DeformableMirror(tel_ext; n_act=4, influence_width=0.3)
    wfs_ext = ShackHartmann(tel_ext; n_subap=4, mode=Diffractive())
    sim_ext = AdaptiveOpticsSim.AOSimulation(tel_ext, atm_ext, src_ext, dm_ext, wfs_ext)
    det_ext = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    null_recon = NullReconstructor()
    @test_throws InvalidConfiguration reconstruct!(similar(dm_ext.state.coefs), null_recon, wfs_ext.state.slopes)

    external_branch = RuntimeBranch(
        :external_branch,
        sim_ext,
        null_recon;
        science_detector=det_ext,
        rng=MersenneTwister(7),
    )
    external_cfg = SingleRuntimeConfig(
        name=:external_runtime_demo,
        branch_label=:external_branch,
        products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
    )
    external_scenario = build_runtime_scenario(external_cfg, external_branch)
    @test external_scenario isa RuntimeScenario
    prepare!(external_scenario)
    external_command = fill(eltype(command(external_scenario))(0.03), length(command(external_scenario)))
    set_command!(external_scenario, external_command)
    sense!(external_scenario)
    @test command(external_scenario) == external_command
    @test slopes(external_scenario) === slopes(external_scenario)
    @test wfs_frame(external_scenario) === wfs_frame(external_scenario)
    @test science_frame(external_scenario) === science_frame(external_scenario)
    @test science_metadata(external_scenario) == science_metadata(external_scenario)
    @test readout(external_scenario).science_frame === science_frame(external_scenario)
    @test_throws InvalidConfiguration step!(external_scenario)

    boundary = SimulationInterface(runtime)
    @test AdaptiveOpticsSim.runtime_export_plan(boundary) isa AdaptiveOpticsSim.DirectRuntimeExportPlan
    boundary_readout = readout(boundary)
    @test length(slopes(boundary)) == length(wfs.state.slopes)
    @test slopes(boundary_readout) === slopes(boundary)
    @test size(science_frame(boundary)) == size(output_frame(det))
    boundary_science_metadata = science_metadata(boundary)
    @test boundary_science_metadata isa DetectorExportMetadata
    @test boundary_science_metadata.output_size == size(output_frame(det))
    @test boundary_science_metadata.frame_size == size(det.state.frame)
    step!(boundary)
    @test command(boundary) == runtime.command

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
    @test wfs_frame(runtime2) === wfs2.state.exported_spot_cube
    @test wfs_frame(runtime2) !== wfs2.state.spot_cube
    @test wfs2.state.spot_cube !== wfs2.state.sampled_spot_cube
    @test wfs_frame(runtime2) !== wfs2.state.sampled_spot_cube
    @test size(wfs_frame(runtime2)) == size(wfs2.state.spot_cube)
    @test all(wfs_frame(runtime2) .>= wfs2.state.spot_cube)
    boundary2 = SimulationInterface(runtime2)
    @test ndims(wfs_frame(boundary2)) == 3
    @test size(wfs_frame(boundary2), 1) == wfs2.params.n_subap^2

    composite = CompositeSimulationInterface(boundary, boundary2)
    @test supports_grouped_execution(composite)
    @test length(command(composite)) == length(command(boundary)) + length(command(boundary2))
    @test length(slopes(composite)) == length(slopes(boundary)) + length(slopes(boundary2))
    @test length(wfs_frame(composite)) == 2
    @test grouped_wfs_stack(composite) === nothing
    step!(composite)
    @test command(composite) == vcat(command(boundary), command(boundary2))

    rng2b = MersenneTwister(12)
    tel2b = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src2b = Source(band=:I, magnitude=0.0)
    atm2b = KolmogorovAtmosphere(tel2b; r0=0.2, L0=25.0)
    dm2b = DeformableMirror(tel2b; n_act=4, influence_width=0.3)
    wfs2b = ShackHartmann(tel2b; n_subap=4, mode=Diffractive())
    sim2b = AdaptiveOpticsSim.AOSimulation(tel2b, atm2b, src2b, dm2b, wfs2b)
    imat2b = interaction_matrix(dm2b, wfs2b, tel2b, src2b; amplitude=0.1)
    recon2b = ModalReconstructor(imat2b; gain=0.5)
    wfs_det_b = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime2b = ClosedLoopRuntime(sim2b, recon2b; rng=rng2b, wfs_detector=wfs_det_b)
    prepare!(runtime2b)
    step!(runtime2b)
    boundary2b = SimulationInterface(runtime2b)
    grouped = CompositeSimulationInterface(
        boundary2,
        boundary2b;
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
    )
    @test length(wfs_frame(grouped)) == 2
    @test !isnothing(grouped_wfs_stack(grouped))
    @test size(grouped_wfs_stack(grouped)) == (size(wfs_frame(boundary2))..., 2)
    @test grouped_science_stack(grouped) === nothing
    step!(grouped)
    @test selectdim(grouped_wfs_stack(grouped), ndims(grouped_wfs_stack(grouped)), 1) == wfs_frame(boundary2)
    @test selectdim(grouped_wfs_stack(grouped), ndims(grouped_wfs_stack(grouped)), 2) == wfs_frame(boundary2b)

    grouped_stack_only = CompositeSimulationInterface(
        boundary2,
        boundary2b;
        products=GroupedRuntimeProductRequirements(wfs_frames=false, science_frames=false, wfs_stack=true, science_stack=false),
    )
    @test wfs_frame(grouped_stack_only) === nothing
    @test !isnothing(grouped_wfs_stack(grouped_stack_only))

    single_cfg = SingleRuntimeConfig(
        name=:single_runtime_demo,
        branch_label=:science_branch,
        products=RuntimeProductRequirements(slopes=true, wfs_pixels=false, science_pixels=true),
    )
    single_branch = RuntimeBranch(
        :science_branch,
        sim,
        recon;
        science_detector=det,
        rng=MersenneTwister(31),
    )
    single_scenario = build_runtime_scenario(single_cfg, single_branch)
    @test single_scenario isa RuntimeScenario
    @test platform_name(single_scenario) == :single_runtime_demo
    @test platform_branch_labels(single_scenario) == (:science_branch,)
    @test simulation_interface(single_scenario) isa SimulationInterface
    @test supports_detector_output(single_scenario)
    single_layout = command_layout(single_scenario)
    @test single_layout isa RuntimeCommandLayout
    @test command_segment_labels(single_layout) == (:dm,)
    @test command_segment_range(command_segments(single_layout)[1]) == 1:length(command(single_scenario))
    branch_command = fill(eltype(command(single_scenario))(0.02), length(command(single_scenario)))
    set_command!(single_scenario, (; dm=branch_command))
    sense!(single_scenario)
    @test command(single_scenario) == branch_command
    step!(single_scenario)
    @test science_frame(single_scenario) !== nothing
    @test wfs_frame(single_scenario) === nothing
    @test readout(single_scenario).science_frame === science_frame(single_scenario)

    split_layout = RuntimeCommandLayout(:woofer => 8, :tweeter => 8)
    split_scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:split_runtime_demo, branch_label=:science_branch,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=false, science_pixels=true)),
        RuntimeBranch(:science_branch, sim, NullReconstructor();
            science_detector=det,
            rng=MersenneTwister(32),
            command_layout=split_layout),
    )
    @test command_segment_labels(command_layout(split_scenario)) == (:woofer, :tweeter)
    woofer_cmd = fill(eltype(command(split_scenario))(0.03), 8)
    tweeter_cmd = fill(eltype(command(split_scenario))(0.04), 8)
    set_command!(split_scenario, (; woofer=woofer_cmd, tweeter=tweeter_cmd))
    sense!(split_scenario)
    @test command(split_scenario) == vcat(woofer_cmd, tweeter_cmd)
    @test science_frame(split_scenario) !== nothing

    tiptilt = TipTiltMirror(tel_ext; scale=0.1, label=:tiptilt)
    dm_combo = DeformableMirror(tel_ext; n_act=4, influence_width=0.3)
    combo_optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm_combo)
    combo_sim = AdaptiveOpticsSim.AOSimulation(tel_ext, atm_ext, src_ext, combo_optic, wfs_ext)
    combo_scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:combo_runtime_demo, branch_label=:external_branch,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true)),
        RuntimeBranch(:external_branch, combo_sim, NullReconstructor();
            wfs_detector=det_ext,
            science_detector=det_ext,
            rng=MersenneTwister(8)),
    )
    prepare!(combo_scenario)
    @test command_segment_labels(command_layout(combo_scenario)) == (:tiptilt, :dm)
    set_command!(combo_scenario, (; tiptilt=fill(eltype(command(combo_scenario))(0.01), 2), dm=fill(eltype(command(combo_scenario))(0.02), 16)))
    sense!(combo_scenario)
    @test command(combo_scenario) == vcat(fill(eltype(command(combo_scenario))(0.01), 2), fill(eltype(command(combo_scenario))(0.02), 16))
    update_command!(combo_scenario, (; tiptilt=fill(eltype(command(combo_scenario))(0.03), 2)))
    @test command(combo_scenario)[1:2] == fill(eltype(command(combo_scenario))(0.03), 2)
    @test command(combo_scenario)[3:end] == fill(eltype(command(combo_scenario))(0.02), 16)
    sense!(combo_scenario)
    @test science_frame(combo_scenario) !== nothing
    @test wfs_frame(combo_scenario) !== nothing

    focus = FocusStage(tel_ext; scale=0.1, label=:focus)
    dm_focus = DeformableMirror(tel_ext; n_act=4, influence_width=0.3)
    focus_optic = CompositeControllableOptic(:focus => focus, :dm => dm_focus)
    focus_sim = AdaptiveOpticsSim.AOSimulation(tel_ext, atm_ext, src_ext, focus_optic, wfs_ext)
    focus_scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:focus_runtime_demo, branch_label=:external_branch,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:external_branch, focus_sim, NullReconstructor();
            wfs_detector=det_ext,
            rng=MersenneTwister(9)),
    )
    prepare!(focus_scenario)
    @test command_segment_labels(command_layout(focus_scenario)) == (:focus, :dm)
    set_command!(focus_scenario, (; focus=fill(eltype(command(focus_scenario))(0.04), 1), dm=fill(eltype(command(focus_scenario))(0.01), 16)))
    sense!(focus_scenario)
    @test wfs_frame(focus_scenario) !== nothing

    grouped_cfg = GroupedRuntimeConfig(
        (:branch_a, :branch_b);
        name=:grouped_runtime_demo,
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
    )
    grouped_scenario = build_runtime_scenario(
        grouped_cfg,
        RuntimeBranch(:branch_a, sim2, recon2; wfs_detector=wfs_det, rng=MersenneTwister(41)),
        RuntimeBranch(:branch_b, sim2b, recon2b; wfs_detector=wfs_det_b, rng=MersenneTwister(42)),
    )
    @test grouped_scenario isa RuntimeScenario
    @test platform_name(grouped_scenario) == :grouped_runtime_demo
    @test platform_branch_labels(grouped_scenario) == (:branch_a, :branch_b)
    @test simulation_interface(grouped_scenario) isa CompositeSimulationInterface
    @test supports_grouped_execution(grouped_scenario)
    grouped_layout = command_layout(grouped_scenario)
    @test grouped_layout isa RuntimeCommandLayout
    @test command_segment_labels(grouped_layout) == (:branch_a, :branch_b)
    @test command_segment_range(command_segments(grouped_layout)[1]) == 1:length(simulation_interface(grouped_scenario).interfaces[1].command)
    @test first(command_segment_range(command_segments(grouped_layout)[2])) == length(simulation_interface(grouped_scenario).interfaces[1].command) + 1
    @test command_segment_labels(branch_command_layout(grouped_scenario, :branch_a)) == (:dm,)
    @test branch_command_layouts(grouped_scenario).branch_b == branch_command_layout(grouped_scenario, :branch_b)
    prepare!(grouped_scenario)
    branch_a_cmd = fill(eltype(command(grouped_scenario))(0.01), length(simulation_interface(grouped_scenario).interfaces[1].command))
    branch_b_cmd = fill(eltype(command(grouped_scenario))(0.02), length(simulation_interface(grouped_scenario).interfaces[2].command))
    set_command!(grouped_scenario, (; branch_a=branch_a_cmd, branch_b=branch_b_cmd))
    sense!(grouped_scenario)
    @test command(grouped_scenario) == vcat(branch_a_cmd, branch_b_cmd)
    step!(grouped_scenario)
    @test !isnothing(grouped_wfs_stack(grouped_scenario))
    @test size(grouped_wfs_stack(grouped_scenario)) == (size(wfs_frame(grouped_scenario)[1])..., 2)
    @test grouped_science_stack(grouped_scenario) === nothing
    grouped_readout = readout(grouped_scenario)
    @test grouped_readout.grouped_wfs_stack === grouped_wfs_stack(grouped_scenario)

    split_grouped_scenario = build_runtime_scenario(
        grouped_cfg,
        RuntimeBranch(:branch_a, sim2, NullReconstructor();
            wfs_detector=wfs_det,
            rng=MersenneTwister(43),
            command_layout=RuntimeCommandLayout(:woofer => 8, :tweeter => 8)),
        RuntimeBranch(:branch_b, sim2b, NullReconstructor();
            wfs_detector=wfs_det_b,
            rng=MersenneTwister(44),
            command_layout=RuntimeCommandLayout(:steering => 2, :dm => 14)),
    )
    prepare!(split_grouped_scenario)
    @test command_segment_labels(command_layout(split_grouped_scenario)) == (:branch_a, :branch_b)
    @test command_segment_labels(branch_command_layout(split_grouped_scenario, :branch_a)) == (:woofer, :tweeter)
    @test command_segment_labels(branch_command_layout(split_grouped_scenario, :branch_b)) == (:steering, :dm)
    set_command!(split_grouped_scenario, (
        branch_a=(; woofer=fill(eltype(command(split_grouped_scenario))(0.05), 8), tweeter=fill(eltype(command(split_grouped_scenario))(0.06), 8)),
        branch_b=(; steering=fill(eltype(command(split_grouped_scenario))(0.07), 2), dm=fill(eltype(command(split_grouped_scenario))(0.08), 14)),
    ))
    sense!(split_grouped_scenario)
    @test command(split_grouped_scenario) == vcat(fill(eltype(command(split_grouped_scenario))(0.05), 8), fill(eltype(command(split_grouped_scenario))(0.06), 8), fill(eltype(command(split_grouped_scenario))(0.07), 2), fill(eltype(command(split_grouped_scenario))(0.08), 14))
    @test !isnothing(grouped_wfs_stack(split_grouped_scenario))

    grouped_tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    grouped_pyr = PyramidWFS(grouped_tel; n_subap=4, modulation=1.0, mode=Diffractive())
    grouped_bio = BioEdgeWFS(grouped_tel; n_subap=4, modulation=1.0, mode=Diffractive())
    @test @inferred(AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(grouped_pyr.state.intensity), grouped_pyr)) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test @inferred(AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(grouped_bio.state.intensity), grouped_bio)) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.reduction_execution_plan(grouped_pyr.state.intensity) isa AdaptiveOpticsSim.DirectReductionPlan
    @test AdaptiveOpticsSim.runtime_export_plan(grouped) isa AdaptiveOpticsSim.CompositeRuntimeExportPlan

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
    @test wfs_frame(runtime2_slopes_only) === nothing
    slopes_only_boundary = SimulationInterface(runtime2_slopes_only)
    @test wfs_frame(slopes_only_boundary) === nothing

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
    @test wfs_frame(runtime3) === wfs3.state.camera_frame
    @test wfs_frame(runtime3) == output_frame(wfs_det3)
    @test size(wfs_frame(runtime3)) == size(wfs3.state.camera_frame)
    @test all(isfinite, slopes(runtime3))

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
        push!(slope_norms, norm(slopes(runtime4)))
        push!(command_norms, norm(command(runtime4)))
        push!(dm_norms, norm(command_storage(runtime4.optic)))
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
