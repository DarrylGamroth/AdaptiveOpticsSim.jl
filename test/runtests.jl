using Test
using AdaptiveOpticsSim
using Random
using SpecialFunctions
using Statistics
using Tables
using TOML

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
using .SubaruAO188Simulation

include("reference_harness.jl")
include("ka_cpu_matrix.jl")
include("tomography.jl")

function run_tutorial_example(name::AbstractString)
    path = joinpath(dirname(@__DIR__), "examples", "tutorials", name)
    mod = Module(Symbol("Tutorial_", replace(basename(name), "." => "_")))
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Base.include(mod, path)
    main_fn = Core.eval(mod, :main)
    return Base.invokelatest(main_fn)
end

function assert_source_interface(src)
    @test hasmethod(wavelength, Tuple{typeof(src)})
end

function assert_atmosphere_interface(atm, tel)
    @test applicable(advance!, atm, tel)
    @test applicable(propagate!, atm, tel)
end

function assert_wfs_interface(wfs, tel)
    @test applicable(update_valid_mask!, wfs, tel)
    @test applicable(measure!, wfs, tel)
end

function assert_detector_interface(det, psf)
    @test applicable(capture!, det, psf)
end

function assert_dm_interface(dm, tel)
    @test applicable(build_influence_functions!, dm, tel)
    @test applicable(apply!, dm, tel, DMAdditive())
end

function subharmonic_tiptilt_power(phs::AbstractMatrix)
    n = size(phs, 1)
    coords = collect(LinRange(-1.0, 1.0, n))
    x = repeat(reshape(coords, 1, n), n, 1)
    y = repeat(reshape(coords, n, 1), 1, n)
    A = hcat(vec(x), vec(y), ones(n * n))
    coeffs = A \ vec(phs)
    return coeffs[1]^2 + coeffs[2]^2
end

function subharmonic_metrics(atm::KolmogorovAtmosphere, D::Real;
    n::Int=24, nsamp::Int=8, kwargs...)
    delta = D / n
    shift = n ÷ 2
    ws = PhaseStatsWorkspace(n; T=Float64)
    tiptilt = Float64[]
    structure = Float64[]
    for seed in 1:nsamp
        phs = ft_sh_phase_screen(atm, n, delta; rng=MersenneTwister(seed), ws=ws, kwargs...)
        push!(tiptilt, subharmonic_tiptilt_power(phs))
        diff = @views phs[:, 1:end-shift] .- phs[:, 1+shift:end]
        push!(structure, Statistics.mean(abs2, diff))
    end
    return (; tiptilt=Statistics.mean(tiptilt), structure=Statistics.mean(structure))
end

@test AdaptiveOpticsSim.PROJECT_STATUS == :in_development

@testset "GPU backend registry" begin
    @test !gpu_backend_loaded(CUDABackendTag)
    @test !gpu_backend_loaded(MetalBackendTag)
    @test !gpu_backend_loaded(AMDGPUBackendTag)
    @test gpu_backend_array_type(CUDABackendTag) === nothing
    @test gpu_backend_array_type(MetalBackendTag) === nothing
    @test gpu_backend_array_type(AMDGPUBackendTag) === nothing
    @test gpu_backend_name(Matrix{Float64}) === nothing
    @test available_gpu_backends() == ()
    @test GPUArrayBuildBackend(CUDABackendTag) isa GPUArrayBuildBackend{CUDABackendTag}
end

@testset "Telescope and PSF" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    psf = compute_psf!(tel, src; zero_padding=2)
    @test size(psf) == (64, 64)
    @test maximum(psf) > 0
    @test isfinite(sum(psf))
    @test tel.state.psf_workspace !== nothing
    cached_ws = tel.state.psf_workspace
    compute_psf!(tel, src; zero_padding=2)
    @test tel.state.psf_workspace === cached_ws

    tel_dim = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2,
        pupil_reflectivity=0.25)
    psf_dim = compute_psf!(tel_dim, src; zero_padding=2)
    @test sum(psf_dim) ≈ 0.25 * sum(psf)
    fmap = flux_map(tel_dim, src)
    @test size(fmap) == size(tel_dim.state.pupil)
    @test maximum(fmap) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"
end

@testset "Zernike basis" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end

@testset "Atmosphere propagation" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (32, 32)
    @test sum(abs.(tel.state.opd)) > 0

    delta = tel.params.diameter / tel.params.resolution
    ensure_psd!(atm, delta)
    psd_snapshot = copy(atm.state.psd)
    ensure_psd!(atm, delta)
    @test psd_snapshot == atm.state.psd

    atm = KolmogorovAtmosphere(tel; r0=0.1, L0=25.0)
    ensure_psd!(atm, delta)
    @test psd_snapshot != atm.state.psd
end

@testset "Sub-harmonic phase screens" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    delta = tel.params.diameter / tel.params.resolution
    rng = MersenneTwister(11)
    ws = PhaseStatsWorkspace(32; T=Float64)
    phs_base = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=false)
    rng = MersenneTwister(11)
    phs_sh = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=true,
        mode=FidelitySubharmonics())
    @test sum(abs.(phs_sh .- phs_base)) > 0
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(25.0, 8.0) == 4
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(200.0, 8.0) > 4

    atm_large = KolmogorovAtmosphere(tel; r0=0.2, L0=200.0)
    rng = MersenneTwister(11)
    phs_default = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FidelitySubharmonics())
    rng = MersenneTwister(11)
    phs_legacy = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FastSubharmonics())
    @test sum(abs.(phs_default .- phs_legacy)) > 0

    metrics_25_none = subharmonic_metrics(atm, tel.params.diameter; subharmonics=false)
    metrics_25_legacy = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_25_default = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
    @test metrics_25_none.tiptilt < metrics_25_legacy.tiptilt < metrics_25_default.tiptilt
    @test metrics_25_none.structure < metrics_25_legacy.structure < metrics_25_default.structure

    metrics_200_none = subharmonic_metrics(atm_large, tel.params.diameter; subharmonics=false)
    metrics_200_legacy = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_200_default = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
@test metrics_200_none.tiptilt < metrics_200_legacy.tiptilt < metrics_200_default.tiptilt
@test metrics_200_none.structure < metrics_200_legacy.structure < metrics_200_default.structure
end

@testset "Fidelity profiles" begin
    @test default_fidelity_profile() isa ScientificProfile
    @test default_subharmonic_mode(ScientificProfile()) isa FidelitySubharmonics
    @test default_subharmonic_mode(FastProfile()) isa FastSubharmonics
    @test default_ncpa_basis(ScientificProfile()).method isa KLHHtPSD
    @test default_ncpa_basis(FastProfile()).method isa KLDMModes

    mixed = ProfileBundle(ScientificProfile(); lift=FastProfile(), tomography=FastProfile())
    @test atmosphere_profile(mixed) isa ScientificProfile
    @test calibration_profile(mixed) isa ScientificProfile
    @test detector_profile(mixed) isa ScientificProfile
    @test lift_profile(mixed) isa FastProfile
    @test tomography_profile(mixed) isa FastProfile
    @test default_subharmonic_mode(mixed) isa FidelitySubharmonics
    @test default_ncpa_basis(mixed).method isa KLHHtPSD

    fast_cal = ProfileBundle(ScientificProfile(); calibration=FastProfile())
    @test default_ncpa_basis(fast_cal).method isa KLDMModes

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.15, L0=200.0)
    scientific = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=ScientificProfile())
    fast = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=FastProfile())
    @test !all(scientific .== fast)

    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    coeffs = zeros(4)
    ncpa_fast = NCPA(tel, dm, atm; profile=FastProfile(), coefficients=coeffs)
    ncpa_scientific = NCPA(tel, dm, atm; profile=ScientificProfile(), coefficients=coeffs)
    @test size(ncpa_fast.opd) == size(tel.state.opd)
    @test size(ncpa_scientific.opd) == size(tel.state.opd)
end

@testset "Multi-layer atmosphere" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_r0=[0.5, 0.5],
        wind_speed=[5.0, 10.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (16, 16)
    @test sum(abs.(tel.state.opd)) > 0
end

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
end

function closed_loop_runtime_allocations()
    rng = MersenneTwister(0)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmann(tel; n_subap=4)
    sim = AOSimulation(tel, atm, src, dm, wfs)
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
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime = ClosedLoopRuntime(sim, recon; rng=rng, science_detector=det)

    step!(runtime)
    @test length(runtime.command) == length(dm.state.coefs)
    @test size(output_frame(det)) == (32, 32)
    @test closed_loop_runtime_allocations() == 0

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
    sim2 = AOSimulation(tel2, atm2, src2, dm2, wfs2)
    imat2 = interaction_matrix(dm2, wfs2, tel2, src2; amplitude=0.1)
    recon2 = ModalReconstructor(imat2; gain=0.5)
    wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime2 = ClosedLoopRuntime(sim2, recon2; rng=rng2, wfs_detector=wfs_det)
    step!(runtime2)
    boundary2 = SimulationInterface(runtime2)
    @test ndims(simulation_wfs_frame(boundary2)) == 3
    @test size(simulation_wfs_frame(boundary2), 1) == wfs2.params.n_subap^2

    composite = CompositeSimulationInterface(boundary, boundary2)
    @test length(simulation_command(composite)) == length(simulation_command(boundary)) + length(simulation_command(boundary2))
    @test length(simulation_slopes(composite)) == length(simulation_slopes(boundary)) + length(simulation_slopes(boundary2))
    @test length(simulation_wfs_frame(composite)) == 2
    step!(composite)
    @test simulation_command(composite) == vcat(simulation_command(boundary), simulation_command(boundary2))

    timing = runtime_timing(runtime; warmup=1, samples=5, gc_before=false)
    @test timing.samples == 5
    @test timing.min_ns >= 0
    @test timing.max_ns >= timing.min_ns
end

@testset "Detector" begin
    psf = fill(1.0, 8, 8)
    det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=2)
    frame = capture!(det, psf; rng=MersenneTwister(2))
    @test size(frame) == (4, 4)
    @test sum(frame) == sum(psf)

    det_sampling = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_sampling = capture!(det_sampling, psf; rng=MersenneTwister(2))
    @test size(frame_sampling) == (2, 2)
    @test sum(frame_sampling) == sum(psf)

    det_tuple = Detector(integration_time=1.0, noise=(NoisePhoton(), NoiseReadout(0.5)),
        qe=1.0, binning=1)
    @test det_tuple.noise isa NoisePhotonReadout

    det_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, full_well=5.0)
    frame_sat = capture!(det_sat, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test maximum(frame_sat) == 5.0

    det_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, bits=8, full_well=10.0)
    frame_adc = capture!(det_adc, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc isa Matrix{UInt8}
    @test output_frame(det_adc) === frame_adc
    @test maximum(frame_adc) == 0xff
    @test minimum(frame_adc) >= 0x00
    @test eltype(det_adc.state.frame) == Float64
    metadata_adc = detector_export_metadata(det_adc)
    @test metadata_adc.noise == :none
    @test metadata_adc.sensor == :ccd
    @test metadata_adc.output_precision == UInt8
    @test metadata_adc.frame_size == (4, 4)
    @test metadata_adc.output_size == (4, 4)

    det_adc_float = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_precision=Float32)
    frame_adc_float = capture!(det_adc_float, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc_float isa Matrix{Float32}
    @test maximum(frame_adc_float) == Float32(255.0)

    det_dark = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, dark_current=100.0)
    frame_dark = capture!(det_dark, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_dark) > 0

    zero_psf = zeros(4, 4)
    rng_ccd = MersenneTwister(7)
    rng_emccd = MersenneTwister(7)
    det_ccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=CCDSensor())
    det_emccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor())
    frame_ccd = copy(capture!(det_ccd, zero_psf; rng=rng_ccd))
    frame_emccd = copy(capture!(det_emccd, zero_psf; rng=rng_emccd))
    @test frame_ccd ≈ 10 .* frame_emccd

    det_buffered = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0, binning=1)
    frame_partial = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test !readout_ready(det_buffered)
    @test sum(frame_partial) == 16.0
    frame_buffered = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test readout_ready(det_buffered)
    @test sum(frame_buffered) == 32.0
    reset_integration!(det_buffered)
    @test readout_ready(det_buffered)
    @test det_buffered.state.integrated_time == 0.0
    metadata_buffered = detector_export_metadata(det_buffered)
    @test metadata_buffered.output_size == size(output_frame(det_buffered))
    @test metadata_buffered.psf_sampling == 1
    @test metadata_buffered.binning == 1

    det_background_flux = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_flux=2.0)
    frame_background_flux = capture!(det_background_flux, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_background_flux) > 0

    det_background_map = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_map=fill(1.0, 4, 4))
    frame_background_map = capture!(det_background_map, zeros(4, 4); rng=MersenneTwister(2))
    @test frame_background_map == fill(-1.0, 4, 4)

    cube = cat(fill(1.0, 4, 4), fill(2.0, 4, 4); dims=3)
    scratch = similar(cube)
    det_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=0.5, binning=1)
    AdaptiveOpticsSim.capture_stack!(det_stack, cube, scratch; rng=MersenneTwister(10))
    @test cube[:, :, 1] ≈ fill(0.5, 4, 4)
    @test cube[:, :, 2] ≈ fill(1.0, 4, 4)

    psf = reshape(Float64.(1:256), 16, 16)
    det_fused = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_fused = copy(AdaptiveOpticsSim.fill_frame!(det_fused, psf, 1.0))
    manual_mid = zeros(Float64, 8, 8)
    manual_out = zeros(Float64, 4, 4)
    AdaptiveOpticsSim.bin2d!(manual_mid, psf, 2)
    AdaptiveOpticsSim.bin2d!(manual_out, manual_mid, 2)
    @test frame_fused == manual_out
end

@testset "Asterism PSF" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    ast = Asterism([src1, src2])
    psf = compute_psf!(tel, ast; zero_padding=2)
    @test size(tel.state.psf_stack, 3) == 2
    @test size(psf) == (32, 32)
    @test sum(psf) >= sum(@view tel.state.psf_stack[:, :, 1])
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    base_sum = sum(tel.state.pupil)
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(tel.state.pupil) < base_sum

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(tel.state.pupil) == sum(custom)
    @test tel.state.pupil_reflectivity == Float64.(custom)

    reflectivity = fill(0.5, 16, 16)
    set_pupil_reflectivity!(tel, reflectivity)
    @test tel.state.pupil_reflectivity[:, 1:8] == fill(0.5, 16, 8)
    @test tel.state.pupil_reflectivity[:, 9:end] == fill(0.0, 16, 8)

    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; n_subap=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    sh = ShackHartmann(tel; n_subap=4)
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = sh.params.n_subap * sh.params.n_subap
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5)

    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_profile = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    sh_sampled = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_sampled_slopes = measure!(sh_sampled, tel, ngs)
    @test length(sh_sampled_slopes) == 2 * 4 * 4

    pyr_sampled = PyramidWFS(tel; n_subap=4, mode=Diffractive(), n_pix_separation=4, binning=2)
    pyr_sampled_slopes = measure!(pyr_sampled, tel, ngs)
    @test length(pyr_sampled_slopes) == 2 * count(pyr_sampled.state.valid_i4q)
    pyr_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    pyr_frame = copy(AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_sampled, tel, pyr_intensity))
    pyr_camera = zeros(Float64, 4, 4)
    pyr_manual = zeros(Float64, 2, 2)
    AdaptiveOpticsSim.bin2d!(pyr_camera, pyr_intensity, 8)
    AdaptiveOpticsSim.bin2d!(pyr_manual, pyr_camera, 2)
    @test pyr_frame == pyr_manual

    bio_sampled = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), binning=2)
    bio_sampled_slopes = measure!(bio_sampled, tel, ngs)
    @test length(bio_sampled_slopes) == 2 * count(bio_sampled.state.valid_i4q)
    bio_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    bio_frame = copy(AdaptiveOpticsSim.sample_bioedge_intensity!(bio_sampled, tel, bio_intensity))
    bio_camera = zeros(Float64, 4, 4)
    bio_manual = similar(bio_frame)
    AdaptiveOpticsSim.bin2d!(bio_camera, bio_intensity, 8)
    AdaptiveOpticsSim.bin2d!(bio_manual, bio_camera, div(size(bio_camera, 1), size(bio_frame, 1)))
    @test bio_frame == bio_manual

    pyr_profile = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_profile_slopes = measure!(pyr_profile, tel, lgs_profile)
    @test all(isfinite, pyr_profile_slopes)

    bio_profile = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_profile_slopes = measure!(bio_profile, tel, lgs_profile)
    @test all(isfinite, bio_profile_slopes)

    pyr = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_det_slopes = measure!(sh_det, tel, ngs, det)
    @test length(sh_det_slopes) == 2 * 4 * 4
    pyr_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_det_slopes = measure!(pyr_det, tel, ngs, det)
    @test length(pyr_det_slopes) == 2 * 4 * 4

    ast = Asterism([ngs, Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))])
    sh_ast = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_ast_slopes = copy(measure!(sh_ast, tel, ast))
    @test length(sh_ast_slopes) == 2 * 4 * 4
    sh_ast_serial = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_serial, tel, ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_serial, tel, ast.sources[1])
    fill!(sh_ast_serial.state.detector_noise_cube, zero(eltype(sh_ast_serial.state.detector_noise_cube)))
    for src in ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_ast_serial, tel, src)
        sh_ast_serial.state.detector_noise_cube .+= sh_ast_serial.state.spot_cube
    end
    copyto!(sh_ast_serial.state.spot_cube, sh_ast_serial.state.detector_noise_cube)
    sh_ast_serial_peak = maximum(sh_ast_serial.state.spot_cube)
    AdaptiveOpticsSim.sh_signal_from_spots!(sh_ast_serial, sh_ast_serial_peak, sh_ast_serial.params.threshold_cog)
    AdaptiveOpticsSim.subtract_reference_and_scale!(sh_ast_serial)
    sh_ast_serial_slopes = copy(sh_ast_serial.state.slopes)
    @test sh_ast_slopes ≈ sh_ast_serial_slopes
    mixed_ngs = Source(wavelength=wavelength(lgs_profile), magnitude=0.0, coordinates=(0.0, 0.0))
    mixed_ast = Asterism([mixed_ngs, lgs_profile])
    sh_mixed_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_mixed_det_slopes = copy(measure!(sh_mixed_det, tel, mixed_ast, det; rng=MersenneTwister(14)))
    sh_mixed_det_frame = copy(sh_mixed_det.state.spot_cube)
    sh_mixed_det_manual = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    fill!(sh_mixed_det_manual.state.detector_noise_cube, zero(eltype(sh_mixed_det_manual.state.detector_noise_cube)))
    for src in mixed_ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_mixed_det_manual, tel, src, det, MersenneTwister(14))
        sh_mixed_det_manual.state.detector_noise_cube .+= sh_mixed_det_manual.state.spot_cube
    end
    copyto!(sh_mixed_det_manual.state.spot_cube, sh_mixed_det_manual.state.detector_noise_cube)
    @test sh_mixed_det_frame ≈ sh_mixed_det_manual.state.spot_cube
    @test length(sh_mixed_det_slopes) == 2 * 4 * 4
    pyr_ast = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_slopes = copy(measure!(pyr_ast, tel, ast))
    @test length(pyr_ast_slopes) == 2 * 4 * 4
    pyr_ast_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_serial, tel, ast.sources[1])
    fill!(pyr_ast_serial.state.intensity, zero(eltype(pyr_ast_serial.state.intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.pyramid_intensity!(pyr_ast_serial.state.temp, pyr_ast_serial, tel, src)
        pyr_ast_serial.state.intensity .+= pyr_ast_serial.state.temp
    end
    pyr_ast_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_serial, tel, pyr_ast_serial.state.intensity)
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_serial, tel, pyr_ast_intensity)
    @. pyr_ast_serial.state.slopes *= pyr_ast_serial.state.optical_gain
    @test pyr_ast_slopes ≈ pyr_ast_serial.state.slopes
    bio_ast = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_slopes = copy(measure!(bio_ast, tel, ast))
    @test length(bio_ast_slopes) == 2 * 4 * 4
    bio_ast_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_serial, tel, ast.sources[1])
    fill!(bio_ast_serial.state.binned_intensity, zero(eltype(bio_ast_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_serial.state.intensity, bio_ast_serial, tel, src)
        bio_ast_serial.state.binned_intensity .+= bio_ast_serial.state.intensity
    end
    bio_ast_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_serial, tel, bio_ast_serial.state.binned_intensity)
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_serial, tel, bio_ast_intensity)
    @test bio_ast_slopes ≈ bio_ast_serial.state.slopes

    pyr_ast_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_det_slopes = copy(measure!(pyr_ast_det, tel, ast, det))
    pyr_ast_det_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_det_serial, tel, ast.sources[1])
    fill!(pyr_ast_det_serial.state.intensity, zero(eltype(pyr_ast_det_serial.state.intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.pyramid_intensity!(pyr_ast_det_serial.state.temp, pyr_ast_det_serial, tel, src)
        pyr_ast_det_serial.state.intensity .+= pyr_ast_det_serial.state.temp
    end
    pyr_ast_det_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_det_serial, tel, pyr_ast_det_serial.state.intensity)
    pyr_ast_det_frame = capture!(det, pyr_ast_det_intensity; rng=MersenneTwister(12))
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_ast_det_serial, size(pyr_ast_det_frame, 1))
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_det_serial, tel, pyr_ast_det_frame)
    @. pyr_ast_det_serial.state.slopes *= pyr_ast_det_serial.state.optical_gain
    @test pyr_ast_det_slopes ≈ pyr_ast_det_serial.state.slopes

    bio_ast_det = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_det_slopes = copy(measure!(bio_ast_det, tel, ast, det; rng=MersenneTwister(13)))
    bio_ast_det_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_det_serial, tel, ast.sources[1])
    fill!(bio_ast_det_serial.state.binned_intensity, zero(eltype(bio_ast_det_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_det_serial.state.intensity, bio_ast_det_serial, tel, src)
        bio_ast_det_serial.state.binned_intensity .+= bio_ast_det_serial.state.intensity
    end
    bio_ast_det_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_det_serial, tel, bio_ast_det_serial.state.binned_intensity)
    bio_ast_det_frame = capture!(det, bio_ast_det_intensity; rng=MersenneTwister(13))
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_ast_det_serial, size(bio_ast_det_frame, 1))
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_det_serial, tel, bio_ast_det_frame)
    @test bio_ast_det_slopes ≈ bio_ast_det_serial.state.slopes
end

@testset "OOPAO parity knobs" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i + j / 10
    end

    sh_plain = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_shift = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        half_pixel_shift=true)
    sh_thresh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        threshold_cog=0.2)
    @test measure!(sh_plain, tel, src) != measure!(sh_shift, tel, src)
    @test measure!(sh_plain, tel, src) != measure!(sh_thresh, tel, src)

    pyr_auto = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    @test size(pyr_auto.state.modulation_phases, 3) == 8

    pyr_path = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=0.0,
        user_modulation_path=((1.0, 0.0), (0.0, 1.0)))
    @test size(pyr_path.state.modulation_phases, 3) == 2

    pyr_default = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    pyr_rooftop = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        rooftop=0.5, theta_rotation=0.2)
    pyr_old = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0, old_mask=true)
    @test pyr_default.state.pyramid_mask != pyr_rooftop.state.pyramid_mask
    @test pyr_default.state.pyramid_mask != pyr_old.state.pyramid_mask

    bio_plain = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    bio_gray = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bio_gray.state.bioedge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bio_plain.state.bioedge_masks != bio_gray.state.bioedge_masks

    bio_gray_slopes = measure!(bio_gray, tel, src)
    @test length(bio_gray_slopes) == 2 * 4 * 4
    @test all(isfinite, bio_gray_slopes)
end
@testset "Calibration vault and modal basis" begin
    D = rand(4, 3)
    vault = CalibrationVault(D)
    @test size(vault.M) == (3, 4)
    @test vault.cond > 0
    @test vault.effective_rank == 3
    vault_trunc = with_truncation(vault, 1)
    @test vault_trunc.n_trunc == 1

    D_sing = [1.0 0.0; 0.0 1e-12]
    vault_exact = CalibrationVault(D_sing; policy=ExactPseudoInverse())
    vault_tsvd = CalibrationVault(D_sing; policy=TSVDInverse(rtol=1e-9))
    @test vault_exact.effective_rank == 2
    @test vault_tsvd.effective_rank == 1
    @test maximum(abs, vault_exact.M .- vault_tsvd.M) > 0
    @test CalibrationVault(D_sing).policy isa TSVDInverse

    imat = InteractionMatrix(D_sing, 0.1)
    recon_exact = ModalReconstructor(imat; policy=ExactPseudoInverse())
    recon_tsvd = ModalReconstructor(imat; policy=TSVDInverse(rtol=1e-9))
    @test recon_exact.effective_rank == 2
    @test recon_tsvd.effective_rank == 1
    @test ModalReconstructor(imat).policy isa TSVDInverse

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    basis = modal_basis(dm, tel; n_modes=2)
    @test size(basis.M2C, 2) == 2
    opd = rand(16, 16)
    fit, corr, turb = fitting_error(opd, basis.projector, basis.basis)
    @test size(fit) == size(opd)
    @test size(corr) == size(opd)
    @test size(turb) == size(opd)

    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    M2C, basis_hht = kl_modal_basis(KLHHtPSD(), dm, tel, atm; n_modes=2)
    @test size(M2C, 2) == 2
    @test size(basis_hht, 3) == 2
    basis2 = modal_basis(dm, tel; n_modes=2, method=KLHHtPSD(), atm=atm)
    @test size(basis2.M2C, 2) == 2
end

@testset "OPD maps and NCPA" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    map = OPDMap(fill(1.0, 8, 8))
    apply!(map, tel, DMReplace())
    @test sum(tel.state.opd) ≈ 64.0

    basis_default = AdaptiveOpticsSim.ncpa_basis(KLBasis(), tel, dm, atm; n_modes=2)
    basis_hht = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLHHtPSD()), tel, dm, atm; n_modes=2)
    basis_dm = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLDMModes()), tel, dm, atm; n_modes=2)
    @test basis_default ≈ basis_hht
    @test sum(abs.(basis_default .- basis_dm)) > 0

    coeffs = [1e-9, 2e-9]
    ncpa_default_kl = NCPA(tel, dm, atm; basis=KLBasis(), coefficients=coeffs)
    ncpa_hht = NCPA(tel, dm, atm; basis=KLBasis(KLHHtPSD()), coefficients=coeffs)
    ncpa_dm = NCPA(tel, dm, atm; basis=KLBasis(KLDMModes()), coefficients=coeffs)
    @test ncpa_default_kl.opd ≈ ncpa_hht.opd
    @test sum(abs.(ncpa_default_kl.opd .- ncpa_dm.opd)) > 0

    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 1e-9, 2e-9])
    @test size(ncpa.opd) == (8, 8)
end

@testset "Spatial filter" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    sf = SpatialFilter(tel; shape=CircularFilter(), diameter=4, zero_padding=2)
    phase, amp = filter!(sf, tel, src)
    @test size(phase) == (8, 8)
    @test size(amp) == (8, 8)
end

@testset "Gain sensing camera" begin
    mask = ones(8, 8)
    basis = rand(8, 8, 3)
    gsc = GainSensingCamera(mask, basis)
    frame = abs.(randn(8, 8))
    calibrate!(gsc, frame)
    og = compute_optical_gains!(gsc, frame)
    @test length(og) == 3
    @test length(weak_mode_mask(gsc)) == 3
    @test all(isfinite, og)
    @test detector_metadata(gsc) === nothing

    weak_gsc = GainSensingCamera(mask, zeros(8, 8, 2); sensitivity_floor=1e-6)
    calibrate!(weak_gsc, frame)
    weak_og = compute_optical_gains!(weak_gsc, frame)
    @test all(weak_mode_mask(weak_gsc))
    @test weak_og == ones(2)

    det = Detector(noise=NoiseReadout(1e-3), integration_time=2.0, qe=0.8, psf_sampling=2, binning=4)
    gsc_with_det = GainSensingCamera(mask, basis; detector=det)
    metadata = detector_metadata(gsc_with_det)
    @test metadata isa GSCDetectorMetadata
    @test metadata.integration_time == 2.0
    @test metadata.qe == 0.8
    @test metadata.psf_sampling == 2
    @test metadata.binning == 4
    @test metadata.noise == :readout
    @test metadata.readout_sigma == 1e-3
    @test occursin("psf_sampling=2", sprint(show, MIME"text/plain"(), gsc_with_det))

    detach_detector!(gsc_with_det)
    @test detector_metadata(gsc_with_det) === nothing
    attach_detector!(gsc_with_det, det)
    @test detector_metadata(gsc_with_det) isa GSCDetectorMetadata
end

@testset "LiFT" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(8, 8, 3)
    diversity = zeros(8, 8)
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=true)
    H = lift_interaction_matrix(lift, zeros(3), [1, 2])
    @test size(H) == (64, 2)
    kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    lift_analytic = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=kernel)
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3), [1, 2])
    @test size(H_analytic) == (64, 2)
    psf = compute_psf!(tel, src; zero_padding=1)
    coeffs = reconstruct(lift, psf, [1, 2])
    @test length(coeffs) == 2
    @test AdaptiveOpticsSim.effective_solve_mode(AdaptiveOpticsSim.ScalarCPUStyle(), LiFTSolveAuto()) isa LiFTSolveQR
    diag = diagnostics(lift)
    @test diag.used_qr isa Bool
    @test isfinite(diag.residual_norm)
    @test isfinite(diag.weighted_residual_norm)
    @test isfinite(diag.update_norm)
    @test isfinite(diag.condition_ratio) || isinf(diag.condition_ratio)
    lift_normal = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, solve_mode=LiFTSolveNormalEquations())
    coeffs_normal = reconstruct(lift_normal, psf, [1, 2])
    @test length(coeffs_normal) == 2
    @test all(isfinite, coeffs_normal)
    @test !diagnostics(lift_normal).used_qr
    lift_damped = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTLevenbergMarquardt())
    coeffs_damped = reconstruct(lift_damped, psf, [1, 2])
    @test length(coeffs_damped) == 2
    @test all(isfinite, coeffs_damped)
    @test diagnostics(lift_damped).regularization >= 0
    lift_adaptive = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTAdaptiveLevenbergMarquardt())
    coeffs_adaptive = reconstruct(lift_adaptive, psf, [1, 2])
    @test length(coeffs_adaptive) == 2
    @test all(isfinite, coeffs_adaptive)
    @test diagnostics(lift_adaptive).regularization >= 0
    det_binned = Detector(noise=NoiseNone(), psf_sampling=1, binning=2)
    frame_binned = capture!(det_binned, psf; rng=MersenneTwister(3))
    lift_binned = LiFT(tel, src, basis, det_binned; diversity_opd=diversity, iterations=2,
        img_resolution=size(frame_binned, 1), numerical=true)
    coeffs_binned = reconstruct(lift_binned, frame_binned, [1, 2])
    @test length(coeffs_binned) == 2
    @test all(isfinite, coeffs_binned)
    det_readout = Detector(noise=NoiseReadout(1e-3), psf_sampling=1)
    lift_readout = LiFT(tel, src, basis, det_readout; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true)
    coeffs_readout = reconstruct(lift_readout, psf, [1, 2])
    @test length(coeffs_readout) == 2
    @test all(isfinite, coeffs_readout)

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    lift_sep = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=sep_kernel)
    @test lift_sep.params.object_kernel isa AdaptiveOpticsSim.LiFTSeparableObjectKernel
    dense_conv = similar(psf)
    sep_conv = similar(psf)
    tmp_conv = similar(psf)
    AdaptiveOpticsSim.conv2d_same!(dense_conv, psf, sep_kernel)
    AdaptiveOpticsSim.conv2d_same_separable!(
        sep_conv,
        tmp_conv,
        psf,
        lift_sep.params.object_kernel.row,
        lift_sep.params.object_kernel.col,
    )
    @test isapprox(sep_conv, dense_conv; rtol=1e-6, atol=1e-6)
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    rho = [0.0, 1e-6, 0.1, 1.0]
    cov = phase_covariance(rho, atm)
    @test length(cov) == length(rho)
    @test all(isfinite, cov)
    @test cov[1] >= cov[2]
    @test cov[2] <= cov[1]
    @test cov[3] <= cov[2]
    @test cov[4] <= cov[3]
    @test abs(cov[1] - cov[2]) / cov[1] < 1e-3
    var = phase_variance(atm)
    @test var > 0
    psd = phase_spectrum([0.1], atm)
    @test length(psd) == 1
    @test psd[1] > 0
    screen = ft_phase_screen(atm, 8, 0.1; rng=MersenneTwister(1))
    @test size(screen) == (8, 8)

    for z in (1e-6, 1e-3, 0.1, 1.0, 4.0, 10.0)
        ref = SpecialFunctions.besselk(5 / 6, z)
        approx = AdaptiveOpticsSim._kv56_scalar(z)
        @test isapprox(real(approx), ref; rtol=2e-4, atol=1e-10)
    end
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    basis = modal_basis(dm, tel; n_modes=2)
    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    est = estimate_misregistration(meta, meta.calib0.D; misregistration_zero=Misregistration())
    @test est.shift_x ≈ 0.0
    @test est.shift_y ≈ 0.0
end

@testset "Interface conformance" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    lgs = LGSSource()
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    wfs = ShackHartmann(tel; n_subap=2)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    det = Detector(noise=NoiseNone())
    psf = fill(1.0, 8, 8)

    assert_source_interface(src)
    assert_source_interface(lgs)
    assert_atmosphere_interface(atm, tel)
    assert_wfs_interface(wfs, tel)
    assert_dm_interface(dm, tel)
    assert_detector_interface(det, psf)
end

@testset "Telemetry and config" begin
    telemetry = Telemetry()
    record!(telemetry, 1, 0.0; wfe_rms=1.0, strehl=0.5, loop_gain=0.2)
    @test length(telemetry) == 1
    row = telemetry[1]
    @test row.iter == 1
    @test row.strehl ≈ 0.5

    @test Tables.istable(Telemetry)
    @test Tables.rowaccess(Telemetry)
    rows = collect(Tables.rows(telemetry))
    @test length(rows) == 1
    @test rows[1].loop_gain ≈ 0.2

    closed_loop = ClosedLoopTrace(Float32[
        100 80 0.50 3 4
        120 90 0.45 5 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(closed_loop) == 2
    @test eltype(typeof(closed_loop)) == ClosedLoopTraceRow{Float32}
    @test closed_loop[2].iter == 2
    @test closed_loop[2].t ≈ 0.102f0
    @test closed_loop[2].residual_rms_nm ≈ 90.0f0
    @test Tables.istable(typeof(closed_loop))
    @test Tables.rowaccess(typeof(closed_loop))
    closed_rows = collect(Tables.rows(closed_loop))
    @test length(closed_rows) == 2
    @test closed_rows[1].command_norm ≈ 4.0f0

    gsc_closed_loop = GSCClosedLoopTrace(Float32[
        100 80 0.50 3 0.9 4
        120 90 0.45 5 0.8 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(gsc_closed_loop) == 2
    @test eltype(typeof(gsc_closed_loop)) == GSCClosedLoopTraceRow{Float32}
    @test gsc_closed_loop[2].mean_optical_gain ≈ 0.8f0
    @test Tables.istable(typeof(gsc_closed_loop))
    @test Tables.rowaccess(typeof(gsc_closed_loop))
    gsc_closed_rows = collect(Tables.rows(gsc_closed_loop))
    @test length(gsc_closed_rows) == 2
    @test gsc_closed_rows[1].slope_norm ≈ 3.0f0

    replay = GSCAtmosphereReplayTrace(Float32[
        140 100 110 0.50 0.45 3 0.9
        150 105 115 0.48 0.42 5 0.8
    ]; dt=0.002f0, t0=0.1f0)
    @test length(replay) == 2
    @test eltype(typeof(replay)) == GSCAtmosphereReplayTraceRow{Float32}
    @test replay[2].sci_strehl ≈ 0.42f0
    @test replay[2].slope_norm ≈ 5.0f0
    @test Tables.istable(typeof(replay))
    @test Tables.rowaccess(typeof(replay))
    replay_rows = collect(Tables.rows(replay))
    @test length(replay_rows) == 2
    @test replay_rows[1].ngs_forcing_rms_nm ≈ 140.0f0

    @test_throws DimensionMismatchError ClosedLoopTrace(zeros(2, 4))
    @test_throws DimensionMismatchError GSCClosedLoopTrace(zeros(2, 5))
    @test_throws DimensionMismatchError GSCAtmosphereReplayTrace(zeros(2, 6))

    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    cfg = snapshot_config(tel=tel, src=src)
    @test haskey(cfg, "tel")
    @test haskey(cfg, "src")
    path, io = mktemp()
    close(io)
    write_config_toml(path, cfg)
    parsed = TOML.parsefile(path)
    @test parsed["tel"]["resolution"] == tel.params.resolution
end

@testset "Reference compare conventions" begin
    convention = parse_reference_compare_convention(Dict(
        "convention" => "oopao_geometric_sh_signal_2d",
    ))
    raw = Float64[1, 2, 3, 4]
    adapted = adapt_compare_convention(convention, raw)
    @test adapted ≈ OOPAO_GEOMETRIC_SH_SLOPE_SCALE .* Float64[3, 4, 1, 2]

    legacy = parse_reference_compare_convention(Dict(
        "swap_halves" => true,
        "scale" => OOPAO_GEOMETRIC_SH_SLOPE_SCALE,
    ))
    @test adapt_compare_convention(legacy, raw) == adapted
    @test parse_reference_compare_convention(nothing) isa IdentityCompareConvention
    @test adapt_compare_convention(IdentityCompareConvention(), raw) === raw
    @test_throws InvalidConfiguration parse_reference_compare_convention(Dict(
        "swap_halves" => true,
    ))
end

@testset "Reference storage and detector conventions" begin
    col_data = Float64[1, 2, 3, 4]
    row_data = Float64[1, 2, 3, 4, 5, 6]
    @test reshape_reference_data(col_data, (2, 2), JuliaColumnMajorStorage()) == [1 3; 2 4]
    @test reshape_reference_data(row_data, (2, 3), NumPyRowMajorStorage()) == [1 2 3; 4 5 6]
    @test parse_reference_storage_convention("F") isa JuliaColumnMajorStorage
    @test parse_reference_storage_convention("numpy_row_major") isa NumPyRowMajorStorage
    @test_throws InvalidConfiguration parse_reference_storage_convention("weird")

    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=2, binning=1)
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}()) == 16
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}("img_resolution" => 12)) == 12
end

@testset "Reference harness fixture" begin
    root = mktempdir()
    create_reference_fixture(root)
    bundle = load_reference_bundle(root)
    @test length(bundle.cases) == 6
    for case in bundle.cases
        result = validate_reference_case(case)
        @test size(result.actual) == size(result.expected)
        @test result.ok
    end
end

@testset "OOPAO reference regression" begin
    root = default_reference_root()
    if has_reference_bundle(root)
        bundle = load_reference_bundle(root)
        @test !isempty(bundle.cases)
        for case in bundle.cases
            result = validate_reference_case(case)
            @test size(result.actual) == size(result.expected)
            @test result.ok
        end
    else
        @info "Skipping OOPAO reference regression; no manifest found" root=root
        @test true
    end
end

@testset "Tutorial examples" begin
    image = run_tutorial_example("image_formation.jl")
    @test size(image.psf_nominal) == size(image.psf_aberrated)
    @test maximum(image.psf_nominal) > maximum(image.psf_aberrated)

    detector = run_tutorial_example("detector.jl")
    @test sum(detector.frame_native) ≈ sum(detector.psf)
    @test size(detector.frame_sampled, 1) < size(detector.frame_native, 1)

    asterism = run_tutorial_example("asterism.jl")
    @test size(asterism.per_source_psf, 3) == 4
    @test sum(asterism.combined_psf) >= sum(@view asterism.per_source_psf[:, :, 1])

    spatial = run_tutorial_example("spatial_filter.jl")
    @test size(spatial.filtered_phase) == size(spatial.filtered_amplitude)
    @test all(isfinite, spatial.filtered_phase)

    ncpa = run_tutorial_example("ncpa.jl")
    @test size(ncpa.psf, 1) == 2 * size(ncpa.ncpa_opd, 1)
    @test size(ncpa.psf, 2) == 2 * size(ncpa.ncpa_opd, 2)
    @test maximum(ncpa.psf) > 0

    lift = run_tutorial_example("lift.jl")
    @test length(lift.coeffs_true) == length(lift.coeffs_fit)
    @test all(isfinite, lift.coeffs_fit)

    sprint = run_tutorial_example("sprint.jl")
    @test isfinite(sprint.estimate.shift_x)
    @test isfinite(sprint.estimate.shift_y)

    gsc = run_tutorial_example("gain_sensing_camera.jl")
    @test length(gsc.optical_gains) == 4
    @test all(isfinite, gsc.optical_gains)
    @test sum(gsc.calibration_frame) > 0
    @test sum(gsc.frame) > 0
    @test size(gsc.atmosphere_trace, 2) == 7
    @test all(isfinite, gsc.atmosphere_trace)

    transfer = run_tutorial_example("transfer_function.jl")
    @test size(transfer.rejection_db, 2) == length(transfer.loop_gains)
    @test all(isfinite, transfer.rejection_db)
    @test all(isfinite, transfer.closed_loop_db)

    tomography = run_tutorial_example("tomography.jl")
    @test all(isfinite, tomography.wavefront[.!isnan.(tomography.wavefront)])
    @test length(tomography.commands) == 4
    @test all(isfinite, tomography.commands)

    for name in ("closed_loop_shack_hartmann.jl", "closed_loop_pyramid.jl", "closed_loop_bioedge.jl")
        loop = run_tutorial_example(name)
        @test length(loop.residual_before) == length(loop.residual_after)
        @test all(isfinite, loop.residual_before)
        @test all(isfinite, loop.residual_after)
        @test maximum(loop.final_psf) > 0
    end
end
