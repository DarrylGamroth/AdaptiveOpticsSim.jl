using AdaptiveOpticsSim
using BenchmarkTools
using LinearAlgebra

const CPU_HOTPATH_CARDS = (
    ("CPU-PERF-01", "extended-source cached asterism", :extended_source_cache),
    ("CPU-PERF-02", "Shack-Hartmann reference subtraction", :sh_reference),
    ("CPU-PERF-03", "subaperture valid-index reuse", :subaperture_layout),
    ("CPU-PERF-04", "VectorDelayLine ring buffer", :delay_line),
    ("CPU-PERF-05", "composite selected apply without Set allocation", :composite_apply),
    ("CPU-PERF-06", "SAPHIRA sampled frame response", :sampled_frame_response),
    ("CPU-PERF-07", "batched SAPHIRA sampled frame response", :batched_sampled_frame_response),
    ("CPU-PERF-08", "detector frame binning", :detector_binning),
    ("CPU-PERF-09", "separable Gaussian frame response", :gaussian_frame_response),
    ("CPU-PERF-10", "batched EMCCD capture", :batched_emccd_capture),
    ("CPU-PERF-11", "lazy Gaussian DM operator application", :gaussian_dm_operator),
    ("CPU-PERF-12", "shared multi-arm optical runtime", :shared_optical_runtime),
)

function configure_cpu_hotpath_benchmarks!()
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
    BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
    BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
    BenchmarkTools.DEFAULT_PARAMETERS.gcsample = false
    return nothing
end

function extended_source_cache_probe()
    src = Source(band=:I, magnitude=0.0)
    model = SampledImageSourceModel(
        [0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0],
        pixel_scale_arcsec=0.2,
    )
    ext = with_extended_source(src, model)
    AdaptiveOpticsSim._cached_extended_source_asterism(ext)
    return () -> AdaptiveOpticsSim._cached_extended_source_asterism(ext)
end

function sh_reference_probe()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    wfs.state.reference_signal_2d .= 0.25
    wfs.state.slopes_units = 2.0
    return function ()
        fill!(wfs.state.slopes, 1.0)
        AdaptiveOpticsSim.subtract_reference_and_scale!(wfs)
        return wfs.state.slopes
    end
end

function subaperture_layout_probes()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    layout = subaperture_layout(wfs)
    geometry_policy = AdaptiveOpticsSim.GeometryValidSubapertures(threshold=0.1)
    flux_policy = FluxThresholdValidSubapertures(light_ratio=0.5)
    support = Float64.(pupil_mask(tel))
    geometry_probe = () -> AdaptiveOpticsSim.update_subaperture_layout!(
        layout, pupil_mask(tel), geometry_policy)
    flux_probe = () -> AdaptiveOpticsSim.update_subaperture_layout!(
        layout, support, flux_policy)
    return geometry_probe, flux_probe
end

function delay_line_probe()
    ref = zeros(Float64, 32)
    line = VectorDelayLine(ref, 4)
    sample = collect(range(0.0, 1.0; length=length(ref)))
    return () -> shift_delay!(line, sample)
end

function composite_apply_probe()
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    optic = CompositeControllableOptic(
        :tiptilt => TipTiltMirror(tel; scale=1.0),
        :dm => DeformableMirror(tel; n_act=4, influence_width=0.3),
    )
    command = vcat([1e-8, -2e-8], fill(1e-8, 16))
    set_command!(optic, command)
    return () -> AdaptiveOpticsSim._apply_selected!(optic, tel, DMReplace(), (:tiptilt,))
end

function sampled_frame_response_probe()
    model = SampledFrameResponse(Float32[
        0.00 0.01 0.00
        0.01 0.96 0.01
        0.00 0.01 0.00
    ]; T=Float32)
    frame = rand(Float32, 96, 96)
    scratch = similar(frame)
    return () -> AdaptiveOpticsSim.apply_response!(
        AdaptiveOpticsSim.ScalarCPUStyle(), model, frame, scratch)
end

function batched_sampled_frame_response_probe()
    model = SampledFrameResponse(Float32[
        0.00 0.01 0.00
        0.01 0.96 0.01
        0.00 0.01 0.00
    ]; T=Float32)
    cube = rand(Float32, 256, 4, 4)
    scratch = similar(cube)
    return () -> AdaptiveOpticsSim._batched_apply_response!(
        AdaptiveOpticsSim.ScalarCPUStyle(), model, cube, scratch)
end

function detector_binning_probe()
    input = rand(Float32, 512, 512)
    output = similar(input, 256, 256)
    return () -> AdaptiveOpticsSim._bin2d!(
        AdaptiveOpticsSim.ScalarCPUStyle(), output, input, 2)
end

function gaussian_frame_response_probe()
    model = GaussianPixelResponse(response_width_px=1.5, T=Float32)
    frame = rand(Float32, 96, 96)
    scratch = similar(frame)
    return () -> AdaptiveOpticsSim.apply_response!(
        AdaptiveOpticsSim.ScalarCPUStyle(), model, frame, scratch)
end

function batched_emccd_capture_probe()
    detector = Detector(
        noise=NoiseNone(),
        sensor=EMCCDSensor(excess_noise_factor=1.4, T=Float32),
        gain=20,
        T=Float32,
    )
    cube = fill(Float32(10), 256, 4, 4)
    scratch = similar(cube)
    rng = runtime_rng(17)
    return function ()
        fill!(cube, Float32(10))
        return AdaptiveOpticsSim.capture_stack!(detector, cube, scratch, rng)
    end
end

function gaussian_dm_operator_probe()
    tel = Telescope(resolution=128, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=32, influence_width=0.3)
    dm.state.coefs .= range(-1e-7, 1e-7; length=length(dm.state.coefs))
    return () -> AdaptiveOpticsSim.apply_opd!(dm, tel)
end

function shared_optical_runtime_probe()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    guide = Source(band=:I, magnitude=0.0)
    science = Source(band=:K, magnitude=1.0, coordinates=(4.0, 90.0))
    atmosphere = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmannWFS(tel; n_lenslets=4)
    simulation = AOSimulation(tel, guide, atmosphere, dm, wfs)
    reconstructor = ModalReconstructor(
        interaction_matrix(dm, wfs, tel; amplitude=0.1);
        gain=0.5,
    )
    primary = AdaptiveOpticsSim.ClosedLoopRuntime(simulation, reconstructor;
        rng=runtime_rng(17))
    arm = SharedOpticalArm(
        :science,
        science;
        wfs_channels=OpticalWFSChannel(
            ShackHartmannWFS(tel; n_lenslets=4)),
        science_detectors=(
            Detector(noise=NoiseNone()),
            Detector(noise=NoiseNone()),
        ),
        science_zero_padding=1,
    )
    runtime = SharedOpticalRuntime(primary, arm)
    prepare!(runtime)
    sense!(runtime)
    return () -> sense!(runtime)
end

function run_probe(card_id::AbstractString, label::AbstractString, f)
    f()
    alloc = @allocated f()
    trial = @benchmark $f()
    median_ns = BenchmarkTools.median(trial).time
    println(card_id, " ", label)
    println("  allocations_after_warmup_bytes: ", alloc)
    println("  median_ns: ", median_ns)
    println("  benchmark_memory_bytes: ", BenchmarkTools.memory(trial))
    println("  benchmark_allocs: ", BenchmarkTools.allocs(trial))
    return (; card_id, label, alloc, median_ns, trial)
end

function run_cpu_hotpath_card_benchmarks()
    configure_cpu_hotpath_benchmarks!()
    geometry_probe, flux_probe = subaperture_layout_probes()
    probes = (
        ("CPU-PERF-01", "extended_source_cache", extended_source_cache_probe()),
        ("CPU-PERF-02", "sh_reference", sh_reference_probe()),
        ("CPU-PERF-03a", "subaperture_geometry", geometry_probe),
        ("CPU-PERF-03b", "subaperture_flux", flux_probe),
        ("CPU-PERF-04", "delay_line", delay_line_probe()),
        ("CPU-PERF-05", "composite_apply", composite_apply_probe()),
        ("CPU-PERF-06", "sampled_frame_response", sampled_frame_response_probe()),
        ("CPU-PERF-07", "batched_sampled_frame_response", batched_sampled_frame_response_probe()),
        ("CPU-PERF-08", "detector_binning", detector_binning_probe()),
        ("CPU-PERF-09", "gaussian_frame_response", gaussian_frame_response_probe()),
        ("CPU-PERF-10", "batched_emccd_capture", batched_emccd_capture_probe()),
        ("CPU-PERF-11", "gaussian_dm_operator", gaussian_dm_operator_probe()),
        ("CPU-PERF-12", "shared_optical_runtime", shared_optical_runtime_probe()),
    )
    results = map(probe -> run_probe(probe...), probes)
    println("recommendation")
    println("  Keep the benchmarked CPU SPRINT AD path as default.")
    println("  The standalone Gaussian influence AD probe was removed because it allocated without a demonstrated benefit.")
    return results
end

run_cpu_hotpath_card_benchmarks()
