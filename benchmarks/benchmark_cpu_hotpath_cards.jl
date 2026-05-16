const ROOT_ENV = normpath(joinpath(@__DIR__, ".."))
const BENCHMARK_ENV = @__DIR__
ROOT_ENV in LOAD_PATH || pushfirst!(LOAD_PATH, ROOT_ENV)
BENCHMARK_ENV in LOAD_PATH || push!(LOAD_PATH, BENCHMARK_ENV)

using AdaptiveOpticsSim
using BenchmarkTools

const CPU_HOTPATH_CARDS = (
    ("CPU-PERF-01", "extended-source cached asterism", :extended_source_cache),
    ("CPU-PERF-02", "ForwardDiff Gaussian influence Jacobian workspace", :gaussian_ad),
    ("CPU-PERF-03", "Shack-Hartmann reference subtraction", :sh_reference),
    ("CPU-PERF-04", "subaperture valid-index reuse", :subaperture_layout),
    ("CPU-PERF-05", "VectorDelayLine ring buffer", :delay_line),
    ("CPU-PERF-06", "composite selected apply without Set allocation", :composite_apply),
)

function configure_cpu_hotpath_benchmarks!()
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

function gaussian_ad_probe()
    mis = Misregistration(shift_x=0.03, shift_y=-0.02, rotation_deg=1.5,
        radial_scaling=1.01, tangential_scaling=0.98)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3,
        central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.35, misregistration=mis)
    fields = (:influence_width, :actuator_x, :actuator_y, :shift_x, :shift_y,
        :rotation_deg, :radial_scaling, :tangential_scaling)
    actuator_index = 6
    return () -> AdaptiveOpticsSim.gaussian_dm_influence_parameter_jacobian(
        tel, dm, actuator_index; field_order=fields)
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
    support = Float64.(tel.state.pupil)
    geometry_probe = () -> AdaptiveOpticsSim.update_subaperture_layout!(
        layout, tel.state.pupil, geometry_policy)
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
        ("CPU-PERF-02", "gaussian_ad", gaussian_ad_probe()),
        ("CPU-PERF-03", "sh_reference", sh_reference_probe()),
        ("CPU-PERF-04a", "subaperture_geometry", geometry_probe),
        ("CPU-PERF-04b", "subaperture_flux", flux_probe),
        ("CPU-PERF-05", "delay_line", delay_line_probe()),
        ("CPU-PERF-06", "composite_apply", composite_apply_probe()),
    )
    results = map(probe -> run_probe(probe...), probes)
    println("recommendation")
    println("  Keep the benchmarked CPU SPRINT AD path as default.")
    println("  Keep the allocating Gaussian influence AD probe explicit unless accuracy justifies it.")
    return results
end

run_cpu_hotpath_card_benchmarks()
