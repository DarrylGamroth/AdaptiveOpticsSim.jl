using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control
using BenchmarkTools
using Random

include(joinpath(@__DIR__, "support", "closed_loop_workload.jl"))

function bench_direct_imaging()
    tel = Telescope(resolution=64, diameter=8.0, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    prepared = prepare_direct_imaging(PupilFunction(tel), src; zero_padding=2)
    return @benchmark form_direct_image!($prepared)
end

function bench_wfs()
    tel = Telescope(resolution=64, diameter=8.0, central_obstruction=0.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=8)
    pupil = PupilFunction(tel)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i
    end
    return @benchmark measure!($wfs, $pupil)
end

function bench_wfs_lgs()
    tel = Telescope(resolution=48, diameter=8.0, central_obstruction=0.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive())
    lgs = LGSSource(elongation_factor=1.3, photon_irradiance=1.0)
    pupil = PupilFunction(tel)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i - j
    end
    return @benchmark measure!($wfs, $pupil, $lgs)
end

function bench_pyramid()
    tel = Telescope(resolution=48, diameter=8.0, central_obstruction=0.0)
    wfs = PyramidWFS(tel; pupil_samples=6, modulation=3.0, modulation_points=4, mode=Diffractive())
    src = Source(band=:I, magnitude=0.0)
    pupil = PupilFunction(tel)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i + j
    end
    return @benchmark measure!($wfs, $pupil, $src)
end

function bench_reconstructor()
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4)
    wfs = ShackHartmannWFS(tel; n_lenslets=4)
    imat = interaction_matrix(dm, wfs, PupilFunction(tel); amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = AdaptiveOpticsSim.WavefrontSensors.slopes(wfs)
    return @benchmark reconstruct($recon, $slopes)
end

function bench_reconstructor_inplace()
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4)
    wfs = ShackHartmannWFS(tel; n_lenslets=4)
    imat = interaction_matrix(dm, wfs, PupilFunction(tel); amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = AdaptiveOpticsSim.WavefrontSensors.slopes(wfs)
    out = similar(slopes, size(recon.reconstructor, 1))
    return @benchmark reconstruct!($out, $recon, $slopes)
end

function bench_closed_loop_workload()
    workload = prepare_closed_loop_workload(T=Float64, seed=0)
    return @benchmark step_closed_loop_workload!($workload)
end

function bench_closed_loop_workload_timing()
    workload = prepare_closed_loop_workload(T=Float64, seed=0)
    return runtime_timing(
        () -> step_closed_loop_workload!(workload);
        warmup=10,
        samples=200,
        gc_before=false,
    )
end

function prepare_lift_benchmark(numerical::Bool)
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    basis = rand(16, 16, 6)
    diversity = zeros(16, 16)
    forward = prepare_lift_forward_model(tel, src, basis, diversity;
        diversity_opd=diversity, focal_resolution=16)
    rate = copy(intensity_values(evaluate_lift_forward!(forward)))
    observation = LiFTObservation(forward, rate)
    product = zeros(3)
    definition = LiFT(iterations=1, mode_ids=1:3,
        jacobian_method=numerical ? LiFTNumericalJacobian() :
            LiFTAnalyticJacobian())
    lift = prepare_lift_estimator(definition, forward, observation, product)
    return lift, zeros(6)
end

function bench_lift(numerical::Bool)
    lift, coeffs = prepare_lift_benchmark(numerical)
    return @benchmark AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix(
        $lift, $coeffs)
end

function bench_lift_inplace(numerical::Bool)
    lift, coeffs = prepare_lift_benchmark(numerical)
    H = lift.workspace.H_buffer
    return @benchmark AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix!(
        $H, $lift, $coeffs)
end

function alloc_checks()
    lift_a, coeffs = prepare_lift_benchmark(false)
    lift_n, _ = prepare_lift_benchmark(true)
    H_a = lift_a.workspace.H_buffer
    H_n = lift_n.workspace.H_buffer
    AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix!(
        H_a, lift_a, coeffs)
    AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix!(
        H_n, lift_n, coeffs)
    alloc_lift_a = @allocated AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix!(
        H_a, lift_a, coeffs)
    alloc_lift_n = @allocated AdaptiveOpticsSim.WavefrontSensors.lift_interaction_matrix!(
        H_n, lift_n, coeffs)

    tel_r = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel_r; n_act=4)
    wfs = ShackHartmannWFS(tel_r; n_lenslets=4)
    imat = interaction_matrix(dm, wfs, PupilFunction(tel_r); amplitude=0.1)
    recon = ModalReconstructor(imat; gain=1.0)
    slopes = AdaptiveOpticsSim.WavefrontSensors.slopes(wfs)
    out = similar(slopes, size(recon.reconstructor, 1))
    reconstruct!(out, recon, slopes)
    alloc_recon = @allocated reconstruct!(out, recon, slopes)

    workload = prepare_closed_loop_workload(T=Float64, seed=0)
    alloc_closed_loop =
        @allocated step_closed_loop_workload!(workload)

    println("Allocation checks:")
    println("  LiFT analytic (in-place): $(alloc_lift_a) bytes")
    println("  LiFT numerical (in-place): $(alloc_lift_n) bytes")
    println("  Reconstructor! (in-place): $(alloc_recon) bytes")
    println("  Closed-loop workload: $(alloc_closed_loop) bytes")
end

println("Direct-imaging benchmark:")
display(bench_direct_imaging())

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

println("Closed-loop workload benchmark:")
display(bench_closed_loop_workload())

println("Closed-loop workload timing:")
display(bench_closed_loop_workload_timing())

println("LiFT analytic benchmark:")
display(bench_lift(false))

println("LiFT numerical benchmark:")
display(bench_lift(true))

println("LiFT analytic in-place benchmark:")
display(bench_lift_inplace(false))

println("LiFT numerical in-place benchmark:")
display(bench_lift_inplace(true))

alloc_checks()
