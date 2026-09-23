using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsCalibration
using BenchmarkTools
using Random

const PR = AdaptiveOpticsCalibration.PhaseRetrieval

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
    src = Source(band=:I, magnitude=0.0)
    rate = shack_hartmann_rate_map(wfs, pupil, src)
    optics_plan = prepare_wfs_optics(shack_hartmann_optics(wfs, src),
        pupil, rate)
    return @benchmark form_wfs_optical_products!($rate, $pupil, $optics_plan)
end

function bench_wfs_lgs()
    tel = Telescope(resolution=48, diameter=8.0, central_obstruction=0.0)
    wfs = ShackHartmannWFS(tel; n_lenslets=6)
    lgs = LGSSource(elongation_factor=1.3, photon_irradiance=1.0)
    pupil = PupilFunction(tel)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i - j
    end
    rate = shack_hartmann_rate_map(wfs, pupil, lgs)
    optics_plan = prepare_wfs_optics(shack_hartmann_optics(wfs, lgs),
        pupil, rate)
    return @benchmark form_wfs_optical_products!($rate, $pupil, $optics_plan)
end

function bench_pyramid()
    tel = Telescope(resolution=48, diameter=8.0, central_obstruction=0.0)
    wfs = PyramidWFS(tel; pupil_samples=6, modulation=3.0,
        modulation_points=4)
    src = Source(band=:I, magnitude=0.0)
    pupil = PupilFunction(tel)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        pupil.opd[i, j] = i + j
    end
    front_end = PyramidOpticalFrontEnd(wfs, src)
    rate = pyramid_rate_map(front_end, pupil)
    optics_plan = prepare_wfs_optics(front_end, pupil, rate)
    return @benchmark form_wfs_optical_products!(
        $rate, $pupil, $optics_plan)
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
    specification = PR.LiFTSpecification(forward, observation)
    method = PR.LiFT(iterations=1, mode_indices=1:3,
        jacobian_method=numerical ? PR.LiFTNumericalJacobian() :
            PR.LiFTAnalyticJacobian())
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = PR.LiFTInputs(observation.values)
    return result, workspace, plan, inputs
end

function bench_lift(numerical::Bool)
    result, workspace, plan, inputs = prepare_lift_benchmark(numerical)
    return @benchmark AdaptiveOpticsCalibration.process!(
        $result, $workspace, $plan, $inputs)
end

function alloc_checks()
    result_a, workspace_a, plan_a, inputs_a = prepare_lift_benchmark(false)
    result_n, workspace_n, plan_n, inputs_n = prepare_lift_benchmark(true)
    AdaptiveOpticsCalibration.process!(result_a, workspace_a, plan_a, inputs_a)
    AdaptiveOpticsCalibration.process!(result_n, workspace_n, plan_n, inputs_n)
    alloc_lift_a = @allocated AdaptiveOpticsCalibration.process!(
        result_a, workspace_a, plan_a, inputs_a)
    alloc_lift_n = @allocated AdaptiveOpticsCalibration.process!(
        result_n, workspace_n, plan_n, inputs_n)

    println("Allocation checks:")
    println("  LiFT analytic (in-place): $(alloc_lift_a) bytes")
    println("  LiFT numerical (in-place): $(alloc_lift_n) bytes")
end

println("Direct-imaging benchmark:")
display(bench_direct_imaging())

println("WFS benchmark:")
display(bench_wfs())

println("WFS LGS benchmark:")
display(bench_wfs_lgs())

println("Pyramid benchmark:")
display(bench_pyramid())

println("AOC LiFT analytic complete-call benchmark:")
display(bench_lift(false))

println("AOC LiFT numerical complete-call benchmark:")
display(bench_lift(true))

alloc_checks()
