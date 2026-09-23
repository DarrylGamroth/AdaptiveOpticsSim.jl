using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsCalibration
using KernelAbstractions
using LinearAlgebra
using Statistics

const PR = AdaptiveOpticsCalibration.PhaseRetrieval

function lift_profile_problem(array_backend, selector, resolution)
    T = Float32
    cpu_telescope = Telescope(resolution=resolution, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    device_telescope = Telescope(resolution=resolution, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)
    source = Source(band=:I, magnitude=zero(T), T=T)
    zernike = ZernikeBasis(cpu_telescope, 5; T=T)
    compute_zernike!(zernike, cpu_telescope)
    basis = copy(@view zernike.modes[:, :, 2:4])
    diversity = T(30e-9) .* @view(zernike.modes[:, :, 5])
    truth = T[8e-9, -4e-9, 0]
    opd = copy(diversity)
    @inbounds for mode in eachindex(truth)
        @views @. opd += truth[mode] * basis[:, :, mode]
    end
    cpu_forward = prepare_lift_forward_model(cpu_telescope, source, basis, opd;
        diversity_opd=copy(diversity), focal_resolution=resolution)
    device_forward = prepare_lift_forward_model(device_telescope, source,
        array_backend(basis), array_backend(opd);
        diversity_opd=array_backend(copy(diversity)),
        focal_resolution=resolution)
    observation = copy(intensity_values(evaluate_lift_forward!(cpu_forward)))
    return cpu_forward, device_forward, observation, truth
end

function lift_profile_plan(forward, execution_backend, iterations)
    specification = PR.LiFTSpecification(LiFTForwardModel(forward),
        PR.LiFTPhotonRate())
    method = PR.LiFT(iterations=iterations,
        jacobian_method=PR.LiFTAnalyticJacobian(),
        solve_mode=PR.LiFTSolveNormalEquations(),
        damping=PR.LiFTAdaptiveLevenbergMarquardt(),
        mode_indices=(1, 2),
        model_scaling=PR.LiFTPhysicalRatePreservation(),
        check_convergence=false)
    return execution_backend === nothing ?
        AdaptiveOpticsCalibration.prepare(method, specification) :
        AdaptiveOpticsCalibration.prepare(method,
            AdaptiveOpticsCalibration.KernelExecution(
                specification, execution_backend;
                workgroup_size=128))
end

function quantile_sorted(samples, fraction)
    return samples[clamp(ceil(Int, fraction * length(samples)), 1,
        length(samples))]
end

@noinline function complete_lift!(result, workspace, plan, inputs, synchronize)
    AdaptiveOpticsCalibration.process!(result, workspace, plan, inputs)
    synchronize()
    return result
end

function profile_lift_adapter(label, array_backend, selector, synchronize,
    device_allocated_bytes, used_device_memory, profile_region)
    resolution = parse(Int, get(ENV, "AOS_LIFT_RESOLUTION", "16"))
    iterations = parse(Int, get(ENV, "AOS_LIFT_ITERATIONS", "2"))
    warmups = parse(Int, get(ENV, "AOS_LIFT_WARMUPS", "10"))
    samples = parse(Int, get(ENV, "AOS_LIFT_SAMPLES", "50"))
    profile_calls = parse(Int, get(ENV, "AOS_LIFT_PROFILE_CALLS", "20"))
    cpu_forward, device_forward, observation, truth =
        lift_profile_problem(array_backend, selector, resolution)
    cpu_plan = lift_profile_plan(cpu_forward, nothing, iterations)
    cpu_result = AdaptiveOpticsCalibration.allocate_result(cpu_plan)
    cpu_workspace = AdaptiveOpticsCalibration.allocate_workspace(cpu_plan)
    cpu_inputs = PR.LiFTInputs(copy(observation))
    complete_lift!(cpu_result, cpu_workspace, cpu_plan, cpu_inputs, () -> nothing)
    cpu_coefficients = copy(Array(PR.lift_coefficients(cpu_result)))
    default_fixture = resolution == 16 && iterations == 2
    truth_tolerance = 2f-10
    if default_fixture
        cpu_truth_error = maximum(abs, cpu_coefficients - truth[1:2])
        cpu_truth_error < truth_tolerance || error(
            "AOS/AOC LiFT CPU truth check failed: $(cpu_truth_error)")
    end

    device_observation = array_backend(observation)
    device_observation_snapshot = Array(device_observation)
    execution_backend = label == "cpu" ? nothing :
        KernelAbstractions.get_backend(device_observation)
    plan = lift_profile_plan(device_forward, execution_backend, iterations)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = PR.LiFTInputs(device_observation)
    for _ in 1:warmups
        complete_lift!(result, workspace, plan, inputs, synchronize)
    end
    coefficients = Array(PR.lift_coefficients(result))
    maximum_cpu_difference = maximum(abs, coefficients - cpu_coefficients)
    maximum_cpu_difference < 2f-10 || error(
        "AOS/AOC LiFT CPU parity failed: $(maximum_cpu_difference)")
    maximum_truth_error = maximum(abs, coefficients - truth[1:2])
    if default_fixture
        maximum_truth_error < truth_tolerance || error(
            "AOS/AOC LiFT truth check failed: $(maximum_truth_error)")
    end
    Array(inputs.observation) == device_observation_snapshot || error(
        "AOS/AOC LiFT changed its caller-owned observation")

    host_bytes = @allocated complete_lift!(result, workspace, plan, inputs,
        synchronize)
    accelerator_bytes = device_allocated_bytes() do
        complete_lift!(result, workspace, plan, inputs, synchronize)
    end
    device_memory_before = used_device_memory()
    complete_lift!(result, workspace, plan, inputs, synchronize)
    device_memory_after = used_device_memory()
    times = Vector{UInt64}(undef, samples)
    for i in eachindex(times)
        start = time_ns()
        complete_lift!(result, workspace, plan, inputs, synchronize)
        times[i] = time_ns() - start
    end
    sort!(times)
    if get(ENV, "AOS_LIFT_PROFILE", "0") == "1"
        profile_region() do
            for _ in 1:profile_calls
                complete_lift!(result, workspace, plan, inputs, synchronize)
            end
        end
    end
    post_profile_difference = maximum(abs,
        Array(PR.lift_coefficients(result)) - cpu_coefficients)
    post_profile_difference < 2f-10 || error(
        "AOS/AOC LiFT post-profile CPU parity failed: $(post_profile_difference)")
    post_profile_truth_error = maximum(abs,
        Array(PR.lift_coefficients(result)) - truth[1:2])
    if default_fixture
        post_profile_truth_error < truth_tolerance || error(
            "AOS/AOC LiFT post-profile truth check failed: $(post_profile_truth_error)")
    end
    Array(inputs.observation) == device_observation_snapshot || error(
        "AOS/AOC LiFT changed its observation after profiling")

    println("backend = ", repr(label))
    println("resolution = ", resolution)
    println("iterations = ", iterations)
    println("warmups = ", warmups)
    println("samples = ", samples)
    println("profile_calls = ", profile_calls)
    println("julia_threads = ", Threads.nthreads())
    println("maximum_cpu_difference = ", maximum_cpu_difference)
    println("maximum_truth_error = ", maximum_truth_error)
    println("post_profile_difference = ", post_profile_difference)
    println("post_profile_truth_error = ", post_profile_truth_error)
    println("host_allocated_bytes = ", host_bytes)
    println("device_allocated_bytes = ", accelerator_bytes)
    println("used_device_memory_before = ", device_memory_before)
    println("used_device_memory_after = ", device_memory_after)
    println("p50_ns = ", quantile_sorted(times, 0.5))
    println("p99_ns = ", quantile_sorted(times, 0.99))
    return nothing
end
