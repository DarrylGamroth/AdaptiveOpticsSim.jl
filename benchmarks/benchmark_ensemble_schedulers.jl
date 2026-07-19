using AcceleratedKernels
using AdaptiveOpticsSim
using Dagger
using LinearAlgebra
using Random
using Statistics

BLAS.set_num_threads(1)
AdaptiveOpticsSim.set_fft_provider_threads!(1)

function make_runtime(seed::Integer, resolution::Int)
    tel = Telescope(
        resolution=resolution,
        diameter=8.0,
        central_obstruction=0.0,
    )
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    n_act = max(2, resolution ÷ 4)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.3)
    wfs = ShackHartmannWFS(tel; n_lenslets=n_act)
    simulation = AOSimulation(tel, src, atm, dm, wfs)

    response = zeros(
        eltype(pupil_reflectivity(tel)),
        length(slopes(wfs)),
        length(dm.state.coefs),
    )
    @inbounds for i in 1:min(size(response)...)
        response[i, i] = one(eltype(response))
    end
    reconstructor = ModalReconstructor(InteractionMatrix(response, 0.1);
        gain=0.5)
    runtime = AdaptiveOpticsSim.ClosedLoopRuntime(
        simulation,
        reconstructor;
        atmosphere_step=1e-3,
        rng=MersenneTwister(seed),
        profile=HILRuntimeProfile(),
    )
    prepare!(runtime)
    return runtime
end

function make_ensemble(policy, members::Int, resolution::Int)
    runtimes = ntuple(i -> make_runtime(0x5eed + i, resolution), members)
    return SimulationEnsemble(runtimes; policy=policy)
end

function measure_scheduler(policy, members::Int, resolution::Int;
    warmup::Int, samples::Int)
    ensemble = make_ensemble(policy, members, resolution)
    for _ in 1:warmup
        step!(ensemble)
    end
    GC.gc()
    latencies_ns = Vector{UInt64}(undef, samples)
    @inbounds for i in eachindex(latencies_ns)
        start_ns = time_ns()
        step!(ensemble)
        latencies_ns[i] = time_ns() - start_ns
    end
    allocated_bytes = @allocated step!(ensemble)
    values = Float64.(latencies_ns)
    return (
        policy=string(nameof(typeof(policy))),
        minimum_ns=minimum(values),
        median_ns=median(values),
        p95_ns=quantile(values, 0.95),
        p99_ns=quantile(values, 0.99),
        maximum_ns=maximum(values),
        allocated_bytes=allocated_bytes,
        members_per_second=members * 1e9 / median(values),
    )
end

function print_result(result)
    println(
        result.policy,
        ": min=", round(Int, result.minimum_ns), " ns",
        " median=", round(Int, result.median_ns), " ns",
        " p95=", round(Int, result.p95_ns), " ns",
        " p99=", round(Int, result.p99_ns), " ns",
        " max=", round(Int, result.maximum_ns), " ns",
        " alloc=", result.allocated_bytes, " B",
        " throughput=", round(result.members_per_second; digits=1),
        " member-steps/s",
    )
end

function main()
    members = parse(Int, get(ENV, "AOS_ENSEMBLE_MEMBERS", "8"))
    resolution = parse(Int, get(ENV, "AOS_ENSEMBLE_RESOLUTION", "16"))
    warmup = parse(Int, get(ENV, "AOS_ENSEMBLE_WARMUP", "10"))
    samples = parse(Int, get(ENV, "AOS_ENSEMBLE_SAMPLES", "100"))
    samples >= 100 || error("AOS_ENSEMBLE_SAMPLES must be at least 100 for p99")

    println("AdaptiveOpticsSim ensemble scheduler benchmark")
    println("julia=", VERSION,
        " threads=", Threads.nthreads(),
        " blas_threads=", BLAS.get_num_threads(),
        " cpu_threads=", Sys.CPU_THREADS,
        " members=", members,
        " resolution=", resolution,
        " warmup=", warmup,
        " samples=", samples)
    println("AK=", Base.pkgversion(AcceleratedKernels),
        " Dagger=", Base.pkgversion(Dagger))
    println("Measurement is a warmed closed-loop offline throughput comparison; it is not a HIL latency contract.")

    policies = (
        SequentialExecution(),
        ThreadedExecution(),
        AcceleratedKernelsExecution(
            max_tasks=Threads.nthreads(),
            min_members_per_task=1,
        ),
        DaggerExecution(scope=Dagger.scope(worker=1)),
    )
    for policy in policies
        print_result(measure_scheduler(
            policy,
            members,
            resolution;
            warmup=warmup,
            samples=samples,
        ))
    end
end

main()
