function dm_surface_allocations(dm)
    update_surface!(dm)
    return @allocated update_surface!(dm)
end

@testset "Deformable-mirror surface formation" begin
    tel = Telescope(
        resolution=32,
        diameter=8.0,
        central_obstruction=0.0,
    )
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    @test dm.state.modes isa AdaptiveOpticsSim.GaussianInfluenceOperator
    @test Base.summarysize(dm.state.modes) <
        prod(size(dm.state.modes)) * sizeof(eltype(dm.state.modes)) ÷ 10

    dm.state.coefs .= range(-0.2, 0.3; length=length(dm.state.coefs))
    update_surface!(dm)
    dense_reference = reshape(
        Array(dm.state.modes) * Array(dm.state.coefs),
        size(dm.state.opd),
    )
    @test dm.state.opd ≈ dense_reference rtol=5e-14
    if coverage_instrumented()
        @test_skip "DM allocation assertions are disabled under coverage instrumentation"
    else
        @test dm_surface_allocations(dm) == 0
    end

    sampled_topology = SampledActuatorTopology(
        actuator_coordinates(dm)[:, 1:4];
        metadata=(source=:measured,),
    )
    measured_modes = Array(dm.state.modes[:, 1:4])
    measured = DeformableMirror(
        tel;
        topology=sampled_topology,
        influence_model=MeasuredInfluenceFunctions(measured_modes),
        actuator_model=CompositeDMActuatorModel(
            ActuatorHealthMap([1.0, 0.0, 0.5, 1.0]),
            ClippedActuators(-0.2, 0.2),
        ),
    )
    set_command!(measured, [0.5, 0.1, -1.0, 0.1])
    update_surface!(measured)
    @test measured.state.actuator_coefs ≈ [0.2, 0.0, -0.2, 0.1]
    @test topology_metadata(measured) == (source=:measured,)
end

@testset "Independent controllable optics compose additively" begin
    tel = Telescope(
        resolution=24,
        diameter=8.0,
        central_obstruction=0.0,
    )
    tiptilt = TipTiltMirror(tel; scale=1.0)
    focus = FocusStage(tel; scale=0.5)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    @test command_storage(tiptilt) !== command_storage(focus)
    @test command_storage(tiptilt) !== command_storage(dm)

    set_command!(tiptilt, [5e-9, -2e-9])
    set_command!(focus, [3e-9])
    set_command!(dm, range(-1e-9, 1e-9; length=n_control_dofs(dm)))

    individual_opd = map((tiptilt, focus, dm)) do optic
        pupil = PupilFunction(tel)
        update_surface!(optic)
        apply_surface!(pupil, optic, DMReplace())
        copy(opd_map(pupil))
    end

    combined = PupilFunction(tel)
    reset_opd!(combined)
    for optic in (tiptilt, focus, dm)
        update_surface!(optic)
        apply_surface!(combined, optic, DMAdditive())
    end
    @test opd_map(combined) ≈
        individual_opd[1] .+ individual_opd[2] .+ individual_opd[3]

    original_dm_command = copy(command_storage(dm))
    set_command!(tiptilt, zeros(2))
    @test command_storage(dm) == original_dm_command
    @test_throws DimensionMismatchError set_command!(tiptilt, zeros(3))
end

@testset "Modal controllable-optic descriptors" begin
    tel = Telescope(
        resolution=16,
        diameter=8.0,
        central_obstruction=0.0,
        backend=CPUBackend(),
    )
    function_modal = ModalControllableOptic(
        tel,
        FunctionModalBasis(((x, y) -> x, (x, y) -> y));
        labels=:function_modal,
    )
    zernike = ModalControllableOptic(
        tel,
        ZernikeOpticBasis([2, 3]),
    )
    matrix_modal = ModalControllableOptic(
        tel,
        MatrixModalBasis(reshape(collect(1.0:512.0), 16 * 16, 2)),
    )
    cartesian = ModalControllableOptic(
        tel,
        CartesianTiltBasis(; scale=0.1),
    )

    @test function_modal.params.labels === :function_modal
    @test zernike.params.labels == (:zernike_2, :zernike_3)
    @test matrix_modal.params.labels === :modal_optic
    @test cartesian.params.labels == (:x_tilt, :y_tilt)
    for optic in (function_modal, zernike, matrix_modal, cartesian)
        @test backend(optic) isa CPUBackend
        @test command_storage(optic) isa Vector
    end

    set_command!(zernike, [0.01, -0.02])
    pupil = PupilFunction(tel)
    update_surface!(zernike)
    apply_surface!(pupil, zernike, DMAdditive())
    @test norm(opd_map(pupil)) > 0
    @test_throws MethodError set_command!(
        zernike,
        (zernike_2=[0.01], zernike_3=[-0.02]),
    )
end

@testset "Calibration and reusable control primitives" begin
    tel = Telescope(
        resolution=16,
        diameter=8.0,
        central_obstruction=0.0,
    )
    pupil = PupilFunction(tel)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmannWFS(tel; n_lenslets=2)
    imat = interaction_matrix(dm, wfs, pupil; amplitude=0.1)
    @test size(imat.matrix) ==
        (length(slopes(wfs)), length(command_storage(dm)))

    fill!(command_storage(dm), 0.03)
    fill!(opd_map(pupil), 0.02)
    command_before = copy(command_storage(dm))
    opd_before = copy(opd_map(pupil))
    storage = similar(imat.matrix)
    inplace = interaction_matrix!(
        storage,
        dm,
        wfs,
        pupil;
        amplitude=0.1,
    )
    @test inplace.matrix === storage
    @test inplace.matrix ≈ imat.matrix atol=0 rtol=0
    @test command_storage(dm) == command_before
    @test opd_map(pupil) == opd_before

    reconstructor = FactorizedReconstructor(imat; gain=1.0)
    command = similar(command_storage(dm))
    reconstruct!(command, reconstructor, slopes(wfs))
    @test length(command) == n_control_dofs(dm)

    delay = VectorDelayLine(command, 2)
    first_sample = copy(shift_delay!(delay, fill(1.0, length(command))))
    second_sample = copy(shift_delay!(delay, fill(2.0, length(command))))
    third_sample = copy(shift_delay!(delay, fill(3.0, length(command))))
    @test all(iszero, first_sample)
    @test all(iszero, second_sample)
    @test all(==(1.0), third_sample)
    @test_throws DimensionMismatchError shift_delay!(
        delay,
        zeros(length(command) + 1),
    )
    @test_throws InvalidConfiguration VectorDelayLine(command, -1)

    zero_delay = VectorDelayLine(command, 0)
    @test shift_delay!(zero_delay, command) == command

    controller = DiscreteIntegratorController(
        length(command);
        gain=0.1,
        tau=0.02,
    )
    output = update!(controller, command, 0.01)
    @test output === controller_output(controller)
    @test reset_controller!(controller) === controller
    @test all(iszero, controller_output(controller))
    @test_throws DimensionMismatchError update!(
        controller,
        zeros(length(command) + 1),
        0.01,
    )
    @test_throws InvalidConfiguration update!(controller, command, -0.01)

    controlled = ControlledReconstructor(
        reconstructor,
        DiscreteIntegratorController(
            length(command);
            gain=1.0,
            tau=1e-3,
        );
        dt=1e-3,
    )
    reconstruct!(command, controlled, slopes(wfs))
    reconstruct!(command, controlled, slopes(wfs))
    @test all(isfinite, command)
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(reconstruct!(command, controlled, slopes(wfs))) == 0
    end

    timing = runtime_timing(() -> reconstruct!(
        command,
        reconstructor,
        slopes(wfs),
    ); warmup=1, samples=3, gc_before=false)
    @test timing.samples == 3
    @test timing.min_ns <= timing.p50_ns <= timing.max_ns
    @test timing.p50_ns <= timing.p95_ns <= timing.p99_ns
end

@testset "AO188 and AO3k model-specific simulation" begin
    params = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        high_order_samples=4,
        low_order_lenslets=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
    )
    simulation = subaru_ao188_simulation(
        ;
        params=params,
        rng=MersenneTwister(1),
    )
    @test simulation isa AO188Simulation
    @test count(simulation.active_mask) == params.n_active_actuators
    @test size(simulation.high_M2C) ==
        (params.n_act^2, params.n_control_modes)
    @test size(simulation.low_M2C) ==
        (params.n_act^2, params.n_low_order_modes)
    step!(simulation)
    step!(simulation)
    @test maximum(abs, simulation.command) > 0
    products = readout(simulation)
    @test products.command === simulation.command
    @test length(products.signals) == 2
    @test length(products.wfs_frames) == 2
    @test length(products.wfs_metadata) == 2
    timing = runtime_timing(
        simulation;
        warmup=1,
        samples=2,
        gc_before=false,
    )
    @test timing.samples == 2

    prepared_params = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        high_order_samples=4,
        low_order_lenslets=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
        replay_mode=PreparedReplayMode(),
    )
    prepared = subaru_ao188_simulation(
        ;
        params=prepared_params,
        rng=MersenneTwister(2),
    )
    @test prepared.replay_prepared
    step!(prepared)

    curvature_params = AO188CurvatureSimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        high_order_samples=4,
        low_order_lenslets=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
    )
    curvature = subaru_ao188_curvature_simulation(
        ;
        params=curvature_params,
        rng=MersenneTwister(3),
    )
    @test curvature.high_wfs isa CurvatureWFS
    @test curvature.high_detector isa APDDetector
    step!(curvature)
    @test size(readout(curvature).wfs_frames[1]) ==
        (2, curvature_params.high_order_samples^2)

    ao3k_params = AO3kSimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        high_order_samples=4,
        low_order_lenslets=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
    )
    ao3k = subaru_ao3k_simulation(
        ;
        params=ao3k_params,
        rng=MersenneTwister(4),
    )
    @test ao3k.high_wfs isa PyramidWFS
    @test ao3k.high_detector.params.sensor isa
        HgCdTeAvalancheArraySensor
    step!(ao3k)
    @test length(readout(ao3k).command) == ao3k_params.n_act^2
end

mutable struct EnsembleCounter
    value::Int
    owner::Any
end

struct ImmutableEnsembleMember
    owner::Base.RefValue{Int}
end

EnsembleCounter(value::Int=0) = EnsembleCounter(value, Ref(value))

AdaptiveOpticsSim.ensemble_ownership_roots(counter::EnsembleCounter) =
    (counter.owner,)
AdaptiveOpticsSim.ensemble_ownership_roots(member::ImmutableEnsembleMember) =
    (member.owner,)

function increment_counter!(counter::EnsembleCounter)
    counter.value += 1
    return counter
end

@testset "Generic coarse ensembles" begin
    first_counter = EnsembleCounter()
    second_counter = EnsembleCounter(10)
    ensemble = SimulationEnsemble(first_counter, second_counter)
    @test ensemble_members(ensemble) ===
        (first_counter, second_counter)
    @test execution_policy(ensemble) isa SequentialExecution
    @test AdaptiveOpticsSim.run_ensemble!(
        increment_counter!,
        ensemble,
    ) === ensemble
    @test (first_counter.value, second_counter.value) == (1, 11)

    immutable_values = SimulationEnsemble(1, 1)
    @test ensemble_members(immutable_values) === (1, 1)

    threaded = SimulationEnsemble(
        EnsembleCounter(),
        EnsembleCounter();
        policy=ThreadedExecution(),
    )
    AdaptiveOpticsSim.run_ensemble!(increment_counter!, threaded)
    @test all(counter -> counter.value == 1, ensemble_members(threaded))

    shared_owner = Ref(0)
    @test_throws InvalidConfiguration SimulationEnsemble(
        EnsembleCounter(0, shared_owner),
        EnsembleCounter(0, shared_owner);
        policy=ThreadedExecution(),
    )
    @test_throws InvalidConfiguration SimulationEnsemble(
        ImmutableEnsembleMember(shared_owner),
        ImmutableEnsembleMember(shared_owner);
        policy=ThreadedExecution(),
    )
    @test_throws InvalidConfiguration SimulationEnsemble(())

    for policy in (
        BackendStreamExecution(),
        AcceleratedKernelsExecution(),
        DaggerExecution(),
    )
        unsupported = SimulationEnsemble(
            EnsembleCounter();
            policy=policy,
        )
        @test_throws UnsupportedAlgorithm AdaptiveOpticsSim.run_ensemble!(
            increment_counter!,
            unsupported,
        )
    end

    if Threads.nthreads() == 1
        deterministic = SimulationEnsemble(
            EnsembleCounter();
            policy=DeterministicExecution(),
        )
        AdaptiveOpticsSim.run_ensemble!(
            increment_counter!,
            deterministic,
        )
        @test first(ensemble_members(deterministic)).value == 1
        @test BLAS.get_num_threads() == 1
    else
        @test_throws InvalidConfiguration SimulationEnsemble(
            EnsembleCounter();
            policy=DeterministicExecution(),
        )
    end
end
