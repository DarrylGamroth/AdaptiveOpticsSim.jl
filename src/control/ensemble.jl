"""
    SimulationEnsemble(members...; policy=SequentialExecution())

Own independent simulation boundaries and execute one coarse operation over
them. The scalar CPU HIL path does not use this scheduler; it continues to call
`step!` directly. Ensemble construction rejects shared mutable plant roots so
threaded or distributed policies cannot silently race on the same simulation.
"""
mutable struct SimulationEnsemble{
    M<:Tuple,
    P<:AbstractExecutionPolicy,
    S,
} <: AbstractControlSimulation
    members::M
    policy::P
    scheduler_state::S
end

@inline ensemble_members(ensemble::SimulationEnsemble) = ensemble.members
@inline execution_policy(ensemble::SimulationEnsemble) = ensemble.policy

@inline _validate_ensemble_member_types(::Tuple{}) = nothing
@inline function _validate_ensemble_member_types(
    members::Tuple{<:AbstractControlSimulation,Vararg})
    _validate_ensemble_member_types(Base.tail(members))
    return nothing
end
_validate_ensemble_member_types(::Tuple) =
    throw(InvalidConfiguration("ensemble members must implement AbstractControlSimulation"))

"""
    ensemble_ownership_roots(simulation)

Return mutable state owners that must not be shared by independently scheduled
ensemble members. Custom simulation boundaries should extend this trait when
they wrap another runtime or plant.
"""
@inline ensemble_ownership_roots(sim::AbstractControlSimulation) = (sim,)

@inline ensemble_source_ownership_roots(::AbstractSource) = ()
@inline ensemble_source_ownership_roots(source::ExtendedSource) = (source,)
@inline ensemble_source_ownership_roots(source::Asterism) =
    _tuple_ownership_roots(Tuple(source.sources))

@inline function ensemble_ownership_roots(sim::AOSimulation)
    return (
        sim.tel,
        sim.atm,
        sim.optic,
        sim.wfs,
        ensemble_source_ownership_roots(sim.src)...,
        ensemble_source_ownership_roots(sim.science_src)...,
    )
end

@inline function ensemble_ownership_roots(runtime::ClosedLoopRuntime)
    return (
        ensemble_ownership_roots(runtime.simulation)...,
        runtime.wfs_detector,
        runtime.science_detector,
        runtime.rng,
        runtime_reconstructor_ownership_roots(runtime.reconstructor)...,
    )
end

@inline ensemble_ownership_roots(interface::SimulationInterface) =
    ensemble_ownership_roots(interface.runtime)

@inline ensemble_ownership_roots(source::AbstractSource) =
    ensemble_source_ownership_roots(source)

@inline _tuple_ownership_roots(::Tuple{}) = ()
@inline function _tuple_ownership_roots(values::Tuple)
    return (
        ensemble_ownership_roots(first(values))...,
        _tuple_ownership_roots(Base.tail(values))...,
    )
end

@inline ensemble_ownership_roots(interface::CompositeSimulationInterface) =
    _tuple_ownership_roots(interface.interfaces)

@inline function _shared_arm_ownership_roots(arm::SharedOpticalArm)
    channel_roots = map(channel -> (channel.wfs, channel.detector), arm.wfs_channels)
    return (
        ensemble_source_ownership_roots(arm.source)...,
        Iterators.flatten(channel_roots)...,
        arm.science_detectors...,
    )
end

@inline _shared_arms_ownership_roots(::Tuple{}) = ()
@inline function _shared_arms_ownership_roots(arms::Tuple)
    return (
        _shared_arm_ownership_roots(first(arms))...,
        _shared_arms_ownership_roots(Base.tail(arms))...,
    )
end

@inline function ensemble_ownership_roots(runtime::SharedOpticalRuntime)
    return (
        ensemble_ownership_roots(runtime.runtime)...,
        _shared_arms_ownership_roots(runtime.arms)...,
    )
end

@inline ensemble_ownership_roots(scenario::ControlLoopScenario) =
    ensemble_ownership_roots(control_loop_boundary(scenario))

@inline _same_ownership_root(::Nothing, other) = false
@inline _same_ownership_root(root, ::Nothing) = false
@inline _same_ownership_root(::Nothing, ::Nothing) = false
@inline _same_ownership_root(root, other) = root === other

function _reject_shared_ownership!(members::Tuple)
    roots = map(ensemble_ownership_roots, members)
    @inbounds for i in eachindex(roots)
        @inbounds for j in (i + 1):length(roots)
            @inbounds for left in roots[i]
                @inbounds for right in roots[j]
                    _same_ownership_root(left, right) &&
                        throw(InvalidConfiguration("ensemble members $i and $j share mutable simulation state ($(typeof(left))); construct independent plants before selecting a parallel scheduler"))
                end
            end
        end
    end
    return nothing
end

@inline validate_ensemble_policy!(::AbstractExecutionPolicy) = nothing

function validate_ensemble_policy!(::DeterministicExecution)
    Threads.nthreads() == 1 ||
        throw(InvalidConfiguration("DeterministicExecution requires starting Julia with JULIA_NUM_THREADS=1"))
    set_fft_provider_threads!(1)
    BLAS.set_num_threads(1)
    return nothing
end

init_ensemble_scheduler(::AbstractExecutionPolicy, members::Tuple) = nothing

function SimulationEnsemble(members::Vararg{AbstractControlSimulation,N};
    policy::AbstractExecutionPolicy=SequentialExecution()) where {N}
    N > 0 || throw(InvalidConfiguration("SimulationEnsemble requires at least one member"))
    validate_ensemble_policy!(policy)
    _reject_shared_ownership!(members)
    state = init_ensemble_scheduler(policy, members)
    return SimulationEnsemble{typeof(members),typeof(policy),typeof(state)}(
        members,
        policy,
        state,
    )
end

function SimulationEnsemble(members::Tuple;
    policy::AbstractExecutionPolicy=SequentialExecution())
    _validate_ensemble_member_types(members)
    return SimulationEnsemble(members...; policy=policy)
end

@inline function _run_ensemble_member!(f!, member)
    f!(member)
    return member
end

@inline _run_ensemble_sequential!(f!, ::Tuple{}) = ()
@inline function _run_ensemble_sequential!(f!, members::Tuple)
    return (
        _run_ensemble_member!(f!, first(members)),
        _run_ensemble_sequential!(f!, Base.tail(members))...,
    )
end

@inline execute_ensemble!(::SequentialExecution, state, f!, members::Tuple) =
    _run_ensemble_sequential!(f!, members)
@inline execute_ensemble!(::DeterministicExecution, state, f!, members::Tuple) =
    _run_ensemble_sequential!(f!, members)

function execute_ensemble!(::ThreadedExecution, state, f!, members::Tuple)
    Threads.nthreads() > 1 || return _run_ensemble_sequential!(f!, members)
    Threads.@threads :dynamic for i in eachindex(members)
        @inbounds _run_ensemble_member!(f!, members[i])
    end
    return members
end

function execute_ensemble!(::BackendStreamExecution, state, f!, members::Tuple)
    throw(UnsupportedAlgorithm("BackendStreamExecution is a backend-specific branch policy, not an ensemble scheduler"))
end

function execute_ensemble!(::AcceleratedKernelsExecution, state, f!, members::Tuple)
    throw(UnsupportedAlgorithm("AcceleratedKernelsExecution requires loading AcceleratedKernels.jl"))
end

function execute_ensemble!(::DaggerExecution, state, f!, members::Tuple)
    throw(UnsupportedAlgorithm("DaggerExecution requires loading Dagger.jl"))
end

"""
    run_ensemble!(f!, ensemble)

Apply one mutating coarse operation to every ensemble member according to the
selected policy. The operation must return normally only after its member is
safe for the next ensemble operation.
"""
function run_ensemble!(f!, ensemble::SimulationEnsemble)
    ensemble.members = execute_ensemble!(
        ensemble.policy,
        ensemble.scheduler_state,
        f!,
        ensemble.members,
    )
    return ensemble
end

function prepare!(ensemble::SimulationEnsemble)
    return run_ensemble!(prepare!, ensemble)
end

function sense!(ensemble::SimulationEnsemble)
    return run_ensemble!(sense!, ensemble)
end

function step!(ensemble::SimulationEnsemble)
    return run_ensemble!(step!, ensemble)
end

function synchronize_runtime!(ensemble::SimulationEnsemble)
    return run_ensemble!(synchronize_runtime!, ensemble)
end

@inline simulation_interface(ensemble::SimulationEnsemble) =
    map(simulation_interface, ensemble.members)
@inline simulation_readout(ensemble::SimulationEnsemble) =
    map(simulation_readout, ensemble.members)
@inline ensemble_readouts(ensemble::SimulationEnsemble) =
    simulation_readout(ensemble)

@inline simulation_command(ensemble::SimulationEnsemble) =
    map(simulation_command, ensemble.members)
@inline simulation_slopes(ensemble::SimulationEnsemble) =
    map(simulation_slopes, ensemble.members)
@inline simulation_wfs_frame(ensemble::SimulationEnsemble) =
    map(simulation_wfs_frame, ensemble.members)
@inline simulation_science_frame(ensemble::SimulationEnsemble) =
    map(simulation_science_frame, ensemble.members)
@inline simulation_wfs_metadata(ensemble::SimulationEnsemble) =
    map(simulation_wfs_metadata, ensemble.members)
@inline simulation_science_metadata(ensemble::SimulationEnsemble) =
    map(simulation_science_metadata, ensemble.members)
@inline simulation_grouped_wfs_stack(ensemble::SimulationEnsemble) =
    map(simulation_grouped_wfs_stack, ensemble.members)
@inline simulation_grouped_science_stack(ensemble::SimulationEnsemble) =
    map(simulation_grouped_science_stack, ensemble.members)
