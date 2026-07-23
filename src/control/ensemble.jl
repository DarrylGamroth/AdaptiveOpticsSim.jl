"""
    SimulationEnsemble(members...; policy=SequentialExecution())

Own independent coarse work units and apply a caller-supplied operation to
them according to an execution policy. Members may be simulations, prepared
plants, calibration jobs, or sweep items; the ensemble does not prescribe a
step-wise simulation interface.

Construction rejects shared ownership roots so parallel policies cannot
silently race on mutable state.
"""
mutable struct SimulationEnsemble{
    M<:Tuple,
    P<:AbstractExecutionPolicy,
    S,
}
    members::M
    policy::P
    scheduler_state::S
end

@inline ensemble_members(ensemble::SimulationEnsemble) = ensemble.members
@inline execution_policy(ensemble::SimulationEnsemble) = ensemble.policy

"""
    ensemble_ownership_roots(member)

Return mutable state owners that must not be shared by independently scheduled
ensemble members. Wrapper types should extend this function and return the
underlying mutable roots they own.
"""
@inline ensemble_ownership_roots(member) = (member,)

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
                        throw(InvalidConfiguration(
                            "ensemble members $i and $j share mutable state " *
                            "($(typeof(left))); construct independent work " *
                            "units before selecting a parallel scheduler"))
                end
            end
        end
    end
    return nothing
end

@inline validate_ensemble_policy!(::AbstractExecutionPolicy) = nothing

function validate_ensemble_policy!(::DeterministicExecution)
    Threads.nthreads() == 1 ||
        throw(InvalidConfiguration(
            "DeterministicExecution requires JULIA_NUM_THREADS=1"))
    set_fft_provider_threads!(1)
    BLAS.set_num_threads(1)
    return nothing
end

init_ensemble_scheduler(::AbstractExecutionPolicy, members::Tuple) = nothing

function SimulationEnsemble(members...;
    policy::AbstractExecutionPolicy=SequentialExecution())
    isempty(members) &&
        throw(InvalidConfiguration(
            "SimulationEnsemble requires at least one member"))
    validate_ensemble_policy!(policy)
    _reject_shared_ownership!(members)
    state = init_ensemble_scheduler(policy, members)
    return SimulationEnsemble{
        typeof(members),typeof(policy),typeof(state),
    }(members, policy, state)
end

function SimulationEnsemble(members::Tuple;
    policy::AbstractExecutionPolicy=SequentialExecution())
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
    Threads.nthreads() > 1 ||
        return _run_ensemble_sequential!(f!, members)
    Threads.@threads :dynamic for i in eachindex(members)
        @inbounds _run_ensemble_member!(f!, members[i])
    end
    return members
end

function execute_ensemble!(::BackendStreamExecution, state, f!,
    members::Tuple)
    throw(UnsupportedAlgorithm(
        "BackendStreamExecution is a backend-specific branch policy, " *
        "not an ensemble scheduler"))
end

function execute_ensemble!(::AcceleratedKernelsExecution, state, f!,
    members::Tuple)
    throw(UnsupportedAlgorithm(
        "AcceleratedKernelsExecution requires loading AcceleratedKernels.jl"))
end

function execute_ensemble!(::DaggerExecution, state, f!, members::Tuple)
    throw(UnsupportedAlgorithm(
        "DaggerExecution requires loading Dagger.jl"))
end

"""
    run_ensemble!(f!, ensemble)

Apply one mutating coarse operation to every ensemble member. The operation
must return normally only after the member is safe for the next ensemble
operation.
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
