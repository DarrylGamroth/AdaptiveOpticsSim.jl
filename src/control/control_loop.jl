#
# Julia-native control-loop orchestration layer
#
# This file introduces typed scenario/config objects for runtime composition
# without switching the package to config-file-driven control.
#

abstract type AbstractControlLoopConfig end

struct ControlLoopBranch{SIM,REC,WD,SD,RNG,PR,CL,B<:AbstractArrayBackend}
    label::Symbol
    simulation::SIM
    reconstructor::REC
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    outputs::PR
    command_layout::CL
end

function ControlLoopBranch(label::Symbol, simulation, reconstructor;
    wfs_detector=nothing,
    science_detector=nothing,
    rng::AbstractRNG=Random.default_rng(),
    outputs=nothing,
    command_layout=nothing)
    selector = require_same_backend(simulation, wfs_detector, science_detector)
    return ControlLoopBranch{
        typeof(simulation),
        typeof(reconstructor),
        typeof(wfs_detector),
        typeof(science_detector),
        typeof(rng),
        typeof(outputs),
        typeof(command_layout),
        typeof(selector),
    }(
        label,
        simulation,
        reconstructor,
        wfs_detector,
        science_detector,
        rng,
        outputs,
        command_layout,
    )
end

@inline backend(::ControlLoopBranch{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

struct SingleControlLoopConfig{RP<:AbstractRuntimeProfile,PR<:RuntimeOutputRequirements,T<:AbstractFloat} <: AbstractControlLoopConfig
    name::Symbol
    branch_label::Symbol
    profile::RP
    outputs::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function SingleControlLoopConfig(;
    name::Symbol=:control_loop_runtime,
    branch_label::Symbol=:main,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    outputs::RuntimeOutputRequirements=default_runtime_outputs(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
)
    T = typeof(float(control_sign))
    return SingleControlLoopConfig{
        typeof(profile),
        typeof(outputs),
        T,
    }(
        name,
        branch_label,
        profile,
        outputs,
        latency,
        T(control_sign),
        Int(science_zero_padding),
    )
end

struct GroupedControlLoopConfig{N,RP<:AbstractRuntimeProfile,PR<:GroupedRuntimeOutputRequirements,T<:AbstractFloat} <: AbstractControlLoopConfig
    name::Symbol
    branch_labels::NTuple{N,Symbol}
    profile::RP
    outputs::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function GroupedControlLoopConfig(branch_labels::NTuple{N,Symbol};
    name::Symbol=:grouped_control_loop_runtime,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    outputs::GroupedRuntimeOutputRequirements=default_grouped_runtime_outputs(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
) where {N}
    T = typeof(float(control_sign))
    return GroupedControlLoopConfig{
        N,
        typeof(profile),
        typeof(outputs),
        T,
    }(
        name,
        branch_labels,
        profile,
        outputs,
        latency,
        T(control_sign),
        Int(science_zero_padding),
    )
end

function GroupedControlLoopConfig(branches::Vararg{ControlLoopBranch,N}; kwargs...) where {N}
    return GroupedControlLoopConfig(ntuple(i -> branches[i].label, N); kwargs...)
end

struct ControlLoopScenario{CFG<:AbstractControlLoopConfig,BR,BD,B<:AbstractArrayBackend} <: AbstractControlSimulation
    config::CFG
    branches::BR
    boundary::BD
end


@inline control_loop_config(scenario::ControlLoopScenario) = scenario.config
@inline control_loop_boundary(scenario::ControlLoopScenario) = scenario.boundary
@inline control_loop_name(config::SingleControlLoopConfig) = config.name
@inline control_loop_name(config::GroupedControlLoopConfig) = config.name
@inline control_loop_name(scenario::ControlLoopScenario) = control_loop_name(control_loop_config(scenario))
@inline control_loop_branch_labels(config::SingleControlLoopConfig) = (config.branch_label,)
@inline control_loop_branch_labels(config::GroupedControlLoopConfig) = config.branch_labels
@inline control_loop_branch_labels(scenario::ControlLoopScenario) = control_loop_branch_labels(control_loop_config(scenario))
@inline runtime_profile(scenario::ControlLoopScenario) = control_loop_config(scenario).profile
@inline runtime_latency(scenario::ControlLoopScenario) = control_loop_config(scenario).latency

@inline function _branch_runtime_outputs(config::SingleControlLoopConfig, branch::ControlLoopBranch)
    isnothing(branch.outputs) || return branch.outputs
    return config.outputs
end

@inline function _branch_runtime_outputs(config::GroupedControlLoopConfig, branch::ControlLoopBranch)
    isnothing(branch.outputs) || return branch.outputs
    return RuntimeOutputRequirements(
        slopes=true,
        wfs_pixels=(config.outputs.wfs_frames || config.outputs.wfs_stack) && !isnothing(branch.wfs_detector),
        science_pixels=(config.outputs.science_frames || config.outputs.science_stack) && !isnothing(branch.science_detector),
    )
end

@inline function _build_control_loop_runtime(config::AbstractControlLoopConfig, branch::ControlLoopBranch)
    return ClosedLoopRuntime(
        branch.simulation,
        branch.reconstructor;
        wfs_detector=branch.wfs_detector,
        science_detector=branch.science_detector,
        rng=branch.rng,
        profile=config.profile,
        outputs=_branch_runtime_outputs(config, branch),
        latency=config.latency,
        control_sign=config.control_sign,
        science_zero_padding=config.science_zero_padding,
        command_layout=branch.command_layout,
    )
end

function build_control_loop_scenario(config::SingleControlLoopConfig, branch::ControlLoopBranch)
    branch.label == config.branch_label ||
        throw(InvalidConfiguration("single control-loop branch label $(branch.label) must match config branch_label $(config.branch_label)"))
    runtime = _build_control_loop_runtime(config, branch)
    boundary = SimulationInterface(runtime)
    selector = backend(branch)
    return ControlLoopScenario{
        typeof(config),
        Tuple{typeof(branch)},
        typeof(boundary),
        typeof(selector),
    }(
        config,
        (branch,),
        boundary,
    )
end

function build_control_loop_scenario(config::GroupedControlLoopConfig{N}, branches::Vararg{ControlLoopBranch,N}) where {N}
    selector = require_same_backend(branches...)
    labels = ntuple(i -> branches[i].label, N)
    labels == config.branch_labels ||
        throw(InvalidConfiguration("grouped control-loop branch labels $labels must match config branch_labels $(config.branch_labels)"))
    runtimes = ntuple(i -> _build_control_loop_runtime(config, branches[i]), N)
    boundary = CompositeSimulationInterface(runtimes...; outputs=config.outputs)
    return ControlLoopScenario{
        typeof(config),
        typeof(branches),
        typeof(boundary),
        typeof(selector),
    }(
        config,
        branches,
        boundary,
    )
end

@inline backend(::ControlLoopScenario{<:Any,<:Any,<:Any,B}) where {B} = B()
simulation_interface(scenario::ControlLoopScenario) = control_loop_boundary(scenario)

readout(sim::AbstractControlSimulation) = simulation_readout(sim)
readout(runtime::ClosedLoopRuntime) = simulation_readout(runtime)
readout(interface::SimulationInterface) = simulation_readout(interface)
readout(interface::CompositeSimulationInterface) = simulation_readout(interface)

function prepare!(scenario::ControlLoopScenario)
    prepare!(control_loop_boundary(scenario))
    return scenario
end

function sense!(scenario::ControlLoopScenario)
    sense!(control_loop_boundary(scenario))
    return scenario
end

function step!(scenario::ControlLoopScenario)
    step!(control_loop_boundary(scenario))
    return scenario
end

function snapshot_outputs!(scenario::ControlLoopScenario)
    snapshot_outputs!(control_loop_boundary(scenario))
    return scenario
end

@inline set_command!(scenario::ControlLoopScenario, command::AbstractVector) = set_command!(control_loop_boundary(scenario), command)
@inline runtime_phase_timing(scenario::ControlLoopScenario; kwargs...) = runtime_phase_timing(control_loop_boundary(scenario); kwargs...)

function command_layout(scenario::ControlLoopScenario)
    labels = control_loop_branch_labels(scenario)
    return scenario_command_layout(control_loop_boundary(scenario), labels)
end

scenario_command_layout(boundary::SimulationInterface, labels) = command_layout(boundary)

function scenario_command_layout(boundary::CompositeSimulationInterface, labels)
    segments = ntuple(i -> RuntimeCommandSegment(labels[i], i == 1 ? 1 : 1 + sum(length(boundary.interfaces[j].command) for j in 1:(i - 1)),
        length(boundary.interfaces[i].command)), length(labels))
    return RuntimeCommandLayout(segments)
end

@inline function branch_command_layout(scenario::ControlLoopScenario, label::Symbol)
    boundary = control_loop_boundary(scenario)
    labels = control_loop_branch_labels(scenario)
    idx = findfirst(==(label), labels)
    isnothing(idx) && throw(InvalidConfiguration("runtime scenario does not contain branch $(label)"))
    return branch_command_layout(boundary, idx)
end

branch_command_layout(boundary::SimulationInterface, idx::Int) = command_layout(boundary)
branch_command_layout(boundary::CompositeSimulationInterface, idx::Int) = command_layout(boundary.interfaces[idx])

@inline function branch_command_layouts(scenario::ControlLoopScenario)
    labels = control_loop_branch_labels(scenario)
    layouts = map(label -> branch_command_layout(scenario, label), labels)
    return NamedTuple{labels}(Tuple(layouts))
end

@inline function update_command!(scenario::ControlLoopScenario, command::NamedTuple)
    boundary = control_loop_boundary(scenario)
    labels = control_loop_branch_labels(scenario)
    update_scenario_command!(boundary, labels, command)
    snapshot_outputs!(scenario)
    return simulation_command(scenario)
end

@inline function set_command!(scenario::ControlLoopScenario, command::NamedTuple)
    boundary = control_loop_boundary(scenario)
    labels = control_loop_branch_labels(scenario)
    set_scenario_command!(boundary, labels, command)
    snapshot_outputs!(scenario)
    return simulation_command(scenario)
end

update_scenario_command!(boundary::SimulationInterface, labels, command::NamedTuple) =
    update_command!(boundary, command)

function update_scenario_command!(boundary::CompositeSimulationInterface, labels, command::NamedTuple)
    all(label -> label in labels, Tuple(keys(command))) ||
        throw(InvalidConfiguration("structured scenario command keys $(Tuple(keys(command))) must be a subset of branch labels $labels"))
    @inbounds for (idx, label) in pairs(labels)
        if hasproperty(command, label)
            apply_branch_update_command!(boundary.interfaces[idx], getproperty(command, label), label)
        end
    end
    return boundary
end

set_scenario_command!(boundary::SimulationInterface, labels, command::NamedTuple) =
    set_command!(boundary, command)

function set_scenario_command!(boundary::CompositeSimulationInterface, labels, command::NamedTuple)
    Tuple(keys(command)) == labels ||
        throw(InvalidConfiguration("structured scenario command keys $(Tuple(keys(command))) must match branch labels $labels"))
    @inbounds for (idx, label) in pairs(labels)
        apply_branch_set_command!(boundary.interfaces[idx], getproperty(command, label), label)
    end
    return boundary
end

apply_branch_update_command!(interface, command::NamedTuple, label::Symbol) = update_command!(interface, command)
apply_branch_update_command!(interface, command::AbstractVector, label::Symbol) = set_command!(interface, command)
apply_branch_update_command!(interface, command, label::Symbol) =
    throw(InvalidConfiguration("structured branch command $(label) must be an AbstractVector or NamedTuple"))

apply_branch_set_command!(interface, command::NamedTuple, label::Symbol) = set_command!(interface, command)
apply_branch_set_command!(interface, command::AbstractVector, label::Symbol) = set_command!(interface, command)
apply_branch_set_command!(interface, command, label::Symbol) =
    throw(InvalidConfiguration("structured branch command $(label) must be an AbstractVector or NamedTuple"))

supports_prepared_runtime(scenario::ControlLoopScenario) =
    any(branch -> supports_prepared_runtime(branch.simulation.wfs, branch.simulation.src), scenario.branches)
supports_detector_output(scenario::ControlLoopScenario) = any(branch -> !isnothing(branch.wfs_detector) || !isnothing(branch.science_detector), scenario.branches)
supports_grouped_execution(scenario::ControlLoopScenario) = length(scenario.branches) > 1
